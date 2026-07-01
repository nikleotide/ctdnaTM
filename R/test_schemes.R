# ============================================================================
# ctdnaTM — Automated test-data generation for filtering schemes
# ============================================================================
#
# ctdna_test_scheme(scheme_name) returns a synthetic data frame designed to
# exercise every pass branch and every gate of a registered filtering scheme.
# Rows are auto-generated from the scheme's rule tree -- no hand-authoring is
# required when a new scheme is added.
#
# Each returned row has two extra columns:
#   filter_decision  ("keep" or "remove") -- the EXPECTED outcome
#   filter_reason    -- a short string describing what the row tests
#
# The function evaluates the scheme on the generated rows before returning
# and errors loudly if the engine disagrees with the expected decision.
# That self-check catches generator bugs (or unintended scheme changes).
#
# Sample-level QC columns (Sample_status, cfDNA_ng, Plasma_ml_input) are
# populated with passing values on every scheme-branch row, then a few
# extra rows are appended that specifically fail QC. After QC fires, every
# scheme branch still has at least one row alive.

# ---- Value bank ------------------------------------------------------------
#
# For each known column, provide:
#   any  — a single matching/realistic default value used when no rule
#          forces a specific value
#   pool — a vector of REALISTIC candidate values (real Guardant Health
#          field values, real gene symbols, etc.) used to find a value
#          NOT matching a constraint when the violator needs one.
#
# Adding a new column referenced by a future scheme means adding one
# entry here with `any` and `pool`. Both must be real-data-realistic;
# no synthetic placeholder strings.

.ctdna_test_value_bank <- function() list(
  Gene = list(
    any  = "TP53",
    pool = c("TP53","RB1","PTEN","APC","KRAS","EGFR","BRAF","ERBB2",
              "PIK3CA","BRCA1","BRCA2","ATM","BARD1","CDK12","PALB2",
              "MYC","JAK2","FGFR3","SMAD4","CDH1","NF1","VHL","STK11",
              "ARID1A","ARID2","SETD2","KMT2D","TSC1","TSC2","FBXW7")),
  Variant_type = list(
    any  = "SNV",
    pool = c("SNV","Indel","CNV","Fusion","LGR")),
  Somatic_status = list(
    any  = "somatic",
    pool = c("somatic","germline","somatic_putative_ch")),
  ClinVar = list(
    any  = "Pathogenic",
    pool = c("Pathogenic","Likely_pathogenic","Pathogenic/Likely_pathogenic",
              "Uncertain_significance","Benign","Likely_benign",
              "Benign/Likely_benign",
              "Conflicting_interpretations_of_pathogenicity")),
  Molecular_consequence = list(
    any  = "missense",
    pool = c("missense","nonsense","frameshift","start_lost","stop_lost",
              "splice_acceptor","splice_donor","splice_region","splice_event",
              "synonymous","intron","inframe_indel","inframe_insertion",
              "inframe_deletion","inframe_duplication","promoter")),
  CNV_type = list(
    any  = "homozygous_deletion",
    pool = c("homozygous_deletion","loh_deletion","amplification",
              "aneuploid_amplification","focal_amplification")),
  Functional_impact = list(
    any  = "deleterious",
    pool = c("deleterious","reversion","reversion_cis")),
  VAF_percentage = list(
    any  = 5,
    pool = c(0, 0.5, 5, 25, 80)),
  gnomAD_AF = list(
    any  = 0.0001,
    pool = c(0.00001, 0.0001, 0.001, 0.005, 0.02, 0.1, 0.5)),
  Copy_number = list(
    any  = 10,
    pool = c(0, 1, 2, 3, 4, 6, 10, 20)),
  Direction_a = list(
    any  = -1,
    pool = c(-1, 1)),
  Direction_b = list(
    any  = 1,
    pool = c(-1, 1)),
  Chromosome = list(
    any  = "chr5",
    pool = c("chr1","chr2","chr3","chr5","chr8","chr10","chr13","chr17","chrX")),
  Fusion_chrom_b = list(
    any  = "chr5",
    pool = c("chr1","chr2","chr3","chr5","chr8","chr10","chr13","chr17","chrX")),
  Fusion_gene_b = list(
    any  = "TP53",
    pool = c("TP53","BRAF","BRCA1","BRCA2","APC","RB1","KRAS","MYC",
              "EGFR","PIK3CA")),
  Mut_aa = list(
    any  = "p.R175H",
    pool = c("p.R175H","p.E545K","p.G12D","p.V600E","p.K3326*","p.X100*")),
  Mutant_allele_status = list(
    any  = "biallelic",
    pool = c("biallelic","monoallelic","unknown")))

# Look up a baseline value for a column. Unknown columns get NA.
.ctdna_test_default <- function(col_name) {
  bank <- .ctdna_test_value_bank()
  if (!is.null(bank[[col_name]])) return(bank[[col_name]]$any)
  NA
}

# Find a REALISTIC value in `pool` that does NOT match `forbidden`.
# Used by the violator to pick a near-miss value.
.ctdna_test_pick_miss <- function(col_name, forbidden) {
  bank <- .ctdna_test_value_bank()
  entry <- bank[[col_name]]
  if (is.null(entry) || is.null(entry$pool)) {
    # No bank entry — fall back to a generic but still realistic-looking
    # value. We choose a real-data-like placeholder rather than a
    # synthetic sentinel.
    return("OTHER")
  }
  pool <- entry$pool
  # Coerce forbidden to comparable form
  candidates <- if (is.numeric(pool)) {
    forb_num <- suppressWarnings(as.numeric(as.character(forbidden)))
    setdiff(pool, forb_num)
  } else {
    forb_str <- tolower(trimws(as.character(forbidden)))
    pool[!tolower(trimws(as.character(pool))) %in% forb_str]
  }
  if (length(candidates) == 0L) {
    # Pool is fully covered by forbidden — pick the first pool value
    # anyway (caller must handle the edge case)
    return(pool[[1]])
  }
  candidates[[1]]
}

# ---- Tree walk helpers -----------------------------------------------------

# Return character vector of all column names referenced anywhere in a tree.
.ctdna_collect_cols <- function(node) {
  if (is.null(node)) return(character(0))
  if (inherits(node, "ctdna_filtering_scheme"))
    return(unique(c(.ctdna_collect_cols(node$scope),
                     .ctdna_collect_cols(node$root))))
  if (inherits(node, "ctdna_rule_atom")) {
    cols <- names(node$conditions)
    # column-vs-column constraints reference a partner column via value
    for (k in names(node$conditions)) {
      v <- node$conditions[[k]]
      if (is.list(v) && all(c("op","value") %in% names(v)) &&
          v$op %in% c("==col","!=col"))
        cols <- c(cols, v$value)
    }
    return(unique(cols))
  }
  if (inherits(node, "ctdna_rule_combinator")) {
    if (node$type == "not") return(.ctdna_collect_cols(node$child))
    return(unique(unlist(lapply(node$children, .ctdna_collect_cols))))
  }
  character(0)
}

# Return a flat list of all ctdna_rule_atom nodes in a tree.
# Used by the drop-row generator to enumerate atoms in a branch
# so we can try violating each one until we find a genuine drop.
.ctdna_collect_atoms <- function(node) {
  if (is.null(node)) return(list())
  if (inherits(node, "ctdna_filtering_scheme"))
    return(c(.ctdna_collect_atoms(node$scope),
              .ctdna_collect_atoms(node$root)))
  if (inherits(node, "ctdna_rule_atom"))
    return(list(node))
  if (inherits(node, "ctdna_rule_combinator")) {
    if (node$type == "not") return(.ctdna_collect_atoms(node$child))
    return(unlist(lapply(node$children, .ctdna_collect_atoms),
                    recursive = FALSE))
  }
  list()
}

# Return a one-row data frame with default ("any") values for every column.
.ctdna_test_blank_row <- function(cols) {
  vals <- lapply(cols, .ctdna_test_default)
  names(vals) <- cols
  as.data.frame(vals, stringsAsFactors = FALSE)
}

# ---- Atom satisfiers -------------------------------------------------------
#
# Each function takes a one-row data frame and a constraint, and returns a
# modified row that satisfies (or violates) the constraint. The walker
# composes these for combinators.

# Modify `row` in place to make atom MATCH.
.ctdna_test_satisfy_atom <- function(row, atom) {
  for (col_name in names(atom$conditions)) {
    constraint <- atom$conditions[[col_name]]
    if (is.list(constraint) && all(c("op","value") %in% names(constraint))) {
      op <- constraint$op; val <- constraint$value
      if (op == "==col") {
        # Force `col_name` to equal the value of column `val` on the same row
        row[[col_name]] <- row[[val]]
      } else if (op == "!=col") {
        # Make different — pick a realistic value from the bank that
        # differs from the partner column's current value
        other_val <- row[[val]]
        row[[col_name]] <- .ctdna_test_pick_miss(col_name, other_val)
      } else if (op %in% c("<","<="))   row[[col_name]] <- val - 1e-6
       else if (op %in% c(">",">="))    row[[col_name]] <- val + 1e-6
       else if (op == "==")             row[[col_name]] <- val
       else if (op == "!=")             {
         # Use a realistic non-matching value from the bank
         row[[col_name]] <- .ctdna_test_pick_miss(col_name, val)
       } else if (op == "%in%")          row[[col_name]] <- val[1]
       else if (op == "not_in")          {
         # Need a realistic value NOT in val
         row[[col_name]] <- .ctdna_test_pick_miss(col_name, val)
       } else                            row[[col_name]] <- val[1]
    } else {
      # Plain value(s) — pick the first
      row[[col_name]] <- constraint[[1]]
    }
  }
  row
}

# Modify `row` to make atom NOT MATCH (force a near-miss).
# `protect` = character vector of columns the violator should AVOID
# touching (because higher-level gates already populated them).
# All "miss" values are realistic Guardant-like values from the bank.
.ctdna_test_violate_atom <- function(row, atom, protect = character(0)) {
  # Prefer columns NOT in protect set; fall back to first if all protected
  all_keys <- names(atom$conditions)
  candidates <- setdiff(all_keys, protect)
  k <- if (length(candidates) > 0) candidates[[1]] else all_keys[[1]]
  constraint <- atom$conditions[[k]]
  if (is.list(constraint) && all(c("op","value") %in% names(constraint))) {
    op <- constraint$op; val <- constraint$value
    if (op == "==col") {
      # Make col != other column: pick a realistic value that differs
      other_val <- row[[val]]
      row[[k]] <- .ctdna_test_pick_miss(k, other_val)
    } else if (op == "!=col") {
      row[[k]] <- row[[val]]  # make equal (violates !=)
    } else if (op %in% c("<","<="))   row[[k]] <- val + 1
     else if (op %in% c(">",">="))    row[[k]] <- val - 1
     else if (op == "==")              {
       if (is.numeric(val)) {
         row[[k]] <- val + 999
       } else {
         row[[k]] <- .ctdna_test_pick_miss(k, val)
       }
     }
     else if (op == "!=")              row[[k]] <- val
     else if (op == "%in%")            row[[k]] <- .ctdna_test_pick_miss(k, val)
     else if (op == "not_in")          row[[k]] <- val[1]
     else                              row[[k]] <- .ctdna_test_pick_miss(k, val)
  } else {
    # Plain value(s) — pick a realistic non-matching value
    row[[k]] <- .ctdna_test_pick_miss(k, constraint)
  }
  row
}

# ---- Tree satisfiers (recursive) -------------------------------------------

# Modify `row` to make `node` evaluate to TRUE.
.ctdna_test_satisfy <- function(row, node, protect = character(0)) {
  if (inherits(node, "ctdna_filtering_scheme"))
    return(.ctdna_test_satisfy(row, node$root, protect))
  if (inherits(node, "ctdna_rule_atom"))
    return(.ctdna_test_satisfy_atom(row, node))
  if (inherits(node, "ctdna_rule_combinator")) {
    if (node$type == "not")  return(.ctdna_test_violate(row, node$child, protect))
    if (node$type == "and")  {
      for (ch in node$children) row <- .ctdna_test_satisfy(row, ch, protect)
      return(row)
    }
    if (node$type == "or")   {
      # Satisfy the first child (any branch suffices)
      return(.ctdna_test_satisfy(row, node$children[[1]], protect))
    }
  }
  row
}

# Modify `row` to make `node` evaluate to FALSE.
.ctdna_test_violate <- function(row, node, protect = character(0)) {
  if (inherits(node, "ctdna_filtering_scheme"))
    return(.ctdna_test_violate(row, node$root, protect))
  if (inherits(node, "ctdna_rule_atom"))
    return(.ctdna_test_violate_atom(row, node, protect))
  if (inherits(node, "ctdna_rule_combinator")) {
    if (node$type == "not")  return(.ctdna_test_satisfy(row, node$child, protect))
    if (node$type == "and")  {
      # Violate at least one child — prefer one whose columns aren't protected
      for (ch in node$children) {
        ch_cols <- .ctdna_collect_cols(ch)
        if (length(setdiff(ch_cols, protect)) > 0) {
          return(.ctdna_test_violate(row, ch, protect))
        }
      }
      # Fall back to first
      return(.ctdna_test_violate(row, node$children[[1]], protect))
    }
    if (node$type == "or")   {
      # Violate ALL children
      for (ch in node$children) row <- .ctdna_test_violate(row, ch, protect)
      return(row)
    }
  }
  row
}

# ---- Branch identification -------------------------------------------------
#
# The "pass branches" of a scheme are the children of the top-level OR node
# (if one exists). Each branch is a distinct way for a row to be kept.
# When there is no OR at the top level, the entire root is treated as a
# single pass branch.

.ctdna_test_find_pass_branches <- function(root) {
  if (inherits(root, "ctdna_rule_combinator") && root$type == "and") {
    # Look for an `or` child — that's the disjunction of pass branches
    or_child <- NULL
    gates <- list()
    for (ch in root$children) {
      if (inherits(ch, "ctdna_rule_combinator") && ch$type == "or" &&
          is.null(or_child)) {
        or_child <- ch
      } else {
        gates <- c(gates, list(ch))
      }
    }
    if (!is.null(or_child)) {
      return(list(branches = or_child$children, gates = gates))
    }
    # No OR — every AND child is itself a gate; treat the whole root as one branch
    return(list(branches = list(root), gates = list()))
  }
  if (inherits(root, "ctdna_rule_combinator") && root$type == "or") {
    return(list(branches = root$children, gates = list()))
  }
  # Atomic root or NOT root
  list(branches = list(root), gates = list())
}

# Short label for a node — used to make filter_reason readable.
# Combinators show a brief summary of their first 1-2 children.
.ctdna_test_label <- function(node, max_len = 80) {
  if (is.null(node)) return("(null)")
  if (inherits(node, "ctdna_filtering_scheme")) return(node$name %||% "scheme")
  if (inherits(node, "ctdna_rule_atom")) {
    parts <- character(0)
    for (k in names(node$conditions)) {
      v <- node$conditions[[k]]
      if (is.list(v) && all(c("op","value") %in% names(v))) {
        val_str <- paste(format(v$value), collapse = ",")
        if (nchar(val_str) > 20) val_str <- paste0(substr(val_str, 1, 17), "...")
        parts <- c(parts, sprintf("%s %s %s", k, v$op, val_str))
      } else {
        val_str <- paste(format(v), collapse = ",")
        if (nchar(val_str) > 20) val_str <- paste0(substr(val_str, 1, 17), "...")
        parts <- c(parts, sprintf("%s=%s", k, val_str))
      }
    }
    s <- paste(parts, collapse = " & ")
    if (nchar(s) > max_len) s <- paste0(substr(s, 1, max_len - 3), "...")
    return(s)
  }
  if (inherits(node, "ctdna_rule_combinator")) {
    if (node$type == "not")
      return(paste0("NOT(", .ctdna_test_label(node$child, max_len - 5), ")"))
    # Show first 2-3 children inline
    op <- if (node$type == "and") " AND " else " OR "
    child_strs <- vapply(node$children,
                          function(c) .ctdna_test_label(c, max_len = 30),
                          character(1))
    s <- paste(child_strs, collapse = op)
    if (nchar(s) > max_len) s <- paste0(substr(s, 1, max_len - 3), "...")
    return(s)
  }
  "?"
}

# ---- Per-row explainer (used when df is user-provided) --------------------
#
# Given ONE row of a data frame and a scheme, return a list with $decision
# ("keep" or "remove") and $reason (a short string). The decision mirrors
# what ctdna_variant_filter(..., filter_scheme = scheme_name) would do for
# this row (QC -> scope -> scheme rules; out-of-scope rows route to
# scheme_basic). The reason explains why.

.ctdna_test_explain_row <- function(row_df, scheme, scheme_name) {
  # row_df is expected to be a one-row data frame

  # ---- 1. Sample-level QC ----
  if (isTRUE(.o("sample_qc_filter"))) {
    status_col <- .o("sample_qc_status_col")
    cfdna_col  <- .o("sample_qc_cfdna_col")
    plasma_col <- .o("sample_qc_plasma_col")
    if (!is.null(status_col) && status_col %in% names(row_df)) {
      pass_val <- .o("sample_qc_status_pass")
      v <- as.character(row_df[[status_col]])
      if (is.na(v) || v != pass_val)
        return(list(decision = "remove",
                     reason = sprintf("QC: %s='%s' (expected '%s')",
                                       status_col, v, pass_val)))
    }
    if (!is.null(cfdna_col) && cfdna_col %in% names(row_df)) {
      min_val <- .o("sample_qc_cfdna_min")
      v <- suppressWarnings(as.numeric(row_df[[cfdna_col]]))
      if (is.na(v) || v < min_val)
        return(list(decision = "remove",
                     reason = sprintf("QC: %s=%s (< %s)",
                                       cfdna_col, v, min_val)))
    }
    if (!is.null(plasma_col) && plasma_col %in% names(row_df)) {
      min_val <- .o("sample_qc_plasma_min")
      v <- suppressWarnings(as.numeric(row_df[[plasma_col]]))
      if (is.na(v) || v < min_val)
        return(list(decision = "remove",
                     reason = sprintf("QC: %s=%s (< %s)",
                                       plasma_col, v, min_val)))
    }
  }

  # ---- 2. Scope check ----
  # If the scheme has a scope rule and this row is out of scope, the
  # composed filter_apply routes the row to scheme_basic. Mirror that.
  in_scope <- TRUE
  if (!is.null(scheme$scope)) {
    sm <- tryCatch(.eval_rule(scheme$scope, row_df),
                    error = function(e) FALSE)
    sm[is.na(sm)] <- FALSE
    in_scope <- isTRUE(sm[1])
    if (!in_scope) {
      basic_mask <- tryCatch(.eval_rule(scheme_basic(), row_df),
                              error = function(e) FALSE)
      basic_mask[is.na(basic_mask)] <- FALSE
      if (isTRUE(basic_mask[1])) {
        return(list(decision = "keep",
                     reason = "out of scheme scope; passed scheme_basic"))
      } else {
        return(list(decision = "remove",
                     reason = "out of scheme scope; failed scheme_basic"))
      }
    }
  }

  # ---- 3. In-scope: evaluate full scheme. The engine treats missing
  # columns as FALSE-matching rules. We mirror that and then explain.
  full_mask <- tryCatch(.eval_rule(scheme, row_df),
                         error = function(e) FALSE)
  full_mask[is.na(full_mask)] <- FALSE
  if (isTRUE(full_mask[1])) {
    # Find which branch matched (for the reason)
    parts <- .ctdna_test_find_pass_branches(scheme$root)
    branches <- parts$branches
    for (i in seq_along(branches)) {
      bm <- tryCatch(.eval_rule(branches[[i]], row_df),
                      error = function(e) FALSE)
      bm[is.na(bm)] <- FALSE
      if (isTRUE(bm[1]))
        return(list(decision = "keep",
                     reason = sprintf("matched pass branch %d: %s", i,
                                       .ctdna_test_label(branches[[i]],
                                                         max_len = 80))))
    }
    # Shouldn't reach here if full_mask was TRUE; fallback
    return(list(decision = "keep", reason = "scheme matched"))
  }

  # ---- 4. Engine dropped it. Walk to find which gate/branch caused it.
  parts <- .ctdna_test_find_pass_branches(scheme$root)
  branches <- parts$branches
  gates    <- parts$gates

  # Note any missing columns the scheme would have looked at -- helps
  # diagnose when the drop reason is "column wasn't there to check".
  needed <- .ctdna_collect_cols(scheme)
  missing_cols <- setdiff(needed, names(row_df))
  miss_suffix <- if (length(missing_cols) > 0L)
    sprintf(" [columns missing from input: %s]",
             paste(missing_cols, collapse = ", ")) else ""

  for (g in gates) {
    gm <- tryCatch(.eval_rule(g, row_df), error = function(e) FALSE)
    gm[is.na(gm)] <- FALSE
    if (!isTRUE(gm[1]))
      return(list(decision = "remove",
                   reason = paste0(sprintf("gate failed: %s",
                                     .ctdna_test_label(g, max_len = 80)),
                                     miss_suffix)))
  }

  list(decision = "remove",
        reason = paste0("no pass branch matched", miss_suffix))
}

# ---- Main generator --------------------------------------------------------

#' Generate or label test data for a registered filtering scheme
#'
#' Two modes of operation, dispatched on the \code{df} argument:
#'
#' \strong{Mode 1 -- simulate (df = NULL, default):} auto-generates a
#' small synthetic data frame that exercises every pass branch and every
#' gate of the named scheme. Rows are constructed from the scheme's rule
#' tree using realistic Guardant Health-like values. Self-verifies
#' against the engine before returning.
#'
#' \strong{Mode 2 -- label (df = user data frame):} takes the user's
#' data frame and adds two columns describing what would happen to each
#' row under the named scheme (combined with sample-level QC). No rows
#' are dropped; the original data is returned with two extra columns.
#'
#' In both modes the returned frame has two added columns:
#' \itemize{
#'   \item \code{filter_decision} -- "keep" or "remove"
#'   \item \code{filter_reason} -- short string explaining why
#' }
#'
#' \strong{How the simulator works (mode 1):} the function walks the
#' scheme's rule tree, identifies gates (top-level AND children) and
#' pass branches (children of the top-level OR node, if present), and
#' generates exactly 2N + M + 3 rows where N is the number of pass
#' branches and M is the number of gates:
#' \itemize{
#'   \item N pass rows -- one per branch satisfying all gates + the
#'     branch -- decision "keep".
#'   \item N drop rows -- one per branch satisfying gates but
#'     failing the branch via an atom break -- decision "remove".
#'   \item M gate-violation rows -- one per gate -- decision "remove".
#'   \item 3 QC-failure rows (Sample_status, cfDNA_ng, Plasma_ml_input)
#'     -- decision "remove".
#' }
#'
#' \strong{How the labeler works (mode 2):} for each row of \code{df}
#' the function checks, in order: sample-level QC (per \code{ctdna_opts}
#' configuration), required columns present, scope match (out-of-scope
#' rows are routed to \code{scheme_basic} just like in
#' \code{ctdna_variant_filter}), each gate in turn, then each pass branch.
#' The first hit determines the decision and reason string.
#'
#' @param scheme_name Character. Name of a registered filtering scheme
#'   (e.g. \code{"TSG_Tier2"}, \code{"HRR14"}, \code{"scheme_basic"})
#'   or a user-registered scheme. Use
#'   \code{\link{ctdna_filter_schemes_list}} to see what's available.
#' @param df Optional data frame. If \code{NULL} (default), the function
#'   simulates a synthetic test frame (mode 1). If a data frame is
#'   provided, the function labels each row (mode 2) and returns the
#'   input frame with two added columns.
#' @param include_qc_failures Mode-1 only. If \code{TRUE} (default),
#'   append rows that specifically fail sample-level QC checks so users
#'   can see QC drops. Ignored in mode 2.
#' @param verify Mode-1 only. If \code{TRUE} (default), run the scheme
#'   on the generated rows and error if the engine disagrees with the
#'   expected decisions. Ignored in mode 2.
#' @param sample_qc_enable Optional logical. Override the global
#'   sample-QC stage just for this call: \code{TRUE} forces QC on,
#'   \code{FALSE} forces it off, \code{NULL} (default) defers to the
#'   value of \code{ctdna_opts("sample_qc_enable")}.
#' @return A data frame with the input rows (or generated rows in mode 1)
#'   plus columns \code{filter_decision}, \code{filter_reason}.
#' @seealso \code{\link{ctdna_variant_filter}},
#'   \code{\link{ctdna_filter_schemes_list}},
#'   \code{\link{create_filtering_scheme}}.
#' @examples
#' \dontrun{
#' # Mode 1: simulate
#' td <- ctdna_test_scheme("TSG_Tier2")
#' table(td$filter_decision)
#'
#' # Mode 2: label the user's data
#' sim <- ctdna_make_mock_study(n_patients = 20, seed = 1)
#' labeled <- ctdna_test_scheme("TSG_Tier2", df = sim$infinity_report)
#' table(labeled$filter_decision)
#' head(labeled[, c("Gene","Variant_type","filter_decision","filter_reason")])
#' }
#' @export
ctdna_test_scheme <- function(scheme_name,
                                df = NULL,
                                include_qc_failures = TRUE,
                                sample_qc_enable = NULL,
                                verify = TRUE) {
  if (!is.character(scheme_name) || length(scheme_name) != 1L)
    stop("`scheme_name` must be a single character string.", call. = FALSE)

  # Resolve QC override (per-call > session opt)
  qc_active <- if (!is.null(sample_qc_enable)) isTRUE(sample_qc_enable)
                else isTRUE(.o("sample_qc_filter"))

  # Resolve the scheme via the existing lookup (handles built-ins +
  # user-registered v1.0.0 schemes).
  scheme <- tryCatch(ctdna_filter_scheme_get(scheme_name),
                       error = function(e) NULL)
  if (is.null(scheme) || !inherits(scheme, "ctdna_filtering_scheme"))
    stop(sprintf("Scheme '%s' not found or not a v1.0.0 rule scheme. ",
                  scheme_name),
          "See ctdna_filter_schemes_list() for available scheme names.",
          call. = FALSE)

  # ----------------------------------------------------------------------
  # Mode 2: user data frame supplied -> label rows, return input + 2 cols
  # ----------------------------------------------------------------------
  if (!is.null(df)) {
    if (!is.data.frame(df))
      stop("`df` must be a data frame (or NULL for simulate mode).",
            call. = FALSE)
    if (nrow(df) == 0L) {
      df$filter_decision <- character(0)
      df$filter_reason   <- character(0)
      return(df)
    }
    # Temporarily set the opt so the explainer mirrors filter_apply
    saved_opt <- .o("sample_qc_filter")
    if (!is.null(sample_qc_enable)) {
      .env$cfg$sample_qc_filter <- qc_active
    }
    on.exit(.env$cfg$sample_qc_filter <- saved_opt, add = TRUE)

    decisions <- character(nrow(df))
    reasons   <- character(nrow(df))
    for (i in seq_len(nrow(df))) {
      row_i <- df[i, , drop = FALSE]
      exp_i <- .ctdna_test_explain_row(row_i, scheme, scheme_name)
      decisions[i] <- exp_i$decision
      reasons[i]   <- exp_i$reason
    }
    df$filter_decision <- decisions
    df$filter_reason   <- reasons
    return(df)
  }

  # ----------------------------------------------------------------------
  # Mode 1: simulate (original behavior)
  # ----------------------------------------------------------------------

  cols <- .ctdna_collect_cols(scheme)
  # Ensure QC columns are present too
  qc_cols <- c("Sample_status","cfDNA_ng","Plasma_ml_input")
  cols <- unique(c(cols, qc_cols))

  # QC defaults for the value bank (all passing)
  qc_defaults <- list(Sample_status   = "SUCCESS",
                       cfDNA_ng        = 10,
                       Plasma_ml_input = 10)
  fill_qc <- function(row) {
    for (k in names(qc_defaults))
      if (k %in% names(row)) row[[k]] <- qc_defaults[[k]]
    row
  }

  # Find pass branches and gates from the scheme root
  parts <- .ctdna_test_find_pass_branches(scheme$root)
  branches <- parts$branches
  gates    <- parts$gates

  # Also: the scheme$scope acts as an implicit gate (if present)
  if (!is.null(scheme$scope))
    gates <- c(list(scheme$scope), gates)

  rows <- list()
  row_id <- 0L

  # ---- Pass rows: one per branch, satisfying all gates + the branch ----
  # Columns touched by gates are PROTECTED when satisfying branches (so a
  # branch's internal NOT doesn't accidentally undo a gate match).
  gate_cols <- unique(unlist(lapply(gates, .ctdna_collect_cols)))
  for (i in seq_along(branches)) {
    branch <- branches[[i]]
    row <- .ctdna_test_blank_row(cols)
    # Satisfy gates first
    for (g in gates) row <- .ctdna_test_satisfy(row, g)
    # Then satisfy the branch, protecting gate columns
    row <- .ctdna_test_satisfy(row, branch, protect = gate_cols)
    row <- fill_qc(row)
    row_id <- row_id + 1L
    row$Row_ID <- row_id
    row$filter_decision <- "keep"
    row$filter_reason <- sprintf("pass branch %d: %s", i,
                                   .ctdna_test_label(branch))
    rows[[length(rows) + 1L]] <- row
  }

  # ---- Drop rows: ONE PER BRANCH ----
  # Strategy: take the row that would PASS branch i, then iteratively
  # break atoms until the FULL scheme evaluates FALSE.
  #
  # For each atom in branch i, we try breaking it; if some other branch
  # is then satisfied (because the value we picked happens to satisfy
  # another path through the OR), we ALSO break an atom in that other
  # branch. We accumulate atom-breaks until no branch passes.
  #
  # This produces realistic near-miss rows: the row is "trying to be"
  # branch i but explicitly fails at least one of its sub-rules AND
  # doesn't accidentally pass any other branch.
  if (length(branches) > 1L) {
    for (i in seq_along(branches)) {
      branch <- branches[[i]]
      # Start from a fresh row that would pass branch i
      base_row <- .ctdna_test_blank_row(cols)
      for (g in gates) base_row <- .ctdna_test_satisfy(base_row, g)
      base_row <- .ctdna_test_satisfy(base_row, branch, protect = gate_cols)

      branch_atoms <- .ctdna_collect_atoms(branch)
      drop_row <- NULL
      broken_atom_labels <- character(0)

      # Try each atom in branch i as the "broken" one
      for (a in branch_atoms) {
        candidate <- .ctdna_test_violate_atom(base_row, a, protect = gate_cols)
        broken_labels <- .ctdna_test_label(a, max_len = 30)
        # If the candidate now passes ANOTHER branch, break an atom in
        # that branch too. Limit to a few iterations to avoid loops.
        for (fix_iter in 1:5) {
          full_pass <- isTRUE(.eval_rule(scheme, candidate))
          if (!full_pass) break
          # Find which other branch passes
          offending <- NULL
          for (j in seq_along(branches)) {
            if (j == i) next
            if (isTRUE(.eval_rule(branches[[j]], candidate))) {
              offending <- branches[[j]]
              break
            }
          }
          if (is.null(offending)) break  # gates failing? scope?
          # Break an atom in the offending branch that makes the FULL
          # scheme FALSE (not just the offending branch). This avoids
          # ping-ponging through shared columns.
          offending_atoms <- .ctdna_collect_atoms(offending)
          broken_now <- FALSE
          for (oa in offending_atoms) {
            candidate2 <- .ctdna_test_violate_atom(candidate, oa,
                                                     protect = gate_cols)
            if (!isTRUE(.eval_rule(scheme, candidate2))) {
              candidate <- candidate2
              broken_labels <- c(broken_labels,
                                  sprintf("also: %s",
                                          .ctdna_test_label(oa, max_len = 25)))
              broken_now <- TRUE
              break
            }
          }
          if (!broken_now) break
        }
        if (!isTRUE(.eval_rule(scheme, candidate))) {
          drop_row <- candidate
          broken_atom_labels <- broken_labels
          break
        }
      }

      # If still no luck, fall back to violating all branches
      if (is.null(drop_row)) {
        candidate <- base_row
        candidate <- .ctdna_test_violate(candidate, branch, protect = gate_cols)
        for (j in seq_along(branches)) {
          if (j != i)
            candidate <- .ctdna_test_violate(candidate, branches[[j]],
                                              protect = gate_cols)
        }
        drop_row <- candidate
        broken_atom_labels <- "all branches violated"
      }

      drop_row <- fill_qc(drop_row)
      row_id <- row_id + 1L
      drop_row$Row_ID <- row_id
      drop_row$filter_decision <- "remove"
      drop_row$filter_reason <- sprintf("near-miss for branch %d (%s)",
        i, paste(broken_atom_labels, collapse = "; "))
      rows[[length(rows) + 1L]] <- drop_row
    }
  }

  # ---- Drop rows: one per gate, violating that gate ----
  for (i in seq_along(gates)) {
    g <- gates[[i]]
    row <- .ctdna_test_blank_row(cols)
    # Satisfy other gates
    for (j in seq_along(gates)) {
      if (j != i) row <- .ctdna_test_satisfy(row, gates[[j]])
    }
    # Satisfy a pass branch
    if (length(branches) > 0L)
      row <- .ctdna_test_satisfy(row, branches[[1]])
    # Violate the chosen gate
    row <- .ctdna_test_violate(row, g)
    row <- fill_qc(row)
    row_id <- row_id + 1L
    row$Row_ID <- row_id
    row$filter_decision <- "remove"
    row$filter_reason <- sprintf("violates gate: %s",
                                   .ctdna_test_label(g))
    rows[[length(rows) + 1L]] <- row
  }

  # ---- QC-failing rows ----
  # Compute values below the CURRENT thresholds (from ctdna_opts) so the
  # QC failure rows actually fail. With sample_qc_plasma_min = 0 (v0.38.0
  # default), the plasma failure row uses -1 so it's strictly below.
  if (isTRUE(include_qc_failures) && qc_active) {
    pass_val   <- .o("sample_qc_status_pass") %||% "SUCCESS"
    cfdna_min  <- .o("sample_qc_cfdna_min")   %||% 5
    plasma_min <- .o("sample_qc_plasma_min")  %||% 0
    # Pick failing values: anything different from pass_val for status,
    # and threshold - 1 for numeric (works for 0 too since -1 < 0).
    qc_cases <- list(
      list(col = "Sample_status",   value = paste0(pass_val, "_FAIL"),
            reason = sprintf("QC: Sample_status != %s", pass_val)),
      list(col = "cfDNA_ng",        value = cfdna_min - 1,
            reason = "QC: cfDNA_ng below threshold"),
      list(col = "Plasma_ml_input", value = plasma_min - 1,
            reason = "QC: Plasma_ml_input below threshold"))
    for (qc in qc_cases) {
      row <- .ctdna_test_blank_row(cols)
      # Satisfy gates + first pass branch (so without QC failure this would pass)
      for (g in gates) row <- .ctdna_test_satisfy(row, g)
      if (length(branches) > 0L)
        row <- .ctdna_test_satisfy(row, branches[[1]], protect = gate_cols)
      row <- fill_qc(row)
      row[[qc$col]] <- qc$value
      row_id <- row_id + 1L
      row$Row_ID <- row_id
      row$filter_decision <- "remove"
      row$filter_reason <- qc$reason
      rows[[length(rows) + 1L]] <- row
    }
  }

  # Bind rows. They may have different column subsets if blank-row builder
  # missed something; rbind.fill via Reduce is safest.
  all_cols <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(r) {
    missing_cols <- setdiff(all_cols, names(r))
    for (mc in missing_cols) r[[mc]] <- NA
    r[all_cols]
  })
  df <- do.call(rbind, rows)
  rownames(df) <- NULL

  # Reorder: Row_ID first, then scheme columns, then QC, then decision/reason
  scheme_cols <- setdiff(all_cols, c("Row_ID", qc_cols,
                                       "filter_decision","filter_reason"))
  ordered_cols <- c("Row_ID", scheme_cols, qc_cols,
                     "filter_decision","filter_reason")
  ordered_cols <- intersect(ordered_cols, names(df))
  df <- df[, ordered_cols, drop = FALSE]

  # ---- Self-verify: run the scheme + QC, compare engine vs expected ----
  if (isTRUE(verify)) {
    expected_keep <- df$Row_ID[df$filter_decision == "keep"]
    out <- tryCatch(
      ctdna_variant_filter(df, filter_scheme = scheme_name, verbose = FALSE),
      error = function(e) {
        stop("ctdna_test_scheme: filter engine errored on generated data: ",
              conditionMessage(e), call. = FALSE)
      })
    actual_keep <- out$Row_ID
    missing_keeps <- setdiff(expected_keep, actual_keep)
    extra_keeps   <- setdiff(actual_keep, expected_keep)
    if (length(missing_keeps) > 0 || length(extra_keeps) > 0) {
      diag <- data.frame(
        Row_ID = sort(unique(c(missing_keeps, extra_keeps))),
        stringsAsFactors = FALSE)
      diag$expected <- ifelse(diag$Row_ID %in% expected_keep, "keep", "remove")
      diag$actual   <- ifelse(diag$Row_ID %in% actual_keep,   "keep", "remove")
      diag$reason   <- df$filter_reason[match(diag$Row_ID, df$Row_ID)]
      msg <- paste0(
        "ctdna_test_scheme: engine disagrees with expected decisions on ",
        nrow(diag), " row(s):\n",
        paste(capture.output(print(diag, row.names = FALSE)), collapse = "\n"),
        "\n\nThis indicates either a generator bug or a real scheme issue.")
      stop(msg, call. = FALSE)
    }
  }

  df
}
