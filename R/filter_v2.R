# =============================================================================
# ctdnaTM v0.71.0 -- filtering makeover
#
# Four public functions total:
#   ctdna_sample_qc      -- sample QC, cascades to variants (kept, unchanged)
#   ctdna_apply_filter   -- variant filtering; prep -> prep; drops rows
#   ctdna_explain_filter -- dry-run; adds filtering_criteria column, no drops
#   ctdna_create_scheme  -- thin builder for user-facing scheme construction
#
# Design rules:
#   * Sample QC is NOT run inside ctdna_prepare -- the user calls
#     ctdna_sample_qc(prep) explicitly. ctdna_apply_filter and
#     ctdna_explain_filter both auto-run ctdna_sample_qc first
#     (idempotent -- reruns drop nothing new).
#   * NO sentinel value trick. filter_scheme = NULL means "no filter"
#     (nothing evaluated, nothing dropped). filter_scheme must be given
#     to actually filter anything.
#   * NO auto-scheme_basic sweep. If unclaimed rows exist, they follow
#     the `unmatched` policy: "drop" (default) or "keep".
#   * NO category system, no deprecated aliases ("general", "default").
# =============================================================================

# ---------------------------------------------------------------------------
# Internal engine: evaluate schemes on a variant data.frame.
# Returns either the FILTERED df (record_criteria = FALSE) or the FULL df
# with a `filtering_criteria` column (record_criteria = TRUE, nothing
# dropped). Multi-scheme semantics: for each row, evaluate every in-scope
# scheme; row PASSES if any in-scope scheme passes. Rows in no scheme's
# scope are handled by `unmatched`.
# ---------------------------------------------------------------------------

.filter_engine <- function(df,
                            filter_scheme,
                            unmatched       = c("drop", "keep"),
                            record_criteria = FALSE,
                            fn_label        = "ctdna_apply_filter",
                            verbose         = FALSE) {

  unmatched <- match.arg(unmatched)
  msg <- function(...) if (isTRUE(verbose)) message(...)

  # No-filter shortcut
  if (is.null(filter_scheme)) {
    msg(fn_label, ": filter_scheme = NULL -> no filtering.")
    if (isTRUE(record_criteria)) {
      df$filtering_criteria <- NA_character_
      return(df)
    }
    return(df)
  }

  schemes <- .ctdna_resolve_schemes(filter_scheme)
  msg(sprintf("%s: %d scheme(s): %s",
              fn_label, length(schemes),
              paste(names(schemes), collapse = ", ")))

  n0 <- nrow(df)
  if (n0 == 0L) {
    if (isTRUE(record_criteria)) df$filtering_criteria <- character(0)
    return(df)
  }

  # Per-scheme scope masks + rule pass masks
  scope_masks <- vector("list", length(schemes))
  pass_masks  <- vector("list", length(schemes))
  names(scope_masks) <- names(schemes)
  names(pass_masks)  <- names(schemes)
  for (nm in names(schemes)) {
    sch <- schemes[[nm]]
    if (is.null(sch$scope)) {
      scope_masks[[nm]] <- rep(TRUE, n0)
    } else {
      sm <- .eval_rule(sch$scope, df)
      sm[is.na(sm)] <- FALSE
      scope_masks[[nm]] <- sm
    }
    fm <- .eval_rule(sch, df)
    fm[is.na(fm)] <- FALSE
    pass_masks[[nm]] <- scope_masks[[nm]] & fm
  }

  # Per-row: which schemes had it in scope; which passed
  any_in_scope <- Reduce(`|`, scope_masks)
  any_passed   <- Reduce(`|`, pass_masks)

  # Row-final: passes if ANY in-scope scheme passes.
  # Rows in no scope: follow `unmatched` policy.
  keep <- any_passed | (!any_in_scope & unmatched == "keep")

  if (isTRUE(record_criteria)) {
    # Pre-compute per-scheme, per-branch pass masks (only for schemes that
    # declare `branches`). scope+gate pass is implicit in pass_masks[[nm]].
    branch_masks <- lapply(names(schemes), function(nm) {
      sch <- schemes[[nm]]
      if (is.null(sch$branches) || !length(sch$branches)) return(NULL)
      # For each named branch, mask = branch-fires AND scope-passes AND gates-pass.
      # `pass_masks[[nm]]` already encodes scope + gates + body. But body includes
      # the anyOf across branches -- to know WHICH branch fired we evaluate them
      # individually, still ANDed with the gate/exclusion rules from sch$gates.
      gate_ok <- if (length(sch$gates)) {
        m <- rep(TRUE, n0)
        for (g in sch$gates) {
          gm <- .eval_rule(g, df); gm[is.na(gm)] <- FALSE
          m <- m & gm
        }
        m
      } else rep(TRUE, n0)
      lapply(sch$branches, function(br) {
        bm <- .eval_rule(br, df); bm[is.na(bm)] <- FALSE
        bm & gate_ok
      })
    })
    names(branch_masks) <- names(schemes)

    crit <- rep(NA_character_, n0)
    for (i in seq_len(n0)) {
      if (!any_in_scope[i]) next   # scope-miss: no annotation
      pass_here <- vapply(names(schemes), function(nm)
        pass_masks[[nm]][i], logical(1))
      fail_here <- vapply(names(schemes), function(nm)
        scope_masks[[nm]][i] && !pass_masks[[nm]][i], logical(1))
      if (any(pass_here)) {
        parts <- vapply(names(schemes)[pass_here], function(nm) {
          bm <- branch_masks[[nm]]
          if (is.null(bm)) return(nm)                    # no branch info
          fired <- names(bm)[vapply(bm, function(x) x[i], logical(1))]
          if (!length(fired)) return(nm)                 # shouldn't happen
          sprintf("%s(%s)", nm, paste(fired, collapse = "+"))
        }, character(1))
        crit[i] <- paste0("PASS::", paste(parts, collapse = ", "))
      } else {
        # FAIL: keep scheme name only (branch info doesn't apply to failures)
        crit[i] <- paste0("FAIL::",
                          paste(names(schemes)[fail_here], collapse = ", "))
      }
    }
    df$filtering_criteria <- crit
    n_pass <- sum(!is.na(crit) & startsWith(crit, "PASS::"))
    n_fail <- sum(!is.na(crit) & startsWith(crit, "FAIL::"))
    n_na   <- sum(is.na(crit))
    msg(sprintf("%s: PASS=%d, FAIL=%d, scope-miss(NA)=%d.",
                fn_label, n_pass, n_fail, n_na))
    return(df)
  }

  out <- df[keep, , drop = FALSE]
  rownames(out) <- NULL
  msg(sprintf("%s: %d -> %d rows.", fn_label, n0, nrow(out)))
  out
}


# ---------------------------------------------------------------------------
# ctdna_variant_filter -- public: variant filter with two mandatory flags
# ---------------------------------------------------------------------------

#' Variant-level filtering with mandatory apply / explain flags
#'
#' Evaluates one or more filter schemes against \code{prep$variants} and,
#' according to two required logical flags, either drops failing rows
#' (\code{apply = TRUE}), stores an annotated per-variant snapshot
#' explaining every scheme decision (\code{explain = TRUE}), or both.
#' This function operates strictly on variants; it does not touch
#' samples or run sample-level QC -- call \code{\link{ctdna_sample_qc}}
#' first if you want that.
#'
#' @details
#' \strong{Two mandatory flags, four modes.} Both \code{apply} and
#' \code{explain} must be supplied explicitly (no defaults). The four
#' combinations behave as follows:
#'
#' \tabular{lll}{
#'   \code{apply} \tab \code{explain} \tab effect \cr
#'   \code{TRUE}  \tab \code{TRUE}  \tab \code{prep$variants} filtered
#'     \emph{and} annotated snapshot stored at
#'     \code{prep$filter_explanation$<sig>} \cr
#'   \code{TRUE}  \tab \code{FALSE} \tab \code{prep$variants} filtered;
#'     no snapshot \cr
#'   \code{FALSE} \tab \code{TRUE}  \tab \code{prep$variants} unchanged;
#'     annotated snapshot stored (audit / dry-run mode) \cr
#'   \code{FALSE} \tab \code{FALSE} \tab error -- nothing to do
#' }
#'
#' \strong{Filtering-criteria column values} (present in the annotated
#' snapshot at \code{prep$filter_explanation$<sig>$filtering_criteria},
#' never on the live \code{prep$variants}):
#'
#' \itemize{
#'   \item \code{NA} -- row was in no scheme's scope. With
#'     \code{apply = TRUE, unmatched = "drop"} (default) these rows are
#'     dropped; with \code{unmatched = "keep"} they survive.
#'   \item \code{"PASS::<scheme>"} -- row passed at least one in-scope
#'     scheme.
#'   \item \code{"PASS::<scheme>(<branch>)"} -- for schemes built with
#'     named branches (e.g. \code{create_filtering_scheme(...,
#'     branches = list(somatic_nonsyn = ..., germline_PLP = ...))}),
#'     the branch that fired is reported alongside the scheme name.
#'   \item \code{"FAIL::<scheme>"} -- in the scheme's scope but failed
#'     every branch.
#' }
#'
#' When multiple schemes are passed (e.g. \code{filter_scheme =
#' c("TSG_Tier1","HRR14")}), PASS and FAIL labels concatenate with
#' commas: \code{"PASS::TSG_Tier1(somatic_nonsyn), FAIL::HRR14"}.
#'
#' \strong{Signature name.} \code{<sig>} is a filesystem-safe
#' concatenation of the passed scheme names, e.g.
#' \code{"TSG_Tier1_HRR14"}. This is what appears as the slot name in
#' \code{prep$filter_explanation}.
#'
#' @param prep A \code{ctdna_prep} object, typically after
#'   \code{\link{ctdna_sample_qc}}.
#' @param filter_scheme The scheme (or schemes) to evaluate against
#'   \code{prep$variants}. Accepts any of:
#'   \itemize{
#'     \item A character scalar -- a built-in scheme name
#'       (\code{"TSG_Tier1"}, \code{"TSG_Tier2"}, \code{"TSG_Tier3"},
#'       \code{"HRR14"}, or \code{"scheme_basic"}) or the name of a
#'       user-registered scheme.
#'     \item A character vector -- evaluate each named scheme
#'       independently and concatenate their PASS/FAIL labels.
#'     \item A \code{ctdna_filtering_scheme} object as returned by
#'       \code{\link{create_filtering_scheme}} or
#'       \code{\link{ctdna_create_scheme}}.
#'     \item A list of any of the above.
#'   }
#'   Required (no default). See \code{\link{ctdna_filter_schemes_list}}
#'   for the registered catalog.
#' @param apply Logical. When \code{TRUE}, remove \code{FAIL::*} rows
#'   from \code{prep$variants}; also remove scope-miss (\code{NA}) rows
#'   unless \code{unmatched = "keep"}. When \code{FALSE},
#'   \code{prep$variants} is left untouched (use with
#'   \code{explain = TRUE} to preview what a filter \emph{would} do
#'   without modifying the prep). No default -- pass explicitly.
#' @param explain Logical. When \code{TRUE}, store the per-variant
#'   annotated snapshot (with the \code{filtering_criteria} column) at
#'   \code{prep$filter_explanation$<sig>}. When \code{FALSE}, no
#'   snapshot is retained. No default -- pass explicitly.
#' @param unmatched How to handle rows that no scheme's scope covers
#'   (i.e. \code{filtering_criteria == NA}) when \code{apply = TRUE}.
#'   \code{"drop"} (default) removes them; \code{"keep"} treats
#'   scope-miss as pass-through. Ignored when \code{apply = FALSE}.
#' @param verbose Logical. When \code{TRUE} (default), prints the
#'   scheme(s) evaluated, PASS/FAIL/scope-miss counts, the snapshot
#'   storage location (if \code{explain = TRUE}), and the pre/post row
#'   counts on the live frame (if \code{apply = TRUE}).
#'
#' @return The input \code{ctdna_prep}, modified in place according to
#'   the flag combination:
#'   \itemize{
#'     \item \code{$variants} -- filtered when \code{apply = TRUE},
#'       otherwise untouched. The \code{filtering_criteria} column is
#'       \emph{never} attached to the live frame.
#'     \item \code{$filter_explanation$<sig>} -- annotated snapshot
#'       (all input variants + a \code{filtering_criteria} column)
#'       when \code{explain = TRUE}. Multiple explain calls with
#'       different schemes accumulate in this list.
#'   }
#'
#' @seealso \code{\link{ctdna_sample_qc}} (mandatory previous step);
#'   \code{\link{ctdna_filter_schemes_list}} (view the scheme catalog);
#'   \code{\link{ctdna_filter_scheme_show}} (inspect a specific
#'   scheme's rule tree); \code{\link{create_filtering_scheme}} and
#'   \code{\link{ctdna_create_scheme}} (build your own scheme);
#'   \code{\link{TSG_Tier1}}, \code{\link{TSG_Tier2}},
#'   \code{\link{TSG_Tier3}}, \code{\link{HRR14}} (built-in schemes).
#'
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 20, seed = 1)
#' prep <- ctdna_prepare(
#'   infinity_report = sim$infinity_report,
#'   adam            = list(adsl = sim$clinical),
#'   verbose         = FALSE)
#' prep <- ctdna_sample_qc(prep, verbose = FALSE)
#'
#' # Basic: apply a single built-in scheme, keep the audit snapshot too.
#' prep1 <- ctdna_variant_filter(prep,
#'   filter_scheme = "TSG_Tier1",
#'   apply         = TRUE,
#'   explain       = TRUE,
#'   verbose       = FALSE)
#' print(prep1)
#'
#' # Dry-run mode: see what the filter WOULD do without touching
#' # $variants. The annotated snapshot is stored, the live frame is
#' # unchanged.
#' prep2 <- ctdna_variant_filter(prep,
#'   filter_scheme = "TSG_Tier1",
#'   apply         = FALSE,
#'   explain       = TRUE,
#'   verbose       = FALSE)
#' # Inspect the filtering_criteria distribution on the snapshot:
#' snap <- prep2$filter_explanation$TSG_Tier1
#' print(table(snap$filtering_criteria, useNA = "ifany"))
#'
#' # Combining schemes: PASS if either TSG_Tier1 OR HRR14 accepts.
#' prep3 <- ctdna_variant_filter(prep,
#'   filter_scheme = c("TSG_Tier1", "HRR14"),
#'   apply         = TRUE,
#'   explain       = TRUE,
#'   verbose       = FALSE)
#' # The snapshot is stored under the combined signature:
#' print(names(prep3$filter_explanation))
#' @export
ctdna_variant_filter <- function(prep,
                                 filter_scheme,
                                 apply,
                                 explain,
                                 unmatched = c("drop", "keep"),
                                 verbose   = TRUE) {

  if (!inherits(prep, "ctdna_prep"))
    stop("ctdna_variant_filter: `prep` must be a ctdna_prep ",
         "(output of ctdna_prepare()).", call. = FALSE)
  if (missing(filter_scheme))
    stop("ctdna_variant_filter: `filter_scheme` is required.", call. = FALSE)
  if (missing(apply) || !is.logical(apply) || length(apply) != 1L)
    stop("ctdna_variant_filter: `apply` must be TRUE or FALSE ",
         "(no default; pass explicitly).", call. = FALSE)
  if (missing(explain) || !is.logical(explain) || length(explain) != 1L)
    stop("ctdna_variant_filter: `explain` must be TRUE or FALSE ",
         "(no default; pass explicitly).", call. = FALSE)
  if (!apply && !explain)
    stop("ctdna_variant_filter: `apply = FALSE` and `explain = FALSE` ",
         "leaves nothing to do. Set at least one to TRUE.",
         call. = FALSE)
  unmatched <- match.arg(unmatched)

  # Always evaluate with record_criteria=TRUE so we have a single
  # source of truth. apply/explain then decide what to store.
  annotated <- .filter_engine(prep$variants,
                               filter_scheme   = filter_scheme,
                               unmatched       = unmatched,
                               record_criteria = TRUE,
                               fn_label        = "ctdna_variant_filter",
                               verbose         = verbose)

  if (isTRUE(explain)) {
    sig <- .filter_scheme_sig(filter_scheme)
    if (is.null(prep$filter_explanation))
      prep$filter_explanation <- list()
    prep$filter_explanation[[sig]] <- annotated
    if (isTRUE(verbose))
      message(sprintf(
        "ctdna_variant_filter: annotated snapshot stored at prep$filter_explanation$%s.",
        sig))
  }

  if (isTRUE(apply)) {
    # Drop rows: FAIL::* always dropped; NA (scope-miss) obeys `unmatched`.
    crit <- annotated$filtering_criteria
    fail_mask <- !is.na(crit) & startsWith(crit, "FAIL::")
    miss_mask <- is.na(crit)
    keep <- !fail_mask & !(unmatched == "drop" & miss_mask)
    filtered <- annotated[keep, , drop = FALSE]
    # Strip the criteria column from the LIVE variants frame -- it belongs
    # only on the explanation snapshot.
    filtered$filtering_criteria <- NULL
    rownames(filtered) <- NULL
    prep$variants <- filtered
    if (isTRUE(verbose))
      message(sprintf(
        "ctdna_variant_filter: applied filter, %d -> %d variants.",
        nrow(annotated), nrow(filtered)))
  }

  prep
}

# Signature name for storing an explanation frame. Concatenates scheme names
# with `_`, strips non-safe chars.
.filter_scheme_sig <- function(filter_scheme) {
  if (is.null(filter_scheme)) return("no_filter")
  if (inherits(filter_scheme, "ctdna_filtering_scheme")) {
    nm <- filter_scheme$name %||% "user_scheme"
    return(gsub("[^A-Za-z0-9_]+", "_", nm))
  }
  if (is.list(filter_scheme))
    nm <- vapply(filter_scheme, function(x)
      if (is.character(x)) x[1] else (x$name %||% "user_scheme"), character(1))
  else
    nm <- as.character(filter_scheme)
  paste(gsub("[^A-Za-z0-9_]+", "_", nm), collapse = "_")
}


# ---------------------------------------------------------------------------
# ctdna_create_scheme -- thin builder for user-facing scheme construction
# ---------------------------------------------------------------------------

#' Create a variant filter scheme (user-friendly wrapper)
#'
#' Thin wrapper over \code{\link{create_filtering_scheme}} that hides the
#' internal rule DSL (\code{rule}/\code{allOf}/\code{anyOf}/\code{not})
#' behind common keyword arguments. Advanced users can still call
#' \code{create_filtering_scheme()} directly.
#'
#' @param name Character. Scheme name (also stored in the catalog).
#' @param genes Character vector of gene symbols to scope this scheme to.
#'   \code{NULL} = scheme applies to all genes.
#' @param criteria Character vector of built-in criteria to apply. One or
#'   more of: \code{"truncating"}, \code{"missense"}, \code{"clinvar_path"},
#'   \code{"cnv_del"} (homozygous or LOH deletion), \code{"cnv_amp"},
#'   \code{"lgr"}, \code{"fusion"}, \code{"somatic"}, \code{"germline"},
#'   \code{"rare"} (<1\% gnomAD AF).
#' @param combine \code{"any"} (default; row passes if ANY criterion matches)
#'   or \code{"all"} (row passes only if ALL match).
#' @param description Optional description.
#' @param overwrite Logical. If \code{TRUE}, allow replacing an existing
#'   scheme with the same \code{name}. Default \code{FALSE} (errors on
#'   collision, protecting against accidental clobber).
#' @return A \code{ctdna_filtering_scheme} registered in the catalog.
#' @seealso \code{\link{create_filtering_scheme}} (the underlying rule-DSL
#'   constructor); \code{\link{HRR14}}, \code{\link{TSG_Tier1}},
#'   \code{\link{scheme_basic}} (built-in schemes);
#'   \code{\link{ctdna_variant_filter}} (apply schemes to a prep).
#' @examples
#' # Scheme: TP53 / RB1 / PTEN, somatic + (truncating OR ClinVar path)
#' sch <- ctdna_create_scheme(
#'   name     = "my_short_tsg",
#'   genes    = c("TP53","RB1","PTEN"),
#'   criteria = c("somatic","truncating","clinvar_path"))
#' print(sch)
#'
#' # Same criteria but require ALL to match (stricter).
#' sch_strict <- ctdna_create_scheme(
#'   name      = "my_short_tsg_strict",
#'   genes     = c("TP53","RB1","PTEN"),
#'   criteria  = c("somatic","truncating","clinvar_path"),
#'   combine   = "all")
#' print(sch_strict)
#'
#' # Use the scheme in the standard filter workflow.
#' sim  <- ctdna_make_mock_study(n_patients = 30, seed = 1)
#' prep <- ctdna_prepare(infinity_report = sim$infinity_report,
#'                       adam = list(adsl = sim$clinical), verbose = FALSE)
#' prep <- ctdna_sample_qc(prep, verbose = FALSE)
#' prep <- ctdna_variant_filter(prep, "my_short_tsg",
#'                              apply = TRUE, explain = TRUE, verbose = FALSE)
#' print(head(prep$filter_explanation$my_short_tsg$filtering_criteria))
#' @export
ctdna_create_scheme <- function(name,
                                genes       = NULL,
                                criteria    = c("truncating","missense",
                                                 "clinvar_path"),
                                combine     = c("any","all"),
                                description = NULL,
                                overwrite   = FALSE) {

  if (missing(name) || !is.character(name) || length(name) != 1L)
    stop("ctdna_create_scheme: `name` (single string) is required.",
         call. = FALSE)
  combine <- match.arg(combine)

  # Map criteria keywords -> rule functions
  criterion_map <- list(
    truncating   = quote(rule_truncating()),
    missense     = quote(rule(Molecular_consequence = "missense")),
    clinvar_path = quote(rule_clinvar_path()),
    cnv_del      = quote(anyOf(rule_cnv_homozyg_del(), rule_cnv_loh_del())),
    cnv_amp      = quote(rule(CNV_type = c("amplification",
                                            "focal_amplification",
                                            "aneuploid_amplification"))),
    lgr          = quote(rule_is_LGR()),
    fusion       = quote(rule_is_Fusion()),
    somatic      = quote(rule_somatic()),
    germline     = quote(rule_germline()),
    rare         = quote(rule_rare_001()))
  unknown <- setdiff(criteria, names(criterion_map))
  if (length(unknown))
    stop("ctdna_create_scheme: unknown criteria: ",
         paste(unknown, collapse = ", "),
         "\n  Available: ",
         paste(names(criterion_map), collapse = ", "), call. = FALSE)

  # Build rules from criteria
  crit_rules <- lapply(criteria, function(k) eval(criterion_map[[k]]))
  body_rule <- if (combine == "any")
    do.call(anyOf, crit_rules)
  else
    do.call(allOf, crit_rules)

  # Scope rule (gene list) -- optional
  scope_rule <- if (!is.null(genes) && length(genes))
    rule(Gene = as.character(genes)) else NULL

  if (is.null(scope_rule))
    create_filtering_scheme(body_rule,
                            name = name, description = description,
                            overwrite = overwrite)
  else
    create_filtering_scheme(scope_rule, body_rule,
                            name = name, description = description,
                            overwrite = overwrite)
}
