# =============================================================================
# v0.30.1 — visualized_scheme(): text-mode decision-tree rendering
# =============================================================================
#
# Renders a filtering scheme as an indented, human-readable text layout
# that separates GATES (AND'd preconditions) from SUCCESS PATHS (OR'd
# alternatives) and prints all values inline (no abbreviation).
#
# Example output for HRR14():
#
#   HRR14  (category: gene_set)
#   HRR / BRCAm classification ...
#
#   REQUIRE (gates -- if any fail -> DROP):
#     |- Gene in {BRCA1, BRCA2, ATM, BARD1, ..., RAD54L}
#     |- ClinVar not in {Benign, Likely_benign, Benign/Likely_benign}
#     `- Somatic_status != somatic_putative_ch     (not CHIP)
#
#   KEEP if ANY of these paths succeeds:
#
#     PATH A -- ClinVar
#       `- ClinVar in {Pathogenic, Likely_pathogenic, ...}
#
#     PATH B -- Molecular_consequence + gnomAD_AF + not Gene + Mut_aa
#       |- Molecular_consequence in {nonsense, frameshift, ...}
#       |- gnomAD_AF < 0.01
#       `- NOT (Gene = BRCA2 AND Mut_aa in {...})
#
#     ...
#
#   Otherwise -> DROP.
#
# Design rationale: graphical visualizers (Graphviz / ggplot) were
# brittle and hard to read at the scheme sizes we use. Text is
# deterministic, copy-pastable into SOPs / logs / Rmd chunks, and the
# AND/OR structure is immediately obvious from indentation.


# -----------------------------------------------------------------------------
# Scheme decomposition: split root into gates + success paths
# -----------------------------------------------------------------------------

# Decompose a scheme into (gates, paths). Returns list(gates, paths).
# - If root is allOf(g1, g2, ..., or_node)   -> gates = [g1, g2, ...],
#                                                paths = or_node$children
# - If root is allOf(g1, g2, ...)            -> gates = [g1, ...], paths = []
# - If root is a single non-and node          -> gates = [], paths = [root]
.decompose_scheme <- function(scheme) {
  root <- scheme$root
  if (inherits(root, "ctdna_rule_combinator") && root$type == "and") {
    kids <- root$children
    or_idx <- which(vapply(kids, function(k)
      inherits(k, "ctdna_rule_combinator") && k$type == "or", logical(1)))
    if (length(or_idx) > 0) {
      last_or <- tail(or_idx, 1)
      return(list(gates = kids[-last_or],
                  paths = kids[[last_or]]$children))
    }
    return(list(gates = kids, paths = list()))
  }
  # Non-AND root: treat as a single path with no gates
  list(gates = list(), paths = list(root))
}


# -----------------------------------------------------------------------------
# Atom rendering: one condition string per column reference
# -----------------------------------------------------------------------------

# Render a single rule_atom's conditions as one or more strings.
# Multiple columns inside one atom (AND'd) produce multiple strings.
# Values are NEVER abbreviated -- all listed inline.
.fmt_atom <- function(atom) {
  conds <- atom$conditions
  vapply(names(conds), function(col) {
    v <- conds[[col]]
    # Operator form
    if (is.list(v) && all(c("op","value") %in% names(v))) {
      val_str <- if (length(v$value) == 1) format(v$value) else
        sprintf("{%s}", paste(format(v$value), collapse = ", "))
      return(sprintf("%s %s %s", col, v$op, val_str))
    }
    # Plain value(s)
    vc <- as.character(v)
    if (length(vc) == 1) return(sprintf("%s = %s", col, vc))
    sprintf("%s in {%s}", col, paste(vc, collapse = ", "))
  }, character(1))
}

# Render a NOT subtree as one or more "not X" condition strings.
# - not(atom)        -> ["NOT (Col op val)"]
# - not(allOf(...)) -> ["NOT (cond1 AND cond2 AND ...)"]
# - not(anyOf(...)) -> ["NOT (cond1 OR cond2 OR ...)"]
.fmt_not <- function(node) {
  child <- node$child
  if (inherits(child, "ctdna_rule_atom")) {
    parts <- .fmt_atom(child)
    if (length(parts) == 1L) return(sprintf("NOT (%s)", parts))
    return(sprintf("NOT (%s)", paste(parts, collapse = " AND ")))
  }
  if (inherits(child, "ctdna_rule_combinator")) {
    parts <- .fmt_rule_inline(child)
    return(sprintf("NOT (%s)", parts))
  }
  "NOT (?)"
}

# Render any rule node as a single inline string (no newlines).
# Used for negated sub-trees and for OR-inside-path collapses.
.fmt_rule_inline <- function(node) {
  if (inherits(node, "ctdna_rule_atom")) {
    parts <- .fmt_atom(node)
    return(paste(parts, collapse = " AND "))
  }
  if (inherits(node, "ctdna_rule_combinator")) {
    if (node$type == "not") return(.fmt_not(node))
    if (node$type == "and") {
      parts <- vapply(node$children, .fmt_rule_inline, character(1))
      return(paste(parts, collapse = " AND "))
    }
    if (node$type == "or") {
      # If all OR children are atoms on the same single column, collapse
      common <- .or_common_column(node)
      if (!is.null(common)) {
        vals <- unlist(lapply(node$children, function(ch) {
          v <- ch$conditions[[common]]
          if (is.list(v) && all(c("op","value") %in% names(v))) return(v$value)
          v
        }))
        vals <- unique(as.character(vals))
        return(sprintf("%s in {%s}", common, paste(vals, collapse = ", ")))
      }
      parts <- vapply(node$children, .fmt_rule_inline, character(1))
      return(sprintf("(%s)", paste(parts, collapse = " OR ")))
    }
  }
  "?"
}

# Detect if every child of an or-node is an atom referencing the same
# single column. Returns the column name or NULL.
.or_common_column <- function(or_node) {
  if (!inherits(or_node, "ctdna_rule_combinator") || or_node$type != "or")
    return(NULL)
  cols <- lapply(or_node$children, function(ch) {
    if (!inherits(ch, "ctdna_rule_atom")) return(NULL)
    nms <- names(ch$conditions)
    if (length(nms) != 1L) return(NULL)
    # Reject operator-form constraints (mixing < and = on the same col
    # collapses awkwardly); only pool plain value matches.
    v <- ch$conditions[[nms]]
    if (is.list(v) && all(c("op","value") %in% names(v))) return(NULL)
    nms
  })
  if (any(vapply(cols, is.null, logical(1)))) return(NULL)
  cols <- unlist(cols)
  if (length(unique(cols)) == 1L) return(cols[1])
  NULL
}


# -----------------------------------------------------------------------------
# Expand a rule tree into a list of condition strings for a path / gate.
# Top-level allOf(...) flattens; everything else collapses inline.
# -----------------------------------------------------------------------------
.expand_conditions <- function(node) {
  if (inherits(node, "ctdna_rule_atom")) return(.fmt_atom(node))
  if (inherits(node, "ctdna_rule_combinator")) {
    if (node$type == "and") {
      out <- character(0)
      for (ch in node$children) out <- c(out, .expand_conditions(ch))
      return(out)
    }
    if (node$type == "not") return(.fmt_not(node))
    if (node$type == "or")  return(.fmt_rule_inline(node))
  }
  "?"
}


# -----------------------------------------------------------------------------
# Path label heuristic: list the distinct columns referenced
# -----------------------------------------------------------------------------
.path_hint <- function(node) {
  if (inherits(node, "ctdna_rule_atom"))
    return(paste(names(node$conditions), collapse = " + "))
  if (inherits(node, "ctdna_rule_combinator")) {
    if (node$type == "not") {
      h <- .path_hint(node$child)
      return(if (!is.null(h)) sprintf("not %s", h) else NULL)
    }
    if (node$type == "and") {
      hints <- vapply(node$children,
                       function(ch) .path_hint(ch) %||% "",
                       character(1))
      hints <- unique(hints[nzchar(hints)])
      if (length(hints) == 0) return(NULL)
      return(paste(hints, collapse = " + "))
    }
    if (node$type == "or") {
      common <- .or_common_column(node)
      if (!is.null(common)) return(common)
      return(.path_hint(node$children[[1]]))
    }
  }
  NULL
}


# -----------------------------------------------------------------------------
# Indent helpers using ASCII box-drawing
# -----------------------------------------------------------------------------
.bullets <- function(items, prefix = "  ") {
  n <- length(items)
  if (n == 0) return(character(0))
  out <- character(n)
  for (i in seq_len(n)) {
    marker <- if (i == n) "`- " else "|- "
    out[i] <- paste0(prefix, marker, items[i])
  }
  out
}


# -----------------------------------------------------------------------------
# Main: text-render a scheme
# -----------------------------------------------------------------------------
.render_scheme_text <- function(scheme, width = 80) {
  d <- .decompose_scheme(scheme)
  out <- character(0)

  # Header
  hdr <- sprintf("%s  (category: %s)",
                  scheme$name %||% "<unnamed>",
                  scheme$category %||% "user")
  out <- c(out, hdr)
  if (!is.null(scheme$description) && nzchar(scheme$description)) {
    desc_wrapped <- strwrap(scheme$description, width = width)
    out <- c(out, desc_wrapped)
  }
  out <- c(out, "")

  # Gates
  if (length(d$gates) > 0) {
    gate_strs <- character(0)
    for (g in d$gates) gate_strs <- c(gate_strs, .expand_conditions(g))
    out <- c(out, "REQUIRE (gates -- if any fail -> DROP):")
    out <- c(out, .bullets(gate_strs, prefix = "  "))
    out <- c(out, "")
  }

  # Paths
  if (length(d$paths) > 0) {
    if (length(d$paths) == 1L) {
      # Single "path" -- just list its conditions as required
      conds <- .expand_conditions(d$paths[[1]])
      out <- c(out, "REQUIRE (additionally):")
      out <- c(out, .bullets(conds, prefix = "  "))
    } else {
      out <- c(out, "KEEP if ANY of these paths succeeds:")
      out <- c(out, "")
      for (i in seq_along(d$paths)) {
        p <- d$paths[[i]]
        letter <- LETTERS[i]
        hint <- .path_hint(p) %||% ""
        if (nzchar(hint))
          out <- c(out, sprintf("  PATH %s -- %s", letter, hint))
        else
          out <- c(out, sprintf("  PATH %s", letter))
        conds <- .expand_conditions(p)
        out <- c(out, .bullets(conds, prefix = "    "))
        out <- c(out, "")
      }
    }
    out <- c(out, "Otherwise -> DROP.")
  } else if (length(d$gates) > 0) {
    out <- c(out, "(No alternative paths -- if all gates pass, KEEP.)")
  }

  out
}


# -----------------------------------------------------------------------------
# Exported entry point
# -----------------------------------------------------------------------------

#' Visualize a filtering scheme as a text decision tree
#'
#' Renders a \code{ctdna_filtering_scheme} (or a scheme name resolvable
#' via the catalog) as an indented text layout that separates AND'd
#' gates from OR'd success paths.
#'
#' Output layout:
#' \itemize{
#'   \item \strong{Header} -- scheme name + category + description.
#'   \item \strong{REQUIRE (gates)} -- the AND'd preconditions. If any
#'     fail, the row is DROPped.
#'   \item \strong{KEEP if ANY of these paths succeeds} -- the OR'd
#'     alternatives. Each path is labeled with the columns it tests.
#'   \item \strong{Otherwise -> DROP} -- terminal.
#' }
#'
#' Vector-valued constraints are listed in full (no abbreviation).
#'
#' \strong{v0.30.1:} the previous graphical engines (Graphviz / ggplot)
#' were removed in favor of text-only rendering. Text reads better at
#' the scheme sizes we use, copy-pastes into SOPs and logs cleanly, and
#' the AND/OR structure is obvious from indentation.
#'
#' @param scheme A \code{ctdna_filtering_scheme} object, or the name of
#'   a built-in / user-registered scheme (string).
#' @param width Integer; description / wrap width in characters.
#'   Default 80.
#' @param print If \code{TRUE} (default), \code{cat()} the rendering to
#'   stdout AND invisibly return the character vector. If \code{FALSE},
#'   just return the character vector (one element per line).
#' @return Invisibly: a character vector with one element per output
#'   line (suitable for \code{writeLines()}).
#' @seealso \code{\link{create_filtering_scheme}},
#'   \code{\link{ctdna_filter_scheme_show}}.
#' @examples
#' \dontrun{
#' visualized_scheme("HRR14")
#' visualized_scheme("TSG_Tier2")
#'
#' my <- create_filtering_scheme(
#'   rule(Gene = c("BRCA1","BRCA2","ATM")),
#'   not(rule_clinvar_benign()),
#'   anyOf(rule_clinvar_path(), rule_truncating()),
#'   name = "my_HRR")
#' visualized_scheme(my)
#'
#' # As character vector for writing to a file:
#' lines <- visualized_scheme("HRR14", print = FALSE)
#' writeLines(lines, "HRR14.txt")
#' }
#' @export
visualized_scheme <- function(scheme, width = 80, print = TRUE) {
  # Resolve scheme arg
  if (is.character(scheme) && length(scheme) == 1L) {
    sch_name <- scheme
    builtin  <- tryCatch(.ctdna_builtin_schemes(), error = function(e) list())
    .ctdna_init_schemes_v1()
    user     <- .env[["filter_schemes_v1"]] %||% list()
    if (sch_name %in% names(builtin))   scheme <- builtin[[sch_name]]
    else if (sch_name %in% names(user)) scheme <- user[[sch_name]]
    else stop(sprintf("No scheme named '%s'. ", sch_name),
              "See ctdna_filter_schemes_list() for the catalog.",
              call. = FALSE)
  }
  if (!inherits(scheme, "ctdna_filtering_scheme"))
    stop("`scheme` must be a ctdna_filtering_scheme object or a scheme name.",
          call. = FALSE)

  lines <- .render_scheme_text(scheme, width = width)
  if (isTRUE(print)) cat(paste(lines, collapse = "\n"), "\n", sep = "")
  invisible(lines)
}
