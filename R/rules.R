# =============================================================================
# v1.0.0 — Composable rule system for filtering ctDNA alteration calls
# =============================================================================
#
# Conceptual model
# ----------------
# An "atomic rule" is a constraint on one or more columns of an alteration
# data frame. e.g. `rule(ClinVar = "Pathogenic")` matches rows where the
# ClinVar column equals "Pathogenic".
#
# Rules compose into schemes via three combinators:
#   allOf(r1, r2, ...) — row passes if ALL of r1, r2, ... pass (AND)
#   anyOf(r1, r2, ...) — row passes if ANY  of r1, r2, ... pass (OR)
#   not(r1)             — row passes if r1 does NOT pass
#
# Equivalent infix operators are provided for fluency: `r1 & r2`, `r1 | r2`,
# `!r1`.
#
# A "filtering scheme" is a composed rule tree returned by
# create_filtering_scheme(), which AND's its top-level arguments.
# Evaluation produces a logical keep-mask of length nrow(df).
#
# All rule objects share the S3 class "ctdna_rule" and have an evaluator
# .eval_rule(rule, df, ...) that returns a logical vector.

# -----------------------------------------------------------------------------
# Atomic rule constructor
# -----------------------------------------------------------------------------

#' Create an atomic filtering rule
#'
#' An atomic rule is a constraint on one or more columns. All constraints
#' within a single \code{rule()} call are AND'd together — every named
#' argument must be satisfied for the rule to match a row.
#'
#' Each named argument can be:
#' \itemize{
#'   \item A single character / numeric value -- equality match
#'     (case-insensitive for strings).
#'   \item A character vector -- value-in-set match.
#'   \item A list with \code{op} and \code{value} -- operator match.
#'     Operators supported: \code{">"}, \code{">="}, \code{"<"},
#'     \code{"<="}, \code{"=="}, \code{"!="}, \code{"%in%"},
#'     \code{"not_in"}.
#' }
#'
#' Argument names refer to columns of the data frame the rule will be
#' evaluated against. Column references can be remapped via
#' \code{ctdna_opts(col_<name> = "<actual_column>")} for vendors using
#' different column headers.
#'
#' NA handling: NA values in the target column produce FALSE (NA does
#' not match any value). The scheme evaluator wraps this to keep NA-safe
#' filtering at the top level — see \code{ctdna_variant_filter}.
#'
#' @param ... Named column constraints. See Details for accepted forms.
#' @return A \code{ctdna_rule} object — a small list with class
#'   \code{c("ctdna_rule_atom", "ctdna_rule")} that other combinators
#'   compose into trees.
#' @seealso \code{\link{allOf}}, \code{\link{anyOf}}, \code{\link{not}},
#'   \code{\link{create_filtering_scheme}}.
#' @examples
#' # Match rows where ClinVar is Pathogenic
#' r1 <- rule(ClinVar = "Pathogenic")
#'
#' # Match rows where ClinVar is any of three labels AND Variant_type is SNV
#' r2 <- rule(ClinVar = c("Pathogenic","Likely_pathogenic"),
#'             Variant_type = "SNV")
#'
#' # Match rows where gnomAD_AF < 0.01
#' r3 <- rule(gnomAD_AF = list(op = "<", value = 0.01))
#' @export
rule <- function(...) {
  conds <- list(...)
  if (length(conds) == 0)
    stop("`rule()` requires at least one column = value pair.")
  if (is.null(names(conds)) || any(names(conds) == ""))
    stop("Every argument to `rule()` must be named (the column to match against).")
  structure(list(type = "atom", conditions = conds),
            class = c("ctdna_rule_atom", "ctdna_rule"))
}


# -----------------------------------------------------------------------------
# Combinators: allOf / anyOf / not
# -----------------------------------------------------------------------------

#' Compose rules with logical AND (\code{allOf}) or OR (\code{anyOf})
#'
#' \code{allOf(...)} matches a row only if EVERY argument matches it.
#' \code{anyOf(...)} matches a row if AT LEAST ONE argument matches it.
#' Arguments may be atomic rules (from \code{\link{rule}}) or other
#' combinator nodes — they nest freely.
#'
#' @param ... One or more \code{ctdna_rule} objects (atoms or combinators).
#' @return A \code{ctdna_rule} combinator object.
#' @seealso \code{\link{not}}, \code{\link{rule}},
#'   \code{\link{create_filtering_scheme}}.
#' @examples
#' r_path <- rule(ClinVar = c("Pathogenic","Likely_pathogenic"))
#' r_rare <- rule(gnomAD_AF = list(op = "<", value = 0.01))
#' both   <- allOf(r_path, r_rare)
#' either <- anyOf(r_path, r_rare)
#' @name combinators
#' @export
allOf <- function(...) {
  args <- list(...)
  if (length(args) == 0)
    stop("`allOf()` requires at least one rule.")
  for (a in args)
    if (!inherits(a, "ctdna_rule"))
      stop("All arguments to `allOf()` must be ctdna_rule objects.")
  structure(list(type = "and", children = args),
            class = c("ctdna_rule_combinator", "ctdna_rule"))
}

#' @rdname combinators
#' @examples
#' # Equivalent to rule_or(rule(Gene = "TP53"), rule(Gene = "KRAS"))
#' anyOf(Gene = c("TP53", "KRAS", "EGFR"))
#'
#' # Use inside a scheme
#' create_filtering_scheme(
#'   anyOf(Variant_type = c("SNV", "Indel")),
#'   name = "snv_or_indel")
#' @export
anyOf <- function(...) {
  args <- list(...)
  if (length(args) == 0)
    stop("`anyOf()` requires at least one rule.")
  for (a in args)
    if (!inherits(a, "ctdna_rule"))
      stop("All arguments to `anyOf()` must be ctdna_rule objects.")
  structure(list(type = "or", children = args),
            class = c("ctdna_rule_combinator", "ctdna_rule"))
}

#' Negate a rule
#'
#' \code{not(r)} matches rows where \code{r} does NOT match.
#'
#' NA handling: if \code{r} returns FALSE for a row because the underlying
#' column was NA, \code{not(r)} returns NA (not TRUE) — we don't claim a
#' row is "the negation of pathogenic" when ClinVar is missing. The
#' scheme evaluator coerces these NAs to FALSE at the top level (rows
#' with missing data fall through to general / basic).
#'
#' @param r A \code{ctdna_rule} object.
#' @return A \code{ctdna_rule} combinator object.
#' @seealso \code{\link{allOf}}, \code{\link{anyOf}}, \code{\link{rule}}.
#' @examples
#' not(rule(ClinVar = c("Benign","Likely_benign")))
#' @export
not <- function(r) {
  if (!inherits(r, "ctdna_rule"))
    stop("`not()` requires a ctdna_rule object.")
  structure(list(type = "not", child = r),
            class = c("ctdna_rule_combinator", "ctdna_rule"))
}


# -----------------------------------------------------------------------------
# Infix operators for fluency
# -----------------------------------------------------------------------------

#' Logical operators for ctdna_rule objects
#'
#' \code{`&`}, \code{`|`}, and \code{`!`} are overloaded for
#' \code{ctdna_rule} objects, providing fluent infix syntax equivalent
#' to \code{\link{allOf}}, \code{\link{anyOf}}, and \code{\link{not}}.
#'
#' @param e1,e2 ctdna_rule objects.
#' @return A combinator ctdna_rule object.
#' @examples
#' r_path <- rule(ClinVar = c("Pathogenic","Likely_pathogenic"))
#' r_rare <- rule(gnomAD_AF = list(op = "<", value = 0.01))
#' r_path & r_rare        # equivalent to allOf(r_path, r_rare)
#' r_path | r_rare        # equivalent to anyOf(r_path, r_rare)
#' !r_path                # equivalent to not(r_path)
#' @name rule-operators
NULL

#' @rdname rule-operators
#' @export
`&.ctdna_rule` <- function(e1, e2) allOf(e1, e2)

#' @rdname rule-operators
#' @export
`|.ctdna_rule` <- function(e1, e2) anyOf(e1, e2)

#' @rdname rule-operators
#' @export
`!.ctdna_rule` <- function(e1) not(e1)


# -----------------------------------------------------------------------------
# Scheme constructor — AND's its top-level arguments
# -----------------------------------------------------------------------------

#' Create a filtering scheme from one or more rules
#'
#' A filtering scheme is a composed rule tree. \code{create_filtering_scheme()}
#' AND's its top-level arguments together — every argument must match for
#' a row to pass the scheme. Use \code{anyOf()} inside to express OR'd
#' success paths and \code{not()} to negate.
#'
#' Schemes can also carry a \code{name} (used for catalog listing) and a
#' \code{description} (printed by \code{ctdna_filter_scheme_show}).
#'
#' @param ... One or more \code{ctdna_rule} objects (atoms or combinators).
#'   All are AND'd together.
#' @param name Optional scheme name (string). If supplied, the scheme is
#'   registered in the user-scheme catalog and is referenceable by name
#'   in \code{filter_scheme = "name"}.
#' @param description Optional human description.
#' @param category One of \code{"gene_set"}, \code{"indication"},
#'   \code{"indication_gene_set"}, \code{"mutation_type"}, \code{"user"}.
#'   Affects which step of the apply hierarchy claims rows; see
#'   \code{\link{ctdna_variant_filter}}. Default \code{"user"}.
#' @param scope Optional \code{ctdna_rule} defining the "scope" of the
#'   scheme — which rows the scheme has an opinion about. Rows matching
#'   the scope are CLAIMED by this scheme: if they also pass the full
#'   rule tree they survive, otherwise they are DROPPED (not handed to
#'   general). Rows OUTSIDE the scope are ignored by this scheme and
#'   fall through to general / next scheme. Default: NULL, meaning the
#'   first top-level positional argument is used as the scope (for
#'   built-in HRR/TSG schemes that's the gene-set rule).
#' @return A \code{ctdna_filtering_scheme} object (also class
#'   \code{ctdna_rule} so it can be nested inside larger schemes).
#' @seealso \code{\link{rule}}, \code{\link{allOf}}, \code{\link{anyOf}},
#'   \code{\link{not}}, \code{\link{ctdna_variant_filter}}.
#' @examples
#' # A 2-line scheme: HRR genes only, drop benign, pathogenic OR truncating-rare.
#' \dontrun{
#' my_HRR <- create_filtering_scheme(
#'   rule(Gene = c("BRCA1","BRCA2","ATM")),
#'   not(rule(ClinVar = c("Benign","Likely_benign"))),
#'   anyOf(
#'     rule(ClinVar = c("Pathogenic","Likely_pathogenic")),
#'     allOf(rule(Molecular_consequence = c("nonsense","frameshift")),
#'             rule(gnomAD_AF = list(op = "<", value = 0.01)))),
#'   name = "my_HRR", category = "gene_set",
#'   description = "HRR with my custom thresholds")
#' }
#' @export
#' Create a filtering scheme (unified creator)
#'
#' The one and only filtering-scheme creator in ctdnaTM. Combines any
#' number of rule-tree expressions (\code{rule()}, \code{allOf()},
#' \code{anyOf()}, \code{not()}) into a named scheme, optionally
#' cloning an existing scheme as a starting point.
#'
#' \strong{v0.36.0:} This is now the only entry point for creating
#' schemes. The older \code{ctdna_filter_scheme_create()} (v0.26
#' gene_sets + overrides container) has been removed -- its template
#' machinery and frozen/overwrite/reserved-name handling have been
#' folded into this function.
#'
#' Positional rule arguments are AND'd together to form the root
#' logic. If \code{template} is given, the rules in this call are
#' AND'd ON TOP OF the template's rule tree (template first, new
#' rules second).
#'
#' @section Templates:
#' \code{template} can name any existing scheme. Three sources are
#' searched in this order:
#' \enumerate{
#'   \item v0.26 built-in indication templates (NSCLC, HNSCC, SCLC,
#'     BRCA, CRC, mCRPC, GBM). These contain a curated driver-gene
#'     list per indication. The template's \code{Drivers} list is
#'     converted to a single \code{rule(Gene = c(...))} that becomes
#'     the base. Any per-set overrides (currently always empty in
#'     the built-ins) are dropped.
#'   \item v1.0.0 built-in schemes (\code{HRR14},
#'     \code{TSG_Tier2}). The complete rule tree, including scope and
#'     all branches (germline/somatic/LGR/CNV paths), is taken as the
#'     base.
#'   \item User-defined schemes already registered by an earlier call
#'     to \code{create_filtering_scheme()}.
#' }
#' Unknown template names produce a clear error listing all
#' available templates across all three registries.
#'
#' @section Registration:
#' If \code{name} is supplied (the typical case), the scheme is
#' registered in the user scheme registry under that name and can be
#' retrieved later by name with \code{\link{ctdna_filter_scheme_get}}.
#' If \code{name} is \code{NULL}, the scheme is returned but not
#' stored anywhere -- useful for one-off evaluation.
#'
#' \strong{Reserved names:} \code{basic}, \code{default}, \code{none},
#' and the names of all built-in templates and schemes are protected.
#' To customize a built-in, pick a different name and pass the
#' built-in via \code{template}.
#'
#' \strong{Overwrite:} re-using an existing user-scheme name errors
#' unless \code{overwrite = TRUE} is set. \code{frozen} schemes
#' cannot be overwritten regardless of the \code{overwrite} flag.
#'
#' @param ... Zero or more \code{ctdna_rule} objects. Positional rules
#'   are AND'd together. Combined with the template's rule tree (if
#'   any) by AND. At least one rule must be supplied across the
#'   positional args + template -- a scheme with no constraints is
#'   not allowed.
#' @param name Scheme name (character). When supplied, the scheme is
#'   registered. NULL = build but do not register.
#' @param description Free-text description of what the scheme does.
#' @param category One of \code{"user"} (default), \code{"gene_set"},
#'   \code{"indication"}, \code{"indication_gene_set"},
#'   \code{"mutation_type"}. Free metadata label used for filtering
#'   in \code{\link{ctdna_filter_schemes_list}}.
#' @param scope A single rule that restricts which rows the scheme
#'   applies to (e.g. \code{rule(Gene = c("BRCA1","BRCA2"))}). Rows
#'   outside the scope are handled by the active general filter. When
#'   NULL, the first positional rule (or the template's scope, if a
#'   template is used) becomes the scope. Pass \code{scope = NA} to
#'   force "scope = all rows".
#' @param template Name of an existing scheme to use as the base. See
#'   \dQuote{Templates} above.
#' @param frozen If \code{TRUE}, the scheme is marked read-only and
#'   future calls with the same \code{name} and \code{overwrite = TRUE}
#'   will be refused. Use for SOP-locked or validated schemes.
#' @param overwrite If \code{TRUE}, allow replacing an existing
#'   (non-frozen) scheme registered under \code{name}. Default
#'   \code{FALSE} so accidental name collisions error loudly.
#' @return A \code{ctdna_filtering_scheme} object.
#' @seealso \code{\link{rule}}, \code{\link{allOf}}, \code{\link{anyOf}},
#'   \code{\link{not}}, \code{\link{ctdna_variant_filter}},
#'   \code{\link{ctdna_filter_scheme_get}},
#'   \code{\link{ctdna_filter_schemes_list}},
#'   \code{\link{ctdna_filter_scheme_show}}.
#' @examples
#' \dontrun{
#' # Build a fresh scheme from rules
#' my_scheme <- create_filtering_scheme(
#'   rule(Gene = c("BRCA1","BRCA2")),
#'   rule(ClinVar = c("Pathogenic","Likely_pathogenic")),
#'   name = "my_brca",
#'   description = "BRCA1/2 P/LP only")
#'
#' # Clone NSCLC's curated driver list and add a rarity filter
#' my_nsclc <- create_filtering_scheme(
#'   rule(gnomAD_AF = list(op = "<", value = 0.01)),
#'   name = "my_NSCLC_rare",
#'   template = "NSCLC")
#'
#' # Clone TSG_Tier2 and AND-in a VAF floor
#' my_tsg <- create_filtering_scheme(
#'   rule(VAF = list(op = ">=", value = 0.005)),
#'   name = "my_TSG_vaf_floor",
#'   template = "TSG_Tier2",
#'   frozen = TRUE)
#' }
#' @export
create_filtering_scheme <- function(...,
                                     name        = NULL,
                                     description = NULL,
                                     category    = c("user","gene_set",
                                                      "indication",
                                                      "indication_gene_set",
                                                      "mutation_type"),
                                     scope       = NULL,
                                     branches    = NULL,
                                     template    = NULL,
                                     frozen      = FALSE,
                                     overwrite   = FALSE) {
  category <- match.arg(category)
  args <- list(...)
  for (a in args)
    if (!inherits(a, "ctdna_rule"))
      stop("All positional arguments to `create_filtering_scheme()` must be ctdna_rule objects.",
            call. = FALSE)

  # ---- Validate branches (optional) ---------------------------------------
  if (!is.null(branches)) {
    if (!is.list(branches) || is.null(names(branches)) ||
        any(!nzchar(names(branches))))
      stop("`branches` must be a NAMED list of ctdna_rule objects.",
           call. = FALSE)
    for (b in branches)
      if (!inherits(b, "ctdna_rule"))
        stop("Every element of `branches` must be a ctdna_rule.", call. = FALSE)
  }

  # ---- Resolve template if provided ---------------------------------------
  tmpl_rule       <- NULL
  tmpl_scope      <- NULL
  tmpl_descr      <- NULL
  if (!is.null(template)) {
    if (!is.character(template) || length(template) != 1L)
      stop("`template` must be a single scheme name.", call. = FALSE)
    .ctdna_init_schemes_v1()
    builtin_v026 <- tryCatch(.ctdna_filter_builtin_templates(),
                               error = function(e) list())
    builtin_v1   <- tryCatch(.ctdna_builtin_schemes(),
                               error = function(e) list())
    user_v1      <- .env[["filter_schemes_v1"]] %||% list()

    if (template %in% names(builtin_v026)) {
      # v0.26 indication template — convert gene_sets to a rule tree
      t <- builtin_v026[[template]]
      all_genes <- unlist(t$gene_sets, use.names = FALSE)
      if (length(all_genes) == 0L)
        stop(sprintf("Template '%s' has no genes.", template), call. = FALSE)
      tmpl_rule  <- rule(Gene = unique(all_genes))
      tmpl_scope <- tmpl_rule
      tmpl_descr <- t$description %||%
        sprintf("Cloned from v0.26 template '%s'.", template)
    } else if (template %in% names(builtin_v1)) {
      sch <- builtin_v1[[template]]
      tmpl_rule  <- sch$root
      tmpl_scope <- sch$scope
      tmpl_descr <- sch$description
    } else if (template %in% names(user_v1)) {
      sch <- user_v1[[template]]
      tmpl_rule  <- sch$root
      tmpl_scope <- sch$scope
      tmpl_descr <- sch$description
    } else {
      all_avail <- unique(c(names(builtin_v026), names(builtin_v1),
                             names(user_v1)))
      stop(sprintf("Unknown template '%s'. Available templates: %s",
                    template, paste(all_avail, collapse = ", ")),
            call. = FALSE)
    }
  }

  # ---- Combine template + new rules + branches into root ------------------
  # If `branches` is given, treat it as an anyOf() body ANDed with the other
  # positional rules. The branches list is stored on the scheme for use by
  # the annotation engine.
  branches_rule <- if (!is.null(branches)) do.call(anyOf, unname(branches)) else NULL
  combined <- c(if (!is.null(tmpl_rule)) list(tmpl_rule),
                args,
                if (!is.null(branches_rule)) list(branches_rule))
  if (length(combined) == 0L)
    stop("`create_filtering_scheme()` requires at least one rule (positional ",
          ", `branches`, or `template`).", call. = FALSE)
  root <- if (length(combined) == 1L) combined[[1]] else do.call(allOf, combined)

  # ---- Resolve scope -------------------------------------------------------
  # Precedence: explicit scope arg > template scope > first positional rule
  # `scope = NA` means "force scope = NULL (no row restriction)"
  scope_rule <- if (!is.null(scope) && !identical(scope, NA)) scope
                else if (identical(scope, NA))                NULL
                else if (!is.null(tmpl_scope))                tmpl_scope
                else if (length(args) > 0)                    args[[1]]
                else                                          NULL

  # ---- Description: caller wins, else template description ----------------
  if (is.null(description)) description <- tmpl_descr

  # ---- Build object --------------------------------------------------------
  scheme <- structure(
    list(type        = "scheme",
         name        = name,
         description = description,
         category    = category,
         root        = root,
         scope       = scope_rule,
         branches    = branches,     # named list of ctdna_rule (or NULL)
         gates       = args,         # positional args = REQUIRE gates
         frozen      = isTRUE(frozen)),
    class = c("ctdna_filtering_scheme", "ctdna_rule"))

  # ---- Register if name supplied -------------------------------------------
  if (!is.null(name)) {
    .ctdna_register_user_scheme(name, scheme, overwrite = isTRUE(overwrite))
  }
  scheme
}


# -----------------------------------------------------------------------------
# Evaluator: turn a rule (tree) into a logical keep-mask over a data frame
# -----------------------------------------------------------------------------

# .eval_rule(): returns a logical vector of length nrow(df). NA propagates
# through atom evaluation (NA in target column -> NA in mask); combinators
# follow standard logic-with-NA rules (FALSE & NA = FALSE, TRUE | NA = TRUE).
# The top-level caller (.ctdna_apply_scheme) coerces remaining NAs to FALSE.
.eval_rule <- function(node, df) {
  if (!inherits(node, "ctdna_rule"))
    stop(".eval_rule: not a ctdna_rule object.")

  if (inherits(node, "ctdna_filtering_scheme"))
    return(.eval_rule(node$root, df))

  switch(node$type,
    "atom" = .eval_atom(node, df),
    "and"  = .eval_and(node, df),
    "or"   = .eval_or(node, df),
    "not"  = .eval_not(node, df),
    "true" = rep(TRUE,  nrow(df)),
    "false"= rep(FALSE, nrow(df)),
    stop(sprintf(".eval_rule: unknown node type '%s'", node$type)))
}

# Single-atom evaluation: AND all named conditions together.
.eval_atom <- function(node, df) {
  conds <- node$conditions
  n <- nrow(df)
  if (n == 0) return(logical(0))
  acc <- rep(TRUE, n)
  for (col_name in names(conds)) {
    actual_col <- .ctdna_resolve_col(col_name, df)
    if (is.null(actual_col)) {
      # Column not present -> rule silently doesn't match (FALSE).
      # This is the right behavior for atoms; the scheme evaluator
      # handles "column missing" at a higher level if needed.
      acc <- acc & FALSE
      next
    }
    col_vals <- df[[actual_col]]
    constraint <- conds[[col_name]]
    mask <- .eval_constraint(col_vals, constraint, df = df)
    acc <- acc & mask
  }
  acc
}

# Resolve a column reference: first check ctdna_opts(col_<name>) mapping,
# then fall back to the literal name.
.ctdna_resolve_col <- function(col_name, df) {
  # Check opts mapping first (e.g. ClinVar -> ctdna_opts$col_clinvar)
  opt_key <- paste0("col_", tolower(col_name))
  mapped <- tryCatch(.o(opt_key), error = function(e) NULL)
  if (!is.null(mapped) && length(mapped) == 1L && is.character(mapped) &&
      mapped %in% names(df)) {
    return(mapped)
  }
  # Fall back to literal column name
  if (col_name %in% names(df)) return(col_name)
  NULL
}

# Apply one constraint to a vector. Returns a logical vector of same length.
.eval_constraint <- function(vec, constraint, df = NULL) {
  # Operator form: list(op, value)
  if (is.list(constraint) && all(c("op","value") %in% names(constraint))) {
    op  <- constraint$op
    val <- constraint$value
    # Column-vs-column equality: op == "==col" or "!=col"
    # `value` is interpreted as the name of another column in df.
    # Used for intragenic LGR checks (e.g. Gene == Fusion_gene_b).
    if (op %in% c("==col","!=col")) {
      if (is.null(df))
        stop(".eval_constraint: column-vs-column op requires df.")
      other_col <- .ctdna_resolve_col(val, df)
      if (is.null(other_col)) {
        # Partner column absent — row cannot satisfy the comparison.
        return(rep(FALSE, length(vec)))
      }
      other_vec <- df[[other_col]]
      if (is.character(vec) || is.factor(vec) ||
           is.character(other_vec) || is.factor(other_vec)) {
        a <- tolower(trimws(as.character(vec)))
        b <- tolower(trimws(as.character(other_vec)))
        hit <- a == b
      } else {
        hit <- vec == other_vec
      }
      hit[is.na(hit)] <- FALSE
      return(if (op == "==col") hit else !hit)
    }
    return(.eval_op(vec, op, val))
  }
  # Plain value(s): case-insensitive %in% match for strings; equality for nums
  if (is.character(vec) || is.factor(vec)) {
    # If the constraint is numeric and the vec is character-like, try
    # numeric coercion first (handles e.g. Direction_a stored as "-1"/"1"
    # vs constraint integer -1/1). Fall back to string match if coercion
    # produces all NAs.
    if (is.numeric(constraint)) {
      vec_num <- suppressWarnings(as.numeric(as.character(vec)))
      if (!all(is.na(vec_num))) {
        mask <- vec_num %in% constraint
        mask[is.na(vec_num)] <- FALSE
        return(mask)
      }
    }
    vec_lower <- tolower(trimws(as.character(vec)))
    val_lower <- tolower(trimws(as.character(constraint)))
    return(vec_lower %in% val_lower)
  }
  # Numeric: %in% (or equality if single value)
  vec %in% constraint
}

# Apply an operator. Numeric ops work elementwise; %in% / not_in work on
# any vector type with case-insensitive matching for strings.
.eval_op <- function(vec, op, val) {
  if (op == "%in%" || op == "not_in") {
    if (is.character(vec) || is.factor(vec)) {
      hit <- tolower(trimws(as.character(vec))) %in%
              tolower(trimws(as.character(val)))
    } else {
      hit <- vec %in% val
    }
    return(if (op == "%in%") hit else !hit)
  }
  # Numeric / comparison ops
  switch(op,
    "==" = vec == val,
    "!=" = vec != val,
    ">"  = vec >  val,
    ">=" = vec >= val,
    "<"  = vec <  val,
    "<=" = vec <= val,
    stop(sprintf(".eval_op: unknown operator '%s'", op)))
}

.eval_and <- function(node, df) {
  n <- nrow(df)
  if (n == 0) return(logical(0))
  if (length(node$children) == 0) return(rep(TRUE, n))
  m <- .eval_rule(node$children[[1]], df)
  for (i in seq_along(node$children)[-1]) {
    m <- m & .eval_rule(node$children[[i]], df)
  }
  m
}

.eval_or <- function(node, df) {
  n <- nrow(df)
  if (n == 0) return(logical(0))
  if (length(node$children) == 0) return(rep(FALSE, n))
  m <- .eval_rule(node$children[[1]], df)
  for (i in seq_along(node$children)[-1]) {
    m <- m | .eval_rule(node$children[[i]], df)
  }
  m
}

.eval_not <- function(node, df) !.eval_rule(node$child, df)


# -----------------------------------------------------------------------------
# Printing / inspection
# -----------------------------------------------------------------------------

#' Print method for ctdna_rule objects
#'
#' Renders a rule tree as indented text, with combinator types
#' (\code{allOf}, \code{anyOf}, \code{not}) shown as parent nodes and
#' atomic rules as leaves with their column/value/operator details.
#'
#' @param x A \code{ctdna_rule} object.
#' @param indent Integer. Current indent depth in spaces (used during
#'   recursion). Callers usually leave as the default \code{0}.
#' @param ... Unused; present for S3 method consistency.
#' @return Invisibly returns \code{x}.
#' @export
print.ctdna_rule <- function(x, indent = 0, ...) {
  pad <- strrep("  ", indent)
  if (inherits(x, "ctdna_filtering_scheme")) {
    cat(sprintf("%s<filtering_scheme%s%s>\n", pad,
                if (!is.null(x$name)) sprintf(" name='%s'", x$name) else "",
                if (!is.null(x$category)) sprintf(" category='%s'", x$category) else ""))
    if (!is.null(x$description)) cat(sprintf("%s  desc: %s\n", pad, x$description))
    print(x$root, indent = indent + 1)
    return(invisible(x))
  }
  if (inherits(x, "ctdna_rule_atom")) {
    parts <- vapply(names(x$conditions), function(nm) {
      v <- x$conditions[[nm]]
      val_str <- if (is.list(v) && all(c("op","value") %in% names(v)))
        sprintf("%s %s", v$op, paste(v$value, collapse = ",")) else
        paste(v, collapse = ",")
      sprintf("%s = %s", nm, val_str)
    }, character(1))
    cat(sprintf("%srule: %s\n", pad, paste(parts, collapse = "; ")))
    return(invisible(x))
  }
  if (inherits(x, "ctdna_rule_combinator")) {
    cat(sprintf("%s%s:\n", pad, switch(x$type,
       "and" = "allOf",
       "or"  = "anyOf",
       "not" = "not")))
    if (x$type == "not") {
      print(x$child, indent = indent + 1)
    } else {
      for (c in x$children) print(c, indent = indent + 1)
    }
    return(invisible(x))
  }
  cat(pad, "<unknown ctdna_rule>\n")
  invisible(x)
}


# -----------------------------------------------------------------------------
# User scheme catalog registration
# -----------------------------------------------------------------------------

.ctdna_register_user_scheme <- function(name, scheme, overwrite = FALSE) {
  if (!nzchar(name))
    stop("Scheme `name` cannot be empty.", call. = FALSE)
  reserved <- tryCatch(.ctdna_reserved_scheme_names(),
                        error = function(e) c("basic","none"))
  if (name %in% reserved)
    stop(sprintf("Scheme name '%s' is reserved. Pick a different name; ", name),
          sprintf("to clone the built-in, use template = '%s'.", name),
          call. = FALSE)
  .ctdna_init_schemes_v1()
  existing <- .env[["filter_schemes_v1"]][[name]]
  if (!is.null(existing)) {
    if (isTRUE(existing$frozen))
      stop(sprintf("Scheme '%s' is frozen and cannot be overwritten. ", name),
            "Pick a different name.", call. = FALSE)
    if (!isTRUE(overwrite))
      stop(sprintf("Scheme '%s' already exists. ", name),
            "Pass `overwrite = TRUE` to replace it.", call. = FALSE)
  }
  .env[["filter_schemes_v1"]][[name]] <- scheme
  invisible(name)
}

# Init the v1.0.0 scheme registry (separate from legacy filter_schemes)
.ctdna_init_schemes_v1 <- function() {
  if (is.null(.env[["filter_schemes_v1"]]))
    .env[["filter_schemes_v1"]] <- list()
}
