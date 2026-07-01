# ============================================================================
# ctdnaTM — configuration, helpers, theme, colors, example data
# ============================================================================

# Private env holding the live config
.env <- new.env(parent = emptyenv())

# Dose palette. Matches the source deck's baseline / reduction plots
# (blue / yellow / gray for 2.8 / 4.8 / 5.8 mg/kg). Same primaries as
# the RECIST palette by design — when plots are stratified by dose OR
# RECIST (not both), the colors are consistent across the deck.
.dose_pal <- c("#1F77B4", "#E6B800", "#7F7F7F", "#4A4A4A", "#D62728")

# Factory defaults
.defaults <- list(
  # ---- v0.41.0: Clinical / subject-level columns ----
  # Defaults match Guardant Health Infinity column names. Legacy short
  # keys (subject, dose, recist, mr, best_change) are unchanged in spelling
  # but their default values now reflect vendor names.
  subject       = "Patient_ID",            # v0.41.0: was "subject_id"
  time          = "Visit_name",            # v0.41.0: was "time_point"; matches vendor
  baseline      = "Baseline",
  dose          = "Dose",                  # v0.41.0: was "dose"
  cohort        = "cohort",
  recist        = "RECIST",
  mr            = "MR_score_percentage",   # v0.44.0: GH header, percent 0-100
  best_change   = "best_pct_change",
  lesion        = "Lesion_sum_mm",         # v0.50.0: per-visit lesion size (SLD)
  lesion_pct    = "Tumor_size_pct_change",  # v0.51.0: per-visit tumour-size %change (ADTR PCHG)
  cancertype    = "Cancertype",            # v0.41.0: new — vendor indication column

  # ---- ctDNA per-sample columns (mostly already vendor-named) ----
  tf            = "Methylation_tumor_fraction_percentage",  # v0.44.0: GH header, percent 0-100
  maxvaf        = "Genomic_max_VAF_percentage",             # v0.44.0: GH header, percent 0-100
  meanvaf       = "Mean_VAF_percentage",                    # v0.44.0: derived per-sample mean, percent
  cfdna_conc    = "cfDNA_ng",                               # v0.44.0: GH header
  n_somatic     = "n.somatic",                              # derived per-sample count
  sample        = "Sample_ID",             # v0.41.0: new — vendor sample identifier
  visit         = "Visit_name",            # v0.41.0: new — vendor visit identifier

  # ---- Genomic-modality columns (vendor names) ----
  gene          = "Gene",                  # v0.41.0: was "gene"
  alteration    = "Variant_type",          # v0.41.0: was "alteration_type"; now matches vendor
  mutation      = "mutation_status",       # derived by prepare (mut / wt)
  panel         = "panel",                 # derived by prepare
  variant_type  = "Variant_type",          # v0.41.0: new — explicit vendor key
  molecular_consequence = "Molecular_consequence",   # v0.41.0: new
  cnv_type      = "CNV_type",              # v0.41.0: new
  copy_number   = "Copy_number",           # v0.41.0: new
  vaf           = "VAF_percentage",        # v0.44.0: GH per-variant header, percent 0-100

  # ---- TMB / MSI (vendor names) ----
  tmb           = "TMB_score",             # v0.41.0: was "tmb"
  tmb_category  = "TMB_category",          # v0.41.0: new
  msi           = "MSI_High",              # v0.41.0: new

  # ---- IHC (vendor names) ----
  marker        = "Marker",                # v0.41.0: was "marker"
  ihc_score     = "Score",                 # v0.41.0: was "score"
  ihc_score_unit = "Score_unit",           # v0.41.0: new

  # ---- Expression (derived by prepare from RNA-seq matrix) ----
  expression    = "expression",            # value column (derived)
  expression_unit = "log2_tpm_plus_one",

  # ---- Other ----
  ctdna_metric  = "best_ratio",            # derived

  # ---- Numeric parameters ----
  display_floor = 0.1,                     # v0.44.0: percent scale (was 0.001 fraction)
  loq           = 1,                       # v0.44.0: percent scale (was 0.01 fraction)
  discord_log10 = 1,

  # ---- v0.44.0: per-sample QC thresholds (used by Sample_QC in prepare,
  #      and by ctdna_sample_qc()). plasma/cfDNA mins drop samples below
  #      threshold; status keeps only matching Sample_status (NULL = off). ----
  sample_qc_status      = c("SUCCESS", "Success", "success", "PASS", "Pass", "pass"),  # acceptable Sample_status; NULL = no status filter
  sample_qc_plasma_min  = 0,                          # min Plasma_ml_input  (0 = off)
  sample_qc_cfdna_min   = 0,                           # min cfDNA_ng         (0 = off)
  col_plasma_input      = "Plasma_ml_input",
  col_cfdna_ng          = "cfDNA_ng",

  # ---- RECIST value mapping ----
  recist_crpr   = c("CR", "PR", "CR/PR"),
  recist_ucrpr  = c("uCR", "uPR", "uCR/PR"),
  recist_sd     = c("SD"),
  recist_pdnena = c("PD", "NE", "NA"),

  # ---- v1.0.0: canonical column-name mapping for the rule system ----
  # Looked up by .ctdna_resolve_col() at rule-evaluation time. These ARE
  # vendor names and don't change in v0.41.0.
  col_clinvar              = "ClinVar",
  col_origin               = "Somatic_status",
  col_variant_type         = "Variant_type",
  col_cnv_type             = "CNV_type",
  col_splice_effect        = "Splice_effect",
  col_molecular_consequence= "Molecular_consequence",
  col_functional_impact    = "Functional_impact",
  col_mutant_allele_status = "Mutant_allele_status",
  col_gnomad_af            = "gnomAD_AF",
  col_gene                 = "Gene",
  col_mut_aa               = "Mut_aa",
  col_cancertype           = "Cancertype",
  col_sample_status        = "Sample_status",

  cancertype_dictionary = list(
    NSCLC = "NSCLC",
    BRCA  = "BRCA",
    SCLC  = "SCLC",
    HNSCC = "HNSCC",
    CRC   = "CRC",
    mCRPC = "mCRPC",
    GBM   = "GBM"),

  default_scheme = "basic",

  # ---- v0.37.0: sample-level QC filtering -------------------------------
  # Sample-level QC is applied INSIDE ctdna_variant_filter() in all three
  # filtering modes (including filter_scheme = NULL). The master switch
  # is sample_qc_filter; when TRUE, each configured QC column is checked
  # against its threshold and rows that fail are dropped.
  #
  # To disable a single check, set its column to NULL:
  #   ctdna_opts(sample_qc_cfdna_col = NULL)
  # To use different column names in your data:
  #   ctdna_opts(sample_qc_status_col = "QC_FLAG")
  # To skip QC entirely:
  #   ctdna_opts(sample_qc_filter = FALSE)
  #
  # If a configured QC column is not present in the input data frame,
  # that check is silently skipped (one-time message per session).
  sample_qc_filter      = TRUE,           # master switch
  sample_qc_status_col  = "Sample_status",
  sample_qc_status_pass = "SUCCESS",
  sample_qc_cfdna_col   = "cfDNA_ng",
  sample_qc_cfdna_min   = 5,
  sample_qc_plasma_col  = "Plasma_ml_input",
  sample_qc_plasma_min  = 0,                # v0.38.0: was 5; lowered to 0

  # ---- v0.38.0: rare-variant filtering switch -------------------------
  # When TRUE (default), rules that check gnomAD_AF (e.g. rule_rare_001
  # in HRR14 path B + TSG_Tier2 rare-pathogenic-synonymous carve-
  # out) use rare_variant_threshold as the cutoff. Set FALSE to disable
  # the rare-variant check entirely (the rule effectively always passes).
  rare_variant_filter    = TRUE,
  rare_variant_threshold = 0.01,

  # ---- v0.38.0: ADSL indication column priority -----------------------
  # When ctdna_make_patient_data() is called without an explicit
  # Cancer_Type frame, it derives Cancer_Type from input frame(s) by
  # trying these columns in order. ADaM names are listed first
  # (PRTUMTY -> COHORT -> STUDYTRT); non-ADaM canonical names
  # (Cancer_Type, Cancertype, Indication, Clin_Indication) come next
  # so user-supplied data frames work too.
  col_indication_adam = c("PRTUMTY","COHORT","STUDYTRT",
                            "Cancer_Type","Cancertype","Indication",
                            "Clin_Indication"),

  # ---- DEPRECATED in v0.26.0 (kept for backward compat through one cycle) ----
  # default_indication_schemes was the v0.22-0.25 mechanism for indication-
  # aware filtering. It's still read by the deprecated `filter_scheme =
  # "default"` alias but new code should use a vector of block names
  # instead, e.g. filter_scheme = c("NSCLC","HNSCC","SCLC",...).
  default_indication_schemes = list(
    NSCLC = "NSCLC", HNSCC = "HNSCC", SCLC = "SCLC",
    BRCA  = "BRCA",  CRC   = "CRC",   mCRPC = "mCRPC",
    GBM   = "GBM"
  ),

  # ---- Statistics defaults ----
  pairwise_test = "wilcox",
  cor_method    = "spearman",
  alpha         = 0.05,

  # ---- Plot display defaults ----
  point_size      = 2,
  point_alpha     = 0.85,
  box_width       = 0.6,            # consistent boxplot width across plots
  bar_width       = 0.55,           # consistent barplot width across plots
  show_stats      = TRUE,
  show_n          = TRUE,
  drop_na         = TRUE,           # v0.44.0: drop rows with NA grouping (e.g. Dose=NA) in plots
  median_color    = "grey20",
  clearance_color = "grey40",
  stat_position   = "subtitle",   # where stats labels appear:
                                  # "subtitle" / "caption" / "on_plot" / "none"
  legend_position = "right"       # ggplot legend position:
                                  # "right" / "left" / "top" / "bottom" / "none"

  # ---- v0.41.0: expression_unit, gene, expression, mutation, panel, marker,
  # ihc_score, alteration, ctdna_metric are defined in the vendor-named
  # defaults section at the top of this list.
)

.env$cfg <- .defaults

# ---------------------------------------------------------------------------
# Restricted-choice + numeric-range validators for ctdna_opts() keys.
# Each entry is either:
#   list(choices = c("a","b","c"))        - allowed string values
#   list(range   = c(min, max))           - allowed numeric range (inclusive)
#   list(type    = "logical")             - must be a TRUE/FALSE flag
# Keys not listed here accept any value (the column-name keys etc.).
# ---------------------------------------------------------------------------
.opts_validators <- list(
  pairwise_test    = list(choices = c("wilcox","t","kruskal","anova","auto")),
  cor_method       = list(choices = c("spearman","pearson","kendall")),
  stat_position    = list(choices = c("subtitle","caption","on_plot","none")),
  legend_position  = list(choices = c("right","left","top","bottom","none")),
  expression_unit  = list(choices = c("log2_tpm_plus_one","tpm")),
  show_stats       = list(type    = "logical"),
  show_n           = list(type    = "logical"),
  display_floor    = list(range   = c(1e-10, 1)),
  loq              = list(range   = c(1e-10, 1)),
  alpha            = list(range   = c(0, 1)),
  point_size       = list(range   = c(0.1, 20)),
  point_alpha      = list(range   = c(0, 1)),
  box_width        = list(range   = c(0.05, 1)),
  bar_width        = list(range   = c(0.05, 1)),
  discord_log10    = list(range   = c(0, 10))
)

# Levenshtein-like nudge for typos.
.opts_did_you_mean <- function(bad, allowed) {
  if (length(allowed) == 0) return(NULL)
  d <- adist(tolower(bad), tolower(allowed))[1, ]
  best <- allowed[which.min(d)]
  if (min(d) <= max(2, nchar(bad) %/% 2)) best else NULL
}

# Validate one key/value pair; return TRUE or stop()
.opts_validate <- function(k, v) {
  rule <- .opts_validators[[k]]
  if (is.null(rule)) return(TRUE)
  if (!is.null(rule$choices)) {
    if (!is.character(v) || length(v) != 1 || !(v %in% rule$choices)) {
      hint <- if (is.character(v) && length(v) == 1)
        .opts_did_you_mean(v, rule$choices) else NULL
      hint_str <- if (!is.null(hint)) sprintf("  Did you mean '%s'?", hint) else ""
      stop(sprintf("Invalid value for '%s'. Allowed: %s.%s",
                   k, paste(shQuote(rule$choices), collapse = ", "),
                   hint_str), call. = FALSE)
    }
  }
  if (!is.null(rule$range)) {
    if (!is.numeric(v) || length(v) != 1 ||
        v < rule$range[1] || v > rule$range[2])
      stop(sprintf("Invalid value for '%s'. Must be a single numeric in [%g, %g].",
                   k, rule$range[1], rule$range[2]), call. = FALSE)
  }
  if (identical(rule$type, "logical")) {
    if (!is.logical(v) || length(v) != 1 || is.na(v))
      stop(sprintf("Invalid value for '%s'. Must be TRUE or FALSE.", k),
           call. = FALSE)
  }
  TRUE
}

#' Get, set, print, reset, or bulk-load ctdnaTM options
#'
#' Single entry point for configuration. Works like base R `options()`.
#'
#' Modes (mutually exclusive):
#' \itemize{
#'   \item No args → print and invisibly return current config
#'   \item One or more unnamed character keys → get values
#'   \item Named arguments → set values
#'   \item \code{.reset = TRUE} → restore factory defaults
#'   \item \code{load = list(...)} → bulk-load a saved config (replaces
#'         the older \code{do.call(ctdna_opts, cfg)} idiom)
#' }
#'
#' @param ... Unnamed character key(s) to get, or named values to set.
#' @param load Optional named list of options to apply in one call
#'   (typically read back from \code{readRDS()}).
#' @param domain Optional character scalar restricting the printout to
#'   one family of options. Accepts: \code{"clinical"} (subject /
#'   cohort / recist / dose / mr / best_change / baseline /
#'   cancertype), \code{"ctdna"} (per-sample ctDNA columns and the
#'   methylation TF), \code{"genomic"} (gene + alteration columns),
#'   \code{"tmb"} (TMB and MSI), \code{"ihc"} (marker + score),
#'   \code{"expression"} (RNA-seq long-form), \code{"qc"} (sample-level
#'   QC), \code{"filter"} (default filter scheme + indication
#'   dictionary), \code{"rules"} (the \code{col_*} mappings used by the
#'   rule engine), \code{"stats"}, \code{"plot"}, \code{"recist_values"}
#'   (the CR/PR/SD/PD vocabulary), or \code{"*"} for the full dump.
#'   Default \code{NULL} prints every key (and emits a tip about this
#'   argument).
#' @param .reset If TRUE, reset to factory defaults.
#' @return Invisibly, the current config (or the requested value).
#' @examples
#' \dontrun{
#' # Configure for a study and save the config
#' ctdna_opts(tf = "TF_methyl_v3", maxvaf = "MAX_VAF")
#' saveRDS(ctdna_opts(), "study123_config.rds")
#'
#' # In a fresh session — load the config in one call
#' cfg <- readRDS("study123_config.rds")
#' ctdna_opts(load = cfg)
#' }
#' @export
ctdna_opts <- function(..., load = NULL, domain = NULL, .reset = FALSE) {
  if (isTRUE(.reset)) { .env$cfg <- .defaults; return(invisible(.env$cfg)) }
  if (!is.null(load)) {
    if (!is.list(load) || is.null(names(load)))
      stop("`load` must be a named list of option values.")
    bad <- setdiff(names(load), names(.defaults))
    if (length(bad))
      stop("Unknown key(s) in load: ", paste(bad, collapse = ", "))
    for (k in names(load)) .opts_validate(k, load[[k]])
    for (k in names(load)) .env$cfg[[k]] <- load[[k]]
    # v0.41.0: hint about domain if a "busy" set call
    if (length(load) >= 3L) .opts_hint()
    return(invisible(.env$cfg))
  }
  args <- list(...)
  if (length(args) == 0) {
    # v0.41.0: domain= argument filters the print to one prefix family
    .print_opts(domain = domain)
    if (is.null(domain)) .opts_hint()
    return(invisible(.env$cfg))
  }
  if (is.null(names(args)) || all(names(args) == "")) {
    keys <- unlist(args)
    bad <- setdiff(keys, names(.defaults))
    if (length(bad)) stop("Unknown key(s): ", paste(bad, collapse = ", "))
    if (length(keys) == 1) return(.env$cfg[[keys]])
    return(.env$cfg[keys])
  }
  bad <- setdiff(names(args), names(.defaults))
  if (length(bad))
    stop("Unknown key(s): ", paste(bad, collapse = ", "),
         "\nValid: ", paste(names(.defaults), collapse = ", "))
  # Validate every key/value before any mutation, so a bad value
  # rejects the whole call rather than half-applying.
  for (k in names(args)) .opts_validate(k, args[[k]])
  for (k in names(args)) .env$cfg[[k]] <- args[[k]]
  # v0.41.0: hint about domain if a "busy" set call
  if (length(args) >= 3L) .opts_hint()
  invisible(.env$cfg)
}

# v0.41.0: short hint pointing users at the domain= argument when the
# output / input is at risk of being overwhelming.
.opts_hint <- function() {
  message("Tip: ctdna_opts(domain = \"X\") shows or sets only one family. ",
          "Domains: clinical, ctdna, genomic, expression, ihc, tmb, qc, ",
          "filter, plot, stats, rules.")
}

.print_opts <- function(domain = NULL) {
  cfg <- .env$cfg
  fmt <- function(v) {
    if (is.null(v)) return("NULL")
    if (is.list(v))
      return(sprintf("[list of %d]", length(v)))
    if (length(v) > 1) paste0("[", paste(v, collapse = ", "), "]")
    else as.character(v)
  }
  mark <- function(k) if (!identical(cfg[[k]], .defaults[[k]])) " *" else ""
  allowed_str <- function(k) {
    rule <- .opts_validators[[k]]
    if (is.null(rule)) return("")
    if (!is.null(rule$choices))
      return(sprintf("  (choices: %s)",
                     paste(shQuote(rule$choices), collapse = ", ")))
    if (!is.null(rule$range))
      return(sprintf("  (range: [%g, %g])", rule$range[1], rule$range[2]))
    if (identical(rule$type, "logical"))
      return("  (TRUE / FALSE)")
    ""
  }

  # v0.41.0: each domain has a curated key list. Keys missing from cfg
  # are silently dropped so the printer never crashes on a stale spec.
  all_keys <- names(cfg)
  domain_specs <- list(
    clinical   = c("subject","cancertype","cohort","recist","dose","mr",
                    "best_change","baseline"),
    ctdna      = c("tf","maxvaf","meanvaf","cfdna_conc","n_somatic",
                    "sample","visit","time","ctdna_metric"),
    tumor      = c("lesion","lesion_pct"),
    genomic    = c("gene","alteration","variant_type","molecular_consequence",
                    "cnv_type","copy_number","vaf","mutation","panel"),
    tmb        = c("tmb","tmb_category","msi"),
    ihc        = c("marker","ihc_score","ihc_score_unit"),
    expression = c("expression","expression_unit"),
    qc         = grep("^sample_qc_", all_keys, value = TRUE),
    filter     = c("default_scheme","cancertype_dictionary",
                    "default_indication_schemes",
                    "rare_variant_filter","rare_variant_threshold",
                    "col_indication_adam"),
    rules      = grep("^col_", all_keys, value = TRUE),
    stats      = c("pairwise_test","cor_method","alpha","loq","discord_log10",
                    "display_floor"),
    plot       = c("point_size","point_alpha","box_width","bar_width",
                    "show_stats","show_n","median_color","clearance_color",
                    "stat_position","legend_position"),
    recist_values = c("recist_crpr","recist_ucrpr","recist_sd","recist_pdnena")
  )
  # Filter each group to only keys actually in cfg
  domain_specs <- lapply(domain_specs, function(ks) ks[ks %in% all_keys])
  used <- unique(unlist(domain_specs))
  leftover <- setdiff(all_keys, used)
  if (length(leftover) > 0L) domain_specs[["other"]] <- leftover

  # v0.41.2: capture the full list of domain names BEFORE we filter, so
  # the printer can surface "other domains: …" at the bottom even when
  # the caller restricts the output to a single domain.
  all_domain_names <- names(domain_specs)

  # v0.41.0: domain filter
  if (!is.null(domain)) {
    if (identical(domain, "*")) {
      # show everything; same as NULL
    } else {
      if (!domain %in% names(domain_specs))
        stop("Unknown domain '", domain, "'. Available: ",
             paste(names(domain_specs), collapse = ", "), call. = FALSE)
      domain_specs <- domain_specs[domain]
    }
  }

  cat("ctdnaTM options (* = modified)")
  if (!is.null(domain) && !identical(domain, "*"))
    cat("  [domain: ", domain, "]", sep = "")
  cat("\n")
  for (g in names(domain_specs)) {
    keys <- domain_specs[[g]]
    if (length(keys) == 0L) next
    cat(sprintf("\n  %s\n", g))
    for (k in keys)
      cat(sprintf("    %-26s : %s%s%s\n",
                   k, fmt(cfg[[k]]), mark(k), allowed_str(k)))
  }
  # v0.41.2: when a domain filter is active, show what other domains
  # exist so the user can navigate without consulting docs.
  if (!is.null(domain) && !identical(domain, "*")) {
    other_domains <- setdiff(all_domain_names, domain)
    if (length(other_domains) > 0L)
      cat(sprintf("\n  other domains: %s\n",
                  paste(other_domains, collapse = ", ")))
  } else {
    # Full dump: tell the user the domains exist
    cat(sprintf("\n  domains available: %s\n",
                paste(all_domain_names, collapse = ", ")))
  }
  cat("\n  set:   ctdna_opts(tf = 'TF_v3')   get: ctdna_opts('tf')   reset: ctdna_opts(.reset = TRUE)\n")
}

#' Show allowed values for restricted ctdna_opts keys
#'
#' Companion to \code{\link{ctdna_opts}}: prints (and invisibly returns)
#' the allowed values or numeric range for every restricted key. Use it
#' when you can't remember which strings a setting accepts.
#'
#' \code{ctdna_opts_choices()} with no argument prints the full table.
#' \code{ctdna_opts_choices("stat_position")} returns the allowed
#' values for that one key as a character vector (or a numeric range
#' for numeric keys).
#'
#' @param key Optional key name. If omitted, prints the full table.
#' @return Invisibly, the allowed values (or NULL for unrestricted keys).
#' @examples
#' ctdna_opts_choices()
#' ctdna_opts_choices("stat_position")
#' @export
ctdna_opts_choices <- function(key = NULL) {
  if (is.null(key)) {
    cat("Restricted ctdna_opts() keys:\n")
    cat(strrep("-", 70), "\n", sep = "")
    for (k in names(.opts_validators)) {
      rule <- .opts_validators[[k]]
      if (!is.null(rule$choices))
        cat(sprintf("  %-18s : %s\n", k,
                    paste(shQuote(rule$choices), collapse = ", ")))
      else if (!is.null(rule$range))
        cat(sprintf("  %-18s : numeric in [%g, %g]\n", k,
                    rule$range[1], rule$range[2]))
      else if (identical(rule$type, "logical"))
        cat(sprintf("  %-18s : TRUE / FALSE\n", k))
    }
    cat(strrep("-", 70), "\n", sep = "")
    cat("Unrestricted keys (column names, palette colors, etc.) accept any value.\n")
    return(invisible(.opts_validators))
  }
  if (!key %in% names(.opts_validators)) {
    message(sprintf("'%s' is not a restricted key. It accepts any value.", key))
    return(invisible(NULL))
  }
  rule <- .opts_validators[[key]]
  out <- rule$choices %||% rule$range %||% NA
  if (identical(rule$type, "logical")) out <- c(TRUE, FALSE)
  out
}

# Internal shortcut
.o <- function(k) .env$cfg[[k]]


# ---- Utilities --------------------------------------------------------------

# Null-coalescing helper (avoid taking on rlang dep)
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Pin tiny ctDNA values for log-scale display
#'
#' Replaces any value smaller than \code{at} with \code{at}. Used before
#' a log scale so that undetectable / sub-floor measurements still appear
#' on the plot instead of becoming \code{-Inf}. Default is
#' \code{ctdna_opts("display_floor")} (0.001 by default).
#' @param x Numeric vector (ctDNA TF or MR values, in fraction units).
#' @param at Floor value (default \code{ctdna_opts("display_floor")}).
#' @return Numeric vector the same length as x, with values \eqn{\geq} at.
#' @examples
#' # Default: clamp to ctdna_opts("display_floor") (0.001)
#' ctdna_floor_tf(c(0, 1e-5, 0.001, 0.05, 0.5))
#'
#' # Custom floor
#' ctdna_floor_tf(c(0, 0.0001, 0.01), at = 0.005)
#' @export
ctdna_floor_tf <- function(x, at = .o("display_floor")) pmax(as.numeric(x), at)

#' Stratify RECIST into 2/3/4 groups
#'
#' Mapping of raw labels to canonical groups is taken from ctdna_opts()
#' so non-standard labels can be remapped via the config.
#'
#' @param x Character vector of raw RECIST values.
#' @param scheme One of `"two"`, `"three"`, `"four"`.
#' Stratify RECIST values into one of six grouping schemes
#'
#' v0.38.0 expands the grouping options from 3 to 6. The argument name
#' has also changed from \code{scheme} to \code{recist_grouping}; the
#' old name is kept as a deprecated alias for backward compat.
#'
#' Grouping options (and their level orderings):
#' \itemize{
#'   \item \code{"two"}      -- CR+PR -> \strong{R}; SD+PD+NE+NA -> \strong{NR}
#'   \item \code{"three"}    -- CR+PR / SD / PD+NE+NA  (default)
#'   \item \code{"four"}     -- CR+PR / SD / PD / NE+NA
#'   \item \code{"four_alt"} -- CR / PR / SD / PD+NE+NA  (splits responders)
#'   \item \code{"five"}     -- CR / PR / SD / PD / NE+NA
#'   \item \code{"six"}      -- CR / PR / SD / PD / NE / NA
#' }
#'
#' @param x Character/factor vector of RECIST values.
#' @param recist_grouping One of \code{"two"}, \code{"three"} (default),
#'   \code{"four"}, \code{"four_alt"}, \code{"five"}, \code{"six"}.
#' @param scheme \strong{Deprecated.} Old name for \code{recist_grouping}.
#'   Still accepted for backward compatibility; emits a one-time warning.
#' @return Ordered factor.
#' @examples
#' x <- c("CR","PR","uCR","uPR","SD","PD","NE","NA")
#' ctdna_stratify_recist(x, "two")        # R / NR
#' ctdna_stratify_recist(x, "three")      # CR+PR / SD / PD+NE+NA
#' ctdna_stratify_recist(x, "four")       # CR+PR / SD / PD / NE+NA
#' ctdna_stratify_recist(x, "six")        # CR / PR / SD / PD / NE / NA
#' @export
ctdna_stratify_recist <- function(x,
                                    recist_grouping = c("three","two","four",
                                                          "four_alt","five","six"),
                                    scheme = NULL) {
  # Backward compat: scheme is deprecated alias for recist_grouping
  if (!is.null(scheme)) {
    .ctdna_warn_once("recist_grouping_scheme",
      "ctdna_stratify_recist: argument `scheme` is DEPRECATED; use `recist_grouping`.")
    recist_grouping <- scheme
  }
  recist_grouping <- match.arg(recist_grouping)
  x <- as.character(x); x[is.na(x)] <- "NA"
  is_cr_only <- toupper(trimws(x)) == "CR"
  is_pr_only <- toupper(trimws(x)) == "PR"
  is_crpr    <- x %in% .o("recist_crpr") | is_cr_only | is_pr_only
  is_sd      <- x %in% .o("recist_sd")
  is_pd      <- toupper(trimws(x)) == "PD"
  is_ne      <- toupper(trimws(x)) == "NE"
  is_na      <- toupper(trimws(x)) %in% c("NA","")

  out <- switch(recist_grouping,
    two = ifelse(is_crpr, "R", "NR"),
    three = ifelse(is_crpr, "CR+PR",
                    ifelse(is_sd, "SD", "PD+NE+NA")),
    four = ifelse(is_crpr, "CR+PR",
                    ifelse(is_sd, "SD",
                            ifelse(is_pd, "PD", "NE+NA"))),
    four_alt = ifelse(is_cr_only, "CR",
                       ifelse(is_pr_only, "PR",
                               ifelse(is_sd, "SD", "PD+NE+NA"))),
    five = ifelse(is_cr_only, "CR",
                   ifelse(is_pr_only, "PR",
                           ifelse(is_sd, "SD",
                                   ifelse(is_pd, "PD", "NE+NA")))),
    six = ifelse(is_cr_only, "CR",
                  ifelse(is_pr_only, "PR",
                          ifelse(is_sd, "SD",
                                  ifelse(is_pd, "PD",
                                          ifelse(is_ne, "NE", "NA")))))
  )
  lvls <- switch(recist_grouping,
    two       = c("R","NR"),
    three     = c("CR+PR","SD","PD+NE+NA"),
    four      = c("CR+PR","SD","PD","NE+NA"),
    four_alt  = c("CR","PR","SD","PD+NE+NA"),
    five      = c("CR","PR","SD","PD","NE+NA"),
    six       = c("CR","PR","SD","PD","NE","NA"))
  factor(out, levels = lvls)
}

#' On-treatment / baseline ratio with LOQ exclusion
#'
#' For every subject, computes `tf / baseline_tf` on each on-treatment
#' visit, where `baseline_tf` is the value at the baseline visit
#' (identified via `ctdna_opts("baseline")`). Subjects whose
#' baseline value is below the limit of quantification
#' (`ctdna_opts("loq")`) are excluded from the output — their ratios
#' are meaningless. Returns a long frame with the original rows that
#' survived plus two added columns (`baseline_tf`, `ratio`); baseline
#' rows themselves are dropped (ratio is undefined there).
#'
#' This is the building block for `ctdna_compute_mr()` and the family
#' of "fold change from baseline" plots.
#'
#' @param df Long data frame.
#' @param tf_col Column name holding the tumor-fraction (or analogous)
#'   per-sample value to ratio. Default \code{NULL}, which resolves
#'   via \code{ctdna_opts("tf")}.
#' @return Filtered data with `baseline_tf` and `ratio` added.
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 10, seed = 1)
#' prep <- ctdna_prepare(sim, verbose = FALSE)
#' r    <- ctdna_ratio(prep$samples)
#' head(r[, c("Patient_ID","Visit_name","methylTF","baseline_tf","ratio")])
#' @export
ctdna_ratio <- function(df, tf_col = NULL) {
  s <- .o("subject"); t <- .o("time")
  if ("Record_type" %in% names(df))
    df <- df[!(as.character(df$Record_type) %in% "tumor"), , drop = FALSE]
  tf <- tf_col %||% .o("tf")
  bl <- .o("baseline"); loq <- .o("loq")
  stopifnot(all(c(s, t, tf) %in% names(df)))
  base <- df[df[[t]] == bl, c(s, tf), drop = FALSE]
  names(base)[2] <- "baseline_tf"
  out <- merge(df, base, by = s, all.x = TRUE)
  out <- out[!is.na(out$baseline_tf) & out$baseline_tf > 0, , drop = FALSE]
  drop_row <- out$baseline_tf < loq & out[[tf]] < loq & out[[t]] != bl
  out <- out[!drop_row, , drop = FALSE]
  out$ratio <- out[[tf]] / out$baseline_tf
  out
}

#' Per-subject ctDNA summary for multi-modal integration
#'
#' Collapses a longitudinal ctDNA dataset into one row per subject with
#' canonical metrics (baseline TF, best ratio, last TF) and metadata.
#' Bridge between ctDNA and other modalities.
#'
#' @param df Longitudinal ctDNA data frame.
#' @return Data frame: one row per subject.
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 30, seed = 1)
#' prep <- ctdna_prepare(sim, verbose = FALSE)
#' ctdna_summary(prep$samples)
#' @export
ctdna_summary <- function(df) {
  s <- .o("subject"); t <- .o("time"); tf <- .o("tf"); bl <- .o("baseline")
  if ("Record_type" %in% names(df))
    df <- df[!(as.character(df$Record_type) %in% "tumor"), , drop = FALSE]
  ratio_df <- ctdna_ratio(df)
  rdf <- ratio_df[ratio_df[[t]] != bl, , drop = FALSE]
  best <- stats::aggregate(ratio ~ get(s), data = rdf, FUN = min)
  names(best) <- c(s, "best_ratio")
  last_tp <- utils::tail(sort(unique(df[[t]][df[[t]] != bl])), 1)
  last_df <- df[df[[t]] == last_tp, c(s, tf), drop = FALSE]
  names(last_df)[2] <- "last_tf"
  base_df <- df[df[[t]] == bl, c(s, tf), drop = FALSE]
  names(base_df)[2] <- "baseline_tf"
  meta_cols <- intersect(c(.o("recist"), .o("dose"), .o("cohort"),
                            .o("mr"), .o("best_change")), names(df))
  meta <- df[!duplicated(as.character(df[[s]])), c(s, meta_cols), drop = FALSE]
  # v0.50.0/v0.51.0/v0.53.0: best lesion-size % change. Prefer an already-present
  # (broadcast) best_lesion_change column, else a direct per-visit %change column
  # (ADTR PCHG), else derive from the per-visit lesion-size (SLD) column.
  les <- NULL; L <- .o("lesion"); LP <- .o("lesion_pct")
  if ("best_lesion_change" %in% names(df)) {
    les <- unique(df[, c(s, "best_lesion_change"), drop = FALSE])
    les <- les[!duplicated(les[[s]]), , drop = FALSE]
  } else if (LP %in% names(df)) {
    lm <- df[df[[t]] != bl, c(s, LP), drop = FALSE]
    lm <- lm[is.finite(lm[[LP]]), , drop = FALSE]
    if (nrow(lm)) {
      les <- stats::aggregate(lm[[LP]], by = list(lm[[s]]), FUN = min)
      names(les) <- c(s, "best_lesion_change")
    }
  } else if (L %in% names(df)) {
    lb <- df[df[[t]] == bl, c(s, L), drop = FALSE]; names(lb)[2] <- ".lb"
    lm <- merge(df[df[[t]] != bl, c(s, t, L), drop = FALSE], lb, by = s, all.x = TRUE)
    lm$.lp <- (lm[[L]] - lm$.lb) / lm$.lb * 100
    lm <- lm[is.finite(lm$.lp), , drop = FALSE]
    if (nrow(lm)) {
      les <- stats::aggregate(.lp ~ get(s), data = lm, FUN = min)
      names(les) <- c(s, "best_lesion_change")
    }
  }
  parts <- list(base_df, best, last_df, meta)
  if (!is.null(les)) parts <- c(parts, list(les))
  Reduce(function(a, b) merge(a, b, by = s, all = TRUE), parts)
}


#' List or register named measurement methods
#'
#' Several deliverables expose a \code{method} argument that selects
#' which underlying measurement to use. Methods are grouped into three
#' DOMAINS:
#'
#' \itemize{
#'   \item \code{"tf"} - ctDNA tumor-fraction signal used by D1, D2,
#'     D3, D5, D7. Choose between methylTF, maxVAF, meanVAF, or a
#'     custom column.
#'   \item \code{"mr"} - molecular-response classifier rule used by
#'     \code{ctdna_compute_mr}. Choose between ratio_drop (TF/baseline
#'     threshold) and absolute_floor (raw TF threshold).
#'   \item \code{"tumor"} - tumor-burden signal used by D7 for the
#'     y-axis. Choose between %-change in target lesions and SLD-based.
#' }
#'
#' By default \code{ctdna_methods()} prints ALL three domains so you
#' can see what's available for any function. Pass \code{domain = "tf"}
#' (etc.) to filter.
#'
#' Built-in TF methods:
#' \itemize{
#'   \item \code{"methyl"}   - methylation-based TF (default; \code{methylTF})
#'   \item \code{"maxvaf"}    - max VAF as TF proxy (\code{maxVAF})
#'   \item \code{"meanvaf"} - mean VAF (\code{meanVAF})
#'   \item \code{"custom"}   - escape hatch; user supplies \code{value_col}
#' }
#'
#' Built-in MR methods:
#' \itemize{
#'   \item \code{"ratio_drop"} (default) - ratio TF/baseline <= threshold
#'   \item \code{"absolute_floor"}       - TF <= absolute value
#' }
#'
#' Built-in tumor methods:
#' \itemize{
#'   \item \code{"recist_pct"} (default) - % change in target lesion SLD
#'   \item \code{"sld_pct"}              - generic SLD %-change
#' }
#'
#' @param domain \code{"all"} (default, shows all three), \code{"tf"},
#'   \code{"mr"}, or \code{"tumor"}.
#' @return Data frame with columns: \code{domain}, \code{method},
#'   \code{value_col}, \code{label}, \code{unit}, \code{floor},
#'   \code{used_by} (which package functions accept this domain).
#' @examples
#' ctdna_methods()             # all three domains
#' ctdna_methods("mr")         # just MR methods
#' @export
ctdna_methods <- function(domain = c("all","tf","mr","tumor")) {
  domain <- match.arg(domain)
  if (is.null(.env[["method_registry"]]))
    .env[["method_registry"]] <- .ctdna_method_registry()
  used_by <- list(
    tf    = "D1, D2, D3, D5, D7",
    mr    = "ctdna_compute_mr",
    tumor = "D7"
  )
  domains <- if (domain == "all") c("tf","mr","tumor") else domain
  out <- do.call(rbind, lapply(domains, function(d) {
    reg <- .env[["method_registry"]][[d]]
    do.call(rbind, lapply(names(reg), function(nm) {
      e <- reg[[nm]]
      data.frame(
        domain       = d,
        method       = nm,
        value_col    = e$value_col %||% NA_character_,
        label        = e$label %||% NA_character_,
        unit         = e$unit %||% NA_character_,
        display_unit = e$display_unit %||% NA_character_,
        floor        = e$floor %||% NA_real_,
        used_by      = used_by[[d]],
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(out) <- NULL
  out
}

#' Register a study-specific measurement method
#'
#' Adds a method to the in-session registry. Becomes usable everywhere
#' the relevant deliverable accepts a \code{method} argument.
#'
#' @param name Method name (single non-empty string).
#' @param domain \code{"tf"} / \code{"mr"} / \code{"tumor"}.
#' @param value_col Column in the data frame to read.
#' @param label Axis label.
#' @param unit Storage unit: \code{"fraction"} (default) or \code{"percent"}.
#' @param display_unit How to render the method on plot axes:
#'   \code{"percent"} (axis ticks formatted as percentages, e.g. 0.1\%),
#'   \code{"fraction"} (axis ticks as-is, e.g. 0.001), or
#'   \code{"actual"} (no log floor, no unit suffix; for count or
#'   concentration variables like TMB or cfDNA).
#' @param floor Display floor for log-scale plots.
#' @return Invisibly, the updated registry.
#' @examples
#' ctdna_register_method("my_tf", domain = "tf",
#'                         value_col = "my_special_TF",
#'                         label = "Special TF", unit = "fraction",
#'                         display_unit = "percent", floor = 1e-4)
#' @export
ctdna_register_method <- function(name, domain = c("tf","mr","tumor"),
                                    value_col = NULL,
                                    label = NULL,
                                    unit = c("fraction","percent"),
                                    display_unit = c("percent","fraction","actual"),
                                    floor = NULL) {
  domain <- match.arg(domain)
  unit <- match.arg(unit)
  display_unit <- match.arg(display_unit)
  if (!is.character(name) || length(name) != 1 || !nzchar(name))
    stop("`name` must be a single non-empty string.")
  if (is.null(.env[["method_registry"]]))
    .env[["method_registry"]] <- .ctdna_method_registry()
  .env[["method_registry"]][[domain]][[name]] <- list(
    value_col = value_col, label = label,
    unit = unit, floor = floor,
    display_unit = display_unit)
  message(sprintf("Registered method '%s' (domain = %s, display = %s).",
                   name, domain, display_unit))
  invisible(.env[["method_registry"]])
}

# Built-in method registry. Three domains: tf, mr, tumor.
.ctdna_method_registry <- function() list(
  tf = list(
    methyl  = list(value_col = .o("tf"),
                    label = "Methylation TF",
                    unit = "percent", floor = .o("display_floor"),
                    display_unit = "percent"),
    maxvaf  = list(value_col = .o("maxvaf"),
                    label = "Max VAF",
                    unit = "percent", floor = .o("display_floor"),
                    display_unit = "percent"),
    meanvaf = list(value_col = .o("meanvaf"),
                    label = "Mean VAF",
                    unit = "percent", floor = .o("display_floor"),
                    display_unit = "percent"),
    custom  = list(value_col = NULL, label = NULL,
                    unit = "percent", floor = .o("display_floor"),
                    display_unit = "percent")
  ),
  mr = list(
    ratio_drop     = list(rule = "ratio_le", default_threshold = 0.5,
                           label = "Ratio drop"),
    absolute_floor = list(rule = "value_le", default_threshold = 0.001,
                           label = "Absolute TF floor"),
    custom         = list(rule = NULL, default_threshold = NA,
                           label = "Custom rule")
  ),
  tumor = list(
    recist_pct = list(value_col = "best_pct_change",
                       label = "% change in target lesion SLD",
                       unit = "percent", floor = NA),
    sld_pct    = list(value_col = "sld_pct_change",
                       label = "% change in SLD",
                       unit = "percent", floor = NA),
    custom     = list(value_col = NULL, label = NULL,
                       unit = "percent", floor = NA)
  )
)

# Resolve a method name to its lookup entry
.resolve_method <- function(method, domain, value_col = NULL, label = NULL) {
  if (is.null(.env[["method_registry"]]))
    .env[["method_registry"]] <- .ctdna_method_registry()
  reg <- .env[["method_registry"]][[domain]]
  if (!method %in% names(reg))
    stop(sprintf("Unknown method '%s' for domain '%s'. ",
                  method, domain),
         "Available: ", paste(names(reg), collapse = ", "),
         ". See ctdna_methods().")
  e <- reg[[method]]
  if (method == "custom") {
    if (is.null(value_col))
      stop(sprintf("method = 'custom' (domain '%s') requires `value_col`.",
                    domain))
    e$value_col <- value_col
    if (!is.null(label)) e$label <- label
    if (is.null(e$label)) e$label <- value_col
  }
  e
}

# ----------------------------------------------------------------------------
# .resolve_methods(): plural / vector form. Accepts:
#   - a single string from the registry (e.g. "methyl")
#   - the literal "all" -> expands to c("methyl","maxvaf","meanvaf") in that
#     canonical order
#   - a length>=2 vector of registry names -> preserves user order, deduped
# Returns a character vector of resolved method keys, length >= 1.
#
# Errors on:
#   - mixing "all" with other names in a vector (ambiguous)
#   - "custom" in a multi-element vector (custom needs value_col + label;
#     not supported in vector mode)
#   - unknown names
# ----------------------------------------------------------------------------
.resolve_methods <- function(method, domain = "tf") {
  if (is.null(method) || length(method) == 0)
    stop("`method` cannot be NULL or empty.")
  m_in <- as.character(method)

  if (is.null(.env[["method_registry"]]))
    .env[["method_registry"]] <- .ctdna_method_registry()
  reg <- .env[["method_registry"]][[domain]]
  if (is.null(reg))
    stop(sprintf("No method registry for domain '%s'.", domain),
         call. = FALSE)

  # v0.41.2: accept three input styles per element, in order:
  #   1. exact registry name (the legacy: "methyl", "maxvaf", "meanvaf", "custom")
  #      — also "all" (a meta-name handled below)
  #   2. case-insensitive variant ("Methyl", "MaxVAF", "MeanVAF", "All")
  #   3. an actual column name that the method would produce
  #      ("methylTF" -> "methyl", "maxVAF" -> "maxvaf", "meanVAF" -> "meanvaf").
  #      Resolved via ctdna_opts() so non-default vendor names also work.
  # The error message lists ALL three forms when something fails to resolve,
  # so users see the actual columns they could also have typed.

  # Build column -> method back-map for the current opts settings
  col_to_method <- list()
  if (domain == "tf") {
    col_to_method[[as.character(.o("tf"))]]      <- "methyl"
    col_to_method[[as.character(.o("maxvaf"))]]  <- "maxvaf"
    col_to_method[[as.character(.o("meanvaf"))]] <- "meanvaf"
  }
  # Plus a lowercase index of registry names for case-insensitive match
  reg_names <- names(reg)
  reg_lc    <- setNames(reg_names, tolower(reg_names))
  meta_lc   <- c(all = "all")

  resolve_one <- function(x) {
    if (x %in% reg_names)                  return(x)
    if (identical(x, "all"))               return("all")
    lc <- tolower(x)
    if (lc %in% names(reg_lc))             return(reg_lc[[lc]])
    if (lc %in% names(meta_lc))            return(meta_lc[[lc]])
    # Column-name back-map (case-sensitive on column names, since vendor
    # casing is meaningful)
    if (x %in% names(col_to_method))       return(col_to_method[[x]])
    NA_character_
  }
  resolved <- vapply(m_in, resolve_one, character(1), USE.NAMES = FALSE)
  bad_idx <- is.na(resolved)
  if (any(bad_idx)) {
    bad <- m_in[bad_idx]
    avail_cols <- if (domain == "tf")
      paste(unique(names(col_to_method)), collapse = ", ") else ""
    stop(sprintf(
      "Unknown method(s) for domain '%s': %s.\n  Available method names: %s\n%s",
      domain, paste(bad, collapse = ", "),
      paste(reg_names, collapse = ", "),
      if (nzchar(avail_cols))
        sprintf("  Or pass an actual column name: %s\n", avail_cols)
      else ""),
      call. = FALSE)
  }
  m <- resolved

  if ("all" %in% m) {
    if (length(m) > 1L)
      stop("`method = \"all\"` cannot be mixed with other method names. ",
            "Pass `\"all\"` alone, or list specific methods.",
           call. = FALSE)
    return(c("methyl","maxvaf","meanvaf"))
  }
  if ("custom" %in% m && length(m) > 1L)
    stop("`method = \"custom\"` cannot be combined with other methods in a vector. ",
          "Use it alone with `value_col`.", call. = FALSE)
  unique(m)  # preserve user-specified order
}


#' Compute molecular response (MR) from longitudinal ctDNA data
#'
#' Derives a per-subject Response / NonResponse call from a longitudinal
#' ctDNA frame, by comparing the on-treatment tumor fraction (TF) to
#' baseline. Two modes are available:
#' \itemize{
#'   \item \code{at = "best"} (default): per subject, take the BEST
#'     (deepest) on-treatment / baseline ratio. The subject is called
#'     \code{"Response"} if that best ratio is \eqn{\le} \code{threshold},
#'     else \code{"NonResponse"}.
#'   \item \code{at = "<time_point>"} (e.g. \code{"C5D1"}): use the TF
#'     at that specific time-point. The subject is called
#'     \code{"Response"} if \eqn{TF_t / TF_{baseline} \le} \code{threshold}.
#' }
#'
#' The threshold is a fold-change ratio (not a percent). The default
#' \code{0.5} corresponds to "TF dropped by at least 50%". Subjects
#' missing the requested timepoint return \code{NA}.
#'
#' @param df Longitudinal ctDNA data with canonical columns
#'   (subject, time, tf — names per \code{ctdna_opts()}).
#' @param at \code{"best"} (default) or a string matching one of the
#'   on-treatment time-point labels in the data (e.g. \code{"C5D1"}).
#' @param threshold Ratio (TF / baseline TF) at or below which the
#'   subject is called a molecular responder. Default 0.5.
#' @param method MR-calling method. Default \code{"ratio_drop"} (the
#'   built-in TF / baseline-TF ratio comparison). User-registered
#'   methods can be referenced by name; see
#'   \code{ctdna_register_mr_method()}.
#' @return Data frame: one row per subject with columns
#'   \code{subject}, \code{baseline_tf}, \code{ratio}, \code{MR}, plus
#'   the time-point used (for traceability).
#' @examples
#' sim <- ctdna_make_mock_study(n_patients = 30)
#' mm  <- ctdna_prepare(sim, verbose = FALSE)
#' ctdna_compute_mr(mm$ctdna)                       # best-ratio MR
#' ctdna_compute_mr(mm$ctdna, at = "C5D1")          # MR at a specific visit
#' ctdna_compute_mr(mm$ctdna, threshold = 0.3)      # stricter (70% drop required)
#' @export
ctdna_compute_mr <- function(df, at = "best", threshold = NULL,
                              method = "ratio_drop") {
  s_col  <- .o("subject"); t_col <- .o("time")
  tf_col <- .o("tf");      bl    <- .o("baseline")

  # Resolve method (built-in or user-registered)
  m <- .resolve_method(method, "mr")
  if (is.null(threshold))
    threshold <- m$default_threshold %||% 0.5

  base <- df[df[[t_col]] == bl, c(s_col, tf_col), drop = FALSE]
  names(base) <- c(s_col, "baseline_tf")

  if (identical(at, "best") || identical(at, "BEST")) {
    rdf <- ctdna_ratio(df, tf_col = tf_col)
    rdf <- rdf[rdf[[t_col]] != bl, , drop = FALSE]
    if (nrow(rdf) == 0) {
      out <- merge(base, data.frame(setNames(list(character(0), numeric(0)),
                   c(s_col, "ratio"))), by = s_col, all.x = TRUE)
      out$MR <- NA_character_
      out$time_used <- "best"
      out$method <- method
      return(out)
    }
    best <- stats::aggregate(ratio ~ get(s_col), data = rdf, FUN = min)
    names(best) <- c(s_col, "ratio")
    # Also keep TF at the best visit for "absolute_floor" rule
    rdf2 <- merge(rdf, best, by = c(s_col, "ratio"))
    rdf2 <- rdf2[!duplicated(rdf2[[s_col]]), c(s_col, tf_col), drop = FALSE]
    names(rdf2) <- c(s_col, "tf_at")
    out <- merge(base, best, by = s_col, all.x = TRUE)
    out <- merge(out, rdf2, by = s_col, all.x = TRUE)
    out$time_used <- "best"
  } else {
    if (!at %in% df[[t_col]])
      stop("Time-point '", at, "' not found in column '", t_col,
           "'. Available: ",
           paste(unique(df[[t_col]]), collapse = ", "))
    tp_df <- df[df[[t_col]] == at, c(s_col, tf_col), drop = FALSE]
    names(tp_df) <- c(s_col, "tf_at")
    out <- merge(base, tp_df, by = s_col, all.x = TRUE)
    out$ratio    <- out$tf_at / out$baseline_tf
    out$time_used <- at
  }
  # Apply the rule
  if (identical(m$rule, "value_le")) {
    out$MR <- ifelse(is.na(out$tf_at), NA_character_,
                     ifelse(out$tf_at <= threshold, "Response", "NonResponse"))
  } else {
    # Default: ratio rule
    out$MR <- ifelse(is.na(out$ratio), NA_character_,
                     ifelse(out$ratio <= threshold, "Response", "NonResponse"))
  }
  out$method <- method
  out[, c(s_col, "baseline_tf", "ratio", "MR", "time_used", "method")]
}


# v0.42.0: tolerant visit-label resolver. Real datasets vary in visit
# labeling (`"Cycle 1 Day 1"`, `"C1D1"`, `"C01D01"`, `"Cycle1Day1"`,
# etc). Strict matching forces users to inspect data first. This:
#   - exact match first (cheap, common case)
#   - then normalized: case-fold + "Cycle X Day Y" -> "CXDY", strip
#     leading zeros, drop whitespace/punctuation
#   - returns the un-normalized df_visits value (so the caller can
#     filter against the actual column values)
#   - errors with edit-distance suggestion on no match
#   - refuses to guess when normalization gives multiple candidates
.resolve_visit <- function(at, df_visits) {
  if (is.null(at) || !is.character(at) || length(at) != 1L)
    return(at)
  if (at %in% df_visits) return(at)
  df_visits <- unique(df_visits)
  normalize <- function(x) {
    x <- toupper(trimws(x))
    x <- gsub("CYCLE\\s*0*([0-9]+)\\s*DAY\\s*0*([0-9]+)", "C\\1D\\2", x)
    x <- gsub("C0+([0-9]+)D0+([0-9]+)", "C\\1D\\2", x)
    x <- gsub("[[:space:]_\\-]+", "", x)
    x
  }
  at_n <- normalize(at)
  df_n <- normalize(df_visits)
  hits_idx <- which(df_n == at_n)
  if (length(hits_idx) == 1L) return(df_visits[hits_idx])
  if (length(hits_idx) > 1L)
    stop(sprintf("Time-point '%s' is ambiguous: matches %s. Use the exact label.",
                  at, paste(shQuote(df_visits[hits_idx]), collapse = ", ")),
         call. = FALSE)
  d <- utils::adist(toupper(at), toupper(df_visits))[1, ]
  nearest <- df_visits[which.min(d)]
  hint <- if (min(d) <= max(3L, nchar(at) %/% 3L))
            sprintf(" Did you mean '%s'?", nearest) else ""
  stop(sprintf("Time-point '%s' not found. Available: %s.%s",
                at, paste(shQuote(df_visits), collapse = ", "), hint),
       call. = FALSE)
}


#' Per-subject ctDNA value at "best" or a specific time-point
#'
#' Computes a single per-subject value from a longitudinal ctDNA frame,
#' on any of three scales, at either the best on-treatment time-point
#' (deepest response) or a specific named visit. Used by the
#' cross-modality plots and \code{ctdna_plot_reduction()} to give the
#' user flexible control over what "ctDNA reduction" means.
#'
#' \strong{Metric definitions:}
#' \itemize{
#'   \item \code{"ratio"} = \eqn{TF_t / TF_{baseline}}. Values in [0, ~1].
#'         0.1 means TF dropped to 10\% of baseline. This is the
#'         classical "molecular reduction" ratio. A threshold of 0.5
#'         (default for MR classification) means "TF dropped by at
#'         least 50\%".
#'   \item \code{"tf"} = raw tumor fraction at the chosen visit
#'         (\eqn{TF_t}). Useful when the absolute value matters.
#'   \item \code{"pct_change"} = \eqn{100 \times (TF_t - TF_{baseline}) /
#'         TF_{baseline}}. Negative values = reduction. -90\% is a
#'         strong response.
#' }
#'
#' \strong{Time-point definitions:}
#' \itemize{
#'   \item \code{at = "best"} - per subject, take the on-treatment
#'         visit with the LOWEST ratio (the deepest molecular response).
#'   \item \code{at = "<visit>"} (e.g. \code{"C5D1"}) - use the value at
#'         that specific visit. Subjects missing the visit return NA.
#'   \item \code{at = "last"} - use the chronologically last on-treatment
#'         visit per subject.
#' }
#'
#' @param df Longitudinal ctDNA data with canonical columns
#'   (subject, time, tf - per current \code{ctdna_opts()}).
#' @param metric \code{"ratio"} (default) / \code{"tf"} / \code{"pct_change"}.
#' @param at \code{"best"} (default) / \code{"last"} / named visit.
#' @return Data frame: one row per subject, with columns
#'   \code{subject}, \code{baseline_tf}, \code{value} (the requested
#'   metric), and \code{time_used} (for traceability).
#' @examples
#' sim <- ctdna_make_mock_study(n_patients = 20, seed = 1)
#' mm  <- ctdna_prepare(sim, verbose = FALSE)
#'
#' # Default — best-ratio per subject
#' head(ctdna_metric_at(mm$ctdna))
#'
#' # TF at C5D1
#' head(ctdna_metric_at(mm$ctdna, metric = "tf", at = "C5D1"))
#'
#' # Percentage change at the last on-treatment visit
#' head(ctdna_metric_at(mm$ctdna, metric = "pct_change", at = "last"))
#' @export
ctdna_metric_at <- function(df,
                             metric = c("ratio","tf","pct_change"),
                             at     = "best",
                             tf_col = NULL) {
  metric <- match.arg(metric)
  s_col <- .o("subject"); t_col <- .o("time")
  if (is.null(tf_col)) tf_col <- .o("tf")
  bl <- .o("baseline")

  base <- df[df[[t_col]] == bl, c(s_col, tf_col), drop = FALSE]
  names(base) <- c(s_col, "baseline_tf")

  if (identical(at, "best") || identical(at, "BEST")) {
    rdf <- ctdna_ratio(df, tf_col = tf_col)
    rdf <- rdf[rdf[[t_col]] != bl, , drop = FALSE]
    if (nrow(rdf) == 0) {
      out <- merge(base,
                    data.frame(setNames(list(character(0), numeric(0),
                                              numeric(0), character(0)),
                                         c(s_col, "ratio", "tf_at",
                                           "time_used"))),
                    by = s_col, all.x = TRUE)
      out$value <- NA_real_
      return(out[, c(s_col, "baseline_tf", "value", "time_used")])
    }
    # Per subject, take the visit with the lowest ratio
    rdf <- rdf[order(rdf[[s_col]], rdf$ratio), ]
    rdf <- rdf[!duplicated(rdf[[s_col]]), ]
    rdf$time_used <- rdf[[t_col]]
    # Get the TF at that best visit
    keys <- paste(rdf[[s_col]], rdf$time_used)
    tf_at <- df[paste(df[[s_col]], df[[t_col]]) %in% keys,
                 c(s_col, t_col, tf_col), drop = FALSE]
    names(tf_at)[3] <- "tf_at"
    rdf <- merge(rdf[, c(s_col, "ratio", "time_used")],
                  tf_at, by.x = c(s_col, "time_used"),
                  by.y = c(s_col, t_col), all.x = TRUE)
    out <- merge(base, rdf, by = s_col, all.x = TRUE)
  } else if (identical(at, "last")) {
    on_tx <- df[df[[t_col]] != bl, , drop = FALSE]
    on_tx <- on_tx[order(on_tx[[s_col]], on_tx[[t_col]]), ]
    on_tx <- on_tx[!duplicated(on_tx[[s_col]], fromLast = TRUE), ]
    last_df <- on_tx[, c(s_col, t_col, tf_col), drop = FALSE]
    names(last_df) <- c(s_col, "time_used", "tf_at")
    out <- merge(base, last_df, by = s_col, all.x = TRUE)
    out$ratio <- out$tf_at / out$baseline_tf
  } else {
    # v0.42.0: tolerant visit-label resolver. "Cycle 1 Day 1" / "C01D01"
    # / "Cycle1Day1" all normalize to "C1D1". On no match, error message
    # includes a fuzzy "Did you mean ...?" suggestion.
    at <- .resolve_visit(at, unique(df[[t_col]]))
    tp_df <- df[df[[t_col]] == at, c(s_col, tf_col), drop = FALSE]
    names(tp_df) <- c(s_col, "tf_at")
    out <- merge(base, tp_df, by = s_col, all.x = TRUE)
    out$ratio <- out$tf_at / out$baseline_tf
    out$time_used <- at
  }

  out$value <- switch(metric,
    ratio      = out$ratio,
    tf         = out$tf_at,
    pct_change = 100 * (out$tf_at - out$baseline_tf) / out$baseline_tf)
  out[, c(s_col, "baseline_tf", "value", "time_used")]
}


# ---- Color scheme + theme ---------------------------------------------------

#' Unified RECIST color scheme
#'
#' Named vector mapping every label that
#' \code{\link{ctdna_stratify_recist}} can emit (across all six schemes:
#' `"two"`, `"three"`, `"four"`, `"four_alt"`, `"five"`, `"six"`) to a
#' hex color. Encodes a green->red response gradient — responders
#' green, stable disease amber, progression red — with non-gradient
#' greys for NE (light) and NA (dark). The returned vector is the
#' source of truth used by [ctdna_scale_recist()] and
#' [ctdna_scale_recist_fill()]; consume it directly if you build your
#' own plot.
#'
#' Legacy `/`-separator keys (e.g. `"CR/PR"`, `"PD/NE/NA"`) are kept
#' as aliases of their `+`-separator equivalents for back-compat with
#' pre-v0.41.4 user code.
#'
#' @return Named character vector of hex codes.
#' @examples
#' pal <- ctdna_colors()
#' pal[["CR"]]; pal[["PD"]]; pal[["NE"]]
#'
#' # Use directly in ggplot
#' # ggplot(df, aes(group, value, color = group)) +
#' #   geom_point() +
#' #   scale_color_manual(values = ctdna_colors())
#' @export
ctdna_colors <- function() c(
  # v0.41.5: green -> red response gradient; NE / NA are non-gradient
  # greys per user spec (NE = light, NA = dark). Keys cover every
  # ctdna_stratify_recist() output across all 6 schemes. Collapsed
  # labels take the color of the worse-response component, since the
  # box represents a group whose clinical worst case is what
  # interpreters react to.
  #
  # ---- responders (greens) ----
  "CR"          = "#1B7837",   # deep green
  "PR"          = "#5AAE61",   # medium green
  "CR+PR"       = "#5AAE61",   # collapsed responders (three/four/four_alt top)
  "R"           = "#5AAE61",   # two-scheme top
  # ---- unconfirmed responders (yellow-greens) ----
  "uCR"         = "#A6DBA0",   # light green
  "uPR"         = "#A6DBA0",
  "uCR+uPR"     = "#A6DBA0",
  # ---- stable disease (amber, gradient midpoint) ----
  "SD"          = "#F1C232",   # amber
  # ---- progression (red) ----
  "PD"          = "#D62728",   # red
  # ---- NE / NA (non-gradient greys per user spec) ----
  "NE"          = "#BDBDBD",   # light grey
  "NA"          = "#424242",   # dark grey
  # ---- collapsed labels containing NE / NA ----
  # Mixed-class collapsed labels: dominated by PD (clinical signal)
  # when PD is present, by the worse grey when only NE / NA present.
  "PD+NE+NA"    = "#D62728",   # PD dominates → red
  "NE+NA"       = "#757575",   # midpoint of light + dark grey
  "NR"          = "#757575",   # two-scheme bottom (SD+PD+NE+NA mix) → mid grey
  # ---- legacy "/" aliases for back-compat (pre-v0.41.4 external
  #      code) -- track their "+" counterparts ----
  "CR/PR"       = "#5AAE61",
  "uCR/PR"      = "#A6DBA0",
  "uCR/uPR"     = "#A6DBA0",
  "SD/PD/NE/NA" = "#757575",
  "PD/NE/NA"    = "#D62728"
)

#' Dose color palette
#'
#' A 5-color palette deliberately chosen so dose and RECIST never
#' visually collide (RECIST uses the green->red gradient from
#' [ctdna_colors()]; dose uses blue / gold / grey / red). Beyond 5
#' levels, the palette is interpolated via
#' `grDevices::colorRampPalette()`.
#'
#' @param n Number of levels to return. Default 3.
#' @return Character vector of `n` hex codes.
#' @examples
#' ctdna_dose_palette(3)
#' ctdna_dose_palette(5)
#'
#' # Beyond five levels, palette is interpolated
#' ctdna_dose_palette(8)
#' @export
ctdna_dose_palette <- function(n = 3) {
  if (n <= length(.dose_pal)) .dose_pal[seq_len(n)] else
    grDevices::colorRampPalette(.dose_pal)(n)
}

#' ggplot2 color scale for RECIST (outlines, points, lines)
#'
#' Wraps `ggplot2::scale_color_manual()` with the value vector from
#' [ctdna_colors()] and a `"grey70"` na.value. Use on layers where
#' RECIST is mapped to `color =` (outlines, points, lines, error
#' bars). For box fills and ribbon fills, use [ctdna_scale_recist_fill()].
#'
#' @param ... Passed to `ggplot2::scale_color_manual()` — typically
#'   `name = "RECIST"`, `breaks = ...`, or `guide = ...`.
#' @return A ggplot2 scale object.
#' @examples
#' # Pair with the deck's light-fill, strong-outline boxplot style
#' # ggplot(df, aes(dose, value, color = RECIST, fill = RECIST)) +
#' #   geom_boxplot(alpha = 0.3) +
#' #   ctdna_scale_recist() +
#' #   ctdna_scale_recist_fill()
#' @export
ctdna_scale_recist <- function(...) {
  ggplot2::scale_color_manual(values = ctdna_colors(), na.value = "grey70", ...)
}

#' ggplot2 fill scale for RECIST (box fills, ribbon fills)
#'
#' Uses the same hex values as [ctdna_scale_recist()]. Pair with
#' `geom_boxplot(alpha = 0.3, ...)` to get the deck's light-fill +
#' strong-outline boxplot style.
#'
#' @param ... Passed to ggplot2::scale_fill_manual.
#' @return A ggplot2 scale object.
#' @examples
#' # See ctdna_scale_recist() for a paired example.
#' ctdna_scale_recist_fill()
#' @export
ctdna_scale_recist_fill <- function(...) {
  ggplot2::scale_fill_manual(values = ctdna_colors(), na.value = "grey85", ...)
}

#' ggplot2 color scale for dose
#'
#' Wraps `ggplot2::scale_color_manual()` with values from
#' [ctdna_dose_palette()] (up to 8 levels). Use on layers where dose
#' is mapped to `color =`; for box fills use [ctdna_scale_dose_fill()].
#'
#' @param ... Passed to `ggplot2::scale_color_manual()`.
#' @return A ggplot2 scale object.
#' @examples
#' # Pair with the deck's light-fill, strong-outline boxplot style
#' # ggplot(df, aes(dose, value, color = dose, fill = dose)) +
#' #   geom_boxplot(alpha = 0.3) +
#' #   ctdna_scale_dose() +
#' #   ctdna_scale_dose_fill()
#' @export
ctdna_scale_dose <- function(...) {
  ggplot2::scale_color_manual(values = ctdna_dose_palette(8),
                              na.value = "grey70", ...)
}

#' ggplot2 fill scale for dose
#'
#' Same hexes as [ctdna_scale_dose()]. Pair with `geom_boxplot(alpha = 0.3)`
#' for the deck's light-fill + strong-outline style.
#'
#' @param ... Passed to ggplot2::scale_fill_manual.
#' @return A ggplot2 scale object.
#' @examples
#' # Paired with ctdna_scale_dose() — outline + fill share hexes
#' # ggplot(df, aes(dose, value, color = dose, fill = dose)) +
#' #   geom_boxplot(alpha = 0.3) +
#' #   ctdna_scale_dose() + ctdna_scale_dose_fill()
#' @export
ctdna_scale_dose_fill <- function(...) {
  ggplot2::scale_fill_manual(values = ctdna_dose_palette(8),
                             na.value = "grey85", ...)
}

#' Canonical shape scale for dose
#'
#' Use this when RECIST already occupies the colour aesthetic and dose
#' needs a second visual channel. Canonical shapes per dose tier:
#' \itemize{
#'   \item Low  = circle (16)
#'   \item Mid  = triangle (17)
#'   \item High = square (15)
#'   \item additional doses cycle through diamond, plus, cross, asterisk
#' }
#' Combined with [ctdna_scale_recist()] for RECIST colour, this gives a
#' clear two-channel encoding (colour = RECIST, shape = dose).
#'
#' @param ... Passed to ggplot2::scale_shape_manual.
#' @return A ggplot2 scale object.
#' @examples
#' # When dose has too many levels for color, use shape
#' # ggplot(df, aes(time, value, shape = dose)) +
#' #   geom_point() + ctdna_scale_dose_shape()
#' @export
ctdna_scale_dose_shape <- function(...) {
  values <- c(Low = 16, Mid = 17, High = 15,
              "Very Low" = 1, "Very High" = 18,
              # numeric ordering fall-throughs
              `1` = 16, `2` = 17, `3` = 15, `4` = 18,
              `5` = 1,  `6` = 2,  `7` = 0,  `8` = 5)
  ggplot2::scale_shape_manual(values = values, na.value = 4, ...)
}

#' ggplot2 theme for ctDNA deliverables
#'
#' Clean theme matching the source deck plots: white panel, light gray
#' strip backgrounds, subtle grid, thin panel border.
#'
#' @param base_size Base font size.
#' @examples
#' # Apply on its own
#' library(ggplot2)
#' ggplot(mtcars, aes(mpg, wt)) + geom_point() + ctdna_theme()
#'
#' # All ctdnaTM plots already use ctdna_theme() internally
#' @export
ctdna_theme <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.background  = ggplot2::element_rect(fill = "white", color = NA),
      panel.border      = ggplot2::element_rect(color = "grey50",
                                                fill  = NA,
                                                linewidth = 0.4),
      panel.grid.major  = ggplot2::element_line(color = "grey92",
                                                linewidth = 0.3),
      panel.grid.minor  = ggplot2::element_blank(),
      strip.background  = ggplot2::element_rect(fill = "grey85", color = NA),
      strip.text        = ggplot2::element_text(face = "bold",
                                                color = "grey15"),
      axis.title        = ggplot2::element_text(face = "bold"),
      axis.text         = ggplot2::element_text(color = "grey20"),
      plot.title        = ggplot2::element_text(face = "bold",
                                                size = base_size + 1),
      plot.subtitle     = ggplot2::element_text(face = "italic",
                                                color = "grey30",
                                                size = base_size - 1),
      plot.caption      = ggplot2::element_text(color = "grey40",
                                                size = base_size - 2),
      legend.background = ggplot2::element_blank(),
      legend.key        = ggplot2::element_blank()
    )
}


# ---- Example data -----------------------------------------------------------

# Internal helper: rename internal canonical names in a data frame to whatever
# the user has configured via ctdna_opts(). For example, after
#   ctdna_opts(tf = "TF_methyl_v3")
# .apply_opts_names(df, c(methylTF = "tf")) will rename df$methylTF -> df$TF_methyl_v3.
.apply_opts_names <- function(df, name_map) {
  # name_map: named character vector where names are current df column names
  # (canonical) and values are ctdna_opts keys whose live value should be used.
  current <- names(df)
  for (i in seq_along(current)) {
    if (current[i] %in% names(name_map)) {
      current[i] <- .o(name_map[[current[i]]])
    }
  }
  names(df) <- current
  df
}

# v0.22.1: harmonize common spellings of canonical TF-domain columns.
# Many vendor pipelines and user scripts use varying casing/separators
# for max VAF, mean VAF, and methyl TF. Before any column-name
# validation runs, walk df and rename any recognised variant to the
# canonical name. Emits one message per renamed column unless quiet.
#
# Canonical names: methylTF, maxVAF, meanVAF.
.ctdna_harmonize_tf_columns <- function(df, quiet = TRUE) {
  if (!is.data.frame(df) || ncol(df) == 0) return(df)
  variants <- list(
    methylTF = c("methylTF","methyltf","methyl_TF","methyl_tf","Methyl_TF",
                  "Methyl_tf","MethylTF","methylTumorFraction",
                  "methyl_tumor_fraction"),
    maxVAF   = c("maxVAF","maxvaf","maxAF","maxaf","max_AF","max_af",
                  "max.af","max.vaf","max_vaf","MaxVAF","Max_VAF","MAX_VAF",
                  "MAX_AF","maxAf"),
    meanVAF  = c("meanVAF","meanvaf","mean_vaf","mean_VAF","mean.vaf",
                  "mean.AF","MeanVAF","Mean_VAF","MEAN_VAF","meanAF","mean_AF")
  )
  current <- names(df)
  changed <- character()
  for (canon in names(variants)) {
    # All hits — including the canonical name itself if present
    hit <- intersect(current, variants[[canon]])
    if (length(hit) <= 1) {
      # 0 hits -> nothing to do; 1 hit -> may need rename if not canonical
      if (length(hit) == 1 && hit != canon) {
        current[current == hit] <- canon
        changed <- c(changed, sprintf("%s -> %s", hit, canon))
      }
      next
    }
    # 2+ hits = ambiguous: canonical + a variant, or two variants
    stop(sprintf(
      "Multiple aliases of canonical column '%s' present in df: %s. ",
      canon, paste(hit, collapse = ", ")),
      "Drop the duplicates and keep one (named '", canon, "').")
  }
  names(df) <- current
  if (length(changed) > 0 && !isTRUE(quiet))
    message(".ctdna_harmonize_tf_columns: renamed ",
             paste(changed, collapse = "; "))
  df
}



.alt_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit)) hit[1] else NULL
}

# Normalize a vector of patient IDs for robust matching:
#   - convert to character
#   - trim leading/trailing whitespace
#   - replace empty strings with NA
.norm_pid <- function(x) {
  out <- trimws(as.character(x))
  out[!nzchar(out)] <- NA_character_
  out
}

# Resolve top_annotation values for an ordered vector of patient IDs.
# Looks in patient_data first (one row per patient, joined by ID);
# falls back to df (alteration frame -- aggregates one value per
# patient via first non-NA, with a tie-breaker warning if multiple
# distinct values exist for a single patient).
#
# Returns a list with $values (named list: ann -> character vector
# aligned to pat_order) and $diag (character vector of one-line
# diagnostic messages).
.resolve_top_annotations <- function(top_annotations, pat_order,
                                       patient_data, df,
                                       patient_col_pd,    # column in pd, or NULL
                                       patient_col_df) {  # column in df
  diag <- character(0)
  values <- list()

  # Normalize patient orders for matching
  po_norm <- .norm_pid(pat_order)

  # Build pd lookup (if patient_data given)
  pd_norm_ids <- NULL
  if (!is.null(patient_data) && !is.null(patient_col_pd) &&
      patient_col_pd %in% names(patient_data)) {
    pd_norm_ids <- .norm_pid(patient_data[[patient_col_pd]])
  }

  # Build df lookup (collapse to one row per patient by first non-NA)
  df_lookup <- NULL
  if (!is.null(df) && !is.null(patient_col_df) &&
      patient_col_df %in% names(df)) {
    df_norm_ids <- .norm_pid(df[[patient_col_df]])
    df_lookup <- list(ids = df_norm_ids, df = df)
  }

  for (ann in top_annotations) {
    src <- NULL
    aligned <- rep(NA_character_, length(po_norm))

    # 1) Try patient_data first
    if (!is.null(pd_norm_ids) && ann %in% names(patient_data)) {
      idx <- match(po_norm, pd_norm_ids)
      aligned <- as.character(patient_data[[ann]])[idx]
      src <- "patient_data"
    } else if (!is.null(df_lookup) && ann %in% names(df_lookup$df)) {
      # 2) Fall back to df: for each patient, take first non-NA value
      vals <- as.character(df_lookup$df[[ann]])
      ids  <- df_lookup$ids
      multi_value_patients <- character(0)
      for (i in seq_along(po_norm)) {
        rows <- which(ids == po_norm[i] & !is.na(vals) & nzchar(vals))
        if (length(rows) == 0L) next
        uniq <- unique(vals[rows])
        aligned[i] <- uniq[1]
        if (length(uniq) > 1L)
          multi_value_patients <- c(multi_value_patients, po_norm[i])
      }
      if (length(multi_value_patients) > 0L)
        diag <- c(diag, sprintf(
          "annotation '%s' had multiple distinct values per patient for %d patient(s); used the first.",
          ann, length(multi_value_patients)))
      src <- "df"
    } else {
      diag <- c(diag, sprintf(
        "annotation '%s' not found in patient_data or df; skipped.", ann))
      next
    }

    # Diagnose unmatched patients (NA after join)
    n_total   <- length(po_norm)
    n_missing <- sum(is.na(aligned))
    if (n_missing > 0L) {
      # Try recovery for patient_data: case-insensitive + numeric-coerce
      if (src == "patient_data" && n_missing > 0L &&
          !is.null(pd_norm_ids)) {
        unmatched <- which(is.na(aligned))
        # Try case-insensitive match
        idx2 <- match(tolower(po_norm[unmatched]), tolower(pd_norm_ids))
        recovered <- !is.na(idx2)
        if (any(recovered)) {
          aligned[unmatched[recovered]] <-
            as.character(patient_data[[ann]])[idx2[recovered]]
          diag <- c(diag, sprintf(
            "annotation '%s': recovered %d patient(s) via case-insensitive match.",
            ann, sum(recovered)))
        }
      }
      n_missing <- sum(is.na(aligned))
      if (n_missing > 0L) {
        unmatched_ids <- pat_order[is.na(aligned)]
        preview <- paste(utils::head(unmatched_ids, 5), collapse = ", ")
        if (length(unmatched_ids) > 5L)
          preview <- paste0(preview, ", ...")
        diag <- c(diag, sprintf(
          "annotation '%s': %d of %d patients had no value in %s (rendered grey). IDs: %s",
          ann, n_missing, n_total, src, preview))
      }
    }
    values[[ann]] <- aligned
  }

  list(values = values, diag = diag)
}

# Apply a "keep these values" filter; FALSE disables; NULL disables; auto
# column resolution. Returns possibly-shrunk df + a 1-row diagnostic.
.apply_keep_filter <- function(df, col, values, label, verbose) {
  if (isFALSE(values) || is.null(values) || is.null(col)) {
    if (verbose && (isFALSE(values) || is.null(values)))
      message(sprintf("  %-22s : skipped", label))
    else if (verbose)
      message(sprintf("  %-22s : column not present, skipped", label))
    return(df)
  }
  n0 <- nrow(df)
  keep <- df[[col]] %in% values | is.na(df[[col]])
  df <- df[keep, , drop = FALSE]
  if (verbose)
    message(sprintf("  %-22s : %d -> %d rows (keep %s in '%s')",
                    label, n0, nrow(df),
                    paste(values, collapse = "|"), col))
  df
}

# Apply a numeric-threshold filter (max value allowed); FALSE disables.
.apply_max_filter <- function(df, col, max_value, label, verbose) {
  if (isFALSE(max_value) || is.null(max_value) || is.null(col)) {
    if (verbose && (isFALSE(max_value) || is.null(max_value)))
      message(sprintf("  %-22s : skipped", label))
    else if (verbose)
      message(sprintf("  %-22s : column not present, skipped", label))
    return(df)
  }
  n0 <- nrow(df)
  vals <- suppressWarnings(as.numeric(df[[col]]))
  keep <- is.na(vals) | vals <= max_value
  df <- df[keep, , drop = FALSE]
  if (verbose)
    message(sprintf("  %-22s : %d -> %d rows ('%s' <= %g)",
                    label, n0, nrow(df), col, max_value))
  df
}

# Apply a numeric-threshold filter (min value required); FALSE disables.
.apply_min_filter <- function(df, col, min_value, label, verbose) {
  if (isFALSE(min_value) || is.null(min_value) || is.null(col)) {
    if (verbose && (isFALSE(min_value) || is.null(min_value)))
      message(sprintf("  %-22s : skipped", label))
    else if (verbose)
      message(sprintf("  %-22s : column not present, skipped", label))
    return(df)
  }
  n0 <- nrow(df)
  vals <- suppressWarnings(as.numeric(df[[col]]))
  keep <- is.na(vals) | vals >= min_value
  df <- df[keep, , drop = FALSE]
  if (verbose)
    message(sprintf("  %-22s : %d -> %d rows ('%s' >= %g)",
                    label, n0, nrow(df), col, min_value))
  df
}




# ---- Column-by-column filter registry (Guardant Infinity) -------------------
#
# A registry-driven, per-column filter for genomic alterations. Each
# Guardant Infinity column has an entry of the form:
#
#   list(enabled = TRUE/FALSE,
#        op      = "%in%" / ">=" / "<=" / "==",
#        value   = allowlist or threshold,
#        desc    = "short description")
#
# Users get the default registry, modify it (turn columns on/off, change
# thresholds), and apply it to a data frame. Per-gene-set overrides are
# merged on top of the general registry so each set can have its own
# rules where needed.
# -----------------------------------------------------------------------------

# Default registry — every column in the Guardant Infinity report that
# is sensibly filterable. enabled = FALSE by default for columns that
# would over-filter in typical analyses; thresholds set to permissive
# defaults so turning them on doesn't remove everything.
.ctdna_filter_registry_defaults <- function() list(

  # ===========================================================================
  # v1.0.0 unified filter registry.
  #
  # Each entry has $tier = "basic" or "general".
  #   tier="basic"   -> applies to ALL surviving rows as the QC safety net
  #                     (step 6 of ctdna_variant_filter hierarchy). Pure QC.
  #   tier="general" -> cohort-tunable rules that apply to rows NOT claimed
  #                     by any user/built-in scheme (step 5 of hierarchy).
  #
  # ctdna_filter_opts() lets the user toggle entries on/off and change
  # thresholds. Basic entries CAN be tuned by the user but the defaults
  # are deliberately minimal (QC only).
  #
  # Values match canonical Guardant Infinity report values (e.g.
  # Sample_status="SUCCESS", Somatic_status="somatic"). For non-Guardant
  # data, users either rename their column values OR set
  # ctdna_opts(col_<name>="<their_col>") and update the entry value.
  # ===========================================================================

  # ---- BASIC tier (pure QC; applies to all surviving rows) ----
  Sample_status = list(enabled = TRUE, op = "%in%",
                        value = "SUCCESS",
                        tier = "basic",
                        desc = "QC: keep only samples with SUCCESS status."),
  ctDNA_detection_status = list(enabled = TRUE, op = "%in%",
                                  value = "Detected",
                                  tier = "basic",
                                  desc = "QC: keep only ctDNA-detected samples."),
  cfDNA_ng = list(enabled = TRUE, op = ">=", value = 5,
                   tier = "basic",
                   desc = "QC: minimum cfDNA input (ng)."),
  Plasma_ml_input = list(enabled = TRUE, op = ">=", value = 4,
                          tier = "basic",
                          desc = "QC: minimum plasma volume input (mL)."),

  # ---- GENERAL tier (cohort-tunable; applies to unclaimed rows) ----

  # --- Identity / context (informational; off by default) -------------------
  Patient_ID = list(enabled = FALSE, op = "%in%", value = character(0),
                     tier = "general",
                     desc = "Restrict to specific patient IDs."),
  Visit_name = list(enabled = FALSE, op = "%in%", value = character(0),
                     tier = "general",
                     desc = "Restrict to specific visit names (e.g. C1D1)."),
  Cancertype = list(enabled = FALSE, op = "%in%", value = character(0),
                     tier = "general",
                     desc = "Restrict to specific indications."),

  # --- Variant identity ----------------------------------------------------
  Variant_type = list(enabled = TRUE, op = "%in%",
                       value = c("SNV","Indel","CNV","Fusion","LGR"),
                       tier = "general",
                       desc = "Allowed variant types (Guardant Infinity values)."),
  Gene = list(enabled = FALSE, op = "%in%", value = character(0),
              tier = "general",
              desc = "Restrict to specific genes."),
  Chromosome = list(enabled = FALSE, op = "%in%", value = character(0),
                    tier = "general",
                    desc = "Restrict to specific chromosomes."),

  # --- Functional / consequence -------------------------------------------
  Functional_impact = list(enabled = FALSE, op = "%in%",
                            value = "deleterious",
                            tier = "general",
                            desc = paste("Guardant Functional_impact values:",
                                          "'deleterious', 'reversion_cis', 'reversion'.",
                                          "Off by default in v1.0.0 — enable when desired.")),
  Molecular_consequence = list(enabled = FALSE, op = "%in%",
                                value = c("missense","nonsense","frameshift",
                                          "splice_acceptor","splice_donor",
                                          "start_lost"),
                                tier = "general",
                                desc = "Guardant Molecular_consequence values."),
  Splice_effect = list(enabled = FALSE, op = "%in%", value = character(0),
                        tier = "general",
                        desc = "Splice effect category."),

  # --- Somatic / germline / clinical interpretation ------------------------
  Somatic_status = list(enabled = TRUE, op = "not_in",
                         value = "somatic_putative_ch",
                         tier = "general",
                         desc = "Drop CHIP (somatic_putative_ch) variants by default."),
  ClinVar = list(enabled = FALSE, op = "%in%",
                  value = c("Pathogenic","Likely_pathogenic",
                            "Pathogenic/Likely_pathogenic"),
                  tier = "general",
                  desc = "ClinVar clinical significance (Guardant labels)."),
  Mutant_allele_status = list(enabled = FALSE, op = "%in%",
                               value = "biallelic",
                               tier = "general",
                               desc = "Zygosity / biallelic status."),

  # --- Numeric thresholds -------------------------------------------------
  VAF_percentage = list(enabled = FALSE, op = ">=", value = 0.5,
                         tier = "general",
                         desc = "Minimum variant allele frequency (%)."),
  Genomic_max_VAF_percentage = list(enabled = FALSE, op = ">=", value = 0.5,
                                     tier = "general",
                                     desc = "Minimum max-VAF across genomic calls (%)."),
  Mol_count = list(enabled = FALSE, op = ">=", value = 50,
                    tier = "general",
                    desc = "Minimum supporting molecule count."),
  Copy_number = list(enabled = FALSE, op = ">=", value = 4,
                     tier = "general",
                     desc = "Min copy number for amplifications (CNV-only filter)."),

  # --- TMB / MSI / HRD -----------------------------------------------------
  TMB_score = list(enabled = FALSE, op = ">=", value = 10,
                    tier = "general",
                    desc = "Minimum TMB score (mut/Mb)."),
  TMB_category = list(enabled = FALSE, op = "%in%",
                       value = c("Intermediate","High"),
                       tier = "general",
                       desc = "TMB category."),
  MSI_High = list(enabled = FALSE, op = "==", value = "Yes",
                   tier = "general",
                   desc = "Require MSI-High calls only."),
  HRD_score = list(enabled = FALSE, op = ">=", value = 42,
                    tier = "general",
                    desc = "Minimum HRD score."),

  # --- Methylation TF --------------------------------------------------------
  Methylation_tumor_fraction_percentage = list(enabled = FALSE, op = ">=",
                                                  value = 0.001,
                                                  tier = "general",
                                                  desc = "Minimum methyl TF (%)."),

  # --- External annotations -------------------------------------------------
  COSMIC = list(enabled = FALSE, op = "has_value", value = TRUE,
                 tier = "general",
                 desc = "Keep rows with any COSMIC ID (annotated). FALSE = novel only."),
  dbSNP = list(enabled = FALSE, op = "has_value", value = FALSE,
                tier = "general",
                desc = paste("FALSE keeps rows missing dbSNP (likely-somatic novel calls).",
                              "TRUE keeps annotated.")),

  # --- gnomAD population AF (basic-tier QC floor + general-tier tunable) ---
  # The gnomAD rule sits in the basic tier as a hard QC floor at 0.001
  # (0.1%). Power users can tighten the cohort-wide threshold via a
  # general-tier rule by enabling a tighter override.
  gnomAD_AF = list(enabled = TRUE, op = "<=", value = 0.001,
                    tier = "basic",
                    desc = paste("Max population AF from gnomAD (0..1).",
                                  "Rows above are dropped as common variants.",
                                  "Basic-tier QC floor. Silently skipped if the",
                                  "gnomAD_AF column is absent — run",
                                  "ctdna_annotate_population_freq() first."))
)



# Print one entry compactly for display
.ctdna_filter_entry_str <- function(name, e) {
  status <- if (isTRUE(e$enabled)) "ON " else "off"
  val <- if (is.character(e$value) || is.factor(e$value)) {
    if (length(e$value) == 0) "<empty - SET TO USE>" else
      paste(shQuote(e$value), collapse = ", ")
  } else if (is.numeric(e$value)) {
    format(e$value, trim = TRUE)
  } else if (is.logical(e$value)) {
    as.character(e$value)
  } else {
    paste(format(e$value), collapse = ", ")
  }
  desc <- if (is.null(e$desc)) "" else paste0("  - ", e$desc)
  sprintf("  [%s]  %-40s  %s %s%s", status, name, e$op, val, desc)
}

#' Get / set / reset the alteration filter registry
#'
#' A column-by-column filter registry that names every filterable
#' Guardant Health Infinity column and pairs it with an
#' \code{enabled} flag, an operator, and a threshold / allowlist.
#' Apply the registry to a data frame with
#' \code{\link{ctdna_variant_filter}} (or pass via \code{filter_set} to
#' \code{\link{ctdna_oncoprint}}).
#'
#' Each entry has the form
#' \preformatted{list(enabled = TRUE/FALSE,
#'      op      = "\%in\%" / ">=" / "<=" / "==",
#'      value   = <allowlist or numeric threshold>,
#'      desc    = "<one-line description>")}
#'
#' @section Modes:
#' \itemize{
#'   \item No arguments: print the current registry.
#'   \item Single column name: return that one entry.
#'   \item Named arguments: update specific entries; pass either a
#'     full list, or an abbreviated list (e.g.
#'     \code{VAF_percentage = list(enabled = TRUE, value = 1.0)})
#'     - missing fields keep their current value.
#'   \item \code{.reset = TRUE}: restore the default registry.
#'   \item \code{load = list(...)}: bulk-load a snapshot.
#' }
#'
#' @param ... Column names (to get) or \code{name = list(...)} pairs
#'   (to set).
#' @param .reset If TRUE, restore the default registry.
#' @param load Optional snapshot list to load.
#' @return Invisibly the current registry (a named list), or the
#'   requested subset.
#' @examples
#' # Print the full registry
#' ctdna_filter_opts()
#'
#' # Get one entry
#' ctdna_filter_opts("VAF_percentage")
#'
#' # Turn on VAF filtering at 1.0%; restrict to High impact
#' ctdna_filter_opts(
#'   VAF_percentage    = list(enabled = TRUE, value = 1.0),
#'   Functional_impact = list(value = "High"))
#'
#' # Save & restore a snapshot
#' cfg <- ctdna_filter_opts()
#' ctdna_filter_opts(.reset = TRUE)
#' ctdna_filter_opts(load = cfg)
#' @export
ctdna_filter_opts <- function(..., .reset = FALSE, load = NULL) {

  # Initialise on first call
  if (is.null(.env[["filter_registry"]])) {
    .env[["filter_registry"]] <- .ctdna_filter_registry_defaults()
  }

  if (isTRUE(.reset)) {
    .env[["filter_registry"]] <- .ctdna_filter_registry_defaults()
    message("ctdna_filter_opts: registry reset to defaults.")
    return(invisible(.env[["filter_registry"]]))
  }
  if (!is.null(load)) {
    if (!is.list(load) || is.null(names(load)))
      stop("`load` must be a named list as returned by ctdna_filter_opts().")
    .env[["filter_registry"]] <- load
    message("ctdna_filter_opts: registry loaded (", length(load),
             " entries).")
    return(invisible(.env[["filter_registry"]]))
  }

  args <- list(...)
  reg <- .env[["filter_registry"]]

  # No args: print full registry
  if (length(args) == 0) {
    on_n  <- sum(vapply(reg, function(e) isTRUE(e$enabled), logical(1)))
    cat(sprintf("ctdna_filter_opts registry  (%d entries, %d enabled)\n",
                 length(reg), on_n))
    cat(strrep("=", 78), "\n", sep = "")
    for (nm in names(reg))
      cat(.ctdna_filter_entry_str(nm, reg[[nm]]), "\n", sep = "")
    cat(strrep("=", 78), "\n", sep = "")
    cat("Edit:   ctdna_filter_opts(<name> = list(enabled = TRUE, value = ...))\n")
    cat("Reset:  ctdna_filter_opts(.reset = TRUE)\n")
    return(invisible(reg))
  }

  # Positional getter: ctdna_filter_opts("VAF_percentage")
  if (length(names(args)) == 0 || all(names(args) == "")) {
    keys <- unlist(args)
    if (length(keys) == 1) return(reg[[keys]])
    return(reg[keys])
  }

  # Setter: each named arg becomes/replaces an entry
  for (nm in names(args)) {
    if (nm == "") next
    new_entry <- args[[nm]]
    if (!is.list(new_entry))
      stop(sprintf("Value for '%s' must be a list with enabled / op / value / desc fields.",
                    nm))
    cur <- if (nm %in% names(reg)) reg[[nm]] else
      list(enabled = TRUE, op = "%in%", value = character(0),
           desc = "user-added column")
    # Merge field by field — partial update OK
    for (k in c("enabled","op","value","desc")) {
      if (k %in% names(new_entry)) cur[[k]] <- new_entry[[k]]
    }
    reg[[nm]] <- cur
  }
  .env[["filter_registry"]] <- reg
  invisible(reg)
}


# Apply ONE filter entry to a data frame; return logical vector of keep-rows.
# Supported operators:
#   %in%       value is a character vector (allowlist)
#   >= / <=    value is a numeric threshold
#   ==         value is a single value to match
#   has_value  value is TRUE/FALSE: TRUE keeps rows where the column is
#              non-NA and non-empty; FALSE keeps rows where it IS missing.
#              Useful e.g. for dbSNP / COSMIC to keep only annotated
#              (has_value = TRUE) or only novel (has_value = FALSE) variants.
.ctdna_apply_entry <- function(df, col, entry, na_keep = TRUE) {
  if (!isTRUE(entry$enabled)) return(rep(TRUE, nrow(df)))
  if (!col %in% names(df))    return(rep(TRUE, nrow(df)))
  x <- df[[col]]

  op  <- entry$op  %||% "%in%"
  val <- entry$value

  keep <- switch(op,
    `%in%` = {
      if (length(val) == 0) rep(TRUE, length(x))
      else tolower(trimws(as.character(x))) %in%
            tolower(trimws(as.character(val)))
    },
    `>=`   = suppressWarnings(as.numeric(x) >= as.numeric(val)),
    `<=`   = suppressWarnings(as.numeric(x) <= as.numeric(val)),
    `==`   = tolower(trimws(as.character(x))) ==
              tolower(trimws(as.character(val))),
    has_value = {
      # TRUE -> keep rows with any value; FALSE -> keep rows that are missing
      has <- !(is.na(x) | as.character(x) %in% c("", "NA", "."))
      if (isTRUE(val)) has else !has
    },
    rep(TRUE, length(x))
  )
  # has_value is NOT subject to na_keep: NAs are exactly what it's about.
  if (na_keep && op != "has_value") keep[is.na(keep)] <- TRUE
  keep
}


# Internal: apply the alteration filter registry to a data frame.
# Walks the active filter registry (from ctdna_filter_opts) and applies
# every enabled entry. Used by ctdna_variant_filter(); not exported
# as a user-facing function (ctdna_variant_filter is the entry point).
#
# v0.26.0: REMOVED inline gnomAD annotation. Annotation is now a
# completely separate step — users call ctdna_annotate_population_freq()
# explicitly before filtering. If the gnomAD rule is enabled and the
# `gnomAD_AF` column is absent on `df`, the gnomAD rule is silently
# skipped with a verbose message.
.ctdna_apply_filter_internal <- function(df, overrides = NULL, na_keep = TRUE,
                                          gnomad_af_col = "gnomAD_AF",
                                          tier = NULL,
                                          verbose = FALSE) {
  if (is.null(.env[["filter_registry"]]))
    .env[["filter_registry"]] <- .ctdna_filter_registry_defaults()
  reg <- .env[["filter_registry"]]

  # v1.0.0: optional tier filter — only apply entries whose $tier matches.
  # tier = NULL -> apply ALL entries (backward compatible).
  # tier = "basic" -> apply only basic-tier entries (QC safety net step).
  # tier = "general" -> apply only general-tier entries (cohort step).
  if (!is.null(tier)) {
    keep_keys <- vapply(names(reg), function(k) {
      t <- reg[[k]]$tier %||% "general"
      t == tier
    }, logical(1))
    reg <- reg[keep_keys]
  }

  # Merge overrides on top of registry for this call
  if (!is.null(overrides)) {
    if (!is.list(overrides) || (length(overrides) > 0 && is.null(names(overrides))))
      stop("`overrides` must be a named list.")
    for (nm in names(overrides)) {
      cur <- if (nm %in% names(reg)) reg[[nm]] else
        list(enabled = TRUE, op = "%in%", value = character(0),
             desc = "override-added")
      new_entry <- overrides[[nm]]
      if (!is.list(new_entry))
        stop(sprintf("Override '%s' must be a list.", nm))
      for (k in c("enabled","op","value","desc"))
        if (k %in% names(new_entry)) cur[[k]] <- new_entry[[k]]
      reg[[nm]] <- cur
    }
  }

  n0 <- nrow(df)
  keep <- rep(TRUE, n0)
  applied <- character(0)
  # v0.23.0: gnomAD_AF is filtered AFTER annotation (which itself runs
  # on the survivors of the other filters, to minimize API payload).
  # Skip it in the first pass.
  filter_names <- setdiff(names(reg), "gnomAD_AF")
  for (nm in filter_names) {
    e <- reg[[nm]]
    if (!isTRUE(e$enabled)) next
    if (!nm %in% names(df)) next
    new_keep <- .ctdna_apply_entry(df, nm, e, na_keep = na_keep)
    keep <- keep & new_keep
    applied <- c(applied, nm)
    if (verbose) {
      n_drop <- sum(!new_keep, na.rm = TRUE)
      message(sprintf("  filter '%s' (%s) dropped %d rows",
                       nm, e$op, n_drop))
    }
  }
  out <- df[keep, , drop = FALSE]

  # v0.23.0/v0.24.0: gnomAD step. Decoupled into:
  #   (a) Annotation: need Chromosome / Position / Mut_nt + annotate=TRUE.
  #       Skipped if AF column is already fully populated (no NAs).
  #   (b) Filtering: need AF column present after step (a). Always
  #       runs if the filter is enabled.
  #   NA rows are NEVER dropped (NA-safe; only rows with AF > threshold
  #   are removed).
  # v0.26.0: gnomAD filter — column-read-only, no fetching.
  # If the gnomAD rule is enabled but the gnomAD_AF column is absent on
  # `df`, the rule is silently skipped (with a verbose message). To
  # populate gnomAD AF, run ctdna_annotate_population_freq() FIRST as a
  # separate step.
  gnomad_entry <- reg[["gnomAD_AF"]]
  gnomad_filter_on <- !is.null(gnomad_entry) &&
                       isTRUE(gnomad_entry$enabled) &&
                       nrow(out) > 0
  if (gnomad_filter_on) {
    has_freq_col <- gnomad_af_col %in% names(out)
    if (has_freq_col) {
      # Apply the gnomAD_AF filter — but route the operation through the
      # user-configurable column name. We temporarily mirror the column
      # so .ctdna_apply_entry can use its standard "gnomAD_AF" key.
      probe_df <- out
      probe_df[["gnomAD_AF"]] <- out[[gnomad_af_col]]
      new_keep <- .ctdna_apply_entry(probe_df, "gnomAD_AF", gnomad_entry,
                                       na_keep = TRUE)  # NA always kept
      n_drop <- sum(!new_keep, na.rm = TRUE)
      out <- out[new_keep, , drop = FALSE]
      applied <- c(applied, "gnomAD_AF")
      if (verbose)
        message(sprintf("  filter 'gnomAD_AF' (%s) dropped %d rows",
                         gnomad_entry$op, n_drop))
    } else if (verbose) {
      message("v0.26.0: gnomAD_AF filter enabled but the '", gnomad_af_col,
               "' column is not on this data frame. Skipped. ",
               "To populate gnomAD AF, run ctdna_annotate_population_freq() ",
               "BEFORE filtering — annotation and filtering are now separate steps.")
    }
  }

  if (verbose)
    message(sprintf("filter: %d -> %d rows (%d filters applied)",
                     n0, nrow(out), length(applied)))
  out
}


# ---- Coerce: schema-driven type / NA / transform cleanup --------------------
#
# Concept
# -------
# Inputs arrive in inconsistent shapes:
#   - patient IDs as integers when functions expect characters
#   - VAFs read as character strings ("0.5", "1.2%") when stats want numeric
#   - "NA" written as a string instead of a real NA
#   - expression as raw TPM when log2(TPM+1) is canonical
#
# ctdna_coerce(df) walks a small per-column registry and applies the
# declared type, NA rule, and transform. Called automatically inside
# every plot/helper that takes a data frame. Decoupled from
# ctdna_opts() — opts sets *intent* (column names), coerce acts on
# *data* (whatever frame is in hand).
# -----------------------------------------------------------------------------

# Per-column schema declaration.
#   type        - target R type. Coerced if it doesn't match.
#                 One of: "character", "numeric", "integer", "factor",
#                         "logical", "date".
#   na_action   - "keep" (default), "drop_row", "replace" (with na_value).
#   na_value    - replacement value when na_action = "replace".
#   transform   - "none" (default), "log2", "log2p1", "zscore", "minmax",
#                 "floor_loq", "percent_to_fraction".
# Anything not listed here passes through untouched.
.ctdna_coerce_registry <- function() list(
  # Patient / sample identity - always character to dodge integer / factor traps
  Patient_ID       = list(type = "character"),
  subject_id       = list(type = "character"),
  Customer_SampleID = list(type = "character"),
  GHSampleID       = list(type = "character"),
  Study_ID         = list(type = "character"),
  Visit_name       = list(type = "character"),
  time_point       = list(type = "character"),

  # Categorical clinical / biology - character so case-insensitive match works
  RECIST           = list(type = "character"),
  Cancertype       = list(type = "character"),
  dose             = list(type = "character"),
  cohort           = list(type = "character"),
  Variant_type     = list(type = "character"),
  Indel_type       = list(type = "character"),
  Somatic_status   = list(type = "character"),
  Functional_impact = list(type = "character"),
  Mutant_allele_status = list(type = "character"),
  Allelic_status   = list(type = "character"),
  Sample_status    = list(type = "character"),
  ClinVar          = list(type = "character"),
  CNV_type         = list(type = "character"),
  ctDNA_detection_status = list(type = "character"),
  TMB_category     = list(type = "character"),
  MSI_High         = list(type = "character"),
  panel            = list(type = "character"),
  gene             = list(type = "character"),
  Gene             = list(type = "character"),
  mutation_status  = list(type = "character"),
  alteration_type  = list(type = "character"),

  # Numeric assay / QC values
  methylTF         = list(type = "numeric"),
  maxVAF           = list(type = "numeric"),
  meanVAF          = list(type = "numeric"),
  gnomAD_AF        = list(type = "numeric"),
  gnomAD_version   = list(type = "character"),
  best_ratio       = list(type = "numeric"),
  ratio            = list(type = "numeric"),
  best_pct_change  = list(type = "numeric"),
  VAF_percentage   = list(type = "numeric"),
  Genomic_max_VAF_percentage = list(type = "numeric"),
  Mol_count        = list(type = "numeric"),
  Copy_number      = list(type = "numeric"),
  TMB_score        = list(type = "numeric"),
  HRD_score        = list(type = "numeric"),
  Methylation_tumor_fraction_percentage = list(type = "numeric"),
  cfDNA_ng         = list(type = "numeric"),
  Plasma_ml_input  = list(type = "numeric"),
  Plasma_ml_remaining = list(type = "numeric"),
  tmb              = list(type = "numeric"),
  expression       = list(type = "numeric"),
  score            = list(type = "numeric"),
  n.somatic        = list(type = "numeric", na_action = "replace", na_value = 0),
  cfdna.conc       = list(type = "numeric"),

  # Dates
  Bloodcoll_date   = list(type = "date"),
  Received_date    = list(type = "date"),
  Reported_date    = list(type = "date")
)

# Apply one coercion step; return the new column + a one-line change message
# (or NULL if no change).
.ctdna_coerce_one <- function(name, x, spec) {
  msgs <- character(0)
  orig_class <- class(x)[1]
  n_total <- length(x)

  # ---- 1. Type ----
  if (!is.null(spec$type)) {
    target <- spec$type
    same_type <- switch(target,
      character = is.character(x),
      numeric   = is.numeric(x) && !is.integer(x),
      integer   = is.integer(x),
      factor    = is.factor(x),
      logical   = is.logical(x),
      date      = inherits(x, "Date"),
      TRUE)
    if (!same_type) {
      n_na_before <- sum(is.na(x))
      x_new <- tryCatch(
        switch(target,
          character = as.character(x),
          numeric   = {
            # Strip common dressings: percent signs, commas, leading/trailing ws
            if (is.character(x) || is.factor(x)) {
              s <- as.character(x)
              s <- gsub("%$", "", trimws(s))
              s <- gsub(",", "", s)
              suppressWarnings(as.numeric(s))
            } else suppressWarnings(as.numeric(x))
          },
          integer   = suppressWarnings(as.integer(x)),
          factor    = as.factor(x),
          logical   = as.logical(x),
          date      = suppressWarnings(as.Date(x)),
          x),
        error = function(e) x)
      n_na_after <- sum(is.na(x_new))
      new_nas <- max(0, n_na_after - n_na_before)
      x <- x_new
      msgs <- c(msgs,
        sprintf("%s: %s -> %s%s", name, orig_class, target,
                if (new_nas > 0)
                  sprintf(" (%d NAs introduced)", new_nas) else
                  sprintf(" (%d rows)", n_total)))
    }
  }

  # ---- 2. NA handling ----
  if (!is.null(spec$na_action) && spec$na_action != "keep") {
    if (spec$na_action == "replace") {
      n_na <- sum(is.na(x))
      if (n_na > 0) {
        x[is.na(x)] <- spec$na_value
        msgs <- c(msgs,
          sprintf("%s: replaced %d NA values with %s",
                  name, n_na, format(spec$na_value)))
      }
    } else if (spec$na_action == "drop_row") {
      # Signal to caller via attribute; row-dropping happens at top level
      attr(x, ".ctdna_drop_na") <- TRUE
    }
  }

  # ---- 3. Transform ----
  if (!is.null(spec$transform) && spec$transform != "none") {
    x_new <- switch(spec$transform,
      log2     = log2(pmax(x, .Machine$double.eps)),
      log2p1   = log2(x + 1),
      zscore   = { mu <- mean(x, na.rm = TRUE); sd <- stats::sd(x, na.rm = TRUE)
                   if (is.na(sd) || sd == 0) x else (x - mu) / sd },
      minmax   = { rng <- range(x, na.rm = TRUE)
                   if (diff(rng) == 0) x else (x - rng[1]) / diff(rng) },
      floor_loq = pmax(x, .o("loq") %||% 0.01),
      percent_to_fraction = x / 100,
      x)
    if (!identical(x_new, x)) {
      msgs <- c(msgs,
        sprintf("%s: applied transform '%s'", name, spec$transform))
      x <- x_new
    }
  }

  list(value = x, messages = msgs)
}

#' Coerce a data frame to the canonical schema
#'
#' Walks a per-column schema registry and coerces each known column to
#' the declared type, applies the declared NA rule, and applies the
#' declared transform. Columns not in the registry pass through
#' untouched.
#'
#' This is called automatically by every plot / helper that takes a
#' data frame, so users rarely need to call it directly. It exists as
#' a user-facing function for ad-hoc cleaning and audit.
#'
#' @section Output log:
#' Prints a single clean block per call: a header naming the data frame
#' (if it has one) and one indented line per change. Use
#' \code{quiet = TRUE} to suppress.
#'
#' @param df A data frame.
#' @param quiet If TRUE, suppress the per-change log.
#' @param df_name Optional label used in the log header (e.g.
#'   "infinity_report").
#' @return The coerced data frame.
#' @examples
#' df <- data.frame(Patient_ID = 1:3,
#'                  VAF_percentage = c("0.5%","1.2%","3.0%"),
#'                  stringsAsFactors = FALSE)
#' ctdna_coerce(df, df_name = "demo")
#' @export
ctdna_coerce <- function(df, quiet = FALSE, df_name = NULL) {
  if (!is.data.frame(df)) return(df)
  reg <- .ctdna_coerce_registry()
  all_msgs <- character(0)
  drop_cols <- character(0)
  for (nm in intersect(names(df), names(reg))) {
    spec <- reg[[nm]]
    r <- .ctdna_coerce_one(nm, df[[nm]], spec)
    df[[nm]] <- r$value
    if (length(r$messages) > 0) all_msgs <- c(all_msgs, r$messages)
    if (isTRUE(attr(r$value, ".ctdna_drop_na"))) drop_cols <- c(drop_cols, nm)
  }
  # Apply any drop_row directives
  if (length(drop_cols) > 0) {
    n_before <- nrow(df)
    keep <- Reduce("&", lapply(drop_cols, function(nm) !is.na(df[[nm]])))
    df <- df[keep, , drop = FALSE]
    n_drop <- n_before - nrow(df)
    if (n_drop > 0)
      all_msgs <- c(all_msgs,
        sprintf("dropped %d row%s due to NA in: %s",
                n_drop, if (n_drop == 1) "" else "s",
                paste(drop_cols, collapse = ", ")))
  }
  # One clean log block
  if (!quiet && length(all_msgs) > 0) {
    hdr <- if (is.null(df_name)) "ctdna_coerce:" else
      sprintf("ctdna_coerce on %s:", df_name)
    message(hdr)
    for (m in all_msgs) message("  ", m)
  }
  df
}


# ---- Cleanup: sample/sequencing-level QC ------------------------------------
#
# Cleanup is the SAMPLE-level counterpart to filter_scheme (which acts on
# alterations). It is applied ONCE at load time inside ctdna_prepare() and
# ctdna_assemble_modalities() so downstream analysis works on a clean cohort
# without re-checking QC in every plot.
#
# What it does:
#   - For alteration frames: defers to the active filter system.
#   - For longitudinal ctDNA frames: drops samples failing sample-level QC
#     (cfDNA_ng, Plasma_ml_input, Mol_count, ctDNA_detection_status).
#   - For expression frames: drops samples whose total library size is below
#     a low-coverage floor.
#   - For TMB / IHC / clinical: pass-through (no sequencing QC applies).
#
# What it does NOT do:
#   - Per-variant filtering (that's filter_scheme's job).
#   - Type coercion (that's ctdna_coerce's job; called inside cleanup too).
# -----------------------------------------------------------------------------

#' Sample-level QC cleanup for ctDNA / expression frames
#'
#' Applied once at load time inside \code{\link{ctdna_prepare}} and
#' \code{\link{ctdna_assemble_modalities}} so the downstream pipeline
#' works on a QC-passed cohort. Distinct from \code{filter_scheme},
#' which acts on alteration calls; cleanup acts on the SAMPLE level
#' (which patients / visits / libraries to keep).
#'
#' @section What cleanup applies:
#' \itemize{
#'   \item ctDNA longitudinal frame - drops samples failing the
#'     standard QC thresholds: \code{cfDNA_ng >= cfdna_min},
#'     \code{Plasma_ml_input >= plasma_min}, \code{Mol_count >= mol_min},
#'     \code{ctDNA_detection_status} required if \code{require_detected = TRUE}.
#'     Columns that are not in the frame are silently skipped.
#'   \item Expression frame - drops samples whose total expression (sum of
#'     non-NA values) is below \code{expr_lib_min}. Useful for catching
#'     RNA-seq libraries with insufficient depth.
#'   \item Other frames (TMB / IHC / clinical) - pass through unchanged.
#' }
#'
#' @param df A data frame. The function auto-detects which kind of frame
#'   it is based on which canonical columns are present.
#' @param kind Optional hint - "ctdna", "expression", "alteration",
#'   "tmb", "ihc", "clinical", or "auto" (default).
#' @param cfdna_min Minimum cfDNA input in ng (default 5).
#' @param plasma_min Minimum plasma volume in mL (default 4).
#' @param mol_min Minimum unique molecule count (default 50).
#' @param require_detected If TRUE, only \code{ctDNA_detection_status = "Detected"}
#'   samples are kept. Default TRUE.
#' @param expr_lib_min Minimum library size for expression frames
#'   (default 1e5; samples with sum of expression values below this
#'   are dropped).
#' @param quiet If TRUE, suppress the per-frame change log.
#' @return The cleaned data frame.
#' @examples
#' sim <- ctdna_make_mock_study(n_patients = 20, seed = 1)
#' mm  <- ctdna_prepare(sim, verbose = FALSE)
#' ctdna_cleanup(mm$ctdna)
#' @export
ctdna_cleanup <- function(df, kind = c("auto","ctdna","expression",
                                          "alteration","tmb","ihc","clinical"),
                            cfdna_min = 5, plasma_min = 4, mol_min = 50,
                            require_detected = TRUE,
                            expr_lib_min = 1e5,
                            quiet = FALSE) {
  kind <- match.arg(kind)
  if (!is.data.frame(df)) return(df)

  if (kind == "auto") {
    cn <- names(df)
    # v0.41.0: derive names from ctdna_opts() so non-default
    # configurations still detect correctly.
    sub_col   <- .o("subject")
    time_col  <- .o("time")
    expr_col  <- .o("expression")
    tmb_col   <- .o("tmb")
    ihc_col   <- .o("ihc_score")
    kind <- if (all(c(sub_col, time_col) %in% cn) &&
                 !"Variant_type" %in% cn) "ctdna"
            else if (expr_col %in% cn) "expression"
            else if ("Variant_type" %in% cn || "Gene" %in% cn) "alteration"
            else if (tmb_col %in% cn) "tmb"
            else if (ihc_col %in% cn) "ihc"
            else "clinical"
  }

  n0 <- nrow(df)
  msgs <- character(0)

  if (kind == "ctdna") {
    # Sample-level sequencing QC
    if (require_detected && "ctDNA_detection_status" %in% names(df)) {
      keep <- as.character(df$ctDNA_detection_status) == "Detected" |
              is.na(df$ctDNA_detection_status)
      n_drop <- sum(!keep)
      df <- df[keep, , drop = FALSE]
      if (n_drop > 0) msgs <- c(msgs,
        sprintf("dropped %d sample(s) (ctDNA_detection_status != 'Detected')",
                n_drop))
    }
    for (col_min in list(c("cfDNA_ng", cfdna_min, "cfDNA_ng"),
                          c("Plasma_ml_input", plasma_min, "plasma volume"),
                          c("Mol_count", mol_min, "Mol_count"))) {
      col <- col_min[[1]]; thresh <- as.numeric(col_min[[2]])
      label <- col_min[[3]]
      if (col %in% names(df)) {
        x <- suppressWarnings(as.numeric(df[[col]]))
        keep <- is.na(x) | x >= thresh
        n_drop <- sum(!keep)
        df <- df[keep, , drop = FALSE]
        if (n_drop > 0) msgs <- c(msgs,
          sprintf("dropped %d sample(s) (%s < %g)", n_drop, label, thresh))
      }
    }
  } else if (kind == "expression") {
    # Drop low-library samples; assumes "subject_id" is per-row
    if (all(c("subject_id","expression") %in% names(df))) {
      lib <- stats::aggregate(expression ~ subject_id, data = df,
                                FUN = function(v) sum(v, na.rm = TRUE))
      bad <- lib$subject_id[lib$expression < expr_lib_min]
      if (length(bad) > 0) {
        n_drop <- sum(df$subject_id %in% bad)
        df <- df[!df$subject_id %in% bad, , drop = FALSE]
        msgs <- c(msgs,
          sprintf("dropped %d row(s) for %d low-library sample(s) (lib < %g)",
                  n_drop, length(bad), expr_lib_min))
      }
    }
  }
  # alteration / tmb / ihc / clinical: cleanup is a no-op
  # (alteration filtering belongs to filter_scheme)

  if (!quiet && length(msgs) > 0) {
    message(sprintf("ctdna_cleanup (%s): %d -> %d rows", kind, n0, nrow(df)))
    for (m in msgs) message("  ", m)
  }
  df
}


# ---- Three-level filter system: schemes -------------------------------------
#
# Concept
# -------
# Filtering has three levels:
#
# 1. DEFAULT GENERAL  -- shipped with the package. Always preserved as a
#    hard reset target so any modification can be undone.
#
# 2. MODIFIED GENERAL -- the user edits the general filter via
#    ctdna_filter_opts(). This becomes the active general filter and is
#    applied uniformly to all genes by every function.
#
# 3. NAMED FILTERING SCHEMES -- created via create_filtering_scheme().
#    A scheme is a composable rule tree (rule() / allOf() / anyOf() /
#    not()) optionally cloned from a template (built-in indication
#    template, built-in v1.0.0 scheme like TSG_Tier2/HRR14, or
#    a user scheme). Schemes registered with `frozen = TRUE` cannot be
#    overwritten until they are explicitly deleted.
#
# Use site
# --------
# Every function that touches alterations takes `filter_scheme`:
#   filter_scheme = NULL          -> current general filter (level 1 or 2)
#   filter_scheme = "myScheme"    -> resolve the named specific scheme
# -----------------------------------------------------------------------------

# Initialise scheme storage on first reference
.ctdna_init_schemes <- function() {
  if (is.null(.env[["filter_schemes"]]))
    .env[["filter_schemes"]] <- list()
  # Default-general is the factory registry, kept untouched as a reset target
  if (is.null(.env[["filter_general_default"]]))
    .env[["filter_general_default"]] <- .ctdna_filter_registry_defaults()
  # Active general is what ctdna_filter_opts() exposes
  if (is.null(.env[["filter_registry"]]))
    .env[["filter_registry"]] <- .ctdna_filter_registry_defaults()
  # Built-in indication templates (read-only, reserved names)
  if (is.null(.env[["filter_builtin_templates"]]))
    .env[["filter_builtin_templates"]] <- .ctdna_filter_builtin_templates()
  invisible(NULL)
}

# Built-in indication-specific filter scheme templates.
# Each template bundles indication-relevant gene sets + per-set overrides
# representing typical clinical-grade analysis defaults.
#
# These are READ-ONLY (their names are reserved); users clone them via
#   create_filtering_scheme(name = "my_NSCLC", template = "NSCLC", ...)
# to get an editable copy.
.ctdna_filter_builtin_templates <- function() list(

  # ---------------------------------------------------------------------------
  # v0.22.0 built-in indication templates.
  # Each template defines ONE gene set, "Drivers", based on a
  # literature-curated list of recurrent driver alterations in that
  # indication (TCGA / IntOGen / COSMIC CGC / NCCN / ESMO guidelines).
  # No overrides — every Drivers row inherits the active general filter.
  # Users should clone (template = "...") and add their own gene_sets /
  # overrides if they need finer-grained per-set rules.
  # ---------------------------------------------------------------------------

  NSCLC = list(
    description = "Non-small cell lung cancer (NSCLC) - recurrent drivers.",
    gene_sets = list(
      Drivers = c("EGFR","KRAS","ALK","ROS1","BRAF","MET","RET","ERBB2",
                  "NTRK1","NTRK2","NTRK3","NRG1","FGFR1","FGFR3",
                  "TP53","STK11","KEAP1","RB1","PTEN","CDKN2A",
                  "PIK3CA","NF1","RBM10","MGA","U2AF1")
    ),
    overrides = list()
  ),

  HNSCC = list(
    description = "Head and neck squamous cell carcinoma - recurrent drivers.",
    gene_sets = list(
      Drivers = c("TP53","CDKN2A","FAT1","NOTCH1","PIK3CA","KMT2D",
                  "NSD1","CASP8","HRAS","NFE2L2","TGFBR2","EPHA2",
                  "AJUBA","FBXW7","RAC1","CDKN2B","CCND1","EGFR",
                  "PTEN","HLA-A","HLA-B")
    ),
    overrides = list()
  ),

  SCLC = list(
    description = "Small cell lung cancer - recurrent drivers.",
    gene_sets = list(
      Drivers = c("TP53","RB1","MYC","MYCL","MYCN","PTEN","CREBBP",
                  "EP300","KMT2D","NOTCH1","NOTCH2","NOTCH3",
                  "NFIB","ASCL1","NEUROD1")
    ),
    overrides = list()
  ),

  BRCA = list(
    description = "Breast cancer - recurrent drivers.",
    gene_sets = list(
      Drivers = c("TP53","PIK3CA","PTEN","AKT1","ESR1","ERBB2",
                  "CDH1","GATA3","MAP3K1","MAP2K4","RUNX1","FOXA1",
                  "NF1","KMT2C","ARID1A","CBFB","CTCF",
                  "BRCA1","BRCA2","PALB2")
    ),
    overrides = list()
  ),

  CRC = list(
    description = "Colorectal cancer - recurrent drivers.",
    gene_sets = list(
      Drivers = c("APC","KRAS","NRAS","BRAF","TP53","SMAD4","PIK3CA",
                  "PTEN","FBXW7","TCF7L2","SMAD2","AMER1","ARID1A",
                  "SOX9","RNF43","ERBB2","ERBB3",
                  "MLH1","MSH2","MSH6","PMS2")
    ),
    overrides = list()
  ),

  mCRPC = list(
    description = "Metastatic castration-resistant prostate cancer - recurrent drivers.",
    gene_sets = list(
      Drivers = c("AR","TP53","PTEN","RB1","FOXA1","SPOP","CDK12",
                  "CTNNB1","BRCA2","BRCA1","ATM","CHEK2","MYC",
                  "ZBTB16","NCOR1","NCOR2")
    ),
    overrides = list()
  ),

  GBM = list(
    description = "Glioblastoma - recurrent drivers.",
    gene_sets = list(
      Drivers = c("EGFR","TP53","PTEN","NF1","CDKN2A","CDKN2B","RB1",
                  "PIK3CA","PIK3R1","ATRX","IDH1","IDH2","TERT",
                  "TP73","PDGFRA","MDM2","MDM4","CDK4","MET")
    ),
    overrides = list()
  )
)

# ---------------------------------------------------------------------------
# v0.26.0 GENE-SET-ONLY block catalog.
# Indication-blind. Each block restricts to genes in the listed pathway/set.
# Combined with the basic ruleset under AND semantics.
# All blocks use the basic rule template (overrides = list()); customize
# per-block rules in later iterations as needed.
# ---------------------------------------------------------------------------
.ctdna_filter_geneset_only_blocks <- function() list(
  HRR = list(
    description = "Homologous recombination repair (HRR) pathway genes.",
    genes = c("BRCA1","BRCA2","PALB2","ATM","CHEK2","RAD51","RAD51B","RAD51C",
              "RAD51D","RAD54L","BARD1","BRIP1","FANCA","FANCC","FANCD2",
              "FANCE","FANCF","FANCG","FANCI","FANCL","FANCM","NBN",
              "MRE11A","RAD50","BLM","WRN","XRCC2","XRCC3"),
    overrides = list()
  ),
  TSGs = list(
    description = "Common tumor suppressor genes (TSGs).",
    genes = c("TP53","RB1","PTEN","APC","NF1","NF2","VHL","CDKN2A","CDKN2B",
              "STK11","BRCA1","BRCA2","SMAD4","CDH1","WT1","BAP1","ARID1A",
              "ARID1B","ARID2","SETD2","KMT2D","KMT2C","TSC1","TSC2","FBXW7"),
    overrides = list()
  ),
  RTK = list(
    description = "Receptor tyrosine kinases (RTK) / MAPK pathway.",
    genes = c("EGFR","ERBB2","ERBB3","ERBB4","MET","ALK","ROS1","RET","KIT",
              "PDGFRA","PDGFRB","FGFR1","FGFR2","FGFR3","FGFR4","NTRK1",
              "NTRK2","NTRK3","KRAS","NRAS","HRAS","BRAF","RAF1","MAP2K1",
              "MAP2K2","NF1","PTPN11"),
    overrides = list()
  ),
  Cell_Cycle = list(
    description = "Cell cycle regulators.",
    genes = c("CDKN2A","CDKN2B","CDKN1A","CDKN1B","CDK4","CDK6","CCND1",
              "CCND2","CCND3","CCNE1","RB1","E2F1","E2F3","MDM2","MDM4",
              "TP53","TP73"),
    overrides = list()
  ),
  TP53 = list(
    description = "TP53 pathway (TP53 + co-regulators).",
    genes = c("TP53","TP63","TP73","MDM2","MDM4","ATM","ATR","CHEK1","CHEK2"),
    overrides = list()
  ),
  MMR = list(
    description = "Mismatch repair (MMR) pathway.",
    genes = c("MLH1","MSH2","MSH6","PMS2","MSH3","MLH3","POLE","POLD1",
              "EXO1","RFC1","PCNA"),
    overrides = list()
  ),
  PI3K = list(
    description = "PI3K / AKT / mTOR pathway.",
    genes = c("PIK3CA","PIK3CB","PIK3R1","PIK3R2","PTEN","AKT1","AKT2","AKT3",
              "MTOR","TSC1","TSC2","RHEB","RPTOR","RICTOR","STK11","INPP4B",
              "PHLPP1","PHLPP2")
  )
)

# ---------------------------------------------------------------------------
# v0.26.0 INDICATION-ONLY block catalog.
# Gene-blind. Each block matches rows by Cancertype only, applying its
# rule overrides (or, with overrides = list(), the basic ruleset).
# All blocks use the basic rule template for now; tune in later iterations.
# ---------------------------------------------------------------------------
.ctdna_filter_indication_only_blocks <- function() list(
  NSCLC_only  = list(description = "NSCLC rows (any gene), basic rules.",
                     indication = "NSCLC",  overrides = list()),
  HNSCC_only  = list(description = "HNSCC rows (any gene), basic rules.",
                     indication = "HNSCC",  overrides = list()),
  SCLC_only   = list(description = "SCLC rows (any gene), basic rules.",
                     indication = "SCLC",   overrides = list()),
  BRCA_only   = list(description = "Breast cancer rows (any gene), basic rules.",
                     indication = "BRCA",   overrides = list()),
  CRC_only    = list(description = "Colorectal cancer rows (any gene), basic rules.",
                     indication = "CRC",    overrides = list()),
  mCRPC_only  = list(description = "mCRPC rows (any gene), basic rules.",
                     indication = "mCRPC",  overrides = list()),
  GBM_only    = list(description = "Glioblastoma rows (any gene), basic rules.",
                     indication = "GBM",    overrides = list())
)

# ---------------------------------------------------------------------------
# v0.26.0 MUTATION-TYPE block catalog.
# Matches rows by Variant_type and (for CNV blocks) CNV_type. Uses the
# basic ruleset for now; tune per-mutation-type in later iterations.
# ---------------------------------------------------------------------------
.ctdna_filter_mutation_type_blocks <- function() list(
  SNV             = list(description = "Single-nucleotide variants.",
                          variant_type = "SNV", overrides = list()),
  Indel           = list(description = "Small insertions and deletions.",
                          variant_type = "Indel", overrides = list()),
  CNV             = list(description = "Copy-number variants (all subtypes).",
                          variant_type = "CNV", overrides = list()),
  Fusion          = list(description = "Gene fusions.",
                          variant_type = "Fusion", overrides = list()),
  LGR             = list(description = "Large genomic rearrangements.",
                          variant_type = "LGR", overrides = list()),
  Focal_Amp       = list(description = "Focal amplifications (CNV subtype).",
                          variant_type = "CNV", cnv_type = "focal_amplification",
                          overrides = list()),
  Homozygous_Del  = list(description = "Homozygous deletions (CNV subtype).",
                          variant_type = "CNV", cnv_type = "homozygous_deletion",
                          overrides = list()),
  LoH             = list(description = "Loss of heterozygosity (CNV subtype).",
                          variant_type = "CNV", cnv_type = "loh_deletion",
                          overrides = list())
)

# Unified catalog combining all four block categories.
# Returns a named list keyed by block name; each entry has a `category`
# slot indicating which group it belongs to.
.ctdna_filter_block_catalog <- function() {
  out <- list()
  # Indication × gene_set: derive from the existing builtin templates
  ind_x_gs <- .ctdna_filter_builtin_templates()
  for (nm in names(ind_x_gs)) {
    e <- ind_x_gs[[nm]]
    out[[nm]] <- list(category   = "indication_gene_set",
                       indication = nm,
                       gene_set   = e$gene_sets$Drivers,
                       overrides  = e$overrides %||% list(),
                       description = e$description)
  }
  # Gene-set-only
  gs <- .ctdna_filter_geneset_only_blocks()
  for (nm in names(gs)) {
    e <- gs[[nm]]
    out[[nm]] <- list(category    = "gene_set_only",
                       gene_set    = e$genes,
                       overrides   = e$overrides %||% list(),
                       description = e$description)
  }
  # Indication-only
  ind <- .ctdna_filter_indication_only_blocks()
  for (nm in names(ind)) {
    e <- ind[[nm]]
    out[[nm]] <- list(category    = "indication_only",
                       indication  = e$indication,
                       overrides   = e$overrides %||% list(),
                       description = e$description)
  }
  # Mutation-type
  mt <- .ctdna_filter_mutation_type_blocks()
  for (nm in names(mt)) {
    e <- mt[[nm]]
    out[[nm]] <- list(category     = "mutation_type",
                       variant_type = e$variant_type,
                       cnv_type     = e$cnv_type %||% NA_character_,
                       overrides    = e$overrides %||% list(),
                       description  = e$description)
  }
  out
}

# Reserved scheme/block names — users can't shadow these.
# v0.26.0: includes the full block catalog plus the reserved literals
# "basic", "none", and the deprecated aliases "general" and "default".
.ctdna_reserved_scheme_names <- function()
  c("basic","none","general","default",
    names(.ctdna_filter_block_catalog()))

# Merge a single override onto a registry entry (partial update OK)
.ctdna_merge_entry <- function(base, override) {
  for (k in c("enabled","op","value","desc"))
    if (k %in% names(override)) base[[k]] <- override[[k]]
  base
}


#' Reset the general filter back to the package default
#'
#' Restores every entry of the active general filter to its factory
#' default. Specific schemes are unaffected (dynamic schemes will
#' re-resolve against the new general state on their next use; frozen
#' schemes are independent by design).
#'
#' @return Invisibly, the restored general registry.
#' @examples
#' ctdna_filter_opts(VAF_percentage = list(enabled = TRUE, value = 5))
#' ctdna_filter_reset_general()
#' ctdna_filter_opts("VAF_percentage")$enabled  # FALSE again
#' @seealso \code{\link{ctdna_filter_opts}}
#' @export
ctdna_filter_reset_general <- function() {
  .ctdna_init_schemes()
  .env[["filter_registry"]] <- .env[["filter_general_default"]]
  message("ctdna_filter_reset_general: general filter restored to factory defaults.")
  invisible(.env[["filter_registry"]])
}



#' List, get, delete, save, or load filter schemes
#'
#' \itemize{
#'   \item \code{ctdna_filter_schemes_list()} -- returns a data frame with
#'     one row per scheme.
#'   \item \code{ctdna_filter_scheme_get(name)} -- returns the scheme object.
#'   \item \code{ctdna_filter_scheme_delete(name)} -- removes the scheme.
#'   \item \code{ctdna_filter_scheme_save(file)} -- saves all schemes
#'     (and the general default) to an RDS file.
#'   \item \code{ctdna_filter_scheme_load(file)} -- restores from an RDS
#'     file written by \code{ctdna_filter_scheme_save}.
#' }
#' To create a scheme, use \code{\link{create_filtering_scheme}}.
#' @param name Scheme name.
#' @param file RDS file path.
#' @return List or data frame depending on the function.
#' @examples
#' create_filtering_scheme(
#'   rule(Gene = "TP53"),
#'   name = "demo")
#' ctdna_filter_schemes_list()
#' ctdna_filter_scheme_delete("demo")
#' @name ctdna_filter_scheme
NULL


#' List all filter blocks in the v0.26.0 catalog
#'
#' Returns a data frame describing every block in the catalog: the four
#' built-in categories (indication x gene-set, gene-set-only,
#' indication-only, mutation-type) plus the special reserved names and
#' any user-defined schemes.
#'
#' This is the v0.26.0 replacement for the older
#' \code{ctdna_filter_schemes_list()} (which is still available and shows
#' only legacy template-style schemes).
#'
#' @param category Optional category filter. One of: \code{"all"}
#'   (default), \code{"indication_gene_set"}, \code{"gene_set_only"},
#'   \code{"indication_only"}, \code{"mutation_type"}, \code{"reserved"},
#'   \code{"user"}.
#' @return A data frame with one row per block. Columns:
#'   \code{name}, \code{category}, \code{indication}, \code{n_genes},
#'   \code{variant_type}, \code{description}.
#' @seealso \code{\link{ctdna_filter_scheme_show}} to inspect a single
#'   block's rules; \code{\link{ctdna_variant_filter}} to apply blocks.
#' @examples
#' ctdna_filter_schemes_list()                       # everything
#' ctdna_filter_schemes_list("gene_set_only")        # just gene-set blocks
#' ctdna_filter_schemes_list("indication_gene_set")  # just indication x gene-set
#' @export
ctdna_filter_schemes_list <- function(category = c("all","indication_gene_set",
                                                     "gene_set_only",
                                                     "indication_only",
                                                     "mutation_type",
                                                     "reserved",
                                                     "user")) {
  category <- match.arg(category)
  .ctdna_init_schemes()
  cat0 <- .ctdna_filter_block_catalog()

  # Catalog rows
  cat_rows <- if (length(cat0) > 0)
    do.call(rbind, lapply(names(cat0), function(nm) {
      e <- cat0[[nm]]
      data.frame(
        name        = nm,
        category    = e$category,
        indication  = e$indication %||% NA_character_,
        n_genes     = if (!is.null(e$gene_set)) length(e$gene_set) else NA_integer_,
        variant_type= e$variant_type %||% NA_character_,
        description = e$description %||% NA_character_,
        stringsAsFactors = FALSE)
    })) else NULL

  # Reserved rows (basic, none, deprecated aliases)
  reserved_rows <- data.frame(
    name = c("basic","none","general","default"),
    category    = c("reserved","reserved","deprecated","deprecated"),
    indication  = NA_character_,
    n_genes     = NA_integer_,
    variant_type= NA_character_,
    description = c(
      "Floor: the baseline rules. Always applied beneath other blocks.",
      "Explicit 'no filter' — same effect as filter_scheme = NULL.",
      "DEPRECATED alias for 'basic'. Will be removed in a future version.",
      paste("DEPRECATED alias resolving to all 7 indication x gene-set blocks.",
            "Replace with an explicit vector of indication names.")),
    stringsAsFactors = FALSE)

  # User-defined schemes (legacy v0.x)
  s <- .env[["filter_schemes"]]
  user_rows <- if (length(s) > 0)
    do.call(rbind, lapply(s, function(x) data.frame(
      name        = x$name,
      category    = "user",
      indication  = NA_character_,
      n_genes     = length(unique(unlist(x$gene_sets, use.names = FALSE))),
      variant_type= NA_character_,
      description = x$description %||% NA_character_,
      stringsAsFactors = FALSE
    ))) else NULL

  # v1.0.0 built-in schemes (HRR14 / TSG_Tier2)
  builtin_v1 <- tryCatch(.ctdna_builtin_schemes(), error = function(e) list())
  v1_builtin_rows <- if (length(builtin_v1) > 0)
    do.call(rbind, lapply(names(builtin_v1), function(nm) {
      sch <- builtin_v1[[nm]]
      data.frame(
        name        = nm,
        category    = sch$category %||% "gene_set",
        indication  = NA_character_,
        n_genes     = NA_integer_,
        variant_type= NA_character_,
        description = sch$description %||% NA_character_,
        stringsAsFactors = FALSE)
    })) else NULL

  # v1.0.0 user-registered schemes (exclude any that match built-in names)
  .ctdna_init_schemes_v1()
  user_v1 <- .env[["filter_schemes_v1"]] %||% list()
  # Drop entries that duplicate a built-in name (the builtin scheme
  # constructors auto-register themselves under the same name).
  user_v1 <- user_v1[setdiff(names(user_v1), names(builtin_v1))]
  v1_user_rows <- if (length(user_v1) > 0)
    do.call(rbind, lapply(names(user_v1), function(nm) {
      sch <- user_v1[[nm]]
      data.frame(
        name        = nm,
        category    = sch$category %||% "user",
        indication  = NA_character_,
        n_genes     = NA_integer_,
        variant_type= NA_character_,
        description = sch$description %||% NA_character_,
        stringsAsFactors = FALSE)
    })) else NULL

  all_rows <- rbind(cat_rows, reserved_rows, user_rows, v1_builtin_rows, v1_user_rows)

  out <- switch(category,
    all                 = all_rows,
    indication_gene_set = cat_rows[cat_rows$category == "indication_gene_set", , drop = FALSE],
    gene_set_only       = cat_rows[cat_rows$category == "gene_set_only", , drop = FALSE],
    indication_only     = cat_rows[cat_rows$category == "indication_only", , drop = FALSE],
    mutation_type       = cat_rows[cat_rows$category == "mutation_type", , drop = FALSE],
    reserved            = reserved_rows,
    user                = rbind(user_rows, v1_user_rows) %||% reserved_rows[0, , drop = FALSE])
  rownames(out) <- NULL
  out
}


#' Show the full rule set for a single filter block
#'
#' Prints (and invisibly returns) the contents of one block from the
#' v0.26.0 catalog: its category, scope (indication / gene set / mutation
#' type), and any rule overrides on top of the basic ruleset. Use it
#' to audit exactly what a block does before applying it.
#'
#' @param name Name of the block, as listed in
#'   \code{\link{ctdna_filter_schemes_list}}.
#' @return The block list (invisibly).
#' @seealso \code{\link{ctdna_filter_schemes_list}},
#'   \code{\link{ctdna_variant_filter}}.
#' @examples
#' ctdna_filter_scheme_show("NSCLC")
#' ctdna_filter_scheme_show("HRR")
#' @export
ctdna_filter_scheme_show <- function(name) {
  if (!is.character(name) || length(name) != 1L)
    stop("`name` must be a single block/scheme name.")

  # v1.0.0: built-in schemes (HRR14 / TSG_Tier2) and user-registered
  # v1.0.0 schemes — print the rule tree via the ctdna_rule print method.
  builtin_v1 <- tryCatch(.ctdna_builtin_schemes(), error = function(e) list())
  .ctdna_init_schemes_v1()
  user_v1 <- .env[["filter_schemes_v1"]] %||% list()
  if (name %in% names(builtin_v1)) {
    cat(sprintf("Built-in v1.0.0 scheme: %s\n", name))
    print(builtin_v1[[name]])
    return(invisible(builtin_v1[[name]]))
  }
  if (name %in% names(user_v1)) {
    cat(sprintf("User-defined v1.0.0 scheme: %s\n", name))
    print(user_v1[[name]])
    return(invisible(user_v1[[name]]))
  }

  # v0.26 block catalog (the richer metadata view of indication / gene-set /
  # mutation-type blocks). Take precedence over the v0.26 scheme container
  # for same-named entries so callers see block-level metadata (category,
  # n_genes, variant_type) on names like "NSCLC".
  cat0 <- .ctdna_filter_block_catalog()
  if (name %in% names(cat0)) {
    e <- cat0[[name]]
    cat(sprintf("Filter block: %s\n", name))
    cat(sprintf("  category    : %s\n", e$category))
    cat(sprintf("  description : %s\n", e$description %||% "(none)"))
    if (!is.null(e$indication))
      cat(sprintf("  indication  : %s\n", e$indication))
    if (!is.null(e$gene_set))
      cat(sprintf("  n_genes     : %d (%s%s)\n",
                  length(e$gene_set),
                  paste(utils::head(e$gene_set, 10), collapse = ", "),
                  if (length(e$gene_set) > 10) ", ..." else ""))
    if (!is.null(e$variant_type))
      cat(sprintf("  variant_type: %s\n", e$variant_type))
    if (!is.null(e$cnv_type) && !is.na(e$cnv_type))
      cat(sprintf("  cnv_type    : %s\n", e$cnv_type))
    cat(sprintf("  overrides   : %s\n",
                if (length(e$overrides) == 0) "(none — inherits basic ruleset)"
                else paste(names(e$overrides), collapse = ", ")))
    return(invisible(e))
  }

  # v0.26 user-defined schemes and built-in templates
  # (gene_sets + per-set overrides container — fallback for names that
  # don't appear in the richer block catalog).
  .ctdna_init_schemes()
  user_schemes  <- .env[["filter_schemes"]]
  builtin_tmpls <- .env[["filter_builtin_templates"]]
  if (name %in% names(user_schemes) || name %in% names(builtin_tmpls)) {
    if (name %in% names(user_schemes)) {
      s <- user_schemes[[name]]
      mode_str <- if (isTRUE(s$frozen)) "frozen" else "dynamic"
      desc_str <- "(user-defined)"
    } else {
      s <- builtin_tmpls[[name]]
      mode_str <- "builtin"
      desc_str <- s$description %||% "(built-in template)"
    }
    cat(sprintf("Filter scheme: %s  [%s]\n", name, mode_str))
    cat(strrep("=", 70), "\n", sep = "")
    cat(desc_str, "\n\n")
    if (length(s$gene_sets) == 0)
      cat("(no gene sets)\n")
    for (set_nm in names(s$gene_sets)) {
      genes <- s$gene_sets[[set_nm]]
      cat(sprintf("  [%s] (%d gene%s):  %s\n",
                  set_nm, length(genes), if (length(genes) == 1) "" else "s",
                  paste(genes, collapse = ", ")))
      ovr <- s$overrides[[set_nm]]
      if (length(ovr) > 0) {
        for (col_nm in names(ovr)) {
          e <- ovr[[col_nm]]
          val_str <- if (is.null(e$value)) "<unset>" else
            paste(format(e$value), collapse = ", ")
          cat(sprintf("       override: %s %s %s\n",
                      col_nm, e$op %||% "<inherit>", val_str))
        }
      } else {
        cat("       (no overrides - inherits general)\n")
      }
      cat("\n")
    }
    cat("All genes NOT in any of these sets use the active general filter.\n")
    return(invisible(s))
  }

  stop(sprintf("No scheme or block named '%s'. ", name),
       "See ctdna_filter_schemes_list() for the catalog.")
}


#' @rdname ctdna_filter_scheme
#' @export
ctdna_filter_scheme_get <- function(name) {
  .ctdna_init_schemes()
  .ctdna_init_schemes_v1()
  # v1.0.0 schemes (where create_filtering_scheme() registers since v0.36.0)
  if (name %in% names(.env[["filter_schemes_v1"]]))
    return(.env[["filter_schemes_v1"]][[name]])
  # v1.0.0 built-in schemes (HRR14, TSG_Tier2)
  builtin_v1 <- tryCatch(.ctdna_builtin_schemes(), error = function(e) list())
  if (name %in% names(builtin_v1))
    return(builtin_v1[[name]])
  # v0.26 user schemes (legacy)
  if (name %in% names(.env[["filter_schemes"]]))
    return(.env[["filter_schemes"]][[name]])
  # v0.26 built-in indication templates
  if (name %in% names(.env[["filter_builtin_templates"]])) {
    s <- .env[["filter_builtin_templates"]][[name]]
    s$builtin <- TRUE; s$name <- name
    return(s)
  }
  stop(sprintf("No filter scheme named '%s'. See ctdna_filter_schemes_list().",
                name))
}

#' @rdname ctdna_filter_scheme
#' @export
ctdna_filter_scheme_delete <- function(name) {
  .ctdna_init_schemes()
  .ctdna_init_schemes_v1()
  found <- FALSE
  # v1.0.0 registry (primary since v0.36.0). Frozen schemes can still be
  # deleted -- frozen only blocks overwrite, not explicit removal.
  if (name %in% names(.env[["filter_schemes_v1"]])) {
    .env[["filter_schemes_v1"]][[name]] <- NULL
    found <- TRUE
  }
  # v0.26 registry (legacy)
  if (name %in% names(.env[["filter_schemes"]])) {
    .env[["filter_schemes"]][[name]] <- NULL
    found <- TRUE
  }
  if (!found) {
    message(sprintf("Scheme '%s' does not exist; nothing to delete.", name))
    return(invisible(NULL))
  }
  message(sprintf("Scheme '%s' deleted.", name))
  invisible(NULL)
}

#' @rdname ctdna_filter_scheme
#' @param name For \code{save}: optionally the name of a single user
#'   scheme to save (NULL = save all user schemes plus the active
#'   general filter).
#' @export
ctdna_filter_scheme_save <- function(file, name = NULL) {
  .ctdna_init_schemes()
  all_schemes <- .env[["filter_schemes"]]

  if (!is.null(name)) {
    if (!name %in% names(all_schemes))
      stop(sprintf("No user scheme named '%s'. See ctdna_filter_schemes_list(\"user\").",
                    name))
    payload <- list(
      tag             = "ctdnaTM_scheme_payload",
      schemes         = all_schemes[name],
      general_active  = .env[["filter_registry"]],
      general_default = .env[["filter_general_default"]]
    )
    saveRDS(payload, file = file)
    message(sprintf("Saved scheme '%s' to '%s'.", name, file))
  } else {
    payload <- list(
      tag             = "ctdnaTM_scheme_payload",
      schemes         = all_schemes,
      general_active  = .env[["filter_registry"]],
      general_default = .env[["filter_general_default"]]
    )
    saveRDS(payload, file = file)
    message(sprintf("Saved %d scheme(s) + general filter to '%s'.",
                     length(payload$schemes), file))
  }
  invisible(file)
}

#' @rdname ctdna_filter_scheme
#' @param only Character vector of scheme names to load from the file
#'   (NULL = all). Errors clearly if any name is absent from the file.
#' @param overwrite If TRUE, silently replace any same-name scheme
#'   already in the session. Default FALSE (errors instead, listing
#'   the collisions).
#' @param prefix Optional string prepended to every loaded scheme's
#'   name. Lets you load schemes alongside existing same-name ones.
#'   Only \code{[A-Za-z0-9_]} characters are allowed.
#' @export
ctdna_filter_scheme_load <- function(file, only = NULL, overwrite = FALSE,
                                       prefix = NULL) {
  .ctdna_init_schemes()
  payload <- readRDS(file)
  # File-format sanity check
  is_scheme_payload <- is.list(payload) && (
    identical(payload$tag, "ctdnaTM_scheme_payload") ||
      "schemes" %in% names(payload))
  is_opts_payload <- is.list(payload) &&
    identical(payload$tag, "ctdnaTM_opts_payload")
  if (is_opts_payload)
    stop("This file is an opts payload, not a scheme payload. ",
          "Use ctdna_opts_load() instead.", call. = FALSE)
  if (!is_scheme_payload)
    stop("File does not look like a ctdna_filter_scheme_save() payload.",
          call. = FALSE)

  schemes <- payload$schemes %||% list()
  if (length(schemes) == 0) {
    message(sprintf("File '%s' contains no schemes.", file))
    return(invisible(list()))
  }

  # Filter by `only`
  if (!is.null(only)) {
    missing <- setdiff(only, names(schemes))
    if (length(missing) > 0)
      stop(sprintf("Schemes not found in file: %s. Available: %s",
                    paste(shQuote(missing), collapse = ", "),
                    paste(shQuote(names(schemes)), collapse = ", ")),
            call. = FALSE)
    schemes <- schemes[only]
  }

  # Optional prefix
  if (!is.null(prefix)) {
    if (!is.character(prefix) || length(prefix) != 1 ||
        !grepl("^[A-Za-z0-9_]+$", prefix))
      stop("`prefix` must be a single string of letters / digits / underscores.",
            call. = FALSE)
    new_names <- paste0(prefix, names(schemes))
    # Renamed schemes must not collide with reserved built-ins
    reserved <- .ctdna_reserved_scheme_names()
    bad <- intersect(new_names, reserved)
    if (length(bad) > 0)
      stop(sprintf("After prefixing, these names collide with reserved scheme names: %s. ",
                    paste(shQuote(bad), collapse = ", ")),
            "Pick a different prefix.", call. = FALSE)
    # Also update the embedded scheme$name fields
    for (i in seq_along(schemes)) schemes[[i]]$name <- new_names[i]
    names(schemes) <- new_names
  }

  # Collision check
  existing <- intersect(names(schemes), names(.env[["filter_schemes"]]))
  if (length(existing) > 0 && !overwrite) {
    stop(sprintf("Loading would overwrite existing scheme(s): %s. ",
                  paste(shQuote(existing), collapse = ", ")),
          "Pass overwrite = TRUE to replace, or prefix = '...' to load alongside.",
          call. = FALSE)
  }

  # Apply
  for (nm in names(schemes))
    .env[["filter_schemes"]][[nm]] <- schemes[[nm]]
  message(sprintf("Loaded %d scheme(s) from '%s': %s",
                   length(schemes), file,
                   paste(names(schemes), collapse = ", ")))
  invisible(schemes)
}


#' Save and load ctdna_opts() configuration
#'
#' \code{ctdna_opts_save(file)} writes the current option registry to
#' an RDS file. \code{ctdna_opts_load(file)} restores it.
#'
#' Opts are session-global - there is only one set at a time -
#' so neither function takes a name argument.
#'
#' @param file RDS file path.
#' @return Invisibly, the file path (for save) or the loaded opts list.
#' @examples
#' ctdna_opts(stat_position = "caption", legend_position = "bottom")
#' f <- tempfile(fileext = ".rds")
#' ctdna_opts_save(f)
#' ctdna_opts(.reset = TRUE)
#' ctdna_opts_load(f)
#' ctdna_opts("stat_position")
#' @name ctdna_opts_io
NULL

#' @rdname ctdna_opts_io
#' @export
ctdna_opts_save <- function(file) {
  payload <- list(
    tag = "ctdnaTM_opts_payload",
    opts = .env$cfg
  )
  saveRDS(payload, file = file)
  message(sprintf("Saved ctdna_opts to '%s'.", file))
  invisible(file)
}

#' @rdname ctdna_opts_io
#' @export
ctdna_opts_load <- function(file) {
  payload <- readRDS(file)
  is_opts_payload <- is.list(payload) &&
    identical(payload$tag, "ctdnaTM_opts_payload")
  is_scheme_payload <- is.list(payload) &&
    identical(payload$tag, "ctdnaTM_scheme_payload")
  if (is_scheme_payload)
    stop("This file is a scheme payload, not an opts payload. ",
          "Use ctdna_filter_scheme_load() instead.", call. = FALSE)
  if (!is_opts_payload)
    stop("File does not look like a ctdna_opts_save() payload.",
          call. = FALSE)
  bad <- setdiff(names(payload$opts), names(.defaults))
  if (length(bad) > 0)
    stop("Unknown option key(s) in file: ", paste(bad, collapse = ", "),
          call. = FALSE)
  for (k in names(payload$opts)) .opts_validate(k, payload$opts[[k]])
  for (k in names(payload$opts)) .env$cfg[[k]] <- payload$opts[[k]]
  message(sprintf("Loaded ctdna_opts from '%s'.", file))
  invisible(.env$cfg)
}


# Resolve a scheme (or NULL = general) against the active general filter
# into a per-gene-set table of effective filter entries.
# Returns a list with: $general (resolved registry), $gene_sets (named
# list of gene-vectors), $resolved (named list keyed by gene-set name,
# each a full registry merged from general + that set's overrides).
.ctdna_resolve_scheme <- function(scheme_name = NULL) {
  .ctdna_init_schemes()
  active_general <- .env[["filter_registry"]]
  if (is.null(scheme_name) || identical(scheme_name, "general")) {
    # No specific scheme: gene-agnostic, one set "all"
    return(list(general  = active_general,
                gene_sets = NULL,
                resolved  = NULL,
                scheme    = NULL))
  }
  # User-defined scheme first, then fall through to built-ins
  user_schemes  <- .env[["filter_schemes"]]
  builtin_tmpls <- .env[["filter_builtin_templates"]]
  if (scheme_name %in% names(user_schemes)) {
    s <- user_schemes[[scheme_name]]
    # Frozen schemes use the snapshot, dynamic schemes use the active general
    general <- if (isTRUE(s$frozen)) s$general_snapshot else active_general
  } else if (scheme_name %in% names(builtin_tmpls)) {
    # Built-in template: treat as dynamic against active general
    tmpl <- builtin_tmpls[[scheme_name]]
    s <- list(name      = scheme_name,
              gene_sets = tmpl$gene_sets,
              overrides = tmpl$overrides,
              frozen    = FALSE,
              builtin   = TRUE,
              description = tmpl$description)
    general <- active_general
  } else {
    stop(sprintf("No filter scheme named '%s'. See ctdna_filter_schemes_list().",
                  scheme_name))
  }
  # Build resolved per-set registries
  resolved <- lapply(names(s$gene_sets), function(set_nm) {
    reg <- general
    ovr <- s$overrides[[set_nm]]
    if (is.list(ovr)) {
      for (nm in names(ovr)) {
        cur <- if (nm %in% names(reg)) reg[[nm]] else
          list(enabled = TRUE, op = "%in%", value = character(0),
               desc = "override-added")
        reg[[nm]] <- .ctdna_merge_entry(cur, ovr[[nm]])
      }
    }
    reg
  })
  names(resolved) <- names(s$gene_sets)
  list(general  = general,
       gene_sets = s$gene_sets,
       resolved  = resolved,
       scheme    = s)
}




# ===========================================================================
# v1.0.0 filter_apply — RULE-SYSTEM ARCHITECTURE
# ===========================================================================
# Replaces the v0.26.0 block catalog with the composable rule system in
# rules.R + rules_library.R. Hierarchy (single-claim semantics):
#
#   1. user-defined schemes  (claim their rows)
#   2. indication x gene_set schemes (claim their rows)
#   3. gene-set schemes  (claim their rows, includes HRR14, TSG_Tier2)
#   4. indication schemes  (claim their rows)
#   5. general  (applies ONLY to rows NOT claimed above)
#   6. basic  (applies to ALL surviving rows as the QC safety net)
#
# A row is "claimed" by the FIRST step whose scheme it matches. Once
# claimed, it skips steps 2-5 and goes directly to step 6 (basic).
# Unclaimed rows go through general (step 5) then basic (step 6).
# ===========================================================================

# Deprecation warning cache (retained for other callers that need one-shot warnings)
.ctdna_filter_dep_warned <- new.env(parent = emptyenv())
.ctdna_warn_once <- function(key, msg) {
  if (!isTRUE(.ctdna_filter_dep_warned[[key]])) {
    warning(msg, call. = FALSE)
    .ctdna_filter_dep_warned[[key]] <- TRUE
  }
}

# -----------------------------------------------------------------------------
# Cancertype dictionary check (LOUD warning on mismatch)
# -----------------------------------------------------------------------------
#
# Scans the Cancertype column against the dictionary in ctdna_opts and
# warns LOUDLY if any value isn't recognized. The warning instructs the
# user to update ctdna_opts(cancertype_dictionary = ...) to map their
# vendor strings to canonical codes.
.ctdna_check_cancertypes <- function(df) {
  col <- .o("col_cancertype") %||% "Cancertype"
  if (!col %in% names(df)) return(invisible())
  vals <- unique(stats::na.omit(as.character(df[[col]])))
  if (length(vals) == 0) return(invisible())
  dict <- .o("cancertype_dictionary") %||% list()
  # All known cancertype strings (canonical names + vendor synonyms)
  known <- unique(c(names(dict), unlist(dict, use.names = FALSE)))
  unknown <- setdiff(vals, known)
  if (length(unknown) > 0) {
    warning(
      "\n",
      "  *** ctdna_variant_filter: UNRECOGNIZED CANCERTYPE VALUES ***\n",
      "  Found in column '", col, "': ",
      paste(shQuote(unknown), collapse = ", "), "\n",
      "  These will NOT match any indication rule (rule_ind_NSCLC, etc.).\n",
      "  To fix, update the cancertype dictionary in ctdna_opts(), e.g.:\n",
      "    ctdna_opts(cancertype_dictionary = list(\n",
      "      NSCLC = c('NSCLC', '", unknown[1], "'),\n",
      "      ...))\n",
      "  Or rename the values in df$", col, " to canonical codes.\n",
      call. = FALSE, immediate. = TRUE)
  }
  invisible()
}

# -----------------------------------------------------------------------------
# Resolve a filter_scheme argument to a list of scheme objects
# -----------------------------------------------------------------------------
#
# Input: filter_scheme can be a character vector of names, a single
# ctdna_filtering_scheme object, or a list of such objects.
# Output: a named list of ctdna_filtering_scheme objects (some may
#   carry $category = "user" if user-defined).
.ctdna_resolve_schemes <- function(filter_scheme) {
  # Single scheme object -> wrap in list
  if (inherits(filter_scheme, "ctdna_filtering_scheme")) {
    nm <- filter_scheme$name %||% "<unnamed>"
    return(setNames(list(filter_scheme), nm))
  }
  # Already a list of schemes?
  if (is.list(filter_scheme) && !is.character(filter_scheme)) {
    out <- list()
    for (i in seq_along(filter_scheme)) {
      x <- filter_scheme[[i]]
      if (!inherits(x, "ctdna_filtering_scheme")) {
        elem_class <- paste(class(x), collapse = "/")
        stop("List element ", i, " of `filter_scheme` is not a ",
              "ctdna_filtering_scheme (got class: ", elem_class, "). ",
              "If you meant to pass extra arguments to the plot ",
              "function, note that some args partial-match ",
              "`filter_scheme` (R's argument-matching rules); use the ",
              "full argument name to avoid ambiguity.",
              call. = FALSE)
      }
      nm <- x$name %||% names(filter_scheme)[i] %||% paste0("scheme_", i)
      out[[nm]] <- x
    }
    return(out)
  }
  # Character vector of names -> resolve each
  if (!is.character(filter_scheme))
    stop("`filter_scheme` must be NULL, a character vector of names, ",
          "or one or more ctdna_filtering_scheme objects.")
  builtin <- .ctdna_builtin_schemes()
  .ctdna_init_schemes_v1()
  user <- .env[["filter_schemes_v1"]] %||% list()
  legacy <- .env[["filter_schemes"]] %||% list()
  out <- list()
  for (nm in filter_scheme) {
    if (nm %in% names(builtin)) {
      out[[nm]] <- builtin[[nm]]
    } else if (nm %in% names(user)) {
      out[[nm]] <- user[[nm]]
    } else if (nm %in% names(legacy)) {
      # v0.26 gene-sets+overrides scheme loaded from an old RDS file.
      # Wrap a thin pseudo-scheme so ctdna_variant_filter() can still
      # find and report it. The unified create_filtering_scheme() no
      # longer produces this shape.
      sch <- legacy[[nm]]
      out[[nm]] <- structure(
        list(type = "scheme", name = nm,
             description = sch$description %||% "legacy scheme",
             category = "user", root = NULL, legacy = sch),
        class = c("ctdna_filtering_scheme","ctdna_rule"))
    } else {
      stop(sprintf("Unknown filter scheme name: '%s'. ", nm),
            "Available built-ins: ", paste(names(builtin), collapse = ", "),
            ". See ctdna_filter_schemes_list() for the full catalog.",
            call. = FALSE)
    }
  }
  out
}

# -----------------------------------------------------------------------------
# Hierarchy ordering: bucket each scheme by category
# -----------------------------------------------------------------------------
#
# Returns a list with 4 slots in apply order:
#   $user, $indication_gene_set, $gene_set, $indication
# Each slot is a named list of ctdna_filtering_scheme objects.
.ctdna_bucket_schemes_by_hierarchy <- function(schemes) {
  out <- list(user = list(),
              indication_gene_set = list(),
              gene_set = list(),
              indication = list())
  for (nm in names(schemes)) {
    sch <- schemes[[nm]]
    cat0 <- sch$category %||% "user"
    if (cat0 == "user")                  out$user[[nm]] <- sch
    else if (cat0 == "indication_gene_set") out$indication_gene_set[[nm]] <- sch
    else if (cat0 == "gene_set")         out$gene_set[[nm]] <- sch
    else if (cat0 == "indication")       out$indication[[nm]] <- sch
    else                                  out$user[[nm]] <- sch
  }
  out
}

# -----------------------------------------------------------------------------
# Evaluate a single scheme against df, returning a keep-mask.
# NA in the mask is coerced to FALSE (NA means "rule couldn't decide" —
# rows fall through to general/basic if no other scheme claims them).
# -----------------------------------------------------------------------------
.ctdna_apply_scheme <- function(scheme, df, gnomad_af_col, verbose) {
  if (!inherits(scheme, "ctdna_filtering_scheme"))
    stop(".ctdna_apply_scheme: not a filtering scheme.")
  n <- nrow(df)
  if (n == 0) return(logical(0))
  # Legacy scheme wrapper -> defer to old internal applier
  if (!is.null(scheme$legacy)) {
    return(rep(TRUE, n))  # legacy schemes treated as no-op in v1.0.0;
    # tests using them get the deprecation warning elsewhere.
  }
  mask <- .eval_rule(scheme, df)
  mask[is.na(mask)] <- FALSE
  mask
}

# -----------------------------------------------------------------------------
# Apply the basic-or-general flat registry to df (used for general + basic
# steps in the hierarchy). Wraps the existing .ctdna_apply_filter_internal
# which already implements the registry semantics.
# -----------------------------------------------------------------------------
.ctdna_apply_registry <- function(df, gnomad_af_col, verbose, tier = NULL) {
  if (nrow(df) == 0) return(df)
  .ctdna_apply_filter_internal(df, overrides = NULL,
                                 gnomad_af_col = gnomad_af_col,
                                 tier = tier,
                                 verbose = verbose)
}


#' Apply filtering schemes to an alteration data frame
#'
#' v1.0.0 filtering: routes rows through a 6-step hierarchy combining
#' user-defined schemes, built-in schemes, the general cohort-wide
#' filter, and basic QC. Schemes are composable rule trees built with
#' \code{\link{rule}}, \code{\link{allOf}}, \code{\link{anyOf}},
#' \code{\link{not}}, and \code{\link{create_filtering_scheme}}.
#'
#' \strong{Hierarchy (in apply order):}
#' \enumerate{
#'   \item User-defined schemes (category \code{"user"})
#'   \item Indication x gene-set schemes (category \code{"indication_gene_set"})
#'   \item Gene-set schemes (category \code{"gene_set"}) -- includes the
#'     built-in \code{HRR14} and \code{TSG_Tier2}
#'   \item Indication schemes (category \code{"indication"})
#'   \item General -- the cohort-wide tunable filter (gnomAD, COSMIC,
#'     ClinVar, dbSNP, etc.). Applies ONLY to rows NOT claimed by any
#'     scheme in steps 1-4.
#'   \item Basic -- the QC safety net (Sample_status, cfDNA, Plasma_ml,
#'     ctDNA_detection_status, very-low gnomAD floor). Applies to ALL
#'     surviving rows.
#' }
#'
#' \strong{Single-claim semantics:} once a scheme in steps 1-4 claims a
#' row, that row skips remaining steps 2-5 and goes straight to basic.
#' This avoids double-jeopardy from competing cohort rules. The basic
#' safety net always fires regardless of which scheme claimed the row.
#'
#' \strong{Cancertype dictionary:} when any scheme uses
#' \code{rule_ind_*}, this function first scans the \code{Cancertype}
#' column and emits a LOUD warning if any value isn't recognized in
#' \code{ctdna_opts("cancertype_dictionary")}. Unrecognized values
#' won't match any indication rule. Fix by updating the dictionary or
#' by renaming column values to canonical codes.
#'
#' \strong{Annotation:} this function never makes network calls. If a
#' rule references \code{gnomAD_AF} but the column is absent, that rule
#' silently evaluates to FALSE for the row (no error). Run
#' \code{\link{ctdna_annotate_population_freq}} BEFORE filtering to
#' populate the column.
#'
#' @param df Alteration data frame.
#' @param filter_scheme One of:
#'   \itemize{
#'     \item \emph{missing} -- the session default from
#'       \code{ctdna_opts("default_scheme")}, typically \code{"basic"}.
#'     \item \code{NULL} or \code{"none"} -- no filtering at all.
#'     \item \code{"basic"} -- skip steps 1-5, apply basic only.
#'     \item Character vector of scheme names -- e.g.
#'       \code{c("HRR14","TSG_Tier2")} or \code{c("HRR","NSCLC")}.
#'     \item A \code{ctdna_filtering_scheme} object from
#'       \code{\link{create_filtering_scheme}}.
#'     \item A list of such objects.
#'   }
#' @param gene_col Column holding gene symbols. Auto-detected.
#' @param cancer_type_col Column holding indication labels.
#'   Default \code{"Cancertype"}.
#' @param gnomad_af_col Column holding gnomAD allele frequency values.
#'   Default \code{"gnomAD_AF"}.
#' @param sample_qc_enable Per-call override for sample-level QC. One of
#'   \code{NULL} (default; respect \code{ctdna_opts("sample_qc_filter")}),
#'   \code{TRUE} (force QC on for this call), or \code{FALSE} (force
#'   QC off for this call). Use \code{FALSE} for quick "what would my
#'   filter return if QC weren't dropping these rows?" tests without
#'   touching the session option.
#' @param verbose If \code{TRUE}, print per-step summaries.
#' @return Filtered data frame with attribute \code{"ctdna_filter_used"}
#'   recording which schemes claimed which rows.
#' @seealso \code{\link{rule}}, \code{\link{create_filtering_scheme}},
#'   \code{\link{HRR14}}, \code{\link{TSG_Tier2}},
#'   \code{\link{ctdna_filter_schemes_list}}.
#' @examples
#' \dontrun{
#' # Apply built-in HRR scheme alone
#' df_hrr <- ctdna_variant_filter(df, filter_scheme = "HRR14")
#'
#' # Multiple schemes — each claims its rows, unclaimed go to general
#' df_combined <- ctdna_variant_filter(df,
#'   filter_scheme = c("HRR14","TSG_Tier2"))
#'
#' # User scheme passed directly
#' my <- create_filtering_scheme(
#'   rule(Gene = c("EGFR","KRAS")),
#'   rule(Variant_type = "SNV"),
#'   name = "my_egfr_kras")
#' df_my <- ctdna_variant_filter(df, filter_scheme = my)
#' }
#' @export
.filter_apply_df <- function(df,
                                 filter_scheme = NULL,
                                 gene_col = NULL,
                                 cancer_type_col = "Cancertype",
                                 gnomad_af_col = "gnomAD_AF",
                                 sample_qc_enable = NULL,
                                 verbose = FALSE) {
  if (!is.data.frame(df)) stop("`df` must be a data frame.")
  msg <- function(...) if (isTRUE(verbose)) message(...)

  # ----------------------------------------------------------------------
  # STAGE 0: sample-level QC (runs in ALL modes including filter_scheme=NULL)
  # ----------------------------------------------------------------------
  # v0.38.0: per-call override via `sample_qc_enable` (TRUE/FALSE) takes
  # precedence over ctdna_opts("sample_qc_filter"). When NULL (default),
  # falls back to the opts setting. This makes quick "QC on/off" testing
  # a one-arg change without touching session options.
  qc_active <- if (!is.null(sample_qc_enable)) {
    isTRUE(sample_qc_enable)
  } else {
    isTRUE(.o("sample_qc_filter"))
  }
  n_pre_qc <- nrow(df)
  qc_dropped <- 0L
  if (qc_active && n_pre_qc > 0L) {
    qc_mask <- rep(TRUE, n_pre_qc)
    status_col <- .o("sample_qc_status_col")
    cfdna_col  <- .o("sample_qc_cfdna_col")
    plasma_col <- .o("sample_qc_plasma_col")
    if (!is.null(status_col) && status_col %in% names(df)) {
      pass_val <- .o("sample_qc_status_pass")
      ok <- as.character(df[[status_col]]) == pass_val
      ok[is.na(ok)] <- FALSE
      qc_mask <- qc_mask & ok
    } else if (!is.null(status_col)) {
      .ctdna_warn_once(paste0("qc_missing_", status_col),
        sprintf("QC: column '%s' not in input -- skipping that QC check.",
                 status_col))
    }
    if (!is.null(cfdna_col) && cfdna_col %in% names(df)) {
      min_val <- .o("sample_qc_cfdna_min")
      v <- suppressWarnings(as.numeric(df[[cfdna_col]]))
      ok <- !is.na(v) & v >= min_val
      qc_mask <- qc_mask & ok
    } else if (!is.null(cfdna_col)) {
      .ctdna_warn_once(paste0("qc_missing_", cfdna_col),
        sprintf("QC: column '%s' not in input -- skipping that QC check.",
                 cfdna_col))
    }
    if (!is.null(plasma_col) && plasma_col %in% names(df)) {
      min_val <- .o("sample_qc_plasma_min")
      v <- suppressWarnings(as.numeric(df[[plasma_col]]))
      ok <- !is.na(v) & v >= min_val
      qc_mask <- qc_mask & ok
    } else if (!is.null(plasma_col)) {
      .ctdna_warn_once(paste0("qc_missing_", plasma_col),
        sprintf("QC: column '%s' not in input -- skipping that QC check.",
                 plasma_col))
    }
    qc_dropped <- sum(!qc_mask)
    df <- df[qc_mask, , drop = FALSE]
    rownames(df) <- NULL
    if (qc_dropped > 0L)
      msg(sprintf("ctdna_variant_filter: QC dropped %d / %d row(s).",
                   qc_dropped, n_pre_qc))
  }

  # ----------------------------------------------------------------------
  # Mode 1: filter_scheme = NULL  ->  QC-only (no variant filtering)
  # ----------------------------------------------------------------------
  if (is.null(filter_scheme)) {
    msg(".filter_apply_df: filter_scheme = NULL -> QC only, no variant filtering.")
    attr(df, "ctdna_filter_used") <- list(mode = "qc_only",
                                            scheme = NULL,
                                            n_in = n_pre_qc,
                                            n_out = nrow(df),
                                            qc_dropped = qc_dropped)
    return(df)
  }

  # --- 'none' still supported as an explicit no-op string --------------
  if (is.character(filter_scheme) && length(filter_scheme) == 1L &&
      filter_scheme == "none") {
    msg(".filter_apply_df: filter_scheme = 'none' -> no filtering.")
    attr(df, "ctdna_filter_used") <- list(mode = "none", scheme = NULL,
                                            n_in = nrow(df), n_out = nrow(df))
    return(df)
  }

  # --- Auto-coerce on entry ---------------------------------------------
  df <- ctdna_coerce(df, quiet = !isTRUE(verbose))
  n0 <- nrow(df)
  if (n0 == 0L) {
    attr(df, "ctdna_filter_used") <- list(mode = "empty", n_in = 0, n_out = 0)
    return(df)
  }

  # ----------------------------------------------------------------------
  # Mode 2: filter_scheme = "basic"  ->  scheme_basic applied to ALL rows
  # ----------------------------------------------------------------------
  if (is.character(filter_scheme) && length(filter_scheme) == 1L &&
      filter_scheme == "basic") {
    msg("ctdna_variant_filter: BASIC mode (scheme_basic on all rows).")
    basic_sch <- scheme_basic()
    mask <- .eval_rule(basic_sch, df)
    mask[is.na(mask)] <- FALSE
    out <- df[mask, , drop = FALSE]
    rownames(out) <- NULL
    attr(out, "ctdna_filter_used") <- list(
      mode = "basic", scheme = "basic",
      n_in = n0, n_out = nrow(out))
    return(out)
  }

  # ----------------------------------------------------------------------
  # Mode 3: filter_scheme = vector of named schemes (each has a gene scope)
  # ----------------------------------------------------------------------
  # Per-row routing:
  #   * row's Gene is in scheme S's gene scope  ->  apply scheme S's rules
  #   * row's Gene matches multiple schemes      ->  first match wins
  #     (order = order passed by user)
  #   * row's Gene matches no listed scheme      ->  apply scheme_basic
  schemes <- .ctdna_resolve_schemes(filter_scheme)
  msg(sprintf("ctdna_variant_filter: %d scheme(s): %s",
                length(schemes), paste(names(schemes), collapse = ", ")))

  # Pre-compute scope membership for each scheme.
  # A scheme's "scope" rule is what defines which rows it claims
  # (gene-set membership for HRR14/TSG23/etc.). If a scheme has NO scope,
  # the scheme applies to ALL genes -- it claims every still-unclaimed
  # row. This is the documented semantics: "if there is no gene set
  # indicated in the scheme definition, it means the scheme rules
  # apply to all genes" (no gene will be subject to scheme_basic).
  scope_masks <- vector("list", length(schemes))
  names(scope_masks) <- names(schemes)
  for (nm in names(schemes)) {
    sch <- schemes[[nm]]
    if (is.null(sch$scope)) {
      scope_masks[[nm]] <- rep(TRUE, nrow(df))
    } else {
      sm <- .eval_rule(sch$scope, df)
      sm[is.na(sm)] <- FALSE
      scope_masks[[nm]] <- sm
    }
  }

  # Walk schemes in order; first-match-wins for the claim mask.
  claimed   <- rep(FALSE, n0)
  claimed_by <- rep(NA_character_, n0)
  pass      <- rep(FALSE, n0)

  for (nm in names(schemes)) {
    sch <- schemes[[nm]]
    to_claim <- scope_masks[[nm]] & !claimed
    n_claim  <- sum(to_claim)
    if (n_claim == 0) {
      msg(sprintf("  '%s': in-scope = %d, new claims = 0", nm,
                   sum(scope_masks[[nm]])))
      next
    }
    # Evaluate the scheme's full rule tree on the claimed rows
    full_mask <- .eval_rule(sch, df)
    full_mask[is.na(full_mask)] <- FALSE
    to_pass <- to_claim & full_mask
    claimed[to_claim]    <- TRUE
    claimed_by[to_claim] <- nm
    pass[to_pass]        <- TRUE
    msg(sprintf("  '%s': in-scope = %d, new claims = %d, passed = %d, dropped = %d",
                 nm, sum(scope_masks[[nm]]), n_claim, sum(to_pass),
                 n_claim - sum(to_pass)))
  }

  # Apply scheme_basic to all unclaimed rows
  unclaimed <- !claimed
  n_unc <- sum(unclaimed)
  if (n_unc > 0) {
    msg(sprintf("  unclaimed -> scheme_basic: %d rows", n_unc))
    basic_sch <- scheme_basic()
    basic_mask <- rep(FALSE, n0)
    basic_mask_sub <- .eval_rule(basic_sch, df[unclaimed, , drop = FALSE])
    basic_mask_sub[is.na(basic_mask_sub)] <- FALSE
    basic_mask[unclaimed] <- basic_mask_sub
    pass[basic_mask & unclaimed] <- TRUE
    claimed_by[basic_mask & unclaimed] <- "scheme_basic"
    msg(sprintf("  unclaimed: passed = %d, dropped = %d",
                 sum(basic_mask & unclaimed),
                 n_unc - sum(basic_mask & unclaimed)))
  }

  out <- df[pass, , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "ctdna_filter_used") <- list(
    mode = "composed", schemes = names(schemes),
    n_in = n0, n_out = nrow(out),
    claimed_counts = if (any(!is.na(claimed_by)))
      table(claimed_by, useNA = "ifany") else NULL)
  out
}


# ---- Oncoprint --------------------------------------------------------------
#
# v0.39.0 — clean rewrite. The public function below is just orchestration;
# all heavy lifting lives in R/oncoprint.R. Both engines (ComplexHeatmap and
# ggplot fallback) share the same input pipeline, palette, scheme, wrap,
# annotation-label, and cumulative-legend behaviour.
#
# Alteration palette (single ground truth -- .oncoprint_alt_palette in
# R/oncoprint.R):
#   Focal_Amp      red          focal_amplification
#   Amp            pink         amplification / aneuploid_amplification
#   Homozygous_Del blue         homozygous_deletion
#   LOH            light blue   loh_deletion
#   Missense       green        missense
#   Truncating     orange       splice_donor / splice_acceptor / nonsense /
#                                  frameshift / stop_lost / start_lost
#   InFrame        pale green   inframe_insertion / _duplication / _deletion /
#                                  _indel
#   Promoter       ochre        promoter
#   LGR            purple       Variant_type = "LGR"
#   Fusion         brown        Variant_type = "Fusion"
#   Other          dark blue    catch-all (synonymous, splice_region, ...)
# -----------------------------------------------------------------------------


#' Print method for ctdna_oncoprint
#' @param x a ctdna_oncoprint object.
#' @param ... unused.
#' @export
print.ctdna_oncoprint <- function(x, ...) {
  if (identical(x$engine, "ComplexHeatmap")) {
    side <- if (is.null(x$legend_side)) "bottom" else x$legend_side
    has_title    <- !is.null(x$title)    && nzchar(as.character(x$title))
    has_subtitle <- !is.null(x$subtitle) && nzchar(as.character(x$subtitle))
    has_caption  <- !is.null(x$caption)  && nzchar(as.character(x$caption))

    # Base draw args (shared with both rendering paths)
    draw_args <- list(
      object                 = x$heatmap,
      heatmap_legend_side    = side,
      annotation_legend_side = side,
      merge_legend           = TRUE,
      padding                = grid::unit(c(8, 8, 8, 8), "mm"))
    if (!is.null(x$gene_set_legend))
      draw_args$annotation_legend_list <- list(x$gene_set_legend)

    # v0.42.2 bugfix: previously the CH branch ignored title/subtitle/
    # caption entirely. Now we render them in a grid layout around the
    # heatmap. Simple path (title only, no subtitle/caption): use CH's
    # native column_title — cheap, no extra viewport gymnastics.
    if (!has_subtitle && !has_caption) {
      if (has_title) {
        draw_args$column_title    <- x$title
        draw_args$column_title_gp <- grid::gpar(fontsize = 14, fontface = "bold")
      }
      do.call(ComplexHeatmap::draw, draw_args)
    } else {
      # Custom layout: title (bold) / subtitle (grey) / heatmap (CH draw
      # called with newpage = FALSE) / caption (small italic grey).
      .ctdna_oncoprint_draw_with_extras(x$title, x$subtitle, x$caption,
                                         has_title, has_subtitle,
                                         has_caption, draw_args)
    }
  } else {
    print(x$plot)
  }
  invisible(x)
}


# Internal helper for print.ctdna_oncoprint when subtitle or caption
# is non-NULL. Builds a vertical grid layout, draws each piece into
# its slot. The CH heatmap goes into the only "null"-sized row so it
# absorbs the remaining vertical space.
.ctdna_oncoprint_draw_with_extras <- function(title, subtitle, caption,
                                                 has_title, has_subtitle,
                                                 has_caption, draw_args) {
  heights <- list()
  if (has_title)    heights <- c(heights, list(grid::unit(2.0, "lines")))
  if (has_subtitle) heights <- c(heights, list(grid::unit(1.5, "lines")))
  heights <- c(heights, list(grid::unit(1, "null")))   # heatmap (absorbs)
  if (has_caption)  heights <- c(heights, list(grid::unit(1.5, "lines")))

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(nrow    = length(heights),
                                ncol    = 1,
                                heights = do.call(grid::unit.c, heights))))
  row <- 1L
  if (has_title) {
    grid::pushViewport(grid::viewport(layout.pos.row = row))
    grid::grid.text(title, y = grid::unit(0.4, "npc"),
                    gp = grid::gpar(fontsize = 14, fontface = "bold"))
    grid::popViewport()
    row <- row + 1L
  }
  if (has_subtitle) {
    grid::pushViewport(grid::viewport(layout.pos.row = row))
    grid::grid.text(subtitle, y = grid::unit(0.5, "npc"),
                    gp = grid::gpar(fontsize = 11, col = "grey40"))
    grid::popViewport()
    row <- row + 1L
  }
  # Heatmap
  grid::pushViewport(grid::viewport(layout.pos.row = row))
  draw_args$newpage <- FALSE
  do.call(ComplexHeatmap::draw, draw_args)
  grid::popViewport()
  row <- row + 1L
  if (has_caption) {
    grid::pushViewport(grid::viewport(layout.pos.row = row))
    grid::grid.text(caption, x = grid::unit(2, "mm"), just = "left",
                    gp = grid::gpar(fontsize = 9, col = "grey40",
                                      fontface = "italic"))
    grid::popViewport()
  }
  grid::popViewport()
}


# ---- Publication-ready tables -----------------------------------------------

#' Publication-quality table as a ggplot
#'
#' Renders a data frame, matrix, or table as a styled, publication-ready
#' table panel using ggplot2. The result is a real ggplot, so you can:
#' \itemize{
#'   \item save it with \code{ggplot2::ggsave()} at any DPI
#'   \item combine it with plots via \code{patchwork}
#'   \item include it directly in an R Markdown report
#' }
#'
#' Style features:
#' \itemize{
#'   \item Title (bold) / subtitle (grey) / caption (small grey, left-aligned)
#'   \item Bold header row with light grey fill and a black underline rule
#'   \item Alternating row stripes for readability
#'   \item Top and bottom horizontal rules in black
#'   \item Optional first-column header (row labels) bolded
#'   \item Optional diagonal-cell tint for contingency / concordance tables
#'   \item Auto-formatted numeric cells
#' }
#'
#' @param x A data frame, matrix, or \code{table}.
#' @param title,subtitle,caption Display labels (NULL to omit).
#' @param row_names Optional character vector of row labels. Default
#'   \code{NULL} uses \code{rownames(x)} if meaningful; set \code{FALSE}
#'   to suppress the row-label column.
#' @param highlight_diag If TRUE and the data portion is square, the
#'   diagonal cells get a light blue tint (useful for contingency
#'   tables to emphasise agreement).
#' @param number_format \code{sprintf} format string for numeric cells
#'   (default \code{"\%g"}). Set to \code{NULL} to keep values as-is.
#' @param font_family Font family (default \code{""} = the device default).
#' @param base_size Base font size for cell text.
#' @return A ggplot object.
#' @examples
#' tab <- as.data.frame(matrix(c(12, 3, 1, 9), 2, 2,
#'                              dimnames = list(MR = c("Response","NonResponse"),
#'                                              RECIST = c("CR/PR","SD/PD/NE"))))
#' ctdna_table(tab,
#'             title    = "MR x RECIST contingency",
#'             subtitle = "n = 25; Cohen's kappa = 0.62",
#'             caption  = "MR computed at C5D1; cutoff = 50% TF drop",
#'             highlight_diag = TRUE)
#' @export
ctdna_table <- function(x,
                         title          = NULL,
                         subtitle       = NULL,
                         caption        = NULL,
                         row_names      = NULL,
                         highlight_diag = FALSE,
                         number_format  = "%g",
                         font_family    = "",
                         base_size      = 11) {

  # --- Coerce input
  if (inherits(x, "table"))  x <- as.data.frame.matrix(x)
  if (is.matrix(x))          x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!is.data.frame(x))
    stop("`x` must be a data frame, matrix, or table.")

  # --- Row-name column handling
  has_rn <- TRUE
  if (isFALSE(row_names)) {
    has_rn <- FALSE
  } else if (is.null(row_names)) {
    if (is.null(rownames(x)) ||
        all(rownames(x) == as.character(seq_len(nrow(x)))))
      has_rn <- FALSE
    else
      rn <- rownames(x)
  } else {
    rn <- as.character(row_names)
  }
  if (has_rn) {
    new_x <- data.frame(rn_col = rn, check.names = FALSE,
                         stringsAsFactors = FALSE)
    names(new_x) <- " "                                     # one-space header
    for (cn in names(x)) new_x[[cn]] <- x[[cn]]
    x <- new_x
  }

  # --- Format numeric cells
  fmt_cell <- function(v) {
    if (is.numeric(v) && !is.null(number_format)) {
      vapply(v, function(z) if (is.na(z)) "" else sprintf(number_format, z),
             character(1))
    } else as.character(v)
  }
  fmt_x <- as.data.frame(lapply(x, fmt_cell), stringsAsFactors = FALSE,
                          check.names = FALSE)
  names(fmt_x) <- names(x)

  nrows <- nrow(fmt_x); ncols <- ncol(fmt_x)

  # --- Data layers (cells: row 1..nrows; header: row 0)
  cells <- data.frame(
    row    = rep(seq_len(nrows), ncols),
    col    = rep(seq_len(ncols), each = nrows),
    text   = unlist(fmt_x, use.names = FALSE),
    stripe = rep(ifelse(seq_len(nrows) %% 2 == 0, "even", "odd"), ncols),
    is_rowname_col = rep(if (has_rn) c(TRUE, rep(FALSE, ncols - 1))
                          else rep(FALSE, ncols),
                          each = nrows),
    stringsAsFactors = FALSE)

  header <- data.frame(row = 0, col = seq_len(ncols),
                        text = names(fmt_x),
                        stringsAsFactors = FALSE)

  # --- Diagonal highlight (for square data portion)
  diag_offset <- if (has_rn) 1 else 0
  n_data_rows <- nrows
  n_data_cols <- ncols - diag_offset
  diag_df <- if (isTRUE(highlight_diag) &&
                   n_data_rows == n_data_cols && n_data_rows > 0) {
    data.frame(row = seq_len(n_data_rows),
                col = seq_len(n_data_rows) + diag_offset)
  } else NULL

  # --- Build plot
  p <- ggplot2::ggplot() +
    # Cell backgrounds (alternating stripes)
    ggplot2::geom_tile(data = cells,
                        ggplot2::aes(x = .data$col, y = -.data$row,
                                     fill = .data$stripe),
                        color = "white", linewidth = 0.4,
                        show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = c(odd = "white", even = "grey97"))

  if (!is.null(diag_df)) {
    p <- p + ggplot2::geom_tile(data = diag_df,
                                  ggplot2::aes(x = .data$col, y = -.data$row),
                                  fill = "#D6EAF8", color = "white",
                                  linewidth = 0.4, inherit.aes = FALSE)
  }

  # Header background
  p <- p + ggplot2::geom_tile(data = header,
                                ggplot2::aes(x = .data$col, y = -.data$row),
                                fill = "grey90", color = "white",
                                linewidth = 0.4, inherit.aes = FALSE)

  # --- Top rule, header underline, bottom rule
  rule_top    <- data.frame(x = 0.5, xend = ncols + 0.5,
                              y = 0.5, yend = 0.5)
  rule_under  <- data.frame(x = 0.5, xend = ncols + 0.5,
                              y = -0.5, yend = -0.5)
  rule_bottom <- data.frame(x = 0.5, xend = ncols + 0.5,
                              y = -nrows - 0.5, yend = -nrows - 0.5)
  rules <- rbind(rule_top, rule_bottom)
  p <- p +
    ggplot2::geom_segment(data = rules,
                            ggplot2::aes(x = .data$x, xend = .data$xend,
                                          y = .data$y, yend = .data$yend),
                            color = "black", linewidth = 0.7,
                            inherit.aes = FALSE) +
    ggplot2::geom_segment(data = rule_under,
                            ggplot2::aes(x = .data$x, xend = .data$xend,
                                          y = .data$y, yend = .data$yend),
                            color = "grey50", linewidth = 0.5,
                            inherit.aes = FALSE)

  # --- Cell text
  body_size   <- base_size / 2.85
  header_size <- base_size / 2.7

  # Row-name column text is bold
  if (has_rn) {
    rn_cells <- cells[cells$is_rowname_col, , drop = FALSE]
    body_cells <- cells[!cells$is_rowname_col, , drop = FALSE]
    p <- p +
      ggplot2::geom_text(data = rn_cells,
                          ggplot2::aes(x = .data$col, y = -.data$row,
                                       label = .data$text),
                          size = body_size, fontface = "bold",
                          family = font_family, color = "grey15",
                          inherit.aes = FALSE) +
      ggplot2::geom_text(data = body_cells,
                          ggplot2::aes(x = .data$col, y = -.data$row,
                                       label = .data$text),
                          size = body_size, family = font_family,
                          color = "grey15", inherit.aes = FALSE)
  } else {
    p <- p +
      ggplot2::geom_text(data = cells,
                          ggplot2::aes(x = .data$col, y = -.data$row,
                                       label = .data$text),
                          size = body_size, family = font_family,
                          color = "grey15", inherit.aes = FALSE)
  }

  p <- p +
    ggplot2::geom_text(data = header,
                        ggplot2::aes(x = .data$col, y = -.data$row,
                                     label = .data$text),
                        size = header_size, fontface = "bold",
                        family = font_family, color = "grey5",
                        inherit.aes = FALSE)

  # --- Theme: clean / table-like
  p +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = c(0.5, 0.5))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = c(1.0, 1.0))) +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption) +
    ggplot2::theme_void(base_size = base_size, base_family = font_family) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", hjust = 0,
                                             size = base_size + 1,
                                             margin = ggplot2::margin(b = 4)),
      plot.subtitle = ggplot2::element_text(color = "grey40", hjust = 0,
                                             size = base_size - 1,
                                             margin = ggplot2::margin(b = 8)),
      plot.caption  = ggplot2::element_text(color = "grey40", hjust = 0,
                                             size = base_size - 2,
                                             margin = ggplot2::margin(t = 8)),
      plot.margin   = ggplot2::margin(10, 10, 10, 10)
    )
}


# ---- Multi-modal assembly ---------------------------------------------------

#' Build a validated multi-modal list from separate data frames
#'
#' Takes one data frame per modality and assembles the standard list
#' consumed by [ctdna_merge_modalities()] and the `plot_ctdna_vs_*()` family.
#' Validates required columns based on the current `ctdna_opts()`
#' configuration.
#'
#' Required columns per modality (using the current config):
#'
#' - **ctdna**: subject, time, RECIST, tf
#' - **genomic_74** / **genomic_500**: subject, gene, mutation_status (+ optional panel, alteration_type, vaf)
#' - **expression**: subject, gene, expression
#' - **mutations** (legacy/generic): subject, gene, mutation_status (+ optional panel)
#' - **tmb**: subject, tmb (+ optional panel)
#' - **ihc**: subject, marker, ihc_score
#'
#' Modalities passed as `NULL` are omitted.
#'
#' @param ctdna Longitudinal ctDNA data frame.
#' @param genomic_74 74-gene panel data frame (e.g. GH Infinity).
#' @param genomic_500 500-gene panel data frame (e.g. GH Infinity).
#' @param expression Long-format gene expression data frame.
#' @param mutations Generic / legacy mutation calls (used when neither
#'   `genomic_74` nor `genomic_500` is supplied, or for non-Infinity data).
#' @param tmb Per-subject (per-panel) TMB data.
#' @param ihc Long-format IHC scores.
#' @param strict If `TRUE`, missing required columns raise an error;
#'   otherwise a warning is issued.
#' @param cleanup If \code{TRUE} (default), run light janitorial steps
#'   on each modality (trim whitespace in ID columns, coerce factor
#'   IDs to character, drop empty rows). Set \code{FALSE} to skip.
#' @param filter_scheme Optional scheme name (string) or list passed to
#'   \code{\link{ctdna_variant_filter}} as a pre-filter on the genomic
#'   modalities (\code{genomic_74}, \code{genomic_500}, \code{mutations}).
#'   \code{NULL} (default) leaves the genomic frames untouched.
#' @return Named list with one element per supplied modality.
#' @examples
#' sim <- ctdna_make_mock_study(n_patients = 15, seed = 1)
#' ms  <- ctdna_assemble_modalities(sim)
#' names(ms)  # available modality frames
#' @export
ctdna_assemble_modalities <- function(ctdna = NULL,
                                genomic_74 = NULL,
                                genomic_500 = NULL,
                                expression = NULL,
                                mutations = NULL,
                                tmb = NULL,
                                ihc = NULL,
                                strict = TRUE,
                                cleanup = TRUE,
                                filter_scheme = NULL) {
  # v0.26.0: filter_scheme default changed from "default" to NULL.
  # Filtering is now a separate user step — call ctdna_variant_filter()
  # explicitly on the inputs you want to filter. The argument is retained
  # for backward compatibility: if non-NULL, the filter is still applied,
  # but a one-time deprecation message is emitted.
  required <- list(
    ctdna       = c(.o("subject"), .o("time"), .o("recist"), .o("tf")),
    genomic_74  = c(.o("subject"), .o("gene"), .o("mutation")),
    genomic_500 = c(.o("subject"), .o("gene"), .o("mutation")),
    expression  = c(.o("subject"), .o("gene"), .o("expression")),
    mutations   = c(.o("subject"), .o("gene"), .o("mutation")),
    tmb         = c(.o("subject"), .o("tmb")),
    ihc         = c(.o("subject"), .o("marker"), .o("ihc_score"))
  )
  inputs <- list(ctdna = ctdna,
                 genomic_74 = genomic_74, genomic_500 = genomic_500,
                 expression = expression, mutations = mutations,
                 tmb = tmb, ihc = ihc)
  inputs <- inputs[!vapply(inputs, is.null, logical(1))]
  if (length(inputs) == 0) stop("Pass at least one modality.")

  # v0.22.1: harmonize TF-domain column names BEFORE validation. Any
  # variant spelling of methylTF / maxVAF / meanVAF in the ctdna frame
  # is renamed to the canonical name. Users can name their input
  # columns however they want; downstream code sees only the canonical
  # spelling.
  if (!is.null(inputs$ctdna))
    inputs$ctdna <- .ctdna_harmonize_tf_columns(inputs$ctdna, quiet = FALSE)

  # Canonicalize expression unit to ctdna_opts("expression_unit") — default
  # log2(TPM + 1). If the input already carries an attr 'units' = 'tpm', it
  # gets converted in-place. If 'units' is missing, we auto-detect TPM vs
  # log2(TPM+1) by value range and warn (so the user is never confused).
  if ("expression" %in% names(inputs)) {
    inputs$expression <- .canonicalize_expression(inputs$expression)
  }

  for (m in names(inputs)) {
    df <- inputs[[m]]
    if (!is.data.frame(df))
      stop(sprintf("'%s' must be a data frame.", m))
    miss <- setdiff(required[[m]], names(df))
    if (length(miss)) {
      msg <- sprintf("'%s' is missing required column(s): %s.\n  Either rename in your data or set with ctdna_opts(...).",
                     m, paste(miss, collapse = ", "))
      if (strict) stop(msg) else warning(msg)
    }
  }

  # Cleanup is sample-level QC and stays here (it's not "filtering" in
  # the alteration-level sense). Alteration filtering moved OUT of this
  # function in v0.26.0 — call ctdna_variant_filter() explicitly afterward.
  if (isTRUE(cleanup)) {
    for (m in names(inputs))
      inputs[[m]] <- ctdna_cleanup(inputs[[m]], quiet = TRUE)
  }
  if (!is.null(filter_scheme)) {
    .ctdna_warn_once("assemble_modalities_filter_scheme", paste0(
      "v0.26.0: ctdna_assemble_modalities() no longer auto-filters by default ",
      "(filtering is now a separate user step). The `filter_scheme` argument ",
      "is retained for backward compatibility but will be removed in a ",
      "future version. Recommended: leave `filter_scheme = NULL` here and ",
      "call ctdna_variant_filter() explicitly on the modalities you need."))
    for (m in c("genomic_74","genomic_500","mutations")) {
      if (m %in% names(inputs))
        inputs[[m]] <- .filter_apply_df(inputs[[m]],
                                            filter_scheme = filter_scheme)
    }
  }
  inputs
}


# ---- Expression-unit canonicalization ---------------------------------------
# Canonical unit: ctdna_opts("expression_unit") — default "log2_tpm_plus_one".
# Two flows are supported:
#   1) The data frame carries attr(df, "units") = "tpm" or "log2_tpm_plus_one".
#      We honour the tag (convert if needed; tag the output).
#   2) No `units` attribute. We auto-detect from the value range:
#        - max value > 30  → almost certainly TPM (log2(TPM+1) rarely exceeds ~25)
#        - else            → assume already log2(TPM+1)
#      A message is printed so the user is never silently confused.

# Distribution-based TPM vs log2(TPM+1) detection.
# Old heuristic used max(values) > 30 — a single huge housekeeping gene
# could fool it. New rule: sample up to 10,000 non-NA, non-zero values
# across genes x samples; use the 99th percentile.
#
# Reasoning:
#  - TPM is heavily right-skewed: median ~1-5, 99th percentile in the
#    hundreds to thousands.
#  - log2(TPM+1) is bounded ~0-15 in practice: 99th percentile rarely > 12.
#  - 50 is a comfortable margin between the two.
.detect_expression_unit <- function(values, sample_size = 10000L) {
  x <- as.numeric(values)
  x <- x[!is.na(x) & x > 0]
  if (length(x) == 0) {
    message("expression: no non-zero values found; defaulting to log2(TPM+1).")
    return("log2_tpm_plus_one")
  }
  if (length(x) > sample_size)
    x <- sample(x, sample_size)
  p99 <- stats::quantile(x, 0.99, na.rm = TRUE, names = FALSE)
  med <- stats::median(x, na.rm = TRUE)
  unit <- if (p99 > 50) "tpm" else "log2_tpm_plus_one"
  message(sprintf(
    "expression: detected '%s' (median=%.2f, 99th pct=%.2f, n=%d sampled). Set attr(df, 'units') explicitly to skip auto-detection.",
    unit, med, p99, length(x)))
  unit
}
.canonicalize_expression <- function(df) {
  target <- .o("expression_unit")
  expr_col <- .o("expression")
  if (!expr_col %in% names(df)) return(df)
  src <- attr(df, "units")
  if (is.null(src)) {
    # The new .detect_expression_unit already prints a one-line summary.
    src <- .detect_expression_unit(df[[expr_col]])
  }
  if (identical(src, target)) {
    attr(df, "units") <- target
    return(df)
  }
  if (src == "tpm" && target == "log2_tpm_plus_one") {
    df[[expr_col]] <- log2(pmax(as.numeric(df[[expr_col]]), 0) + 1)
    message("expression: converted TPM -> log2(TPM + 1)")
  } else if (src == "log2_tpm_plus_one" && target == "tpm") {
    df[[expr_col]] <- 2 ^ as.numeric(df[[expr_col]]) - 1
    message("expression: converted log2(TPM + 1) -> TPM")
  } else {
    warning(sprintf("expression: don't know how to convert '%s' to '%s'; leaving values unchanged.", src, target))
  }
  attr(df, "units") <- target
  df
}

#' Canonicalize gene expression values to the package unit
#'
#' Converts a long-format expression frame between \code{TPM} and
#' \code{log2(TPM + 1)} so every plotting / statistical function sees
#' the same unit. The target unit is taken from
#' \code{ctdna_opts("expression_unit")} (default
#' \code{"log2_tpm_plus_one"}).
#'
#' The function honours \code{attr(df, "units")} when set (recommended);
#' otherwise it auto-detects from the value range (max > 30 → assumed
#' TPM) and emits a message so the assumption is visible.
#'
#' @param df Long-format expression frame (\code{subject}, \code{gene},
#'   \code{expression} per current \code{ctdna_opts()}).
#' @return The same frame with values in the canonical unit and
#'   \code{attr(.,"units")} set accordingly.
#' @examples
#' # Build a tiny TPM frame and convert it
#' df <- data.frame(subject_id = "S001", gene = c("A","B","C"),
#'                  expression = c(0, 5.4, 412))
#' attr(df, "units") <- "tpm"
#' ctdna_to_log2tpm(df)
#' @export
ctdna_to_log2tpm <- function(df) {
  .canonicalize_expression(df)
}

# ============================================================================
# ctdnaTM — ctdna_make_mock_study()
# ============================================================================
# Comprehensive study simulator. Produces realistic mock data for every
# data type a ctDNA translational-medicine analysis touches, using the
# RAW VENDOR COLUMN NAMES as delivered by Guardant Health Infinity (plus
# IHC, tissue WES/WGS alteration calls, and RNA-seq count + TPM matrices).
#
# The 8 possible outputs are (each gated by its modality flag in v0.40.0):
#   $infinity_report      — File 1: per-sample per-alteration master report
#   $tf_change            — File 2: paired TF change, MR_panel = "TF.methyl"
#   $panel_74_response    — File 3: paired 74-gene response report
#   $panel_500_response   — File 4: paired 500-gene response report
#   $methylation_response — File 5: paired methylation panel response
#                                    (MR_panel = "legacy.methyl")
#   $ihc                  — per-patient IHC marker scores
#   $dnaseq               — tissue DNA-seq (WES/WGS) alteration calls
#                           (v0.40.0; previously $wes)
#   $rnaseq               — list($counts, $tpm, $sample_metadata)
#
# Conventions:
#   - All patients share the same Study_ID
#   - Each (Patient_ID, Visit_name) has a unique GHSampleID
#   - Visit A is always C1D1 (baseline); Visit B is each later visit
#   - Cancertype is fixed per patient (the indication)
#   - Methylation TF stored in PERCENTAGE units (matches vendor column name)
#   - Dates: real Date objects with realistic offsets per visit cycle
# ============================================================================


# ---- Internal helpers --------------------------------------------------------

# Gene panels (representative subsets — real Infinity panels are larger)
.gh_genes_74 <- c("TP53","KRAS","EGFR","BRAF","PIK3CA","ERBB2","MET","RET",
                  "ROS1","ALK","NRAS","PTEN","CDKN2A","FGFR1","FGFR3",
                  "IDH1","IDH2","KIT","CCND1","CDK4")
.gh_genes_500_only <- c("BRCA1","BRCA2","ATM","CHEK2","RB1","NF1","APC",
                        "MSH2","MLH1","STK11","PALB2","BARD1","BAP1","BLM",
                        "MYC","CTNNB1","SMAD4","ARID1A","SETD2","KMT2D")
.gh_fusion_partners <- list(ALK = "EML4", ROS1 = "CD74", RET = "KIF5B",
                            NTRK1 = "TPM3", FGFR3 = "TACC3", BRAF = "AGK")

# Build a single alteration's variant-specific fields (NA for fields that
# don't apply to that variant_type). Returns a named list of 21 fields.
.gh_build_variant <- function(variant_type, sample_tf_pct) {
  g <- sample(.gh_genes_74, 1)
  row <- list(
    Variant_type          = variant_type,
    Indel_type            = NA_character_,
    Gene                  = g,
    Chromosome            = paste0("chr", sample(1:22, 1)),
    Position              = sample(1e6:2e8, 1),
    Exon                  = NA_integer_,
    Mut_aa                = NA_character_,
    Mut_nt                = NA_character_,
    Mut_cdna              = NA_character_,
    Transcript            = NA_character_,
    VAF_percentage        = NA_real_,
    Splice_effect         = NA_character_,
    Somatic_status        = NA_character_,
    Molecular_consequence = NA_character_,
    Fusion_chrom_b        = NA_character_,
    Fusion_gene_b         = NA_character_,
    Fusion_position_a     = NA_real_,
    Fusion_position_b     = NA_real_,
    Direction_a           = NA_character_,
    Direction_b           = NA_character_,
    Downstream_gene       = NA_character_,
    Copy_number           = NA_real_,
    CNV_type              = NA_character_,
    COSMIC                = NA_character_,
    dbSNP                 = NA_character_,
    ClinVar               = NA_character_,
    ClinVarID             = NA_character_,
    Functional_impact     = NA_character_,
    Mutant_allele_status  = NA_character_,
    Mol_count             = NA_integer_,
    Alleletype            = NA_character_
  )

  aa <- c("A","R","N","D","C","E","Q","G","H","I","L","K","M","F",
          "P","S","T","W","Y","V")
  nt <- c("A","C","G","T")

  if (variant_type %in% c("SNV", "Indel")) {
    row$Exon       <- sample(1:25, 1)
    row$Transcript <- sprintf("NM_%06d.%d", sample(1000:9999, 1),
                              sample(1:5, 1))
    if (variant_type == "SNV") {
      row$Mut_aa <- sprintf("p.%s%d%s", sample(aa, 1),
                            sample(50:600, 1), sample(aa, 1))
      # v0.24.0: Mut_nt as explicit ref>alt (Guardant Infinity convention)
      .ref_b <- sample(nt, 1)
      .alt_b <- sample(setdiff(nt, .ref_b), 1)
      row$Mut_nt <- sprintf("%s>%s", .ref_b, .alt_b)
      row$Molecular_consequence <- sample(
        c("missense","nonsense","synonymous","splice_donor","splice_acceptor"),
        1, prob = c(0.55, 0.15, 0.10, 0.10, 0.10))
      row$Splice_effect <- if (grepl("splice", row$Molecular_consequence))
        sample(c("Donor lost","Acceptor lost","Cryptic donor"), 1) else NA
    } else {
      row$Indel_type <- sample(c("Insertion", "Deletion"), 1)
      row$Mut_aa <- sprintf("p.%s%dfs", sample(aa, 1), sample(50:600, 1))
      # v0.24.0: Mut_nt for indels also uses explicit ref>alt with
      # a left-anchoring base, the VCF-normalized form gnomAD wants.
      .anchor <- sample(nt, 1)
      .ins_seq <- paste(sample(nt, sample(2:6, 1), replace = TRUE),
                         collapse = "")
      row$Mut_nt <- if (row$Indel_type == "Deletion")
                       sprintf("%s%s>%s", .anchor, .ins_seq, .anchor)
                    else
                       sprintf("%s>%s%s", .anchor, .anchor, .ins_seq)
      row$Molecular_consequence <- "frameshift"
    }
    row$Mut_cdna   <- row$Mut_nt
    # VAF roughly tracks sample TF, with skewed distribution
    base_vaf <- stats::rbeta(1, 2, 50) * 100
    row$VAF_percentage <- round(min(base_vaf + sample_tf_pct / 4, 50), 3)
    row$Somatic_status <- sample(c("somatic","germline","somatic_putative_ch"),
                                  1, prob = c(0.80, 0.15, 0.05))
    row$COSMIC    <- if (stats::runif(1) < 0.7)
      sprintf("COSM%d", sample(1e4:9.99e5, 1)) else NA
    row$dbSNP     <- if (stats::runif(1) < 0.4)
      sprintf("rs%d", sample(1e5:9e8, 1)) else NA
    row$ClinVar   <- sample(c(NA, "Pathogenic", "Likely_pathogenic",
                              "Uncertain_significance","Benign","Likely_benign"),
                            1, prob = c(0.45, 0.18, 0.12, 0.15, 0.05, 0.05))
    row$ClinVarID <- if (!is.na(row$ClinVar))
      sprintf("VCV%07d", sample(1e5:1e7, 1)) else NA
    row$Functional_impact <- sample(c("deleterious",NA),
                                     1, prob = c(0.45, 0.55))
    row$Mutant_allele_status <- sample(c(NA,"biallelic"),
                                        1, prob = c(0.85, 0.15))
    row$Mol_count   <- sample(20:5000, 1)
    row$Alleletype  <- "Somatic"
  } else if (variant_type == "CNV") {
    row$CNV_type    <- sample(c("amplification","focal_amplification","aneuploid_amplification",
                                            "homozygous_deletion","loh_deletion"),
                              1, prob = c(0.25, 0.25, 0.15, 0.20, 0.15))
    row$Copy_number <- if (grepl("amplification", row$CNV_type))
      sample(4:25, 1) else sample(0:1, 1)
    row$Functional_impact <- "deleterious"
    row$Somatic_status <- "somatic"
    row$Alleletype  <- "Somatic"
  } else if (variant_type == "Fusion") {
    partner <- if (g %in% names(.gh_fusion_partners))
      .gh_fusion_partners[[g]] else sample(.gh_genes_74, 1)
    row$Fusion_gene_b     <- partner
    row$Fusion_chrom_b    <- paste0("chr", sample(1:22, 1))
    row$Fusion_position_a <- sample(1e6:2e8, 1)
    row$Fusion_position_b <- sample(1e6:2e8, 1)
    row$Direction_a       <- sample(c("+","-"), 1)
    row$Direction_b       <- sample(c("+","-"), 1)
    row$Downstream_gene   <- partner
    row$Functional_impact <- "deleterious"
    row$Somatic_status    <- "somatic"
  } else if (variant_type == "LGR") {
    row$Functional_impact <- "deleterious"
    row$Somatic_status    <- "somatic"
  }
  row
}


# ---- Main API ---------------------------------------------------------------

#' Simulate a complete ctDNA translational-medicine study
#'
#' Produces realistic mock data for every data type a ctDNA analysis
#' touches, using the \strong{raw vendor column names} as delivered by
#' Guardant Health Infinity, plus IHC, tissue DNA-seq (WES/WGS)
#' alteration calls, and RNA-seq (count + TPM matrices with sample
#' metadata).
#'
#' Modality flags let you generate only the modalities you need. The
#' returned list contains exactly the modalities you asked for — there
#' are no \code{NULL} placeholders for skipped modalities.
#'
#' Possible outputs:
#' \itemize{
#'   \item \strong{ctdna = TRUE} (default) — adds \code{infinity_report},
#'     \code{tf_change}, \code{panel_74_response},
#'     \code{panel_500_response}, \code{methylation_response}
#'   \item \strong{clinical = TRUE} (default) — adds \code{clinical}
#'     (per-patient RECIST / Dose / Cancertype / best_pct_change / MR)
#'   \item \strong{ihc = TRUE} (default) — adds \code{ihc} (per-patient
#'     IHC marker scores)
#'   \item \strong{dnaseq = TRUE} (default) — adds \code{dnaseq} (tissue
#'     WES/WGS alteration calls)
#'   \item \strong{rnaseq = TRUE} (default) — adds \code{rnaseq} (a list
#'     of \code{counts}, \code{tpm}, and \code{sample_metadata})
#' }
#'
#' Patients, visits, sample IDs, dates and cancer types are consistent
#' across all outputs (a single patient appears in all enabled outputs,
#' joinable by \code{Patient_ID}).
#'
#' @param n_patients Number of patients. Default 30.
#' @param visits Visit labels in chronological order. The first must be
#'   the baseline (e.g., \code{"C1D1"}).
#' @param cancer_types Indications to sample from.
#' @param doses Dose levels to sample from.
#' @param n_rnaseq_genes Number of genes in the RNA-seq matrices. Used
#'   only when \code{rnaseq = TRUE}.
#' @param study_id Study identifier used as the prefix throughout.
#' @param seed Random seed for reproducibility. Default 42.
#' @param ctdna Logical. If \code{TRUE} (default), include the GH
#'   Infinity ctDNA bundle (\code{infinity_report}, \code{tf_change},
#'   \code{panel_74_response}, \code{panel_500_response},
#'   \code{methylation_response}).
#' @param clinical Logical. If \code{TRUE} (default), include the
#'   per-patient \code{clinical} frame.
#' @param ihc Logical. If \code{TRUE} (default), include the \code{ihc}
#'   frame.
#' @param dnaseq Logical. If \code{TRUE} (default), include the
#'   \code{dnaseq} (tissue WES/WGS alteration calls) frame.
#' @param rnaseq Logical. If \code{TRUE} (default), include the
#'   \code{rnaseq} list.
#' @return Named list containing only the requested modalities. The
#'   list carries \code{platform}, \code{study_id} and
#'   \code{n_patients} attributes.
#' @examples
#' sim <- ctdna_make_mock_study(n_patients = 15, seed = 1)
#' names(sim)   # all modality frames + clinical
#' head(sim$infinity_report)
#' head(sim$tf_change)
#' @export
ctdna_make_mock_study <- function(n_patients   = 30,
                            visits         = c("C1D1","C2D1","C3D1","C5D1","EOT"),
                            cancer_types   = c("NSCLC","CRC","BRCA","mCRPC","GBM"),
                            doses          = c("Low","Mid","High"),
                            n_rnaseq_genes = 2000,
                            study_id       = "STUDY01",
                            seed           = 42,
                            ctdna          = TRUE,
                            clinical       = TRUE,
                            ihc            = TRUE,
                            dnaseq         = TRUE,
                            rnaseq         = TRUE) {
  set.seed(seed)

  # v0.40.0: capture flag values before any local var named the same
  # gets bound to a constructed data frame.
  .keep_ctdna    <- isTRUE(ctdna)
  .keep_clinical <- isTRUE(clinical)
  .keep_ihc      <- isTRUE(ihc)
  .keep_dnaseq   <- isTRUE(dnaseq)
  .keep_rnaseq   <- isTRUE(rnaseq)

  # ============================================================
  # Patient master table (internal — has biology fields with .)
  # ============================================================
  rec_levels <- c("CR","PR","uCR","SD","PD","NE")
  pts <- data.frame(
    Patient_ID = sprintf("PT%03d", seq_len(n_patients)),
    Cancertype = sample(cancer_types, n_patients, replace = TRUE),
    Dose       = sample(doses, n_patients, replace = TRUE),
    RECIST     = sample(rec_levels, n_patients, replace = TRUE,
                        prob = c(0.10, 0.20, 0.10, 0.30, 0.25, 0.05)),
    stringsAsFactors = FALSE
  )
  best_eff <- c(CR=-100, PR=-50, uCR=-25, SD=-10, PD=25, NE=0)
  pts$best_pct_change <- best_eff[pts$RECIST] +
    stats::rnorm(n_patients, 0, 10)
  # MR derived from RECIST (CR/PR/uCR → Response; SD/PD/NE → NonResponse)
  pts$MR <- ifelse(pts$RECIST %in% c("CR","PR","uCR"),
                   "Response", "NonResponse")

  # Internal biology fields (prefixed . for clarity)
  rec_drop <- c(CR=1.5, PR=1.0, uCR=0.7, SD=0.3, PD=-0.5, NE=0)
  pts$.drop_per_cycle <- rec_drop[pts$RECIST]
  pts$.tf_base_pct    <- 10 ^ stats::runif(n_patients, -2, 0.5)  # 0.01-3 %
  pts$.sld_base       <- round(stats::runif(n_patients, 40, 120), 1)  # baseline SLD (mm)
  pts$.tmb_base       <- pmax(stats::rgamma(n_patients, 2, 0.5), 0.1)
  pts$.msi_high       <- sample(c("Yes","No"), n_patients,
                                replace = TRUE, prob = c(0.05, 0.95))
  pts$.hrd            <- ifelse(stats::runif(n_patients) < 0.15,
                                round(stats::runif(n_patients, 30, 80), 1),
                                NA_real_)
  pts$.baseline_date  <- as.Date("2024-01-01") +
    sample(0:365, n_patients, replace = TRUE)

  # ============================================================
  # Sample master table (one row per Patient × Visit)
  # ============================================================
  visit_days <- c(C1D1=0, C2D1=21, C3D1=42, C5D1=84, EOT=168)
  sample_rows <- list()
  sidx <- 1
  for (i in seq_len(nrow(pts))) {
    pat <- pts[i, ]
    for (v in visits) {
      vidx <- which(visits == v) - 1
      day_off <- as.integer(visit_days[[v]] + sample(-2:2, 1))
      bloodcoll <- pat$.baseline_date + day_off

      if (vidx == 0) {
        tf <- pat$.tf_base_pct
      } else {
        log10_red <- pat$.drop_per_cycle * vidx + stats::rnorm(1, 0, 0.3)
        tf <- max(pat$.tf_base_pct * 10 ^ (-log10_red), 0.001)
      }
      # Biological cap: methylation TF rarely exceeds 50-80% even in
      # advanced disease. Clip to a realistic ceiling.
      tf <- min(tf, 80)
      qc <- sample(c("SUCCESS","SUCCESS","SUCCESS","SUCCESS","FAIL"), 1)
      sample_rows[[sidx]] <- data.frame(
        sample_seq          = sidx,
        Study_ID            = study_id,
        Patient_ID          = pat$Patient_ID,
        Visit_name          = v,
        Cancertype          = pat$Cancertype,
        Sample_status       = qc,
        Sample_comment      = if (qc == "SUCCESS") NA_character_
                              else "Sample below input requirements",
        Customer_SampleID   = sprintf("%s_%s_%s", study_id, pat$Patient_ID, v),
        GHSampleID          = sprintf("GHS%07d", 5000000 + sidx),
        GHRequestID         = sprintf("GHR%05d", 10000 + sidx),
        Bloodcoll_date      = bloodcoll,
        Received_date       = bloodcoll + sample(1:3, 1),
        Reported_date       = bloodcoll + sample(7:14, 1),
        TF_pct              = tf,
        Lesion_sum_mm       = {
          nvis <- length(visits); frac <- if (nvis > 1) vidx / (nvis - 1) else 0
          round(max(pat$.sld_base * (1 + (pat$best_pct_change/100) * frac) + stats::rnorm(1, 0, 2), 1), 1)
        },
        TMB                 = max(pat$.tmb_base + stats::rnorm(1, 0, 0.5), 0.1),
        MSI_High            = pat$.msi_high,
        HRD_score           = pat$.hrd,
        cfDNA_ng            = round(stats::rgamma(1, 2, 0.5) * 20, 1),
        Plasma_ml_input     = sample(c(4, 6, 8, 10), 1),
        Plasma_ml_remaining = round(stats::runif(1, 0, 5), 1),
        stringsAsFactors    = FALSE
      )
      sidx <- sidx + 1
    }
  }
  samples <- do.call(rbind, sample_rows)
  samples$ctDNA_detection_status <- ifelse(
    samples$TF_pct < 0.1, "Not Detected",
    ifelse(samples$TF_pct < 1, "Below LOQ", "Detected"))
  samples$TMB_category <- ifelse(
    samples$TMB < 6, "Low",
    ifelse(samples$TMB >= 10, "High", "Intermediate"))

  # ============================================================
  # File 1: InfinityReport (per-sample, per-alteration rows)
  # ============================================================
  infinity_report <- tf_change <- panel_74_response <-
    panel_500_response <- methylation_response <- NULL
  if (.keep_ctdna) {
  inf_rows <- list()
  ridx <- 1
  for (si in seq_len(nrow(samples))) {
    sr <- samples[si, ]

    common_id <- list(
      Study_ID          = sr$Study_ID,
      Customer_SampleID = sr$Customer_SampleID,
      GHRequestID       = sr$GHRequestID,
      GHSampleID        = sr$GHSampleID,
      Patient_ID        = sr$Patient_ID,
      Visit_name        = sr$Visit_name
    )
    common_bio <- list(
      Genomic_max_VAF_percentage            = NA_real_,
      HRD_score                             = sr$HRD_score,
      ctDNA_detection_status                = sr$ctDNA_detection_status,
      Methylation_tumor_fraction_percentage = sr$TF_pct,
      TMB_score                             = round(sr$TMB, 2),
      TMB_category                          = sr$TMB_category,
      MSI_High                              = sr$MSI_High,
      cfDNA_ng                              = sr$cfDNA_ng,
      Plasma_ml_input                       = sr$Plasma_ml_input,
      Plasma_ml_remaining                   = sr$Plasma_ml_remaining,
      Received_date                         = sr$Received_date,
      Bloodcoll_date                        = sr$Bloodcoll_date,
      Reported_date                         = sr$Reported_date,
      Cancertype                            = sr$Cancertype,
      Practice_name                         = "Sample Oncology Practice",
      Physician_name                        = "Dr. Sample"
    )
    na_variant <- .gh_build_variant("SNV", 0)  # template
    na_variant[] <- lapply(na_variant, function(x) x[NA])

    if (sr$Sample_status != "SUCCESS") {
      inf_rows[[ridx]] <- as.data.frame(c(
        common_id,
        list(Genomic_somatic_alteration_detected = NA_character_,
             Sample_status  = sr$Sample_status,
             Sample_comment = sr$Sample_comment),
        na_variant,
        lapply(common_bio, function(x) if (length(x)) x else NA)
      ), stringsAsFactors = FALSE)
      ridx <- ridx + 1
      next
    }

    expected_n <- 1 + (sr$TF_pct > 0.5) * 2 + (sr$TF_pct > 5) * 3
    n_alt <- stats::rpois(1, expected_n)

    if (n_alt == 0) {
      inf_rows[[ridx]] <- as.data.frame(c(
        common_id,
        list(Genomic_somatic_alteration_detected = "No",
             Sample_status  = sr$Sample_status,
             Sample_comment = NA_character_),
        na_variant,
        common_bio
      ), stringsAsFactors = FALSE)
      ridx <- ridx + 1
      next
    }

    alt_types <- sample(c("SNV","Indel","CNV","Fusion","LGR"),
                        n_alt, replace = TRUE,
                        prob = c(0.60, 0.10, 0.15, 0.10, 0.05))
    for (alt_t in alt_types) {
      variant <- .gh_build_variant(alt_t, sr$TF_pct)
      inf_rows[[ridx]] <- as.data.frame(c(
        common_id,
        list(Genomic_somatic_alteration_detected = "Yes",
             Sample_status  = sr$Sample_status,
             Sample_comment = NA_character_),
        variant,
        common_bio
      ), stringsAsFactors = FALSE)
      ridx <- ridx + 1
    }
  }
  infinity_report <- do.call(rbind, inf_rows)
  rownames(infinity_report) <- NULL

  # Fill in Genomic_max_VAF_percentage per sample
  for (gh in unique(infinity_report$GHSampleID)) {
    rows <- which(infinity_report$GHSampleID == gh)
    vafs <- infinity_report$VAF_percentage[rows]
    if (any(!is.na(vafs)))
      infinity_report$Genomic_max_VAF_percentage[rows] <-
        max(vafs, na.rm = TRUE)
  }

  # Enforce exact column order from spec
  inf_col_order <- c(
    "Study_ID","Customer_SampleID","GHRequestID","GHSampleID","Patient_ID",
    "Visit_name","Genomic_somatic_alteration_detected","Sample_status",
    "Sample_comment","Variant_type","Indel_type","Gene","Chromosome","Position",
    "Exon","Mut_aa","Mut_nt","Mut_cdna","Transcript","VAF_percentage",
    "Splice_effect","Somatic_status","Molecular_consequence","Fusion_chrom_b",
    "Fusion_gene_b","Fusion_position_a","Fusion_position_b","Direction_a",
    "Direction_b","Downstream_gene","Copy_number","CNV_type","COSMIC","dbSNP",
    "ClinVar","ClinVarID","Functional_impact","Mutant_allele_status",
    "Mol_count","Genomic_max_VAF_percentage","Alleletype","HRD_score",
    "ctDNA_detection_status","Methylation_tumor_fraction_percentage",
    "TMB_score","TMB_category","MSI_High","cfDNA_ng","Plasma_ml_input",
    "Plasma_ml_remaining","Received_date","Bloodcoll_date","Reported_date",
    "Cancertype","Practice_name","Physician_name")
  infinity_report <- infinity_report[, inf_col_order]

  # ============================================================
  # Build paired metadata (Visit A = C1D1, Visit B = each later visit)
  # ============================================================
  pair_list <- list()
  pidx <- 1
  for (i in seq_len(nrow(pts))) {
    pid <- pts$Patient_ID[i]
    a <- samples[samples$Patient_ID == pid & samples$Visit_name == "C1D1", ]
    if (nrow(a) == 0 || a$Sample_status != "SUCCESS") next
    for (v in setdiff(visits, "C1D1")) {
      b <- samples[samples$Patient_ID == pid & samples$Visit_name == v, ]
      if (nrow(b) == 0 || b$Sample_status != "SUCCESS") next
      pair_list[[pidx]] <- data.frame(
        Study_ID            = study_id,
        Patient_ID          = pid,
        Visit_name_A        = "C1D1",
        Customer_SampleID_A = a$Customer_SampleID,
        GHSampleID_A        = a$GHSampleID,
        GHRequestID_A       = a$GHRequestID,
        Cancertype_A        = a$Cancertype,
        Sample_status_A     = a$Sample_status,
        Sample_comment_A    = a$Sample_comment,
        Bloodcoll_date_A    = a$Bloodcoll_date,
        Visit_name_B        = v,
        Customer_SampleID_B = b$Customer_SampleID,
        GHSampleID_B        = b$GHSampleID,
        GHRequestID_B       = b$GHRequestID,
        Cancertype_B        = b$Cancertype,
        Sample_status_B     = b$Sample_status,
        Sample_comment_B    = b$Sample_comment,
        Bloodcoll_date_B    = b$Bloodcoll_date,
        .tf_a               = a$TF_pct,
        .tf_b               = b$TF_pct,
        .vaf_a              = NA_real_,
        .vaf_b              = NA_real_,
        stringsAsFactors    = FALSE
      )
      pidx <- pidx + 1
    }
  }
  paired <- do.call(rbind, pair_list)

  # Per-pair derived stats (shared)
  paired$days_between <- as.integer(paired$Bloodcoll_date_B -
                                     paired$Bloodcoll_date_A)
  paired$Ctdna_pct_change <- (paired$.tf_b - paired$.tf_a) /
    pmax(paired$.tf_a, 0.001) * 100
  paired$Ctdna_lvl_change <- ifelse(paired$.tf_b < paired$.tf_a * 0.5,
                                     "Decrease",
                                     ifelse(paired$.tf_b > paired$.tf_a * 2,
                                            "Increase","No change"))
  paired$MR_QC_status <- sample(c("PASS","PASS","PASS","FAIL"),
                                 nrow(paired), replace = TRUE)
  paired$MR_quantifiable <- sample(c("Yes","No"), nrow(paired),
                                    replace = TRUE, prob = c(0.85, 0.15))
  # Per-panel mean/max VAF based on TF
  paired$.meanvaf_a <- pmax(paired$.tf_a / 3 +
                              stats::rnorm(nrow(paired), 0, 0.5), 0.01)
  paired$.meanvaf_b <- pmax(paired$.tf_b / 3 +
                              stats::rnorm(nrow(paired), 0, 0.5), 0.01)
  paired$.maxvaf_a  <- paired$.meanvaf_a * stats::runif(nrow(paired), 1.3, 2.0)
  paired$.maxvaf_b  <- paired$.meanvaf_b * stats::runif(nrow(paired), 1.3, 2.0)
  paired$.mr_variant_count <- pmax(round(paired$.tf_a / 2), 0L)

  # ============================================================
  # File 2: tf_change (MR_panel = "TF.methyl")
  # ============================================================
  tf_change <- data.frame(
    Study_ID            = paired$Study_ID,
    Patient_ID          = paired$Patient_ID,
    Visit_name_A        = paired$Visit_name_A,
    Customer_SampleID_A = paired$Customer_SampleID_A,
    GHSampleID_A        = paired$GHSampleID_A,
    GHRequestID_A       = paired$GHRequestID_A,
    Cancertype_A        = paired$Cancertype_A,
    Sample_status_A     = paired$Sample_status_A,
    Sample_comment_A    = paired$Sample_comment_A,
    Bloodcoll_date_A    = paired$Bloodcoll_date_A,
    Visit_name_B        = paired$Visit_name_B,
    Customer_SampleID_B = paired$Customer_SampleID_B,
    GHSampleID_B        = paired$GHSampleID_B,
    GHRequestID_B       = paired$GHRequestID_B,
    Cancertype_B        = paired$Cancertype_B,
    Sample_status_B     = paired$Sample_status_B,
    Sample_comment_B    = paired$Sample_comment_B,
    Bloodcoll_date_B    = paired$Bloodcoll_date_B,
    Days_between        = paired$days_between,
    MR_panel            = "TF.methyl",
    MR_QC_status        = paired$MR_QC_status,
    MR_quantifiable     = paired$MR_quantifiable,
    MR_score_percentage = round(paired$.tf_b, 4),
    Ctdna_percentage_change = round(paired$Ctdna_pct_change, 2),
    Ctdna_level_change      = paired$Ctdna_lvl_change,
    Methylation_tumor_fraction_percentage_A = round(paired$.tf_a, 4),
    Methylation_tumor_fraction_percentage_B = round(paired$.tf_b, 4),
    stringsAsFactors    = FALSE
  )

  # ============================================================
  # File 3 & 4: panel_74_response / panel_500_response
  #   (Shared schema with MR_variant_count + mean/max VAF columns;
  #   500-gene generally has more variants → slightly different stats.)
  # ============================================================
  build_panel_response <- function(panel_label, vaf_scale) {
    data.frame(
      Study_ID            = paired$Study_ID,
      Patient_ID          = paired$Patient_ID,
      Visit_name_A        = paired$Visit_name_A,
      Customer_SampleID_A = paired$Customer_SampleID_A,
      GHSampleID_A        = paired$GHSampleID_A,
      GHRequestID_A       = paired$GHRequestID_A,
      Cancertype_A        = paired$Cancertype_A,
      Sample_status_A     = paired$Sample_status_A,
      Sample_comment_A    = paired$Sample_comment_A,
      Visit_name_B        = paired$Visit_name_B,
      Customer_SampleID_B = paired$Customer_SampleID_B,
      GHSampleID_B        = paired$GHSampleID_B,
      GHRequestID_B       = paired$GHRequestID_B,
      Cancertype_B        = paired$Cancertype_B,
      Sample_status_B     = paired$Sample_status_B,
      Sample_comment_B    = paired$Sample_comment_B,
      Days_between_samples = paired$days_between,
      MR_panel            = panel_label,
      MR_variant_count    = if (panel_label == "74genes")
                              paired$.mr_variant_count
                            else paired$.mr_variant_count + sample(0:3,
                              nrow(paired), replace = TRUE),
      MR_QC_status        = paired$MR_QC_status,
      Mean_VAF_percentage_A = round(paired$.meanvaf_a * vaf_scale, 4),
      Mean_VAF_percentage_B = round(paired$.meanvaf_b * vaf_scale, 4),
      Max_VAF_percentage_A  = round(paired$.maxvaf_a  * vaf_scale, 4),
      Max_VAF_percentage_B  = round(paired$.maxvaf_b  * vaf_scale, 4),
      MR_quantifiable     = paired$MR_quantifiable,
      MR_score_percentage = round(paired$.meanvaf_b * vaf_scale, 4),
      Ctdna_percentage_change = round(paired$Ctdna_pct_change, 2),
      Ctdna_level_change      = paired$Ctdna_lvl_change,
      stringsAsFactors    = FALSE
    )
  }
  panel_74_response  <- build_panel_response("74genes",  1.0)
  panel_500_response <- build_panel_response("500genes", 1.05)

  # ============================================================
  # File 5: methylation_response (MR_panel = "legacy.methyl")
  # ============================================================
  methylation_response <- data.frame(
    Study_ID            = paired$Study_ID,
    Patient_ID          = paired$Patient_ID,
    Visit_name_A        = paired$Visit_name_A,
    Customer_SampleID_A = paired$Customer_SampleID_A,
    GHSampleID_A        = paired$GHSampleID_A,
    GHRequestID_A       = paired$GHRequestID_A,
    Cancertype_A        = paired$Cancertype_A,
    Sample_status_A     = paired$Sample_status_A,
    Sample_comment_A    = paired$Sample_comment_A,
    Visit_name_B        = paired$Visit_name_B,
    Customer_SampleID_B = paired$Customer_SampleID_B,
    GHSampleID_B        = paired$GHSampleID_B,
    GHRequestID_B       = paired$GHRequestID_B,
    Cancertype_B        = paired$Cancertype_B,
    Sample_status_B     = paired$Sample_status_B,
    Sample_comment_B    = paired$Sample_comment_B,
    Days_between_samples = paired$days_between,
    MR_panel            = "legacy.methyl",
    MR_QC_status        = paired$MR_QC_status,
    MR_quantifiable     = paired$MR_quantifiable,
    MR_score_percentage = round(paired$.tf_b * stats::runif(nrow(paired),
                                                            0.9, 1.1), 4),
    Ctdna_percentage_change = round(paired$Ctdna_pct_change *
                                     stats::runif(nrow(paired), 0.8, 1.2), 2),
    Ctdna_level_change      = paired$Ctdna_lvl_change,
    stringsAsFactors    = FALSE
  )
  } # end if (.keep_ctdna)

  # ============================================================
  # IHC (per patient × marker, single baseline tissue sample)
  # ============================================================
  ihc <- NULL
  if (.keep_ihc) {
  ihc_rows <- list()
  iidx <- 1
  for (i in seq_len(nrow(pts))) {
    pat <- pts[i, ]
    rec_effect <- rec_drop[[pat$RECIST]]
    tissue_id  <- sprintf("%s_%s_TISSUE", study_id, pat$Patient_ID)
    collected  <- pat$.baseline_date - sample(7:30, 1)
    for (mk in c("PDL1","CD8","KI67")) {
      raw <- if (mk == "KI67") 40 + stats::rnorm(1, 0, 15)
             else 20 * rec_effect + stats::rnorm(1, 30, 15)
      score <- min(max(round(raw, 1), 0), 100)
      ihc_rows[[iidx]] <- data.frame(
        Study_ID         = study_id,
        Patient_ID       = pat$Patient_ID,
        Tissue_sample_ID = tissue_id,
        Cancertype       = pat$Cancertype,
        Visit_name       = "Baseline_tissue",
        Marker           = mk,
        Method           = ifelse(mk == "PDL1", "TPS",
                                  ifelse(mk == "CD8",
                                         "cells_per_mm2", "percent_positive")),
        Score            = score,
        Score_unit       = ifelse(mk == "CD8", "cells/mm2", "percent"),
        Collection_date  = collected,
        Reported_date    = collected + sample(7:21, 1),
        stringsAsFactors = FALSE
      )
      iidx <- iidx + 1
    }
  }
  ihc <- do.call(rbind, ihc_rows)
  rownames(ihc) <- NULL
  } # end if (.keep_ihc)

  # ============================================================
  # DNA-seq (tissue WES/WGS alteration calls)
  # ============================================================
  dnaseq <- NULL
  if (.keep_dnaseq) {
  wes_rows <- list()
  widx <- 1
  for (i in seq_len(nrow(pts))) {
    pat <- pts[i, ]
    tissue_id <- sprintf("%s_%s_TISSUE", study_id, pat$Patient_ID)
    collected <- pat$.baseline_date - sample(7:30, 1)
    n_var <- stats::rpois(1, lambda = 6) + 1
    var_types <- sample(c("SNV","Indel","CNV","SV"),
                        n_var, replace = TRUE,
                        prob = c(0.65, 0.15, 0.15, 0.05))
    for (vt in var_types) {
      g <- sample(c(.gh_genes_74, .gh_genes_500_only), 1)
      tum_vaf <- if (vt %in% c("SNV","Indel"))
        round(stats::rbeta(1, 5, 10) * 100, 2) else NA_real_
      total_d <- if (vt %in% c("SNV","Indel"))
        sample(100:500, 1) else NA_integer_
      alt_d <- if (!is.na(tum_vaf))
        round(tum_vaf * total_d / 100) else NA_integer_
      wes_rows[[widx]] <- data.frame(
        Study_ID         = study_id,
        Patient_ID       = pat$Patient_ID,
        Tissue_sample_ID = tissue_id,
        Tissue_type      = "FFPE biopsy",
        Sequencing_panel = sample(c("WES","WGS"), 1, prob = c(0.85, 0.15)),
        Cancertype       = pat$Cancertype,
        Collection_date  = collected,
        Gene             = g,
        Chromosome       = paste0("chr", sample(1:22, 1)),
        Position         = sample(1e6:2e8, 1),
        Variant_type     = vt,
        Mut_aa           = if (vt == "SNV")
          sprintf("p.%s%d%s",
                  sample(c("A","R","D","E","K","L","N","P","Q","S","T","V"), 1),
                  sample(50:600, 1),
                  sample(c("A","R","D","E","K","L","N","P","Q","S","T","V"), 1))
          else if (vt == "Indel") sprintf("p.%s%dfs", "L", sample(50:600, 1))
          else NA_character_,
        Mut_nt           = if (vt %in% c("SNV","Indel")) {
                              # v0.24.0: explicit ref>alt
                              .nb <- c("A","C","G","T")
                              .r <- sample(.nb, 1)
                              .a <- sample(setdiff(.nb, .r), 1)
                              if (vt == "SNV") sprintf("%s>%s", .r, .a)
                              else {
                                # Indel: anchored ref>alt
                                .ins <- paste(sample(.nb, sample(2:5, 1),
                                                       replace = TRUE),
                                                collapse = "")
                                if (sample(c(TRUE,FALSE), 1))
                                  sprintf("%s%s>%s", .r, .ins, .r)  # del
                                else sprintf("%s>%s%s", .r, .r, .ins)  # ins
                              }
                            } else NA_character_,
        Transcript       = sprintf("NM_%06d.%d", sample(1000:9999, 1),
                                    sample(1:5, 1)),
        Tumor_VAF        = tum_vaf,
        Normal_VAF       = if (!is.na(tum_vaf))
                             round(stats::rbeta(1, 1, 100) * 100, 3)
                           else NA_real_,
        Total_depth      = total_d,
        Alt_depth        = alt_d,
        Copy_number      = if (vt == "CNV")
                             sample(c(0, 1, 4, 6, 8, 12), 1) else NA_real_,
        Functional_impact = sample(c("deleterious",NA), 1,
                                    prob = c(0.45, 0.55)),
        Somatic_status   = sample(c("somatic","germline","somatic_putative_ch"), 1,
                                   prob = c(0.80, 0.15, 0.05)),
        COSMIC           = if (stats::runif(1) < 0.5)
                              sprintf("COSM%d", sample(1e4:9.99e5, 1))
                            else NA_character_,
        ClinVar          = sample(c(NA, "Pathogenic","Likely_pathogenic",
                                     "Uncertain_significance","Benign","Likely_benign"),
                                   1, prob = c(0.50, 0.18, 0.12, 0.10, 0.05, 0.05)),
                gnomAD_AF        = if (stats::runif(1) < 0.7) 0
                            else 10 ^ stats::runif(1, -5, -1),
        KG_AF            = if (stats::runif(1) < 0.7) 0
                            else 10 ^ stats::runif(1, -5, -1),
        dbSNP            = if (stats::runif(1) < 0.4)
                              sprintf("rs%d", sample(1e5:9e8, 1))
                            else NA_character_,
        Mutant_allele_status = sample(c(NA,"biallelic"),
                                       1, prob = c(0.85, 0.15)),
        Molecular_consequence = if (vt == "SNV") sample(
          c("missense","nonsense","synonymous","splice_donor","splice_acceptor"),
          1, prob = c(0.55, 0.15, 0.10, 0.10, 0.10))
          else if (vt == "Indel") "frameshift"
          else NA_character_,
        stringsAsFactors = FALSE
      )
      widx <- widx + 1
    }
  }
  wes <- do.call(rbind, wes_rows)
  rownames(wes) <- NULL
  dnaseq <- wes      # v0.40.0: renamed; same content, new name in output
  } # end if (.keep_dnaseq)

  # ============================================================
  # RNA-seq — count matrix + TPM + sample metadata
  # ============================================================
  rnaseq <- NULL
  if (.keep_rnaseq) {
  # Sample metadata: one tissue RNA-seq sample per patient (baseline)
  rna_meta <- data.frame(
    Sample_ID        = sprintf("%s_%s_RNA", study_id, pts$Patient_ID),
    Patient_ID       = pts$Patient_ID,
    Cancertype       = pts$Cancertype,
    Visit_name       = "Baseline_tissue",
    Tissue_type      = "FFPE biopsy",
    Library_type     = "polyA",
    Sequencing_date  = pts$.baseline_date - sample(7:30, n_patients,
                                                    replace = TRUE),
    Reads_M          = round(stats::runif(n_patients, 25, 60), 1),
    RIN_score        = round(stats::runif(n_patients, 5.5, 9.0), 1),
    stringsAsFactors = FALSE
  )

  # Gene set: real-looking symbols + numbered placeholders for the rest
  real_genes <- c(.gh_genes_74, .gh_genes_500_only,
                  "CD8A","PDL1","PDCD1","IFNG","MKI67","TGFB1","FOXP3",
                  "CD274","CTLA4","LAG3","TIGIT","HAVCR2","GZMB","PRF1",
                  "CXCL9","CXCL10","CXCL11","IDO1","CD68","CD163","FCGR3A",
                  "VIM","CDH1","CDH2","SNAI1","TWIST1","ZEB1","ESR1","PGR",
                  "AR","FOLH1","KLK3","TMPRSS2","ERG")
  n_extra <- n_rnaseq_genes - length(unique(real_genes))
  if (n_extra > 0) {
    gene_set <- c(unique(real_genes), sprintf("GENE_%05d", seq_len(n_extra)))
  } else {
    gene_set <- unique(real_genes)[seq_len(n_rnaseq_genes)]
  }

  # Per-gene baseline expression (log-scale)
  gene_mu <- stats::rnorm(length(gene_set), mean = 3, sd = 1.5)
  names(gene_mu) <- gene_set

  # Patient × gene effects driven by RECIST for immune genes
  immune_genes <- intersect(c("CD8A","PDL1","PDCD1","IFNG","GZMB","CXCL9",
                              "CXCL10","CTLA4","LAG3","TIGIT","HAVCR2"),
                            gene_set)
  prolif_genes <- intersect(c("MKI67","CDK4","CCND1"), gene_set)
  rec_score    <- rec_drop[pts$RECIST]   # higher = responder

  counts <- matrix(0L,
                    nrow = length(gene_set),
                    ncol = nrow(rna_meta),
                    dimnames = list(gene_set, rna_meta$Sample_ID))
  for (s in seq_len(ncol(counts))) {
    base <- gene_mu + stats::rnorm(length(gene_set), 0, 0.5)
    if (length(immune_genes))
      base[immune_genes] <- base[immune_genes] + 0.7 * rec_score[s]
    if (length(prolif_genes))
      base[prolif_genes] <- base[prolif_genes] - 0.5 * rec_score[s]
    mu_lin <- exp(base) * rna_meta$Reads_M[s] / 30
    counts[, s] <- stats::rnbinom(length(gene_set), mu = mu_lin, size = 5)
  }

  # TPM: scale counts by gene length proxy, then normalize per sample to 1e6
  gene_len_kb <- pmax(stats::rgamma(length(gene_set), 2, 1), 0.3)
  rpk <- counts / gene_len_kb
  per_sample_scaling <- colSums(rpk) / 1e6
  tpm <- sweep(rpk, 2, per_sample_scaling, "/")
  tpm <- round(tpm, 3)

  rnaseq <- list(
    counts          = counts,
    tpm             = tpm,
    sample_metadata = rna_meta
  )
  } # end if (.keep_rnaseq)

  # ============================================================
  # Clinical / efficacy frame (per patient, joined into ctDNA
  # downstream via ctdna_prepare()). RECIST, dose, best_pct_change
  # and MR live here because real studies deliver them in a
  # separate clinical CRF file — not on the vendor reports.
  # ============================================================
  clinical <- NULL
  if (.keep_clinical) {
  clinical <- data.frame(
    Study_ID        = study_id,
    Patient_ID      = pts$Patient_ID,
    Cancertype      = pts$Cancertype,
    Dose            = pts$Dose,
    RECIST          = pts$RECIST,
    best_pct_change = round(pts$best_pct_change, 2),
    MR              = pts$MR,
    stringsAsFactors = FALSE
  )
  } # end if (.keep_clinical)

  # ============================================================
  # ADaM ADTR (tumour assessments) — on their OWN imaging schedule,
  # deliberately DIFFERENT visit labels from the ctDNA collection visits
  # (imaging is q9w, ctDNA is per cycle). Some patients have fewer scans.
  # ============================================================
  adtr <- NULL
  {
    tvis  <- c("Screening", "Week 9", "Week 18", "Week 27", "Week 36")
    tfrac <- c(0, 0.45, 0.75, 0.92, 1.0)        # fraction toward best by scan
    arows <- list(); a <- 1
    for (i in seq_len(n_patients)) {
      bp <- pts$best_pct_change[i]; s0 <- pts$.sld_base[i]
      nv <- sample(3:5, 1)
      for (j in seq_len(nv)) {
        pchg <- if (j == 1) NA_real_ else round(bp * tfrac[j] + stats::rnorm(1, 0, 5), 1)
        sld  <- round(max(s0 * (1 + (if (is.na(pchg)) 0 else pchg) / 100) + stats::rnorm(1, 0, 1.5), 1), 1)
        arows[[a]] <- data.frame(
          SUBJID = pts$Patient_ID[i], AVISIT = tvis[j],
          PARCAT1 = "Target Lesion", PARAM = "Sum of Diameters", PARAMCD = "SUMDIAM",
          AVAL = sld, PCHG = pchg, stringsAsFactors = FALSE)
        a <- a + 1
      }
    }
    adtr <- do.call(rbind, arows)
  }

  # ============================================================
  # Assemble final output (v0.40.0: includes only requested modalities)
  # ============================================================
  out <- list()
  if (.keep_ctdna) {
    out$infinity_report      <- infinity_report
    out$tf_change            <- tf_change
    out$panel_74_response    <- panel_74_response
    out$panel_500_response   <- panel_500_response
    out$methylation_response <- methylation_response
  }
  if (.keep_clinical) out$clinical <- clinical
  out$adtr <- adtr
  if (.keep_ihc)      out$ihc      <- ihc
  if (.keep_dnaseq)   out$dnaseq   <- dnaseq
  if (.keep_rnaseq)   out$rnaseq   <- rnaseq
  attr(out, "platform")    <- "guardant_health_infinity"
  attr(out, "study_id")    <- study_id
  attr(out, "n_patients")  <- n_patients
  out
}


# ---- Pipeline ---------------------------------------------------------------

#' Run the standard ctdnaTM plotting pipeline
#'
#' Builds every standard deliverable (D1-D8, the cross-modality
#' scatters, oncoprint and alteration-grid) from a prepared dataset
#' in a single call. Each deliverable is opt-in via the available
#' modalities on the input list: plots that need data the user
#' didn't provide are silently skipped (and reported in the
#' \code{$skipped} slot).
#'
#' This is the v0.42.2 rewrite — replaces the v0.39-era
#' implementation that drifted out of sync with the v0.40-v0.42 API
#' changes. The new pipeline:
#' \itemize{
#'   \item takes one argument: the output of
#'     \code{\link{ctdna_prepare}()}. No need to disassemble the
#'     prep list and pass individual modality frames.
#'   \item dispatches via a single registry so every plot uses the
#'     current API (no stale \code{vars =}, no invalid
#'     \code{filter = list(...)}, no dead \code{stat_position}).
#'   \item collects failures and skips into named slots
#'     (\code{$failures} / \code{$skipped}) instead of swallowing
#'     them silently. A failed plot tells you which function,
#'     which arguments, and the exact error.
#'   \item exposes per-plot timing in \code{$timing}.
#'   \item lets the user pick a subset via \code{which =}, with
#'     standard short names (\code{"baseline"},
#'     \code{"reduction"}, \code{"oncoprint"}, ...).
#'   \item supports per-plot kwarg overrides via \code{overrides =}.
#' }
#'
#' @param prep A \code{ctdna_prep} object from \code{\link{ctdna_prepare}()}.
#'   Steps in the built-in registry read the standard prep frames
#'   (\code{$samples}, \code{$variants}, \code{$clinical},
#'   \code{$assessments}) plus any optional cross-modality frames the
#'   user attaches (\code{$tmb}, \code{$ihc}, \code{$expression},
#'   \code{$genomic_74}, \code{$genomic_500}). Steps whose required
#'   frame is missing are skipped with a clear reason.
#' @param which Character vector of pipeline keys to run (see
#'   \code{ctdna_pipeline_steps()} for the full list). Default
#'   \code{NULL} runs every deliverable whose modalities are
#'   present.
#' @param scheme Default RECIST stratification scheme used by all
#'   deliverables that take a \code{scheme =} argument. One of
#'   \code{"three"} (default), \code{"two"}, \code{"four"},
#'   \code{"four_alt"}, \code{"five"}, \code{"six"}.
#' @param visit Default visit passed to plots that take a
#'   \code{visit =} argument (e.g. oncoprint, alteration_grid).
#'   Default \code{"baseline"} (case-insensitive; accepts also
#'   \code{"base"}, \code{"screening"}, or a canonical visit like
#'   \code{"C1D1"}).
#' @param wrap Default column name used as the row-facet variable by
#'   plots that take a \code{wrap =} argument. \code{NULL} (default)
#'   = single panel / no row facet.
#' @param gene_sets Default gene sets passed to landscape plots
#'   (oncoprint, alteration_grid). Named list of gene-symbol vectors
#'   or combine wrappers (\code{alt_any()}, \code{alt_all()},
#'   \code{alt_any_n()}). \code{NULL} = built-in default
#'   (\code{list(HRR14 = ctdna_gene_set("HRR14"), TSG = <core TSGs>)}).
#' @param filter_scheme Default variant filter scheme(s) applied
#'   inside plotting steps. Scheme name(s), a scheme object, or a
#'   list of these. \code{NULL} = no filter.
#' @param stat Default statistical test passed to steps that offer
#'   one (e.g. contingency comparisons). \code{"fisher"} (default) or
#'   \code{"chisq"}.
#' @param engine Rendering engine hint for the landscape plots.
#'   \code{"auto"} (ComplexHeatmap if installed, else ggplot),
#'   \code{"complexheatmap"}, or \code{"ggplot"}.
#' @param overrides Named list of per-plot kwarg lists. Keys are
#'   pipeline keys; values are arg lists that override the
#'   pipeline's defaults for that plot. Example:
#'   \code{list(baseline = list(method = "maxvaf"),
#'              reduction = list(metric = "tf"))}.
#' @param compose If \code{TRUE}, also build a patchwork composite
#'   of every plot that succeeded.
#' @param verbose Print step-by-step progress messages.
#' @param dry_run If \code{TRUE}, do not build any plots — just
#'   return a list reporting which steps would run, which would be
#'   skipped (and why), per the available modalities.
#' @return A list with elements:
#'   \itemize{
#'     \item \code{plots}     — named list of ggplot objects
#'       (one per successful step).
#'     \item \code{failures}  — named list of error messages for
#'       steps that errored (empty if all succeeded).
#'     \item \code{skipped}   — named list of "reason" strings for
#'       steps skipped because a required modality was missing.
#'     \item \code{timing}    — named numeric of seconds per step.
#'     \item \code{composite} — patchwork composite (only if
#'       \code{compose = TRUE}).
#'   }
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 30, seed = 1)
#' prep <- ctdna_prepare(sim, verbose = FALSE)
#'
#' # Run everything the prep list supports
#' out <- ctdna_pipeline(prep)
#' names(out$plots)             # all generated plot keys
#' out$skipped                  # what was skipped and why
#' out$timing                   # seconds per step
#'
#' # Inspect a single plot
#' # print(out$plots$baseline)
#'
#' # Run only three plots, with custom overrides
#' out <- ctdna_pipeline(prep,
#'                       which = c("baseline","reduction","oncoprint"),
#'                       overrides = list(
#'                         baseline  = list(method = "maxvaf"),
#'                         reduction = list(metric = "tf")))
#'
#' # See what would run on this dataset without building anything
#' dry <- ctdna_pipeline(prep, dry_run = TRUE)
#' @export
ctdna_pipeline <- function(prep,
                            which         = NULL,
                            scheme        = "three",
                            visit         = "baseline",
                            wrap          = NULL,
                            gene_sets     = NULL,
                            filter_scheme = NULL,
                            stat          = "fisher",
                            engine        = "auto",
                            overrides     = list(),
                            compose       = FALSE,
                            verbose       = TRUE,
                            dry_run       = FALSE) {
  if (!is.list(prep))
    stop("`prep` must be a list (typically the output of ",
          "ctdna_prepare()).", call. = FALSE)
  registry  <- .pipeline_registry()
  all_keys  <- names(registry)
  if (is.null(which)) {
    which <- all_keys
  } else {
    bad <- setdiff(which, all_keys)
    if (length(bad))
      stop("Unknown pipeline step(s): ",
           paste(bad, collapse = ", "),
           ".\n  Available: ", paste(all_keys, collapse = ", "),
           call. = FALSE)
  }
  msg <- function(...) if (isTRUE(verbose)) message("[pipeline] ", ...)

  # ---- default gene_sets (used only if the user didn't pass one) ----
  gene_sets_is_user <- !is.null(gene_sets)
  if (!gene_sets_is_user)
    gene_sets <- list(
      HRR14 = tryCatch(ctdna_gene_set("HRR14"),
                       error = function(e) c("BRCA1","BRCA2","ATM","PALB2")),
      TSG   = c("TP53","RB1","PTEN","APC","STK11","CDKN2A","SMAD4","VHL"))

  # ---- opts bag passed to each step's arg-builder ----
  user_flags <- list(
    scheme        = !missing(scheme),
    visit         = !missing(visit),
    wrap          = !missing(wrap) && !is.null(wrap),
    gene_sets     = gene_sets_is_user,
    filter_scheme = !missing(filter_scheme) && !is.null(filter_scheme),
    stat          = !missing(stat),
    engine        = !missing(engine))
  opts <- list(scheme = scheme, visit = visit, wrap = wrap,
               gene_sets = gene_sets, filter_scheme = filter_scheme,
               stat = stat, engine = engine,
               .user = user_flags)

  out <- list(plots = list(), failures = list(), skipped = list(),
               timing = numeric())

  for (key in which) {
    spec   <- registry[[key]]
    needs  <- spec$needs %||% character(0)
    have   <- vapply(needs, function(n) !is.null(prep[[n]]) &&
                       (is.data.frame(prep[[n]]) ||
                        is.list(prep[[n]])), logical(1))
    miss   <- needs[!have]
    if (length(miss)) {
      hint <- spec$rename_hint %||% NULL
      hint_txt <- if (!is.null(hint) && any(names(hint) %in% miss)) {
        pieces <- vapply(miss, function(m) {
          if (m %in% names(hint))
            sprintf("%s -> overrides=list(%s=list(%s=prep$<yourname>))",
                    m, key, hint[[m]])
          else m
        }, character(1))
        paste0("  Rename hint: ", paste(pieces, collapse = "; "))
      } else ""
      reason <- sprintf("needs prep$%s (not found in prep).%s",
                        paste(miss, collapse = ", prep$"), hint_txt)
      out$skipped[[key]] <- reason
      msg("skip ", key, " -- ", reason)
      next
    }
    if (isTRUE(dry_run)) {
      out$plots[[key]] <- "(dry run)"
      msg("dry: ", key)
      next
    }
    # Build args: spec defaults + user-level overrides
    args <- tryCatch(
      if (length(formals(spec$args)) >= 2L && "opts" %in%
            names(formals(spec$args))) spec$args(prep, opts)
      else spec$args(prep, scheme),
      error = function(e) NULL)
    if (is.null(args)) {
      reason <- "argument builder failed"
      out$failures[[key]] <- reason
      msg("FAIL ", key, ": ", reason)
      next
    }
    over <- overrides[[key]] %||% list()
    if (length(over)) args <- utils::modifyList(args, over)
    fn <- get(spec$fn, mode = "function")
    # Drop kwargs the target function doesn't accept.
    fn_formals <- names(formals(fn))
    args <- args[intersect(names(args), fn_formals)]

    # ---- verbose defaults report (per step) ----
    if (isTRUE(verbose)) {
      tracked <- intersect(c("scheme","visit","wrap","gene_sets",
                              "filter_scheme","stat","engine"), names(args))
      if (length(tracked)) {
        tag <- function(nm) {
          u <- isTRUE(user_flags[[nm]]) || nm %in% names(over)
          v <- args[[nm]]
          v_txt <- if (is.null(v)) "NULL"
                   else if (is.list(v)) sprintf("<%d set(s)>", length(v))
                   else if (length(v) > 1L)
                     sprintf("c(%s)", paste(shQuote(as.character(v)),
                                             collapse=","))
                   else as.character(v)
          sprintf("%s=%s<%s>", nm, v_txt, if (u) "user" else "default")
        }
        msg(key, ": ", paste(vapply(tracked, tag, character(1)),
                              collapse = ", "))
      }
    }

    t0 <- Sys.time()
    res <- tryCatch(do.call(fn, args),
                     error = function(e) structure(list(err = e),
                                                    class = "pipeline_error"))
    out$timing[[key]] <- as.numeric(Sys.time() - t0, units = "secs")
    if (inherits(res, "pipeline_error")) {
      err <- conditionMessage(res$err)
      out$failures[[key]] <- err
      msg("FAIL ", key, ": ", err)
      next
    }
    # ---- normalise storage so `result$<key>$plot` always works ----
    # For a bare ggplot return, wrap as list(plot = ...). For an
    # oncoprint / grid return that already has $plot, keep as-is.
    normalise <- function(x) {
      if (inherits(x, "ctdna_oncoprint")) return(x)  # already has $plot
      if (is.list(x) && "plot" %in% names(x))    return(x)
      list(plot = x)
    }
    val <- normalise(res)
    out[[key]]        <- val
    out$plots[[key]]  <- val$plot   # back-compat: `res$plots$<key>` = the ggplot
    msg("ok   ", key, sprintf(" (%.2fs)", out$timing[[key]]))
  }

  if (isTRUE(compose) && length(out$plots)) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      msg("compose = TRUE but patchwork not installed; skipping composite")
    } else {
      msg("composing layout")
      keep <- out$plots[vapply(out$plots,
                                  function(p) inherits(p, "gg"), logical(1))]
      if (length(keep)) {
        n <- length(keep)
        ncol_ <- if (n <= 4L) 2L else 3L
        out$composite <-
          patchwork::wrap_plots(keep, ncol = ncol_) +
          patchwork::plot_annotation(
            title    = "ctdnaTM pipeline output",
            subtitle = sprintf("%d plots", n))
      }
    }
  }

  class(out) <- c("ctdna_pipeline_result", "list")
  out
}

#' List the steps the pipeline knows how to run
#'
#' Each row of the returned data frame describes one pipeline step:
#' its key (what you pass in \code{which =}), the underlying plot
#' function, and the prep-list modalities it needs.
#'
#' @return A data frame with columns \code{key}, \code{fn},
#'   \code{needs}.
#' @examples
#' ctdna_pipeline_steps()
#' @export
ctdna_pipeline_steps <- function() {
  reg <- .pipeline_registry()
  data.frame(
    key   = names(reg),
    fn    = vapply(reg, `[[`, character(1), "fn"),
    needs = vapply(reg, function(s)
                     paste(s$needs %||% character(0), collapse = ", "),
                     character(1)),
    row.names = NULL,
    stringsAsFactors = FALSE)
}

#' Print method for pipeline output
#' @param x A \code{ctdna_pipeline_result} object.
#' @param ... Ignored.
#' @return Invisibly, \code{x}.
#' @export
print.ctdna_pipeline_result <- function(x, ...) {
  cat("ctdnaTM pipeline output\n")
  cat(sprintf("  Generated: %d plots\n", length(x$plots)))
  cat(sprintf("  Failed   : %d\n",       length(x$failures)))
  cat(sprintf("  Skipped  : %d\n",       length(x$skipped)))
  cat(sprintf("  Total    : %.2fs\n",    sum(x$timing %||% 0)))
  if (length(x$failures)) {
    cat("\n  Failures:\n")
    for (k in names(x$failures))
      cat(sprintf("    - %-24s %s\n", k, x$failures[[k]]))
  }
  if (length(x$skipped)) {
    cat("\n  Skipped:\n")
    for (k in names(x$skipped))
      cat(sprintf("    - %-24s %s\n", k, x$skipped[[k]]))
  }
  cat("\n  Access:  result$<key>$plot   (or result$<key>$stats, etc.)\n",
      "         result$plots$<key>  (legacy: just the ggplot)\n",
      "  Available keys:\n   ",
       paste(setdiff(intersect(names(x),
                                names(.pipeline_registry())),
                     c("plots","failures","skipped","timing","composite")),
             collapse = ", "), "\n")
  invisible(x)
}

# Internal: pipeline registry. One entry per deliverable; keeps the
# stale-args drift problem from recurring because every plot's
# defaults live in exactly one place.
.pipeline_registry <- function() {
  list(
    # --- Sample-level per-visit plots -----------------------------------
    # All read prep$samples (long per-sample-per-visit frame with TF, MR,
    # RECIST, Dose, Indication, etc.). Plot functions take prep first-arg
    # and pull what they need internally.
    baseline = list(
      fn     = "ctdna_plot_baseline",
      needs  = "samples",
      args   = function(prep, opts) list(
        prep     = prep,
        group_by = "RECIST",
        scheme   = opts$scheme,
        title    = "D1 - Baseline by RECIST")),
    baseline_dose = list(
      fn     = "ctdna_plot_baseline",
      needs  = "samples",
      args   = function(prep, opts) list(
        prep     = prep,
        group_by = "Dose",
        scheme   = opts$scheme,
        title    = "D1b - Baseline by dose")),
    reduction = list(
      fn     = "ctdna_plot_reduction",
      needs  = "samples",
      args   = function(prep, opts) list(
        prep   = prep,
        scheme = opts$scheme,
        title  = "D3 - Best ctDNA reduction")),
    longitudinal = list(
      fn     = "ctdna_plot_longitudinal",
      needs  = "samples",
      args   = function(prep, opts) list(
        prep   = prep,
        scheme = opts$scheme,
        title  = "D5 - Longitudinal by dose")),
    pct_change_by_dose_time = list(
      fn     = "ctdna_plot_pct_change_by_dose_time",
      needs  = "samples",
      args   = function(prep, opts) list(
        df              = prep$samples,
        recist_grouping = opts$scheme,
        title           = "Percent change by dose and time")),
    # --- Cross-modality (samples + external frame) ----------------------
    vs_tmb = list(
      fn     = "ctdna_plot_vs_tmb",
      needs  = "samples",
      args   = function(prep, opts) list(
        ctdna_df = prep$samples,
        tmb_df   = prep$samples,
        title    = "ctDNA vs TMB")),
    vs_ihc = list(
      fn     = "ctdna_plot_vs_ihc",
      needs  = c("samples","ihc"),
      args   = function(prep, opts) list(
        ctdna_df = prep$samples,
        ihc_df   = prep$ihc,
        marker   = "PDL1",
        title    = "ctDNA vs PDL1 IHC")),
    vs_expression = list(
      fn     = "ctdna_plot_vs_expression",
      needs  = c("samples","expression"),
      args   = function(prep, opts) list(
        ctdna_df = prep$samples,
        expr_df  = prep$expression,
        gene     = "PDL1",
        title    = "ctDNA vs PDL1 expression")),
    vs_mutation = list(
      fn     = "ctdna_plot_vs_mutation",
      needs  = c("samples","variants"),
      args   = function(prep, opts) list(
        ctdna_df = prep$samples,
        mut_df   = prep$variants,
        gene     = "TP53",
        title    = "ctDNA vs TP53 mutation")),
    # --- Genomic landscape (variant-level) ------------------------------
    # Take prep first-arg; every user-level flag (visit/wrap/gene_sets/
    # filter_scheme/scheme/engine) flows through opts.
    oncoprint = list(
      fn     = "ctdna_oncoprint",
      needs  = "variants",
      args   = function(prep, opts) list(
        prep          = prep,
        gene_sets     = opts$gene_sets,
        visit         = opts$visit,
        wrap          = opts$wrap,
        filter_scheme = opts$filter_scheme,
        scheme        = opts$scheme,
        engine        = opts$engine,
        title         = "D9 - Genomic landscape (oncoprint)")),
    alteration_grid = list(
      fn      = "ctdna_alteration_grid",
      needs   = "variants",
      args    = function(prep, opts) list(
        prep          = prep,
        gene_sets     = opts$gene_sets,
        visit         = opts$visit,
        wrap          = opts$wrap,
        filter_scheme = opts$filter_scheme,
        scheme        = opts$scheme,
        stat          = opts$stat,
        title         = "D10 - Alteration grid")))
}


# v0.42.2: ctdna_run_pipeline() is kept as a thin wrapper around
# ctdna_pipeline() so external code that called the old function
# keeps working. The old function emitted a deprecation warning
# pointing at ctdna_pipeline.

#' Run the full ctdnaTM plotting pipeline (deprecated wrapper)
#'
#' Deprecated in v0.42.2 — call \code{\link{ctdna_pipeline}} directly.
#' This wrapper exists only so existing user code that calls
#' \code{ctdna_run_pipeline()} keeps working. It silently translates
#' the old positional layout into the new
#' \code{ctdna_pipeline(prep, ...)} contract.
#'
#' @param ctdna_df,expr_df,tmb_df,ihc_df,g74_df,g500_df,sim Legacy args.
#'   Assembled into a temporary \code{prep}-style list and passed to
#'   \code{ctdna_pipeline}.
#' @param stat_position Unused since v0.42.0 (kept for back-compat).
#' @param legend_position Default legend position, passed to plots
#'   that accept it.
#' @param scheme,plot_overrides,compose,verbose Passed through to
#'   \code{ctdna_pipeline}; \code{plot_overrides} is renamed to
#'   \code{overrides}.
#' @return Same shape as \code{ctdna_pipeline} output (a
#'   \code{ctdna_pipeline_result}).
#' @export
ctdna_run_pipeline <- function(ctdna_df,
                                expr_df  = NULL,
                                tmb_df   = NULL,
                                ihc_df   = NULL,
                                g74_df   = NULL,
                                g500_df  = NULL,
                                sim      = NULL,
                                stat_position   = NULL,
                                legend_position = .o("legend_position"),
                                scheme          = "three",
                                plot_overrides  = list(),
                                compose         = FALSE,
                                verbose         = TRUE) {
  warning("ctdna_run_pipeline() is deprecated in v0.42.2; ",
          "use ctdna_pipeline() with a prep list from ctdna_prepare().",
          call. = FALSE)
  prep <- list(ctdna = ctdna_df,
                expression = expr_df,
                tmb        = tmb_df,
                ihc        = ihc_df,
                genomic_74 = g74_df,
                genomic_500 = g500_df,
                infinity_report = if (!is.null(sim)) sim$infinity_report,
                clinical        = if (!is.null(sim)) sim$clinical)
  prep <- prep[!vapply(prep, is.null, logical(1))]
  ctdna_pipeline(prep = prep,
                  scheme = scheme,
                  overrides = plot_overrides,
                  compose = compose,
                  verbose = verbose)
}


# ---- Help / function catalog ------------------------------------------------

# The catalog itself is a constant. Each entry is a list with:
#   category   - one of the section names below
#   summary    - one-line description, shown in compact mode
#   args       - character vector of key argument names (or NULL = all)
.ctdna_function_catalog <- function() list(

  # --- Configuration ----------------------------------------------------------
  ctdna_opts = list(
    category = "Configuration",
    summary  = "Get / set / reset all package configuration (column names, plot defaults, statistical defaults).",
    args     = c("...", "load", ".reset")
  ),
  ctdna_platforms = list(
    category = "Configuration",
    summary  = "Lists supported vendor / technology combinations that ctdna_prepare() handles."
  ),

  # --- Mock data --------------------------------------------------------------
  ctdna_make_mock_study = list(
    category = "Mock data",
    summary  = "Simulate a complete Guardant Health Infinity vendor bundle (5 files + clinical + IHC + DNA-seq + RNA-seq) with realistic biology. Modality flags select which modalities to generate."
  ),

  # --- Vendor -> canonical -----------------------------------------------------
  ctdna_prepare = list(
    category = "Vendor -> canonical",
    summary  = "Transform raw vendor files (Infinity report, TF change, panel response, clinical) into canonical long-format frames."
  ),
  ctdna_assemble_modalities = list(
    category = "Vendor -> canonical",
    summary  = "Validate required columns and assemble a multi-modal database. Auto-converts expression to log2(TPM+1)."
  ),
  ctdna_to_log2tpm = list(
    category = "Vendor -> canonical",
    summary  = "Convert expression values from TPM to log2(TPM+1). Auto-detects when units attribute is missing."
  ),

  # --- Helpers ----------------------------------------------------------------
  ctdna_floor_tf = list(
    category = "Helpers",
    summary  = "Replace zeros (and values below the LOQ) with a display floor so they render on log plots."
  ),
  ctdna_stratify_recist = list(
    category = "Helpers",
    summary  = "Collapse RECIST values into clinically-useful strata. scheme = 'two' / 'three' / 'four'."
  ),
  ctdna_ratio = list(
    category = "Helpers",
    summary  = "Add per-subject baseline_tf and per-visit ratio = TF_visit / TF_baseline."
  ),
  ctdna_summary = list(
    category = "Helpers",
    summary  = "Collapse a longitudinal frame to one row per subject (baseline, best ratio, last TF, metadata)."
  ),
  ctdna_compute_mr = list(
    category = "Helpers",
    summary  = "Per-subject Molecular Response / Non-Response classifier. at = 'best' or named visit; threshold = ratio cutoff (default 0.5)."
  ),
  ctdna_metric_at = list(
    category = "Helpers",
    summary  = "Flexible per-subject ctDNA summary: metric = ratio / tf / pct_change at best / last / a named visit."
  ),
  ctdna_merge_modalities = list(
    category = "Helpers",
    summary  = "Per-subject merge of ctDNA summary with another modality (TMB, IHC, expression)."
  ),
  ctdna_filter_opts = list(
    category = "Filtering",
    summary  = "Get / set / reset the per-column alteration filter registry. Level 1 default + level 2 modified general."
  ),
  ctdna_opts_choices = list(
    category = "Configuration",
    summary  = "Show allowed values for restricted ctdna_opts keys (typo-proofing aid)."
  ),
  ctdna_opts_save = list(
    category = "Configuration",
    summary  = "Save current ctdna_opts() to an RDS file (counterpart to ctdna_opts_load)."
  ),
  ctdna_opts_load = list(
    category = "Configuration",
    summary  = "Load ctdna_opts() from an RDS file written by ctdna_opts_save."
  ),
  ctdna_coerce = list(
    category = "Helpers",
    summary  = "Coerce a data frame to canonical types (auto-strips %, fixes integer/character IDs, applies declared NA/transform rules)."
  ),
  ctdna_cleanup = list(
    category = "Helpers",
    summary  = "Sample-level QC cleanup: drop QC-failed ctDNA samples and low-library expression samples. Applied at load time by ctdna_prepare / ctdna_assemble_modalities."
  ),
  ctdna_methods = list(
    category = "Helpers",
    summary  = "List measurement methods across all three domains (tf / mr / tumor) with a used_by column."
  ),
  ctdna_register_method = list(
    category = "Helpers",
    summary  = "Register a study-specific measurement method (custom TF column, custom MR rule, etc.)."
  ),
  ctdna_help = list(
    category = "Helpers",
    summary  = "Searchable function catalog grouped by category."
  ),
  ctdna_logo = list(
    category = "Helpers",
    summary  = "Path to the bundled package logo PNG; convenient for Rmd headers via knitr::include_graphics()."
  ),
  ctdna_filter_reset_general = list(
    category = "Filtering",
    summary  = "Restore the general filter to factory defaults (level 1)."
  ),
  ctdna_annotate_population_freq = list(
    category = "Filtering",
    summary  = "Annotate alterations with gnomAD population AF. v0.25.0: routes rows to one of three gnomAD datasets (short variants -> gnomad_r2_1, CNVs -> gnomad_cnv_r4, fusions/LGR -> gnomad_sv_r2_1) by Variant_type. Cache-then-annotate: builds per-gene caches, then looks up rows from cache. Used automatically by ctdna_variant_filter() when the gnomAD_AF filter is enabled and the column is missing."
  ),
  ctdna_save_popfreq_cache = list(
    category = "Filtering",
    summary  = "v0.25.0: write the three in-session gnomAD caches (short-variant, CNV, SV) plus metadata to an RDS file. Useful for sharing the cache or persisting across R sessions."
  ),
  ctdna_load_popfreq_cache = list(
    category = "Filtering",
    summary  = "v0.25.0: restore gnomAD caches from an RDS file written by ctdna_save_popfreq_cache(). merge=TRUE (default) preserves existing in-session keys; merge=FALSE replaces them."
  ),
  ctdna_filter_scheme_get = list(
    category = "Filtering",
    summary  = "Retrieve a filtering scheme by name."
  ),
  create_filtering_scheme = list(
    category = "Filtering",
    summary  = "Build a named filtering scheme from rule-tree expressions (rule/allOf/anyOf/not). Optionally clones an existing scheme via template = '...'. Registers under `name` unless name = NULL."
  ),
  ctdna_filter_scheme_delete = list(
    category = "Filtering",
    summary  = "Delete a filtering scheme by name."
  ),
  ctdna_filter_scheme_save = list(
    category = "Filtering",
    summary  = "Save all schemes + general filter to an RDS file."
  ),
  ctdna_filter_scheme_load = list(
    category = "Filtering",
    summary  = "Load schemes + general filter from an RDS file written by ctdna_filter_scheme_save."
  ),
  ctdna_variant_filter = list(
    category = "Filtering",
    summary  = "Route a data frame through the three-level filter system. filter_scheme=NULL uses general; a name uses that specific scheme."
  ),

  # --- Statistics -------------------------------------------------------------
  ctdna_pval_label = list(
    category = "Statistics",
    summary  = "Vectorized p-value formatting: 'p=0.012', 'p<0.001', etc."
  ),
  ctdna_run_group_test = list(
    category = "Statistics",
    summary  = "Auto / Wilcoxon / Kruskal-Wallis / t-test / ANOVA depending on the number of groups."
  ),
  ctdna_pairwise_ranksum = list(
    category = "Statistics",
    summary  = "All pairwise Wilcoxon / t-tests with multiplicity reporting."
  ),
  ctdna_cor_test = list(
    category = "Statistics",
    summary  = "Spearman / Pearson / Kendall correlation with rho + p-value + n."
  ),
  ctdna_paired_test = list(
    category = "Statistics",
    summary  = "Paired Wilcoxon / paired t-test (baseline vs on-treatment)."
  ),
  ctdna_cohen_kappa = list(
    category = "Statistics",
    summary  = "Cohen's kappa for a 2x2 contingency table (chance-adjusted agreement)."
  ),
  ctdna_concordance_test = list(
    category = "Statistics",
    summary  = "Agreement % + Cohen's kappa + (2x2 only) McNemar p, in one call."
  ),

  # --- Theme & scales ---------------------------------------------------------
  ctdna_theme = list(
    category = "Theme & scales",
    summary  = "Package ggplot2 theme: clean white panels, grey borders, bold titles, italic subtitles."
  ),
  ctdna_colors = list(
    category = "Theme & scales",
    summary  = "Named character vector with RECIST colours (both collapsed labels and raw values)."
  ),
  ctdna_dose_palette = list(
    category = "Theme & scales",
    summary  = "Dose-level colour palette (distinct from the RECIST palette so they never visually collide)."
  ),
  ctdna_scale_recist = list(
    category = "Theme & scales",
    summary  = "ggplot2 colour scale for RECIST (boxes, points, lines)."
  ),
  ctdna_scale_recist_fill = list(
    category = "Theme & scales",
    summary  = "ggplot2 fill scale for RECIST (box fills, ribbon fills) using the same hues at alpha = 0.3."
  ),
  ctdna_scale_dose = list(
    category = "Theme & scales",
    summary  = "ggplot2 colour scale for dose (used when dose is the only stratifier)."
  ),
  ctdna_scale_dose_fill = list(
    category = "Theme & scales",
    summary  = "ggplot2 fill scale for dose."
  ),
  ctdna_scale_dose_shape = list(
    category = "Theme & scales",
    summary  = "ggplot2 shape scale for dose (used when RECIST is on colour and dose needs another aesthetic)."
  ),

  # --- 8 deliverable plots ----------------------------------------------------
  ctdna_plot_baseline = list(
    category = "Deliverable plots (D1-D8)",
    summary  = "D1: Boxplot of baseline value(s) grouped by any clinical factor. scheme + recist_shape supported."
  ),
  ctdna_plot_baseline_dose_recist = list(
    category = "Deliverable plots (D1-D8)",
    summary  = "D2: Baseline methylTF dodged by dose x RECIST. Pairwise within each dose, scheme-aware."
  ),
  ctdna_plot_reduction = list(
    category = "Deliverable plots (D1-D8)",
    summary  = "D3: Per-patient ctDNA reduction. metric in {ratio, tf, pct_change}; at in {best, last, named visit}; mr_threshold drawn."
  ),
  ctdna_plot_mr_recist = list(
    category = "Deliverable plots (D1-D8)",
    summary  = "D4: MR <-> RECIST 2x2 concordance heatmap with agreement / kappa / McNemar. Returns plot + table_plot + numbers."
  ),
  ctdna_plot_longitudinal = list(
    category = "Deliverable plots (D1-D8)",
    summary  = "D5: Per-subject spaghetti trajectory + thick median per dose. Optional recist_color via scheme."
  ),
  ctdna_plot_longitudinal_recist = list(
    category = "Deliverable plots (D1-D8)",
    summary  = "D6: Longitudinal trajectory faceted by dose x RECIST stratum, scheme-aware."
  ),
  ctdna_plot_reduction_vs_tumor = list(
    category = "Deliverable plots (D1-D8)",
    summary  = "D7: Scatter of ctDNA best-ratio vs % tumor change. RECIST = colour, dose = shape."
  ),
  ctdna_plot_mut_methyl = list(
    category = "Deliverable plots (D1-D8)",
    summary  = "D8: Per-baseline maxVAF vs methylTF concordance, with discordance threshold. Optional recist_shape."
  ),
  ctdna_plot_pct_change_by_dose_time = list(
    category = "Deliverable plots (D1-D8)",
    summary  = "v0.38.0: Wrapped boxplots; cols=time, rows=RECIST, boxes=dose. Y=Ctdna_percentage_change. Empty cells preserved. Optional jittered points + per-cell stats."
  ),

  # --- Cross-modality ---------------------------------------------------------
  ctdna_plot_vs_expression = list(
    category = "Cross-modality",
    summary  = "ctDNA metric vs expression of one gene. Optional recist_color via scheme."
  ),
  ctdna_plot_vs_ihc = list(
    category = "Cross-modality",
    summary  = "ctDNA metric vs IHC score of one marker. Optional recist_color via scheme."
  ),
  ctdna_plot_vs_mutation = list(
    category = "Cross-modality",
    summary  = "Boxplot of ctDNA metric by mutation status of one gene. Optional recist_facet via scheme."
  ),
  ctdna_plot_vs_tmb = list(
    category = "Cross-modality",
    summary  = "ctDNA metric vs TMB. Optional recist_color via scheme."
  ),

  # --- Genomic landscape ------------------------------------------------------
  ctdna_oncoprint = list(
    category = "Genomic landscape",
    summary  = "Oncoprint (ComplexHeatmap or ggplot engine). Patient base = variant-matrix patients at baseline filtered by filter_scheme; drops patients without canonical Indication."
  ),
  ctdna_alteration_grid = list(
    category = "Genomic landscape",
    summary  = "Stacked-bar grid: Indication x gene-set, x = response category, fill = alt vs wt, Fisher p per panel. Same patient base as ctdna_oncoprint."
  ),

  # --- Tables -----------------------------------------------------------------
  ctdna_table = list(
    category = "Tables",
    summary  = "Publication-ready ggplot table from any data frame / matrix / table. Title, subtitle, caption, alternating row stripes, optional diagonal tint."
  ),

  # --- Pipeline ---------------------------------------------------------------
  ctdna_run_pipeline = list(
    category = "Pipeline",
    summary  = "Run every standard deliverable (D1-D8 + cross-modality + oncoprint + alteration grid) in one call with consistent display options."
  )
)


#' Print a categorized listing of every function in ctdnaTM
#'
#' Prints the catalog of all 48 exported functions in this package,
#' grouped by purpose, with their key arguments and a one-line
#' description of what each one does. Useful when you can't remember
#' which function does what.
#'
#' @param category Optional. Restrict to one category. One of:
#'   `"Configuration"`, `"Mock data"`, `"Vendor -> canonical"`,
#'   `"Helpers"`, `"Statistics"`, `"Theme & scales"`,
#'   `"Deliverable plots (D1-D8)"`, `"Cross-modality"`,
#'   `"Genomic landscape"`, `"Tables"`, `"Pipeline"`.
#'   Partial matching supported (`"plot"`, `"stat"`, `"deliv"`, ...).
#' @param search Optional. A regex; only functions whose name OR
#'   summary matches are shown.
#' @param show_args If TRUE (default) print key argument names under
#'   each function.
#' @param width Maximum line width for the printed summary (default 78).
#' @return Invisibly returns a data frame (`name`, `category`, `summary`,
#'   `args`) of the matching rows, suitable for further programmatic use.
#' @examples
#' ctdna_help()                           # everything, grouped
#' ctdna_help("plot")                     # all plot functions
#' ctdna_help("stat")                     # the statistics group
#' ctdna_help(search = "RECIST")          # name or summary contains "RECIST"
#' ctdna_help(show_args = FALSE)          # compact mode: one line per fn
#'
#' # Programmatic use: get the catalog as a data frame
#' df <- ctdna_help()
#' subset(df, category == "Helpers")
#' @export
ctdna_help <- function(category = NULL,
                        search = NULL,
                        show_args = TRUE,
                        width = 78) {

  cat_tbl <- .ctdna_function_catalog()
  # Augment each entry with argument names from formals()
  rows <- lapply(names(cat_tbl), function(nm) {
    entry <- cat_tbl[[nm]]
    fn <- tryCatch(get(nm, envir = asNamespace("ctdnaTM")),
                    error = function(e) NULL)
    actual_args <- if (!is.null(fn) && is.function(fn))
      names(formals(fn)) else character(0)
    args_use <- if (is.null(entry$args)) actual_args else entry$args
    data.frame(name     = nm,
                category = entry$category,
                summary  = entry$summary,
                args     = paste(actual_args, collapse = ", "),
                stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)

  # Filter by category
  if (!is.null(category)) {
    cats <- unique(df$category)
    hit <- grep(category, cats, ignore.case = TRUE, value = TRUE)
    if (length(hit) == 0)
      stop("No category matches '", category, "'. Available: ",
           paste(cats, collapse = "; "))
    df <- df[df$category %in% hit, , drop = FALSE]
  }

  # Filter by search
  if (!is.null(search)) {
    keep <- grepl(search, df$name, ignore.case = TRUE) |
            grepl(search, df$summary, ignore.case = TRUE)
    df <- df[keep, , drop = FALSE]
    if (nrow(df) == 0) {
      message("No matches for '", search, "'.")
      return(invisible(df))
    }
  }

  # Render
  v <- utils::packageVersion("ctdnaTM")
  cat(sprintf("ctdnaTM v%s\n", v))
  cat(strrep("=", min(width, 78)), "\n", sep = "")

  for (cat_name in unique(df$category)) {
    sub <- df[df$category == cat_name, , drop = FALSE]
    cat(sprintf("\n[%s]  (%d functions)\n",
                cat_name, nrow(sub)))
    cat(strrep("-", min(width, 78)), "\n", sep = "")
    for (i in seq_len(nrow(sub))) {
      cat(sprintf("  %s()\n", sub$name[i]))
      # Wrap summary
      summary_lines <- strwrap(sub$summary[i], width = width - 6,
                                exdent = 0)
      for (line in summary_lines)
        cat(sprintf("      %s\n", line))
      if (show_args && nzchar(sub$args[i])) {
        arg_text <- paste0("args: ", sub$args[i])
        arg_lines <- strwrap(arg_text, width = width - 6, exdent = 12)
        for (line in arg_lines)
          cat(sprintf("      %s\n",
                       gsub("\u00a0", " ", line, fixed = TRUE)))
      }
    }
  }
  cat("\n", strrep("=", min(width, 78)), "\n", sep = "")
  cat(sprintf("Total: %d functions shown",
              nrow(df)))
  if (!is.null(category) || !is.null(search))
    cat(sprintf("  (filtered from %d total)",
                 length(.ctdna_function_catalog())))
  cat("\nFor full docs: ?function_name. Vignette: vignette('ctdnaTM').\n")

  invisible(df)
}


# ============================================================================
# v0.42.2: Update-check
# ============================================================================

#' Check the project URL for a newer ctdnaTM release
#'
#' Fetches \code{LATEST_VERSION.txt} from the project hosting URL
#' (\code{https://nikleotide.com/wp-content/products/ctdnatm/}) and
#' compares with the installed package version. Reports if an update
#' is available and gives the download URL.
#'
#' If \code{LATEST_VERSION.txt} cannot be reached, falls back to
#' parsing the directory listing for \code{ctdnaTM_X.Y.Z.tar.gz}
#' filenames.
#'
#' Auto-check on package load is opt-in. Set
#' \code{Sys.setenv(CTDNATM_AUTO_UPDATE_CHECK = "TRUE")} in
#' \code{~/.Rprofile} to enable. The manual call always works.
#'
#' @param quiet If \code{TRUE}, do not emit messages; the return value
#'   is still informative.
#' @return Invisibly, a list with \code{installed}, \code{latest},
#'   \code{up_to_date} (logical), and \code{download_url}. Returns
#'   \code{NULL} if the URL is unreachable or the response is
#'   malformed.
#' @examples
#' \dontrun{
#' res <- ctdna_check_for_updates()
#' if (!is.null(res) && !res$up_to_date)
#'   message("update available; see ", res$download_url)
#' }
#' @export
ctdna_check_for_updates <- function(quiet = FALSE) {
  installed <- as.character(utils::packageVersion("ctdnaTM"))
  base_url  <- "https://nikleotide.com/wp-content/products/ctdnatm"
  ver_url   <- paste0(base_url, "/LATEST_VERSION.txt")

  # Try the explicit version file first (cheap, one HTTP GET)
  latest <- tryCatch({
    con <- url(ver_url, open = "r")
    on.exit(close(con))
    readLines(con, n = 1L, warn = FALSE)
  }, error = function(e) NULL, warning = function(w) NULL)

  # Fallback: parse the directory listing for tarball filenames
  if (is.null(latest) || !length(latest) || !nzchar(latest[1L])) {
    listing <- tryCatch({
      con <- url(paste0(base_url, "/"), open = "r")
      on.exit(close(con))
      paste(readLines(con, warn = FALSE), collapse = "\n")
    }, error = function(e) NULL, warning = function(w) NULL)
    if (!is.null(listing)) {
      matches <- regmatches(
        listing,
        gregexpr("ctdnaTM_[0-9]+\\.[0-9]+\\.[0-9]+\\.tar\\.gz", listing))[[1L]]
      if (length(matches)) {
        versions <- gsub("ctdnaTM_|\\.tar\\.gz", "", matches)
        versions <- versions[order(numeric_version(versions))]
        latest   <- versions[length(versions)]
      }
    }
  }

  if (is.null(latest) || !length(latest) || !nzchar(latest[1L])) {
    if (!quiet)
      message("ctdnaTM: could not check for updates ",
              "(URL unreachable or no version metadata).")
    return(invisible(NULL))
  }
  latest <- trimws(latest[1L])
  if (!grepl("^[0-9]+\\.[0-9]+\\.[0-9]+$", latest)) {
    if (!quiet)
      message("ctdnaTM: unexpected version format from update check: ",
              latest)
    return(invisible(NULL))
  }

  up_to_date <- utils::compareVersion(latest, installed) <= 0
  download_url <- sprintf("%s/ctdnaTM_%s.tar.gz", base_url, latest)
  if (!quiet) {
    if (up_to_date) {
      message(sprintf("ctdnaTM: up to date (v%s).", installed))
    } else {
      message(sprintf(
        "ctdnaTM: a newer version is available (%s; you have %s).\n  Download: %s",
        latest, installed, download_url))
    }
  }
  invisible(list(installed    = installed,
                  latest       = latest,
                  up_to_date   = up_to_date,
                  download_url = download_url))
}
