# ============================================================================
# ctdnaTM — deliverable plots (D1-D8)
# ============================================================================
# Every plot accepts:
#   stat        choice of statistical test (see each function's docs)
#   title, subtitle, caption, xlab, ylab
#   show_stats  print test results in subtitle (default TRUE)
#   show_n      annotate group sizes (default TRUE)
#   point_size, point_alpha
#   y_scale     "linear" | "log2" | "log10" | "log" (natural) | "log_<N>"
# Defaults match the pipeline conventions doc.

# ----------------------------------------------------------------------------
# .ctdna_resolve_y_scale(): build a ggplot scale_y_* layer for the requested
# transform.
#
# y_scale options:
#   "linear" -> scale_y_continuous(...)
#   "log2"   -> scale_y_continuous(transform = "log2", ...)
#   "log10"  -> scale_y_continuous(transform = "log10", ...)
#   "log"    -> natural log (base e)
#   "log_N"  -> log base N (any positive number; e.g. "log_5")
#
# `labels` is optional and passed straight to scale_y_continuous(). If NULL,
# ggplot uses its default labeller.
# `floor_value` (optional, numeric) is the small positive value substituted
# for zeros / negatives before any log transform, to avoid -Inf. Linear
# scales ignore it.
# Returns either a ggplot scale layer, or NULL if y_scale is invalid (caller
# should warn and fall back to linear).
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# v0.43.5: auto-detect baseline visit. Used by ctdna_plot_baseline and
# ctdna_plot_baseline_dose_recist to filter df to baseline rows without
# the user having to subset manually.
#
# `baseline_visit` arg:
#   - "auto" (default): scan unique Visit_name values for common baseline
#     patterns (C1D1, C01D01, Baseline, BL, Screening, Pre-dose, etc.).
#     Pick the first match; if multiple visit labels qualify, take them all.
#   - character vector of literal visit names: use those exactly
#     (case-insensitive, trimmed).
#   - NULL: skip the filter entirely.
# ----------------------------------------------------------------------------
.ctdna_detect_baseline <- function(df, baseline_visit = "auto",
                                     visit_col = "Visit_name",
                                     verbose = FALSE) {
  if (is.null(baseline_visit)) return(df)
  if (!visit_col %in% names(df)) {
    if (verbose)
      message(sprintf("  no `%s` column found; baseline filter skipped",
                      visit_col))
    return(df)
  }

  visits     <- as.character(df[[visit_col]])
  visits_tr  <- trimws(visits)

  # Resolve baseline_visit -> a character vector of visit values to keep
  if (length(baseline_visit) == 1L &&
      tolower(as.character(baseline_visit)) == "auto") {
    fuzzy_patterns <- c(
      "^c0*1d0*1$",                    # C1D1, C01D01, c1d1
      "^cycle\\s*0*1\\s*day\\s*0*1$",  # Cycle 1 Day 1, Cycle1Day1
      "^baseline$",
      "^bl$",
      "^screen(ing)?$",
      "^pre.?treatment$",
      "^pre.?dose$",
      "^visit\\s*0*1$"
    )
    chosen <- character(0)
    for (p in fuzzy_patterns) {
      hits <- unique(visits_tr[grepl(p, visits_tr, ignore.case = TRUE)])
      if (length(hits) > 0L) {
        chosen <- hits
        break
      }
    }
    if (length(chosen) == 0L) {
      stop(sprintf(
        paste0("Auto-baseline detection found no candidate visits in `%s`. ",
               "Unique values present: %s. Pass `baseline_visit = '<your_baseline>'` ",
               "to specify (or `baseline_visit = NULL` to skip the filter)."),
        visit_col, paste(unique(visits_tr), collapse = ", ")),
        call. = FALSE)
    }
    baseline_visit <- chosen
    if (verbose)
      message(sprintf("  auto-detected baseline visit(s): %s",
                      paste(chosen, collapse = ", ")))
  }

  mask <- tolower(visits_tr) %in% tolower(trimws(as.character(baseline_visit)))
  out  <- df[mask, , drop = FALSE]
  if (nrow(out) == 0L) {
    warning(sprintf(
      paste0("ctdna baseline filter: no rows match baseline_visit = c(%s) ",
             "in column `%s`. Unique values present: %s. Returning input ",
             "unfiltered."),
      paste(shQuote(as.character(baseline_visit)), collapse=", "),
      visit_col,
      paste(unique(visits_tr), collapse=", ")), call. = FALSE)
    return(df)
  }
  if (verbose)
    message(sprintf("  baseline filter: kept %d of %d rows",
                    nrow(out), nrow(df)))
  out
}


.ctdna_resolve_y_scale <- function(y_scale = "linear", labels = NULL,
                                     floor_value = NULL, ...) {
  if (is.null(y_scale)) y_scale <- "linear"
  y_scale <- tolower(trimws(as.character(y_scale)[1]))
  args <- list(...)
  if (!is.null(labels)) args$labels <- labels
  # v0.42.5: top expansion 0.03 (was 0.00 in v0.42.4). Pure zero-top
  # cuts off bracket text glyphs that overshoot their y coordinate
  # because vjust=-0.1 puts the baseline just above the y_text value,
  # and the actual text ascender extends slightly above. 3% of span
  # gives the glyphs room to render fully without pushing the axis
  # into a different order of magnitude. Bottom expansion preserved.
  if (is.null(args$expand)) args$expand <- ggplot2::expansion(mult = c(0.05, 0.03))

  if (y_scale == "linear") {
    return(do.call(ggplot2::scale_y_continuous, args))
  }
  # v0.38.0: log1p and symlog support added
  if (y_scale == "log1p") {
    args$trans <- scales::trans_new(
      "log1p",
      transform = function(x) log1p(x),
      inverse   = function(x) expm1(x))
    return(do.call(ggplot2::scale_y_continuous, args))
  }
  if (y_scale == "symlog") {
    # Symmetric log: log of absolute value, signed. Handles negatives.
    args$trans <- scales::trans_new(
      "symlog",
      transform = function(x) sign(x) * log1p(abs(x)),
      inverse   = function(x) sign(x) * expm1(abs(x)))
    return(do.call(ggplot2::scale_y_continuous, args))
  }
  trans <- switch(y_scale,
                   "log2"  = "log2",
                   "log10" = "log10",
                   "log"   = scales::log_trans(base = exp(1)),
                   NULL)
  if (is.null(trans)) {
    # Try "log_N" pattern
    m <- regmatches(y_scale, regexec("^log[_\\.]?([0-9]+(?:\\.[0-9]+)?)$",
                                       y_scale, perl = TRUE))[[1]]
    if (length(m) == 2) {
      base <- as.numeric(m[2])
      if (!is.na(base) && base > 0 && base != 1) {
        trans <- scales::log_trans(base = base)
      }
    }
  }
  if (is.null(trans)) {
    warning(sprintf(
      "y_scale='%s' is not a recognized transform. Falling back to linear. ",
      y_scale, call. = FALSE),
      "Valid: 'linear', 'log2', 'log10', 'log', or 'log_<N>'.")
    return(do.call(ggplot2::scale_y_continuous, args))
  }
  args$trans <- trans
  do.call(ggplot2::scale_y_continuous, args)
}

# ----------------------------------------------------------------------------
# .ctdna_apply_facet(): build a ggplot facet layer from a `facet` argument.
#
# v0.38.0: universal `facet` arg semantics for deliverable plots.
#   facet = NULL                -> no faceting (default)
#   facet = "indication"        -> facet_wrap(~ Cancer_Type) (or chosen col)
#   facet = "recist"            -> facet_wrap(~ RECIST_factor)
#   facet = c("indication","recist") -> facet_grid(recist ~ indication)
#   facet = c("recist","indication") -> facet_grid(indication ~ recist)
#
# Caller must pass `df` so we can verify the column(s) exist, and a list
# `col_map` mapping facet keyword -> actual column name to use.
# Returns a ggplot facet layer, or NULL when facet is NULL/empty.
# Stops with an informative error if requested keyword has no matching col.
# ----------------------------------------------------------------------------
.ctdna_apply_facet <- function(facet, df, col_map = NULL) {
  if (is.null(facet)) return(NULL)
  if (length(facet) == 0L) return(NULL)
  if (!is.character(facet))
    stop(".ctdna_apply_facet: `facet` must be a character vector (or NULL).",
          call. = FALSE)
  if (length(facet) > 2L)
    stop(".ctdna_apply_facet: `facet` can have at most 2 elements ",
          "(got ", length(facet), "). Use one for facet_wrap, two for facet_grid.",
          call. = FALSE)

  # Default keyword -> column mapping
  default_map <- list(
    indication   = c("Cancer_Type","Cancertype","Indication","Clin_Indication"),
    cancer_type  = c("Cancer_Type","Cancertype","Indication","Clin_Indication"),
    recist       = c(".r","RECIST","Response_Subcategory","Response","BCR"),
    response     = c("Response","Response_Subcategory"),
    dose         = c("Dose","DOSE","dose"),
    arm          = c("ARM","ACTARM","TRT01A","treatment_arm"),
    cohort       = c("Cohort","COHORT","cohort"),
    sex          = c("Sex","SEX"),
    time_point   = c("time_point","AVISIT","VISIT","Time","time"))
  if (!is.null(col_map))
    default_map[names(col_map)] <- col_map

  resolve_col <- function(keyword) {
    keyword <- tolower(trimws(keyword))
    # First: maybe user passed a literal column name
    if (keyword %in% tolower(names(df))) {
      hit <- names(df)[tolower(names(df)) == keyword][1]
      return(hit)
    }
    candidates <- default_map[[keyword]]
    if (is.null(candidates))
      stop(sprintf(
        ".ctdna_apply_facet: facet keyword '%s' is not recognized. ",
        keyword),
        "Known keywords: ", paste(names(default_map), collapse = ", "),
        ". Or pass a literal column name present in the data frame.",
        call. = FALSE)
    hit <- intersect(candidates, names(df))
    if (length(hit) == 0L)
      stop(sprintf(
        ".ctdna_apply_facet: no column found in `df` for facet '%s'. ",
        keyword),
        "Looked for: ", paste(candidates, collapse = ", "), ".",
        call. = FALSE)
    hit[1]
  }

  cols <- vapply(facet, resolve_col, character(1))
  if (length(cols) == 1L) {
    return(ggplot2::facet_wrap(stats::as.formula(paste0("~ ", cols[1]))))
  }
  # length 2: facet_grid(rows ~ cols). User said:
  #   c("indication","recist") -> indication cols, recist rows
  # So second element is the column variable, first is the row variable.
  # Actually per spec: c("indication","recist") -> cols=indication, rows=recist
  # We interpret first = column facet (x), second = row facet (y).
  ggplot2::facet_grid(stats::as.formula(paste0(cols[2], " ~ ", cols[1])))
}

# ---- D1: Baseline characteristics -------------------------------------------

#' D1: Baseline characteristics — imbalance check across a grouping variable
#'
#' Faceted boxplots of baseline ctDNA-related variables across a stratifying
#' factor (dose, cohort, subpopulation). Runs a per-variable group-difference
#' test and annotates the p-value in each facet.
#'
#' @section Recommended statistic:
#' Use `"auto"` (default): Kruskal-Wallis when >2 groups, Wilcoxon rank-sum
#' for 2 groups. Use `"anova"`/`"t"` only if the variable is approximately
#' normal — uncommon for ctDNA metrics, which are typically log-distributed.
#'
#' @param df Baseline data (one row per subject).
#' @param group Stratifying column. Accepts either an opts key
#'   (e.g. `"dose"`, `"recist"`) or a literal column name on `df`
#'   (e.g. `"Dose"`, `"RECIST"`). Default: `ctdna_opts("dose")`.
#' @param specific_columns Optional character vector of column names on
#'   `df` to plot directly. When supplied, overrides `method` entirely
#'   — only these columns are used as panels. Strict: column names are
#'   case-sensitive and must exist on `df`. Use this when you need
#'   non-TF columns (e.g. `cfDNA_ng`, `n_somatic`, `TMB_score`) in
#'   addition to or instead of the TF panels. Default `NULL`, in which
#'   case `method` drives the panel set.
#' @param method TF-domain panel selector. Resolved via the method
#'   registry to actual column names of `df`. Accepts (case-insensitive
#'   as of v0.41.2):
#'   `"methyl"` (methylation TF), `"maxvaf"`, `"meanvaf"`, the meta
#'   `"all"` (expands to all three TF panels), `"custom"` (uses
#'   `value_col`), or a length-`n` vector ordering several methods.
#'   Actual column names (e.g. `"methylTF"`, `"maxVAF"`, `"meanVAF"`)
#'   are also accepted and back-mapped. Ignored if `specific_columns`
#'   is supplied. Default `"all"`.
#' @param value_col Column name to plot when `method = "custom"`.
#' @param recist_shape Optional column to encode RECIST via point
#'   shape. Accepts an opts key or a literal column name.
#' @param scheme RECIST grouping scheme passed to
#'   \code{\link{ctdna_stratify_recist}} when grouping by RECIST.
#'   One of `"four"` (default), `"three"`, `"two"`, `"four_alt"`,
#'   `"five"`, `"six"`.
#' @param scales Facet scales. One of `"free_y"` (default), `"fixed"`,
#'   `"free_x"`, `"free"`.
#' @param stat One of `"auto"` (default — Kruskal-Wallis if >2 groups,
#'   Wilcoxon if 2), `"kruskal"`, `"wilcox"`, `"anova"`, `"t"`,
#'   `"none"`. Recommended: `"auto"`.
#' @param show_stats,show_n Toggle stats subtitle and n labels.
#' @param title,subtitle,caption,xlab,ylab Plot labels.
#' @param point_size,point_alpha Display.
#' @param stat_position,legend_position See `ctdna_opts()`.
#' @param facet Optional faceting spec (length 1 or 2). Keywords:
#'   `"indication"`, `"recist"`, `"dose"`, `"arm"`, `"cohort"`, `"sex"`,
#'   `"time_point"`, or pass a literal column name. Length-2 produces
#'   facet_grid (first = columns, second = rows). NULL = no faceting.
#' @param y_scale One of `"linear"` (default), `"log10"`, `"log2"`,
#'   `"log"` (natural), `"log_<N>"` for any base, `"log1p"`, or
#'   `"symlog"` (handles signed values). Linear unless noted otherwise.
#' @param vars DEPRECATED in v0.41.3 — use `specific_columns` instead.
#'   Accepted with a one-time warning per session. Removed in a future
#'   release.
#' @param indications Optional character vector of Indication values (e.g.
#'   c("NSCLC","mCRPC")). When set, `df` is restricted to patients whose
#'   Indication is in this set. Lookup falls through Indication -> indication
#'   -> Cancertype in the frame itself; if not found in the frame, a Patient_ID
#'   join against prep$clinical is attempted. NULL (default) = no filter.
#' @return A ggplot object.
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 30, seed = 1)
#' prep <- ctdna_prepare(sim, verbose = FALSE)
#'
#' # Default: all three TF panels grouped by dose
#' ctdna_plot_baseline(prep$ctdna, group = "Dose")
#'
#' # By RECIST with the three-group scheme; bracketed pairwise stats
#' ctdna_plot_baseline(prep$ctdna, group = "RECIST",
#'                      scheme = "three", show_stats = TRUE)
#'
#' # A single TF measure on the y-axis, no inline stats
#' ctdna_plot_baseline(prep$ctdna, group = "Dose",
#'                      method = "maxVAF", show_stats = FALSE)
#'
#' # Bring in a non-TF column (cfDNA) using specific_columns =
#' ctdna_plot_baseline(prep$ctdna, group = "Dose",
#'                      specific_columns = c("methylTF", "cfDNA_ng"))
#' @export
ctdna_plot_baseline <- function(df,
                          group = .o("dose"),
                          specific_columns = NULL,
                          method = c("all","methyl","maxvaf","meanvaf","custom"),
                          value_col = NULL,
                          recist_shape = NULL,
                          scheme = c("four","three","two","four_alt","five","six"),
                          scales = c("free_y","fixed","free_x","free"),
                          stat = c("auto","kruskal","wilcox","anova","t","none"),
                          show_stats = .o("show_stats"),
                          show_n = .o("show_n"),
                          title = "Baseline characteristics",
                          subtitle = NULL,
                          caption = NULL,
                          xlab = NULL, ylab = NULL,
                          point_size = .o("point_size"),
                          point_alpha = .o("point_alpha"),
                          stat_position = .o("stat_position"),
    legend_position = .o("legend_position"),
    facet = NULL,
    y_scale = "log10",
    ylim = NULL,
    vars = NULL,
    baseline_visit = "auto",
                          indications = NULL) {
  # ---- cohort restriction ------------------------
  df <- .ctdna_filter_by_indication(df, indications)
  # -----------------------------------------------
  # v0.41.3: `vars` is the deprecated alias for `specific_columns`. Accept
  # it for back-compat but warn once per session. `specific_columns` is
  # the new explicit name — it makes clear this arg takes column names of
  # `df` directly, in contrast to `method` which takes category names
  # ("methyl", "maxvaf", "meanvaf", "all", "custom") and resolves to
  # columns via ctdna_opts(). When both are passed, specific_columns
  # wins; vars is ignored.
  if (!is.null(vars) && is.null(specific_columns)) {
    .ctdna_warn_once("ctdna_plot_baseline_vars",
      "ctdna_plot_baseline: argument `vars` is DEPRECATED; use `specific_columns`.")
    specific_columns <- vars
  } else if (!is.null(vars) && !is.null(specific_columns)) {
    .ctdna_warn_once("ctdna_plot_baseline_both",
      "ctdna_plot_baseline: both `vars` and `specific_columns` passed; using `specific_columns`.")
  }
  if (identical(method, c("all","methyl","maxvaf","meanvaf","custom")))
    method <- "all"
  stat <- match.arg(stat)
  scheme <- match.arg(scheme)
  scales <- match.arg(scales)

  # v0.43.5: auto-filter df to baseline rows. `baseline_visit = "auto"`
  # (default) scans Visit_name for common baseline labels (C1D1, Baseline,
  # Screening, etc.). Pass a literal value (e.g. "Cycle 1 Day 1") to
  # override, or NULL to skip the filter.
  df <- .ctdna_detect_baseline(df, baseline_visit = baseline_visit,
                                  visit_col = .o("visit"))

  # v0.41.1: resolve `group` and `recist_shape` to actual column names.
  # Accepts either an opts key ("dose", "recist") or a literal column name
  # ("Dose", "RECIST"). Lets pre-v0.41.0 code that used canonical lowercase
  # keys keep working after the vendor-name default flip.
  group <- .resolve_col(group, df, what = "group")
  if (!is.null(recist_shape))
    recist_shape <- .resolve_col(recist_shape, df, what = "recist_shape")

  # Resolve methods (handles "all" expansion + vector input; case-
  # insensitive per v0.41.2). `specific_columns` does NOT get the same
  # treatment — it's intentionally strict, since the user is asserting
  # exact column names.
  if (length(method) == 1L && method == "custom") {
    methods_resolved <- "custom"
  } else {
    methods_resolved <- .resolve_methods(method, domain = "tf")
  }

  # ---- v0.22 / v0.25 method dispatch ------------------------------------
  # method picks which TF column(s) appear as boxplot panels. The non-TF
  # baseline characteristics (cfDNA, n_somatic, TMB) are ALWAYS shown —
  # they are independent of TF method choice.
  #
  # method examples:
  #   "methyl"                 -> single methylTF panel
  #   "all"                    -> methyl + maxvaf + meanvaf (in that order)
  #   c("maxvaf","methyl")     -> two panels in user-specified order
  #   c("methyl","meanvaf","maxvaf") -> same as "all" but reordered
  #   "custom"                 -> user value_col only
  #
  # Per-panel y-axis labels are set via a labeller so each facet shows
  # its own units / scale description (methylTF panel says "Methylation
  # TF", TMB panel says "TMB (mut/Mb)", etc.).
  tf_panels <- setNames(
    lapply(methods_resolved, function(m)
      .resolve_method(m, "tf", value_col = if (m == "custom") value_col else NULL)),
    methods_resolved)
  tf_cols <- vapply(tf_panels, `[[`, character(1), "value_col")

  # Resolve panel labels keyed by COLUMN NAME for facet_wrap labeller
  panel_labels <- setNames(
    vapply(tf_panels, `[[`, character(1), "label"),
    tf_cols)
  panel_labels <- c(panel_labels,
    setNames(c("Max VAF"),  .o("maxvaf")),
    setNames(c("Mean VAF"), .o("meanvaf")),
    setNames(c("cfDNA concentration (ng/mL)"),  .o("cfdna_conc")),
    setNames(c("# somatic alterations"),        .o("n_somatic")),
    setNames(c("TMB (mut/Mb)"),                 .o("tmb")))
  # Dedupe (e.g. method = "maxvaf" already supplies the maxVAF label)
  panel_labels <- panel_labels[!duplicated(names(panel_labels))]

  # v0.23.0: display_unit drives both the y-axis formatter and the
  # floor caption. For the standard plotted TF panels (methyl / maxvaf
  # / meanvaf), display_unit is "percent" by default. Resolved from
  # the method registry.
  display_units <- vapply(tf_panels, function(e)
                            e$display_unit %||% "fraction",
                            character(1))
  # Dominant display unit for the plot. If all panels agree (the
  # common case), use that. If mixed (custom + something else), fall
  # back to "fraction" — the most generic.
  dominant_unit <- if (length(unique(display_units)) == 1)
                      unique(display_units) else "fraction"

  # When showing percent, augment per-panel labels with a "(%)" suffix
  # so the strip text matches the axis.
  if (dominant_unit == "percent") {
    for (col in tf_cols) {
      if (col %in% names(panel_labels) &&
          !grepl("%", panel_labels[[col]], fixed = TRUE))
        panel_labels[[col]] <- paste0(panel_labels[[col]], " (%)")
    }
  }

  # Caption: if user provided one, keep it. Otherwise compute from
  # the dominant display unit.
  if (is.null(caption)) {
    floor_val <- .o("display_floor")
    caption <- switch(dominant_unit,
      percent  = sprintf("Floor=%g%%; log y-scale", floor_val * 100),
      fraction = sprintf("Floor=%g; log y-scale", floor_val),
      actual   = "log y-scale")
  }

  if (is.null(specific_columns)) {
    # The user's method choice determines what is plotted:
    #   method = "methyl"  -> single methylTF panel
    #   method = "maxvaf"  -> single maxVAF panel
    #   method = "meanvaf" -> single meanVAF panel
    #   method = "all"     -> all three TF panels
    #   method = "custom"  -> single user-supplied value_col panel
    # No other variables sneak in. cfDNA / n_somatic / TMB are not
    # TF-domain quantities and must be requested explicitly via
    # specific_columns = .
    specific_columns <- tf_cols
    # Per-method y-axis label. For multi-method (length>1), the panels are
    # mixed (e.g. methylTF vs maxVAF), so we drop the global y-axis title
    # and let each facet's strip label identify what it shows.
    if (is.null(ylab)) {
      if (length(methods_resolved) > 1L) {
        ylab <- ""                              # suppress global title; strips disambiguate
      } else {
        base_lbl <- tf_panels[[1]]$label
        ylab <- if (dominant_unit == "percent")
                   paste0(base_lbl, " (%)") else base_lbl
      }
    }
  } else if (is.null(ylab)) {
    # User passed an explicit specific_columns =. If they passed a SINGLE
    # column, use that column's label (looked up from panel_labels if
    # known, else the bare column name). Otherwise generic.
    if (length(specific_columns) == 1 &&
        specific_columns %in% names(panel_labels)) {
      ylab <- panel_labels[[specific_columns]]
    } else {
      ylab <- "Value (log scale)"
    }
  }
  # v0.41.3: track what the user/method requested vs what's actually on df
  # so we can give an informative error when nothing matches, AND warn the
  # user when some requested columns are silently dropped (e.g. method =
  # "all" but the platform handler didn't produce one of the panels).
  requested_columns <- specific_columns
  specific_columns  <- intersect(requested_columns, names(df))
  missing_cols      <- setdiff(requested_columns, specific_columns)
  if (length(missing_cols) > 0L && length(specific_columns) > 0L) {
    .ctdna_warn_once(paste0("ctdna_plot_baseline_missing_",
                              paste(missing_cols, collapse = "_")),
      sprintf("ctdna_plot_baseline: dropping requested column(s) not on df: %s. Kept: %s.",
              paste(missing_cols, collapse = ", "),
              paste(specific_columns, collapse = ", ")))
  }
  if (length(specific_columns) == 0L)
    stop("None of the requested column(s) are in `df`.\n",
         "  Requested: ", paste(requested_columns, collapse = ", "), "\n",
         "  df columns: ", paste(names(df), collapse = ", "),
         call. = FALSE)
  # Keep the local name `vars` for the rest of the function body
  # (pivot_longer call, factor_canonical, etc.) so the diff stays small.
  vars <- specific_columns

  # If the user is grouping by RECIST, apply the chosen scheme so that
  # boxes, colors, AND the statistical groups are all in the same strata.
  recist_col <- .o("recist")
  if (identical(group, recist_col) || identical(group, "RECIST")) {
    df[[group]] <- ctdna_stratify_recist(df[[group]], scheme)
  }

  cols <- c(group, vars, if (!is.null(recist_shape)) recist_shape)
  long <- tidyr::pivot_longer(df[, cols, drop = FALSE],
                              cols = tidyr::all_of(vars),
                              names_to = "variable", values_to = "value")
  long <- long[!is.na(long$value), ]
  # v0.22.0.1: floor TF-domain values to display_floor before log-scaling.
  # Only TF-domain columns are floored — n_somatic, tmb, cfdna_conc keep
  # their natural values.
  tf_floor_vars <- unique(c(tf_cols, .o("maxvaf"), .o("meanvaf"),
                              "methylTF", "maxVAF", "meanVAF"))
  long$value <- ifelse(long$variable %in% tf_floor_vars,
                        ctdna_floor_tf(long$value),
                        long$value)
  long[[group]] <- .ctdna_factor_canonical(long[[group]], kind = group)

  long[[group]] <- .ctdna_factor_canonical(long[[group]], kind = group)

  # v0.42.0: hand-rolled comparison brackets + stats caption. show_stats
  # is now strictly TRUE/FALSE. When TRUE, brackets are drawn above
  # boxes per facet (one set per `variable`), overall test name + N go
  # in the caption.
  facet_stats <- NULL
  if (isTRUE(show_stats)) {
    facet_stats <- lapply(split(long, long$variable), function(s) {
      .compute_box_stats(s, group_col = group, value_col = "value")
    })
    facet_stats <- facet_stats[!vapply(facet_stats, is.null, logical(1))]
  }

  p <- ggplot2::ggplot(long, ggplot2::aes(.data[[group]], .data$value,
                                          color = .data[[group]],
                                          fill  = .data[[group]])) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.3,
                          linewidth = 0.7, fatten = 2.2,
                          width = .o("box_width"))

  if (!is.null(recist_shape) && recist_shape %in% names(long)) {
    long$.r <- ctdna_stratify_recist(long[[recist_shape]], scheme)
    p <- p + ggplot2::geom_jitter(data = long,
                                  ggplot2::aes(shape = .data$.r),
                                  width = 0.15, size = point_size,
                                  alpha = point_alpha,
                                  color = "grey20", inherit.aes = TRUE) +
      ggplot2::scale_shape_manual(values = c("CR/PR" = 16, "uCR/PR" = 17,
                                             "SD" = 15, "PD/NE/NA" = 18,
                                             "R" = 16, "NR" = 15,
                                             "SD/PD/NE/NA" = 15),
                                  na.value = 4)
  } else {
    p <- p + ggplot2::geom_jitter(width = 0.15, size = point_size,
                                  alpha = point_alpha,
                                  shape = 21, stroke = 0.3, color = "grey20")
  }

  recist_col <- .o("recist")
  group_is_recist <- identical(group, recist_col) || identical(group, "RECIST")
  scale_fn <- if (group_is_recist) ctdna_scale_recist else ctdna_scale_dose

  pct_labeller <- function(x) paste0(formatC(x * 100, digits = 3,
                                               format = "g"), "%")
  scale_y_fn <- if (dominant_unit == "percent") {
    function(...) .ctdna_resolve_y_scale(y_scale, labels = pct_labeller, ...)
  } else {
    function(...) .ctdna_resolve_y_scale(y_scale, ...)
  }

  p <- p + ggplot2::facet_wrap(~ variable, scales = scales,
                                labeller = ggplot2::labeller(
                                  variable = panel_labels)) +
    scale_y_fn() +
    scale_fn(aesthetics = c("color", "fill")) +
    ggplot2::labs(x = group, color = group, fill = group, shape = "RECIST") +
    ctdna_theme() +
    NULL  # legend_position handled by .finalize

  if (isTRUE(show_stats) && length(facet_stats)) {
    # Per-facet bracket layers. Each facet is a panel for a single
    # `variable`; we filter `long` to that variable and add bracketed
    # significance markers above the data in transformed y-space.
    # v0.42.4: pass facet_data so the bracket geoms get the panel's
    # `variable` value attached, keeping them in the right facet.
    # v0.42.5: new .add_brackets defaults (expand_top=0.28, base_off=0.04)
    # kick in automatically.
    for (var_name in names(facet_stats)) {
      sub <- long[long$variable == var_name, , drop = FALSE]
      p <- .add_brackets(p, facet_stats[[var_name]], sub,
                          group_col  = group, value_col = "value",
                          y_scale    = y_scale,
                          facet_data = list(variable = var_name))
    }
    # v0.42.5 caption: report the OVERALL p-value PER FACET rather
    # than broadcasting only the first facet's stats (which was
    # misleading when methods diverged). Single-facet case keeps
    # the legacy shape with N per group.
    scale_str <- if (grepl("^log", y_scale)) sprintf("%s y-scale", y_scale)
                  else "linear y-scale"
    floor_str <- switch(dominant_unit,
        percent  = sprintf("Floor=%g%%", .o("display_floor") * 100),
        fraction = sprintf("Floor=%g",   .o("display_floor")),
        actual   = NULL)
    test_names <- unique(vapply(facet_stats,
                                  function(s) s$overall$name,
                                  character(1)))
    hdr_parts <- c(floor_str, scale_str,
                    sprintf("Test: %s", paste(test_names, collapse = " / ")))
    hdr_parts <- hdr_parts[!vapply(hdr_parts, is.null, logical(1))]
    header <- paste(hdr_parts, collapse = "; ")
    if (length(facet_stats) == 1L) {
      e <- facet_stats[[1L]]
      n_str <- paste(sprintf("%s=%d",
                                names(e$n_per_group), e$n_per_group),
                       collapse = ", ")
      caption <- sprintf("%s; overall p = %s; N per group: %s",
                          header,
                          .format_p_value(e$overall$p), n_str)
    } else {
      pvals_line <- paste(
        vapply(names(facet_stats), function(nm)
          sprintf("%s: p=%s",
                    panel_labels[[nm]] %||% nm,
                    .format_p_value(facet_stats[[nm]]$overall$p)),
          character(1)),
        collapse = " | ")
      caption <- sprintf("%s\nOverall p (per panel): %s",
                          header, pvals_line)
    }
  }
  if (show_n) {
    # v0.43.8 n-labels: two complementary renderings —
    # (a) appended to x-axis tick labels as "<group>\n(n=N)" so they
    #     are ALWAYS visible, regardless of plot height, clipping, or
    #     overlap with bracket annotations. N is per-group overall
    #     (across all variables / facets), counting distinct samples.
    # (b) per-facet inside-panel `n=N` labels, rendered as filled white
    #     boxes with bold text via geom_label for high contrast. This
    #     captures the case where a particular variable has NAs in some
    #     samples (so the facet's effective n differs from the overall).
    pfc <- as.data.frame(table(long$variable, long[[group]]))
    names(pfc) <- c("variable", "g", "n")
    pfc <- pfc[pfc$n > 0, , drop = FALSE]
    pfc$g <- factor(pfc$g, levels = levels(long[[group]]))

    # (a) X-axis labels with overall N per group from the input df
    if (group %in% names(df)) {
      df_nz <- df[!is.na(df[[group]]), , drop = FALSE]
      n_overall <- setNames(as.integer(table(df_nz[[group]])),
                             names(table(df_nz[[group]])))
      p <- p + ggplot2::scale_x_discrete(labels = function(x) {
        n <- n_overall[as.character(x)]
        ifelse(is.na(n), as.character(x),
                sprintf("%s\n(n=%d)", as.character(x), n))
      })
    }

    # (b) Per-facet inside-panel labels with visible white-box styling
    if (nrow(pfc) > 0L) {
      p <- p + ggplot2::geom_label(
        data = pfc,
        mapping = ggplot2::aes(x = .data$g, y = -Inf,
                                 label = paste0("n=", .data$n)),
        inherit.aes   = FALSE,
        vjust         = -0.2,
        size          = 2.8,
        label.padding = ggplot2::unit(0.12, "lines"),
        label.size    = 0,
        color         = "grey15",
        fill          = "white",
        fontface      = "bold")
    }
  }
  # v0.42.5: coord_cartesian(clip="off") always; ylim if provided.
  cc_args <- list(clip = "off")
  if (!is.null(ylim)) cc_args$ylim <- ylim
  p <- p + do.call(ggplot2::coord_cartesian, cc_args)

  .finalize(p, title, subtitle, caption, xlab, ylab, legend_position = legend_position)
}


# ---- D2: Baseline methylTF by dose & RECIST ---------------------------------

#' D2: Baseline methylTF by dose and RECIST
#'
#' Boxplot of baseline methylTF grouped by dose and colored by RECIST.
#' Computes pairwise tests within each dose across RECIST strata.
#'
#' @section Recommended statistic:
#' `"wilcox"` (default): Wilcoxon rank-sum for two-group RECIST pairs.
#' Switch to `"t"` only if values are approximately normal on the log
#' scale.
#'
#' @param df Baseline data.
#' @param scheme RECIST stratification.
#' @param stat One of `"wilcox"` (default), `"t"`, `"none"`. Recommended: `"wilcox"`.
#' @param show_stats,show_n Toggles.
#' @param title,subtitle,caption,xlab,ylab Labels.
#' @param point_size,point_alpha Display.
#' @param facet Optional faceting spec (length 1 or 2). Keywords:
#'   `"indication"`, `"recist"`, `"dose"`, `"arm"`, `"cohort"`, `"sex"`,
#'   `"time_point"`, or pass a literal column name. Length-2 produces
#'   facet_grid (first = columns, second = rows). NULL = no faceting.
#' @param y_scale One of `"linear"` (default), `"log10"`, `"log2"`,
#'   `"log"` (natural), `"log_<N>"` for any base, `"log1p"`, or
#'   `"symlog"` (handles signed values). Linear unless noted otherwise.
#' @param indications Optional character vector of Indication values (e.g.
#'   c("NSCLC","mCRPC")). When set, `df` is restricted to patients whose
#'   Indication is in this set. Lookup falls through Indication -> indication
#'   -> Cancertype in the frame itself; if not found in the frame, a Patient_ID
#'   join against prep$clinical is attempted. NULL (default) = no filter.
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 30, seed = 1)
#' prep <- ctdna_prepare(sim, verbose = FALSE)
#'
#' # Baseline TF by dose, RECIST shown via box colour, three-group scheme
#' ctdna_plot_baseline_dose_recist(prep$ctdna, scheme = "three")
#'
#' # Use maxVAF as the y-axis measurement
#' ctdna_plot_baseline_dose_recist(prep$ctdna,
#'                                  method = "maxvaf", scheme = "three")
#' @export
ctdna_plot_baseline_dose_recist <- function(df,
                                      scheme = c("three","two","four","four_alt","five","six"),
                                      method = "methyl",
                                      value_col = NULL,
                                      stat = c("wilcox","t","none"),
                                      show_stats = .o("show_stats"),
                                      show_n = .o("show_n"),
                                      title = NULL,
                                      subtitle = NULL,
                                      caption = sprintf("Floor=%g; log y", .o("display_floor")),
                                      xlab = NULL,
                                      ylab = NULL,
                                      point_size = .o("point_size"),
                                      point_alpha = .o("point_alpha"),
    stat_position = .o("stat_position"),
    legend_position = .o("legend_position"),
    facet = NULL,
    y_scale = "log10",
    ylim = NULL,
    baseline_visit = "auto",
                          indications = NULL) {
  # ---- cohort restriction ------------------------
  df <- .ctdna_filter_by_indication(df, indications)
  # -----------------------------------------------
  # v0.43.5: auto-filter df to baseline rows.
  df <- .ctdna_detect_baseline(df, baseline_visit = baseline_visit,
                                  visit_col = .o("visit"))

  # v0.25: method accepts a vector for side-by-side comparison. When
  # length > 1, we recursively call this function once per method and
  # stack the resulting plots vertically via patchwork.
  if (length(method) > 1L || identical(method, "all")) {
    methods_resolved <- .resolve_methods(method, domain = "tf")
    plots <- lapply(methods_resolved, function(m) {
      ctdna_plot_baseline_dose_recist(df, scheme = scheme, method = m,
        value_col = if (m == "custom") value_col else NULL,
        stat = stat, show_stats = show_stats, show_n = show_n,
        title = if (!is.null(title)) sprintf("%s — %s", title,
                  .resolve_method(m, "tf")$label) else NULL,
        subtitle = subtitle, caption = caption, xlab = xlab, ylab = ylab,
        point_size = point_size, point_alpha = point_alpha,
        stat_position = stat_position, legend_position = legend_position,
        y_scale = y_scale, ylim = ylim,
        baseline_visit = NULL)   # already filtered above; don't refilter
    })
    return(patchwork::wrap_plots(plots, ncol = 1L))
  }

  .stats_label <- NULL

  scheme <- match.arg(scheme); stat <- match.arg(stat)
  m <- .resolve_method(method, "tf", value_col = value_col)
  tf_col <- m$value_col
  if (is.null(title)) title <- paste0("Baseline ", m$label, " by dose and RECIST")
  if (is.null(ylab))  ylab  <- paste0(m$label, " (log)")
  df$.r <- ctdna_stratify_recist(df[[.o("recist")]], scheme)
  df$.tf <- ctdna_floor_tf(df[[tf_col]])
  df[[.o("dose")]] <- .ctdna_factor_canonical(df[[.o("dose")]], kind = "dose")

  # v0.42.0: between-dose stats with bracket rendering above the boxes.
  # show_stats is now strictly TRUE/FALSE. The RECIST colour split is
  # retained as visual context; brackets compare the dose distributions
  # marginal of RECIST. Per-RECIST-within-dose detail is summarised in
  # the caption when more than one RECIST group is present.
  stats_res <- NULL
  if (isTRUE(show_stats))
    stats_res <- .compute_box_stats(df, group_col = .o("dose"),
                                      value_col = ".tf")

  p <- ggplot2::ggplot(df, ggplot2::aes(.data[[.o("dose")]], .data$.tf,
                                        color = .data$.r, fill = .data$.r)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.3,
                          linewidth = 0.7, fatten = 2.2,
                          width = .o("box_width"),
                          position = ggplot2::position_dodge(0.8)) +
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(
      jitter.width = 0.15, dodge.width = 0.8),
      size = point_size, alpha = point_alpha,
      shape = 21, stroke = 0.3, color = "grey20") +
    ctdna_scale_recist(aesthetics = c("color", "fill")) + .log_clearance(y_scale) +
    ggplot2::labs(x = .o("dose"), color = "RECIST", fill = "RECIST") +
    ctdna_theme()

  if (isTRUE(show_stats) && !is.null(stats_res)) {
    # v0.42.5: .add_brackets gets new defaults (expand_top=0.28,
    # base_off=0.04) and linear-anchor-at-v_data_max fix automatically.
    p <- .add_brackets(p, stats_res, df,
                        group_col = .o("dose"), value_col = ".tf",
                        y_scale = y_scale)
    caption <- .compose_stats_caption(stats_res, "box", base_caption = caption)
  }
  if (show_n) {
    # v0.43.8: same two-rendering approach as D1 — x-axis label gets
    # "(n=N)" appended for always-on visibility; per-group label
    # rendered inside the panel as a white box with bold text.
    counts <- as.data.frame(table(df[[.o("dose")]]))
    names(counts) <- c("g", "n")
    counts <- counts[counts$n > 0, , drop = FALSE]
    counts$g <- factor(counts$g, levels = levels(df[[.o("dose")]]))

    # (a) X-axis labels
    n_overall <- setNames(as.integer(counts$n), as.character(counts$g))
    p <- p + ggplot2::scale_x_discrete(labels = function(x) {
      n <- n_overall[as.character(x)]
      ifelse(is.na(n), as.character(x),
              sprintf("%s\n(n=%d)", as.character(x), n))
    })

    # (b) Inside-panel labels
    if (nrow(counts) > 0L) {
      p <- p + ggplot2::geom_label(
        data = counts,
        mapping = ggplot2::aes(x = .data$g, y = -Inf,
                                 label = paste0("n=", .data$n)),
        inherit.aes   = FALSE,
        vjust         = -0.2,
        size          = 2.8,
        label.padding = ggplot2::unit(0.12, "lines"),
        label.size    = 0,
        color         = "grey15",
        fill          = "white",
        fontface      = "bold")
    }
  }
  # v0.42.5: coord_cartesian(clip="off") always; ylim if provided.
  cc_args <- list(clip = "off")
  if (!is.null(ylim)) cc_args$ylim <- ylim
  p <- p + do.call(ggplot2::coord_cartesian, cc_args)

  .finalize(p, title, subtitle, caption, xlab, ylab,
             legend_position = legend_position)
}


# ---- D3: Reduction distribution ---------------------------------------------

#' D3: ctDNA reduction distribution across dose and RECIST
#'
#' Boxplot of on-treatment / baseline methylTF ratio, faceted by dose and
#' colored by RECIST. LOQ exclusion applied.
#'
#' @section Recommended statistic:
#' `"wilcox"` (default): Wilcoxon rank-sum within each dose. Reduction
#' ratios are heavily skewed even on log scale, so non-parametric is
#' the right call.
#'
#' @inheritParams ctdna_plot_baseline_dose_recist
#' @param indications Optional character vector of Indication values (e.g.
#'   c("NSCLC","mCRPC")). When set, `df` is restricted to patients whose
#'   Indication is in this set. Lookup falls through Indication -> indication
#'   -> Cancertype in the frame itself; if not found in the frame, a Patient_ID
#'   join against prep$clinical is attempted. NULL (default) = no filter.
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 30, seed = 1)
#' prep <- ctdna_prepare(sim, verbose = FALSE)
#'
#' # Default — best-ratio (deepest drop from baseline), grouped by dose
#' ctdna_plot_reduction(prep$ctdna, scheme = "three", metric = "ratio",
#'                       at = "best", mr_threshold = 0.5)
#'
#' # Percent-change at a specific visit (tolerant matching: "C2D1" works
#' # but so does "Cycle 2 Day 1")
#' ctdna_plot_reduction(prep$ctdna, scheme = "three",
#'                       metric = "pct_change", at = "Cycle 2 Day 1")
#'
#' # Raw TF at the best on-treatment visit
#' ctdna_plot_reduction(prep$ctdna, scheme = "three",
#'                       metric = "tf", at = "best")
#' @export
ctdna_plot_reduction <- function(df,
                           scheme = c("three","two","four","four_alt","five","six"),
                           at     = "best",
                           metric = c("ratio","tf","pct_change"),
                           method = "methyl",
                           value_col = NULL,
                           mr_threshold = 0.5,
                           stat = c("wilcox","t","none"),
                           show_stats = .o("show_stats"),
                           show_n = .o("show_n"),
                           title = "ctDNA reduction by dose and RECIST",
                           subtitle = NULL,
                           caption = NULL,
                           xlab = NULL,
                           ylab = NULL,
                           point_size = .o("point_size"),
                           point_alpha = .o("point_alpha"),
    stat_position = .o("stat_position"),
    legend_position = .o("legend_position"),
    facet = NULL,
    y_scale = "log10",
    ylim = NULL,
    y_clip = NULL,
                          indications = NULL) {
  # ---- cohort restriction ------------------------
  df <- .ctdna_filter_by_indication(df, indications)
  # -----------------------------------------------
  # v0.42.5: ylim is the canonical name (matches D1/D2). y_clip kept
  # as a deprecated alias for back-compat — if user passes y_clip and
  # not ylim, promote y_clip into ylim.
  if (is.null(ylim) && !is.null(y_clip)) ylim <- y_clip

  # v0.42.5 BUG-FIX: match.arg(metric) MUST run before the multi-method
  # recursion. Inside the lapply closure, match.arg() looks up the
  # formals of the calling function, which is the anonymous function
  # (not ctdna_plot_reduction), and errors with "arg must be of
  # length 1" when metric is the default vector.
  metric <- match.arg(metric)

  # v0.25: method accepts a vector for side-by-side comparison.
  if (length(method) > 1L || identical(method, "all")) {
    methods_resolved <- .resolve_methods(method, domain = "tf")
    plots <- lapply(methods_resolved, function(m) {
      ctdna_plot_reduction(df, scheme = scheme, at = at, metric = metric,
        method = m,
        value_col = if (m == "custom") value_col else NULL,
        mr_threshold = mr_threshold, stat = stat,
        show_stats = show_stats, show_n = show_n,
        title = sprintf("%s — %s", title, .resolve_method(m, "tf")$label),
        subtitle = subtitle, caption = caption, xlab = xlab, ylab = ylab,
        point_size = point_size, point_alpha = point_alpha,
        stat_position = stat_position, legend_position = legend_position,
        y_scale = y_scale, ylim = ylim)
    })
    return(patchwork::wrap_plots(plots, ncol = 1L))
  }

  .stats_label <- NULL

  scheme <- match.arg(scheme); stat <- match.arg(stat)
  m <- .resolve_method(method, "tf", value_col = value_col)
  tf_col <- m$value_col

  # v0.42.5 BUG-FIX: ctdna_metric_at computes pct_change as
  #   100 * (tf_at - baseline_tf) / baseline_tf.
  # When baseline_tf is near LOQ (e.g. 0.001), even a moderate tf_at
  # blows up the metric (0.5 / 0.001 = +49900%). Floor baseline at LOQ
  # before dividing, consistent with how the rest of the library
  # treats sub-LOQ values (display_floor). Only relevant for pct_change;
  # ratio / tf metrics use log scale + ctdna_floor_tf on the output.
  if (metric == "pct_change") {
    s_col   <- .o("subject")
    t_col   <- .o("time")
    bl      <- .o("baseline")
    floor_v <- .o("display_floor")
    base <- df[df[[t_col]] == bl, c(s_col, tf_col), drop = FALSE]
    names(base) <- c(s_col, "baseline_tf_raw")
    base$baseline_tf <- pmax(base$baseline_tf_raw, floor_v, na.rm = TRUE)
    ma_raw <- ctdna_metric_at(df, metric = "tf", at = at, tf_col = tf_col)
    ma_raw <- merge(ma_raw[, c(s_col, "value", "time_used")],
                     base[, c(s_col, "baseline_tf")],
                     by = s_col, all.x = TRUE)
    names(ma_raw)[names(ma_raw) == "value"] <- "tf_at"
    ma <- data.frame(
      stats::setNames(list(ma_raw[[s_col]]), s_col),
      baseline_tf = ma_raw$baseline_tf,
      value       = 100 * (ma_raw$tf_at - ma_raw$baseline_tf) /
                         ma_raw$baseline_tf,
      time_used   = ma_raw$time_used,
      check.names = FALSE,
      stringsAsFactors = FALSE)
  } else {
    ma <- ctdna_metric_at(df, metric = metric, at = at, tf_col = tf_col)
  }

  # Join the per-subject RECIST and dose
  meta_cols <- intersect(c(.o("recist"), .o("dose")), names(df))
  meta <- unique(df[, c(.o("subject"), meta_cols), drop = FALSE])
  d <- merge(ma, meta, by = .o("subject"), all.x = TRUE)

  d$.r <- ctdna_stratify_recist(d[[.o("recist")]], scheme)
  d[[.o("dose")]] <- .ctdna_factor_canonical(d[[.o("dose")]], kind = "dose")

  # Decide y-axis & log scaling based on metric
  use_log <- metric %in% c("ratio","tf")
  if (use_log) {
    d$.y <- ctdna_floor_tf(d$value)
  } else {
    d$.y <- d$value          # pct_change can be negative; linear scale
  }
  d_plotted <- d[!is.na(d$.y) & is.finite(d$.y), , drop = FALSE]

  # v0.42.5 BUG-FIX: resolve ylim EARLY so we can cap the data BEFORE
  # stats / brackets compute. Default for pct_change = c(-100, 200).
  # Then build d_capped (same data as d_plotted but .y clipped to
  # ylim) and use d_capped for boxplot, stats AND brackets. With
  # clip="off", uncapped outliers would render WAY above the panel,
  # making the y-axis look like it goes to 10000+%.
  if (is.null(ylim) && metric == "pct_change") ylim <- c(-100, 200)
  d_capped <- d_plotted
  if (!is.null(ylim)) {
    d_capped$.y <- pmin(pmax(d_plotted$.y, ylim[1]), ylim[2])
  }

  # MR classification (only meaningful for the ratio metric)
  if (metric == "ratio") {
    d_capped$MR <- ifelse(is.na(d_capped$value), NA_character_,
                    ifelse(d_capped$value <= mr_threshold,
                           "Response", "NonResponse"))
  }

  # v0.42.5 BUG-FIX: ylab uses m$label (e.g. "Max VAF") instead of
  # hardcoded "TF", so multi-method patchwork stacks display each
  # method's actual name on its y-axis.
  if (is.null(ylab)) ylab <- switch(metric,
    ratio      = sprintf("%s(%s) / %s(baseline) [log]",
                          m$label,
                          if (identical(at, "best")) "best" else at,
                          m$label),
    tf         = sprintf("%s at %s [log]",
                          m$label,
                          if (identical(at, "best")) "best" else at),
    pct_change = sprintf("Percent change in %s at %s",
                          m$label,
                          if (identical(at, "best")) "best" else at))
  if (is.null(caption)) caption <- switch(metric,
    ratio      = sprintf("metric = ratio at %s; floor=%g, log y; MR threshold = %g (ratio<=thr -> Response)",
                          if (identical(at, "best")) "best" else at,
                          .o("display_floor"), mr_threshold),
    tf         = sprintf("metric = TF at %s; floor=%g, log y",
                          if (identical(at, "best")) "best" else at,
                          .o("display_floor")),
    pct_change = sprintf("metric = pct_change at %s; negative = reduction",
                          if (identical(at, "best")) "best" else at))

  # v0.42.0: between-dose stats with bracket rendering above the boxes.
  stats_res <- NULL
  if (isTRUE(show_stats))
    stats_res <- .compute_box_stats(d_capped, group_col = .o("dose"),
                                      value_col = ".y")

  # v0.42.5: build the plot from d_capped (NOT d_plotted) so outlier
  # points don't render via clip="off" above the panel.
  p <- ggplot2::ggplot(d_capped, ggplot2::aes(.data[[.o("dose")]], .data$.y,
                                       color = .data$.r, fill = .data$.r)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.3,
                          linewidth = 0.7, fatten = 2.2,
                          width = .o("box_width"),
                          position = ggplot2::position_dodge(0.8)) +
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(
      jitter.width = 0.15, dodge.width = 0.8),
      size = point_size, alpha = point_alpha,
      shape = 21, stroke = 0.3, color = "grey20") +
    ctdna_scale_recist(aesthetics = c("color", "fill")) +
    ggplot2::labs(x = .o("dose"), color = "RECIST", fill = "RECIST") +
    ctdna_theme()

  # Log scale + clearance line only for ratio / tf
  if (use_log) p <- p + .log_clearance(y_scale)

  # Dashed line at the MR threshold so MR responders are visually
  # obvious. Label sits at the LEFT edge so it never overlaps boxplot
  # data on the right.
  if (metric == "ratio" && !is.null(mr_threshold) &&
      is.finite(mr_threshold) && mr_threshold > 0) {
    p <- p +
      ggplot2::geom_hline(yintercept = mr_threshold, linetype = "dashed",
                            color = "grey25", linewidth = 0.4) +
      ggplot2::annotate("text", x = -Inf, y = mr_threshold,
                         label = sprintf("MR threshold = %g", mr_threshold),
                         hjust = -0.05, vjust = -0.5,
                         size = 3, color = "grey25")
  }

  if (show_n) {
    # v0.43.8: x-axis labels + visible inside-panel n (same pattern as D1/D2)
    counts <- as.data.frame(table(d_capped[[.o("dose")]]))
    names(counts) <- c("g", "n")
    counts <- counts[counts$n > 0, , drop = FALSE]
    counts$g <- factor(counts$g, levels = levels(d_capped[[.o("dose")]]))

    n_overall <- setNames(as.integer(counts$n), as.character(counts$g))
    p <- p + ggplot2::scale_x_discrete(labels = function(x) {
      n <- n_overall[as.character(x)]
      ifelse(is.na(n), as.character(x),
              sprintf("%s\n(n=%d)", as.character(x), n))
    })
    if (nrow(counts) > 0L) {
      p <- p + ggplot2::geom_label(
        data = counts,
        mapping = ggplot2::aes(x = .data$g, y = -Inf,
                                 label = paste0("n=", .data$n)),
        inherit.aes   = FALSE,
        vjust         = -0.2,
        size          = 2.8,
        label.padding = ggplot2::unit(0.12, "lines"),
        label.size    = 0,
        color         = "grey15",
        fill          = "white",
        fontface      = "bold")
    }
  }
  if (isTRUE(show_stats) && !is.null(stats_res)) {
    # v0.42.5: brackets compute on d_capped (NOT d_plotted) so bracket
    # math stays within the visible range; new .add_brackets defaults
    # (expand_top=0.28, base_off=0.04) + linear-anchor fix kick in.
    p <- .add_brackets(p, stats_res, d_capped,
                        group_col = .o("dose"), value_col = ".y",
                        y_scale = if (use_log) "log10" else "linear")
    caption <- .compose_stats_caption(stats_res, "box", base_caption = caption)
  }
  # v0.42.5: coord_cartesian(clip="off") always; ylim if provided.
  cc_args <- list(clip = "off")
  if (!is.null(ylim)) cc_args$ylim <- ylim
  p <- p + do.call(ggplot2::coord_cartesian, cc_args)

  .finalize(p, title, subtitle, caption, xlab, ylab,
             legend_position = legend_position)
}


# ---- D4: MR-RECIST concordance ----------------------------------------------

#' D4: MR-RECIST concordance with kappa and agreement
#'
#' Heatmap of the MR x RECIST contingency table with diagonal-emphasis
#' coloring. Reports overall % agreement, Cohen's kappa, and (for 2x2)
#' McNemar's test.
#'
#' @section Recommended statistic:
#' `"both"` (default): shows agreement % and Cohen's kappa together.
#' Use `"mcnemar"` if you specifically want a marginal-homogeneity test
#' (only valid for 2x2 tables).
#'
#' @section Source of MR values:
#' By default MR is **recomputed from longitudinal ctDNA dynamics** via
#' [ctdna_compute_mr()] — taking the deepest on-treatment / baseline TF ratio
#' per subject and calling Response if that best ratio is at or below
#' `mr_threshold` (default 0.5).
#'
#' Set `mr_timepoint = "C5D1"` (or any other on-treatment label) to use
#' the MR call at that specific visit instead of the best.
#'
#' Set `mr_timepoint = "use_column"` to fall back to a pre-computed MR
#' column already present in `df` (the legacy behaviour) — useful when
#' the MR call was made upstream (e.g. by the vendor or by a separate
#' clinical script).
#'
#' @param df Longitudinal ctDNA data (subject × time × TF × RECIST).
#' @param scheme RECIST stratification.
#' @param mr_timepoint One of `"best"` (default — deepest ratio across
#'   visits), a specific on-treatment time-point label (e.g. `"C5D1"`),
#'   or `"use_column"` (use the MR column already on `df`).
#' @param mr_threshold Ratio (TF / baseline TF) at or below which a
#'   subject is called a molecular responder. Default 0.5.
#' @param stat One of `"both"` (default), `"kappa"`, `"agreement"`,
#'   `"mcnemar"`, `"none"`. Recommended: `"both"`.
#' @param show_stats,title,subtitle,caption,xlab,ylab Display options.
#' @param indications Optional character vector of Indication values (e.g.
#'   c("NSCLC","mCRPC")). When set, `df` is restricted to patients whose
#'   Indication is in this set. Lookup falls through Indication -> indication
#'   -> Cancertype in the frame itself; if not found in the frame, a Patient_ID
#'   join against prep$clinical is attempted. NULL (default) = no filter.
#' @return list(plot, table, agreement, kappa, mcnemar_p, mr_source).
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 30, seed = 1)
#' prep <- ctdna_prepare(sim, verbose = FALSE)
#'
#' # Concordance heatmap with overall statistics in caption
#' ctdna_plot_mr_recist(prep$ctdna, scheme = "three", stat = "both")
#' @export
ctdna_plot_mr_recist <- function(df,
                           scheme = c("three","two","four","four_alt","five","six"),
                           mr_timepoint = "best",
                           mr_threshold = 0.5,
                           stat = c("both","kappa","agreement","mcnemar","none"),
                           show_stats = .o("show_stats"),
                           title = "MR-RECIST concordance",
                           subtitle = NULL,
                           caption = NULL,
                           xlab = "RECIST",
                           ylab = "Molecular Response",
    stat_position = .o("stat_position"),
    legend_position = .o("legend_position"),
    facet = NULL,
    y_scale = "linear",
                          indications = NULL) {  # ---- cohort restriction ------------------------
  df <- .ctdna_filter_by_indication(df, indications)
  # -----------------------------------------------
  # categorical y; arg present for API consistency
  .stats_label <- NULL

  scheme <- match.arg(scheme); stat <- match.arg(stat)

  # Resolve MR values: compute fresh, or use the column already on df.
  if (identical(mr_timepoint, "use_column")) {
    if (!.o("mr") %in% names(df))
      stop("mr_timepoint='use_column' requested but column '", .o("mr"),
           "' not present in df.")
    per_subj <- unique(df[, c(.o("subject"), .o("mr"), .o("recist"))])
    mr_vals  <- per_subj[[.o("mr")]]
    rec_vals <- per_subj[[.o("recist")]]
    mr_source <- sprintf("MR taken from df$%s", .o("mr"))
  } else {
    mr_df <- ctdna_compute_mr(df, at = mr_timepoint, threshold = mr_threshold)
    rec_lookup <- unique(df[, c(.o("subject"), .o("recist"))])
    merged <- merge(mr_df, rec_lookup, by = .o("subject"))
    mr_vals  <- merged$MR
    rec_vals <- merged[[.o("recist")]]
    mr_source <- sprintf(
      "MR recomputed at %s; ratio threshold = %g",
      if (identical(mr_timepoint, "best")) "best on-treatment visit"
      else paste0("visit ", mr_timepoint),
      mr_threshold)
  }

  if (is.null(caption)) caption <- mr_source

  r <- ctdna_stratify_recist(rec_vals, scheme)
  m <- factor(mr_vals)
  tab <- table(MR = m, RECIST = r)
  ct <- ctdna_concordance_test(tab, method = stat)

  # v0.42.0: concordance stats are single statistics (kappa, agreement,
  # McNemar p, N) so they fit naturally in the caption. show_stats is
  # strictly TRUE/FALSE.
  if (isTRUE(show_stats) && stat != "none") {
    parts <- c()
    if (stat %in% c("both","agreement"))
      parts <- c(parts, sprintf("Agreement = %.0f%%",
                                100 * ct$agreement))
    if (stat %in% c("both","kappa"))
      parts <- c(parts, sprintf("kappa = %.2f", ct$kappa))
    if (stat == "mcnemar" && !is.na(ct$mcnemar_p))
      parts <- c(parts, sprintf("McNemar %s", ctdna_pval_label(ct$mcnemar_p)))
    parts <- c(parts, sprintf("N = %d", ct$n))
    extra <- paste(parts, collapse = "  |  ")
    caption <- if (is.null(caption) || !nzchar(caption)) extra
                else paste(caption, extra, sep = "  |  ")
  }

  long <- as.data.frame(tab)
  p <- ggplot2::ggplot(long, ggplot2::aes(.data$RECIST, .data$MR,
                                          fill = .data$Freq)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = .data$Freq)) +
    ggplot2::scale_fill_gradient(low = "white", high = "#3182BD") +
    ctdna_theme()
  p <- .finalize(p, title, subtitle, caption, xlab, ylab,
                  legend_position = legend_position)

  # Publication-style table panel for the same contingency
  tbl_caption <- paste(
    "Agreement = % of patients whose MR class matches their RECIST class on the diagonal.",
    if (!is.null(mr_source)) mr_source else "",
    sep = "  ")
  table_plot <- ctdna_table(
    as.data.frame.matrix(tab),
    title    = if (is.null(title)) NULL else paste(title, "(table)"),
    subtitle = sprintf("Agreement %.0f%% | Cohen's kappa = %.2f | N = %d",
                       100 * ct$agreement, ct$kappa, ct$n),
    caption  = tbl_caption,
    highlight_diag = TRUE)

  list(plot = p, table = tab, table_plot = table_plot,
       agreement = ct$agreement,
       kappa = ct$kappa, mcnemar_p = ct$mcnemar_p,
       mr_source = mr_source)
}


# ---- D5: Longitudinal by dose -----------------------------------------------

#' D5: Longitudinal ctDNA reduction by dose (RECIST aggregated)
#'
#' Per-subject spaghetti + thick median per dose. Computes a paired test
#' between baseline and each on-treatment time point, per dose.
#'
#' @section Recommended statistic:
#' `"paired_wilcox"` (default): paired Wilcoxon — appropriate because
#' the same subject is observed at baseline and each on-treatment time
#' point. Switch to `"paired_t"` only if log-ratios appear normal.
#'
#' @param df Longitudinal data.
#' @param stat One of `"paired_wilcox"` (default), `"paired_t"`,
#'   `"wilcox"`, `"t"`, `"none"`. Recommended: `"paired_wilcox"`.
#' @param show_stats,show_n,title,subtitle,caption,xlab,ylab Display options.
#' @param point_size,point_alpha Display.
#' @param facet Optional faceting spec (length 1 or 2). Keywords:
#'   `"indication"`, `"recist"`, `"dose"`, `"arm"`, `"cohort"`, `"sex"`,
#'   `"time_point"`, or pass a literal column name. Length-2 produces
#'   facet_grid (first = columns, second = rows). NULL = no faceting.
#' @param y_scale One of `"linear"` (default), `"log10"`, `"log2"`,
#'   `"log"` (natural), `"log_<N>"` for any base, `"log1p"`, or
#'   `"symlog"` (handles signed values). Linear unless noted otherwise.
#' @param indications Optional character vector of Indication values (e.g.
#'   c("NSCLC","mCRPC")). When set, `df` is restricted to patients whose
#'   Indication is in this set. Lookup falls through Indication -> indication
#'   -> Cancertype in the frame itself; if not found in the frame, a Patient_ID
#'   join against prep$clinical is attempted. NULL (default) = no filter.
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 15, seed = 1)
#' prep <- ctdna_prepare(sim, verbose = FALSE)
#'
#' # Per-subject TF trajectories
#' ctdna_plot_longitudinal(prep$ctdna)
#'
#' # Colored by dose
#' ctdna_plot_longitudinal(prep$ctdna, color_by = "Dose")
#' @export
ctdna_plot_longitudinal <- function(df,
                              stat = c("paired_wilcox","paired_t","wilcox","t","none"),
                              method = "methyl",
                              value_col = NULL,
                              scales = c("fixed","free_y","free_x","free"),
                              recist_color = FALSE,
                              scheme = c("three","two","four","four_alt","five","six"),
                              show_stats = .o("show_stats"),
                              show_n = .o("show_n"),
                              title = "Longitudinal ctDNA reduction by dose",
                              subtitle = "Thick line = median",
                              caption = sprintf("Floor=%g; LOQ=%g; baseline vs each timepoint per dose",
                                                .o("display_floor"), .o("loq")),
                              xlab = "Time point",
                              ylab = "Reduction ratio (log)",
                              point_size = .o("point_size"),
                              point_alpha = .o("point_alpha"),
    stat_position = .o("stat_position"),
    legend_position = .o("legend_position"),
    facet = NULL,
    y_scale = "log10",
                          indications = NULL) {
  # ---- cohort restriction ------------------------
  df <- .ctdna_filter_by_indication(df, indications)
  # -----------------------------------------------
  # v0.25: method accepts a vector for side-by-side comparison.
  if (length(method) > 1L || identical(method, "all")) {
    methods_resolved <- .resolve_methods(method, domain = "tf")
    plots <- lapply(methods_resolved, function(m) {
      ctdna_plot_longitudinal(df, stat = stat, method = m,
        value_col = if (m == "custom") value_col else NULL,
        scales = scales, recist_color = recist_color, scheme = scheme,
        show_stats = show_stats, show_n = show_n,
        title = sprintf("%s — %s", title, .resolve_method(m, "tf")$label),
        subtitle = subtitle, caption = caption, xlab = xlab, ylab = ylab,
        point_size = point_size, point_alpha = point_alpha,
        stat_position = stat_position, legend_position = legend_position,
        y_scale = y_scale)
    })
    return(patchwork::wrap_plots(plots, ncol = 1L))
  }
  stat <- match.arg(stat)
  scheme <- match.arg(scheme)
  scales <- match.arg(scales)
  m <- .resolve_method(method, "tf", value_col = value_col)
  d <- ctdna_ratio(df, tf_col = m$value_col)
  d$.ratio <- ctdna_floor_tf(d$ratio)
  d[[.o("dose")]] <- .ctdna_factor_canonical(d[[.o("dose")]], kind = "dose")
  d[[.o("time")]] <- .ctdna_factor_canonical(d[[.o("time")]], kind = "time_point")
  if (recist_color && .o("recist") %in% names(d)) {
    d$.r <- ctdna_stratify_recist(d[[.o("recist")]], scheme)
  }
  med <- stats::aggregate(.ratio ~ get(.o("dose")) + get(.o("time")),
                          data = d, FUN = stats::median)
  names(med) <- c(.o("dose"), .o("time"), ".m")

  pw <- NULL
  if (show_stats && stat != "none") {
    bl_lbl <- .o("baseline"); rows <- list()
    for (ds in levels(d[[.o("dose")]])) {
      sub <- d[d[[.o("dose")]] == ds, ]
      # subject-aligned: take only subjects with both BL and the timepoint
      bl_df <- sub[sub[[.o("time")]] == bl_lbl, c(.o("subject"), ".ratio")]
      names(bl_df)[2] <- ".bl"
      for (tp in setdiff(levels(d[[.o("time")]]), bl_lbl)) {
        ot_df <- sub[sub[[.o("time")]] == tp, c(.o("subject"), ".ratio")]
        names(ot_df)[2] <- ".ot"
        joined <- merge(bl_df, ot_df, by = .o("subject"))
        if (nrow(joined) < 2) next
        r <- ctdna_paired_test(joined$.bl, joined$.ot, method = stat)
        rows[[length(rows) + 1]] <- data.frame(
          dose = ds, time = tp, label = ctdna_pval_label(r$p),
          stringsAsFactors = FALSE)
      }
    }
    pw <- if (length(rows)) do.call(rbind, rows) else NULL
    if (!is.null(pw)) names(pw) <- c(.o("dose"), .o("time"), "label")
  }

  p <- if (recist_color && ".r" %in% names(d)) {
    ggplot2::ggplot(d, ggplot2::aes(.data[[.o("time")]], .data$.ratio,
                                     group = .data[[.o("subject")]],
                                     color = .data$.r)) +
      ggplot2::geom_line(alpha = 0.45) +
      ggplot2::geom_line(data = med,
                          ggplot2::aes(.data[[.o("time")]], .data$.m,
                                       group = .data[[.o("dose")]]),
                          color = .o("median_color"), linewidth = 1.5,
                          inherit.aes = FALSE) +
      ctdna_scale_recist() +
      ggplot2::labs(color = "RECIST")
  } else {
    ggplot2::ggplot(d, ggplot2::aes(.data[[.o("time")]], .data$.ratio,
                                     group = .data[[.o("subject")]])) +
      ggplot2::geom_line(alpha = 0.25, color = "grey50") +
      ggplot2::geom_line(data = med,
                          ggplot2::aes(.data[[.o("time")]], .data$.m,
                                       group = .data[[.o("dose")]]),
                          color = .o("median_color"), linewidth = 1.5,
                          inherit.aes = FALSE)
  }
  p <- p +
    ggplot2::facet_wrap(stats::as.formula(paste("~", .o("dose"))),
                         scales = scales) +
    .log_clearance(y_scale) + ctdna_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, vjust = 1, hjust = 1,
                                            size = 7))

  if (show_stats && !is.null(pw)) {
    p <- p + ggplot2::geom_text(data = pw,
                                ggplot2::aes(x = .data[[.o("time")]], y = -Inf,
                                             label = .data$label),
                                vjust = -1.5, size = 2.7,
                                color = "grey20",
                                inherit.aes = FALSE)
  }
  if (show_n) {
    # v0.43.8: per-(dose, time) counts as a visible white-box label inside
    # each panel. Anchored to the top edge (y = Inf, vjust = 1.2). x-axis
    # labels NOT modified here since x is time (varies per panel).
    counts <- stats::aggregate(get(.o("subject")) ~ get(.o("dose")) +
                                 get(.o("time")), data = d,
                               FUN = function(x) length(unique(x)))
    names(counts) <- c(.o("dose"), .o("time"), "n")
    counts$label <- paste0("n=", counts$n)
    p <- p + ggplot2::geom_label(
      data = counts,
      mapping = ggplot2::aes(x = .data[[.o("time")]], y = Inf,
                              label = .data$label),
      inherit.aes   = FALSE,
      vjust         = 1.2,
      size          = 2.8,
      label.padding = ggplot2::unit(0.12, "lines"),
      label.size    = 0,
      color         = "grey15",
      fill          = "white",
      fontface      = "bold")
  }
  .finalize(p, title, subtitle, caption, xlab, ylab, legend_position = legend_position)
}


# ---- D6: Longitudinal by dose AND RECIST ------------------------------------

#' D6: Longitudinal ctDNA reduction by dose AND RECIST
#'
#' Subject-level lines faceted by both dose (columns) and RECIST (rows).
#' Reports a per-RECIST overall dose-effect test in the subtitle.
#'
#' @section Recommended statistic:
#' `"kruskal"` (default): non-parametric test of dose effect within each
#' RECIST stratum. Use `"wilcox"` if dose has exactly two levels;
#' `"anova"` only if log-ratios are approximately normal.
#'
#' @param df Longitudinal data.
#' @param scheme RECIST stratification.
#' @param stat One of `"kruskal"` (default), `"wilcox"`, `"anova"`, `"none"`.
#'   Recommended: `"kruskal"`.
#' @param show_stats,show_n,title,subtitle,caption,xlab,ylab,point_size,point_alpha
#'   Display.
#' @param indications Optional character vector of Indication values (e.g.
#'   c("NSCLC","mCRPC")). When set, `df` is restricted to patients whose
#'   Indication is in this set. Lookup falls through Indication -> indication
#'   -> Cancertype in the frame itself; if not found in the frame, a Patient_ID
#'   join against prep$clinical is attempted. NULL (default) = no filter.
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 15, seed = 1)
#' prep <- ctdna_prepare(sim, verbose = FALSE)
#'
#' # Per-subject TF trajectories colored by RECIST (three-group scheme)
#' ctdna_plot_longitudinal_recist(prep$ctdna, scheme = "three")
#' @export
ctdna_plot_longitudinal_recist <- function(df,
                                     scheme = c("three","two","four","four_alt","five","six"),
                                     method = "methyl",
                                     value_col = NULL,
                                     stat = c("kruskal","wilcox","anova","none"),
                                     scales = c("fixed","free_y","free_x","free"),
                                     show_stats = .o("show_stats"),
                                     show_n = .o("show_n"),
                                     title = "ctDNA reduction by dose x RECIST",
                                     subtitle = NULL,
                                     caption = sprintf("Floor=%g; LOQ=%g",
                                                       .o("display_floor"), .o("loq")),
                                     xlab = "Time point",
                                     ylab = "Reduction ratio (log)",
                                     point_size = .o("point_size"),
                                     point_alpha = .o("point_alpha"),
    stat_position = .o("stat_position"),
    legend_position = .o("legend_position"),
    facet = NULL,
    y_scale = "log10",
                          indications = NULL) {
  # ---- cohort restriction ------------------------
  df <- .ctdna_filter_by_indication(df, indications)
  # -----------------------------------------------
  # v0.25: method accepts a vector for side-by-side comparison.
  if (length(method) > 1L || identical(method, "all")) {
    methods_resolved <- .resolve_methods(method, domain = "tf")
    plots <- lapply(methods_resolved, function(m) {
      ctdna_plot_longitudinal_recist(df, scheme = scheme, method = m,
        value_col = if (m == "custom") value_col else NULL,
        stat = stat, scales = scales,
        show_stats = show_stats, show_n = show_n,
        title = sprintf("%s — %s", title, .resolve_method(m, "tf")$label),
        subtitle = subtitle, caption = caption, xlab = xlab, ylab = ylab,
        point_size = point_size, point_alpha = point_alpha,
        stat_position = stat_position, legend_position = legend_position,
        y_scale = y_scale)
    })
    return(patchwork::wrap_plots(plots, ncol = 1L))
  }
  .stats_label <- NULL

  scheme <- match.arg(scheme); stat <- match.arg(stat)
  scales <- match.arg(scales)
  m <- .resolve_method(method, "tf", value_col = value_col)
  d <- ctdna_ratio(df, tf_col = m$value_col)
  d$.ratio <- ctdna_floor_tf(d$ratio)
  d$.r <- ctdna_stratify_recist(d[[.o("recist")]], scheme)
  d[[.o("dose")]] <- .ctdna_factor_canonical(d[[.o("dose")]], kind = "dose")
  d[[.o("time")]] <- .ctdna_factor_canonical(d[[.o("time")]], kind = "time_point")

  if (show_stats && stat != "none" ) {
    parts <- lapply(split(d, d$.r), function(s) {
      r <- ctdna_run_group_test(s, ".ratio", .o("dose"), method = stat)
      # Abbreviate test names (KW = Kruskal-Wallis, WX = Wilcoxon)
      tname <- switch(r$name,
                      "Kruskal-Wallis" = "KW",
                      "Wilcoxon"       = "WX",
                      "ANOVA"          = "ANOVA",
                      "t-test"         = "t",
                      r$name)
      sprintf("%s: %s %s", unique(s$.r), tname, ctdna_pval_label(r$p))
    })
    .stats_label <- paste(unlist(parts), collapse = " | ")
  }

  ggplot2::ggplot(d, ggplot2::aes(.data[[.o("time")]], .data$.ratio,
                                  group = .data[[.o("subject")]],
                                  color = .data$.r)) +
    ggplot2::geom_line(alpha = 0.45) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::facet_grid(stats::as.formula(paste(".r ~", .o("dose"))),
                         scales = scales) +
    ctdna_scale_recist() + .log_clearance(y_scale) + ctdna_theme() +
    ggplot2::labs(color = "RECIST", x = "Time point",
                   y = "Reduction ratio (log)") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, vjust = 1, hjust = 1,
                                            size = 7)) -> p
    # Route stats label per stat_position option
  if (!is.null(.stats_label)) {
    .routed <- .route_stats(.stats_label, stat_position, subtitle, caption)
    subtitle <- .routed$subtitle; caption <- .routed$caption
  
    on_plot_label <- .routed$on_plot
  } else on_plot_label <- NULL
  if (is.null(caption))
    caption <- sprintf("Dotted line = display floor (%g)", .o("display_floor"))
  # (on_plot_label annotated below after .finalize)
  .annotate_stat(.finalize(p, title, subtitle, caption, xlab, ylab, legend_position = legend_position), on_plot_label)
}


# ---- D7: Best reduction vs best % tumor change ------------------------------

#' D7: Best ctDNA reduction vs best % change in target lesions
#'
#' Scatter of per-subject best log-scale reduction ratio (x) vs the best
#' % change in sum of target lesions (y), with the -30% RECIST threshold.
#'
#' @section Recommended statistic:
#' `"spearman"` (default): preferred because ratios at clearance produce
#' many ties; Spearman handles ties robustly. Use `"pearson"` only if
#' both axes are approximately linear-Gaussian. `"kendall"` is an
#' alternative rank-based option.
#'
#' @param df Longitudinal data.
#' @param scheme RECIST stratification (used for shape encoding).
#' @param stat One of `"spearman"` (default), `"pearson"`, `"kendall"`,
#'   `"none"`. Recommended: `"spearman"`.
#' @param show_stats,title,subtitle,caption,xlab,ylab,point_size,point_alpha
#'   Display.
#' @param indications Optional character vector of Indication values (e.g.
#'   c("NSCLC","mCRPC")). When set, `df` is restricted to patients whose
#'   Indication is in this set. Lookup falls through Indication -> indication
#'   -> Cancertype in the frame itself; if not found in the frame, a Patient_ID
#'   join against prep$clinical is attempted. NULL (default) = no filter.
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 30, seed = 1)
#' prep <- ctdna_prepare(sim, verbose = FALSE)
#' # ctDNA reduction vs tumor-size change (RECIST sum of diameters)
#' ctdna_plot_reduction_vs_tumor(prep$ctdna, prep$tumor)
#' @export
ctdna_plot_reduction_vs_tumor <- function(df,
                                    scheme = c("three","two","four","four_alt","five","six"),
                                    method = "methyl",
                                    tumor_method = "recist_pct",
                                    value_col = NULL,
                                    stat = c("spearman","pearson","kendall","none"),
                                    show_stats = .o("show_stats"),
                                    title = "ctDNA reduction vs best % tumor change",
                                    subtitle = NULL,
                                    caption = "Red dashed = -30% (RECIST PR threshold)",
                                    xlab = NULL,
                                    ylab = NULL,
                                    point_size = .o("point_size") + 1,
                                    point_alpha = .o("point_alpha"),
    stat_position = .o("stat_position"),
    legend_position = .o("legend_position"),
    facet = NULL,
    y_scale = "linear",
                          indications = NULL) {
  # ---- cohort restriction ------------------------
  df <- .ctdna_filter_by_indication(df, indications)
  # -----------------------------------------------
  # v0.25: method accepts a vector for side-by-side comparison.
  if (length(method) > 1L || identical(method, "all")) {
    methods_resolved <- .resolve_methods(method, domain = "tf")
    plots <- lapply(methods_resolved, function(m) {
      ctdna_plot_reduction_vs_tumor(df, scheme = scheme, method = m,
        tumor_method = tumor_method,
        value_col = if (m == "custom") value_col else NULL,
        stat = stat, show_stats = show_stats,
        title = sprintf("%s — %s", title, .resolve_method(m, "tf")$label),
        subtitle = subtitle, caption = caption, xlab = xlab, ylab = ylab,
        point_size = point_size, point_alpha = point_alpha,
        stat_position = stat_position, legend_position = legend_position,
        y_scale = y_scale)
    })
    return(patchwork::wrap_plots(plots, ncol = 1L))
  }
  .stats_label <- NULL

  scheme <- match.arg(scheme); stat <- match.arg(stat)
  m_tf <- .resolve_method(method, "tf", value_col = value_col)
  m_t  <- .resolve_method(tumor_method, "tumor")
  if (is.null(xlab)) xlab <- paste0("Best reduction ratio (", m_tf$label, ", log)")
  if (is.null(ylab)) ylab <- m_t$label
  d <- ctdna_summary(df)
  d$.ratio <- ctdna_floor_tf(d$best_ratio)
  d$.r <- ctdna_stratify_recist(d[[.o("recist")]], scheme)
  d[[.o("dose")]] <- .ctdna_factor_canonical(d[[.o("dose")]], kind = "dose")
  if (show_stats && stat != "none" ) {
    ct <- ctdna_cor_test(d$.ratio, d[[.o("best_change")]], method = stat)
    .stats_label <- sprintf("%s rho = %.2f, %s, n = %d",
                        ct$name, ct$rho, ctdna_pval_label(ct$p), ct$n)
  }

  ggplot2::ggplot(d, ggplot2::aes(.data$.ratio,
                                  .data[[.o("best_change")]],
                                  color = .data$.r,
                                  shape = .data[[.o("dose")]])) +
    ggplot2::geom_hline(yintercept = -30, color = "red", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = .o("display_floor"), color = .o("clearance_color"),
                        linetype = "dotted") +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::scale_x_log10() +
    ctdna_scale_recist() +
    ctdna_scale_dose_shape() +
    ggplot2::labs(color = "RECIST", shape = .o("dose")) +
    .ctdna_resolve_y_scale(y_scale) +
    ctdna_theme() -> p
    # Route stats label per stat_position option
  if (!is.null(.stats_label)) {
    .routed <- .route_stats(.stats_label, stat_position, subtitle, caption)
    subtitle <- .routed$subtitle; caption <- .routed$caption
  
    on_plot_label <- .routed$on_plot
  } else on_plot_label <- NULL
  # (on_plot_label annotated below after .finalize)
  .annotate_stat(.finalize(p, title, subtitle, caption, xlab, ylab, legend_position = legend_position), on_plot_label)
}


# ---- D8 supplemental: maxVAF vs methylTF -------------------------------------

#' D8: maxVAF vs methylTF concordance
#'
#' Log-log scatter with identity line. Flags points beyond a configurable
#' log10 distance (default >1, i.e. >10-fold) as discordant.
#'
#' @section Recommended statistic:
#' `"spearman"` (default): robust to ties caused by flooring.
#' `"pearson"` works if both axes are normal on the log scale.
#'
#' @param df Data with maxVAF and methylTF.
#' @param stat One of `"spearman"` (default), `"pearson"`, `"kendall"`,
#'   `"none"`. Recommended: `"spearman"`.
#' @param show_stats,title,subtitle,caption,xlab,ylab,point_size,point_alpha
#'   Display.
#' @param indications Optional character vector of Indication values (e.g.
#'   c("NSCLC","mCRPC")). When set, `df` is restricted to patients whose
#'   Indication is in this set. Lookup falls through Indication -> indication
#'   -> Cancertype in the frame itself; if not found in the frame, a Patient_ID
#'   join against prep$clinical is attempted. NULL (default) = no filter.
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 30, seed = 1)
#' prep <- ctdna_prepare(sim, verbose = FALSE)
#' # Per-subject mutation TF vs methylation TF
#' ctdna_plot_mut_methyl(prep$ctdna)
#' @export
ctdna_plot_mut_methyl <- function(df,
                            stat = c("spearman","pearson","kendall","none"),
                            recist_shape = FALSE,
                            scheme = c("three","two","four","four_alt","five","six"),
                            show_stats = .o("show_stats"),
                            title = "maxVAF vs methylTF concordance",
                            subtitle = NULL,
                            caption = sprintf("Discordant if log10 distance > %g",
                                              .o("discord_log10")),
                            xlab = "maxVAF (log)",
                            ylab = "methylTF (log)",
                            point_size = .o("point_size"),
                            point_alpha = .o("point_alpha"),
    stat_position = .o("stat_position"),
    legend_position = .o("legend_position"),
    facet = NULL,
    y_scale = "log10",
                          indications = NULL) {
  # ---- cohort restriction ------------------------
  df <- .ctdna_filter_by_indication(df, indications)
  # -----------------------------------------------
  .stats_label <- NULL

  stat <- match.arg(stat)
  scheme <- match.arg(scheme)
  d <- data.frame(maxvaf = ctdna_floor_tf(df[[.o("maxvaf")]]),
                  tf    = ctdna_floor_tf(df[[.o("tf")]]))
  d$.disc <- abs(log10(d$maxvaf) - log10(d$tf)) > .o("discord_log10")
  if (recist_shape && .o("recist") %in% names(df)) {
    d$.r <- ctdna_stratify_recist(df[[.o("recist")]], scheme)
  }
  n_disc <- sum(d$.disc, na.rm = TRUE)
  if (show_stats && stat != "none" ) {
    ct <- ctdna_cor_test(d$maxvaf, d$tf, method = stat)
    .stats_label <- sprintf("%s rho = %.2f, %s, n = %d; %d discordant",
                        ct$name, ct$rho, ctdna_pval_label(ct$p), ct$n, n_disc)
  }
  p <- if (recist_shape && ".r" %in% names(d)) {
    ggplot2::ggplot(d, ggplot2::aes(.data$maxvaf, .data$tf,
                                     color = .data$.disc,
                                     shape = .data$.r)) +
      ggplot2::geom_abline(slope = 1, color = "grey50", linetype = "dashed") +
      ggplot2::geom_point(size = point_size + 0.5, alpha = point_alpha) +
      ggplot2::scale_color_manual(values = c(`FALSE` = "grey20",
                                              `TRUE` = "#E31A1C"),
                                   labels = c("concordant", "discordant"),
                                   name = NULL) +
      ggplot2::scale_shape_manual(values = c("CR/PR" = 16, "uCR/PR" = 17,
                                              "SD" = 15, "PD/NE/NA" = 18,
                                              "R" = 16, "NR" = 15,
                                              "SD/PD/NE/NA" = 15),
                                   na.value = 4, name = "RECIST")
  } else {
    ggplot2::ggplot(d, ggplot2::aes(.data$maxvaf, .data$tf,
                                     color = .data$.disc)) +
      ggplot2::geom_abline(slope = 1, color = "grey50", linetype = "dashed") +
      ggplot2::geom_point(size = point_size, alpha = point_alpha) +
      ggplot2::scale_color_manual(values = c(`FALSE` = "grey20",
                                              `TRUE` = "#E31A1C"),
                                   labels = c("concordant", "discordant"),
                                   name = NULL)
  }
  p <- p + ggplot2::scale_x_log10() + .ctdna_resolve_y_scale(y_scale) + ctdna_theme()
    # Route stats label per stat_position option
  if (!is.null(.stats_label)) {
    .routed <- .route_stats(.stats_label, stat_position, subtitle, caption)
    subtitle <- .routed$subtitle; caption <- .routed$caption
  
    on_plot_label <- .routed$on_plot
  } else on_plot_label <- NULL
  # (on_plot_label annotated below after .finalize)
  .annotate_stat(.finalize(p, title, subtitle, caption, xlab, ylab, legend_position = legend_position), on_plot_label)
}

# ----------------------------------------------------------------------------
# v0.38.0: ctdna_plot_pct_change_by_dose_time
# ----------------------------------------------------------------------------

#' Percent-change boxplots by dose and time, grouped by RECIST
#'
#' Wrapped boxplots showing percent change in tumor burden (default
#' column \code{Ctdna_percentage_change}) stratified by dose (x within
#' each panel), time point (panel columns), and RECIST grouping (panel
#' rows). Empty time x RECIST cells are kept as empty rectangles via
#' \code{drop = FALSE} so the grid stays rectangular.
#'
#' @param df Long data frame with required columns (see arguments).
#' @param time_col Column for the time axis. Default \code{"time_point"}.
#' @param dose_col Column for dose (x within each panel). Default
#'   \code{ctdna_opts("dose")}.
#' @param recist_col Column for RECIST values. Default
#'   \code{ctdna_opts("recist")}.
#' @param value_col Column for the y-axis (percent change). Default
#'   \code{"Ctdna_percentage_change"}.
#' @param recist_grouping One of \code{"three"} (default), \code{"two"}
#'   (R/NR), \code{"four"}, \code{"four_alt"}, \code{"five"},
#'   \code{"six"}; see \code{\link{ctdna_stratify_recist}}.
#' @param show_points Logical. Overlay individual subject points jittered
#'   on top of each box. Default \code{TRUE}.
#' @param show_stats Logical. Annotate per-cell p-values (Kruskal-Wallis
#'   across doses within each cell). Default \code{FALSE} (opt-in).
#' @param y_scale One of \code{"linear"} (default), \code{"log10"},
#'   \code{"log2"}, \code{"log"}, or \code{"log_<N>"}.
#' @param title,subtitle,caption,xlab,ylab Plot labels (auto-derived if NULL).
#' @param point_size,point_alpha Display for jittered points.
#' @param legend_position One of \code{"right"} (default), \code{"left"},
#'   \code{"top"}, \code{"bottom"}, \code{"none"}.
#' @param indications Optional character vector of Indication values (e.g.
#'   c("NSCLC","mCRPC")). When set, `df` is restricted to patients whose
#'   Indication is in this set. Lookup falls through Indication -> indication
#'   -> Cancertype in the frame itself; if not found in the frame, a Patient_ID
#'   join against prep$clinical is attempted. NULL (default) = no filter.
#' @return A ggplot object.
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 30, seed = 1)
#' prep <- ctdna_prepare(sim, verbose = FALSE)
#' d <- prep$ctdna
#' d$Ctdna_percentage_change <- 100 * (d$methylTF / d$methylTF[1] - 1)
#'
#' # Percent change by dose, faceted by visit x RECIST scheme
#' ctdna_plot_pct_change_by_dose_time(d, recist_grouping = "three",
#'                                     show_stats = TRUE)
#' @export
ctdna_plot_pct_change_by_dose_time <- function(
    df,
    time_col = .o("time"),
    dose_col = .o("dose"),
    recist_col = .o("recist"),
    value_col = "Ctdna_percentage_change",
    recist_grouping = c("three","two","four","four_alt","five","six"),
    show_points = TRUE,
    show_stats = FALSE,
    y_scale = "linear",
    y_clip = c(-100, 200),
    title = NULL, subtitle = NULL, caption = NULL,
    xlab = NULL, ylab = NULL,
    point_size = .o("point_size"),
    point_alpha = .o("point_alpha"),
    legend_position = .o("legend_position"),
                          indications = NULL) {
  # ---- cohort restriction ------------------------
  df <- .ctdna_filter_by_indication(df, indications)
  # -----------------------------------------------
  if (!is.data.frame(df))
    stop("ctdna_plot_pct_change_by_dose_time: `df` must be a data frame.",
          call. = FALSE)
  recist_grouping <- match.arg(recist_grouping)

  req <- c(time_col, dose_col, recist_col, value_col)
  miss <- setdiff(req, names(df))
  if (length(miss) > 0L)
    stop("ctdna_plot_pct_change_by_dose_time: missing column(s): ",
          paste(miss, collapse = ", "), ".", call. = FALSE)

  d <- df
  d$.r <- ctdna_stratify_recist(d[[recist_col]],
                                  recist_grouping = recist_grouping)
  d$.t <- factor(d[[time_col]])
  d$.d <- factor(d[[dose_col]])
  d$.y <- suppressWarnings(as.numeric(d[[value_col]]))
  d <- d[!is.na(d$.y), , drop = FALSE]

  p <- ggplot2::ggplot(d,
                        ggplot2::aes(x = .data$.d, y = .data$.y,
                                      fill = .data$.d)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.75)

  if (isTRUE(show_points)) {
    p <- p + ggplot2::geom_jitter(width = 0.15, height = 0,
                                    size = point_size,
                                    alpha = point_alpha,
                                    color = "grey30",
                                    na.rm = TRUE)
  }

  p <- p +
    ggplot2::facet_grid(.r ~ .t, drop = FALSE) +
    .ctdna_resolve_y_scale(y_scale) +
    ctdna_theme() +
    ggplot2::theme(legend.position = legend_position,
                    axis.text.x = ggplot2::element_text(angle = 45,
                                                          hjust = 1))

  # v0.42.0: per-facet (time x recist) bracket renderer. show_stats is
  # strictly TRUE/FALSE; brackets compare dose distributions within each
  # cell; overall test type + first-cell N go in the caption.
  if (isTRUE(show_stats)) {
    cells <- unique(d[, c(".t",".r")])
    facet_stats <- vector("list", nrow(cells))
    for (i in seq_len(nrow(cells))) {
      sub <- d[d$.t == cells$.t[i] & d$.r == cells$.r[i], , drop = FALSE]
      sub$.d <- droplevels(sub$.d)
      facet_stats[[i]] <- list(
        cell  = cells[i, ],
        sub   = sub,
        stats = .compute_box_stats(sub, group_col = ".d", value_col = ".y"))
    }
    # Add brackets per facet
    # v0.42.4: pass facet_data so each bracket layer lives in its
    # own (.r, .t) panel — no cross-panel leak from facet_grid.
    for (fs in facet_stats) {
      if (is.null(fs$stats)) next
      p <- .add_brackets(p, fs$stats, fs$sub,
                          group_col  = ".d", value_col = ".y",
                          y_scale    = y_scale,
                          facet_data = list(.r = fs$cell$.r,
                                              .t = fs$cell$.t))
    }
    first <- Find(function(x) !is.null(x$stats), facet_stats)
    if (!is.null(first))
      caption <- .compose_stats_caption(first$stats, "box",
                                          base_caption = caption)
  }

  # v0.42.2: cap the visible y range. Percent-change plots can have
  # extreme treatment-failure outliers (+5000%, +9900% etc) that
  # otherwise compress the data of interest into invisibility AND
  # cause the y-axis to read 10,000% or worse. `y_clip = c(lo, hi)`
  # is applied via coord_cartesian (clipping is visual only — the
  # data and stats are computed on the unclipped values). Pass
  # `y_clip = NULL` to render the full unclipped range.
  if (!is.null(y_clip) && length(y_clip) == 2L)
    p <- p + ggplot2::coord_cartesian(ylim = y_clip)

  .finalize(p,
             title  %||% "Percent change by dose and time",
             subtitle, caption,
             xlab %||% dose_col,
             ylab %||% value_col,
             legend_position = legend_position)
}


# ============================================================================
# ctdna_boxplot — UNIFIED grouped-boxplot function (added v0.42.5)
# ============================================================================
# Combines the functionality of D1 + D2 + D3 into a single configurable
# function. The existing D1/D2/D3 stay as named entry points (their
# defaults match the canonical deliverables exactly); ctdna_boxplot is
# the flexible alternative for any other Dose/RECIST/method/metric
# combination users want.

#' Unified ctDNA boxplot: configurable x / colour / facet / compare_by.
#'
#' Single grouped-boxplot function covering the structural space spanned
#' by \code{\link{ctdna_plot_baseline}} (D1),
#' \code{\link{ctdna_plot_baseline_dose_recist}} (D2), and the
#' between-dose comparisons used in \code{\link{ctdna_plot_reduction}}
#' (D3). Use the dedicated D1/D2/D3 functions when their defaults match
#' your needs; use \code{ctdna_boxplot} when you need to vary
#' \code{x} / \code{color} / \code{facet} / \code{compare_by}
#' independently.
#'
#' @section Reaching D1, D2, D3 from \code{ctdna_boxplot}:
#' \itemize{
#'   \item \code{ctdna_boxplot(d, x = "RECIST")} — D1 default (RECIST on
#'     x, one facet per TF method).
#'   \item \code{ctdna_boxplot(d, x = "Dose", color = "RECIST", facet = NULL, method = "methyl")} —
#'     D2 default (Dose on x, RECIST as colour, single panel).
#'   \item \code{ctdna_boxplot(d, x = "Dose", color = "RECIST", facet = NULL, value = "best_pct_change")} —
#'     D3-shape on a pre-computed \code{best_pct_change} column.
#' }
#'
#' @param df Long ctDNA data frame (typically from
#'   \code{ctdna_prepare()$ctdna}).
#'
#' @param x Column name for the x-axis. Defaults to \code{"Dose"}. Can
#'   be any column on \code{df}. Use \code{"RECIST"} for D1's default
#'   shape.
#'
#' @param color Column name to split boxes by (colour + fill + dodge).
#'   \code{NULL} (default) draws a single set of boxes. Use
#'   \code{"RECIST"} for D2's colour split.
#'
#' @param facet Faceting strategy. Three forms:
#' \itemize{
#'   \item \code{"method"} (default): when \code{value = "tf"} and
#'     multiple TF methods are requested via \code{method}, produces
#'     one panel per method.
#'   \item \code{NULL}: single panel.
#'   \item Any column name on \code{df}: facets by that column.
#' }
#'
#' @param value What to plot on the y-axis. Two forms:
#' \itemize{
#'   \item \code{"tf"} (default): use TF columns selected via
#'     \code{method}.
#'   \item A literal column name on \code{df}, e.g.
#'     \code{"best_pct_change"} or \code{"best_ratio"}, for
#'     reduction-style plots.
#' }
#'
#' @param method For \code{value = "tf"} only: which TF metric(s) to
#'   plot. Accepts \code{"all"} (default; expands to methyl + maxvaf +
#'   meanvaf), a single method (\code{"methyl"}, \code{"maxvaf"},
#'   \code{"meanvaf"}, \code{"custom"}), or a character vector. Ignored
#'   when \code{value} is a literal column.
#'
#' @param scheme RECIST collapse scheme; applied when \code{x} or
#'   \code{color} resolves to the RECIST column. One of \code{"two"},
#'   \code{"three"} (default), \code{"four"}, \code{"four_alt"},
#'   \code{"five"}, \code{"six"}. See
#'   \code{\link{ctdna_stratify_recist}}.
#'
#' @param compare_by Statistical comparison strategy. Determines what
#'   the p-value brackets compare:
#' \itemize{
#'   \item \code{"x"} (default): pairwise comparisons between
#'     \code{x} groups, pooling colour groups within each \code{x}.
#'     ONE bracket stack per panel. D1's and D2's default behaviour.
#'   \item \code{"color"}: within each \code{x} level, pairwise
#'     comparisons between colour groups. ONE stack per \code{x}
#'     level. Requires \code{color} to be set.
#'   \item \code{"x_within_color"}: turns \code{color} into a facet
#'     dimension; pairwise comparisons between \code{x} groups within
#'     each colour value. ONE stack per colour-facet. If \code{facet}
#'     is also set, you get a \code{facet_grid(facet ~ color)} layout.
#'     Requires \code{color} to be set.
#' }
#'
#' @param show_stats Logical. Compute and render p-value brackets +
#'   caption. Default \code{TRUE}.
#'
#' @param show_n Logical. Render \code{n=X} labels inside the bottom
#'   of each panel. Default \code{TRUE}.
#'
#' @param title,subtitle,caption,xlab,ylab Plot labels. \code{NULL} =
#'   auto-derived from the data.
#'
#' @param y_scale One of \code{"log10"} (default), \code{"log2"},
#'   \code{"log"} (natural), or \code{"linear"}.
#'
#' @param ylim Numeric length-2 vector for
#'   \code{coord_cartesian(ylim = ...)} to clip the visible y range
#'   without dropping data. Default \code{NULL} (auto-fit).
#'
#' @param indications Optional character vector of Indication values (e.g.
#'   c("NSCLC","mCRPC")). When set, `df` is restricted to patients whose
#'   Indication is in this set. Lookup falls through Indication -> indication
#'   -> Cancertype in the frame itself; if not found in the frame, a Patient_ID
#'   join against prep$clinical is attempted. NULL (default) = no filter.
#' @return A \code{ggplot} object.
#'
#' @section Input validation:
#' Reuse of the same column in two structural roles produces broken or
#' redundant plots and is rejected upfront:
#' \itemize{
#'   \item \code{color == x}: colour split would be redundant with the
#'     x-axis.
#'   \item \code{facet == x}: each facet would contain only one x
#'     value (empty boxplots at other positions).
#'   \item \code{facet == color}: facets would be monochromatic; the
#'     colour split is redundant.
#' }
#'
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 30, seed = 1)
#' prep <- ctdna_prepare(sim, verbose = FALSE)
#' d    <- prep$ctdna[prep$ctdna$Visit_name == ctdna_opts("baseline"), ]
#'
#' # D1-default: by RECIST, one panel per TF method
#' ctdna_boxplot(d, x = "RECIST")
#'
#' # D2-default: Dose on x, RECIST as colour, single panel
#' ctdna_boxplot(d, x = "Dose", color = "RECIST", facet = NULL,
#'               method = "methyl")
#'
#' # D3-shape on a pre-computed column
#' ctdna_boxplot(d, value = "best_pct_change", color = "RECIST",
#'               facet = NULL, y_scale = "linear", ylim = c(-100, 200))
#'
#' # Compare colours within each x rather than across x's
#' ctdna_boxplot(d, x = "Dose", color = "RECIST", facet = NULL,
#'               compare_by = "color")
#'
#' # Auto-facet by colour, compare x's within each facet
#' ctdna_boxplot(d, x = "Dose", color = "RECIST", facet = NULL,
#'               compare_by = "x_within_color")
#'
#' @export
ctdna_boxplot <- function(df,
                           x          = "Dose",
                           color      = NULL,
                           facet      = "method",
                           value      = "tf",
                           method     = "all",
                           scheme     = "three",
                           compare_by = c("x", "color", "x_within_color"),
                           show_stats = .o("show_stats"),
                           show_n     = .o("show_n"),
                           title      = NULL,
                           subtitle   = NULL,
                           caption    = NULL,
                           xlab       = NULL,
                           ylab       = NULL,
                           y_scale    = "log10",
                           ylim       = NULL,
                           point_size = .o("point_size"),
                           point_alpha = .o("point_alpha"),
                           legend_position = .o("legend_position"),
                          indications = NULL) {

  # ---- cohort restriction ------------------------
  df <- .ctdna_filter_by_indication(df, indications)
  # -----------------------------------------------
  compare_by <- match.arg(compare_by)
  if (compare_by != "x" && is.null(color))
    stop("compare_by = '", compare_by,
         "' requires `color` to be set.", call. = FALSE)

  # ----- input validation: reuse-of-column guards ------------------
  if (!is.null(color) && identical(color, x))
    stop("`color` and `x` are both '", x, "'. The colour split would be ",
          "redundant with the x-axis. Use a different column for `color`, ",
          "or set `color = NULL`.", call. = FALSE)
  if (!is.null(facet) && !identical(facet, "method") &&
      identical(facet, x))
    stop("`facet` and `x` are both '", x, "'. Each facet would contain ",
          "only one x value, leaving empty boxplots. Pick a different ",
          "column for `facet`, set `facet = NULL`, or use ",
          "`facet = \"method\"` for per-method panels.",
          call. = FALSE)
  if (!is.null(color) && !is.null(facet) && !identical(facet, "method") &&
      identical(facet, color))
    stop("`facet` and `color` are both '", color, "'. The colour split ",
          "would be redundant with the facets. Pick a different column ",
          "for one of them.", call. = FALSE)

  # ----- resolve x, color columns ----------------------------------
  x <- .resolve_col(x, df, what = "x")
  if (!is.null(color))
    color <- .resolve_col(color, df, what = "color")

  # ----- resolve value & method -----------------------------------
  if (value == "tf") {
    methods_resolved <- .resolve_methods(method, domain = "tf")
    tf_panels <- setNames(
      lapply(methods_resolved, function(mm) .resolve_method(mm, "tf")),
      methods_resolved)
    value_cols    <- unname(vapply(tf_panels, `[[`, character(1), "value_col"))
    display_units <- unname(vapply(tf_panels,
                                      function(e) e$display_unit %||% "fraction",
                                      character(1)))
    dominant_unit <- if (length(unique(display_units)) == 1L)
                          unique(display_units) else "fraction"
    panel_labels  <- setNames(vapply(tf_panels, `[[`, character(1), "label"),
                                value_cols)
  } else {
    if (!value %in% names(df))
      stop("`value = \"", value, "\"` is not a column on df.", call. = FALSE)
    value_cols <- value
    panel_labels <- setNames(value, value)
    # raw_pct = column already in percent units (e.g. best_pct_change).
    # No multiplication by 100 in the y-axis labeller.
    dominant_unit <- if (grepl("pct_change|pct_chg|percent_change|^pct$|_pct$|_pct_",
                                  value, ignore.case = TRUE)) "raw_pct"
                      else if (grepl("ratio", value, ignore.case = TRUE)) "actual"
                      else "fraction"
    tf_panels <- NULL
  }

  # ----- RECIST stratification ------------------------------------
  if (identical(x, "RECIST"))
    df$RECIST <- ctdna_stratify_recist(df$RECIST, recist_grouping = scheme)
  if (identical(color, "RECIST"))
    df$RECIST <- ctdna_stratify_recist(df$RECIST, recist_grouping = scheme)

  # ----- resolve facet --------------------------------------------
  if (is.null(facet)) {
    facet_col <- NULL
  } else if (identical(facet, "method")) {
    facet_col <- if (length(value_cols) > 1L) "variable" else NULL
  } else {
    facet_col <- .resolve_col(facet, df, what = "facet")
  }

  # ----- compare_by = "x_within_color": color becomes a facet -----
  using_facet_grid <- FALSE
  facet_col_2     <- NULL
  if (compare_by == "x_within_color") {
    if (is.null(facet_col)) {
      facet_col <- color
    } else {
      facet_col_2 <- color
      using_facet_grid <- TRUE
    }
  }

  # ----- reshape to long ------------------------------------------
  needed <- unique(c(x, color, facet_col, facet_col_2, value_cols))
  needed <- intersect(needed, names(df))
  if (length(value_cols) > 1L) {
    long <- tidyr::pivot_longer(df[, needed, drop = FALSE],
                                  cols = tidyr::all_of(value_cols),
                                  names_to = "variable",
                                  values_to = ".value")
  } else {
    long <- df[, needed, drop = FALSE]
    long$.value <- df[[value_cols]]
    long$variable <- value_cols
  }
  long <- long[!is.na(long$.value), ]

  # ----- floor TF values ------------------------------------------
  if (value == "tf") {
    long$.value <- ifelse(long$variable %in% c(value_cols, "methylTF",
                                                   "maxVAF", "meanVAF"),
                            ctdna_floor_tf(long$.value),
                            long$.value)
  }

  # ----- canonical factor orders ----------------------------------
  long[[x]] <- .ctdna_factor_canonical(long[[x]], kind = x)
  if (!is.null(color))
    long[[color]] <- .ctdna_factor_canonical(long[[color]], kind = color)
  if (!is.null(facet_col) && !identical(facet_col, "variable") &&
      !identical(facet_col, color))
    long[[facet_col]] <- .ctdna_factor_canonical(long[[facet_col]],
                                                    kind = facet_col)

  # ----- stats per compare_by mode --------------------------------
  iter_facets <- function(long_df) {
    if (is.null(facet_col) && is.null(facet_col_2)) {
      list(list(key = NULL, sub = long_df))
    } else if (using_facet_grid) {
      combos <- unique(long_df[, c(facet_col, facet_col_2), drop = FALSE])
      lapply(seq_len(nrow(combos)), function(i) {
        row <- combos[i, ]
        sub <- long_df[long_df[[facet_col]]   == row[[facet_col]] &
                        long_df[[facet_col_2]] == row[[facet_col_2]], ]
        k <- setNames(list(row[[facet_col]], row[[facet_col_2]]),
                        c(facet_col, facet_col_2))
        list(key = k, sub = sub)
      })
    } else {
      lapply(unique(long_df[[facet_col]]), function(fv) {
        sub <- long_df[long_df[[facet_col]] == fv, ]
        list(key = setNames(list(fv), facet_col), sub = sub)
      })
    }
  }

  stats_entries <- list()
  if (isTRUE(show_stats)) {
    if (compare_by == "x") {
      for (fc in iter_facets(long)) {
        s <- .compute_box_stats(fc$sub, group_col = x, value_col = ".value")
        if (!is.null(s))
          stats_entries[[length(stats_entries) + 1L]] <- list(
            stats = s, facet_key = fc$key, x_lookup = NULL, x_val = NULL)
      }
    } else if (compare_by == "x_within_color") {
      for (fc in iter_facets(long)) {
        s <- .compute_box_stats(fc$sub, group_col = x, value_col = ".value")
        if (!is.null(s))
          stats_entries[[length(stats_entries) + 1L]] <- list(
            stats = s, facet_key = fc$key, x_lookup = NULL, x_val = NULL)
      }
    } else if (compare_by == "color") {
      K       <- length(levels(long[[color]]))
      dodge_w <- 0.8
      offsets <- (seq_len(K) - (K + 1)/2) * (dodge_w / K)
      names(offsets) <- levels(long[[color]])
      for (fc in iter_facets(long)) {
        x_levels <- levels(droplevels(fc$sub[[x]]))
        for (xv in x_levels) {
          sub <- fc$sub[fc$sub[[x]] == xv, ]
          if (!nrow(sub)) next
          s <- .compute_box_stats(sub, group_col = color, value_col = ".value")
          if (is.null(s)) next
          x_pos <- which(levels(long[[x]]) == xv)
          xl    <- x_pos + offsets[s$group_levels]
          stats_entries[[length(stats_entries) + 1L]] <- list(
            stats = s, facet_key = fc$key, x_lookup = xl, x_val = xv)
        }
      }
    }
  }

  # ----- build base plot -----------------------------------------
  if (is.null(color)) {
    p <- ggplot2::ggplot(long, ggplot2::aes(.data[[x]], .data$.value,
                                              color = .data[[x]],
                                              fill  = .data[[x]])) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.3,
                              linewidth = 0.7, fatten = 2.2,
                              width = .o("box_width")) +
      ggplot2::geom_jitter(width = 0.15, size = point_size,
                             alpha = point_alpha,
                             shape = 21, stroke = 0.3, color = "grey20")
  } else if (compare_by == "x_within_color") {
    p <- ggplot2::ggplot(long, ggplot2::aes(.data[[x]], .data$.value,
                                              color = .data[[color]],
                                              fill  = .data[[color]])) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.3,
                              linewidth = 0.7, fatten = 2.2,
                              width = .o("box_width")) +
      ggplot2::geom_jitter(width = 0.15, size = point_size,
                             alpha = point_alpha,
                             shape = 21, stroke = 0.3, color = "grey20")
  } else {
    p <- ggplot2::ggplot(long, ggplot2::aes(.data[[x]], .data$.value,
                                              color = .data[[color]],
                                              fill  = .data[[color]])) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.3,
                              linewidth = 0.7, fatten = 2.2,
                              width = .o("box_width"),
                              position = ggplot2::position_dodge(0.8)) +
      ggplot2::geom_point(position = ggplot2::position_jitterdodge(
                              jitter.width = 0.15, dodge.width = 0.8),
                            size = point_size, alpha = point_alpha,
                            shape = 21, stroke = 0.3, color = "grey20")
  }

  # ----- y-scale with appropriate labeller ------------------------
  pct_x100 <- function(x) paste0(formatC(x * 100, digits = 3, format = "g"), "%")
  pct_raw  <- function(x) paste0(formatC(x, digits = 3, format = "g"), "%")
  scale_y_args <- list(y_scale, expand = ggplot2::expansion(mult = c(0.05, 0.03)))
  if (dominant_unit == "percent")      scale_y_args$labels <- pct_x100
  else if (dominant_unit == "raw_pct") scale_y_args$labels <- pct_raw
  p <- p + do.call(.ctdna_resolve_y_scale, scale_y_args)

  # ----- colour scale --------------------------------------------
  pal_key <- color %||% x
  if (identical(pal_key, "RECIST")) {
    p <- p + ctdna_scale_recist(aesthetics = c("color", "fill"))
  } else {
    p <- p + ctdna_scale_dose(aesthetics = c("color", "fill"))
  }

  # ----- LOQ floor line for TF data -------------------------------
  if (dominant_unit %in% c("percent", "fraction") && value == "tf")
    p <- p + ggplot2::geom_hline(yintercept = .o("display_floor"),
                                    linetype = "dotted",
                                    color = .o("clearance_color"),
                                    linewidth = 0.4)

  # ----- facet with deterministic labeller ------------------------
  facet_labeller_fn <- function(labels) {
    lapply(labels, function(col) {
      vals   <- as.character(col)
      mapped <- panel_labels[vals]
      ifelse(is.na(mapped), vals, unname(mapped))
    })
  }
  if (using_facet_grid) {
    p <- p + ggplot2::facet_grid(rows = ggplot2::vars(.data[[facet_col]]),
                                    cols = ggplot2::vars(.data[[facet_col_2]]),
                                    scales = "free_y",
                                    labeller = facet_labeller_fn)
  } else if (!is.null(facet_col)) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[[facet_col]]),
                                    scales = "free_y",
                                    labeller = facet_labeller_fn)
  }

  p <- p + ctdna_theme() +
    ggplot2::labs(x     = x,
                    color = if (is.null(color)) x else color,
                    fill  = if (is.null(color)) x else color)

  # ----- brackets -------------------------------------------------
  if (length(stats_entries)) {
    for (e in stats_entries) {
      sub_df <- long
      if (!is.null(e$facet_key))
        for (k in names(e$facet_key))
          sub_df <- sub_df[sub_df[[k]] == e$facet_key[[k]], , drop = FALSE]
      if (compare_by == "color" && !is.null(e$x_val)) {
        sub_df <- sub_df[sub_df[[x]] == e$x_val, , drop = FALSE]
        group_col_b <- color
      } else {
        group_col_b <- x
      }
      p <- .add_brackets(p, e$stats, sub_df,
                          group_col  = group_col_b,
                          value_col  = ".value",
                          y_scale    = y_scale,
                          x_lookup   = e$x_lookup,
                          facet_data = e$facet_key)
    }
  }

  # ----- n-labels: x-axis tick labels + visible inside-panel boxes -----
  if (show_n) {
    # v0.43.8: per-(facet, x) counts as inside-panel white-box labels,
    # plus overall n per x-group appended to x-axis tick labels.
    if (using_facet_grid) {
      pfc <- as.data.frame(table(long[[facet_col]],
                                   long[[facet_col_2]],
                                   long[[x]]))
      names(pfc) <- c(facet_col, facet_col_2, "g", "n")
    } else if (!is.null(facet_col)) {
      pfc <- as.data.frame(table(long[[facet_col]], long[[x]]))
      names(pfc) <- c(facet_col, "g", "n")
    } else {
      pfc <- as.data.frame(table(long[[x]]))
      names(pfc) <- c("g", "n")
    }
    pfc <- pfc[pfc$n > 0, , drop = FALSE]
    pfc$g <- factor(pfc$g, levels = levels(long[[x]]))

    # (a) X-axis tick labels — overall n per x-group from `long`
    n_overall_tab <- table(long[[x]][!is.na(long[[x]])])
    n_overall <- setNames(as.integer(n_overall_tab), names(n_overall_tab))
    p <- p + ggplot2::scale_x_discrete(labels = function(z) {
      n <- n_overall[as.character(z)]
      ifelse(is.na(n), as.character(z),
              sprintf("%s\n(n=%d)", as.character(z), n))
    })

    # (b) Per-facet inside-panel labels with visible styling
    if (nrow(pfc) > 0L) {
      p <- p + ggplot2::geom_label(
        data = pfc,
        mapping = ggplot2::aes(x = .data$g, y = -Inf,
                                 label = paste0("n=", .data$n)),
        inherit.aes   = FALSE,
        vjust         = -0.2,
        size          = 2.8,
        label.padding = ggplot2::unit(0.12, "lines"),
        label.size    = 0,
        color         = "grey15",
        fill          = "white",
        fontface      = "bold")
    }
  }

  # ----- coord_cartesian(clip="off") + optional ylim --------------
  cc_args <- list(clip = "off")
  if (!is.null(ylim)) cc_args$ylim <- ylim
  p <- p + do.call(ggplot2::coord_cartesian, cc_args)

  # ----- caption --------------------------------------------------
  if (isTRUE(show_stats) && length(stats_entries) && is.null(caption)) {
    scale_str <- if (grepl("^log", y_scale)) sprintf("%s y-scale", y_scale)
                  else "linear y-scale"
    floor_str <- if (value == "tf")
                    sprintf("Floor=%g%%", .o("display_floor") * 100)
                  else NULL
    test_names <- unique(vapply(stats_entries,
                                  function(e) e$stats$overall$name,
                                  character(1)))
    cmp_desc <- switch(compare_by,
        "x"              = sprintf("comparing %s groups", x),
        "color"          = sprintf("comparing %s within each %s", color, x),
        "x_within_color" = sprintf("comparing %s within each %s", x, color))
    hdr_parts <- c(floor_str, scale_str,
                    sprintf("Test: %s (%s)",
                              paste(test_names, collapse = " / "), cmp_desc))
    hdr_parts <- hdr_parts[!vapply(hdr_parts, is.null, logical(1))]
    header <- paste(hdr_parts, collapse = "; ")
    if (length(stats_entries) == 1L) {
      e <- stats_entries[[1L]]
      n_str <- paste(sprintf("%s=%d",
                                names(e$stats$n_per_group),
                                e$stats$n_per_group),
                       collapse = ", ")
      caption <- sprintf("%s; overall p = %s; N: %s",
                          header,
                          .format_p_value(e$stats$overall$p), n_str)
    } else {
      caption <- sprintf("%s\n%d comparisons; see per-panel brackets",
                          header, length(stats_entries))
    }
  }

  # ----- defaults + finalize --------------------------------------
  if (is.null(title)) {
    title <- if (value == "tf")
                paste0("ctDNA boxplot: ",
                        if (length(value_cols) > 1L) "all TF methods"
                          else panel_labels[[value_cols]],
                        " by ", x,
                        if (!is.null(color)) paste0(" \u00d7 ", color) else "")
              else
                paste0("ctDNA boxplot: ", panel_labels[[value]],
                        " by ", x,
                        if (!is.null(color)) paste0(" \u00d7 ", color) else "")
  }
  if (is.null(ylab)) {
    ylab <- if (length(value_cols) == 1L)
                paste0(panel_labels[[value_cols]],
                        if (dominant_unit %in% c("percent", "raw_pct")) " (%)"
                        else "")
              else "value"
  }

  .finalize(p, title, subtitle, caption, xlab, ylab,
             legend_position = legend_position)
}
