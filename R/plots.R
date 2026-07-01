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


.ctdna_log_minor_breaks <- function(limits) {
  limits <- limits[is.finite(limits) & limits > 0]
  if (length(limits) < 1) return(numeric(0))
  rng <- log10(range(limits))
  p   <- seq(floor(rng[1]) - 1, ceiling(rng[2]) + 1)
  as.vector(outer(c(2,3,4,5,6,7,8,9), 10^p))
}

# v0.48.0: x-axis log scale layers with minor gridlines + bottom logticks
# (mirror of the y resolver, so a logarithmic x-axis also gets minor ticks).
.ctdna_x_log_layers <- function(x_scale = "linear") {
  x_scale <- tolower(trimws(as.character(x_scale)[1]))
  if (!x_scale %in% c("log10","log","log2"))
    return(list(ggplot2::scale_x_continuous()))
  trans <- switch(x_scale, log2 = "log2", log10 = "log10",
                  log = scales::log_trans(base = exp(1)))
  sc <- ggplot2::scale_x_continuous(trans = trans, minor_breaks = .ctdna_log_minor_breaks)
  if (x_scale %in% c("log10","log"))
    return(list(sc, ggplot2::annotation_logticks(
      sides = "b", outside = FALSE, colour = "grey45",
      short = ggplot2::unit(0.5, "mm"), mid = ggplot2::unit(1.2, "mm"),
      long = ggplot2::unit(2, "mm"))))
  list(sc)
}

.ctdna_resolve_y_scale <- function(y_scale = "linear", labels = NULL,
                                     floor_value = NULL, ...) {
  if (is.null(y_scale)) y_scale <- "linear"
  y_scale <- tolower(trimws(as.character(y_scale)[1]))
  args <- list(...)
  if (!is.null(labels)) args$labels <- labels
  if (is.null(args$expand)) args$expand <- ggplot2::expansion(mult = c(0.05, 0.03))

  if (y_scale == "linear")
    return(list(do.call(ggplot2::scale_y_continuous, args)))
  if (y_scale == "log1p") {
    args$trans <- scales::trans_new("log1p", function(x) log1p(x), function(x) expm1(x))
    return(list(do.call(ggplot2::scale_y_continuous, args)))
  }
  if (y_scale == "symlog") {
    args$trans <- scales::trans_new("symlog",
      function(x) sign(x) * log1p(abs(x)), function(x) sign(x) * expm1(abs(x)))
    return(list(do.call(ggplot2::scale_y_continuous, args)))
  }
  trans <- switch(y_scale, "log2" = "log2", "log10" = "log10",
                   "log" = scales::log_trans(base = exp(1)), NULL)
  if (is.null(trans)) {
    m <- regmatches(y_scale, regexec("^log[_\\.]?([0-9]+(?:\\.[0-9]+)?)$", y_scale, perl = TRUE))[[1]]
    if (length(m) == 2) {
      base <- as.numeric(m[2])
      if (!is.na(base) && base > 0 && base != 1) trans <- scales::log_trans(base = base)
    }
  }
  if (is.null(trans)) {
    warning(sprintf("y_scale='%s' is not a recognized transform. Falling back to linear.", y_scale),
            call. = FALSE)
    return(list(do.call(ggplot2::scale_y_continuous, args)))
  }
  # v0.44.3: minor gridlines + log tick marks between major ticks
  if (is.null(args$minor_breaks)) args$minor_breaks <- .ctdna_log_minor_breaks
  args$trans <- trans
  sc <- do.call(ggplot2::scale_y_continuous, args)
  if (y_scale %in% c("log10", "log"))
    return(list(sc, ggplot2::annotation_logticks(
      sides = "l", outside = FALSE, colour = "grey45",
      short = ggplot2::unit(0.5, "mm"), mid = ggplot2::unit(1.2, "mm"),
      long = ggplot2::unit(2, "mm"))))
  list(sc)
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
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 30, seed = 1)
#' prep <- ctdna_prepare(sim, verbose = FALSE)
#'
#' # Default — best-ratio (deepest drop from baseline), grouped by dose
#' ctdna_plot_reduction(prep$samples, scheme = "three", metric = "ratio",
#'                       at = "best", mr_threshold = 0.5)
#'
#' # Percent-change at a specific visit (tolerant matching: "C2D1" works
#' # but so does "Cycle 2 Day 1")
#' ctdna_plot_reduction(prep$samples, scheme = "three",
#'                       metric = "pct_change", at = "Cycle 2 Day 1")
#'
#' # Raw TF at the best on-treatment visit
#' ctdna_plot_reduction(prep$samples, scheme = "three",
#'                       metric = "tf", at = "best")
#' @export
ctdna_plot_reduction <- function(prep,
                           group_by    = .o("dose"),
                           subgroup_by = .o("recist"),
                           compare_by  = c("group","subgroup","subgroup_within_group","group_within_subgroup"),
                           pairwise    = TRUE,
                           p_adjust    = FALSE,
                           scheme = c("three","two","four","four_alt","five","six"),
                           at     = "all",
                           metric = c("ratio","tf","pct_change"),
                           measure = "methyl",
                           value_col = NULL,
                           mr_threshold = 0.5,
                           stat = c("wilcox","t","none"),
                           show_stats = .o("show_stats"),
                           show_n = .o("show_n"),
                           drop_na = .o("drop_na"),
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
    xlim = NULL) {
  compare_by <- match.arg(compare_by)
  # v0.44.0: accept a ctdna_prep (uses $samples) or a per-sample frame.
  df <- .ctdna_grain(prep, "samples")
  if ("Record_type" %in% names(df)) df <- df[!(as.character(df$Record_type) %in% "tumor"), , drop = FALSE]
  group_by <- .resolve_col(group_by, df, what = "group_by")
  if (!is.null(subgroup_by)) {
    subgroup_by <- tryCatch(.resolve_col(subgroup_by, df, what = "subgroup_by"),
                            error = function(e) NULL)
  }
  if (compare_by != "group" && is.null(subgroup_by))
    stop("ctdna_plot_reduction: compare_by = '", compare_by,
         "' requires subgroup_by.", call. = FALSE)

  # v0.42.5 BUG-FIX: match.arg(metric) MUST run before the multi-method
  # recursion. Inside the lapply closure, match.arg() looks up the
  # formals of the calling function, which is the anonymous function
  # (not ctdna_plot_reduction), and errors with "arg must be of
  # length 1" when metric is the default vector.
  metric <- match.arg(metric)

  # v0.25: measure accepts a vector for side-by-side comparison.
  if (length(measure) > 1L || identical(measure, "all")) {
    methods_resolved <- .resolve_methods(measure, domain = "tf")
    plots <- lapply(methods_resolved, function(mm) {
      ctdna_plot_reduction(df, group_by = group_by, subgroup_by = subgroup_by,
        compare_by = compare_by, pairwise = pairwise, p_adjust = p_adjust, scheme = scheme, at = at, metric = metric,
        measure = mm,
        value_col = if (mm == "custom") value_col else NULL,
        mr_threshold = mr_threshold, stat = stat,
        show_stats = show_stats, show_n = show_n, drop_na = drop_na,
        title = sprintf("%s — %s", title, .resolve_method(mm, "tf")$label),
        subtitle = subtitle, caption = caption, xlab = xlab, ylab = ylab,
        point_size = point_size, point_alpha = point_alpha,
        stat_position = stat_position, legend_position = legend_position,
        y_scale = y_scale, ylim = ylim)
    })
    return(patchwork::wrap_plots(plots, ncol = 1L))
  }

  .stats_label <- NULL

  scheme <- match.arg(scheme); stat <- match.arg(stat)
  m <- .resolve_method(measure, "tf", value_col = value_col)
  tf_col <- m$value_col

  # v0.42.5 BUG-FIX: ctdna_metric_at computes pct_change as
  #   100 * (tf_at - baseline_tf) / baseline_tf.
  # When baseline_tf is near LOQ (e.g. 0.001), even a moderate tf_at
  # blows up the metric (0.5 / 0.001 = +49900%). Floor baseline at LOQ
  # before dividing, consistent with how the rest of the library
  # treats sub-LOQ values (display_floor). Only relevant for pct_change;
  # ratio / tf metrics use log scale + ctdna_floor_tf on the output.
  s_col <- .o("subject"); t_col <- .o("time"); bl <- .o("baseline")
  floor_v <- .o("display_floor")

  # one-`at` metric for a subject set (handles pct_change baseline floor)
  compute_one <- function(at_one) {
    if (metric == "pct_change") {
      base <- df[df[[t_col]] == bl, c(s_col, tf_col), drop = FALSE]
      names(base) <- c(s_col, "baseline_tf_raw")
      base$baseline_tf <- pmax(base$baseline_tf_raw, floor_v, na.rm = TRUE)
      ma_raw <- ctdna_metric_at(df, metric = "tf", at = at_one, tf_col = tf_col)
      ma_raw <- merge(ma_raw[, c(s_col, "value", "time_used")],
                       base[, c(s_col, "baseline_tf")], by = s_col, all.x = TRUE)
      names(ma_raw)[names(ma_raw) == "value"] <- "tf_at"
      data.frame(stats::setNames(list(ma_raw[[s_col]]), s_col),
                 baseline_tf = ma_raw$baseline_tf,
                 value = 100 * (ma_raw$tf_at - ma_raw$baseline_tf) / ma_raw$baseline_tf,
                 time_used = ma_raw$time_used, check.names = FALSE,
                 stringsAsFactors = FALSE)
    } else {
      ctdna_metric_at(df, metric = metric, at = at_one, tf_col = tf_col)
    }
  }

  # v0.44.3: `at` controls the visit facets.
  #   "best"            -> single "Best" panel
  #   "all" (default)   -> one panel per on-treatment visit + a "Best" panel
  #   single visit/vec  -> a panel per requested visit (no Best)
  onx_visits <- setdiff(unique(as.character(df[[t_col]])), bl)
  onx_visits <- as.character(.ctdna_factor_canonical(onx_visits, kind = "visit"))
  onx_visits <- onx_visits[!is.na(onx_visits)]
  if (identical(at, "all")) {
    facet_ats <- c(stats::setNames(onx_visits, onx_visits), Best = "best")
  } else if (identical(at, "best") || identical(at, "BEST")) {
    facet_ats <- c(Best = "best")
  } else {
    av <- as.character(at); facet_ats <- stats::setNames(av, av)
  }

  ma <- do.call(rbind, lapply(seq_along(facet_ats), function(i) {
    one <- compute_one(unname(facet_ats[i]))
    one$.vfacet <- names(facet_ats)[i]
    one
  }))
  vlev <- names(facet_ats)              # facet order (visits then Best)
  ma$.vfacet <- factor(ma$.vfacet, levels = vlev)

  # Join per-subject grouping columns (group_by + optional subgroup_by)
  meta_cols <- intersect(unique(c(group_by, subgroup_by, .o("recist"), .o("dose"))), names(df))
  meta <- unique(df[, c(.o("subject"), meta_cols), drop = FALSE])
  d <- merge(ma, meta, by = .o("subject"), all.x = TRUE)

  recist_col <- .o("recist")
  # stratify whichever grouping column is RECIST
  if (identical(group_by, recist_col))    d[[group_by]]    <- ctdna_stratify_recist(d[[group_by]], scheme)
  if (!is.null(subgroup_by) && identical(subgroup_by, recist_col))
    d[[subgroup_by]] <- ctdna_stratify_recist(d[[subgroup_by]], scheme)
  d[[group_by]] <- .ctdna_factor_canonical(d[[group_by]], kind = group_by)
  if (!is.null(subgroup_by))
    d[[subgroup_by]] <- .ctdna_factor_canonical(d[[subgroup_by]], kind = subgroup_by)

  # v0.44.0: drop rows with NA grouping
  if (isTRUE(drop_na)) {
    keep <- !is.na(d[[group_by]])
    if (!is.null(subgroup_by)) keep <- keep & !is.na(d[[subgroup_by]])
    d <- d[keep, , drop = FALSE]
  }

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

  # ylab (faceting carries the visit, so it's not embedded here)
  if (is.null(ylab)) ylab <- switch(metric,
    ratio      = sprintf("%s / %s(baseline) [log]", m$label, m$label),
    tf         = sprintf("%s [log]", m$label),
    pct_change = sprintf("Percent change in %s", m$label))

  # caption: line 1 = settings, line 2 = the exact formula used (#5)
  if (is.null(caption)) {
    line1 <- switch(metric,
      ratio      = sprintf("metric = ratio; floor=%g%%, log y; MR threshold = %g (ratio<=thr -> Response)",
                            .o("display_floor"), mr_threshold),
      tf         = sprintf("metric = TF; floor=%g%%, log y", .o("display_floor")),
      pct_change = "metric = pct_change; negative = reduction")
    line2 <- switch(metric,
      ratio      = sprintf("ratio = %s(visit) / %s(baseline)", m$label, m$label),
      tf         = sprintf("value = %s at visit (raw)", m$label),
      pct_change = sprintf("%% change = 100 \u00d7 (%s(visit) \u2212 %s(baseline)) / %s(baseline)",
                            m$label, m$label, m$label))
    caption <- paste0(line1, "\n", line2)
  }

  # MR threshold line is meaningful only for ratio (handled below)
  stat_col <- switch(compare_by, group = group_by, subgroup = subgroup_by,
                     subgroup_within_group = subgroup_by,
                     group_within_subgroup = group_by)
  # per-facet stats (gate / overall-p source)
  facet_stats <- NULL
  if (isTRUE(show_stats) && !is.null(stat_col)) {
    facet_stats <- lapply(split(d_capped, d_capped$.vfacet), function(s)
      if (nrow(s)) .compute_box_stats(s, group_col = stat_col, value_col = ".y", p_adjust = p_adjust) else NULL)
    facet_stats <- facet_stats[!vapply(facet_stats, is.null, logical(1))]
  }

  # build the plot from d_capped so clipped outliers don't escape the panel
  if (is.null(subgroup_by)) {
    p <- ggplot2::ggplot(d_capped, ggplot2::aes(.data[[group_by]], .data$.y,
                                         color = .data[[group_by]], fill = .data[[group_by]])) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.3, linewidth = 0.7,
                            fatten = 2.2, width = .o("box_width")) +
      ggplot2::geom_jitter(width = 0.15, size = point_size, alpha = point_alpha,
                           shape = 21, stroke = 0.3, color = "grey20")
    color_fn <- if (identical(group_by, recist_col)) ctdna_scale_recist else ctdna_scale_dose
    legend_lab <- group_by
  } else {
    p <- ggplot2::ggplot(d_capped, ggplot2::aes(.data[[group_by]], .data$.y,
                                         color = .data[[subgroup_by]], fill = .data[[subgroup_by]])) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.3, linewidth = 0.7,
                            fatten = 2.2, width = .o("box_width"),
                            position = ggplot2::position_dodge(0.8, preserve = "single")) +
      ggplot2::geom_point(position = ggplot2::position_jitterdodge(
        jitter.width = 0.15, dodge.width = 0.8),
        size = point_size, alpha = point_alpha,
        shape = 21, stroke = 0.3, color = "grey20")
    color_fn <- if (identical(subgroup_by, recist_col)) ctdna_scale_recist else ctdna_scale_dose
    legend_lab <- subgroup_by
  }
  p <- p +
    color_fn(aesthetics = c("color", "fill")) +
    ggplot2::facet_wrap(~ .vfacet, scales = "free_x") +
    ggplot2::labs(x = group_by, color = legend_lab, fill = legend_lab) +
    ctdna_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

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
    # inside-panel n per (facet x group), so each visit panel shows its own n
    counts <- as.data.frame(table(d_capped$.vfacet, d_capped[[group_by]]))
    names(counts) <- c(".vfacet", "g", "n")
    counts <- counts[counts$n > 0, , drop = FALSE]
    counts$g <- factor(counts$g, levels = levels(d_capped[[group_by]]))
    counts$.vfacet <- factor(counts$.vfacet, levels = levels(d_capped$.vfacet))
    if (nrow(counts) > 0L) {
      p <- p + ggplot2::geom_label(
        data = counts,
        mapping = ggplot2::aes(x = .data$g, y = -Inf, label = paste0("n=", .data$n)),
        inherit.aes = FALSE, vjust = -0.2, size = 2.8,
        label.padding = ggplot2::unit(0.12, "lines"),
        label.size = 0, color = "grey15", fill = "white", fontface = "bold")
    }
  }
  if (isTRUE(show_stats) && length(facet_stats)) {
    # x-aligned brackets per facet, only when comparison IS the x grouping
    if (compare_by == "group") {
      for (v in names(facet_stats)) {
        sub <- d_capped[d_capped$.vfacet == v, , drop = FALSE]
        p <- .add_brackets(p, facet_stats[[v]], sub,
                            group_col = group_by, value_col = ".y",
                            y_scale = if (use_log) "log10" else "linear",
                            facet_data = list(.vfacet = v))
      }
      pv <- paste(vapply(names(facet_stats), function(v)
        sprintf("%s: p=%s", v, .format_p_value(facet_stats[[v]]$overall$p)),
        character(1)), collapse = " | ")
      caption <- paste0(caption, "\nCompare by ", group_by, " — ", pv)

    } else if (compare_by == "subgroup_within_group" && !is.null(subgroup_by)) {
      # within each x-group, compare subgroups (pairwise brackets by default)
      for (v in levels(d_capped$.vfacet)) {
        sub <- d_capped[d_capped$.vfacet == v, , drop = FALSE]
        p <- .add_dodge_brackets(p, sub, group_col = group_by, subgroup_col = subgroup_by,
                                 value_col = ".y", mode = "within_group", pairwise = pairwise, p_adjust = p_adjust,
                                 facet_data = list(.vfacet = v),
                                 y_scale = if (use_log) "log10" else "linear")
      }
      caption <- paste0(caption, "\nCompare ", subgroup_by, " within each ", group_by,
                        if (pairwise) " (pairwise)" else " (overall)")

    } else if (compare_by == "group_within_subgroup" && !is.null(subgroup_by)) {
      # within each subgroup, compare x-groups (pairwise brackets by default)
      for (v in levels(d_capped$.vfacet)) {
        sub <- d_capped[d_capped$.vfacet == v, , drop = FALSE]
        p <- .add_dodge_brackets(p, sub, group_col = group_by, subgroup_col = subgroup_by,
                                 value_col = ".y", mode = "within_subgroup", pairwise = pairwise, p_adjust = p_adjust,
                                 facet_data = list(.vfacet = v),
                                 y_scale = if (use_log) "log10" else "linear")
      }
      caption <- paste0(caption, "\nCompare ", group_by, " within each ", subgroup_by,
                        if (pairwise) " (pairwise)" else " (overall)")

    } else {  # compare_by == "subgroup": one pooled p per facet, centered at top
      ng <- nlevels(droplevels(d_capped[[group_by]]))
      labdf <- do.call(rbind, lapply(names(facet_stats), function(v)
        data.frame(.vfacet = factor(v, levels = levels(d_capped$.vfacet)), .cx = (ng + 1) / 2,
                   lab = sprintf("%s p=%s", subgroup_by, .format_p_value(facet_stats[[v]]$overall$p)),
                   stringsAsFactors = FALSE)))
      p <- p + ggplot2::geom_text(data = labdf,
        mapping = ggplot2::aes(x = .data$.cx, y = Inf, label = .data$lab),
        inherit.aes = FALSE, vjust = 1.5, size = 3.0, fontface = "bold", color = "grey15")
      caption <- paste0(caption, "\nCompare by ", subgroup_by, " (pooled)")
    }
    if (.resolve_p_adjust(p_adjust) != "none")
      caption <- paste0(caption, "; pairwise p ", .resolve_p_adjust(p_adjust), "-adjusted")
  }
  # v0.42.5: coord_cartesian(clip="off") always; ylim if provided.
  # v0.44.6: for dodge-bracket modes the data is already capped to ylim, so let
  # the axis autoscale (the stacked brackets define the top) instead of clipping.
  cc_args <- list(clip = "off")
  if (!is.null(ylim) && !(compare_by %in% c("subgroup_within_group","group_within_subgroup")))
    cc_args$ylim <- ylim
  if (!is.null(xlim)) cc_args$xlim <- xlim
  p <- p + do.call(ggplot2::coord_cartesian, cc_args)

  .finalize(p, title, subtitle, caption, xlab, ylab,
             legend_position = legend_position)
}


# ---- D4: Categorical concordance (generalized; replaces mr_recist) ----------

#' Concordance between two categorical/ordinal variables
#'
#' Renders the contingency tile of two categorical (or ordinal) variables with
#' agreement statistics annotated. Defaults reproduce the former
#' \code{ctdna_plot_mr_recist}: \code{x = "RECIST"} (stratified by \code{scheme})
#' versus \code{y = "MR"} (molecular response, derived by thresholding the
#' per-subject best ctDNA reduction at \code{mr_threshold}). Either axis can be
#' any categorical column on the per-sample frame (e.g. \code{"Dose"},
#' \code{"Arm"}, \code{"Cancertype"}, or a second RECIST scheme).
#'
#' @section Statistics (what each number means):
#' \itemize{
#'   \item \strong{\% agreement} — proportion of subjects falling on the
#'     matching diagonal (only meaningful when both axes share categories).
#'   \item \strong{Cohen's kappa} — agreement corrected for chance:
#'     0 = chance, 1 = perfect, <0 = worse than chance.
#'   \item \strong{Weighted kappa} (ordinal) — like kappa but gives partial
#'     credit when the two raters land in adjacent categories.
#'   \item \strong{McNemar} (2x2) — tests whether the two raters' marginal
#'     rates differ (a systematic shift), not overall agreement.
#'   \item \strong{Chi-square} — tests association/independence between the
#'     two variables (used when categories differ, e.g. MR vs RECIST).
#'   \item \strong{Cramer's V} — effect size for that association, 0 (none)
#'     to 1 (perfect).
#'   \item \strong{N} — subjects contributing to the table.
#' }
#' \code{stat = "auto"} shows agreement + kappa when the axes share a category
#' set, otherwise chi-square + Cramer's V. Continuous variables are rejected.
#'
#' @param prep A \code{ctdna_prep} (uses \code{$samples}) or a per-sample frame.
#' @param x,y Axis variables. Special tokens \code{"RECIST"} and \code{"MR"} are
#'   derived; any other value is resolved as a categorical column. Defaults
#'   \code{x = "RECIST"}, \code{y = "MR"}.
#' @param scheme RECIST stratification (applied to a \code{"RECIST"} axis).
#' @param mr_timepoint For an \code{"MR"} axis: \code{"best"} (default), a visit
#'   label, or \code{"use_column"} to use a pre-computed MR column.
#' @param mr_threshold Reduction ratio at/below which a subject is a molecular
#'   responder (default 0.5).
#' @param ordinal Treat the axes as ordinal and report weighted kappa
#'   (default FALSE).
#' @param weights Weighting for ordinal weighted kappa: \code{"quadratic"} or \code{"linear"}.
#' @param stat \code{"auto"} (default), \code{"kappa"}, \code{"agreement"},
#'   \code{"mcnemar"}, \code{"association"}, or \code{"none"}.
#' @param show_stats Annotate the statistics in the caption (default TRUE).
#' @param title,subtitle,caption,xlab,ylab,legend_position Display.
#' @return A ggplot tile (categorical concordance). The contingency table and
#'   statistics are attached as \code{attr(p, "concordance")}.
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 40, seed = 1)
#' prep <- ctdna_prepare(infinity_report = sim$infinity_report,
#'                       tf_change = sim$tf_change, clinical = sim$clinical, verbose = FALSE)
#' ctdna_concordance(prep)                                  # MR vs RECIST (default)
#' ctdna_concordance(prep, x = "Dose", y = "RECIST")        # any two categoricals
#' @export
ctdna_concordance <- function(prep,
                              x = "RECIST",
                              y = "MR",
                              scheme = c("three","two","four","four_alt","five","six"),
                              mr_timepoint = "best",
                              mr_threshold = 0.5,
                              ordinal = FALSE,
                              weights = c("quadratic","linear"),
                              stat = c("auto","kappa","agreement","mcnemar","association","none"),
                              show_stats = .o("show_stats"),
                              title = NULL,
                              subtitle = NULL,
                              caption = NULL,
                              xlab = NULL,
                              ylab = NULL,
                              legend_position = .o("legend_position")) {
  scheme  <- match.arg(scheme)
  weights <- match.arg(weights)
  stat    <- match.arg(stat)
  df    <- .ctdna_grain(prep, "samples")
  if ("Record_type" %in% names(df)) df <- df[!(as.character(df$Record_type) %in% "tumor"), , drop = FALSE]
  subj  <- .o("subject")
  src   <- NULL

  # resolve one axis token to a per-subject factor (+ label, ordinal flag)
  axis_vec <- function(token) {
    # v0.57.0: categorical modality tokens (mutation:GENE / alteration:GENE) ->
    # per-subject factor from prep$variants (subsumes ctdna_plot_vs_mutation x RECIST).
    .mr <- .ctdna_modality_resolve(prep, token, scheme)
    if (!is.null(.mr) && identical(.mr$kind, "cat")) {
      src <<- sprintf("%s from $variants", .mr$label)
      return(list(v = factor(.mr$values), label = .mr$label, ordinal = FALSE))
    }
    if (toupper(token) == "MR") {
      if (identical(mr_timepoint, "use_column")) {
        if (!.o("mr") %in% names(df))
          stop("y/x='MR' with mr_timepoint='use_column' but column '", .o("mr"), "' is absent.", call. = FALSE)
        raw <- tapply(df[[.o("mr")]], df[[subj]], function(z) z[!is.na(z)][1])
        v   <- factor(raw); src <<- sprintf("MR from df$%s", .o("mr"))
      } else {
        mr  <- ctdna_compute_mr(df, at = mr_timepoint, threshold = mr_threshold)
        v   <- stats::setNames(factor(mr$MR), mr[[subj]])
        src <<- sprintf("MR recomputed at %s; ratio threshold = %g",
                        if (identical(mr_timepoint, "best")) "best on-treatment visit"
                        else paste0("visit ", mr_timepoint), mr_threshold)
      }
      list(v = v, label = "Molecular Response", ordinal = FALSE)
    } else if (toupper(token) == "RECIST") {
      raw <- tapply(df[[.o("recist")]], df[[subj]], function(z) z[!is.na(z)][1])
      list(v = stats::setNames(ctdna_stratify_recist(raw, scheme), names(raw)),
           label = "RECIST", ordinal = TRUE)
    } else {
      col <- .resolve_col(token, df, what = "axis")
      raw <- tapply(df[[col]], df[[subj]], function(z) z[!is.na(z)][1])
      if (!.ctdna_is_categorical(raw))
        stop(sprintf(paste0("Axis '%s' looks continuous (%d distinct values). ",
             "ctdna_concordance compares only categorical/ordinal variables. ",
             "For continuous-vs-categorical use ctdna_boxplot(); for ",
             "continuous-vs-continuous use ctdna_plot_scatter()."),
             col, length(unique(raw[is.finite(raw)]))), call. = FALSE)
      list(v = stats::setNames(.ctdna_factor_canonical(raw, kind = col), names(raw)),
           label = col, ordinal = is.ordered(raw))
    }
  }

  xa <- axis_vec(x); ya <- axis_vec(y)
  if (identical(xa$label, ya$label)) {        # same source on both axes
    xa$label <- paste0(xa$label, " (x)"); ya$label <- paste0(ya$label, " (y)")
  }
  keep <- intersect(names(xa$v), names(ya$v))
  if (!length(keep)) stop("No subjects shared between the two axes.", call. = FALSE)
  xf <- droplevels(factor(xa$v[keep])); yf <- droplevels(factor(ya$v[keep]))
  tab <- table(setNames(list(yf, xf), c(ya$label, xa$label)))

  is_ord     <- isTRUE(ordinal) || (xa$ordinal && ya$ordinal)
  same_levs  <- setequal(rownames(tab), colnames(tab))
  ct  <- ctdna_concordance_test(tab, method = "both")
  chisq_p <- tryCatch(suppressWarnings(stats::chisq.test(tab, correct = FALSE)$p.value),
                      error = function(e) NA_real_)
  cv  <- .cramers_v(tab)
  wk  <- if (is_ord && same_levs) ctdna_weighted_kappa(tab, weights) else NA_real_

  if (is.null(src)) src <- sprintf("%s vs %s", ya$label, xa$label)
  if (is.null(caption)) caption <- src
  if (is.null(title))   title <- sprintf("%s vs %s concordance", ya$label, xa$label)
  if (is.null(xlab))    xlab <- xa$label
  if (is.null(ylab))    ylab <- ya$label

  if (isTRUE(show_stats) && stat != "none") {
    show_agree <- stat %in% c("agreement","both") ||
                  (stat == "auto" && same_levs)
    show_kappa <- stat %in% c("kappa","both") ||
                  (stat == "auto" && same_levs)
    show_assoc <- stat == "association" || (stat == "auto" && !same_levs)
    parts <- c()
    if (show_agree) parts <- c(parts, sprintf("Agreement = %.0f%%", 100 * ct$agreement))
    if (show_kappa) parts <- c(parts, sprintf("kappa = %.2f", ct$kappa))
    if (show_kappa && is_ord && !is.na(wk))
      parts <- c(parts, sprintf("weighted kappa (%s) = %.2f", weights, wk))
    if (stat == "mcnemar" && !is.na(ct$mcnemar_p))
      parts <- c(parts, sprintf("McNemar %s", ctdna_pval_label(ct$mcnemar_p)))
    if (show_assoc) {
      if (!is.na(chisq_p)) parts <- c(parts, sprintf("Chi-square %s", ctdna_pval_label(chisq_p)))
      if (!is.na(cv))      parts <- c(parts, sprintf("Cramer's V = %.2f", cv))
    }
    parts <- c(parts, sprintf("N = %d", ct$n))
    caption <- paste(caption, paste(parts, collapse = "  |  "), sep = "  |  ")
    # plain-English gloss for whatever was reported
    gl <- c()
    if (show_agree) gl <- c(gl, "Agreement = % of subjects on the matching diagonal")
    if (show_kappa) gl <- c(gl, "kappa = chance-corrected agreement (0 = chance, 1 = perfect)")
    if (show_kappa && is_ord && !is.na(wk))
      gl <- c(gl, "weighted kappa = kappa that credits near-miss (adjacent) categories")
    if (stat == "mcnemar" && !is.na(ct$mcnemar_p))
      gl <- c(gl, "McNemar = test of a marginal shift between the two raters (2x2)")
    if (show_assoc)
      gl <- c(gl, "Chi-square = test of association (independence); Cramer's V = association strength, 0 (none) to 1 (perfect)")
    if (length(gl)) caption <- paste0(caption, "\n", paste(gl, collapse = "; "))
  }

  long <- as.data.frame(tab, stringsAsFactors = TRUE)
  names(long)[1:2] <- c(ya$label, xa$label)   # dim order is (y, x); keep exact labels
  p <- ggplot2::ggplot(long, ggplot2::aes(.data[[xa$label]], .data[[ya$label]], fill = .data$Freq)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = .data$Freq)) +
    ggplot2::scale_fill_gradient(low = "white", high = "#3182BD") +
    ctdna_theme()
  p <- .finalize(p, title, subtitle, caption, xlab, ylab, legend_position = legend_position)
  attr(p, "concordance") <- list(
    table = tab, agreement = ct$agreement, kappa = ct$kappa,
    weighted_kappa = wk, mcnemar_p = ct$mcnemar_p, chisq_p = chisq_p,
    cramers_v = cv, n = ct$n, source = src)
  p
}


# ---- D5: Longitudinal ctDNA trajectories (merged) ---------------------------

# Whole-trajectory comparison between two panels: LMM logTF ~ panel*visit +
# (1|subject) when a subject id and >=2 repeated measures exist; AUC-Wilcoxon
# fallback otherwise. Returns list(level_p, shape_p, method).
.longitudinal_traj_test <- function(sub, t_col, s_col) {
  sub <- sub[is.finite(sub$.y) & !is.na(sub[[t_col]]), , drop = FALSE]
  if (!nrow(sub)) return(list(level_p = NA_real_, shape_p = NA_real_, method = "none"))
  sub$.logy  <- log10(pmax(sub$.y, .Machine$double.eps))
  sub$.panel <- droplevels(factor(sub$.panel))
  if (nlevels(sub$.panel) < 2L)
    return(list(level_p = NA_real_, shape_p = NA_real_, method = "none"))
  has_subj <- s_col %in% names(sub) && anyDuplicated(sub[[s_col]]) > 0L
  if (has_subj && requireNamespace("nlme", quietly = TRUE)) {
    sub$.subj <- factor(sub[[s_col]]); sub$.vf <- factor(sub[[t_col]])
    fit <- tryCatch({
      m  <- nlme::lme(.logy ~ .panel * .vf, random = ~ 1 | .subj,
                      data = sub, na.action = stats::na.omit, method = "REML")
      an <- as.data.frame(stats::anova(m)); rn <- rownames(an)
      list(level_p = an[rn == ".panel", "p-value"][1],
           shape_p = an[grepl(":", rn), "p-value"][1],
           method  = "LMM logTF ~ panel*visit + (1|subject)")
    }, error = function(e) NULL)
    if (!is.null(fit)) return(fit)
  }
  # AUC fallback: per-subject area under log-TF over visit index, compared by panel
  idx <- as.integer(factor(sub[[t_col]], levels = levels(factor(sub[[t_col]]))))
  sub$.vi <- idx
  key <- if (has_subj) sub[[s_col]] else seq_len(nrow(sub))
  auc <- tapply(seq_len(nrow(sub)), list(key), function(i) {
    o <- order(sub$.vi[i]); x <- sub$.vi[i][o]; y <- sub$.logy[i][o]
    if (length(x) < 2L) return(mean(y))
    sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1)) / 2) / (max(x) - min(x))
  })
  pan <- tapply(seq_len(nrow(sub)), list(key), function(i) as.character(sub$.panel[i][1]))
  d2  <- data.frame(auc = as.numeric(auc), panel = factor(unlist(pan)))
  st  <- .compute_box_stats(d2, group_col = "panel", value_col = "auc")
  list(level_p = if (is.null(st)) NA_real_ else st$overall$p,
       shape_p = NA_real_, method = "per-subject AUC (Wilcoxon)")
}

#' Longitudinal ctDNA reduction trajectories
#'
#' Per-subject reduction-ratio trajectories with the x-axis fixed to Visit.
#' Merges the former \code{ctdna_plot_longitudinal} and
#' \code{ctdna_plot_longitudinal_recist}: \code{group_by} colours the lines
#' (e.g. RECIST, Dose) and \code{facet} lays out panels (one variable ->
#' \code{facet_wrap}, two -> \code{facet_grid}). The thick line is the median.
#'
#' @section Statistics:
#' On the plot, when \code{group_by} is set, a per-visit omnibus p (Kruskal-Wallis
#' / Wilcoxon) compares the groups at each visit. Two tables ride along on the
#' object and are read back with \code{longitudinal_stats(p)}:
#' \code{$within_panel} (overall + pairwise group comparison per panel, pooled
#' across visits) and \code{$between_panels} (whole-trajectory LMM comparison of
#' every panel pair). The median line is never used in any test.
#'
#' @param prep A \code{ctdna_prep} (uses \code{$samples}) or a per-sample frame.
#' @param group_by Line/point colour stratifier (Dose/RECIST/Arm/Indication/...);
#'   \code{NULL} for a single colour.
#' @param facet One variable (\code{facet_wrap}) or a length-2 vector
#'   (\code{facet_grid}). x-axis is always Visit.
#' @param measure ctDNA quantity: \code{"methyl"}/\code{"maxvaf"}/\code{"meanvaf"}/\code{"custom"}, or \code{"all"}.
#' @param value_col Column for \code{measure = "custom"}.
#' @param scheme RECIST stratification scheme.
#' @param stat Logical; compute and show the comparisons (default TRUE).
#' @param p_adjust Multiple-test correction for pairwise families. \code{FALSE}/
#'   \code{"none"} = none (default); \code{TRUE} = BH; or a method name
#'   (case-insensitive: bh/fdr, by, holm, hochberg, hommel, bonferroni).
#' @param show_stats,show_n Display toggles.
#' @param y_scale,ylim,scales Y-axis transform (log gets minor ticks), limits, facet scales.
#' @param title,subtitle,caption,xlab,ylab,point_size,point_alpha,legend_position Display.
#' @return A ggplot (classed \code{ctdna_longitudinal}) carrying comparison
#'   tables; retrieve them with \code{longitudinal_stats()}.
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 40, seed = 1)
#' prep <- ctdna_prepare(infinity_report = sim$infinity_report,
#'                       tf_change = sim$tf_change, clinical = sim$clinical, verbose = FALSE)
#' p <- ctdna_plot_longitudinal(prep, group_by = "RECIST", facet = "Dose")
#' longitudinal_stats(p)$within_panel
#' @export
ctdna_plot_longitudinal <- function(prep,
                              group_by = NULL,
                              facet    = .o("dose"),
                              measure  = "methyl",
                              value_col = NULL,
                              scheme   = c("three","two","four","four_alt","five","six"),
                              stat     = TRUE,
                              p_adjust = FALSE,
                              show_stats = .o("show_stats"),
                              show_n   = .o("show_n"),
                              y_scale  = "log10",
                              ylim     = NULL,
                              xlim     = NULL,
                              scales   = c("fixed","free_y","free_x","free"),
                              title = "Longitudinal ctDNA reduction",
                              subtitle = "Thick line = median",
                              caption = NULL,
                              xlab = "Visit",
                              ylab = NULL,
                              point_size = .o("point_size"),
                              point_alpha = .o("point_alpha"),
                              legend_position = .o("legend_position")) {
  scheme <- match.arg(scheme); scales <- match.arg(scales)
  df <- .ctdna_grain(prep, "samples")
  if ("Record_type" %in% names(df)) df <- df[!(as.character(df$Record_type) %in% "tumor"), , drop = FALSE]

  # multi-measure -> stacked panels
  if (length(measure) > 1L || identical(measure, "all")) {
    mm <- .resolve_methods(measure, domain = "tf")
    plots <- lapply(mm, function(one)
      ctdna_plot_longitudinal(df, group_by = group_by, facet = facet, measure = one,
        value_col = if (one == "custom") value_col else NULL, scheme = scheme,
        stat = stat, p_adjust = p_adjust, show_stats = show_stats, show_n = show_n,
        y_scale = y_scale, ylim = ylim, scales = scales,
        title = sprintf("%s - %s", title, .resolve_method(one, "tf")$label),
        subtitle = subtitle, caption = caption, xlab = xlab, ylab = ylab,
        point_size = point_size, point_alpha = point_alpha, legend_position = legend_position))
    return(patchwork::wrap_plots(plots, ncol = 1L))
  }

  m <- .resolve_method(measure, "tf", value_col = value_col)
  t_col <- .o("time"); s_col <- .o("subject")
  d <- ctdna_ratio(df, tf_col = m$value_col)
  d$.y <- ctdna_floor_tf(d$ratio)
  d[[t_col]] <- .ctdna_factor_canonical(d[[t_col]], kind = "time_point")
  if (is.null(ylab)) ylab <- sprintf("%s reduction ratio (log)", m$label)

  # resolve + canonicalise grouping / facet columns
  if (!is.null(group_by)) group_by <- .resolve_col(group_by, d, what = "group_by")
  facet <- if (is.null(facet) || !length(facet)) character(0)
           else vapply(as.character(facet), function(f)
                  tryCatch(.resolve_col(f, d, what = "facet"), error = function(e) NA_character_),
                  character(1))
  facet <- unique(facet[!is.na(facet) & facet %in% names(d)])
  if (length(facet) > 2L) facet <- facet[1:2]
  done_cols <- character(0)
  canon <- function(col) {
    if (col %in% done_cols) return(invisible())
    if (identical(col, .o("recist"))) d[[col]] <<- ctdna_stratify_recist(d[[col]], scheme)
    d[[col]] <<- .ctdna_factor_canonical(d[[col]], kind = col)
    done_cols <<- c(done_cols, col)
  }
  if (!is.null(group_by)) canon(group_by)
  for (f in facet) canon(f)

  # drop rows with NA in any grouping/facet variable (matches baseline/reduction)
  na_cols <- c(if (!is.null(group_by)) group_by, facet)
  if (length(na_cols)) {
    keep <- stats::complete.cases(d[, na_cols, drop = FALSE])
    d <- d[keep, , drop = FALSE]
    for (col in na_cols) if (is.factor(d[[col]])) d[[col]] <- droplevels(d[[col]])
  }

  # panel id (one cell per facet combination) for stats/tables
  d$.panel <- if (length(facet))
    interaction(d[, facet, drop = FALSE], sep = " | ", drop = TRUE, lex.order = TRUE)
  else factor("all")

  # median trajectory per (facet cell x group)
  by_cols <- c(t_col, facet, if (!is.null(group_by)) group_by)
  med <- stats::aggregate(d$.y, by = d[by_cols], FUN = stats::median, na.rm = TRUE)
  names(med)[ncol(med)] <- ".m"
  med$.mg <- if (!is.null(group_by))
    interaction(med[, c(facet, group_by), drop = FALSE], drop = TRUE)
  else interaction(med[, c(facet), drop = FALSE], drop = TRUE)

  # ---- plot ----
  if (!is.null(group_by)) {
    p <- ggplot2::ggplot(d, ggplot2::aes(.data[[t_col]], .data$.y,
                                         group = .data[[s_col]], color = .data[[group_by]])) +
      ggplot2::geom_line(alpha = 0.40) +
      ggplot2::geom_point(size = point_size, alpha = point_alpha) +
      ggplot2::geom_line(data = med,
        ggplot2::aes(.data[[t_col]], .data$.m, group = .data$.mg, color = .data[[group_by]]),
        linewidth = 1.5, inherit.aes = FALSE) +
      (if (identical(group_by, .o("recist"))) ctdna_scale_recist() else ggplot2::scale_color_viridis_d()) +
      ggplot2::labs(color = group_by)
  } else {
    p <- ggplot2::ggplot(d, ggplot2::aes(.data[[t_col]], .data$.y, group = .data[[s_col]])) +
      ggplot2::geom_line(alpha = 0.25, color = "grey50") +
      ggplot2::geom_line(data = med,
        ggplot2::aes(.data[[t_col]], .data$.m, group = .data$.mg),
        color = .o("median_color"), linewidth = 1.5, inherit.aes = FALSE)
  }
  # facet: 1 var -> wrap, 2 vars -> grid
  if (length(facet) == 1L)
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[[facet]]), scales = scales)
  else if (length(facet) == 2L)
    p <- p + ggplot2::facet_grid(rows = ggplot2::vars(.data[[facet[1]]]),
                                 cols = ggplot2::vars(.data[[facet[2]]]), scales = scales)
  p <- p + .log_clearance(y_scale) + ctdna_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1, size = 7))
  p <- .ctdna_apply_limits(p, xlim, ylim)

  # ---- statistics ----
  within_panel <- NULL; between_panels <- NULL
  do_stats <- isTRUE(stat) && isTRUE(show_stats) && !is.null(group_by)
  if (do_stats) {
    # per-visit omnibus on the plot + within_panel table
    vis_rows <- list(); wp_rows <- list()
    for (pn in split(d, d$.panel, drop = TRUE)) {
      if (!nrow(pn)) next
      fvals <- pn[1, facet, drop = FALSE]
      for (tp in levels(d[[t_col]])) {
        sst <- .compute_box_stats(pn[pn[[t_col]] == tp, , drop = FALSE],
                                  group_col = group_by, value_col = ".y", p_adjust = p_adjust)
        if (is.null(sst)) next
        rr <- data.frame(.tt = tp, label = paste0("p=", .format_p_value(sst$overall$p)),
                         stringsAsFactors = FALSE)
        vis_rows[[length(vis_rows) + 1L]] <- if (length(facet)) cbind(fvals, rr, row.names = NULL) else rr
      }
      pst <- .compute_box_stats(pn, group_col = group_by, value_col = ".y", p_adjust = p_adjust)
      if (!is.null(pst)) {
        plab <- as.character(pn$.panel[1])
        wp_rows[[length(wp_rows) + 1L]] <- rbind(
          data.frame(panel = plab, comparison = paste0(pst$overall$name, " (overall)"),
                     p = pst$overall$p, stringsAsFactors = FALSE),
          data.frame(panel = plab, comparison = paste0(pst$pairwise$left, " vs ", pst$pairwise$right),
                     p = pst$pairwise$p, stringsAsFactors = FALSE))
      }
    }
    if (length(vis_rows)) {
      vdf <- do.call(rbind, vis_rows); names(vdf)[names(vdf) == ".tt"] <- t_col
      vdf[[t_col]] <- factor(vdf[[t_col]], levels = levels(d[[t_col]]))
      p <- p + ggplot2::geom_text(data = vdf,
        mapping = ggplot2::aes(x = .data[[t_col]], y = -Inf, label = .data$label),
        inherit.aes = FALSE, vjust = -0.8, size = 2.5, color = "grey20")
    }
    within_panel <- if (length(wp_rows)) do.call(rbind, wp_rows) else NULL

    # between_panels table (whole-trajectory)
    panels <- levels(droplevels(d$.panel))
    if (length(panels) >= 2L) {
      bp_rows <- list()
      for (pr in utils::combn(panels, 2L, simplify = FALSE)) {
        res <- .longitudinal_traj_test(d[d$.panel %in% pr, , drop = FALSE], t_col, s_col)
        bp_rows[[length(bp_rows) + 1L]] <- data.frame(
          panel_a = pr[1], panel_b = pr[2], level_p = res$level_p,
          shape_p = res$shape_p, method = res$method, stringsAsFactors = FALSE)
      }
      between_panels <- do.call(rbind, bp_rows)
      # between_panels is already pairwise (one row per panel pair); correct the
      # level_p / shape_p families across those pairs.
      padj <- .resolve_p_adjust(p_adjust)
      if (!is.null(between_panels) && padj != "none") {
        between_panels$level_p <- stats::p.adjust(between_panels$level_p, padj)
        between_panels$shape_p <- stats::p.adjust(between_panels$shape_p, padj)
        between_panels$p_adjust <- padj
      }
    }
  }

  if (show_n) {
    counts <- stats::aggregate(d[[s_col]], by = c(d[c(t_col, facet)]),
                               FUN = function(x) length(unique(x)))
    names(counts)[ncol(counts)] <- "n"; counts$label <- paste0("n=", counts$n)
    p <- p + ggplot2::geom_label(data = counts,
      mapping = ggplot2::aes(x = .data[[t_col]], y = Inf, label = .data$label),
      inherit.aes = FALSE, vjust = 1.2, size = 2.6,
      label.padding = ggplot2::unit(0.12, "lines"), label.size = 0,
      color = "grey15", fill = "white", fontface = "bold")
  }

  if (is.null(caption)) {
    caption <- sprintf("Floor=%g%%; thick line = median", .o("display_floor"))
    if (do_stats) {
      padj <- .resolve_p_adjust(p_adjust)
      caption <- paste0(caption,
        "\nP above each visit: groups compared at that visit (Kruskal-Wallis",
        if (padj != "none") paste0(", ", padj, "-adjusted pairwise") else "", ").",
        "\nPooled and between-panel comparisons: longitudinal_stats(p)")
    }
  }

  p <- .finalize(p, title, subtitle, caption, xlab, ylab, legend_position = legend_position)
  attr(p, "longitudinal_stats") <- list(within_panel = within_panel, between_panels = between_panels)
  class(p) <- c("ctdna_longitudinal", class(p))
  if (do_stats) message("Comparison tables available: longitudinal_stats(p)")
  p
}

#' Retrieve longitudinal comparison tables
#'
#' @param p A plot from \code{ctdna_plot_longitudinal()}.
#' @return A list with \code{$within_panel} and \code{$between_panels}.
#' @export
longitudinal_stats <- function(p) {
  s <- attr(p, "longitudinal_stats")
  if (is.null(s)) {
    message("No longitudinal comparison tables on this object.")
    return(invisible(NULL))
  }
  s
}

#' @export
print.ctdna_longitudinal <- function(x, ...) {
  NextMethod()
  if (!is.null(attr(x, "longitudinal_stats")))
    message("Comparison tables available: longitudinal_stats(p)")
  invisible(x)
}

# ---- D7: Continuous scatter (generalized; absorbs reduction_vs_tumor, mut_methyl)

# Correlation for all three methods with confidence intervals. Pearson CI comes
# from cor.test; Spearman/Kendall use a Fisher-z approximation (flagged approx).
.scatter_cor_all <- function(x, y, conf = 0.95) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- sum(ok)
  z <- stats::qnorm((1 + conf) / 2)
  one <- function(m) {
    if (n < 3L) return(data.frame(method = m, estimate = NA_real_, p = NA_real_,
                                  ci_low = NA_real_, ci_high = NA_real_,
                                  ci_approx = NA, n = n, stringsAsFactors = FALSE))
    ct <- suppressWarnings(stats::cor.test(x, y, method = m, exact = FALSE, conf.level = conf))
    r  <- unname(ct$estimate); approx <- FALSE
    ci <- ct$conf.int
    if (is.null(ci)) {                       # spearman/kendall: Fisher-z approx
      approx <- TRUE
      if (is.finite(r) && abs(r) < 1 && n > 4L) {
        se <- if (m == "spearman") sqrt((1 + r^2 / 2) / (n - 3))
              else                 sqrt(0.437 / (n - 4))
        ci <- tanh(atanh(r) + c(-1, 1) * z * se)
      } else ci <- c(NA_real_, NA_real_)
    }
    data.frame(method = m, estimate = r, p = ct$p.value,
               ci_low = ci[1], ci_high = ci[2], ci_approx = approx, n = n,
               stringsAsFactors = FALSE)
  }
  do.call(rbind, lapply(c("spearman","pearson","kendall"), one))
}

#' Scatter of two continuous quantities with correlation
#'
#' One scatter for any two continuous variables, with a rank/linear correlation,
#' optional regression line and confidence-interval shading. Replaces the former
#' `ctdna_plot_reduction_vs_tumor` and `ctdna_plot_mut_methyl`.
#'
#' @details
#' Axis tokens and how derived quantities are computed:
#' \itemize{
#'   \item \code{"maxvaf"}, \code{"meanvaf"}, \code{"methyl"}/\code{"tf"} — the
#'     per-sample ctDNA level columns, floored at \code{ctdna_opts("display_floor")}
#'     so clearance values stay on a log axis.
#'   \item \code{"reduction"} (a.k.a. \code{"reduction_ratio"}) — the per-subject
#'     \strong{best reduction ratio} = the deepest on-treatment value divided by the
#'     subject's baseline (C1D1) value, floored. 1 = no change, <1 = reduction,
#'     0.001 = clearance floor.
#'   \item \code{"tumor_change"} (a.k.a. \code{"best_tumor_change"},
#'     \code{"target_lesion_change"}, \code{"tumor"}) — the per-subject \strong{best
#'     percent change in the sum of target-lesion diameters} from baseline
#'     (RECIST). \eqn{-30\%} is the RECIST partial-response threshold.
#'   \item \code{"mr"} — molecular-response score; or any numeric column by name.
#' }
#' Categorical axes are rejected with a pointer to [ctdna_concordance()].
#'
#' @param prep A `ctdna_prep` (uses `$samples`) or a per-sample data frame.
#' @param x,y Axis tokens (see Details) or numeric column names. Defaults
#'   `x="maxvaf"`, `y="methyl"`.
#' @param grain `"auto"` (default), `"sample"`, or `"subject"`. `"auto"` uses
#'   subject grain when either axis is `reduction`/`tumor_change`, else sample.
#' @param color_by Categorical colour variable (e.g. `"RECIST"`, `"Dose"`);
#'   `NULL` for a single colour.
#' @param shape_by Categorical shape variable (e.g. `"RECIST"`); `NULL` for one shape.
#' @param facet One variable (`facet_wrap`) or a length-2 vector (`facet_grid`),
#'   e.g. `"Visit_name"`; `NULL` for no faceting.
#' @param scheme RECIST stratification used when `color_by`/`shape_by`/`facet` is RECIST.
#' @param stat Correlation reported on the plot: `"spearman"` (default),
#'   `"pearson"`, `"kendall"`, or `"none"`.
#' @param conf Confidence level for the reported correlation interval (default 0.95).
#' @param regression Draw the linear (lm) regression line (default `FALSE`).
#' @param ci_band Shade the regression confidence band (default `FALSE`).
#' @param x_scale,y_scale Axis transform per axis: `"log10"` (default) /
#'   `"linear"` / `"log2"` / `"log"`. Log axes get minor tick marks.
#' @param abline Draw the identity (slope-1) reference line (default `FALSE`).
#' @param hline,vline Optional reference-line intercepts (e.g. `hline = -30`).
#' @param discord Colour points concordant/discordant by log10 distance >
#'   `ctdna_opts("discord_log10")` (default `FALSE`).
#' @param show_stats Annotate the correlation in the caption (default `TRUE`).
#' @param title,subtitle,caption Plot labels (auto-derived if `NULL`).
#' @param xlab,ylab Axis labels (auto-derived, with a `(log)` suffix on log axes, if `NULL`).
#' @param point_size,point_alpha Point aesthetics.
#' @param legend_position One of `"right"`/`"left"`/`"top"`/`"bottom"`/`"none"`.
#' @return A ggplot scatter. Correlations (all three methods, with CIs) are
#'   retrievable with [scatter_stats()].
#' @seealso [scatter_stats()], [ctdna_concordance()]
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 40, seed = 1)
#' prep <- ctdna_prepare(infinity_report = sim$infinity_report,
#'                       tf_change = sim$tf_change, clinical = sim$clinical, verbose = FALSE)
#' ctdna_plot_scatter(prep)                                       # maxVAF vs methylTF
#' p <- ctdna_plot_scatter(prep, x = "reduction", y = "tumor_change",
#'                         x_scale = "log10", y_scale = "linear", hline = -30,
#'                         color_by = "Dose", shape_by = "RECIST",
#'                         regression = TRUE, ci_band = TRUE)
#' scatter_stats(p)                                               # all correlations + CIs
#' @export
ctdna_plot_scatter <- function(prep,
                               x = "maxvaf",
                               y = "methyl",
                               grain = c("auto","sample","subject"),
                               color_by = NULL,
                               shape_by = NULL,
                               facet = NULL,
                               scheme = c("three","two","four","four_alt","five","six"),
                               stat = c("spearman","pearson","kendall","none"),
                               conf = 0.95,
                               regression = FALSE,
                               ci_band = FALSE,
                               x_scale = "log10",
                               y_scale = "log10",
                               abline = FALSE,
                               hline = NULL,
                               vline = NULL,
                               xlim = NULL,
                               ylim = NULL,
                               discord = FALSE,
                               show_stats = .o("show_stats"),
                               title = NULL,
                               subtitle = NULL,
                               caption = NULL,
                               xlab = NULL,
                               ylab = NULL,
                               point_size = .o("point_size"),
                               point_alpha = .o("point_alpha"),
                               legend_position = .o("legend_position")) {
  scheme <- match.arg(scheme); stat <- match.arg(stat); grain <- match.arg(grain)
  df <- .ctdna_augment_frame(.ctdna_grain(prep, "samples"), prep)
  red_tokens  <- c("reduction","reduction_ratio","best_ratio")
  cpc_tokens  <- c("ctdna_change","pct_change","ctdna_pct_change")
  les_tokens  <- unique(c("tumor_size_change","lesion_change","sld_change","lesion_size_change",
                          "best_lesion_change","best_tumor_size_change","pchg",
                          tolower(.o("lesion_pct"))))
  tum_tokens  <- c("tumor","tumor_change","best_tumor_change","target_lesion_change","best_change")
  subj_tokens <- c(red_tokens, cpc_tokens, les_tokens, tum_tokens)
  uses_subj_tok <- tolower(x) %in% subj_tokens || tolower(y) %in% subj_tokens
  # v0.57.0: modality axis tokens (expression:/ihc:/tmb/mutation:/alteration:) are
  # per-subject -> resolved against prep$rnaseq/ihc/dnaseq/variants (subsumes the
  # ctdna_plot_vs_* helpers). They force subject grain.
  mod_x <- .ctdna_modality_resolve(prep, x, scheme)
  mod_y <- .ctdna_modality_resolve(prep, y, scheme)
  uses_modality <- !is.null(mod_x) || !is.null(mod_y)
  # resolve facet names early (against the raw sample frame): a per-visit facet
  # forces per-sample grain (reduction computed at each visit, not "best").
  fct <- if (is.null(facet) || !length(facet)) character(0)
         else unique(vapply(as.character(facet), function(f)
                tryCatch(.resolve_col(f, df, "facet"), error = function(e) NA_character_), character(1)))
  fct <- fct[!is.na(fct) & fct %in% names(df)]; if (length(fct) > 2L) fct <- fct[1:2]
  needs_visit <- .o("time") %in% fct
  if (grain == "auto")
    grain <- if (needs_visit) "sample" else if (uses_subj_tok || uses_modality) "subject" else "sample"
  dat <- if (grain == "subject") ctdna_summary(df)
         else if (uses_subj_tok) ctdna_ratio(df)   # per on-treatment sample: subject/Visit/Dose/RECIST/ratio
         else df
  if (!is.null(mod_x)) { ax <- .ctdna_attach_modality(prep, dat, x, scheme); dat <- ax$frame; x <- ax$token; if (is.null(xlab)) xlab <- ax$label }
  if (!is.null(mod_y)) { ay <- .ctdna_attach_modality(prep, dat, y, scheme); dat <- ay$frame; y <- ay$token; if (is.null(ylab)) ylab <- ay$label }
  if (grain == "sample") {                 # per-sample reduction/tumor for visit faceting
    if (tolower(x) %in% red_tokens || tolower(y) %in% red_tokens)
      dat$.red <- ctdna_floor_tf(dat$ratio)
    if (any(c(tolower(x), tolower(y)) %in% c(tum_tokens, les_tokens))) {
      sm <- ctdna_summary(df)
      if (.o("best_change") %in% names(sm))
        dat$.tum <- sm[[.o("best_change")]][match(dat[[.o("subject")]], sm[[.o("subject")]])]
      L <- .o("lesion"); LP <- .o("lesion_pct"); tt <- .o("time")
      if (LP %in% names(df)) {
        key <- paste(dat[[.o("subject")]], dat[[tt]])
        dat$.lp <- df[[LP]][match(key, paste(df[[.o("subject")]], df[[tt]]))]
      } else if (L %in% names(df)) {
        lb <- df[df[[tt]] == .o("baseline"), c(.o("subject"), L), drop = FALSE]; names(lb)[2] <- ".lb"
        lm <- merge(df[, c(.o("subject"), tt, L), drop = FALSE], lb, by = .o("subject"), all.x = TRUE)
        lm$.lp <- (lm[[L]] - lm$.lb) / lm$.lb * 100
        key <- paste(dat[[.o("subject")]], dat[[tt]])
        dat$.lp <- lm$.lp[match(key, paste(lm[[.o("subject")]], lm[[tt]]))]
      }
    }
  }

  meas_col <- function(tok) switch(tolower(tok),
    maxvaf = .o("maxvaf"), meanvaf = .o("meanvaf"),
    methyl = , tf = .o("tf"), mr = .o("mr"), NA_character_)

  resolve_axis <- function(tok) {
    t <- tolower(tok)
    if (t %in% red_tokens) {
      v   <- if (grain == "subject") ctdna_floor_tf(dat$best_ratio) else dat$.red
      if (is.null(v) || !length(v))
        stop("Axis '", tok, "' needs a reduction ratio, but none is available in this prep ",
             "(no on-treatment vs baseline ctDNA). Check that tf_change / longitudinal samples are present.",
             call. = FALSE)
      lbl <- if (grain == "subject") "Best reduction ratio" else "ctDNA change from baseline"
      return(list(v = v, label = lbl, floored = TRUE))
    }
    if (t %in% tum_tokens || t %in% les_tokens) {
      prefer_reported <- t %in% tum_tokens
      rep_v <- if (grain == "subject") dat[[.o("best_change")]] else dat$.tum
      der_v <- if (grain == "subject") dat$best_lesion_change   else dat$.lp
      has_rep <- !is.null(rep_v) && any(is.finite(rep_v))
      has_der <- !is.null(der_v) && any(is.finite(der_v))
      pick <- if (prefer_reported) (if (has_rep) "rep" else if (has_der) "der" else "none")
              else                 (if (has_der) "der" else if (has_rep) "rep" else "none")
      if (pick == "none")
        stop("Axis '", tok, "' needs tumour-size data, which is absent from this prep. ",
             "Pass your ADaM ADTR frame via ctdna_prepare(adtr = <ADTR>) (per-visit PCHG / ",
             "best = nadir), or supply a reported '", .o("best_change"), "' column. ",
             "(For the ctDNA percent change, use 'ctdna_change'.)", call. = FALSE)
      v   <- if (pick == "rep") rep_v else der_v
      lbl <- if (pick == "rep") "Best % change in target lesions"
             else if (grain == "subject") "Best % change in lesion size"
             else "Lesion size % change from baseline"
      return(list(v = v, label = lbl, floored = FALSE))
    }
    if (t %in% c("ctdna_change","pct_change","ctdna_pct_change")) {
      r <- if (grain == "subject") dat$best_ratio else dat$ratio
      if (is.null(r) || !length(r))
        stop("Axis '", tok, "' needs a reduction ratio to derive ctDNA % change; none available.",
             call. = FALSE)
      return(list(v = (r - 1) * 100,
                  label = if (grain == "subject") "Best ctDNA % change from baseline"
                          else "ctDNA % change from baseline", floored = FALSE))
    }
    mc <- meas_col(tok)
    col <- if (!is.na(mc) && mc %in% names(dat)) mc else .resolve_col(tok, dat, what = "axis")
    v <- dat[[col]]
    if (!is.numeric(v))
      stop(sprintf(paste0("Axis '%s' is categorical. ctdna_plot_scatter handles only ",
           "continuous variables; use ctdna_concordance() for categorical comparisons."),
           col), call. = FALSE)
    floored <- col %in% c(.o("maxvaf"), .o("meanvaf"), .o("tf"))
    list(v = if (floored) ctdna_floor_tf(v) else v, label = col, floored = floored)
  }

  # v0.55.0: classify each axis as continuous/categorical and route mixed or
  # categorical combinations to the correct tool with an emphasised, user-facing
  # message (so it reads as "use the right function", not a library failure).
  axis_kind <- function(tok) {
    t <- tolower(tok)
    if (t %in% c(red_tokens, cpc_tokens, les_tokens, tum_tokens)) return("cont")
    mc <- meas_col(tok)
    col <- if (!is.na(mc) && mc %in% names(dat)) mc
           else tryCatch(.resolve_col(tok, dat, "axis"), error = function(e) NA_character_)
    if (is.na(col)) return("cont")          # unknown -> let resolve_axis raise the precise error
    if (is.numeric(dat[[col]])) "cont" else "cat"
  }
  kx <- axis_kind(x); ky <- axis_kind(y)
  if (kx == "cat" || ky == "cat") {
    if (kx == "cat" && ky == "cat")
      stop(sprintf(paste0(">> ctdna_plot_scatter() plots TWO CONTINUOUS axes - here BOTH '%s' and '%s' ",
                   "are categorical. <<\n   For a category-vs-category comparison use:  ",
                   "ctdna_concordance(prep, x = \"%s\", y = \"%s\")"), x, y, x, y), call. = FALSE)
    catv  <- if (kx == "cat") x else y
    contv <- if (kx == "cat") y else x
    stop(sprintf(paste0(">> ctdna_plot_scatter() needs TWO CONTINUOUS axes - '%s' is categorical. <<\n",
                 "   For a continuous-vs-categorical comparison use a boxplot:  ",
                 "ctdna_boxplot(prep, y = \"%s\", x = \"%s\")"), catv, contv, catv), call. = FALSE)
  }

  xa <- resolve_axis(x); ya <- resolve_axis(y)
  d <- data.frame(.x = xa$v, .y = ya$v)
  # v0.55.0: when both axes carry the same unit (both % change, or both ratios)
  # and the caller gave no explicit limits, share one square range on both axes.
  axis_unit <- function(tok, lbl) {
    t <- tolower(tok)
    if (t %in% c("reduction_ratio", "best_ratio") ||
        grepl("ratio", lbl, ignore.case = TRUE)) return("ratio")
    if (t %in% c("reduction", cpc_tokens, les_tokens, tum_tokens) ||
        grepl("%|percent|pct|change", lbl, ignore.case = TRUE)) return("pct")
    "other"
  }
  ux <- axis_unit(x, xa$label); uy <- axis_unit(y, ya$label)
  eq_limits <- FALSE
  if (is.null(xlim) && is.null(ylim) && ux == uy && ux %in% c("pct", "ratio")) {
    rng <- suppressWarnings(range(c(d$.x, d$.y), na.rm = TRUE))
    log_any <- grepl("log", x_scale) || grepl("log", y_scale)
    if (all(is.finite(rng)) && diff(rng) > 0 && (!log_any || rng[1] > 0)) {
      pad <- 0.04 * diff(rng)
      xlim <- ylim <- c(rng[1] - pad, rng[2] + pad)
      eq_limits <- TRUE
    }
  }
  if (is.null(xlab)) xlab <- if (grepl("log", x_scale)) paste0(xa$label, " (log)") else xa$label
  if (is.null(ylab)) ylab <- if (grepl("log", y_scale)) paste0(ya$label, " (log)") else ya$label

  add_cat <- function(col) {
    cc <- .resolve_col(col, dat, what = "variable")
    v  <- dat[[cc]]
    if (identical(cc, .o("recist"))) v <- ctdna_stratify_recist(v, scheme)
    else v <- .ctdna_factor_canonical(v, kind = cc)
    list(v = v, label = cc)
  }
  ca <- if (!is.null(color_by)) add_cat(color_by) else NULL
  sa <- if (!is.null(shape_by)) add_cat(shape_by) else NULL
  if (!is.null(ca)) d$.col <- ca$v
  if (!is.null(sa)) d$.shp <- sa$v
  for (f in fct) {
    fv <- dat[[f]]
    d[[f]] <- if (identical(f, .o("recist"))) ctdna_stratify_recist(fv, scheme)
              else .ctdna_factor_canonical(fv, kind = f)
  }
  if (isTRUE(discord))
    d$.disc <- abs(log10(pmax(d$.x, .Machine$double.eps)) -
                   log10(pmax(d$.y, .Machine$double.eps))) > .o("discord_log10")

  # correlation: all three methods (+ CIs), plus the chosen one for the caption
  cor_tab <- .scatter_cor_all(d$.x, d$.y, conf = conf)
  if (show_stats && stat != "none") {
    nm  <- c(spearman = "Spearman", pearson = "Pearson", kendall = "Kendall")[stat]
    row <- cor_tab[cor_tab$method == stat, ]
    sym <- if (stat == "kendall") "tau" else if (stat == "pearson") "r" else "rho"
    ci_txt <- if (all(is.finite(c(row$ci_low, row$ci_high))))
      sprintf(" (%g%% CI%s %.2f to %.2f)", 100 * conf,
              if (isTRUE(row$ci_approx)) " approx" else "", row$ci_low, row$ci_high) else ""
    lab <- sprintf("%s %s = %.2f%s, %s, n = %d", nm, sym, row$estimate, ci_txt,
                   ctdna_pval_label(row$p), row$n)
    gl <- switch(stat,
      spearman = "rho: monotonic rank correlation, -1 to 1 (0 = none)",
      pearson  = "r: linear correlation, -1 to 1 (0 = none)",
      kendall  = "tau: rank concordance, -1 to 1 (0 = none)")
    caption <- if (is.null(caption)) paste0(lab, "\n", gl) else paste0(caption, "\n", lab, "\n", gl)
    if (isTRUE(discord))
      caption <- paste0(caption, sprintf("\n%d discordant (log10 distance > %g)",
                                          sum(d$.disc, na.rm = TRUE), .o("discord_log10")))
  }

  # build
  aes_args <- list(x = quote(.data$.x), y = quote(.data$.y))
  if (isTRUE(discord) && is.null(ca)) aes_args$color <- quote(.data$.disc)
  else if (!is.null(ca))              aes_args$color <- quote(.data$.col)
  if (!is.null(sa))                   aes_args$shape <- quote(.data$.shp)
  p <- ggplot2::ggplot(d, do.call(ggplot2::aes, aes_args))
  if (isTRUE(abline)) p <- p + ggplot2::geom_abline(slope = 1, color = "grey50", linetype = "dashed")
  if (!is.null(hline)) p <- p + ggplot2::geom_hline(yintercept = hline, color = "red", linetype = "dashed")
  if (!is.null(vline)) p <- p + ggplot2::geom_vline(xintercept = vline, color = "grey50", linetype = "dotted")
  if (isTRUE(ci_band))                       # CI shading (ribbon), independent of the line
    p <- p + ggplot2::geom_smooth(mapping = ggplot2::aes(x = .data$.x, y = .data$.y),
             inherit.aes = FALSE, method = "lm", formula = y ~ x, se = TRUE,
             color = NA, fill = "grey60", alpha = 0.2)
  if (isTRUE(regression))
    p <- p + ggplot2::geom_smooth(mapping = ggplot2::aes(x = .data$.x, y = .data$.y),
             inherit.aes = FALSE, method = "lm", formula = y ~ x, se = FALSE,
             color = "grey25", linewidth = 0.7)
  p <- p + ggplot2::geom_point(size = point_size, alpha = point_alpha)
  if (isTRUE(discord) && is.null(ca))
    p <- p + ggplot2::scale_color_manual(values = c(`FALSE` = "grey20", `TRUE` = "#E31A1C"),
                                         labels = c("concordant","discordant"), name = NULL)
  else if (!is.null(ca))
    p <- p + (if (identical(ca$label, .o("recist"))) ctdna_scale_recist()
              else ggplot2::scale_color_viridis_d()) + ggplot2::labs(color = ca$label)
  if (!is.null(sa)) p <- p + ggplot2::labs(shape = sa$label)
  if (length(fct) == 1L) p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[[fct]]))
  else if (length(fct) == 2L)
    p <- p + ggplot2::facet_grid(rows = ggplot2::vars(.data[[fct[1]]]), cols = ggplot2::vars(.data[[fct[2]]]))
  p <- p + .ctdna_x_log_layers(x_scale) + .ctdna_resolve_y_scale(y_scale) + ctdna_theme()

  if (is.null(title)) title <- sprintf("%s vs %s", ya$label, xa$label)
  p <- .finalize(p, title, subtitle, caption, xlab, ylab, legend_position = legend_position)
  p <- if (eq_limits)
         p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE, clip = "off")
       else .ctdna_apply_limits(p, xlim, ylim)
  attr(p, "scatter_stats") <- cor_tab
  p
}

#' Retrieve scatter correlation statistics
#'
#' Returns the Spearman, Pearson and Kendall correlations (estimate, p-value,
#' confidence interval, and n) computed for the last [ctdna_plot_scatter()] call —
#' so you do not have to dig into plot attributes.
#'
#' @param p A plot from [ctdna_plot_scatter()].
#' @return A data frame with one row per method (`method`, `estimate`, `p`,
#'   `ci_low`, `ci_high`, `ci_approx`, `n`).
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 30, seed = 1)
#' prep <- ctdna_prepare(infinity_report = sim$infinity_report,
#'                       tf_change = sim$tf_change, clinical = sim$clinical, verbose = FALSE)
#' scatter_stats(ctdna_plot_scatter(prep))
#' @export
scatter_stats <- function(p) {
  s <- attr(p, "scatter_stats")
  if (is.null(s)) { message("No scatter statistics on this object."); return(invisible(NULL)) }
  s
}

#' Retrieve concordance statistics
#'
#' Returns the agreement statistics from the last [ctdna_concordance()] call
#' (contingency table, % agreement, kappa, weighted kappa, McNemar p,
#' chi-square p, Cramer's V, N) — a friendlier alternative to reading plot
#' attributes.
#'
#' @param p A plot from [ctdna_concordance()].
#' @return A named list of the concordance statistics and the contingency table.
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 30, seed = 1)
#' prep <- ctdna_prepare(infinity_report = sim$infinity_report,
#'                       tf_change = sim$tf_change, clinical = sim$clinical, verbose = FALSE)
#' concordance_stats(ctdna_concordance(prep))
#' @export
concordance_stats <- function(p) {
  s <- attr(p, "concordance")
  if (is.null(s)) { message("No concordance statistics on this object."); return(invisible(NULL)) }
  s
}

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
#' @return A ggplot object.
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 30, seed = 1)
#' prep <- ctdna_prepare(sim, verbose = FALSE)
#' d <- prep$samples
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
    xlim = NULL, ylim = NULL,
    legend_position = .o("legend_position")) {
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
  .yl <- if (!is.null(ylim)) ylim else if (!is.null(y_clip) && length(y_clip) == 2L) y_clip else NULL
  p <- .ctdna_apply_limits(p, xlim, .yl)

  .finalize(p,
             title  %||% "Percent change by dose and time",
             subtitle, caption,
             xlab %||% dose_col,
             ylab %||% value_col,
             legend_position = legend_position)
}


# ---- General continuous-by-categorical boxplot --------------------------------

#' Boxplot of any continuous variable across any categorical variable
#'
#' The continuous-by-categorical member of the comparison trio
#' (`ctdna_plot_scatter` = continuous-vs-continuous,
#' `ctdna_concordance` = categorical-vs-categorical, `ctdna_boxplot` =
#' continuous-vs-categorical). `y` (continuous) is boxed across `x`
#' (categorical), with the same Wilcoxon / Kruskal-Wallis + pairwise-bracket
#' engine used elsewhere. Any column from the per-sample (ctDNA) frame **or** the
#' clinical frame can be used for `y`, `x`, `color_by`, or `facet`; set `visit`
#' to a time point (e.g. `"C1D1"`) to reproduce the baseline-characteristics view.
#'
#' @details
#' `y` accepts the ctDNA tokens `"methyl"`/`"tf"`, `"maxvaf"`, `"meanvaf"`,
#' `"mr"`; the derived `"reduction"` (best ratio) and `"ctdna_change"`
#' (= \eqn{(ratio-1)\times100}, the ctDNA % change from baseline, always derivable
#' from the TF); the imaging-derived `"tumor_change"` (best % change in target
#' lesions, needs RECIST sum-of-diameters data); or any numeric column by name.
#' Categorical `y` is rejected (use [ctdna_concordance()]); continuous `x` is
#' rejected (use [ctdna_plot_scatter()]).
#'
#' @param prep A `ctdna_prep` (uses `$samples` + `$clinical`) or a data frame.
#' @param y Continuous variable: a token (see Details) or numeric column name.
#' @param x Categorical grouping variable (e.g. `"RECIST"`, `"Dose"`, `"Arm"`).
#' @param color_by Optional categorical sub-group (dodged); `NULL` for none.
#' @param facet One variable (`facet_wrap`) or a length-2 vector (`facet_grid`).
#' @param visit Optional visit label(s) to subset to (e.g. `"C1D1"`, `"baseline"`).
#' @param compare_by `"group"` (x levels, default), `"subgroup"` (color levels),
#'   `"subgroup_within_group"`, or `"group_within_subgroup"`.
#' @param pairwise Draw pairwise brackets for the within-modes (default `TRUE`).
#' @param scheme RECIST stratification when `x`/`color_by`/`facet` is RECIST.
#' @param stat Reserved for test selection (currently non-parametric throughout).
#' @param p_adjust Multiple-test correction for pairwise families
#'   (`FALSE`/`TRUE`/method name; see [ctdna_plot_baseline]).
#' @param y_scale `"linear"` (default) / `"log10"` / ... (log adds minor ticks).
#' @param ylim Optional y-limits.
#' @param show_stats,show_n Display toggles.
#' @param title,subtitle,caption,xlab,ylab,point_size,point_alpha,legend_position Display.
#' @return A ggplot boxplot.
#' @seealso [ctdna_plot_scatter()], [ctdna_concordance()], [ctdna_plot_baseline()]
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 40, seed = 1)
#' prep <- ctdna_prepare(infinity_report = sim$infinity_report,
#'                       tf_change = sim$tf_change, clinical = sim$clinical, verbose = FALSE)
#' ctdna_boxplot(prep, y = "methyl", x = "RECIST")                 # methylTF by RECIST
#' ctdna_boxplot(prep, y = "TMB_score", x = "Dose")               # any numeric column
#' ctdna_boxplot(prep, y = "methyl", x = "Dose", color_by = "RECIST",
#'               visit = "C1D1")                                    # ~ baseline characteristics
#' ctdna_boxplot(prep, y = "ctdna_change", x = "RECIST", visit = "C2D1")
#' @export
ctdna_boxplot <- function(prep,
                          y, x,
                          color_by = NULL,
                          facet = NULL,
                          visit = NULL,
                          compare_by = c("group","subgroup","subgroup_within_group","group_within_subgroup"),
                          pairwise = TRUE,
                          scheme = c("three","two","four","four_alt","five","six"),
                          stat = c("auto","wilcox","t","none"),
                          p_adjust = FALSE,
                          y_scale = "linear",
                          ylim = NULL,
                          xlim = NULL,
                          show_stats = .o("show_stats"),
                          show_n = .o("show_n"),
                          title = NULL,
                          subtitle = NULL,
                          caption = NULL,
                          xlab = NULL,
                          ylab = NULL,
                          point_size = .o("point_size"),
                          point_alpha = .o("point_alpha"),
                          legend_position = .o("legend_position")) {
  compare_by <- match.arg(compare_by); scheme <- match.arg(scheme); stat <- match.arg(stat)
  s <- .o("subject"); t <- .o("time")
  # v0.56.0: per-visit tumour size lives in $assessments, not $samples. When y is a
  # per-visit tumour token, draw from the tumour assessments (Patient_ID x Visit).
  .pervisit_tumor <- tolower(y) %in% c("tumor_size_change","lesion_change","sld_change",
      "lesion_size_change","pchg", tolower(.o("lesion_pct")))
  .tl <- if (.pervisit_tumor) .ctdna_tumor_long(prep) else NULL
  D <- if (!is.null(.tl) && nrow(.tl)) .tl
       else .ctdna_augment_frame(.ctdna_grain(prep, "samples"), prep)

  # v0.57.0: modality tokens (expression:/ihc:/tmb/mutation:/alteration:) -> columns
  # on D keyed by subject (subsumes ctdna_plot_vs_mutation / _vs_ihc, etc.).
  for (.nm in c("y", "x", "color_by")) {
    .tok <- get(.nm); if (is.null(.tok)) next
    .am <- .ctdna_attach_modality(prep, D, .tok, scheme)
    if (isTRUE(.am$modality)) { D <- .am$frame; assign(.nm, .am$token) }
  }
  if (!is.null(facet) && length(facet))
    facet <- vapply(as.character(facet), function(ft) {
      am <- .ctdna_attach_modality(prep, D, ft, scheme)
      if (isTRUE(am$modality)) { D <<- am$frame; am$token } else ft }, character(1))

  # bring in clinical columns (per subject) so x/y/color/facet can be clinical
  cl <- if (inherits(prep, "ctdna_prep")) prep$clinical else NULL
  ensure_col <- function(col) {
    if (col %in% names(D)) return(invisible())
    if (!is.null(cl) && col %in% names(cl))
      D[[col]] <<- cl[[col]][match(D[[s]], cl[[s]])]
  }

  # derived columns: per-visit (.cpc ctDNA % change, .lp lesion % change) and
  # subject-level best (.best_ratio, .best_cpc, .best_lesion, .tum)
  rr <- tryCatch(ctdna_ratio(D), error = function(e) NULL)
  if (!is.null(rr)) {
    key <- paste(D[[s]], D[[t]])
    D$.ratio <- rr$ratio[match(key, paste(rr[[s]], rr[[t]]))]
    D$.cpc   <- (D$.ratio - 1) * 100
  }
  L <- .o("lesion"); LP <- .o("lesion_pct")
  if (LP %in% names(D)) {
    D$.lp <- suppressWarnings(as.numeric(D[[LP]]))
  } else if (L %in% names(D)) {
    base_map <- tapply(D[[L]][D[[t]] == .o("baseline")],
                       as.character(D[[s]])[D[[t]] == .o("baseline")], function(z) z[1])
    bl_les <- base_map[as.character(D[[s]])]
    D$.lp  <- (D[[L]] - bl_les) / bl_les * 100
  }
  sm <- tryCatch(ctdna_summary(D), error = function(e) NULL)
  if (!is.null(sm)) {
    bm <- function(col) if (col %in% names(sm)) sm[[col]][match(D[[s]], sm[[s]])] else NULL
    D$.best_ratio  <- bm("best_ratio")
    if (!is.null(D$.best_ratio)) D$.best_cpc <- (D$.best_ratio - 1) * 100
    D$.tum         <- bm(.o("best_change"))
    D$.best_lesion <- bm("best_lesion_change")
  }

  meas_col <- function(tok) switch(tolower(tok),
    maxvaf = .o("maxvaf"), meanvaf = .o("meanvaf"),
    methyl = , tf = .o("tf"), mr = .o("mr"), NA_character_)

  resolve_y <- function(tok) {
    tl <- tolower(tok)
    spec <-
      if (tl %in% c("reduction","reduction_ratio","best_ratio")) list(".best_ratio","Best reduction ratio",TRUE)
      else if (tl %in% c("best_ctdna_change")) list(".best_cpc","Best ctDNA % change from baseline",TRUE)
      else if (tl %in% c("ctdna_change","pct_change","ctdna_pct_change")) list(".cpc","ctDNA % change from baseline",FALSE)
      else if (tl %in% c("tumor_size_change","lesion_change","sld_change","lesion_size_change",
                          "pchg", tolower(.o("lesion_pct")))) list(".lp","Lesion size % change from baseline",FALSE)
      else if (tl %in% c("best_tumor_size_change","best_lesion_change")) list(".best_lesion","Best % change in lesion size",TRUE)
      else if (tl %in% c("tumor","tumor_change","best_tumor_change","target_lesion_change","best_change")) {
        if (!is.null(D$.tum) && any(is.finite(D$.tum))) list(".tum","Best % change in target lesions",TRUE)
        else list(".best_lesion","Best % change in lesion size",TRUE)   # fall back to derived
      }
      else NULL
    if (!is.null(spec)) { col <- spec[[1]]; lbl <- spec[[2]]; subj <- spec[[3]] }
    else {
      mc <- meas_col(tok)
      col <- if (!is.na(mc) && mc %in% names(D)) mc else { ensure_col(.resolve_col(tok, D, "y")); .resolve_col(tok, D, "y") }
      lbl <- col; subj <- FALSE
    }
    v <- D[[col]]
    if (is.null(v) || !length(v) || all(is.na(v))) {
      hint <- if (col %in% c(".lp",".best_lesion",".tum"))
        paste0(" Pass your ADaM ADTR frame via ctdna_prepare(adtr = <ADTR>) (per-visit PCHG / ",
               "best = nadir), or supply a reported '", .o("best_change"), "' column.") else ""
      stop("y = '", tok, "' is not available in this prep.", hint, call. = FALSE)
    }
    if (!is.numeric(v))
      stop(sprintf(paste0(">> ctdna_boxplot() needs a CONTINUOUS y - '%s' is categorical. <<\n",
           "   If your x is continuous, swap them so the continuous variable is y:  ",
           "ctdna_boxplot(prep, y = \"<continuous>\", x = \"%s\").\n",
           "   For a category-vs-category comparison use:  ctdna_concordance(prep, x = ..., y = \"%s\")."),
           tok, tok, tok), call. = FALSE)
    list(col = col, label = lbl, subject = subj)
  }

  resolve_cat <- function(col, role) {
    cc <- .resolve_col(col, D, role); ensure_col(cc)
    if (!cc %in% names(D)) ensure_col(cc)
    v <- D[[cc]]
    if (is.numeric(v) && length(unique(v[is.finite(v)])) > 6L)
      stop(role, " = '", col, "' looks CONTINUOUS. ctdna_boxplot needs a categorical ", role,
           "; for continuous-vs-continuous use ctdna_plot_scatter().", call. = FALSE)
    cc
  }

  ya <- resolve_y(y)
  # subject-level y (best_*): one row per subject (best is across visits, so visit= is ignored);
  # per-visit y: apply the visit subset.
  if (isTRUE(ya$subject)) {
    D <- D[!duplicated(as.character(D[[s]])), , drop = FALSE]
  } else if (!is.null(visit)) {
    vt <- as.character(visit); vt[tolower(vt) %in% c("baseline","bl")] <- .o("baseline")
    D <- D[.ctdna_normalize_visit(D[[t]]) %in% .ctdna_normalize_visit(vt), , drop = FALSE]
    if (!nrow(D)) stop("No samples at visit(s): ", paste(visit, collapse = ", "), call. = FALSE)
  }
  xcol <- resolve_cat(x, "x")
  ccol <- if (!is.null(color_by)) resolve_cat(color_by, "color_by") else NULL
  fcols <- if (is.null(facet) || !length(facet)) character(0)
           else vapply(as.character(facet), function(f) resolve_cat(f, "facet"), character(1))
  fcols <- unique(fcols); if (length(fcols) > 2L) fcols <- fcols[1:2]

  canon <- function(col) if (identical(col, .o("recist"))) ctdna_stratify_recist(D[[col]], scheme)
                         else .ctdna_factor_canonical(D[[col]], kind = col)
  pdat <- data.frame(.y = D[[ya$col]], .x = canon(xcol), stringsAsFactors = FALSE)
  if (!is.null(ccol)) pdat$.col <- canon(ccol)
  for (f in fcols) pdat[[f]] <- canon(f)
  pdat <- pdat[is.finite(pdat$.y) & !is.na(pdat$.x), , drop = FALSE]
  use_log <- grepl("log", y_scale)
  if (use_log) pdat <- pdat[pdat$.y > 0, , drop = FALSE]
  if (!nrow(pdat)) stop("No finite data to plot after filtering.", call. = FALSE)

  if (is.null(xlab)) xlab <- xcol
  if (is.null(ylab)) ylab <- if (use_log) paste0(ya$label, " (log)") else ya$label
  if (is.null(title)) title <- sprintf("%s by %s%s", ya$label, xcol,
                                        if (!is.null(visit)) paste0(" (", paste(visit, collapse = ", "), ")") else "")

  # plot
  dodge <- !is.null(ccol)
  fill_by_x <- is.null(ccol)            # v0.55.0: colour boxes by x when no color_by
  aes0 <- if (dodge) ggplot2::aes(.data$.x, .data$.y, fill = .data$.col)
          else if (fill_by_x) ggplot2::aes(.data$.x, .data$.y, fill = .data$.x)
          else ggplot2::aes(.data$.x, .data$.y)
  p <- ggplot2::ggplot(pdat, aes0) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.3, linewidth = 0.7, fatten = 2.2,
                          width = .o("box_width"),
                          position = ggplot2::position_dodge(0.8, preserve = "single"))
  p <- p + (if (dodge)
              ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
                                  size = point_size, alpha = point_alpha)
            else ggplot2::geom_jitter(width = 0.15, size = point_size, alpha = point_alpha))
  if (length(fcols) == 1L) p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[[fcols]]))
  else if (length(fcols) == 2L)
    p <- p + ggplot2::facet_grid(rows = ggplot2::vars(.data[[fcols[1]]]), cols = ggplot2::vars(.data[[fcols[2]]]))
  if (!is.null(ccol)) p <- p + (if (identical(ccol, .o("recist"))) ctdna_scale_recist_fill()
                                else ggplot2::scale_fill_viridis_d()) + ggplot2::labs(fill = ccol)
  else if (fill_by_x) p <- p + (if (identical(xcol, .o("recist"))) ctdna_scale_recist_fill()
                                else ggplot2::scale_fill_viridis_d()) +
                               ggplot2::guides(fill = "none")   # legend redundant with x axis
  p <- p + .ctdna_resolve_y_scale(y_scale) + ctdna_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, vjust = 1, hjust = 1))
  cc_args <- list(clip = "off")
  if (!is.null(ylim) && !(compare_by %in% c("subgroup_within_group","group_within_subgroup"))) cc_args$ylim <- ylim
  if (!is.null(xlim)) cc_args$xlim <- xlim
  p <- p + do.call(ggplot2::coord_cartesian, cc_args)

  # stats (reuse the shared engine; mirrors baseline/reduction)
  if (isTRUE(show_stats)) {
    fdata <- if (length(fcols)) split(pdat, interaction(pdat[fcols], drop = TRUE)) else list(all = pdat)
    if (compare_by == "group") {
      for (nm in names(fdata)) {
        sub <- fdata[[nm]]
        st <- .compute_box_stats(sub, ".x", ".y", p_adjust = p_adjust)
        if (!is.null(st)) p <- .add_brackets(p, st, sub, ".x", ".y",
          y_scale = y_scale, facet_data = if (length(fcols)) sub[1, fcols, drop = FALSE] else NULL)
      }
      caption <- paste0(caption %||% "", "\nCompare ", xcol, if (p_adjust_on(p_adjust)) " (pairwise adjusted)" else "")
    } else if (compare_by == "subgroup" && !is.null(ccol)) {
      ng <- nlevels(droplevels(pdat$.x))
      for (nm in names(fdata)) {
        sub <- fdata[[nm]]; st <- .compute_box_stats(sub, ".col", ".y", p_adjust = p_adjust)
        if (is.null(st)) next
        lab <- data.frame(.x = factor(levels(pdat$.x)[ceiling(ng/2)], levels = levels(pdat$.x)),
                          .y = Inf, lab = sprintf("%s p=%s", ccol, .format_p_value(st$overall$p)),
                          stringsAsFactors = FALSE)
        if (length(fcols)) lab <- cbind(lab, sub[1, fcols, drop = FALSE])
        p <- p + ggplot2::geom_text(data = lab, mapping = ggplot2::aes(x = .data$.x, y = .data$.y, label = .data$lab),
          inherit.aes = FALSE, vjust = 1.5, size = 3.0, fontface = "bold", color = "grey15")
      }
    } else if (compare_by %in% c("subgroup_within_group","group_within_subgroup") && !is.null(ccol)) {
      mode <- if (compare_by == "subgroup_within_group") "within_group" else "within_subgroup"
      for (nm in names(fdata)) {
        sub <- fdata[[nm]]
        p <- .add_dodge_brackets(p, sub, group_col = ".x", subgroup_col = ".col", value_col = ".y",
          mode = mode, pairwise = pairwise, p_adjust = p_adjust, y_scale = y_scale,
          facet_data = if (length(fcols)) sub[1, fcols, drop = FALSE] else NULL)
      }
    }
  }

  if (show_n) {
    nlab <- stats::aggregate(rep(1L, nrow(pdat)), by = c(list(.x = pdat$.x), pdat[fcols]), FUN = length)
    names(nlab)[ncol(nlab)] <- "n"; nlab$lab <- paste0("n=", nlab$n)
    p <- p + ggplot2::geom_text(data = nlab,
      mapping = ggplot2::aes(x = .data$.x, y = -Inf, label = .data$lab),
      inherit.aes = FALSE, vjust = -0.6, size = 2.5, color = "grey30")
  }
  .finalize(p, title, subtitle, caption, xlab, ylab, legend_position = legend_position)
}

p_adjust_on <- function(x) .resolve_p_adjust(x) != "none"

# ---- v0.52.0: shared x/y limit application (data-space zoom; never drops points) ----
# coord_cartesian keeps all data and zooms; safe on log scales (limits in data units).
# Applied as the LAST coord, so it supersedes any internal clip a function set.
.ctdna_apply_limits <- function(p, xlim = NULL, ylim = NULL) {
  if (is.null(xlim) && is.null(ylim)) return(p)
  p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, clip = "off")
}

# ---- v0.53.0: make prep a truly linked object for plotting ------------------
# Per-visit tumour size (ADTR) is folded into $samples by ctdna_prepare, and the
# per-patient nadir/best into $clinical. Here we join $clinical (per patient)
# onto the per-sample frame by subject, so any plotting function can use ANY
# column from either frame as an axis/aesthetic.
.ctdna_augment_frame <- function(df, prep) {
  if (!is.list(prep) || is.data.frame(prep) || is.null(df)) return(df)
  subj <- .o("subject")
  cl <- prep$clinical
  if (!is.null(cl) && subj %in% names(cl) && subj %in% names(df)) {
    add <- setdiff(setdiff(names(cl), names(df)), subj)
    if (length(add)) {
      idx <- match(df[[subj]], cl[[subj]])
      for (cc in add) df[[cc]] <- cl[[cc]][idx]
    }
  }
  df
}

# v0.56.0: the per-visit tumour view now lives in $assessments (Patient_ID x
# Visit_name x param). Return its tumour rows (non-NA Tumor_size_pct_change),
# clinical-augmented, with a Tumor_visit_name alias so x="Tumor_visit_name" works.
.ctdna_tumor_long <- function(prep) {
  if (!is.list(prep) || is.data.frame(prep) || is.null(prep$assessments)) return(NULL)
  a <- prep$assessments; lp <- .o("lesion_pct")
  if (!(lp %in% names(a))) return(NULL)
  tum <- a[is.finite(suppressWarnings(as.numeric(a[[lp]]))), , drop = FALSE]
  if (!nrow(tum)) return(NULL)
  # tumour visit lives in Tumor_visit_name (raw AVISIT). Mirror it onto the generic
  # Visit_name so both x = "Tumor_visit_name" and x = "Visit_name" resolve in plots.
  if ("Tumor_visit_name" %in% names(tum)) tum[[.o("time")]] <- tum[["Tumor_visit_name"]]
  .ctdna_augment_frame(tum, prep)
}

# v0.57.0: resolve a MODALITY axis token against the prep's modality frames,
# WITHOUT changing any data structure. Subsumes ctdna_plot_vs_expression /
# _vs_ihc / _vs_tmb / _vs_mutation. Token grammar (per subject):
#   "expression:GENE" -> numeric, from prep$rnaseq      (cont)
#   "ihc:MARKER"      -> numeric, from prep$ihc          (cont)
#   "tmb"             -> numeric, from prep$dnaseq/$samples TMB (cont)
#   "mutation:GENE"   -> mut/wt,  from prep$variants      (cat)
#   "alteration:GENE" -> alteration type, from prep$variants (cat)
# Returns list(values = vector named by subject, label, kind) or NULL.
.ctdna_modality_resolve <- function(prep, token, scheme = "three") {
  if (!inherits(prep, "ctdna_prep") || is.null(token) || length(token) != 1) return(NULL)
  subj <- .o("subject")
  parts <- strsplit(as.character(token), ":", fixed = TRUE)[[1]]
  key <- tolower(parts[1]); arg <- if (length(parts) >= 2) paste(parts[-1], collapse = ":") else NA
  first_by <- function(val, by) tapply(val, as.character(by), function(z) z[!is.na(z)][1])
  if (key %in% c("expression", "expr") && !is.na(arg)) {
    e <- prep$rnaseq; if (is.null(e)) return(NULL)
    g <- .o("gene"); ec <- .o("expression")
    if (!all(c(g, ec, subj) %in% names(e))) return(NULL)
    sub <- e[as.character(e[[g]]) == arg, , drop = FALSE]
    if (!nrow(sub)) stop("expression: gene not found in $rnaseq: ", arg, call. = FALSE)
    return(list(values = first_by(as.numeric(sub[[ec]]), sub[[subj]]),
                label = sprintf("%s expression", arg), kind = "cont"))
  }
  if (key == "ihc" && !is.na(arg)) {
    h <- prep$ihc; if (is.null(h)) return(NULL)
    m <- .o("marker"); sc <- .o("ihc_score")
    if (!all(c(m, sc, subj) %in% names(h))) return(NULL)
    sub <- h[as.character(h[[m]]) == arg, , drop = FALSE]
    if (!nrow(sub)) stop("ihc: marker not found in $ihc: ", arg, call. = FALSE)
    return(list(values = first_by(as.numeric(sub[[sc]]), sub[[subj]]),
                label = sprintf("%s IHC score", arg), kind = "cont"))
  }
  if (key == "tmb") {
    tc <- .o("tmb"); src <- NULL
    for (fr in list(prep$dnaseq, prep$samples))
      if (!is.null(fr) && tc %in% names(fr) && subj %in% names(fr)) { src <- fr; break }
    if (is.null(src)) return(NULL)
    return(list(values = first_by(suppressWarnings(as.numeric(src[[tc]])), src[[subj]]),
                label = "TMB", kind = "cont"))
  }
  if (key %in% c("mutation", "mut", "alteration", "alt") && !is.na(arg)) {
    vfr <- prep$variants; if (is.null(vfr)) return(NULL)
    g <- .o("gene"); if (!all(c(g, subj) %in% names(vfr))) return(NULL)
    allsubj <- if (!is.null(prep$clinical) && subj %in% names(prep$clinical))
                 unique(as.character(prep$clinical[[subj]]))
               else unique(as.character(prep$samples[[subj]]))
    hit <- vfr[as.character(vfr[[g]]) == arg, , drop = FALSE]
    if (key %in% c("alteration", "alt")) {
      acol <- tryCatch(.o("alteration"), error = function(e) NULL)
      if (is.null(acol) || !nzchar(acol)) acol <- "alteration_type"
      v <- stats::setNames(rep("wt", length(allsubj)), allsubj)
      if (nrow(hit) && acol %in% names(hit)) {
        a <- tapply(as.character(hit[[acol]]), as.character(hit[[subj]]), function(z) z[1])
        v[names(a)] <- a
      }
      return(list(values = v, label = sprintf("%s alteration", arg), kind = "cat"))
    }
    mutsubj <- unique(as.character(hit[[subj]]))
    v <- stats::setNames(ifelse(allsubj %in% mutsubj, "mut", "wt"), allsubj)
    return(list(values = v, label = sprintf("%s status", arg), kind = "cat"))
  }
  NULL
}

# Attach a modality token as a column on `frame` (keyed by subject) and return the
# new column name; pass-through (modality=FALSE) when the token is not a modality.
.ctdna_attach_modality <- function(prep, frame, token, scheme = "three") {
  mr <- .ctdna_modality_resolve(prep, token, scheme)
  if (is.null(mr)) return(list(frame = frame, token = token, modality = FALSE))
  subj <- .o("subject")
  cn <- make.names(gsub("[^A-Za-z0-9]+", "_", token))
  idx <- as.character(frame[[subj]])
  frame[[cn]] <- if (identical(mr$kind, "cat")) factor(unname(mr$values[idx]))
                 else as.numeric(mr$values[idx])
  list(frame = frame, token = cn, modality = TRUE, label = mr$label, kind = mr$kind)
}
