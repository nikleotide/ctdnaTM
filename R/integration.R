# ============================================================================
# ctdnaTM — multi-modal integration
# ============================================================================

# v0.61.0: ONE visit filter shared by ctdna_oncoprint and
# ctdna_alteration_grid, so the SAME `visit` resolves the SAME variant
# rows in both. Format-tolerant: a passed "Cycle 1 Day 1" matches a stored
# "C1D1" and vice versa (both sides normalised through .ctdna_normalize_visit).
.ctdna_apply_visit_filter <- function(df, visit, visit_col = "Visit_name",
                                      fn = "ctdna_plot", verbose = TRUE) {
  if (is.null(visit)) return(df)
  if (!visit_col %in% names(df))
    stop(fn, ": `visit` was set but column '", visit_col, "' is not in the data. ",
         "Set `visit_col=` to the right column or pass `visit = NULL`.", call. = FALSE)
  n0   <- nrow(df)
  keep <- .ctdna_normalize_visit(df[[visit_col]]) %in% .ctdna_normalize_visit(visit)
  df   <- df[!is.na(df[[visit_col]]) & keep, , drop = FALSE]
  if (!nrow(df))
    stop(fn, ": no rows left after filtering to visit = ",
         paste(sQuote(visit), collapse = ", "), ".", call. = FALSE)
  if (isTRUE(verbose) && nrow(df) < n0)
    message(sprintf("%s: visit filter kept %d of %d rows.", fn, nrow(df), n0))
  df
}
# Every cross-modality plot accepts a `stat` parameter with a recommended
# default. See each function's docs for the recommendation.

.ctdna_metric_col <- function(metric) {
  if (!metric %in% c("best_ratio","baseline_tf","last_tf"))
    stop("metric must be one of: best_ratio, baseline_tf, last_tf")
  metric
}

#' Merge per-subject ctDNA summary with another modality
#'
#' Convenience wrapper around [ctdna_summary()] + [merge()].
#'
#' @param ctdna_df Longitudinal ctDNA data.
#' @param other_df Modality data with a subject column.
#' @param by Join key (default `ctdna_opts("subject")`).
#' @return Merged data frame.
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 15, seed = 1)
#' prep <- ctdna_prepare(sim, verbose = FALSE)
#' # Attach a TMB summary to the ctDNA frame
#' merged <- ctdna_merge_modalities(prep$samples, sim$tmb)
#' head(merged)
#' @export
ctdna_merge_modalities <- function(ctdna_df, other_df, by = .o("subject")) {
  cs <- ctdna_summary(ctdna_df)
  merge(cs, other_df, by = by, all = FALSE)
}


# ---- ctDNA x RNA-seq --------------------------------------------------------

#' Plot ctDNA metric vs gene expression
#'
#' Filters a long-format expression frame to a single gene and correlates
#' with a per-subject ctDNA metric.
#'
#' @section Recommended statistic:
#' `"spearman"` (default): expression values are often non-normal
#' (especially without log/voom transformation) and ctDNA metrics span
#' orders of magnitude. Use `"pearson"` only if both axes are clearly
#' normal on the chosen scale.
#'
#' @param ctdna_df Longitudinal ctDNA data.
#' @param expr_df Long-format expression: subject + gene + expression.
#' @param gene Gene symbol.
#' @param metric ctDNA metric: `"best_ratio"` (default), `"baseline_tf"`,
#'   `"last_tf"`.
#' @param log_x Log-scale x-axis (default TRUE).
#' @param recist_color Color points by RECIST.
#' @param scheme RECIST scheme.
#' @param stat One of `"spearman"` (default), `"pearson"`, `"kendall"`,
#'   `"none"`. Recommended: `"spearman"`.
#' @param show_stats,title,subtitle,caption,xlab,ylab,point_size,point_alpha
#'   Display.
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 30, seed = 1)
#' prep <- ctdna_prepare(sim, verbose = FALSE)
#'
#' # ctDNA best-ratio vs PD-L1 expression
#' ctdna_plot_vs_expression(prep$samples, sim$expression, gene = "CD274")
#' @export
ctdna_plot_vs_expression <- function(ctdna_df, expr_df, gene,
                                     metric = .o("ctdna_metric"),
                                     log_x = TRUE,
                                     recist_color = TRUE,
                                     scheme = c("three","two","four","four_alt","five","six"),
                                     stat = c("spearman","pearson","kendall","none"),
                                     show_stats = .o("show_stats"),
                                     title = NULL, subtitle = NULL,
                                     caption = NULL,
                                     xlab = NULL, ylab = NULL,
                                     point_size = .o("point_size") + 1,
                                     point_alpha = .o("point_alpha"),
    stat_position = .o("stat_position"),
    legend_position = .o("legend_position"),
    xlim = NULL, ylim = NULL,
    y_scale = "linear") {
  .stats_label <- NULL

  scheme <- match.arg(scheme); stat <- match.arg(stat)
  metric <- .ctdna_metric_col(metric)
  gcol <- .o("gene"); ecol <- .o("expression")
  # Canonicalize expression unit (no-op if already canonical)
  expr_df <- .canonicalize_expression(expr_df)
  expr <- expr_df[expr_df[[gcol]] == gene, c(.o("subject"), ecol), drop = FALSE]
  if (nrow(expr) == 0) stop("Gene not found in expression frame: ", gene)
  d <- ctdna_merge_modalities(ctdna_df, expr)
  d$.x <- if (log_x) ctdna_floor_tf(d[[metric]]) else d[[metric]]
  d$.y <- d[[ecol]]
  d$.r <- ctdna_stratify_recist(d[[.o("recist")]], scheme)
  if (is.null(title)) title <- sprintf("ctDNA (%s) vs %s expression", metric, gene)
  # v0.42.0: scatter correlation -> caption. show_stats is strictly
  # TRUE/FALSE.
  stats_res <- NULL
  if (isTRUE(show_stats)) {
    d_corr <- d
    d_corr$.logx <- if (log_x) suppressWarnings(log10(d_corr$.x)) else d_corr$.x
    stats_res <- .compute_scatter_stats(d_corr, x_col = ".logx",
                                          y_col = ".y", method = "spearman")
  }
  if (is.null(xlab)) xlab <- if (log_x) sprintf("%s (log)", metric) else metric
  if (is.null(ylab)) {
    unit_lbl <- switch(.o("expression_unit"),
                       log2_tpm_plus_one = "log2(TPM+1)",
                       tpm               = "TPM",
                       .o("expression_unit"))
    ylab <- sprintf("%s expression [%s]", gene, unit_lbl)
  }

  p <- ggplot2::ggplot(d, ggplot2::aes(.data$.x, .data$.y))
  if (recist_color) {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = .data$.r),
                                 size = point_size, alpha = point_alpha) +
      ctdna_scale_recist() + ggplot2::labs(color = "RECIST")
  } else {
    p <- p + ggplot2::geom_point(size = point_size, alpha = point_alpha,
                                 color = "grey20")
  }
  p <- p + ggplot2::geom_smooth(method = "lm", se = TRUE,
                                color = "grey30", fill = "grey80",
                                linetype = "dashed", linewidth = 0.6,
                                formula = y ~ x, alpha = 0.4) +
    ctdna_theme() +
    .ctdna_resolve_y_scale(y_scale)
  if (log_x) p <- p + ggplot2::scale_x_log10()
  if (isTRUE(show_stats) && !is.null(stats_res))
    caption <- .compose_stats_caption(stats_res, "scatter",
                                        base_caption = caption)
  .ctdna_apply_limits(.finalize(p, title, subtitle, caption, xlab, ylab,
             legend_position = legend_position), xlim, ylim)
}


# ---- ctDNA x mutation status ------------------------------------------------

#' Plot ctDNA metric stratified by mutation status of a gene
#'
#' Builds a box-plot of a per-subject ctDNA metric (default: best ratio)
#' across groups defined by mutation status (or full alteration type)
#' of a specified gene in a long-format mutation frame. The mutation
#' frame is joined to the ctDNA frame by subject ID; subjects absent
#' from `mut_df` are treated as wild-type for that gene. When
#' `show_stats = TRUE`, hand-rolled comparison brackets are drawn
#' above the boxes (Wilcoxon for two groups, Kruskal-Wallis overall
#' + pairwise Wilcoxon for 3+), and the overall test type + N per
#' group go in the caption.
#'
#' @section Recommended statistic:
#' `"wilcox"` (default): rank-sum between mutation carriers and
#' non-carriers — handles skewed ctDNA values. Use `"t"` only if the
#' ctDNA metric is normal on the working scale. Use `"kruskal"` when
#' grouping by `alteration_type` (multi-state).
#'
#' @param ctdna_df Longitudinal ctDNA data.
#' @param mut_df Long-format mutations: subject + gene + mutation_status
#'   (+ optional panel, alteration_type).
#' @param gene Gene symbol.
#' @param panel Optional: restrict to a panel (e.g. `"74genes"`, `"500genes"`).
#' @param by Stratify by `"mutation_status"` (default; binary wt vs mut)
#'   or `"alteration_type"` (snv / cnv_amp / cnv_loss / lgr / fusion / wt).
#'   Resolved via `ctdna_opts(mutation)` / `ctdna_opts(alteration)`.
#' @param metric ctDNA metric.
#' @param stat One of `"wilcox"` (default), `"t"`, `"kruskal"`, `"none"`.
#'   Recommended: `"wilcox"` (two groups) or `"kruskal"` (multi-state).
#' @param show_stats,title,subtitle,caption,xlab,ylab,point_size,point_alpha
#'   Display.
#' @return A ggplot.
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 30, seed = 1)
#' prep <- ctdna_prepare(sim, verbose = FALSE)
#'
#' # ctDNA best-ratio in TP53 mutants vs wild-type
#' ctdna_plot_vs_mutation(prep$samples, sim$mutations, gene = "TP53")
#'
#' # Drill down by alteration type (uses Kruskal-Wallis automatically)
#' ctdna_plot_vs_mutation(prep$samples, sim$mutations, gene = "TP53",
#'                         by = "alteration_type")
#' @export
ctdna_plot_vs_mutation <- function(ctdna_df, mut_df, gene,
                                   panel = NULL,
                                   by = c("mutation_status", "alteration_type"),
                                   metric = .o("ctdna_metric"),
                                   recist_facet = FALSE,
                                   scheme = c("three","two","four","four_alt","five","six"),
                                   scales = c("fixed","free_y","free_x","free"),
                                   filter_scheme = NULL,
                                   stat = c("wilcox","t","kruskal","none"),
                                   show_stats = .o("show_stats"),
                                   title = NULL, subtitle = NULL,
                                   caption = NULL,
                                   xlab = "Mutation status",
                                   ylab = NULL,
                                   point_size = .o("point_size"),
                                   point_alpha = .o("point_alpha"),
    stat_position = .o("stat_position"),
    legend_position = .o("legend_position"),
    xlim = NULL, ylim = NULL,
    y_scale = "log10") {
  .stats_label <- NULL

  stat <- match.arg(stat)
  by   <- match.arg(by)
  scheme <- match.arg(scheme)
  scales <- match.arg(scales)
  metric <- .ctdna_metric_col(metric)
  # v1.0.0: filtering is a separate user step. If filter_scheme is
  # non-NULL, apply it for backward compat; otherwise leave mut_df alone.
  if (!is.null(filter_scheme))
    mut_df <- .filter_apply_df(mut_df, filter_scheme = filter_scheme)
  gcol <- .o("gene"); mcol <- .o("mutation"); pcol <- .o("panel")
  acol <- .o("alteration")
  group_col <- if (by == "mutation_status") mcol else acol
  mut <- mut_df[mut_df[[gcol]] == gene, , drop = FALSE]
  if (!is.null(panel) && pcol %in% names(mut))
    mut <- mut[mut[[pcol]] == panel, , drop = FALSE]
  if (nrow(mut) == 0) stop("Gene/panel not found in mutation frame.")
  if (!group_col %in% names(mut))
    stop(sprintf("Column '%s' (selected via by='%s') not found in mutation data.",
                 group_col, by))
  mut <- mut[, c(.o("subject"), group_col), drop = FALSE]
  d <- ctdna_merge_modalities(ctdna_df, mut)
  d$.y <- ctdna_floor_tf(d[[metric]])
  d[[group_col]] <- .ctdna_factor_canonical(d[[group_col]], kind = group_col)
  if (recist_facet && .o("recist") %in% names(d)) {
    d$.r <- ctdna_stratify_recist(d[[.o("recist")]], scheme)
  }

  if (is.null(title))
    title <- sprintf("ctDNA (%s) by %s %s%s", metric, gene, by,
                     if (!is.null(panel)) sprintf(" [%s]", panel) else "")
  # v0.42.0: bracket renderer for between-group comparison. show_stats
  # is strictly TRUE/FALSE.
  stats_res <- NULL
  if (isTRUE(show_stats))
    stats_res <- .compute_box_stats(d, group_col = group_col, value_col = ".y")
  if (is.null(ylab)) ylab <- sprintf("%s (log)", metric)
  if (identical(xlab, "Mutation status") && by == "alteration_type")
    xlab <- "Alteration type"

  p <- ggplot2::ggplot(d, ggplot2::aes(.data[[group_col]], .data$.y,
                                       color = .data[[group_col]],
                                       fill  = .data[[group_col]])) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.3,
                          linewidth = 0.7, fatten = 2.2,
                          width = .o("box_width")) +
    ggplot2::geom_jitter(width = 0.15, size = point_size, alpha = point_alpha,
                         shape = 21, stroke = 0.3, color = "grey20") +
    ctdna_scale_dose(aesthetics = c("color", "fill")) +
    .log_clearance(y_scale) + ctdna_theme() +
    NULL  # legend_position handled by .finalize
  if (recist_facet && ".r" %in% names(d)) {
    p <- p + ggplot2::facet_wrap(~ .data$.r, scales = scales) +
      ggplot2::labs(subtitle = if (is.null(subtitle))
                                  sprintf("Faceted by RECIST (scheme=%s)", scheme)
                                else subtitle)
  }
  counts <- as.data.frame(table(d[[group_col]]))
  names(counts) <- c("g","n"); counts$label <- paste0("n=", counts$n)
  # Put n into the x-tick label to avoid overlap with data
  new_labels <- setNames(
    paste0(as.character(counts$g), "\n(n=", counts$n, ")"),
    as.character(counts$g))
  p <- p + ggplot2::scale_x_discrete(labels = new_labels) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(8, 8, 14, 8),
      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 4)))
  if (isTRUE(show_stats) && !is.null(stats_res)) {
    p <- .add_brackets(p, stats_res, d,
                        group_col = group_col, value_col = ".y",
                        y_scale = y_scale)
    caption <- .compose_stats_caption(stats_res, "box", base_caption = caption)
  }
  .ctdna_apply_limits(.finalize(p, title, subtitle, caption, xlab, ylab,
             legend_position = legend_position), xlim, ylim)
}


# ---- ctDNA x TMB ------------------------------------------------------------

#' Plot ctDNA metric vs TMB from a chosen panel
#'
#' Scatter of a per-subject ctDNA metric against TMB from a long-format
#' TMB frame. ctDNA is plotted on a log10 x-axis (TMB on linear y by
#' default; `y_scale` accepts the usual transforms). A smoothed
#' linear-model fit with confidence band is overlaid. When
#' `show_stats = TRUE`, the chosen correlation (Spearman by default)
#' is computed on `log10(ctDNA)` vs `TMB` and reported in the caption
#' as `Spearman ρ = ...; p = ...; N = ...`.
#'
#' @section Recommended statistic:
#' `"spearman"` (default): TMB is right-skewed and ctDNA metrics span
#' orders of magnitude. Pearson is rarely appropriate.
#'
#' @param ctdna_df Longitudinal ctDNA data.
#' @param tmb_df Data with subject + TMB column (+ optional panel).
#' @param panel Optional panel filter.
#' @param metric ctDNA metric.
#' @param tmb_col Column name in `tmb_df` (default `ctdna_opts("tmb")`).
#' @param stat One of `"spearman"` (default), `"pearson"`, `"kendall"`,
#'   `"none"`. Recommended: `"spearman"`.
#' @param show_stats,title,subtitle,caption,xlab,ylab,point_size,point_alpha
#'   Display.
#' @return A ggplot.
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 30, seed = 1)
#' prep <- ctdna_prepare(sim, verbose = FALSE)
#'
#' # Default scatter with Spearman correlation in the caption
#' ctdna_plot_vs_tmb(prep$samples, sim$tmb)
#'
#' # Color by RECIST scheme
#' ctdna_plot_vs_tmb(prep$samples, sim$tmb,
#'                    recist_color = TRUE, scheme = "three")
#' @export
ctdna_plot_vs_tmb <- function(ctdna_df, tmb_df,
                              panel = NULL,
                              metric = .o("ctdna_metric"),
                              tmb_col = .o("tmb"),
                              recist_color = FALSE,
                              scheme = c("three","two","four","four_alt","five","six"),
                              stat = c("spearman","pearson","kendall","none"),
                              show_stats = .o("show_stats"),
                              title = NULL, subtitle = NULL,
                              caption = NULL,
                              xlab = NULL, ylab = NULL,
                              point_size = .o("point_size") + 1,
                              point_alpha = .o("point_alpha"),
    stat_position = .o("stat_position"),
    legend_position = .o("legend_position"),
    xlim = NULL, ylim = NULL,
    y_scale = "linear") {
  .stats_label <- NULL

  stat <- match.arg(stat)
  scheme <- match.arg(scheme)
  metric <- .ctdna_metric_col(metric)
  pcol <- .o("panel")
  tmb <- tmb_df
  if (!is.null(panel) && pcol %in% names(tmb))
    tmb <- tmb[tmb[[pcol]] == panel, , drop = FALSE]
  tmb <- tmb[, c(.o("subject"), tmb_col), drop = FALSE]
  d <- ctdna_merge_modalities(ctdna_df, tmb)
  d$.x <- ctdna_floor_tf(d[[metric]])
  d$.y <- d[[tmb_col]]
  if (recist_color && .o("recist") %in% names(d)) {
    d$.r <- ctdna_stratify_recist(d[[.o("recist")]], scheme)
  }
  if (is.null(title))
    title <- sprintf("ctDNA (%s) vs TMB%s", metric,
                     if (!is.null(panel)) sprintf(" [%s]", panel) else "")
  # v0.42.0: scatter correlation -> caption. No brackets for scatter.
  # show_stats is strictly TRUE/FALSE.
  stats_res <- NULL
  if (isTRUE(show_stats)) {
    d_corr <- d
    d_corr$.logx <- suppressWarnings(log10(d_corr$.x))
    stats_res <- .compute_scatter_stats(d_corr, x_col = ".logx",
                                          y_col = ".y", method = "spearman")
  }
  if (is.null(xlab)) xlab <- sprintf("%s (log)", metric)
  if (is.null(ylab)) ylab <- "TMB (mut/Mb)"

  p <- if (recist_color && ".r" %in% names(d)) {
    ggplot2::ggplot(d, ggplot2::aes(.data$.x, .data$.y, color = .data$.r)) +
      ggplot2::geom_point(size = point_size, alpha = point_alpha) +
      ggplot2::geom_smooth(method = "lm", se = TRUE, inherit.aes = FALSE,
                           ggplot2::aes(.data$.x, .data$.y),
                           color = "grey30", fill = "grey80",
                           linetype = "dashed", linewidth = 0.6,
                           formula = y ~ x, alpha = 0.4) +
      ctdna_scale_recist() +
      ggplot2::labs(color = "RECIST")
  } else {
    ggplot2::ggplot(d, ggplot2::aes(.data$.x, .data$.y)) +
      ggplot2::geom_point(size = point_size, alpha = point_alpha,
                          color = "grey20") +
      ggplot2::geom_smooth(method = "lm", se = TRUE,
                           color = "grey30", fill = "grey80",
                           linetype = "dashed", linewidth = 0.6,
                           formula = y ~ x, alpha = 0.4)
  }
  p <- p + ggplot2::scale_x_log10() + .ctdna_resolve_y_scale(y_scale) + ctdna_theme()
  if (isTRUE(show_stats) && !is.null(stats_res))
    caption <- .compose_stats_caption(stats_res, "scatter",
                                        base_caption = caption)
  .ctdna_apply_limits(.finalize(p, title, subtitle, caption, xlab, ylab,
             legend_position = legend_position), xlim, ylim)
}


# ---- ctDNA x IHC ------------------------------------------------------------

#' Plot ctDNA metric vs IHC marker score
#'
#' Scatter of a per-subject ctDNA metric against the IHC score for a
#' chosen marker from a long-format IHC frame. ctDNA is on a log10
#' x-axis (`y_scale` adjusts the y-axis); a smoothed linear-model fit
#' with confidence band is overlaid. When `show_stats = TRUE`, the
#' Spearman correlation is reported in the caption.
#'
#' @section Recommended statistic:
#' `"spearman"` (default): IHC scores are often bounded/skewed.
#'
#' @param ctdna_df Longitudinal ctDNA data.
#' @param ihc_df Long-format IHC: subject + marker + score.
#' @param marker Marker name.
#' @param metric ctDNA metric.
#' @param stat One of `"spearman"` (default), `"pearson"`, `"kendall"`,
#'   `"none"`. Recommended: `"spearman"`.
#' @param show_stats,title,subtitle,caption,xlab,ylab,point_size,point_alpha
#'   Display.
#' @return A ggplot.
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 30, seed = 1)
#' prep <- ctdna_prepare(sim, verbose = FALSE)
#'
#' # ctDNA vs PD-L1 IHC score with Spearman correlation in the caption
#' ctdna_plot_vs_ihc(prep$samples, sim$ihc, marker = "PDL1")
#' @export
ctdna_plot_vs_ihc <- function(ctdna_df, ihc_df, marker,
                              metric = .o("ctdna_metric"),
                              recist_color = FALSE,
                              scheme = c("three","two","four","four_alt","five","six"),
                              stat = c("spearman","pearson","kendall","none"),
                              show_stats = .o("show_stats"),
                              title = NULL, subtitle = NULL,
                              caption = NULL,
                              xlab = NULL, ylab = NULL,
                              point_size = .o("point_size") + 1,
                              point_alpha = .o("point_alpha"),
    stat_position = .o("stat_position"),
    legend_position = .o("legend_position"),
    xlim = NULL, ylim = NULL,
    y_scale = "linear") {
  .stats_label <- NULL

  stat <- match.arg(stat)
  scheme <- match.arg(scheme)
  metric <- .ctdna_metric_col(metric)
  mcol <- .o("marker"); scol <- .o("ihc_score")
  ihc <- ihc_df[ihc_df[[mcol]] == marker, c(.o("subject"), scol), drop = FALSE]
  if (nrow(ihc) == 0) stop("Marker not found in IHC frame: ", marker)
  d <- ctdna_merge_modalities(ctdna_df, ihc)
  d$.x <- ctdna_floor_tf(d[[metric]])
  d$.y <- d[[scol]]
  if (recist_color && .o("recist") %in% names(d)) {
    d$.r <- ctdna_stratify_recist(d[[.o("recist")]], scheme)
  }
  if (is.null(title))
    title <- sprintf("ctDNA (%s) vs %s IHC score", metric, marker)
  # v0.42.0: scatter correlation -> caption. show_stats is strictly
  # TRUE/FALSE.
  stats_res <- NULL
  if (isTRUE(show_stats)) {
    d_corr <- d
    d_corr$.logx <- suppressWarnings(log10(d_corr$.x))
    stats_res <- .compute_scatter_stats(d_corr, x_col = ".logx",
                                          y_col = ".y", method = "spearman")
  }
  if (is.null(xlab)) xlab <- sprintf("%s (log)", metric)
  if (is.null(ylab)) ylab <- sprintf("%s IHC score", marker)

  p <- if (recist_color && ".r" %in% names(d)) {
    ggplot2::ggplot(d, ggplot2::aes(.data$.x, .data$.y, color = .data$.r)) +
      ggplot2::geom_point(size = point_size, alpha = point_alpha) +
      ggplot2::geom_smooth(method = "lm", se = TRUE, inherit.aes = FALSE,
                            ggplot2::aes(.data$.x, .data$.y),
                            color = "grey30", fill = "grey80",
                            linetype = "dashed", linewidth = 0.6,
                            formula = y ~ x, alpha = 0.4) +
      ctdna_scale_recist() + ggplot2::labs(color = "RECIST")
  } else {
    ggplot2::ggplot(d, ggplot2::aes(.data$.x, .data$.y)) +
      ggplot2::geom_point(size = point_size, alpha = point_alpha,
                          color = "grey20") +
      ggplot2::geom_smooth(method = "lm", se = TRUE,
                           color = "grey30", fill = "grey80",
                           linetype = "dashed", linewidth = 0.6,
                           formula = y ~ x, alpha = 0.4)
  }
  p <- p + ggplot2::scale_x_log10() + .ctdna_resolve_y_scale(y_scale) + ctdna_theme()
  if (isTRUE(show_stats) && !is.null(stats_res))
    caption <- .compose_stats_caption(stats_res, "scatter",
                                        base_caption = caption)
  .ctdna_apply_limits(.finalize(p, title, subtitle, caption, xlab, ylab,
             legend_position = legend_position), xlim, ylim)
}


