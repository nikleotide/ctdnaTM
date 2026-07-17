# ============================================================================
# ctdnaTM — multi-modal integration
# ============================================================================
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
#' merged <- ctdna_merge_modalities(prep$ctdna, sim$tmb)
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
#' @param indications Optional character vector of Indication values (e.g.
#'   c("NSCLC","mCRPC")). When set, `ctdna_df` is restricted to patients whose
#'   Indication is in this set. Lookup falls through Indication -> indication
#'   -> Cancertype in the frame itself; if not found in the frame, a Patient_ID
#'   join against prep$clinical is attempted. NULL (default) = no filter.
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 30, seed = 1)
#' prep <- ctdna_prepare(sim, verbose = FALSE)
#'
#' # ctDNA best-ratio vs PD-L1 expression
#' ctdna_plot_vs_expression(prep$ctdna, sim$expression, gene = "CD274")
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
    y_scale = "linear",
                          indications = NULL) {
  # ---- cohort restriction ------------------------
  ctdna_df <- .ctdna_filter_by_indication(ctdna_df, indications)
  # -----------------------------------------------
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
  .finalize(p, title, subtitle, caption, xlab, ylab,
             legend_position = legend_position)
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
#' @param indications Optional character vector of Indication values (e.g.
#'   c("NSCLC","mCRPC")). When set, `ctdna_df` is restricted to patients whose
#'   Indication is in this set. Lookup falls through Indication -> indication
#'   -> Cancertype in the frame itself; if not found in the frame, a Patient_ID
#'   join against prep$clinical is attempted. NULL (default) = no filter.
#' @return A ggplot.
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 30, seed = 1)
#' prep <- ctdna_prepare(sim, verbose = FALSE)
#'
#' # ctDNA best-ratio in TP53 mutants vs wild-type
#' ctdna_plot_vs_mutation(prep$ctdna, sim$mutations, gene = "TP53")
#'
#' # Drill down by alteration type (uses Kruskal-Wallis automatically)
#' ctdna_plot_vs_mutation(prep$ctdna, sim$mutations, gene = "TP53",
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
    y_scale = "log10",
                          indications = NULL) {
  # ---- cohort restriction ------------------------
  ctdna_df <- .ctdna_filter_by_indication(ctdna_df, indications)
  # -----------------------------------------------
  .stats_label <- NULL

  stat <- match.arg(stat)
  by   <- match.arg(by)
  scheme <- match.arg(scheme)
  scales <- match.arg(scales)
  metric <- .ctdna_metric_col(metric)
  # v1.0.0: filtering is a separate user step. If filter_scheme is
  # non-NULL, apply it for backward compat; otherwise leave mut_df alone.
  if (!is.null(filter_scheme))
    mut_df <- ctdna_filter_apply(mut_df, filter_scheme = filter_scheme)
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
  .finalize(p, title, subtitle, caption, xlab, ylab,
             legend_position = legend_position)
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
#' @param indications Optional character vector of Indication values (e.g.
#'   c("NSCLC","mCRPC")). When set, `ctdna_df` is restricted to patients whose
#'   Indication is in this set. Lookup falls through Indication -> indication
#'   -> Cancertype in the frame itself; if not found in the frame, a Patient_ID
#'   join against prep$clinical is attempted. NULL (default) = no filter.
#' @return A ggplot.
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 30, seed = 1)
#' prep <- ctdna_prepare(sim, verbose = FALSE)
#'
#' # Default scatter with Spearman correlation in the caption
#' ctdna_plot_vs_tmb(prep$ctdna, sim$tmb)
#'
#' # Color by RECIST scheme
#' ctdna_plot_vs_tmb(prep$ctdna, sim$tmb,
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
    y_scale = "linear",
                          indications = NULL) {
  # ---- cohort restriction ------------------------
  ctdna_df <- .ctdna_filter_by_indication(ctdna_df, indications)
  # -----------------------------------------------
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
  .finalize(p, title, subtitle, caption, xlab, ylab,
             legend_position = legend_position)
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
#' @param indications Optional character vector of Indication values (e.g.
#'   c("NSCLC","mCRPC")). When set, `ctdna_df` is restricted to patients whose
#'   Indication is in this set. Lookup falls through Indication -> indication
#'   -> Cancertype in the frame itself; if not found in the frame, a Patient_ID
#'   join against prep$clinical is attempted. NULL (default) = no filter.
#' @return A ggplot.
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 30, seed = 1)
#' prep <- ctdna_prepare(sim, verbose = FALSE)
#'
#' # ctDNA vs PD-L1 IHC score with Spearman correlation in the caption
#' ctdna_plot_vs_ihc(prep$ctdna, sim$ihc, marker = "PDL1")
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
    y_scale = "linear",
                          indications = NULL) {
  # ---- cohort restriction ------------------------
  ctdna_df <- .ctdna_filter_by_indication(ctdna_df, indications)
  # -----------------------------------------------
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
  .finalize(p, title, subtitle, caption, xlab, ylab,
             legend_position = legend_position)
}



# ---- Combinatorial alteration grid ------------------------------------------
#
# Faceted stacked bar plot showing the rate of "any alteration" in
# pre-defined gene sets across patient subgroups (response category x
# indication / cancer type). Adapted from a common TM deliverable —
# columns are gene sets, rows are indications, x within each panel is
# the response category, fill is alt vs wt with text labels showing
# percentages and counts.
#
# Statistics: per panel, Fisher's exact test on the (response x alt/wt)
# contingency table; reported as p-value above each panel.
# -----------------------------------------------------------------------------

#' Combinatorial alteration grid (indication x gene set)
#'
#' Stacked-bar grid showing the proportion of patients with any
#' alteration in each of several gene sets, broken down by response
#' category and faceted by indication.
#'
#' Mirrors the deck-style "Combinatorial Alteration Status" plot:
#' rows = indication (e.g. cancer type), columns = gene set, x-axis
#' within each panel = response category, fill = "alt" vs "wt", with
#' percentage and count labels in each bar and a Fisher's-exact p-value
#' annotated in the top-left of each panel.
#'
#' "Alt" is defined as: the patient has at least one detected alteration
#' in any gene of that set (after filtering). "Wt" is the complement.
#'
#' @param df Alteration data frame (vendor or canonical schema).
#' @param gene_sets Named list of gene-symbol vectors (one column per set).
#' @param patient_data Per-patient frame with at least the patient ID,
#'   response category, and indication columns. (v0.32.0+) Also accepts
#'   a NAMED list of CDISC ADaM datasets, e.g.
#'   \code{list(adsl = adsl, adrs = adrs, adtr = adtr)}; the package
#'   auto-builds the patient frame via
#'   \code{\link{ctdna_make_patient_data}}.
#' @param response_col Name of the response category column in
#'   \code{patient_data} (default \code{"RECIST"}).
#' @param indication_col Name of the indication / cancer-type column in
#'   \code{patient_data} (default \code{"Cancertype"}).
#' @param response_levels Optional level order for the x-axis (defaults
#'   to alphabetical order of values found in the data).
#' @param filter Optional named list of arguments forwarded to
#'   status. NULL = use \code{df} as given.
#' @param patient_col,gene_col Column overrides (auto-detected when NULL).
#' @param stat \code{"fisher"} (default), \code{"chisq"}, or \code{"none"}.
#' @param show_counts If TRUE, write "n / total" in the bar segments.
#' @param title,subtitle,caption Display options.
#' @param indications Optional character vector of Indication values (e.g.
#'   c("NSCLC","mCRPC")). When set, `df` is restricted to patients whose
#'   Indication is in this set. Lookup falls through Indication -> indication
#'   -> Cancertype in the frame itself; if not found in the frame, a Patient_ID
#'   join against prep$clinical is attempted. NULL (default) = no filter.
#' @return list(plot, summary, stats) where \code{summary} is a long-
#'   format per-cell summary and \code{stats} is one p-value per panel.
#' @examples
#' sim <- ctdna_make_mock_study(n_patients = 100, seed = 1)
#' ctdna_plot_alteration_grid(
#'   sim$infinity_report,
#'   gene_sets = list(
#'     KRAS  = "KRAS",
#'     AV    = c("ALK","ROS1","NTRK1","BRAF","RET"),
#'     Hansh = c("BRCA1","BRCA2","ATM","CHEK2","PALB2"),
#'     TSG   = c("TP53","RB1","MYC"),
#'     TOP1  = c("TOP1","TOP2A","ERCC1")
#'   ),
#'   patient_data   = sim$clinical,
#'   response_col   = "RECIST",
#'   indication_col = "Cancertype")$plot
#' @export
ctdna_plot_alteration_grid <- function(df,
                                        gene_sets,
                                        patient_data,
                                        response_col    = "RECIST",
                                        indication_col  = "Cancertype",
                                        response_levels = NULL,
                                        scheme          = c("raw","two","three","four"),
                                        scales          = c("fixed","free_y","free_x","free"),
                                        filter_scheme   = NULL,
                                        patient_col     = NULL,
                                        gene_col        = NULL,
                                        stat            = c("fisher","chisq","none"),
                                        show_counts     = TRUE,
                                        title    = "Combinatorial alteration status",
                                        subtitle = NULL,
                                        caption  = NULL,
                                        y_scale  = "linear",
                          indications = NULL) {
  # ---- cohort restriction ------------------------
  df <- .ctdna_filter_by_indication(df, indications)
  # -----------------------------------------------
  # y_scale is accepted for API consistency, but ignored — the y axis here
  # is a per-gene-set alteration frequency (0-1 stack) with fixed limits;
  # a log transform doesn't make sense and would break the c(0, 1.18) cap.
  stat <- match.arg(stat)
  scheme <- match.arg(scheme)
  scales <- match.arg(scales)

  if (!is.list(gene_sets) || is.null(names(gene_sets)))
    stop("`gene_sets` must be a NAMED list of gene-symbol vectors.")

  # v0.32.0: accept ADaM list as patient_data and auto-build
  patient_data <- .resolve_patient_data(patient_data)

  if (!is.data.frame(patient_data))
    stop("`patient_data` must be a data frame with patient ID + ",
         "response_col + indication_col, or a list of ADaM datasets ",
         "(adsl/adrs/adtr).")

  # Apply scheme collapse to the response column if requested
  if (scheme != "raw" && response_col %in% names(patient_data)) {
    patient_data[[response_col]] <- ctdna_stratify_recist(
      patient_data[[response_col]], scheme)
  }

  # v1.0.0: filtering is a separate user step. If filter_scheme is
  # non-NULL, apply it for backward compat; otherwise leave df alone.
  if (!is.null(filter_scheme))
    df <- ctdna_filter_apply(df, filter_scheme = filter_scheme,
                                      gene_col = gene_col)

  if (is.null(patient_col))
    patient_col <- .alt_col(df, c("Patient_ID", .o("subject"), "subject_id"))
  if (is.null(gene_col))
    gene_col    <- .alt_col(df, c("Gene", .o("gene"), "gene"))
  pid_p <- intersect(c(patient_col, "Patient_ID", "subject_id"),
                      names(patient_data))[1]
  if (is.na(pid_p))
    stop("Could not find a patient ID column in patient_data.")

  # All patients (from patient_data, so wt patients without alteration
  # rows still show up)
  all_pat <- patient_data[[pid_p]]

  # For each (patient, gene_set), is there any alteration?
  rows <- list()
  for (set_name in names(gene_sets)) {
    set_genes <- as.character(gene_sets[[set_name]])
    has_alt_pats <- unique(df[[patient_col]][
      !is.na(df[[gene_col]]) & df[[gene_col]] %in% set_genes])
    rows[[set_name]] <- data.frame(
      Patient_ID = all_pat,
      gene_set   = set_name,
      value      = ifelse(all_pat %in% has_alt_pats, "alt", "wt"),
      stringsAsFactors = FALSE)
  }
  long <- do.call(rbind, rows)
  names(long)[1] <- pid_p
  long$gene_set <- factor(long$gene_set, levels = names(gene_sets))

  # Join clinical
  pd <- patient_data[, c(pid_p, response_col, indication_col), drop = FALSE]
  long <- merge(long, pd, by = pid_p)

  if (is.null(response_levels)) {
    # Use canonical RECIST ordering when applicable; fall back to sort
    rl_factor <- .ctdna_factor_canonical(long[[response_col]],
                                          kind = if (response_col == "RECIST")
                                                    "RECIST" else NULL)
    response_levels <- levels(rl_factor)
  }
  long[[response_col]]   <- factor(long[[response_col]],   levels = response_levels)
  long[[indication_col]] <- .ctdna_factor_canonical(long[[indication_col]])
  long$value <- factor(long$value, levels = c("alt","wt"))

  # Per-cell summary
  summary <- aggregate(
    stats::as.formula(paste0("rep(1, nrow(long)) ~ ", indication_col, " + ",
                              response_col, " + gene_set + value")),
    data = transform(long, .one = 1),
    FUN  = length)
  names(summary)[5] <- "n"
  totals <- aggregate(
    stats::as.formula(paste0("n ~ ", indication_col, " + ",
                              response_col, " + gene_set")),
    data = summary, FUN = sum)
  names(totals)[4] <- "total"
  summary <- merge(summary, totals, by = c(indication_col, response_col, "gene_set"))
  summary$pct <- summary$n / summary$total

  # Per-panel stats
  stats_df <- NULL
  if (stat != "none") {
    keys <- unique(long[, c(indication_col, "gene_set")])
    stats_rows <- list()
    for (i in seq_len(nrow(keys))) {
      ind <- keys[i, indication_col]; gs <- keys[i, "gene_set"]
      sub <- long[long[[indication_col]] == ind & long$gene_set == gs, ]
      tab <- table(sub[[response_col]], sub$value)
      if (nrow(tab) < 2 || ncol(tab) < 2) {
        p <- NA_real_
      } else {
        if (stat == "fisher") {
          p <- tryCatch(stats::fisher.test(tab, workspace = 2e6,
                                            simulate.p.value = TRUE,
                                            B = 5000)$p.value,
                         error = function(e) NA_real_)
        } else {
          p <- tryCatch(suppressWarnings(
                          stats::chisq.test(tab)$p.value),
                         error = function(e) NA_real_)
        }
      }
      stats_rows[[i]] <- data.frame(
        setNames(list(ind, gs, p),
                  c(indication_col, "gene_set", "p")),
        stringsAsFactors = FALSE)
    }
    stats_df <- do.call(rbind, stats_rows)
    stats_df$label <- ctdna_pval_label(stats_df$p)
  }

  # Plot
  p <- ggplot2::ggplot(summary,
                       ggplot2::aes(x = .data[[response_col]], y = .data$pct,
                                    fill = .data$value)) +
    ggplot2::geom_col(width = .o("bar_width"), position = "stack") +
    ggplot2::scale_fill_manual(
      values = c(alt = "#F2786F", wt = "#3CB6BD"),
      breaks = c("alt","wt"),
      labels = c("alt","wt"),
      name   = "value") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(),
                                 limits = c(0, 1.18),
                                 expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::facet_grid(stats::as.formula(paste0(indication_col,
                                                  " ~ gene_set")),
                         scales = scales, space = "fixed") +
    ctdna_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1,
                                           size = 8),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.spacing = grid::unit(0.4, "lines")) +
    ggplot2::labs(x = "Response category",
                   y = "Percent genetic alteration",
                   title = title, subtitle = subtitle, caption = caption)

  if (show_counts) {
    # Label placement: use position_stack(vjust = 0.5) so each label
    # lands at the mid-height of its OWN segment, in whatever order
    # ggplot stacks the bars. Format: "10% (5)" — percent with count
    # in parentheses. Tiny segments (< 6%) are dropped to avoid text
    # crashes; can be tuned via the threshold below.
    summary$lbl <- sprintf("%.0f%% (%d)", 100 * summary$pct, summary$n)
    label_df <- summary[summary$pct >= 0.06, , drop = FALSE]
    if (nrow(label_df) > 0) {
      p <- p + ggplot2::geom_text(
        data = label_df,
        ggplot2::aes(x = .data[[response_col]],
                      y = .data$pct,
                      label = .data$lbl,
                      group = .data$value),
        position = ggplot2::position_stack(vjust = 0.5),
        size = 2.6, color = "grey10", fontface = "bold",
        inherit.aes = FALSE)
    }
  }

  if (!is.null(stats_df)) {
    # Place p-value at the TOP-LEFT corner inside each panel; clears the
    # alt-segment labels (which sit at any x-position within the panel
    # at y = ymid or y = 1.03 above the bar). The bar at x=1 is rare
    # for small labels in practice.
    p <- p + ggplot2::geom_text(
      data = stats_df,
      ggplot2::aes(x = -Inf, y = 1.15, label = .data$label),
      hjust = -0.08, vjust = 1, size = 2.4, color = "grey30",
      fontface = "italic",
      inherit.aes = FALSE)
  }

  list(plot = p, summary = summary, stats = stats_df)
}
