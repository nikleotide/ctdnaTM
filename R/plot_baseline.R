# ---- D1: Baseline characteristics -------------------------------------------

#' Baseline characteristics across one or two grouping factors
#'
#' Unified baseline boxplot (v0.44.0). Filters a prepared bundle to baseline
#' samples and shows TF panels (and optionally non-TF characteristics) split by
#' \code{group_by} (the x-axis) and, optionally, \code{subgroup_by} (drawn as
#' colour). Replaces the former \code{ctdna_plot_baseline} +
#' \code{ctdna_plot_baseline_dose_recist} pair.
#'
#' @param prep A \code{ctdna_prep} object (uses \code{prep$samples}) or a
#'   per-sample data frame.
#' @param group_by Primary split = x-axis. Opts key (e.g. \code{"dose"}) or a
#'   literal column name. Default \code{ctdna_opts("dose")}.
#' @param subgroup_by Optional second split, drawn as colour (e.g.
#'   \code{"RECIST"}). \code{NULL} (default) = single-factor overview.
#' @param compare_by Statistical comparison strategy:
#'   \code{"group"} (between \code{group_by} levels; the only mode that draws
#'   x-aligned brackets), \code{"subgroup"} (between \code{subgroup_by} levels),
#'   or \code{"subgroup_within_group"} (between subgroups within each group).
#'   The latter two report p-values in the caption. Requires
#'   \code{subgroup_by} unless \code{"group"}.
#' @param measure Which TF panel(s): \code{"all"}, \code{"methyl"},
#'   \code{"maxvaf"}, \code{"meanvaf"}, or \code{"custom"} (+\code{value_col}).
#' @param specific_columns Optional explicit column names to plot as panels
#'   (overrides \code{measure}); use for non-TF columns (cfDNA_ng, n.somatic,
#'   TMB_score).
#' @param value_col Column for \code{measure = "custom"}.
#' @param scheme RECIST stratification scheme.
#' @param scales facet scales. @param stat test family.
#' @param show_stats,show_n Toggles. @param title,subtitle,caption,xlab,ylab Labels.
#' @param point_size,point_alpha,stat_position,legend_position Display.
#' @param facet Optional extra faceting spec. @param y_scale y transform.
#' @param ylim Optional y limits. @param baseline_visit Baseline filter
#'   (\code{"auto"} default; literal visit; or \code{NULL} to skip).
#' @return A ggplot object.
#' @examples
#' prep <- ctdna_prepare(infinity_report = ctdna_make_mock_study(seed = 1)$infinity_report,
#'                       verbose = FALSE)
#' ctdna_plot_baseline(prep)                                       # overview by Dose
#' ctdna_plot_baseline(prep, group_by = "Dose", subgroup_by = "RECIST",
#'                     compare_by = "group", measure = "methyl")    # dose x RECIST
#' @export
ctdna_plot_baseline <- function(prep,
                          group_by    = .o("dose"),
                          subgroup_by = NULL,
                          compare_by  = c("group", "subgroup", "subgroup_within_group", "group_within_subgroup"),
                          pairwise = TRUE,
                          p_adjust = FALSE,
                          measure = c("all","methyl","maxvaf","meanvaf","custom"),
                          specific_columns = NULL,
                          value_col = NULL,
                          scheme = c("four","three","two","four_alt","five","six"),
                          scales = c("free_y","fixed","free_x","free"),
                          stat = c("auto","kruskal","wilcox","anova","t","none"),
                          show_stats = .o("show_stats"),
                          show_n = .o("show_n"),
                          drop_na = .o("drop_na"),
                          title = "Baseline characteristics",
                          subtitle = NULL, caption = NULL,
                          xlab = NULL, ylab = NULL,
                          point_size = .o("point_size"),
                          point_alpha = .o("point_alpha"),
                          stat_position = .o("stat_position"),
                          legend_position = .o("legend_position"),
                          facet = NULL,
                          y_scale = "log10",
                          ylim = NULL,
                          xlim = NULL,
                          baseline_visit = "auto") {

  compare_by <- match.arg(compare_by)
  if (identical(measure, c("all","methyl","maxvaf","meanvaf","custom"))) measure <- "all"
  stat <- match.arg(stat); scheme <- match.arg(scheme); scales <- match.arg(scales)
  # v0.44.5: tolerate y_scale being given a facet-scales keyword (a common
  # mix-up with `scales`): treat it as "fix/free the panels, linear transform".
  if (length(y_scale) == 1L && y_scale %in% c("fixed","free","free_y","free_x")) {
    scales  <- y_scale
    y_scale <- "linear"
  }

  # accept a whole prep or a per-sample frame; pull the per-sample grain
  df <- .ctdna_grain(prep, "samples")

  # filter to baseline rows
  df <- .ctdna_detect_baseline(df, baseline_visit = baseline_visit,
                               visit_col = .o("visit"))

  # resolve grouping columns (opts key or literal column name)
  group_by <- .resolve_col(group_by, df, what = "group_by")
  if (!is.null(subgroup_by))
    subgroup_by <- .resolve_col(subgroup_by, df, what = "subgroup_by")
  if (compare_by != "group" && is.null(subgroup_by))
    stop("ctdna_plot_baseline: compare_by = '", compare_by,
         "' requires subgroup_by to be set.", call. = FALSE)

  # v0.44.0: drop rows with NA in the grouping factor(s) — an NA dose/RECIST
  # group carries no meaning and only adds an empty box.
  if (isTRUE(drop_na)) {
    keep_rows <- !is.na(df[[group_by]])
    if (!is.null(subgroup_by)) keep_rows <- keep_rows & !is.na(df[[subgroup_by]])
    df <- df[keep_rows, , drop = FALSE]
  }

  # resolve TF panels
  if (length(measure) == 1L && measure == "custom") methods_resolved <- "custom"
  else methods_resolved <- .resolve_methods(measure, domain = "tf")
  tf_panels <- stats::setNames(
    lapply(methods_resolved, function(m)
      .resolve_method(m, "tf", value_col = if (m == "custom") value_col else NULL)),
    methods_resolved)
  tf_cols <- vapply(tf_panels, `[[`, character(1), "value_col")

  panel_labels <- stats::setNames(vapply(tf_panels, `[[`, character(1), "label"), tf_cols)
  panel_labels <- c(panel_labels,
    stats::setNames("Max VAF (%)",                .o("maxvaf")),
    stats::setNames("Mean VAF (%)",               .o("meanvaf")),
    stats::setNames("cfDNA concentration (ng/mL)",.o("cfdna_conc")),
    stats::setNames("# somatic alterations",      .o("n_somatic")),
    stats::setNames("TMB (mut/Mb)",               .o("tmb")))
  panel_labels <- panel_labels[!duplicated(names(panel_labels))]

  if (is.null(specific_columns)) {
    specific_columns <- tf_cols
    if (is.null(ylab))
      ylab <- if (length(methods_resolved) > 1L) "" else tf_panels[[1]]$label
  } else if (is.null(ylab)) {
    ylab <- if (length(specific_columns) == 1 &&
                specific_columns %in% names(panel_labels))
              panel_labels[[specific_columns]] else "Value (log scale)"
  }
  requested_columns <- specific_columns
  specific_columns  <- intersect(requested_columns, names(df))
  missing_cols      <- setdiff(requested_columns, specific_columns)
  if (length(missing_cols) && length(specific_columns))
    .ctdna_warn_once(paste0("ctdna_plot_baseline_missing_", paste(missing_cols, collapse="_")),
      sprintf("ctdna_plot_baseline: dropping column(s) not on samples: %s. Kept: %s.",
              paste(missing_cols, collapse=", "), paste(specific_columns, collapse=", ")))
  if (!length(specific_columns))
    stop("ctdna_plot_baseline: none of the requested column(s) are in samples.\n",
         "  Requested: ", paste(requested_columns, collapse=", "), "\n",
         "  samples columns: ", paste(names(df), collapse=", "), call. = FALSE)
  vars <- specific_columns

  # v0.44.0: TF columns are stored as PERCENT (0-100). Floor + "%" labels,
  # NO x100 scaling.
  tf_floor_vars <- unique(c(tf_cols, .o("tf"), .o("maxvaf"), .o("meanvaf"), .o("mr")))
  is_pct <- all(vars %in% tf_floor_vars)
  pct_labeller <- function(x) paste0(formatC(x, digits = 3, format = "g"), "%")
  scale_y_fn <- if (is_pct)
    function(...) .ctdna_resolve_y_scale(y_scale, labels = pct_labeller, ...)
  else function(...) .ctdna_resolve_y_scale(y_scale, ...)
  if (is.null(caption)) {
    fv <- .o("display_floor")
    caption <- if (is_pct) sprintf("Floor=%g%%; %s y-scale", fv, y_scale)
               else sprintf("%s y-scale", y_scale)
  }

  # RECIST stratification on whichever grouping column is RECIST
  recist_col <- .o("recist")
  if (identical(group_by, recist_col))
    df[[group_by]] <- ctdna_stratify_recist(df[[group_by]], scheme)
  if (!is.null(subgroup_by) && identical(subgroup_by, recist_col))
    df[[subgroup_by]] <- ctdna_stratify_recist(df[[subgroup_by]], scheme)

  # long over panels
  cols <- unique(c(group_by, subgroup_by, vars))
  long <- tidyr::pivot_longer(df[, cols, drop = FALSE],
                              cols = tidyr::all_of(vars),
                              names_to = "variable", values_to = "value")
  long <- long[!is.na(long$value), ]
  long$value <- ifelse(long$variable %in% tf_floor_vars,
                       ctdna_floor_tf(long$value), long$value)
  long[[group_by]] <- .ctdna_factor_canonical(long[[group_by]], kind = group_by)
  if (!is.null(subgroup_by))
    long[[subgroup_by]] <- .ctdna_factor_canonical(long[[subgroup_by]], kind = subgroup_by)

  # ---- stats ----
  stat_col <- switch(compare_by, group = group_by, subgroup = subgroup_by,
                     subgroup_within_group = subgroup_by,
                     group_within_subgroup = group_by)
  facet_stats <- NULL
  if (isTRUE(show_stats)) {
    facet_stats <- lapply(split(long, long$variable), function(s)
      .compute_box_stats(s, group_col = stat_col, value_col = "value", p_adjust = p_adjust))
    facet_stats <- facet_stats[!vapply(facet_stats, is.null, logical(1))]
  }

  # ---- plot ----
  if (is.null(subgroup_by)) {
    p <- ggplot2::ggplot(long, ggplot2::aes(.data[[group_by]], .data$value,
            color = .data[[group_by]], fill = .data[[group_by]])) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.3, linewidth = 0.7,
                            fatten = 2.2, width = .o("box_width")) +
      ggplot2::geom_jitter(width = 0.15, size = point_size, alpha = point_alpha,
                           shape = 21, stroke = 0.3, color = "grey20")
    scale_fn   <- if (identical(group_by, recist_col)) ctdna_scale_recist else ctdna_scale_dose
    legend_lab <- group_by
  } else {
    p <- ggplot2::ggplot(long, ggplot2::aes(.data[[group_by]], .data$value,
            color = .data[[subgroup_by]], fill = .data[[subgroup_by]])) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.3, linewidth = 0.7,
                            fatten = 2.2, width = .o("box_width"),
                            position = ggplot2::position_dodge(0.8, preserve = "single")) +
      ggplot2::geom_point(position = ggplot2::position_jitterdodge(
                            jitter.width = 0.15, dodge.width = 0.8),
                          size = point_size, alpha = point_alpha,
                          shape = 21, stroke = 0.3, color = "grey20")
    scale_fn   <- if (identical(subgroup_by, recist_col)) ctdna_scale_recist else ctdna_scale_dose
    legend_lab <- subgroup_by
  }

  p <- p +
    ggplot2::facet_wrap(~ variable, scales = scales,
                        labeller = ggplot2::labeller(variable = panel_labels)) +
    scale_y_fn() +
    scale_fn(aesthetics = c("color", "fill")) +
    ggplot2::labs(x = group_by, color = legend_lab, fill = legend_lab) +
    ctdna_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  # stats display, per compare_by
  if (isTRUE(show_stats) && length(facet_stats)) {
    if (compare_by == "group") {
      # x-aligned brackets per panel (comparison IS the x-axis grouping)
      for (var_name in names(facet_stats)) {
        sub <- long[long$variable == var_name, , drop = FALSE]
        p <- .add_brackets(p, facet_stats[[var_name]], sub,
                           group_col = group_by, value_col = "value",
                           y_scale = y_scale, facet_data = list(variable = var_name))
      }
      pvals_line <- paste(vapply(names(facet_stats), function(nm)
        sprintf("%s: p=%s", panel_labels[[nm]] %||% nm,
                .format_p_value(facet_stats[[nm]]$overall$p)), character(1)),
        collapse = " | ")
      caption <- sprintf("%s\nCompare by %s — overall p: %s", caption, group_by, pvals_line)

    } else if (compare_by == "subgroup_within_group" && !is.null(subgroup_by)) {
      for (var_name in names(facet_stats)) {
        sub <- long[long$variable == var_name, , drop = FALSE]
        p <- .add_dodge_brackets(p, sub, group_col = group_by, subgroup_col = subgroup_by,
                                 value_col = "value", mode = "within_group", pairwise = pairwise, p_adjust = p_adjust,
                                 facet_data = list(variable = var_name), y_scale = y_scale)
      }
      caption <- sprintf("%s\nCompare %s within each %s%s", caption, subgroup_by, group_by,
                         if (pairwise) " (pairwise)" else " (overall)")

    } else if (compare_by == "group_within_subgroup" && !is.null(subgroup_by)) {
      for (var_name in names(facet_stats)) {
        sub <- long[long$variable == var_name, , drop = FALSE]
        p <- .add_dodge_brackets(p, sub, group_col = group_by, subgroup_col = subgroup_by,
                                 value_col = "value", mode = "within_subgroup", pairwise = pairwise, p_adjust = p_adjust,
                                 facet_data = list(variable = var_name), y_scale = y_scale)
      }
      caption <- sprintf("%s\nCompare %s within each %s%s", caption, group_by, subgroup_by,
                         if (pairwise) " (pairwise)" else " (overall)")

    } else {  # subgroup: one pooled p per panel, centered at top
      ng <- nlevels(droplevels(long[[group_by]]))
      labdf <- do.call(rbind, lapply(names(facet_stats), function(nm)
        data.frame(variable = nm, .cx = (ng + 1) / 2,
                   lab = sprintf("%s p=%s", subgroup_by, .format_p_value(facet_stats[[nm]]$overall$p)),
                   stringsAsFactors = FALSE)))
      p <- p + ggplot2::geom_text(data = labdf,
        mapping = ggplot2::aes(x = .data$.cx, y = Inf, label = .data$lab),
        inherit.aes = FALSE, vjust = 1.5, size = 3.2, fontface = "bold", color = "grey15")
      caption <- sprintf("%s\nCompare by %s (pooled)", caption, subgroup_by)
    }
    if (.resolve_p_adjust(p_adjust) != "none")
      caption <- sprintf("%s; pairwise p %s-adjusted", caption, .resolve_p_adjust(p_adjust))
  }

  # n labels (per group_by on the x-axis + inside-panel white boxes)
  if (show_n) {
    if (group_by %in% names(df)) {
      df_nz <- df[!is.na(df[[group_by]]), , drop = FALSE]
      tb <- table(df_nz[[group_by]])
      n_overall <- stats::setNames(as.integer(tb), names(tb))
      p <- p + ggplot2::scale_x_discrete(labels = function(x) {
        n <- n_overall[as.character(x)]
        ifelse(is.na(n), as.character(x), sprintf("%s\n(n=%d)", as.character(x), n))
      })
    }
    pfc <- as.data.frame(table(long$variable, long[[group_by]]))
    names(pfc) <- c("variable", "g", "n")
    pfc <- pfc[pfc$n > 0, , drop = FALSE]
    pfc$g <- factor(pfc$g, levels = levels(long[[group_by]]))
    if (nrow(pfc) > 0L)
      p <- p + ggplot2::geom_label(data = pfc,
        mapping = ggplot2::aes(x = .data$g, y = -Inf, label = paste0("n=", .data$n)),
        inherit.aes = FALSE, vjust = -0.2, size = 2.8,
        label.padding = ggplot2::unit(0.12, "lines"),
        label.size = 0, color = "grey15", fill = "white", fontface = "bold")
  }

  cc_args <- list(clip = "off"); if (!is.null(ylim)) cc_args$ylim <- ylim
  if (!is.null(xlim)) cc_args$xlim <- xlim
  p <- p + do.call(ggplot2::coord_cartesian, cc_args)

  .finalize(p, title, subtitle, caption, xlab, ylab, legend_position = legend_position)
}
