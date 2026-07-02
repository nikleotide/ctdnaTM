# ============================================================================
# ctdnaTM — statistics helpers
# ============================================================================
# Each plotting function exposes a `stat` parameter that selects the test.
# This file provides the test catalog and small executors.

#' Format p-values for plot annotations
#'
#' Vectorized: accepts a scalar or numeric vector.
#'
#' @param p Numeric p-value(s).
#' @return Character vector of the same length, with entries like
#'   \code{"p=0.012"}, \code{"p<0.001"}, or \code{"p=NA"}.
#' @examples
#' ctdna_pval_label(0.0001)   # "<0.001"
#' ctdna_pval_label(0.034)    # "p = 0.034"
#' ctdna_pval_label(0.85)     # "p = 0.85"
#' @export
ctdna_pval_label <- function(p) {
  vapply(p, function(x) {
    if (is.na(x)) return("p=NA")
    if (x < 0.001) return("p<0.001")
    sprintf("p=%.3f", x)
  }, character(1))
}

# ---- catalogs ---------------------------------------------------------------
# Allowed `stat` values for group / correlation / concordance / paired plots.
.STAT_GROUP    <- c("auto","wilcox","t","kruskal","anova","none")
.STAT_TWO      <- c("wilcox","t","none")
.STAT_MULTI    <- c("kruskal","anova","wilcox","none")
.STAT_COR      <- c("spearman","pearson","kendall","none")
.STAT_PAIRED   <- c("paired_wilcox","paired_t","wilcox","t","none")
.STAT_CONCORD  <- c("both","kappa","agreement","mcnemar","none")


# ---- two-group + multi-group ------------------------------------------------

#' Run a single group-comparison test
#'
#' Used inside plot functions. `method = "auto"` selects Kruskal-Wallis
#' for >2 groups and Wilcoxon rank-sum for 2 groups.
#'
#' @param data Data frame.
#' @param value Numeric column name.
#' @param group Grouping column name.
#' @param method One of `"auto","wilcox","t","kruskal","anova","none"`.
#' @return list(method, name, p, n).
#' @examples
#' df <- data.frame(group = rep(c("a","b"), each = 10),
#'                   x     = c(rnorm(10), rnorm(10, 1)))
#' ctdna_run_group_test(df, "x", "group", method = "wilcox")
#' @export
ctdna_run_group_test <- function(data, value, group, method = "auto") {
  g <- data[[group]]; v <- data[[value]]
  groups <- unique(stats::na.omit(g))
  k <- length(groups); n <- sum(!is.na(v) & !is.na(g))
  if (method == "none" || k < 2)
    return(list(method = "none", name = "none", p = NA_real_, n = n))
  if (method == "auto") method <- if (k > 2) "kruskal" else "wilcox"
  p <- tryCatch({
    switch(method,
      wilcox  = stats::wilcox.test(v ~ g, exact = FALSE)$p.value,
      t       = stats::t.test(v ~ g)$p.value,
      kruskal = stats::kruskal.test(v ~ g)$p.value,
      anova   = stats::anova(stats::lm(v ~ g))$`Pr(>F)`[1])
  }, error = function(e) NA_real_)
  name <- switch(method,
    wilcox = "Wilcoxon", t = "t-test",
    kruskal = "Kruskal-Wallis", anova = "ANOVA")
  list(method = method, name = name, p = p, n = n)
}

#' All-pairs comparison test
#' @param data,value,group See [ctdna_run_group_test()].
#' @param method `"wilcox"` (default) or `"t"`.
#' @return Data frame: group1, group2, p, n1, n2.
#' @export
ctdna_pairwise_ranksum <- function(data, value, group, method = "wilcox") {
  g <- data[[group]]; v <- data[[value]]
  groups <- as.character(unique(stats::na.omit(g)))
  if (length(groups) < 2) return(NULL)
  pairs <- utils::combn(groups, 2, simplify = FALSE)
  do.call(rbind, lapply(pairs, function(p) {
    a <- v[g == p[1]]; b <- v[g == p[2]]
    na <- length(stats::na.omit(a)); nb <- length(stats::na.omit(b))
    if (na < 2 || nb < 2)
      return(data.frame(group1 = p[1], group2 = p[2], p = NA_real_,
                        n1 = na, n2 = nb))
    pv <- tryCatch({
      if (method == "wilcox") stats::wilcox.test(a, b, exact = FALSE)$p.value
      else stats::t.test(a, b)$p.value
    }, error = function(e) NA_real_)
    data.frame(group1 = p[1], group2 = p[2], p = pv, n1 = na, n2 = nb)
  }))
}


# ---- correlations -----------------------------------------------------------

#' Correlation test with rho, p, n
#'
#' @param x,y Numeric vectors.
#' @param method One of `"spearman"` (default), `"pearson"`, `"kendall"`.
#' @return list(method, name, rho, p, n).
#' @examples
#' set.seed(1)
#' x <- runif(30)
#' y <- 2 * x + rnorm(30, sd = 0.3)
#' ctdna_cor_test(x, y, method = "spearman")
#' @export
ctdna_cor_test <- function(x, y, method = "spearman") {
  ok <- !is.na(x) & !is.na(y); n <- sum(ok)
  if (method == "none" || n < 3)
    return(list(method = method, name = "none",
                rho = NA, p = NA, n = n))
  ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = method,
                                         exact = FALSE))
  list(method = method,
       name = switch(method, spearman = "Spearman",
                     pearson = "Pearson", kendall = "Kendall"),
       rho = unname(ct$estimate), p = ct$p.value, n = n)
}


# ---- paired tests (baseline vs each timepoint) ------------------------------

#' Paired test
#' @param x,y Numeric vectors of equal length (paired observations).
#' @param method One of `"paired_wilcox"`, `"paired_t"`, `"wilcox"`, `"t"`.
#' @return list(method, name, p, n).
#' @export
ctdna_paired_test <- function(x, y, method = "paired_wilcox") {
  ok <- !is.na(x) & !is.na(y); n <- sum(ok)
  if (method == "none" || n < 2)
    return(list(method = "none", name = "none", p = NA, n = n))
  paired <- method %in% c("paired_wilcox","paired_t")
  base <- if (paired) sub("paired_", "", method) else method
  p <- tryCatch({
    if (base == "wilcox")
      stats::wilcox.test(x[ok], y[ok], paired = paired, exact = FALSE)$p.value
    else
      stats::t.test(x[ok], y[ok], paired = paired)$p.value
  }, error = function(e) NA_real_)
  list(method = method,
       name = paste0(if (paired) "paired " else "",
                     if (base == "wilcox") "Wilcoxon" else "t-test"),
       p = p, n = n)
}


# ---- concordance ------------------------------------------------------------

#' Cohen's kappa from a contingency table
#' @param tab Square table or matrix.
#' @return Numeric kappa.
#' @export
ctdna_cohen_kappa <- function(tab) {
  tab <- as.matrix(tab)
  common <- intersect(rownames(tab), colnames(tab))
  tab <- tab[common, common, drop = FALSE]
  n <- sum(tab); if (n == 0) return(NA_real_)
  po <- sum(diag(tab)) / n
  pe <- sum(rowSums(tab) * colSums(tab)) / n^2
  if (pe == 1) return(NA_real_)
  (po - pe) / (1 - pe)
}

# v0.46.0: is a vector categorical/ordinal rather than continuous?
.ctdna_is_categorical <- function(v) {
  if (is.factor(v) || is.character(v) || is.logical(v)) return(TRUE)
  if (is.numeric(v)) {
    u <- unique(v[is.finite(v)])
    return(length(u) <= 6L && all(u == round(u)))   # small discrete integer set
  }
  FALSE
}

# Cramér's V (association effect size) from a contingency table.
.cramers_v <- function(tab) {
  tab <- as.matrix(tab); n <- sum(tab)
  if (n == 0 || nrow(tab) < 2L || ncol(tab) < 2L) return(NA_real_)
  chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE)$statistic)
  unname(sqrt((chi / n) / min(nrow(tab) - 1L, ncol(tab) - 1L)))
}

#' Weighted Cohen's kappa for ordinal agreement
#'
#' Linear or quadratic weighted kappa on a square contingency table whose
#' row and column categories match (ordered the same way).
#' @param tab A contingency table (matching row/col ordinal categories).
#' @param weights `"quadratic"` (default) or `"linear"`.
#' @return Weighted kappa, or `NA` if not computable.
#' @export
ctdna_weighted_kappa <- function(tab, weights = c("quadratic","linear")) {
  weights <- match.arg(weights)
  tab <- as.matrix(tab)
  common <- intersect(rownames(tab), colnames(tab))
  if (length(common) < 2L) return(NA_real_)
  tab <- tab[common, common, drop = FALSE]
  n <- sum(tab); if (n == 0) return(NA_real_)
  k <- nrow(tab); idx <- seq_len(k)
  w  <- outer(idx, idx, function(i, j)
    if (weights == "quadratic") ((i - j)^2) / ((k - 1)^2) else abs(i - j) / (k - 1))
  po <- sum((1 - w) * tab) / n
  e  <- outer(rowSums(tab), colSums(tab)) / n
  pe <- sum((1 - w) * e) / n
  if (pe == 1) return(NA_real_)
  (po - pe) / (1 - pe)
}

#' Concordance summary for a 2-way table
#' @param tab Contingency table.
#' @param method One of `"both","kappa","agreement","mcnemar","none"`.
#' @return list with named numeric components.
#' @export
ctdna_concordance_test <- function(tab, method = "both") {
  if (method == "none") return(list(method = "none"))
  agree <- {
    cm <- intersect(rownames(tab), colnames(tab))
    sum(diag(as.matrix(tab)[cm, cm, drop = FALSE])) / sum(tab)
  }
  k <- ctdna_cohen_kappa(tab)
  mcn <- if (method %in% c("mcnemar","both") &&
              nrow(tab) == 2 && ncol(tab) == 2) {
    tryCatch(stats::mcnemar.test(tab)$p.value, error = function(e) NA_real_)
  } else NA_real_
  list(method = method, agreement = agree, kappa = k,
       mcnemar_p = mcn, n = sum(tab))
}


# ---- internal: finalize titles and clearance markers ------------------------

.finalize <- function(p, title = NULL, subtitle = NULL, caption = NULL,
                     xlab = NULL, ylab = NULL,
                     legend_position = .o("legend_position")) {
  l <- list()
  if (!is.null(title))    l$title    <- title
  if (!is.null(subtitle)) l$subtitle <- subtitle
  if (!is.null(caption))  l$caption  <- caption
  if (!is.null(xlab))     l$x        <- xlab
  if (!is.null(ylab))     l$y        <- ylab
  if (length(l)) p <- p + do.call(ggplot2::labs, l)
  if (!is.null(legend_position))
    p <- p + ggplot2::theme(legend.position = legend_position)
  p
}

# Routes a stats label string to its display location.
#  position = "subtitle" : returns list(subtitle = label, caption = caption, on_plot = NULL)
#  position = "caption"  : returns list(subtitle = subtitle, caption = label_or_appended, on_plot = NULL)
#  position = "on_plot"  : returns list(subtitle = subtitle, caption = caption, on_plot = label)
#  position = "none"     : returns list(subtitle = subtitle, caption = caption, on_plot = NULL)
.route_stats <- function(label, position,
                          subtitle = NULL, caption = NULL,
                          wrap_width = 70) {
  position <- match.arg(position,
                         c("subtitle","caption","on_plot","none"))
  out <- list(subtitle = subtitle, caption = caption, on_plot = NULL)
  if (position == "none" || is.null(label) || !nzchar(label)) return(out)
  # Helper: wrap a long stats label across multiple lines using existing
  # separators (`|`, `;`, ` - `) as soft break points.
  wrap_label <- function(s, width) {
    if (nchar(s) <= width) return(s)
    # Replace " - " with a sentinel so it's treated as a break, then split
    s_norm <- gsub(" - ", "\u0001", s, fixed = TRUE)
    parts <- strsplit(s_norm, "\\s*(\\||;|\u0001)\\s*")[[1]]
    out_lines <- character(0); cur <- ""
    for (p in parts) {
      cand <- if (nzchar(cur)) paste(cur, p, sep = "; ") else p
      if (nchar(cand) > width && nzchar(cur)) {
        out_lines <- c(out_lines, cur); cur <- p
      } else cur <- cand
    }
    if (nzchar(cur)) out_lines <- c(out_lines, cur)
    paste(out_lines, collapse = "\n")
  }
  label <- wrap_label(label, wrap_width)
  if (position == "subtitle") {
    out$subtitle <- if (is.null(subtitle) || !nzchar(subtitle)) label
                    else paste(subtitle, label, sep = "\n")
  } else if (position == "caption") {
    out$caption  <- if (is.null(caption) || !nzchar(caption)) label
                    else paste(caption, label, sep = "\n")
  } else if (position == "on_plot") {
    out$on_plot  <- label
  }
  out
}

# Annotate a single overall-stat label inside the plot panel. Uses
# `annotation_custom` with grid::textGrob in NPC coordinates so the label
# is anchored to the panel even when scales are log / transformed and
# never falls outside the panel boundaries (which the older -Inf/Inf
# approach did with clip = 'on').
.annotate_stat <- function(p, label, x_npc = 0.02, y_npc = 0.97,
                            hjust = 0, vjust = 1) {
  if (is.null(label) || !nzchar(label)) return(p)
  grob <- grid::textGrob(
    label = label,
    x = grid::unit(x_npc, "npc"),
    y = grid::unit(y_npc, "npc"),
    hjust = hjust, vjust = vjust,
    gp = grid::gpar(col = "grey20", fontsize = 9))
  p + ggplot2::annotation_custom(grob = grob,
                                  xmin = -Inf, xmax = Inf,
                                  ymin = -Inf, ymax = Inf)
}

# Draw pairwise-comparison brackets above a discrete x-axis boxplot.
# pw_df: data frame with columns group1, group2, p (in plotting order)
# x_levels: factor levels of the x-axis variable (used for positions)
# y_top: log10 (if log=TRUE) or linear top-of-data y-value
# log: TRUE if y-axis is log10
# max_show: cap number of brackets to avoid stacking
.add_pairwise_brackets <- function(p, pw_df, x_levels, y_top,
                                    log = TRUE, max_show = 5) {
  if (is.null(pw_df) || nrow(pw_df) == 0) return(p)
  pw_df <- pw_df[!is.na(pw_df$p), ]
  pw_df <- pw_df[order(pw_df$p), ]
  if (nrow(pw_df) > max_show) pw_df <- utils::head(pw_df, max_show)
  if (nrow(pw_df) == 0) return(p)
  pw_df$x1 <- match(as.character(pw_df$group1), x_levels)
  pw_df$x2 <- match(as.character(pw_df$group2), x_levels)
  pw_df    <- pw_df[!is.na(pw_df$x1) & !is.na(pw_df$x2), ]
  if (nrow(pw_df) == 0) return(p)

  # Stack brackets vertically: bigger stack for more pairs
  step <- if (log) 0.12 else (y_top * 0.08)
  for (i in seq_len(nrow(pw_df))) {
    y_level <- if (log) 10 ^ (log10(y_top) + step * i)
               else y_top + step * i
    p <- p +
      ggplot2::annotate("segment",
                         x = pw_df$x1[i], xend = pw_df$x2[i],
                         y = y_level, yend = y_level,
                         color = "grey30", linewidth = 0.4) +
      ggplot2::annotate("text",
                         x = (pw_df$x1[i] + pw_df$x2[i]) / 2,
                         y = y_level,
                         label = ctdna_pval_label(pw_df$p[i]),
                         vjust = -0.3, size = 3, color = "grey20")
  }
  p
}

.log_clearance <- function(y_scale = "log10") {
  # Dotted horizontal line at the display floor. The label that used to
  # be drawn in-panel is now communicated via the plot caption (set by
  # each calling plot), avoiding overlap with data and with facet ticks.
  # y_scale (v0.25.0) lets the caller swap out the log10 transform.
  c(
    list(ggplot2::geom_hline(yintercept = .o("display_floor"), linetype = "dotted",
                        color = .o("clearance_color"), linewidth = 0.4)),
    .ctdna_resolve_y_scale(y_scale)
  )
}

.pairwise_str <- function(pw, max_show = 3) {
  if (is.null(pw) || nrow(pw) == 0) return("")
  pw <- pw[order(pw$p), ]
  shown <- utils::head(pw, max_show)
  paste(sprintf("%s vs %s: %s", shown$group1, shown$group2,
                vapply(shown$p, ctdna_pval_label, character(1))),
        collapse = "; ")
}


# ---- v0.22.0: canonical axis ordering ---------------------------------------
# Many categorical variables have a meaningful, non-alphabetic order:
#   dose:   Low < Mid < High
#   RECIST: CR < PR < uCR < uPR < SD < PD < NE < NA
#   visits: chronological (Baseline first, then sorted by cycle/day)
#   Functional_impact: High > Moderate > Low
#   TMB_category:      Low < Intermediate < High
#   ctDNA_detection_status: Detected > Below LOQ > Not Detected
#
# .ctdna_factor_canonical() returns the value as a factor with the
# semantic order applied (unrecognised levels appended at the end in
# their natural sort order). Used inside plot functions to set axis
# order without forcing each call site to know the canonical sequence.

.CANONICAL_LEVELS <- list(
  dose            = c("Low","Mid","High",
                       "low","mid","high",
                       "LOW","MID","HIGH"),
  # v0.41.1: include stratify_recist outputs from every recist_grouping
  # scheme so .ctdna_factor_canonical(kind="RECIST") preserves order
  # after pivot_longer demotes the factor to character. Order is response
  # quality, best → worst: CR/PR variants first, then uCR/uPR, then SD,
  # then PD/NE/NA variants, then the two-scheme "NR" tail.
  RECIST          = c("CR","PR","CR+PR","R",
                       "uCR","uPR",
                       "SD",
                       "PD","NE","NA","NE+NA","PD+NE+NA","NR"),
  RECIST_collapsed = c("CR/PR","CR/PR + uCR/uPR",
                        "uCR/uPR","SD","PD/NE/NA",
                        "Responders","NonResponders"),
  Functional_impact = c("High","Moderate","Low"),
  TMB_category    = c("Low","Intermediate","High"),
  ctDNA_detection_status = c("Detected","Below LOQ","Not Detected"),
  Mutation_status = c("alt","mut","wt"),
  alteration_type = c("SNV","Indel","CNV","Fusion","LGR")
)

# Visit ordering is data-dependent (Baseline first, then C<n>D<n> by
# cycle, then EOT last). Returns a character vector of levels in the
# right order given the unique values found in `x`.
.canonical_visit_order <- function(x) {
  u <- unique(as.character(stats::na.omit(x)))
  if (length(u) == 0) return(u)
  # Baseline always first
  bl <- grep("^(baseline|Baseline|BASELINE|C1D1)$", u, value = TRUE)
  # EOT-like always last
  eot <- grep("^(EOT|End of treatment|End_of_treatment)$", u, value = TRUE,
              ignore.case = TRUE)
  # Middle: extract numbers and sort
  middle <- setdiff(u, c(bl, eot))
  # Try to sort by extracted (cycle, day) numbers
  cyc <- suppressWarnings(as.integer(sub(".*C(\\d+)D.*", "\\1", middle)))
  day <- suppressWarnings(as.integer(sub(".*D(\\d+).*", "\\1", middle)))
  ord <- order(cyc, day, middle, na.last = TRUE)
  middle <- middle[ord]
  c(unique(bl), middle, unique(eot))
}

# v0.44.0: clean dose labels as comprehensively as safely possible.
# (1) If a value contains a "N mg" / "N mg/kg" token, reduce it to that token
#     (drops compound codes and schedule suffixes).
# (2) Otherwise (categorical doses with no mg pattern), remove the whitespace
#     tokens common to every label, keeping only the distinguishing part,
#     never reducing a label to empty.
.ctdna_strip_common_tokens <- function(x) {
  v <- trimws(as.character(x))
  uq <- unique(v[!is.na(v) & nzchar(v)])
  if (length(uq) < 2) return(v)
  toks <- strsplit(uq, "\\s+")
  common <- Reduce(intersect, toks)
  if (!length(common)) return(v)
  map <- vapply(uq, function(one) {
    keep <- setdiff(strsplit(one, "\\s+")[[1]], common)
    out  <- paste(keep, collapse = " ")
    if (!nzchar(out)) one else out          # never blank out
  }, character(1))
  unname(map[match(v, names(map))])
}

.ctdna_clean_dose <- function(x) {
  v   <- trimws(as.character(x))
  pat <- "[0-9]+(?:\\.[0-9]+)?\\s*mg(?:\\s*/\\s*kg)?"
  m   <- regexpr(pat, v, perl = TRUE, ignore.case = TRUE)
  hit <- !is.na(v) & m > 0
  cleaned <- v
  if (any(hit)) cleaned[hit] <- regmatches(v, m)   # keep only the "N mg(/kg)" token
  if (!any(hit)) cleaned <- .ctdna_strip_common_tokens(v)  # categorical: strip shared tokens
  cleaned <- gsub("\\s*/\\s*", "/", cleaned)   # "mg / kg" -> "mg/kg"
  cleaned <- gsub("\\s+", " ", cleaned)
  trimws(cleaned)
}

.ctdna_factor_canonical <- function(x, kind = NULL) {
  if (is.factor(x)) x <- as.character(x)
  # v0.44.0: clean dose labels as much as possible.
  if (!is.null(kind) && tolower(kind) == "dose") x <- .ctdna_clean_dose(x)
  u <- unique(as.character(stats::na.omit(x)))
  if (length(u) == 0) return(factor(x))
  # Choose canonical level set
  lev <- NULL
  if (!is.null(kind)) {
    if (identical(kind, "visit") || identical(kind, "time_point")) {
      lev <- .canonical_visit_order(x)
    } else if (kind %in% names(.CANONICAL_LEVELS)) {
      lev <- .CANONICAL_LEVELS[[kind]]
    }
  }
  if (is.null(lev)) {
    # Heuristic auto-detection from values
    for (k in names(.CANONICAL_LEVELS)) {
      cand <- .CANONICAL_LEVELS[[k]]
      if (all(u %in% cand)) { lev <- cand; break }
    }
    if (is.null(lev) && any(grepl("^C\\d+D\\d+", u))) {
      lev <- .canonical_visit_order(x)
    }
  }
  if (is.null(lev)) return(factor(x))   # fall back to alphabetical
  # Keep only levels actually present in x, in canonical order, then
  # append any leftover levels in their natural sort order.
  present_canonical <- intersect(lev, u)
  leftover <- setdiff(u, lev)
  if (length(leftover) > 0) leftover <- sort(leftover)
  factor(x, levels = c(present_canonical, leftover))
}

# v0.41.1: resolve a user-supplied identifier to an actual column name on
# `df`. Accepts either (a) a literal column name already present in `df`,
# or (b) an opts key whose value resolves to a column in `df`. Used by
# plot functions so users can pass `group = "dose"` (an opts key) OR
# `group = "Dose"` (the literal column) interchangeably — the case-shift
# from canonical-named to vendor-named defaults in v0.41.0 broke the
# former without a resolver in place.
#
# Errors with a clear message listing the columns actually on `df` when
# neither resolution path hits. `what` is included in the error to point
# the user at which argument failed (e.g. `what = "group"`).
.resolve_col <- function(x, df, what = "column") {
  if (is.null(x) || length(x) != 1L || !is.character(x) || !nzchar(x))
    stop(sprintf("`%s` must be a single non-empty string (got: %s).",
                  what, deparse(x)[1]), call. = FALSE)
  if (x %in% names(df)) return(x)
  # case-insensitive column match
  ci <- which(tolower(names(df)) == tolower(x))
  if (length(ci) == 1L) return(names(df)[ci])
  # opts key (exact or case-insensitive) whose value matches a column
  hit <- tryCatch(.o(x), error = function(e) NULL)
  if (is.null(hit)) {
    kk <- tryCatch(names(ctdna_opts()), error = function(e) NULL)
    if (!is.null(kk)) { m <- kk[tolower(kk) == tolower(x)]; if (length(m) == 1L) hit <- .o(m) }
  }
  if (!is.null(hit) && is.character(hit) && length(hit) == 1L) {
    if (hit %in% names(df)) return(hit)
    ci2 <- which(tolower(names(df)) == tolower(hit)); if (length(ci2) == 1L) return(names(df)[ci2])
  }
  stop(sprintf("`%s = \"%s\"` is neither a column of the data nor an opts key whose value matches a column.\n  data columns: %s",
                what, x, paste(names(df), collapse = ", ")),
       call. = FALSE)
}


# ============================================================
# v0.42.0: hand-rolled comparison brackets + stats helpers
# ============================================================
#
# Replaces the per-function ad-hoc stats code that lived inside each
# plot. All ctdnaTM plot functions (except oncoprint) now share these
# helpers so the rendering is consistent: brackets above the data
# stacked top-down by comparison width, p-values printed above each
# bracket, overall test type + N in the caption.
#
# Three entry points:
#
#   .format_p_value(p)           -> human-readable string ("<0.001", "0.034")
#   .compute_box_stats(df, ...)  -> overall + pairwise tests for box-plot family
#   .compute_scatter_stats(...)  -> correlation for scatter family
#   .add_brackets(p, ...)        -> appends bracket layers to a ggplot
#
# The bracket layout works in TRANSFORMED y-space (log10 for "log10"
# y_scale, linear otherwise). Brackets stack from lowest (narrowest
# comparison) to highest (widest), with a fixed fraction-of-data-span
# gap between them. The plot's y limits are expanded via the helper's
# returned `expand_factor` so the brackets never overflow.

.format_p_value <- function(p) {
  if (length(p) != 1L) return(vapply(p, .format_p_value, character(1)))
  if (is.na(p) || !is.finite(p))                return("NA")
  if (p < 0.001)                                return("<0.001")
  if (p < 0.01)                                 return(format(round(p, 3), nsmall = 3))
  if (p < 0.1)                                  return(format(round(p, 2), nsmall = 2))
  format(round(p, 2), nsmall = 2)
}

# Compute box-family stats: an overall test (Wilcoxon if 2 groups,
# Kruskal-Wallis if 3+) plus pairwise Wilcoxon tests for the bracket
# rendering. Pair selection: all pairs for K<=3; adjacent + extremes
# for larger K (keeps the panel readable).
#
# Returns list(overall, pairwise, n_per_group) or NULL if data too
# sparse to test.
#
# v0.42.3: relaxed the per-group n threshold from >=2 to >=1, and now
# accepts pairings where one group has n=1 as long as the other has
# >=2. wilcox.test(c(1), c(2,3,4)) is a valid test (low power, but
# the user gets a real p-value rather than a silent NULL). Previous
# behaviour dropped these pairings entirely, producing no brackets
# for plots like ctdna_plot_vs_mutation() on rare-mutation genes
# where only one or two carriers are present.
# v0.45.0: resolve a user-facing p_adjust value to a stats::p.adjust method.
# FALSE / NULL / "none" / absent -> "none"; TRUE -> "BH"; a method string is
# matched case-insensitively to the canonical name.
.resolve_p_adjust <- function(x) {
  if (is.null(x) || (is.logical(x) && length(x) && isFALSE(x))) return("none")
  if (is.logical(x) && length(x) && isTRUE(x))                  return("BH")
  s <- tolower(trimws(as.character(x)[1]))
  if (!nzchar(s) || s == "none") return("none")
  switch(s,
    "bh" = , "fdr"          = "BH",
    "by"                    = "BY",
    "holm"                  = "holm",
    "hochberg"              = "hochberg",
    "hommel"                = "hommel",
    "bonferroni" = , "bonf" = "bonferroni",
    { if (s %in% tolower(stats::p.adjust.methods))
        stats::p.adjust.methods[tolower(stats::p.adjust.methods) == s][1]
      else "none" })
}

.compute_box_stats <- function(df, group_col, value_col, test = "auto",
                               p_adjust = "none") {
  if (!group_col %in% names(df) || !value_col %in% names(df)) return(NULL)
  g <- df[[group_col]]; v <- df[[value_col]]
  ok <- !is.na(g) & !is.na(v) & is.finite(v)
  g <- g[ok]; v <- v[ok]
  if (!length(g)) return(NULL)

  # Preserve factor level order if present, else alphabetical
  group_levels <- if (is.factor(df[[group_col]])) levels(droplevels(df[[group_col]][ok]))
                  else sort(unique(as.character(g)))
  g <- factor(as.character(g), levels = group_levels)
  values <- split(v, g, drop = FALSE)
  # v0.42.3: keep groups with >=1 observation (was >=2). Need at
  # least 2 groups overall AND at least one group with n>=2 to make
  # any test meaningful (Wilcoxon with two singletons has no
  # discriminating power).
  values <- values[lengths(values) >= 1L]
  if (length(values) < 2L) return(NULL)
  if (max(lengths(values)) < 2L) return(NULL)
  group_names <- names(values)
  K <- length(values)
  n_per_group <- lengths(values)

  # Overall test
  overall <- if (K == 2L) {
    list(name = "Wilcoxon rank-sum",
         p    = suppressWarnings(stats::wilcox.test(values[[1]], values[[2]])$p.value))
  } else {
    list(name = "Kruskal-Wallis",
         p    = suppressWarnings(stats::kruskal.test(values)$p.value))
  }

  # Pairwise comparisons
  if (K <= 3L) {
    pairs <- combn(group_names, 2L, simplify = FALSE)
  } else {
    # adjacent pairs + extremes (i,i+1) + (1,K)
    adj <- lapply(seq_len(K - 1L), function(i) c(group_names[i], group_names[i + 1L]))
    pairs <- unique(c(adj, list(c(group_names[1L], group_names[K]))))
  }
  pw <- do.call(rbind, lapply(pairs, function(p) {
    t <- suppressWarnings(stats::wilcox.test(values[[p[1L]]], values[[p[2L]]]))
    data.frame(left = p[1L], right = p[2L], p = t$p.value,
                stringsAsFactors = FALSE)
  }))
  # Drop any pairs where wilcox couldn't return a p (e.g. n=1 vs n=1)
  pw <- pw[!is.na(pw$p), , drop = FALSE]
  if (!nrow(pw)) return(NULL)
  # v0.45.0: multiple-test correction across the pairwise family (omnibus raw)
  padj <- .resolve_p_adjust(p_adjust)
  pw$p_raw <- pw$p
  if (padj != "none") pw$p <- stats::p.adjust(pw$p, method = padj)
  list(overall = overall, pairwise = pw, n_per_group = n_per_group,
       group_levels = group_levels, p_adjust = padj)
}

# Compute scatter-family stats: Spearman correlation by default
# (robust to outliers and works on log/linear identically). Returns
# list(method, rho, p, n) or NULL.
.compute_scatter_stats <- function(df, x_col, y_col, method = "spearman") {
  if (!x_col %in% names(df) || !y_col %in% names(df)) return(NULL)
  x <- df[[x_col]]; y <- df[[y_col]]
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 3L) return(NULL)
  ct <- suppressWarnings(stats::cor.test(x, y, method = method))
  list(method = method, rho = unname(ct$estimate), p = ct$p.value,
       n = length(x))
}

# Append comparison brackets to a ggplot. Works on the transformed
# y-scale: brackets sit at log10 positions if use_log = TRUE, linear
# otherwise. Returns the modified ggplot.
#
# Args:
#   p           ggplot already constructed (boxes / points / etc.)
#   stats_res   .compute_box_stats() output
#   df          the same data passed to ggplot (used for y range)
#   group_col   column on df used for the x mapping
#   value_col   column on df used for the y mapping
#   y_scale     "log10" / "log2" / "linear" (default linear)
#   expand_top  upper expand fraction (default 0.18)
# v0.44.6: brackets drawn between DODGED sub-positions (not x-axis ticks).
#   mode "within_group"    -> within each group_col level, compare subgroup_col levels
#   mode "within_subgroup" -> within each subgroup_col level, compare group_col levels
# pairwise = TRUE draws every pair; FALSE draws one overall-p bracket per cluster.
.add_dodge_brackets <- function(p, df, group_col, subgroup_col, value_col,
                                mode = c("within_group","within_subgroup"),
                                pairwise = TRUE, facet_data = NULL,
                                y_scale = "linear", dodge_w = 0.8, test = "auto",
                                p_adjust = "none") {
  mode <- match.arg(mode)
  if (!nrow(df)) return(p)
  gfac <- if (is.factor(df[[group_col]]))    df[[group_col]]    else factor(df[[group_col]])
  sfac <- if (is.factor(df[[subgroup_col]])) df[[subgroup_col]] else factor(df[[subgroup_col]])
  glev <- levels(gfac); G <- length(glev)
  slev <- levels(sfac); K <- length(slev)
  if (G < 1L || K < 2L) return(p)
  off    <- function(k) dodge_w * ((2 * k - 1 - K) / (2 * K))  # subgroup-k center offset
  subpos <- function(gi, k) gi + off(k)
  use_log <- y_scale %in% c("log10","log","log2")
  v <- df[[value_col]]; v <- v[is.finite(v)]
  if (!length(v)) return(p)
  ytop <- max(v, na.rm = TRUE)
  rng  <- range(v, na.rm = TRUE); span <- diff(rng)
  if (!is.finite(span) || span == 0) span <- abs(rng[2]) + 1
  y_at <- function(lvl) if (use_log) ytop * (1.5 ^ (lvl + 1)) else rng[2] + span * 0.09 * (lvl + 1)

  rows <- list()
  push <- function(x1, x2, lvl, p_lab)
    rows[[length(rows) + 1L]] <<- data.frame(
      x = min(x1, x2), xend = max(x1, x2), y = y_at(lvl), label = p_lab,
      stringsAsFactors = FALSE)

  anchors <- if (mode == "within_group") glev else slev
  for (ai in seq_along(anchors)) {
    if (mode == "within_group") {
      d_i <- df[as.character(gfac) == anchors[ai], , drop = FALSE]
      st  <- .compute_box_stats(d_i, group_col = subgroup_col, value_col = value_col, test = test, p_adjust = p_adjust)
      if (is.null(st)) next
      if (pairwise && !is.null(st$pairwise)) {
        pw <- st$pairwise
        for (r in seq_len(nrow(pw))) {
          ka <- match(pw$left[r], slev); kb <- match(pw$right[r], slev)
          if (!is.na(ka) && !is.na(kb))
            push(subpos(ai, ka), subpos(ai, kb), r - 1L, paste0("p=", .format_p_value(pw$p[r])))
        }
      } else push(subpos(ai, 1L), subpos(ai, K), 0L, paste0("p=", .format_p_value(st$overall$p)))
    } else {  # within_subgroup
      d_i <- df[as.character(sfac) == anchors[ai], , drop = FALSE]
      st  <- .compute_box_stats(d_i, group_col = group_col, value_col = value_col, test = test, p_adjust = p_adjust)
      if (is.null(st)) next
      if (pairwise && !is.null(st$pairwise)) {
        pw <- st$pairwise
        for (r in seq_len(nrow(pw))) {
          ia <- match(pw$left[r], glev); ib <- match(pw$right[r], glev)
          if (!is.na(ia) && !is.na(ib))
            push(subpos(ia, ai), subpos(ib, ai), (ai - 1L) + (r - 1L), paste0("p=", .format_p_value(pw$p[r])))
        }
      } else push(subpos(1L, ai), subpos(G, ai), ai - 1L, paste0("p=", .format_p_value(st$overall$p)))
    }
  }
  if (!length(rows)) return(p)
  bdf <- do.call(rbind, rows)
  if (!is.null(facet_data)) for (nm in names(facet_data)) bdf[[nm]] <- facet_data[[nm]]
  p +
    ggplot2::geom_segment(data = bdf,
      mapping = ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$y),
      inherit.aes = FALSE, linewidth = 0.4, color = "grey20") +
    ggplot2::geom_text(data = bdf,
      mapping = ggplot2::aes(x = (.data$x + .data$xend) / 2, y = .data$y, label = .data$label),
      inherit.aes = FALSE, vjust = -0.3, size = 2.7, color = "grey15")
}

.add_brackets <- function(p, stats_res, df, group_col, value_col,
                          y_scale = "linear", expand_top = 0.28,
                          base_off = NULL, x_lookup = NULL,
                          facet_data = NULL) {
  if (is.null(stats_res) || is.null(stats_res$pairwise) ||
      !nrow(stats_res$pairwise)) return(p)
  pw <- stats_res$pairwise

  # v0.42.5: x_lookup lets the caller override the default
  # group-name -> x-position mapping. This is needed when bracket
  # comparisons run between dodged sub-positions (e.g. compare_by =
  # "color" in ctdna_boxplot where K colour groups occupy K dodged
  # positions within each integer x). If NULL, fall back to integer
  # match() against group_levels as before.
  group_levels <- stats_res$group_levels
  if (is.null(x_lookup)) {
    pw$x_left  <- match(pw$left,  group_levels)
    pw$x_right <- match(pw$right, group_levels)
  } else {
    pw$x_left  <- x_lookup[pw$left]
    pw$x_right <- x_lookup[pw$right]
  }
  ok <- !is.na(pw$x_left) & !is.na(pw$x_right)
  pw <- pw[ok, , drop = FALSE]
  if (!nrow(pw)) return(p)
  pw$width <- pw$x_right - pw$x_left
  # Stack: narrower comparisons low, wider ones above
  pw <- pw[order(pw$width, pw$x_left), , drop = FALSE]

  # y-space: work in log10 if log-y, else linear
  use_log <- grepl("^log", y_scale)
  v_obs <- df[[value_col]]
  v_obs <- v_obs[is.finite(v_obs)]
  if (use_log) v_obs <- v_obs[v_obs > 0]
  if (!length(v_obs)) return(p)
  v_obs_t <- if (use_log) log10(v_obs) else v_obs

  q05    <- as.numeric(stats::quantile(v_obs_t, 0.05, na.rm = TRUE))
  q25    <- as.numeric(stats::quantile(v_obs_t, 0.25, na.rm = TRUE))
  q75    <- as.numeric(stats::quantile(v_obs_t, 0.75, na.rm = TRUE))
  q95    <- as.numeric(stats::quantile(v_obs_t, 0.95, na.rm = TRUE))
  iqr    <- q75 - q25
  upper_fence  <- q75 + 1.5 * iqr
  v_top_robust <- min(q95, upper_fence, na.rm = TRUE)
  if (!is.finite(v_top_robust) || v_top_robust <= q25)
    v_top_robust <- max(v_obs_t)
  span_robust <- max(q95 - q05, iqr, abs(v_top_robust) * 0.1, 1e-6)
  if (!is.finite(span_robust) || span_robust <= 0)
    span_robust <- abs(v_top_robust) + 1
  v_data_max <- max(v_obs_t)

  # v0.42.5 — anchor logic depends on scale.
  #   LOG: original heuristic (anchor at data_max if outliers modest,
  #        at robust top otherwise). This handles log's wide span.
  #   LINEAR: always anchor at v_data_max. The robust-top swap is a
  #        log-only protection; on linear scale, capping below data_max
  #        causes brackets to render INSIDE the boxes for any group
  #        with non-skewed data (e.g. RECIST=PD whisker above v_top).
  v_anchor <- if (use_log) {
                outlier_gap <- v_data_max - v_top_robust
                if (outlier_gap <= 1.5 * span_robust) v_data_max
                else v_top_robust
              } else {
                v_data_max
              }

  # v0.42.5 — base_off is now exposed as a parameter (was hardcoded
  # 0.07 in log, 0.05*span in linear). Default = 0.04 (log) or
  # 0.025*span_robust (linear); callers can pass tighter values for
  # bracket stacks that need to fit a constrained ylim.
  if (is.null(base_off)) {
    base_off <- if (use_log) 0.04 else 0.025 * span_robust
  }
  if (use_log) {
    max_y_band <- expand_top   # 0.28 log units default (was 0.40)
  } else {
    max_y_band <- expand_top * span_robust * 0.75
  }
  n_b <- nrow(pw)
  # v0.42.4 hard-cap on bracket-driven y expansion preserved.
  budget_step <- (max_y_band - base_off) / max(n_b - 0.35, 0.65)

  min_step <- if (use_log) 0.10 else max(0.04 * span_robust, 0.025 * abs(v_anchor))
  step <- if (budget_step < min_step) max(budget_step, 0.001)
          else max(min_step, min(budget_step,
                                  (expand_top * 0.80) * span_robust / n_b))

  pw$y_lin   <- v_anchor + base_off + (seq_len(n_b) - 1L) * step
  pw$y_text  <- pw$y_lin + step * 0.45
  tick_drop  <- step * 0.30
  pw$y_tick  <- pw$y_lin - tick_drop

  # Top of stack in transformed space (with a tiny final margin)
  top_in_t <- max(pw$y_text) + step * 0.20

  if (use_log) {
    pw$y      <- 10 ^ pw$y_lin
    pw$y_text <- 10 ^ pw$y_text
    pw$y_tick <- 10 ^ pw$y_tick
    top_in_data <- 10 ^ top_in_t
  } else {
    pw$y <- pw$y_lin
    top_in_data <- top_in_t
  }
  pw$label <- vapply(pw$p, .format_p_value, character(1))

  # Build geom data
  seg <- do.call(rbind, lapply(seq_len(n_b), function(i) {
    data.frame(
      x    = c(pw$x_left[i],  pw$x_left[i],  pw$x_right[i]),
      xend = c(pw$x_left[i],  pw$x_right[i], pw$x_right[i]),
      y    = c(pw$y_tick[i],  pw$y[i],       pw$y[i]),
      yend = c(pw$y[i],       pw$y[i],       pw$y_tick[i]),
      stringsAsFactors = FALSE)
  }))
  txt <- data.frame(
    x = (pw$x_left + pw$x_right) / 2,
    y = pw$y_text,
    label = pw$label,
    stringsAsFactors = FALSE)

  # v0.42.4 — facet awareness. When .add_brackets is called once per
  # facet (e.g. ctdna_plot_baseline, ctdna_plot_pct_change_by_dose_time),
  # the caller passes facet_data with the panel's facet column
  # values. We attach those to the bracket geom data so ggplot's
  # facet_wrap/facet_grid renders each bracket layer in its own
  # panel only — no cross-leak. Without this, every facet shows
  # every facet's brackets.
  if (!is.null(facet_data) && length(facet_data)) {
    for (k in names(facet_data)) {
      seg[[k]] <- facet_data[[k]]
      txt[[k]] <- facet_data[[k]]
    }
  }

  # Expand y limits to fit the bracket stack. Don't extend beyond
  # top_in_data — this is the hard cap that v0.42.4 enforces.
  p <- p + ggplot2::expand_limits(y = c(NA_real_, top_in_data))

  p +
    ggplot2::geom_segment(data = seg,
       mapping     = ggplot2::aes(x = .data$x, xend = .data$xend,
                                    y = .data$y, yend = .data$yend),
       inherit.aes = FALSE,
       color       = "grey25",
       linewidth   = 0.4) +
    ggplot2::geom_text(data = txt,
       mapping     = ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
       inherit.aes = FALSE,
       size        = 2.8,
       color       = "grey15",
       vjust       = -0.1)
}

# Build the caption text from stats results. For box family: "Test:
# Wilcoxon rank-sum; overall p = 0.034; N = 8, 12, 11" or similar.
# For scatter: "Spearman rho = 0.42; p = 0.013; N = 31".
.compose_stats_caption <- function(stats_res, kind = c("box","scatter"),
                                    base_caption = NULL) {
  kind <- match.arg(kind)
  if (is.null(stats_res)) return(base_caption)
  if (kind == "box") {
    n_str <- paste(stats_res$n_per_group, collapse = ", ")
    extra <- sprintf("Test: %s; overall p = %s; N per group = %s",
                      stats_res$overall$name,
                      .format_p_value(stats_res$overall$p),
                      n_str)
  } else {
    extra <- sprintf("Spearman \u03c1 = %.2f; p = %s; N = %d",
                      stats_res$rho, .format_p_value(stats_res$p), stats_res$n)
  }
  if (is.null(base_caption) || !nzchar(base_caption)) extra
  else paste(base_caption, extra, sep = "  |  ")
}
