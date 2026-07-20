# ============================================================================
# ctdna_concordance_oncoprint_core()  --  step 2 (full signature restored)
# ----------------------------------------------------------------------------
# Signature is IDENTICAL to ctdna_oncoprint(), plus one extra `tissue_df` arg.
# Body is a copy of ctdna_oncoprint's body up to the matrix build, then a
# small tissue-expansion block that doubles global_mat's columns to
# [P1|ctDNA, P1|Tissue, P2|ctDNA, P2|Tissue, ...], then continues with
# ctdna_oncoprint's downstream code verbatim.
#
# The only thing added beyond copying ctdna_oncoprint is: (1) the tissue
# expansion, and (2) a temporary override of .oncoprint_sort_patients that
# keeps pairs adjacent even when the panel renderer re-sorts by burden.
# ============================================================================

ctdna_concordance_oncoprint_core <- function(
    prep,
    tissue_df          = NULL,
    gene_sets,
    filter_scheme      = NULL,
    visit              = "C1D1",
    visit_col          = "Visit_name",
    wrap               = NULL,
    top_annotations    = NULL,
    annotation_labels  = NULL,
    group_patients_by  = NULL,
    scheme             = c("three", "two", "four", "raw"),
    alterations        = NULL,
    sort_genes         = c("global", "within_set", "none", "overlap"),
    sort_by            = c("burden", "overlap"),
    show_all_patients  = TRUE,
    show_patient_names = TRUE,
    show_freq_bar      = TRUE,
    show_gene_set_annotation = TRUE,
    engine             = c("auto", "complexheatmap", "ggplot"),
    title              = NULL,
    subtitle           = NULL,
    caption            = NULL,
    patient_col        = "Patient_ID",
    gene_col           = "Gene",
    variant_col        = NULL,
    indications        = NULL) {

  # Default tissue_df to prep$dnaseq if not provided
  if (is.null(tissue_df)) {
    if (is.null(prep$dnaseq) || !is.data.frame(prep$dnaseq) ||
        !nrow(prep$dnaseq))
      stop("tissue_df is NULL and prep$dnaseq is empty. ",
           "Attach tissue data via ctdna_prep_add(prep, dnaseq = <df>) ",
           "or pass tissue_df = <data.frame> explicitly.", call. = FALSE)
    tissue_df <- prep$dnaseq
  }

  # ==== VERBATIM FROM ctdna_oncoprint ======================================
  prep <- ctdnaTM:::.ctdna_filter_prep_by_indication(prep, indications)
  engine     <- match.arg(engine)
  scheme     <- match.arg(scheme)
  sort_genes <- match.arg(sort_genes)
  sort_by    <- match.arg(sort_by)
  ctdnaTM:::.reject_cohort_args(top_annotations   = top_annotations,
                                 group_patients_by = group_patients_by,
                                 wrap              = wrap,
                                 annotation_labels = names(annotation_labels))
  if (missing(gene_sets) || is.null(gene_sets))
    stop("gene_sets is required (named list, or built-in name like \"HRR14\").",
         call. = FALSE)

  base <- ctdnaTM:::.build_landscape_base(prep, filter_scheme, visit, visit_col,
                                            patient_col, gene_col,
                                            "ctdna_concordance_oncoprint_core")
  df           <- base$variants
  patient_data <- base$patient_data
  pd_id_col    <- patient_col

  has_ch <- requireNamespace("ComplexHeatmap", quietly = TRUE)
  if (engine == "auto") engine <- if (has_ch) "complexheatmap" else "ggplot"
  if (engine == "complexheatmap" && !has_ch) {
    message("ComplexHeatmap not installed; falling back to ggplot engine.")
    engine <- "ggplot"
  }

  df <- ctdnaTM:::.oncoprint_classify(df)
  if (!is.null(alterations))
    df <- df[df$Alteration_class %in% alterations, , drop = FALSE]
  if (!nrow(df))
    stop("no variants remain after alteration filter.", call. = FALSE)

  rgs        <- ctdnaTM:::.oncoprint_resolve_gene_sets(df, gene_sets, gene_col)
  resolved   <- rgs$resolved
  gene_union <- rgs$gene_union

  if (is.null(top_annotations) && "RECIST" %in% names(patient_data) &&
      any(!is.na(patient_data$RECIST)))
    top_annotations <- "RECIST"

  sch <- ctdnaTM:::.oncoprint_apply_scheme(patient_data, top_annotations,
                                            group_patients_by, scheme)
  patient_data     <- sch$patient_data
  recist_cols_used <- sch$recist_cols_used

  explicit_cohort <- if (is.character(show_all_patients))
                       show_all_patients else NULL
  keep_zero <- isTRUE(show_all_patients) || !is.null(explicit_cohort)
  # For gene-overlap sort, run build_matrix in "global" mode; reorder rows
  # after tissue_mat is available. Otherwise pass sort_genes through as-is.
  internal_sort_genes <- if (identical(sort_genes, "overlap"))
                            "global" else sort_genes
  built <- ctdnaTM:::.oncoprint_build_matrix(df, patient_data, pd_id_col,
                                              patient_col, gene_col,
                                              variant_col, gene_union,
                                              resolved, internal_sort_genes,
                                              explicit_cohort = explicit_cohort)
  global_mat    <- built$global_mat
  ordered_genes <- built$ordered_genes
  row_split     <- built$row_split

  # ==== TISSUE-COLUMN EXPANSION (the only real addition) ==================
  patients_orig <- colnames(global_mat)

  # Restrict tissue_df to the ctDNA cohort, then classify + optionally filter
  tissue_df <- tissue_df[
    as.character(tissue_df[[patient_col]]) %in% patients_orig, , drop = FALSE]
  tissue_df <- .co_classify_tissue_local(tissue_df)
  if (!is.null(alterations))
    tissue_df <- tissue_df[tissue_df$Alteration_class %in% alterations, ,
                            drop = FALSE]

  # Build tissue_mat parallel to global_mat (same rows, same patient cols)
  tissue_mat <- matrix("", nrow = length(ordered_genes),
                        ncol = length(patients_orig),
                        dimnames = list(ordered_genes, patients_orig))
  if (nrow(tissue_df)) {
    t_kept <- tissue_df[tissue_df[[gene_col]] %in% ordered_genes, , drop = FALSE]
    if (nrow(t_kept)) {
      prio <- c("Focal_Amp","Amp","Homozygous_Del","LOH",
                "Truncating","Missense","InFrame","Promoter",
                "Fusion","LGR","Other")
      t_kept$.p <- match(t_kept$Alteration_class, prio)
      t_kept$.p[is.na(t_kept$.p)] <- length(prio) + 1L
      t_kept <- t_kept[order(t_kept[[patient_col]], t_kept[[gene_col]],
                               t_kept$.p), ]
      t_kept <- t_kept[!duplicated(
        t_kept[, c(patient_col, gene_col)]), ]
      ri <- match(t_kept[[gene_col]],    rownames(tissue_mat))
      ci <- match(t_kept[[patient_col]], colnames(tissue_mat))
      ok <- !is.na(ri) & !is.na(ci)
      tissue_mat[cbind(ri[ok], ci[ok])] <- t_kept$Alteration_class[ok]
    }
  }

  # ---- OVERLAP gene sort ---------------------------------------------------
  # Reorder rows by Jaccard overlap between ctDNA and tissue calls per gene.
  # Ties broken by combined burden. Applied AFTER tissue_mat is built so we
  # have both matrices to compare.
  if (identical(sort_genes, "overlap")) {
    ct_hits <- global_mat[ordered_genes, patients_orig, drop = FALSE] != ""
    ti_hits <- tissue_mat != ""
    gene_uni   <- rowSums(ct_hits | ti_hits)
    gene_inter <- rowSums(ct_hits & ti_hits)
    gene_score <- ifelse(gene_uni > 0, gene_inter / gene_uni, 0)
    gene_burden <- rowSums(ct_hits) + rowSums(ti_hits)
    new_ord <- order(-gene_score, -gene_burden, ordered_genes)
    ordered_genes <- ordered_genes[new_ord]
    global_mat    <- global_mat[ordered_genes, , drop = FALSE]
    tissue_mat    <- tissue_mat[ordered_genes, , drop = FALSE]
    if (!is.null(row_split)) row_split <- row_split[new_ord]
  }

  # Interleave columns
  ctdna_suffix  <- " \u00B7 ctDNA"
  tissue_suffix <- " \u00B7 Tissue"
  paired_cols <- as.vector(rbind(paste0(patients_orig, ctdna_suffix),
                                   paste0(patients_orig, tissue_suffix)))
  paired_mat  <- matrix("", nrow = length(ordered_genes),
                          ncol = length(paired_cols),
                          dimnames = list(ordered_genes, paired_cols))
  paired_mat[, seq(1L, ncol(paired_mat), by = 2L)] <-
    global_mat[ordered_genes, patients_orig]
  paired_mat[, seq(2L, ncol(paired_mat), by = 2L)] <- tissue_mat

  # Expand patient_data: duplicate each row with the two suffixes so that
  # annotation lookup (match on pd_id_col) resolves for both.
  pd_c <- patient_data
  pd_c[[pd_id_col]] <- paste0(patient_data[[pd_id_col]], ctdna_suffix)
  pd_t <- patient_data
  pd_t[[pd_id_col]] <- paste0(patient_data[[pd_id_col]], tissue_suffix)
  patient_data <- rbind(pd_c, pd_t)

  # Replace global_mat + explicit_cohort so panel renderers see the paired set.
  global_mat <- paired_mat
  if (is.character(explicit_cohort))
    explicit_cohort <- paired_cols   # keep the interleave the caller asked for

  # Global alteration-class set now must include tissue-only classes.
  tissue_classes <- setdiff(unique(as.vector(tissue_mat)), "")

  # ==== TEMPORARY OVERRIDE of .oncoprint_sort_patients =====================
  # Keeps ctDNA / Tissue pairs adjacent. Sorts pairs by combined burden
  # (highest first), respecting group_patients_by when provided.
  paired_sort <- function(sub, burden, group_patients_by, patient_data,
                            pd_id_col, recist_cols_used) {
    cn <- colnames(sub)
    is_c <- endsWith(cn, ctdna_suffix)
    is_t <- endsWith(cn, tissue_suffix)
    pair_id <- cn
    pair_id[is_c] <- substr(cn[is_c], 1L,
                              nchar(cn[is_c]) - nchar(ctdna_suffix))
    pair_id[is_t] <- substr(cn[is_t], 1L,
                              nchar(cn[is_t]) - nchar(tissue_suffix))
    modality <- rep(3L, length(cn))
    modality[is_c] <- 1L
    modality[is_t] <- 2L
    pair_burden <- ave(burden, pair_id, FUN = sum)

    # Overlap score per pair (Jaccard: |ctDNA n tissue| / |ctDNA u tissue|).
    # 0 for pairs with no alterations in either modality.
    if (identical(sort_by, "overlap")) {
      unique_pairs <- unique(pair_id)
      pair_score <- vapply(unique_pairs, function(pid) {
        cc <- paste0(pid, ctdna_suffix)
        tc <- paste0(pid, tissue_suffix)
        ca <- if (cc %in% cn) sub[, cc] != "" else rep(FALSE, nrow(sub))
        ta <- if (tc %in% cn) sub[, tc] != "" else rep(FALSE, nrow(sub))
        n_union <- sum(ca | ta)
        if (n_union == 0L) 0 else sum(ca & ta) / n_union
      }, numeric(1))
      names(pair_score) <- unique_pairs
      overlap_per_col <- as.numeric(pair_score[pair_id])
    } else {
      overlap_per_col <- rep(0, length(cn))
    }

    if (!is.null(group_patients_by) && !is.null(patient_data) &&
        !is.na(pd_id_col) && group_patients_by %in% names(patient_data)) {
      raw_grp <- patient_data[[group_patients_by]]
      grp <- as.character(raw_grp)[
        match(cn, as.character(patient_data[[pd_id_col]]))]
      grp[is.na(grp) | !nzchar(grp)] <- "NA"
      grp_pair <- ave(grp, pair_id, FUN = function(x)
        x[!is.na(x) & nzchar(x)][1])
      grp_pair[is.na(grp_pair) | !nzchar(grp_pair)] <- "NA"

      level_order <- if (group_patients_by %in% recist_cols_used) {
        u <- unique(grp_pair)
        c(intersect(ctdnaTM:::.oncoprint_recist_order, u),
          setdiff(u, ctdnaTM:::.oncoprint_recist_order))
      } else if (is.factor(raw_grp)) {
        intersect(levels(raw_grp), unique(grp_pair))
      } else sort(unique(grp_pair))

      if (identical(sort_by, "overlap")) {
        sub[, order(match(grp_pair, level_order),
                      -overlap_per_col, -pair_burden,
                      pair_id, modality), drop = FALSE]
      } else {
        sub[, order(match(grp_pair, level_order),
                      -pair_burden, pair_id, modality), drop = FALSE]
      }
    } else {
      if (identical(sort_by, "overlap")) {
        sub[, order(-overlap_per_col, -pair_burden,
                      pair_id, modality), drop = FALSE]
      } else {
        sub[, order(-pair_burden, pair_id, modality), drop = FALSE]
      }
    }
  }
  old_sort <- ctdnaTM:::.oncoprint_sort_patients
  assignInNamespace(".oncoprint_sort_patients", paired_sort, ns = "ctdnaTM")
  on.exit(assignInNamespace(".oncoprint_sort_patients", old_sort,
                              ns = "ctdnaTM"), add = TRUE)

  # ==== TEMPORARY OVERRIDE of .oncoprint_panel_ch ==========================
  # Returns list(heatmap, freq_anno). The freq bar is built as a STANDALONE
  # HeatmapAnnotation (which = "column"), NOT as bottom_annotation of the
  # heatmap. This is because %v% collapses middle heatmaps' top/bottom
  # annotations to first/last only. Putting the freq bar as its own element
  # in the vertical stack keeps a per-panel freq bar visible.
  paired_panel_ch <- function(global_mat, ordered_genes, row_split,
                                panel_pids, panel_title,
                                is_first, is_last,
                                patient_data, pd_id_col,
                                group_patients_by, recist_cols_used,
                                ann_col_global, ann_at_global,
                                annotation_labels,
                                resolved, sort_genes,
                                keep_zero_burden_cols, show_patient_names,
                                show_freq_bar,
                                alter_fun, col_vec, global_classes) {
    panel_pids <- unique(panel_pids)
    sub <- global_mat[ordered_genes,
                       intersect(colnames(global_mat), panel_pids),
                       drop = FALSE]
    burden <- colSums(sub != "")
    if (!isTRUE(keep_zero_burden_cols)) {
      keep_pat <- names(burden)[burden > 0]
      sub <- sub[, keep_pat, drop = FALSE]; burden <- burden[keep_pat]
    }
    if (ncol(sub) == 0L) return(NULL)

    sub <- ctdnaTM:::.oncoprint_sort_patients(
      sub, burden, group_patients_by,
      patient_data, pd_id_col, recist_cols_used)

    # ---- Patient annotations (Dose, RECIST) as rowAnnotation on the LEFT
    patient_row_anno <- NULL
    if (length(ann_col_global)) {
      idx <- match(colnames(sub), as.character(patient_data[[pd_id_col]]))
      ann_args <- list()
      for (ann in names(ann_col_global))
        ann_args[[ann]] <- as.character(patient_data[[ann]])[idx]
      ann_names <- names(ann_args)
      ann_labels <- ann_names
      if (!is.null(annotation_labels)) {
        rn <- annotation_labels[ann_names]
        ann_labels <- ifelse(is.na(rn) | !nzchar(rn), ann_names, rn)
      }
      patient_row_anno <- do.call(ComplexHeatmap::HeatmapAnnotation,
        c(ann_args,
          list(which                = "row",
               col                  = ann_col_global[ann_names],
               na_col               = "#E0E0E0",
               annotation_label     = ann_labels,
               annotation_name_side = "top",
               annotation_name_gp   = grid::gpar(fontsize = 8),
               simple_anno_size     = grid::unit(4, "mm"),
               gap                  = grid::unit(0.5, "mm"),
               show_legend          = FALSE)))
    }

    # ---- TOP annotation: per-gene stacked alteration count SPLIT by
    # source. Two adjacent stacked bars per gene: ctDNA (left) and Tissue
    # (right). Each stacked by alteration type. Each labeled with a %
    # using its own sample count as denominator.
    freq_anno <- NULL
    if (isTRUE(show_freq_bar)) {
      types_here <- names(col_vec)
      cn         <- colnames(sub)
      is_c       <- endsWith(cn, ctdna_suffix)
      is_t       <- endsWith(cn, tissue_suffix)
      sub_c      <- sub[, is_c,  drop = FALSE]
      sub_t      <- sub[, is_t,  drop = FALSE]

      per_gene_c <- matrix(0L, nrow = nrow(sub), ncol = length(types_here),
                            dimnames = list(rownames(sub), types_here))
      per_gene_t <- per_gene_c
      count_cell <- function(mat, dest) {
        if (!ncol(mat)) return(dest)
        for (gi in seq_len(nrow(mat))) {
          for (pj in seq_len(ncol(mat))) {
            cell <- mat[gi, pj]
            if (is.na(cell) || !nzchar(cell)) next
            for (p in unlist(strsplit(cell, ";", fixed = TRUE))) {
              if (p %in% types_here) dest[gi, p] <- dest[gi, p] + 1L
            }
          }
        }
        dest
      }
      per_gene_c <- count_cell(sub_c, per_gene_c)
      per_gene_t <- count_cell(sub_t, per_gene_t)

      bar_cols <- unname(col_vec[types_here])
      n_c      <- max(ncol(sub_c), 1L)
      n_t      <- max(ncol(sub_t), 1L)
      alt_c    <- !is.na(sub_c) & sub_c != ""
      alt_t    <- !is.na(sub_t) & sub_t != ""
      n_alt_c  <- as.integer(rowSums(alt_c, na.rm = TRUE))
      n_alt_t  <- as.integer(rowSums(alt_t, na.rm = TRUE))
      pct_c    <- as.integer(round(100 * n_alt_c / n_c))
      pct_t    <- as.integer(round(100 * n_alt_t / n_t))
      pct_c[is.na(pct_c)] <- 0L
      pct_t[is.na(pct_t)] <- 0L
      max_h    <- max(c(rowSums(per_gene_c), rowSums(per_gene_t), 1L))

      stacked_pct_fun <- ComplexHeatmap::AnnotationFunction(
        fun = function(index, k, n_slice) {
          n_here <- length(index)
          grid::pushViewport(grid::viewport(
            xscale = c(0, n_here),
            yscale = c(0, max_h_v * 1.6)))
          bar_w <- 0.36
          off   <- 0.20   # horizontal offset of each bar center from cell center
          for (i in seq_along(index)) {
            gi <- index[i]
            # ---- ctDNA stack (left half) ----
            counts_c <- per_gene_c_v[gi, ]
            tot_c    <- sum(counts_c)
            y0 <- 0
            for (t_idx in seq_along(counts_c)) {
              cnt <- counts_c[t_idx]; if (cnt == 0L) next
              grid::grid.rect(
                x     = grid::unit(i - 0.5 - off, "native"),
                y     = grid::unit(y0, "native"),
                width = grid::unit(bar_w, "native"),
                height = grid::unit(cnt, "native"),
                just  = c("center", "bottom"),
                gp    = grid::gpar(fill = colors_v[t_idx], col = NA))
              y0 <- y0 + cnt
            }
            if (tot_c > 0L) {
              grid::grid.text(
                sprintf("%d%%", pct_c_v[gi]),
                x = grid::unit(i - 0.5 - off, "native"),
                y = grid::unit(tot_c, "native") + grid::unit(0.8, "mm"),
                rot = 90,
                just = c("left", "center"),
                gp = grid::gpar(fontsize = 6, fontface = "bold",
                                col = "#2C7FB8"))
            }
            # ---- Tissue stack (right half) ----
            counts_t <- per_gene_t_v[gi, ]
            tot_t    <- sum(counts_t)
            y0 <- 0
            for (t_idx in seq_along(counts_t)) {
              cnt <- counts_t[t_idx]; if (cnt == 0L) next
              grid::grid.rect(
                x     = grid::unit(i - 0.5 + off, "native"),
                y     = grid::unit(y0, "native"),
                width = grid::unit(bar_w, "native"),
                height = grid::unit(cnt, "native"),
                just  = c("center", "bottom"),
                gp    = grid::gpar(fill = colors_v[t_idx], col = NA))
              y0 <- y0 + cnt
            }
            if (tot_t > 0L) {
              grid::grid.text(
                sprintf("%d%%", pct_t_v[gi]),
                x = grid::unit(i - 0.5 + off, "native"),
                y = grid::unit(tot_t, "native") + grid::unit(0.8, "mm"),
                rot = 90,
                just = c("left", "center"),
                gp = grid::gpar(fontsize = 6, fontface = "bold",
                                col = "#D95F0E"))
            }
          }
          # ---- Vertical Y-axis on the LEFT (perpendicular to concordance
          # axis) showing the count scale (actual event numbers).
          tick_at <- pretty(c(0, max_h_v), n = 3)
          tick_at <- tick_at[tick_at <= max_h_v & tick_at >= 0]
          if (!length(tick_at)) tick_at <- c(0, max_h_v)
          axis_x <- grid::unit(0, "native") - grid::unit(1, "mm")
          grid::grid.lines(
            x = grid::unit.c(axis_x, axis_x),
            y = grid::unit(range(tick_at), "native"),
            gp = grid::gpar(col = "black", lwd = 0.8))
          for (tv in tick_at) {
            grid::grid.lines(
              x = grid::unit.c(axis_x, axis_x - grid::unit(1, "mm")),
              y = grid::unit(c(tv, tv), "native"),
              gp = grid::gpar(col = "black", lwd = 0.6))
            grid::grid.text(
              as.character(tv),
              x = axis_x - grid::unit(1.5, "mm"),
              y = grid::unit(tv, "native"),
              just = c("right", "center"),
              gp = grid::gpar(fontsize = 7))
          }
          grid::grid.text(
            "n altered",
            x = axis_x - grid::unit(6, "mm"),
            y = grid::unit(mean(range(tick_at)), "native"),
            rot = 90,
            just = c("center", "center"),
            gp = grid::gpar(fontsize = 8))
          grid::popViewport()
        },
        fun_name    = "split_stacked_pct_bar",
        which       = "column",
        var_import  = list(per_gene_c_v = per_gene_c,
                            per_gene_t_v = per_gene_t,
                            colors_v     = bar_cols,
                            pct_c_v      = pct_c,
                            pct_t_v      = pct_t,
                            max_h_v      = max_h),
        n           = nrow(sub),
        subsettable = TRUE,
        height      = grid::unit(2.5, "cm"))

      freq_anno <- ComplexHeatmap::HeatmapAnnotation(
        "n altered"          = stacked_pct_fun,
        which                = "column",
        show_annotation_name = FALSE,
        height               = grid::unit(2.5, "cm"))
    }

    # ---- RIGHT annotation: per-PATIENT concordance (Jaccard between the
    # ctDNA and Tissue rows for the same patient). Both rows of a pair get
    # the same score so the bar visually spans the pair.
    conc_row_anno <- NULL
    {
      cn <- colnames(sub)   # includes " · ctDNA" / " · Tissue" suffixes
      strip <- function(s) {
        s <- sub(paste0(ctdna_suffix,  "$"), "", s)
        s <- sub(paste0(tissue_suffix, "$"), "", s)
        s
      }
      pat <- vapply(cn, strip, character(1))
      unique_pats <- unique(pat)
      conc_per_pat <- setNames(rep(NA_real_, length(unique_pats)), unique_pats)
      for (pt in unique_pats) {
        c_col <- which(cn == paste0(pt, ctdna_suffix))
        t_col <- which(cn == paste0(pt, tissue_suffix))
        if (length(c_col) == 0L || length(t_col) == 0L) next
        c_alt <- which(sub[, c_col] != "")
        t_alt <- which(sub[, t_col] != "")
        if (length(c_alt) == 0L && length(t_alt) == 0L) {
          conc_per_pat[pt] <- NA_real_
        } else {
          conc_per_pat[pt] <- length(intersect(c_alt, t_alt)) /
                              max(length(union(c_alt, t_alt)), 1L)
        }
      }
      # Value per row of t_sub (same for ctDNA and Tissue rows of same pt)
      conc_per_row <- conc_per_pat[pat]
      # Also compute n_altered per pair (|union(ctDNA, Tissue)|) — same for
      # both rows of a pair. Companion metric next to the concordance ratio.
      n_alt_per_pat <- setNames(rep(0L, length(unique_pats)), unique_pats)
      for (pt in unique_pats) {
        c_col <- which(cn == paste0(pt, ctdna_suffix))
        t_col <- which(cn == paste0(pt, tissue_suffix))
        if (length(c_col) == 0L || length(t_col) == 0L) next
        c_alt <- which(sub[, c_col] != "")
        t_alt <- which(sub[, t_col] != "")
        n_alt_per_pat[pt] <- length(union(c_alt, t_alt))
      }
      n_alt_per_row <- n_alt_per_pat[pat]
      # Anchor NA to 0 for the bar (grey out via a companion column)
      bar_vals <- ifelse(is.na(conc_per_row), 0, conc_per_row)
      conc_row_anno <- ComplexHeatmap::rowAnnotation(
        "concordance\nscore" = ComplexHeatmap::anno_barplot(
          bar_vals,
          gp        = grid::gpar(fill = "#4C72B0", col = NA),
          which     = "row",
          border    = FALSE,
          bar_width = 0.95,
          ylim      = c(0, 1),
          width     = grid::unit(1.8, "cm"),
          axis      = TRUE,
          axis_param = list(side = "top", labels_rot = 0,
                            gp = grid::gpar(fontsize = 7))),
        show_annotation_name = TRUE,
        annotation_name_side = "top",
        annotation_name_gp   = grid::gpar(fontsize = 8),
        width = grid::unit(1.8, "cm"))
    }

    # ---- TRANSPOSE and render. Square cells: 5mm × 5mm.
    t_sub <- t(sub)   # paired_patients × genes
    row_h <- grid::unit(5, "mm")
    col_w <- grid::unit(5, "mm")

    ht <- ComplexHeatmap::oncoPrint(
      t_sub, alter_fun = alter_fun, col = col_vec,
      row_title            = panel_title,
      row_title_gp         = grid::gpar(fontsize = 11, fontface = "bold"),
      row_title_side       = "left",
      row_names_side       = "left",
      row_names_gp         = grid::gpar(fontsize = 6),
      column_names_gp      = grid::gpar(fontsize = 9),
      column_names_side    = "bottom",
      left_annotation      = patient_row_anno,
      top_annotation       = NULL,     # nothing at the top
      bottom_annotation    = NULL,     # freq bar is a stack element, not here
      right_annotation     = conc_row_anno,
      row_order            = seq_len(nrow(t_sub)),
      column_order         = seq_len(ncol(t_sub)),
      column_split         = row_split,
      row_gap              = grid::unit(0.5, "mm"),
      column_gap           = if (!is.null(row_split)) grid::unit(3, "mm")
                              else grid::unit(0.5, "mm"),
      width                = col_w * ncol(t_sub),
      heatmap_height       = row_h * nrow(t_sub),
      show_row_names       = isTRUE(show_patient_names),
      show_column_names    = FALSE,   # v0.76.1 — drawn via explicit stack element below
      show_pct             = FALSE,
      remove_empty_columns = FALSE,
      show_heatmap_legend  = FALSE)

    list(heatmap = ht, freq_anno = freq_anno)
  }
  old_panel_ch <- ctdnaTM:::.oncoprint_panel_ch
  assignInNamespace(".oncoprint_panel_ch", paired_panel_ch, ns = "ctdnaTM")
  on.exit(assignInNamespace(".oncoprint_panel_ch", old_panel_ch,
                              ns = "ctdnaTM"), add = TRUE)

  # ==== VERBATIM FROM ctdna_oncoprint (from legend down) ===================
  leg <- ctdnaTM:::.oncoprint_build_legend_state(top_annotations, patient_data,
                                                  recist_cols_used)
  ann_col_global <- leg$col
  ann_at_global  <- leg$at

  if (is.null(wrap)) {
    panel_values <- list(all = NULL)
  } else {
    if (length(wrap) != 1L || !is.character(wrap))
      stop("`wrap` must be a single column name or NULL.", call. = FALSE)
    if (!(wrap %in% names(patient_data)))
      stop("wrap column '", wrap, "' not found in patient_data. Available: ",
           paste(names(patient_data), collapse = ", "), ".", call. = FALSE)
    pv <- unique(as.character(patient_data[[wrap]]))
    pv <- pv[!is.na(pv) & nzchar(pv)]
    # If the caller passed `indications = c(...)`, honour that ORDER for the
    # panels. Values not in the indications list (e.g. wrap on a different
    # column) get appended alphabetically.
    if (!is.null(indications) && length(indications)) {
      ind_str <- as.character(indications)
      pv <- c(intersect(ind_str, pv), sort(setdiff(pv, ind_str)))
    } else {
      pv <- sort(pv)
    }
    panel_values <- as.list(pv)
    names(panel_values) <- unlist(panel_values)
  }

  # Include tissue-only classes so their color appears in the legend too.
  global_classes <- sort(unique(c(df$Alteration_class, tissue_classes)))
  alter_fun      <- ctdnaTM:::.oncoprint_alter_fun(global_classes)
  col_vec        <- ctdnaTM:::.oncoprint_alt_palette[global_classes]

  panels <- vector("list", length(panel_values))
  names(panels) <- names(panel_values)
  for (i in seq_along(panel_values)) {
    label <- names(panel_values)[i]
    if (is.null(wrap)) {
      pids <- intersect(unique(as.character(patient_data[[pd_id_col]])),
                          colnames(global_mat))
      panel_title <- ""
    } else {
      pids <- patient_data[[pd_id_col]][
        as.character(patient_data[[wrap]]) == label]
      pids <- intersect(unique(as.character(pids)), colnames(global_mat))
      panel_title <- label
    }
    if (!length(pids)) next

    args <- list(global_mat        = global_mat,
                 ordered_genes     = ordered_genes,
                 row_split         = row_split,
                 panel_pids        = pids,
                 panel_title       = panel_title,
                 is_first          = i == 1L,
                 is_last           = i == length(panel_values),
                 patient_data      = patient_data,
                 pd_id_col         = pd_id_col,
                 group_patients_by = group_patients_by,
                 recist_cols_used  = recist_cols_used,
                 ann_col_global    = ann_col_global,
                 ann_at_global     = ann_at_global,
                 annotation_labels = annotation_labels,
                 resolved          = resolved,
                 sort_genes        = sort_genes,
                 keep_zero_burden_cols = keep_zero,
                 show_patient_names = show_patient_names,
                 show_freq_bar      = show_freq_bar,
                 global_classes     = global_classes)
    if (engine == "complexheatmap") {
      args$alter_fun <- alter_fun
      args$col_vec   <- col_vec
      panels[[label]] <- do.call(ctdnaTM:::.oncoprint_panel_ch, args)
    } else {
      args$patient_col <- patient_col
      args$gene_col    <- gene_col
      args$variant_col <- variant_col
      panels[[label]]  <- do.call(ctdnaTM:::.oncoprint_panel_gg, args)
    }
  }
  panels <- panels[!vapply(panels, is.null, logical(1))]
  if (!length(panels))
    stop("no panel had any patient to render.", call. = FALSE)

  if (engine == "complexheatmap") {
    # Each CH panel result is list(heatmap, freq_anno). We build a vertical
    # stack of the form  [H1, freq1, H2, freq2, ..., gene_set_bottom]  so
    # that each panel's freq bar sits directly under that panel's heatmap
    # (as its own HeatmapAnnotation stack element, which %v% preserves),
    # and the gene-set membership strip appears only at the very bottom.
    #
    # `%v%` isn't in base R's operator table, so dispatch through the
    # ComplexHeatmap namespace.
    ch_vcat <- get("%v%", envir = asNamespace("ComplexHeatmap"))

    # -- Gene-set annotation built ONCE, ONLY at the bottom --------------
    gene_set_bottom_anno <- NULL
    set_names <- names(resolved)
    if (isTRUE(show_gene_set_annotation) &&
        length(set_names) >= 1L && !identical(sort_genes, "within_set") &&
        !(length(set_names) == 1L && set_names[1] == "All genes")) {
      set_pal <- c("#377EB8","#4DAF4A","#77EDDD","#E41A1C",
                   "#984EA3","#FF7F00","#A65628","#F781BF",
                   "#999999","#1F78B4")
      nsets <- length(set_names)
      fills <- if (nsets <= length(set_pal)) set_pal[seq_len(nsets)]
                else grDevices::colorRampPalette(set_pal)(nsets)
      names(fills) <- set_names
      gsets_df <- data.frame(row.names = ordered_genes,
                              check.names = FALSE,
                              stringsAsFactors = FALSE)
      gsets_col <- list()
      for (s in set_names) {
        gsets_df[[s]] <- ifelse(ordered_genes %in% resolved[[s]],
                                 s, NA_character_)
        gsets_col[[s]] <- setNames(fills[[s]], s)
      }
      short_label <- function(s) {
        lookup <- c("HRR14" = "HRR14", "TSG" = "TSG", "RTK" = "RTK",
                    "Cell_Cycle" = "Cell", "TP53_pathway" = "TP53",
                    "MMR" = "MMR", "PI3K" = "PI3K")
        if (s %in% names(lookup)) lookup[[s]] else
          substr(gsub("_", " ", s), 1, 8)
      }
      gene_set_bottom_anno <- ComplexHeatmap::HeatmapAnnotation(
        df                   = gsets_df,
        which                = "column",
        col                  = gsets_col,
        na_col               = "transparent",
        show_legend          = rep(FALSE, length(set_names)),
        annotation_label     = vapply(set_names, short_label, character(1)),
        annotation_name_side = "right",
        annotation_name_rot  = 0,
        annotation_name_gp   = grid::gpar(fontsize = 8),
        simple_anno_size     = grid::unit(3, "mm"))
    }

    # -- Stack: freq_anno BEFORE each heatmap so it appears ON TOP of that
    # panel's body (matching original ctdna_oncoprint()'s top_annotation
    # position, just with our new ctDNA/Tissue bars). Gene set at very end.
    #   [freq1, H1, freq2, H2, ..., gene_set]
    # Record what's what for ht_gap computation.
    stack_list <- list()
    kind       <- character(0)   # "freq" or "heatmap" or "geneset"
    for (p in panels) {
      if (is.null(p)) next
      if (!is.null(p$freq_anno)) {
        stack_list[[length(stack_list) + 1L]] <- p$freq_anno
        kind[length(stack_list)] <- "freq"
      }
      stack_list[[length(stack_list) + 1L]] <- p$heatmap
      kind[length(stack_list)] <- "heatmap"
    }
    if (!is.null(gene_set_bottom_anno)) {
      # v0.76.1: draw gene names via an explicit anno_text stack element
      # ABOVE the gene-set strip so the two never collide. Column names
      # via show_column_names on the Heatmap itself have a layout race
      # with the following annotation strip that CH doesn't resolve.
      gene_names_anno <- ComplexHeatmap::HeatmapAnnotation(
        genes = ComplexHeatmap::anno_text(
          ordered_genes,
          gp       = grid::gpar(fontsize = 9),
          rot      = 90,
          just     = "right",
          location = grid::unit(1, "npc")),
        which                = "column",
        show_annotation_name = FALSE,
        annotation_height    = grid::unit(20, "mm"))
      stack_list[[length(stack_list) + 1L]] <- gene_names_anno
      kind[length(stack_list)] <- "genelabel"
      stack_list[[length(stack_list) + 1L]] <- gene_set_bottom_anno
      kind[length(stack_list)] <- "geneset"
    }

    obj <- Reduce(ch_vcat, stack_list)

    # ht_gap between consecutive elements:
    #   freq       -> heatmap   : 0mm   (freq attaches to top of its panel)
    #   heatmap    -> freq      : 5mm   (visual separation between panels)
    #   heatmap    -> genelabel : 2mm   (names hug the bottom of last panel)
    #   genelabel  -> geneset   : 4mm   (small visual separation)
    ht_gap_custom <- NULL
    if (length(stack_list) > 1L) {
      gaps <- vapply(seq_len(length(stack_list) - 1L), function(i) {
        this_k <- kind[i]; next_k <- kind[i + 1L]
        if      (this_k == "freq"      && next_k == "heatmap")   0
        else if (this_k == "heatmap"   && next_k == "freq")      5
        else if (this_k == "heatmap"   && next_k == "genelabel") 2
        else if (this_k == "genelabel" && next_k == "geneset")   4
        else if (this_k == "heatmap"   && next_k == "geneset")   12
        else                                                     2
      }, numeric(1))
      ht_gap_custom <- grid::unit(gaps, "mm")
    }
    # ---- Build the full legend set ourselves so ComplexHeatmap doesn't
    # trim it to whatever appears in the first panel. Every value that
    # exists ANYWHERE in the plot gets a legend key, coloured consistently.
    extra_legends <- list()

    # 1. Alteration legend (all classes present in ctDNA or tissue)
    if (length(global_classes)) {
      alt_pal <- ctdnaTM:::.oncoprint_alt_palette[global_classes]
      extra_legends[[length(extra_legends) + 1L]] <- ComplexHeatmap::Legend(
        labels    = global_classes,
        title     = "Alteration",
        legend_gp = grid::gpar(fill = unname(alt_pal)),
        title_gp  = grid::gpar(fontsize = 10, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = 9))
    }

    # 2. One legend per annotation column (Dose, RECIST, ...) using the
    # global set of values from ann_at_global (built by
    # .oncoprint_build_legend_state over the full patient_data).
    for (ann in names(ann_col_global)) {
      vals <- ann_at_global[[ann]]
      if (!length(vals)) next
      cols <- ann_col_global[[ann]][vals]
      ann_title <- if (!is.null(annotation_labels) &&
                        ann %in% names(annotation_labels) &&
                        nzchar(annotation_labels[[ann]]))
                     annotation_labels[[ann]] else ann
      extra_legends[[length(extra_legends) + 1L]] <- ComplexHeatmap::Legend(
        labels    = vals,
        title     = ann_title,
        legend_gp = grid::gpar(fill = unname(cols)),
        title_gp  = grid::gpar(fontsize = 10, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = 9))
    }

    # 3. Gene-set legend (captured from .oncoprint_left_anno)
    gsl <- ctdnaTM:::.oncoprint_left_anno(ordered_genes, resolved,
                                            internal_sort_genes)
    if (!is.null(gsl) && !is.null(gsl$legend))
      extra_legends[[length(extra_legends) + 1L]] <- gsl$legend

    # 4. Source legend for the ctDNA / Tissue bar colours
    if (isTRUE(show_freq_bar)) {
      extra_legends[[length(extra_legends) + 1L]] <- ComplexHeatmap::Legend(
        labels    = c("ctDNA", "Tissue"),
        title     = "Source",
        legend_gp = grid::gpar(fill = c("#2C7FB8","#D95F0E")),
        title_gp  = grid::gpar(fontsize = 10, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = 9))
    }

    return(structure(list(plot          = obj,
                          engine        = "complexheatmap",
                          subtitle      = subtitle,
                          caption       = caption,
                          title         = title,
                          extra_legends = extra_legends,
                          ht_gap        = ht_gap_custom),
                     class = c("ctdna_concordance_oncoprint",
                                 "ctdna_oncoprint", "list")))
  }

  assembled <- ctdnaTM:::.oncoprint_assemble_gg(panels, title, subtitle, caption,
                                                  legend_position = "right")
  structure(list(plot     = assembled,
                 engine   = "ggplot2",
                 subtitle = subtitle,
                 caption  = caption,
                 title    = title),
            class = c("ctdna_concordance_oncoprint", "ctdna_oncoprint", "list"))
}


# ---------------------------------------------------------------------------
# Tissue classifier -- Guardant schema pass-through, else Personalis translate
# ---------------------------------------------------------------------------
.co_classify_tissue_local <- function(df) {
  if ("Alteration_class" %in% names(df)) return(df)
  if (all(c("Molecular_consequence","CNV_type") %in% names(df)))
    return(ctdnaTM:::.oncoprint_classify(df))

  vt <- as.character(df$Variant_type)
  ve <- if ("Variant_Effect" %in% names(df))
          tolower(as.character(df$Variant_Effect)) else rep("", nrow(df))
  fc <- if ("Functional_Class" %in% names(df))
          toupper(as.character(df$Functional_Class)) else rep("", nrow(df))

  vt_norm <- vt
  vt_norm[tolower(vt) %in% c("snp","snv","point_mutation")]  <- "SNV"
  vt_norm[tolower(vt) %in% c("insertion","deletion","indel",
                              "del","ins","mnv","mnp","complex")] <- "Indel"
  vt_norm[tolower(vt) %in% c("cnv","copy_number","amplification",
                              "deletion_cn","loss","gain")]  <- "CNV"
  vt_norm[tolower(vt) %in% c("fusion","gene_fusion")]         <- "Fusion"
  vt_norm[tolower(vt) %in% c("lgr","large_rearrangement",
                              "structural_variant","sv")]    <- "LGR"

  mc <- rep("", nrow(df))
  mc[grepl("missense", ve)]                            <- "missense"
  mc[grepl("stop_gained|nonsense", ve)]                <- "nonsense"
  mc[grepl("frameshift", ve)]                          <- "frameshift"
  mc[grepl("splice_donor", ve)]                        <- "splice_donor"
  mc[grepl("splice_acceptor", ve)]                     <- "splice_acceptor"
  mc[grepl("stop_lost", ve)]                           <- "stop_lost"
  mc[grepl("start_lost|initiator_codon", ve)]          <- "start_lost"
  mc[grepl("inframe_ins", ve)]                         <- "inframe_insertion"
  mc[grepl("inframe_del", ve)]                         <- "inframe_deletion"
  mc[grepl("^promoter|regulatory_region_variant", ve)] <- "promoter"
  mc[mc == "" & fc == "MISSENSE"] <- "missense"
  mc[mc == "" & fc == "NONSENSE"] <- "nonsense"

  cn <- if ("CNV_type" %in% names(df)) as.character(df$CNV_type)
        else rep("", nrow(df))

  df$Variant_type          <- vt_norm
  df$Molecular_consequence <- mc
  df$CNV_type              <- cn
  ctdnaTM:::.oncoprint_classify(df)
}


# ---------------------------------------------------------------------------
# print method -- adds the Source legend when drawing the CH plot
# ---------------------------------------------------------------------------
print.ctdna_concordance_oncoprint <- function(x, ...) {
  if (identical(x$engine, "complexheatmap") ||
      identical(x$engine, "ComplexHeatmap")) {
    draw_args <- list(object = x$plot,
                      merge_legend = TRUE,
                      auto_adjust = FALSE,
                      annotation_legend_side = "right",
                      heatmap_legend_side    = "right",
                      padding = grid::unit(c(8, 8, 8, 8), "mm"))
    if (!is.null(x$ht_gap)) draw_args$ht_gap <- x$ht_gap
    if (!is.null(x$extra_legends))
      draw_args$annotation_legend_list <- x$extra_legends
    if (!is.null(x$title) && nzchar(as.character(x$title))) {
      draw_args$column_title    <- x$title
      draw_args$column_title_gp <- grid::gpar(fontsize = 14, fontface = "bold")
    }
    do.call(ComplexHeatmap::draw, draw_args)
  } else {
    print(x$plot)
  }
  invisible(x)
}
