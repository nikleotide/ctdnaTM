#' Paired ctDNA / tissue concordance oncoprint
#'
#' Extension of \code{\link{ctdna_oncoprint}} that renders paired
#' columns per patient -- one ctDNA sample and one tissue sample side
#' by side -- so alteration concordance between the two modalities can
#' be read directly off the matrix. Signature mirrors
#' \code{ctdna_oncoprint} exactly, plus \code{tissue_df} for the
#' tissue calls and two extra ordering controls (\code{sort_by},
#' \code{show_gene_set_annotation}).
#'
#' @details
#' The function assembles two matrices with identical row (gene) and
#' column (patient) orders -- one from \code{prep$variants} (ctDNA)
#' and one from \code{tissue_df} -- then interleaves the columns as
#' \eqn{[P_1|ctDNA, P_1|Tissue, P_2|ctDNA, P_2|Tissue, \ldots]}. Each
#' pair stays adjacent even when the panel renderer would otherwise
#' re-sort columns by burden; a temporary internal override of the
#' sort keeps pairs together.
#'
#' \strong{Tissue source.} \code{tissue_df} defaults to
#' \code{prep$dnaseq} when not supplied. If both are absent the
#' function errors and points to \code{\link{ctdna_prep_add}(prep,
#' dnaseq = ...)}.
#'
#' \strong{Patient base and rendering behaviour} are inherited from
#' \code{\link{ctdna_oncoprint}} -- including the "no Cohort" naming
#' rule, the Indication / Dose fall-through, and the
#' ComplexHeatmap -> ggplot fallback when the CH engine is
#' unavailable.
#'
#' @param prep A \code{ctdna_prep} object.
#' @param tissue_df Optional tissue-DNA-seq data frame (WES/WGS
#'   alteration calls) providing the second modality of each patient
#'   pair. \code{NULL} (default) uses \code{prep$dnaseq}; errors if
#'   both are missing.
#' @param gene_sets Named list of gene-symbol vectors (one row split
#'   per set) or a built-in name (e.g. \code{"HRR14"}). Required.
#' @param filter_scheme Optional filter-scheme name(s) applied to
#'   \code{prep$variants} before the matrix is built.
#' @param visit Visit to restrict to (default \code{"C1D1"}). Accepts
#'   \code{"Cycle 1 Day 1"} and \code{"baseline"} (case-insensitive).
#' @param visit_col Column in \code{prep$variants} carrying the visit.
#' @param wrap Optional column name to split panels by (\code{"Dose"},
#'   \code{"Indication"}, ...). \code{NULL} (default) = single panel.
#' @param facet Alias for \code{wrap}; both names are accepted on every
#'   function that supports panelling. Names a column to split the output
#'   into separate panels. If both are supplied and they differ,
#'   \code{facet} wins and a message is emitted.
#' @param top_annotations Character vector of \code{patient_data}
#'   columns to show as coloured tracks above the matrix.
#' @param annotation_labels Optional named character vector mapping
#'   column names to display labels (e.g.
#'   \code{c(RECIST = "Response")}).
#' @param group_patients_by Optional column used to sort patients
#'   within each panel.
#' @param scheme RECIST collapsing scheme used when RECIST feeds a
#'   panel-splitting or annotation column. One of \code{"three"}
#'   (default), \code{"two"}, \code{"four"}, or \code{"raw"}.
#' @param alterations Optional character vector of alteration classes
#'   to keep on display (see \code{\link{ctdna_alteration_types}}).
#' @param sort_genes Row-sort strategy: \code{"global"} (default;
#'   sort by prevalence across the whole cohort),
#'   \code{"within_set"} (sort within each row split),
#'   \code{"none"}, or \code{"overlap"} (rank genes by how often
#'   ctDNA and tissue call the same gene in the same patient).
#' @param sort_by Column-sort strategy: \code{"burden"} (default;
#'   patients with the most alterations first) or \code{"overlap"}
#'   (patients whose ctDNA and tissue calls overlap the most, first).
#' @param show_all_patients Logical. \code{TRUE} (default) shows
#'   every patient in the base -- even patients whose two samples
#'   are both all-wildtype in the displayed genes. \code{FALSE}
#'   drops those pairs.
#' @param show_patient_names Logical. Print patient IDs on the column
#'   axis. Default \code{TRUE}.
#' @param show_freq_bar Logical. Draw the per-gene frequency bar on
#'   the right of the matrix. Default \code{TRUE}.
#' @param show_gene_set_annotation Logical. Draw a left-side legend
#'   showing which gene set each row belongs to. Default \code{TRUE}.
#' @param engine Rendering engine hint: \code{"auto"} (ComplexHeatmap
#'   if installed, else ggplot), \code{"complexheatmap"}, or
#'   \code{"ggplot"}.
#' @param title Main title. \code{NULL} = none.
#' @param subtitle Subtitle text. \code{NULL} = none.
#' @param caption Caption below the plot. \code{NULL} = none.
#' @param patient_col Name of the patient-ID column on
#'   \code{prep$variants} (default \code{"Patient_ID"}). Must also
#'   exist on \code{tissue_df} for the join.
#' @param gene_col Name of the gene-symbol column (default
#'   \code{"Gene"}). Same column name is required on
#'   \code{tissue_df}.
#' @param variant_col Optional override for the variant/alteration
#'   column. \code{NULL} = auto-detect.
#' @param indication Optional indication (cohort) restriction. A single string
#'   or a character vector of Indication values, e.g.
#'   \code{indication = "HNSCC"} or \code{indication = c("HNSCC","ES_SCLC")}.
#'   \code{NULL} (the default) uses every indication present in the data.
#'   Values are matched against the canonical \code{Indication} column,
#'   resolved via \code{ctdna_opts("indication")} with fallbacks to
#'   \code{Indication}, \code{indication}, \code{Cancertype} and
#'   \code{Cancer_Type}. This is the modern spelling of the historical
#'   \code{cancer_type=} argument.
#' @param indications Optional character vector of Indication values
#'   (e.g. \code{c("NSCLC","mCRPC")}). When set, restricts the
#'   analysis to patients whose Indication is in this set. Lookup
#'   falls through \code{Indication -> indication -> Cancertype ->
#'   Cancer_Type}. \code{NULL} (default) = no filter.
#'
#' @return A \code{ctdna_oncoprint} object with paired columns, ready
#'   to \code{print()}. Rendering follows the same engine dispatch as
#'   \code{\link{ctdna_oncoprint}}.
#'
#' @seealso \code{\link{ctdna_oncoprint}} (the single-modality
#'   equivalent this function extends);
#'   \code{\link{ctdna_alteration_grid}} (bar-chart view of the same
#'   patient base); \code{\link{ctdna_prep_add}} (attach a tissue
#'   \code{dnaseq} frame after \code{ctdna_prepare()}).
#'
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 40, seed = 1)
#' # Pass dnaseq = sim$dnaseq so prep$dnaseq is populated (the default
#' # tissue source for the concordance oncoprint).
#' prep <- ctdna_prepare(
#'   infinity_report = sim$infinity_report,
#'   adam            = list(adsl = sim$clinical),
#'   dnaseq          = sim$dnaseq,
#'   verbose         = FALSE)
#' prep <- ctdna_sample_qc(prep, verbose = FALSE)
#'
#' # Baseline paired oncoprint: ctDNA and tissue side by side per patient.
#' op <- ctdna_concordance_oncoprint_core(prep,
#'   gene_sets = list(TSG = c("TP53","PTEN","CDKN2A","PIK3CA")))
#' print(op)
#' @export

ctdna_concordance_oncoprint_core <- function(
    prep,
    tissue_df          = NULL,
    gene_sets,
    filter_scheme      = NULL,
    visit              = "C1D1",
    visit_col          = "Visit_name",
    wrap               = NULL,
    facet              = NULL,
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
    indication         = NULL,
    indications        = NULL) {
  # v0.79.0: `wrap` is an alias for `facet` (and vice versa).
  facet <- .ctdna_resolve_facet(facet, wrap, "ctdna_concordance_oncoprint_core")
  wrap  <- facet


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
  indications <- .ctdna_resolve_indication(indication, indications)
  prep <- .ctdna_filter_prep_by_indication(prep, indications)
  engine     <- match.arg(engine)
  scheme     <- match.arg(scheme)
  sort_genes <- match.arg(sort_genes)
  sort_by    <- match.arg(sort_by)
  .reject_cohort_args(top_annotations   = top_annotations,
                                 group_patients_by = group_patients_by,
                                 wrap              = wrap,
                                 annotation_labels = names(annotation_labels))
  if (missing(gene_sets) || is.null(gene_sets))
    stop("gene_sets is required (named list, or built-in name like \"HRR14\").",
         call. = FALSE)

  base <- .build_landscape_base(prep, filter_scheme, visit, visit_col,
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

  df <- .oncoprint_classify(df)
  if (!is.null(alterations))
    df <- df[df$Alteration_class %in% alterations, , drop = FALSE]
  if (!nrow(df))
    stop("no variants remain after alteration filter.", call. = FALSE)

  rgs        <- .oncoprint_resolve_gene_sets(df, gene_sets, gene_col)
  resolved   <- rgs$resolved
  gene_union <- rgs$gene_union

  if (is.null(top_annotations) && "RECIST" %in% names(patient_data) &&
      any(!is.na(patient_data$RECIST)))
    top_annotations <- "RECIST"

  sch <- .oncoprint_apply_scheme(patient_data, top_annotations,
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
  built <- .oncoprint_build_matrix(df, patient_data, pd_id_col,
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
    # Prevalence first: the number of patients altered in EITHER modality, so
    # the most prevalent gene pairs sit at the top. Overlap and burden only
    # break ties.
    new_ord <- order(-gene_uni, -gene_score, -gene_burden, ordered_genes)
    ordered_genes <- ordered_genes[new_ord]
    global_mat    <- global_mat[ordered_genes, , drop = FALSE]
    tissue_mat    <- tissue_mat[ordered_genes, , drop = FALSE]
    if (!is.null(row_split)) row_split <- row_split[new_ord]
  }

  # ==== v0.81.0 RE-PIVOT ====================================================
  # The ctDNA / Tissue pairing lives on the GENE (row) axis, not the patient
  # (column) axis. Rows are "TP53 . ctDNA" / "TP53 . Tissue" interleaved so
  # each gene's pair is adjacent; columns are one per patient.
  #
  # Before v0.81.0 this was the other way round (2 columns per patient, genes
  # as columns after a t()), which produced a 600 x 15 figure on a 300-patient
  # study. Genes-as-rows / patients-as-columns matches ctdna_oncoprint and is
  # the conventional oncoprint shape.
  ctdna_suffix  <- " ctDNA"
  tissue_suffix <- " Tissue"
  paired_rows <- as.vector(rbind(paste0(ordered_genes, ctdna_suffix),
                                   paste0(ordered_genes, tissue_suffix)))
  paired_mat  <- matrix("", nrow = length(paired_rows),
                          ncol = length(patients_orig),
                          dimnames = list(paired_rows, patients_orig))
  paired_mat[seq(1L, nrow(paired_mat), by = 2L), ] <-
    global_mat[ordered_genes, patients_orig]
  paired_mat[seq(2L, nrow(paired_mat), by = 2L), ] <- tissue_mat

  # Columns are plain patient IDs again, so patient_data needs no suffix
  # duplication -- annotation lookup on pd_id_col resolves directly.

  # Row split follows the paired rows: each gene's two rows inherit the gene's
  # gene-set label, so ComplexHeatmap keeps the pair inside one set block.
  if (!is.null(row_split))
    row_split <- rep(as.character(row_split), each = 2L)

  # Replace global_mat + explicit_cohort so panel renderers see the paired set.
  global_mat     <- paired_mat
  paired_genes   <- paired_rows          # row order used downstream
  # explicit_cohort refers to PATIENT columns and is unchanged by the re-pivot.

  # Global alteration-class set now must include tissue-only classes.
  tissue_classes <- setdiff(unique(as.vector(tissue_mat)), "")

  # ==== TEMPORARY OVERRIDE of .oncoprint_sort_patients =====================
  # v0.81.0: after the re-pivot each patient is a SINGLE column, so there is no
  # pair adjacency to preserve here -- the ctDNA / Tissue pair lives on the row
  # axis. Burden and overlap are computed per patient by reading that patient's
  # ctDNA rows and Tissue rows.
  paired_sort <- function(sub, burden, group_patients_by, patient_data,
                            pd_id_col, recist_cols_used) {
    cn <- colnames(sub)
    rn <- rownames(sub)
    ct_rows <- endsWith(rn, ctdna_suffix)
    ti_rows <- endsWith(rn, tissue_suffix)

    # Combined burden across both modalities, per patient column.
    pair_burden <- colSums(sub != "")

    # Overlap per patient (Jaccard: |ctDNA n tissue| / |ctDNA u tissue|),
    # matching gene-for-gene between the two row blocks. 0 when neither
    # modality has an alteration.
    if (identical(sort_by, "overlap")) {
      g_ct <- substr(rn[ct_rows], 1L, nchar(rn[ct_rows]) - nchar(ctdna_suffix))
      g_ti <- substr(rn[ti_rows], 1L, nchar(rn[ti_rows]) - nchar(tissue_suffix))
      shared <- intersect(g_ct, g_ti)
      if (length(shared)) {
        r_ct <- match(paste0(shared, ctdna_suffix),  rn)
        r_ti <- match(paste0(shared, tissue_suffix), rn)
        ca <- sub[r_ct, , drop = FALSE] != ""
        ta <- sub[r_ti, , drop = FALSE] != ""
        n_int   <- colSums(ca & ta)
        n_union <- colSums(ca | ta)
        overlap_per_col <- ifelse(n_union == 0L, 0, n_int / n_union)
      } else {
        overlap_per_col <- rep(0, length(cn))
      }
    } else {
      overlap_per_col <- rep(0, length(cn))
    }

    if (!is.null(group_patients_by) && !is.null(patient_data) &&
        !is.na(pd_id_col) && group_patients_by %in% names(patient_data)) {
      raw_grp <- patient_data[[group_patients_by]]
      grp <- as.character(raw_grp)[
        match(cn, as.character(patient_data[[pd_id_col]]))]
      grp[is.na(grp) | !nzchar(grp)] <- "NA"

      level_order <- if (group_patients_by %in% recist_cols_used) {
        u <- unique(grp)
        c(intersect(.oncoprint_recist_order, u),
          setdiff(u, .oncoprint_recist_order))
      } else if (is.factor(raw_grp)) {
        intersect(levels(raw_grp), unique(grp))
      } else sort(unique(grp))

      if (identical(sort_by, "overlap")) {
        sub[, order(match(grp, level_order),
                      -overlap_per_col, -pair_burden, cn), drop = FALSE]
      } else {
        sub[, order(match(grp, level_order), -pair_burden, cn), drop = FALSE]
      }
    } else {
      # No group split -> memo sort the columns so each gene's mutated patients
      # form one contiguous block, top gene first (the conventional oncoprint
      # waterfall). Row order as displayed sets the gene priority, so
      # `sort_genes` composes with this.
      alt <- sub != ""
      key <- lapply(seq_len(nrow(alt)), function(i) -as.integer(alt[i, ]))
      memo <- do.call(order, c(key, list(-pair_burden, cn)))
      sub[, memo, drop = FALSE]
    }
  }
  old_sort <- .oncoprint_sort_patients
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
    # v0.81.0: global_mat rows are the PAIRED rows ("GENE . ctDNA" /
    # "GENE . Tissue"), not bare gene names, so subset by rownames() rather
    # than by `ordered_genes`.
    rows_use <- intersect(rownames(global_mat),
                           as.vector(rbind(paste0(ordered_genes, ctdna_suffix),
                                            paste0(ordered_genes, tissue_suffix))))
    sub <- global_mat[rows_use,
                       intersect(colnames(global_mat), panel_pids),
                       drop = FALSE]
    burden <- colSums(sub != "")
    if (!isTRUE(keep_zero_burden_cols)) {
      keep_pat <- names(burden)[burden > 0]
      sub <- sub[, keep_pat, drop = FALSE]; burden <- burden[keep_pat]
    }
    if (ncol(sub) == 0L) return(NULL)

    sub <- .oncoprint_sort_patients(
      sub, burden, group_patients_by,
      patient_data, pd_id_col, recist_cols_used)

    # ---- Patient annotations (Dose, RECIST) as a TOP column annotation.
    # v0.81.0: patients are columns after the re-pivot, so per-patient
    # annotations sit above the matrix (they were a left rowAnnotation while
    # patients were rows). Matches ctdna_oncoprint.
    patient_col_anno <- NULL
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
      patient_col_anno <- do.call(ComplexHeatmap::HeatmapAnnotation,
        c(ann_args,
          list(which                = "column",
               col                  = ann_col_global[ann_names],
               na_col               = "#E0E0E0",
               annotation_label     = ann_labels,
               annotation_name_side = "left",
               annotation_name_gp   = grid::gpar(fontsize = 8),
               simple_anno_size     = grid::unit(4, "mm"),
               gap                  = grid::unit(0.5, "mm"),
               show_legend          = FALSE)))
    }

    # ---- TOP annotation: per-gene stacked alteration count SPLIT by
    # source. Two adjacent stacked bars per gene: ctDNA (left) and Tissue
    # (right). Each stacked by alteration type. Each labeled with a %
    # using its own sample count as denominator.
    # ---- LEFT annotation: per-ROW stacked alteration-type barplot.
    # v0.81.0: after the re-pivot each row is a single gene x source, so each
    # row needs ONE stacked bar (previously two sub-bars shared one gene row,
    # which needed hand-rolled grid geometry). anno_barplot handles the
    # stacking, the colours and the count axis natively; direction="reverse"
    # grows the bars leftward, away from the matrix.
    freq_anno <- NULL
    if (isTRUE(show_freq_bar)) {
      types_here <- intersect(global_classes, unique(as.vector(sub)))
      types_here <- types_here[nzchar(types_here)]
      if (length(types_here)) {
        # counts[row, type] -- alterations of each class per row
        counts <- matrix(0L, nrow = nrow(sub), ncol = length(types_here),
                          dimnames = list(rownames(sub), types_here))
        for (ti in seq_along(types_here)) {
          tt <- types_here[ti]
          counts[, ti] <- rowSums(vapply(seq_len(ncol(sub)), function(j)
            grepl(tt, sub[, j], fixed = TRUE), logical(nrow(sub))))
        }
        n_pat  <- ncol(sub)
        n_alt  <- rowSums(counts)
        pct_lb <- ifelse(n_alt > 0L,
                          sprintf("%d%%", round(100 * rowSums(sub != "") / max(n_pat, 1L))),
                          "")
        freq_anno <- ComplexHeatmap::rowAnnotation(
          "#Altered" = ComplexHeatmap::anno_barplot(
            counts,
            which      = "row",
            gp         = grid::gpar(fill = unname(col_vec[types_here]), col = NA),
            border     = FALSE,
            bar_width  = 0.9,
            width      = grid::unit(2.6, "cm"),
            axis       = TRUE,
            axis_param = list(side = "bottom", direction = "reverse",
                               gp = grid::gpar(fontsize = 7))),
          pct = ComplexHeatmap::anno_text(
            pct_lb, gp = grid::gpar(fontsize = 7, fontface = "bold",
                                     col = "#2C7FB8"),
            just = "right", location = grid::unit(1, "npc")),
          show_annotation_name = c("#Altered" = TRUE, pct = FALSE),
          annotation_name_side = "bottom",
          annotation_name_gp   = grid::gpar(fontsize = 8),
          annotation_name_rot  = 0)
      }
    }

    # ---- Gene-set strip as a LEFT rowAnnotation, placed leftmost.
    # v0.81.0: genes are rows now, so the set membership strip is a row
    # annotation. It was lost when the old bottom column annotation was
    # dropped in the re-pivot.
    gset_anno <- NULL
    if (length(resolved)) {
      set_names <- names(resolved)
      base_g <- sub(paste0("(", ctdna_suffix, "|", tissue_suffix, ")$"), "",
                     rownames(sub))
      gs_df <- data.frame(row.names = rownames(sub), check.names = FALSE,
                           stringsAsFactors = FALSE)
      gs_col <- list()
      # Same palette that drives the "Gene sets" legend, so the strip and the
      # legend agree.
      set_pal <- c("#377EB8","#4DAF4A","#77EDDD","#E41A1C",
                   "#984EA3","#FF7F00","#A65628","#F781BF",
                   "#999999","#1F78B4")
      pal <- if (length(set_names) <= length(set_pal))
        set_pal[seq_along(set_names)]
        else grDevices::colorRampPalette(set_pal)(length(set_names))
      for (k in seq_along(set_names)) {
        nm <- set_names[k]
        gs_df[[nm]] <- ifelse(base_g %in% resolved[[nm]], nm, NA_character_)
        gs_col[[nm]] <- setNames(pal[k], nm)
      }
      gset_anno <- ComplexHeatmap::rowAnnotation(
        df                   = gs_df,
        col                  = gs_col,
        na_col               = "transparent",
        show_legend          = rep(FALSE, length(set_names)),
        annotation_name_side = "bottom",
        annotation_name_rot  = 90,
        annotation_name_gp   = grid::gpar(fontsize = 7),
        simple_anno_size     = grid::unit(3, "mm"))
    }

    # ---- Status strips: mutated / not-mutated per patient, one per source.
    # v0.81.0: replaces the per-patient concordance barplot. Placed adjacent
    # immediately above the matrix so concordance is scannable -- both
    # "Mutated" = concordant positive, both "Not mutated" = concordant
    # negative, mismatched = discordant.
    status_anno <- NULL
    {
      rn   <- rownames(sub)
      r_ct <- endsWith(rn, ctdna_suffix)
      r_ti <- endsWith(rn, tissue_suffix)
      # Returned as a FACTOR with both levels always declared. ComplexHeatmap
      # derives simple-annotation legend keys from the values actually present,
      # so on a cohort where every patient is mutated (or none are) the missing
      # level was silently dropped from the legend. Declaring the levels keeps
      # both keys visible regardless of the data.
      st_levels <- c("Mutated", "Not mutated")
      mut <- function(rows) {
        v <- if (!any(rows)) rep("Not mutated", ncol(sub))
             else ifelse(colSums(sub[rows, , drop = FALSE] != "") > 0L,
                          "Mutated", "Not mutated")
        factor(v, levels = st_levels)
      }
      st_cols <- c("Mutated" = "#5FA88C", "Not mutated" = "#1B3A5C")
      status_args <- list()
      status_args[["Status (Tissue)"]] <- mut(r_ti)
      status_args[["Status (ctDNA)"]]  <- mut(r_ct)
      status_anno <- do.call(ComplexHeatmap::HeatmapAnnotation,
        c(status_args,
          list(which                = "column",
               col                  = list("Status (Tissue)" = st_cols,
                                            "Status (ctDNA)"  = st_cols),
               annotation_name_side = "left",
               annotation_name_gp   = grid::gpar(fontsize = 8),
               simple_anno_size     = grid::unit(4, "mm"),
               gap                  = grid::unit(0.5, "mm"),
               annotation_legend_param = list(
                 "Status (Tissue)" = list(title = "Status")),
               show_legend          = c(TRUE, FALSE))))
    }

    # ---- Render. v0.81.0: no transpose -- genes stay as rows (2 per gene:
    # ctDNA + Tissue), patients are columns (1 each). Narrower columns since a
    # 300-patient cohort is now 300 columns wide rather than 15.
    t_sub <- sub      # (2 x genes) x patients
    row_h <- grid::unit(5, "mm")
    col_w <- grid::unit(2.6, "mm")

    ht <- ComplexHeatmap::oncoPrint(
      t_sub, alter_fun = alter_fun, col = col_vec,
      row_title            = panel_title,
      row_title_gp         = grid::gpar(fontsize = 11, fontface = "bold"),
      row_title_side       = "left",
      row_names_side       = "right",
      row_names_gp         = grid::gpar(fontsize = 7),
      column_names_gp      = grid::gpar(fontsize = 5),
      column_names_side    = "bottom",
      # v0.81.0: Status strips sit directly above the matrix, with the
      # requested per-patient annotations (RECIST / Dose) above them. Both are
      # column annotations now that patients are columns.
      top_annotation       = .cc_stack_col_anno(status_anno, patient_col_anno),
      left_annotation      = .cc_stack_col_anno(gset_anno, freq_anno),
      bottom_annotation    = NULL,
      right_annotation     = NULL,
      row_order            = seq_len(nrow(t_sub)),
      column_order         = seq_len(ncol(t_sub)),
      # v0.81.0: gene sets split ROWS now (each gene's ctDNA/Tissue pair stays
      # inside its set); patient grouping splits COLUMNS.
      row_split            = row_split,
      row_gap              = if (!is.null(row_split)) grid::unit(3, "mm")
                              else grid::unit(0.5, "mm"),
      column_gap           = grid::unit(0.5, "mm"),
      width                = col_w * ncol(t_sub),
      heatmap_height       = row_h * nrow(t_sub),
      show_row_names       = TRUE,
      show_column_names    = isTRUE(show_patient_names),
      show_pct             = FALSE,
      remove_empty_columns = FALSE,
      show_heatmap_legend  = FALSE)

    list(heatmap = ht, freq_anno = NULL)
  }
  old_panel_ch <- .oncoprint_panel_ch
  assignInNamespace(".oncoprint_panel_ch", paired_panel_ch, ns = "ctdnaTM")
  on.exit(assignInNamespace(".oncoprint_panel_ch", old_panel_ch,
                              ns = "ctdnaTM"), add = TRUE)

  # ==== VERBATIM FROM ctdna_oncoprint (from legend down) ===================
  leg <- .oncoprint_build_legend_state(top_annotations, patient_data,
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
  alter_fun      <- .oncoprint_alter_fun(global_classes)
  col_vec        <- .oncoprint_alt_palette[global_classes]

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
      panels[[label]] <- do.call(.oncoprint_panel_ch, args)
    } else {
      args$patient_col <- patient_col
      args$gene_col    <- gene_col
      args$variant_col <- variant_col
      panels[[label]]  <- do.call(.oncoprint_panel_gg, args)
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
    # v0.81.0: the gene-name labels and the gene-set strip were %v% stack
    # elements only because genes used to live on the COLUMN axis. Genes are
    # rows now, so oncoPrint handles both natively -- gene x source labels via
    # show_row_names, gene sets via row_split -- and stacking column
    # annotations sized to n_genes against a heatmap sized to n_patients is
    # exactly the `nobs` mismatch that made the vertical list invalid.

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
      alt_pal <- .oncoprint_alt_palette[global_classes]
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
    gsl <- .oncoprint_left_anno(ordered_genes, resolved,
                                            internal_sort_genes)
    if (!is.null(gsl) && !is.null(gsl$legend))
      extra_legends[[length(extra_legends) + 1L]] <- gsl$legend

    # 4. Source legend -- colours match the per-stack % label text
    #    (ctDNA = #2C7FB8, Tissue = #D95F0E). Source is encoded by bar
    #    POSITION within each gene column (left = ctDNA, right = Tissue);
    #    the label suffix names the position so the swatches aren't read
    #    as bar fills (which are coloured by alteration class instead).
    if (isTRUE(show_freq_bar)) {
      extra_legends[[length(extra_legends) + 1L]] <- ComplexHeatmap::Legend(
        labels    = c("ctDNA (row)", "Tissue (row)"),
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

  assembled <- .oncoprint_assemble_gg(panels, title, subtitle, caption,
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
    return(.oncoprint_classify(df))

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
  .oncoprint_classify(df)
}


# ---------------------------------------------------------------------------
# print method -- adds the Source legend when drawing the CH plot
# ---------------------------------------------------------------------------
#' @export
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


# ============================================================================
# .cc_stack_col_anno(...)
# ----------------------------------------------------------------------------
# v0.81.0: combine several column HeatmapAnnotation objects into one, dropping
# NULLs, so the Status strips and the caller's per-patient annotations can be
# passed as a single `top_annotation`. ComplexHeatmap has no public constructor
# for concatenating built annotations, so the underlying vectors are merged via
# c() on the annotation list.
# ============================================================================
.cc_stack_col_anno <- function(...) {
  parts <- Filter(Negate(is.null), list(...))
  if (!length(parts)) return(NULL)
  if (length(parts) == 1L) return(parts[[1]])
  Reduce(function(a, b) c(a, b), parts)
}
