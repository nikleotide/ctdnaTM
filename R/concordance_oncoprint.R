# ============================================================================
# ctdna_concordance_oncoprint() -- ctDNA vs tissue paired oncoprint
# ----------------------------------------------------------------------------
# Transposed oncoprint: genes in columns, patients in rows with two stacked
# tracks per patient (ctDNA on top, tissue underneath). Glyphs and palette
# match ctdna_plot_oncoprint() byte-for-byte. Full annotation surface:
# patient_data / annotations -> clinical strips at left; show_freq_bar -> per-
# gene prevalence at top; alterations -> alteration-type filter; sort_genes,
# indications, subtitle/caption/legend_position.
# ============================================================================

.co_alt_palette <- c(
  Focal_Amp      = "#FF0000",
  Amp            = "#FFC0CB",
  Homozygous_Del = "#0000FF",
  LOH            = "#ADD8E6",
  Missense       = "#008000",
  Truncating     = "#FFA500",
  InFrame        = "#98FB98",
  Promoter       = "#A28520",
  LGR            = "#800080",
  Fusion         = "#A52A2A",
  Other          = "#224A85")

.co_thin_classes <- c("Missense","Truncating","InFrame","Promoter")

.co_class_priority <- c(
  "Focal_Amp","Amp","Homozygous_Del","LOH",
  "Truncating","Missense","InFrame","Promoter",
  "Fusion","LGR","Other")

.co_annot_palettes <- list(
  Dose       = c("6 mg/kg" = "#4C9EFF", "8 mg/kg" = "#FF6B4C"),
  RECIST     = c(CR = "#1B7837", PR = "#5AAE61", SD = "#F1B6DA",
                  PD = "#C51B7D", NE = "#B2ABD2", "NE/NA" = "#B2ABD2"),
  Indication = c(mCRPC = "#66C2A5", NSCLC = "#FC8D62", HNSCC = "#8DA0CB",
                  SCLC = "#E78AC3", CRC = "#A6D854", GBM = "#FFD92F",
                  BRCA = "#E5C494"))

.co_classify_ctdna <- function(df) {
  vt <- as.character(df$Variant_type)
  mc <- as.character(df$Molecular_consequence)
  cn <- as.character(df$CNV_type)
  truncating <- c("splice_donor","splice_acceptor","nonsense","frameshift",
                  "stop_lost","start_lost")
  inframe    <- c("inframe_insertion","inframe_duplication",
                  "inframe_deletion","inframe_indel")
  amp_nf     <- c("aneuploid_amplification","amplification")
  cls <- rep("Other", nrow(df))
  cls[vt == "CNV" & cn == "focal_amplification"] <- "Focal_Amp"
  cls[vt == "CNV" & cn %in% amp_nf]              <- "Amp"
  cls[vt == "CNV" & cn == "homozygous_deletion"] <- "Homozygous_Del"
  cls[vt == "CNV" & cn == "loh_deletion"]        <- "LOH"
  cls[vt == "Fusion"]                             <- "Fusion"
  cls[vt == "LGR"]                                <- "LGR"
  cls[vt %in% c("SNV","Indel") & mc == "missense"]   <- "Missense"
  cls[vt %in% c("SNV","Indel") & mc %in% truncating] <- "Truncating"
  cls[vt %in% c("SNV","Indel") & mc %in% inframe]    <- "InFrame"
  cls[vt %in% c("SNV","Indel") & mc == "promoter"]   <- "Promoter"
  cls
}

.co_classify_tissue <- function(df) {
  if ("Alteration_class" %in% names(df))
    return(as.character(df$Alteration_class))
  vt <- as.character(df$Variant_type)
  cls <- rep("Other", nrow(df))
  cls[vt %in% c("SNV","Indel")] <- "Missense"
  cls[vt == "CNV"]              <- "Amp"
  cls[vt == "Focal_Amp"]        <- "Focal_Amp"
  cls[vt == "Homozygous_Del"]   <- "Homozygous_Del"
  cls[vt == "LOH"]              <- "LOH"
  cls[vt == "Fusion"]           <- "Fusion"
  cls[vt == "LGR"]              <- "LGR"
  for (nm in c("Missense","Truncating","InFrame","Promoter"))
    cls[vt == nm] <- nm
  cls
}

.co_resolve_baseline <- function(visit) {
  if (is.null(visit)) return(NULL)
  v <- as.character(visit)
  is_bl <- tolower(trimws(v)) %in% c("baseline","base","screening")
  if (any(is_bl))
    v <- unique(c(v[!is_bl],
                  "C1D1","Cycle 1 Day 1","Baseline","Screening"))
  v
}

.co_collapse_per_patient <- function(long_df) {
  long_df <- long_df[!is.na(long_df$Gene) & nzchar(long_df$Gene), , drop = FALSE]
  long_df <- long_df[!is.na(long_df$Patient_ID) & nzchar(long_df$Patient_ID), , drop = FALSE]
  long_df$.priority <- match(long_df$Alteration_class, .co_class_priority)
  long_df$.priority[is.na(long_df$.priority)] <- length(.co_class_priority) + 1L
  long_df <- long_df[order(long_df$Patient_ID, long_df$Gene, long_df$.priority), ]
  long_df <- long_df[!duplicated(long_df[, c("Patient_ID","Gene")]), ]
  long_df$.priority <- NULL
  long_df[, c("Patient_ID","Gene","Alteration_class")]
}

#' Concordance oncoprint (ctDNA vs tissue)
#'
#' Transposed oncoprint that pairs each patient's ctDNA calls with their
#' orthogonal tissue-genomic calls. Rows: patients, two per patient
#' (ctDNA on top, tissue underneath). Columns: genes. Cells: alteration
#' glyph with the same palette and thin-bar / full-cell rules as
#' [ctdna_plot_oncoprint()].
#'
#' Full annotation parity with [ctdna_plot_oncoprint()]: pass
#' `annotations = c("Dose","RECIST")` for clinical strips at the left, keep
#' `show_freq_bar = TRUE` for the per-gene prevalence bar at the top, filter
#' alteration types with `alterations`, restrict cohort with `indications`.
#'
#' @param prep A `ctdna_prep` object from [ctdna_prepare()].
#' @param tissue_df data.frame with at least Patient_ID, Gene, Variant_type
#'   (e.g. `prep$dnaseq` or a WES call frame). Optionally an
#'   `Alteration_class` column overrides the built-in classifier.
#' @param gene_sets Named list of gene vectors (e.g.
#'   `list(TSG = c("TP53","RB1"), HRR14 = ctdna_gene_set("HRR14"))`).
#'   `NULL` (default) uses the union of every gene called in either data
#'   source. Set names are used by `sort_genes = "within_set"`.
#' @param filter_scheme Optional variant-filter scheme(s) to apply to
#'   `prep$variants` before analysis (character vector of scheme names).
#'   `NULL` (default) = no filter.
#' @param visit Visit filter for ctDNA. `"baseline"` (default) expands to
#'   `c("C1D1","Cycle 1 Day 1","Baseline","Screening")`. Use `"all"` or a
#'   specific visit name to override. Tissue is not filtered by visit.
#' @param wrap Optional column in `prep$clinical` to facet the plot by
#'   (e.g. `"Indication"`).
#' @param indications Optional character vector of Indication values to
#'   restrict the cohort to (e.g. `c("NSCLC","mCRPC")`). `NULL` (default) =
#'   no filter.
#' @param patient_data Optional per-patient annotation data.frame keyed by
#'   Patient_ID. Overrides the `annotations` arg if both are supplied.
#' @param annotations Character vector of column names in `prep$clinical` to
#'   render as annotation strips beside each patient. Each column gets its
#'   own legend with its column name as the legend title (requires the
#'   `ggnewscale` package).
#' @param show_freq_bar Show the per-gene prevalence bar at the top of the
#'   plot (dodged bars, ctDNA vs Tissue). Default `TRUE`.
#' @param alterations Optional character vector restricting which alteration
#'   classes are drawn. Must be a subset of the 11 canonical classes
#'   (Focal_Amp, Amp, Homozygous_Del, LOH, Missense, Truncating, InFrame,
#'   Promoter, LGR, Fusion, Other). `NULL` (default) = all.
#' @param sort_genes Gene column ordering. `"global"` = by combined
#'   prevalence across ctDNA + tissue (default). `"within_set"` = sort within
#'   each `gene_sets` group. `"none"` = keep input order.
#' @param max_patients_per_panel Cap the number of patients shown per facet,
#'   keeping top-N by combined burden. Default `25`.
#' @param show_patient_names `TRUE` / `FALSE` / `NA`. `NA` (default) = auto
#'   (show if a facet has <=30 patients, hide otherwise).
#' @param sort_by Patient row ordering within each facet. One of
#'   `"combined_burden"` (default), `"ctdna_burden"`, `"tissue_burden"`,
#'   `"patient"` (alphabetical).
#' @param title,subtitle,caption Plot chrome.
#' @param legend_position Passed to `theme(legend.position = ...)`.
#'
#' @return A composed `ggplot` (via patchwork) with class
#'   `"ctdna_concordance_oncoprint"`. Add layers directly with `+`. Auxiliary
#'   data is attached as attributes:
#'   \describe{
#'     \item{`ctdna_matrix`}{Patient x Gene wide matrix of ctDNA
#'       alteration classes (character; "" = wildtype).}
#'     \item{`tissue_matrix`}{Same for tissue.}
#'     \item{`plot_dim`}{Named numeric vector `c(width, height)` recommended
#'       for ggsave.}
#'     \item{`patients`}{Patient IDs kept (in row order).}
#'     \item{`genes`}{Gene names kept (in column order).}
#'     \item{`long`}{Long-format data.frame that fed the main tile layers.}
#'     \item{`n_dropped_by_topN`}{Count of patients trimmed by
#'       `max_patients_per_panel`.}
#'   }
#'
#' @examples
#' \dontrun{
#' sim <- ctdna_make_mock_study(n_patients = 60, seed = 3)
#' sim$clinical$indication <- sim$clinical$Cancertype
#' prep <- ctdna_prepare(infinity_report = sim$infinity_report,
#'                        clinical        = sim$clinical,
#'                        dnaseq          = sim$dnaseq,
#'                        verbose         = FALSE)
#' prep <- ctdna_qc_filter(prep, verbose = FALSE)
#'
#' res <- ctdna_concordance_oncoprint(
#'   prep,
#'   tissue_df   = sim$dnaseq,
#'   gene_sets   = list(TSG = c("TP53","RB1","PTEN"),
#'                       HRR = c("BRCA1","BRCA2","ATM","PALB2")),
#'   visit       = "baseline",
#'   wrap        = "Indication",
#'   indications = c("mCRPC","NSCLC"),
#'   annotations = c("Dose","RECIST"),
#'   title       = "ctDNA vs tissue concordance")
#'
#' # res is a ggplot -- chain layers directly
#' res + ggplot2::theme(axis.text = ggplot2::element_text(size = 6))
#'
#' # Auxiliary data via attributes
#' dims <- attr(res, "plot_dim")
#' ggplot2::ggsave("concordance.png", res,
#'                  width = dims["width"], height = dims["height"],
#'                  dpi = 120, limitsize = FALSE)
#' }
#'
#' @export
ctdna_concordance_oncoprint <- function(
    prep,
    tissue_df,
    gene_sets              = NULL,
    filter_scheme          = NULL,
    visit                  = "baseline",
    wrap                   = NULL,
    indications            = NULL,
    patient_data           = NULL,
    annotations            = NULL,
    show_freq_bar          = TRUE,
    alterations            = NULL,
    sort_genes             = c("global","within_set","none"),
    max_patients_per_panel = 25,
    show_patient_names     = NA,
    sort_by                = c("combined_burden","ctdna_burden",
                                "tissue_burden","patient"),
    title                  = "ctDNA vs tissue concordance",
    subtitle               = NULL,
    caption                = NULL,
    legend_position        = "right") {

  if (!is.list(prep) || is.data.frame(prep))
    stop("prep must be a ctdna_prep object (a list returned by ctdna_prepare()).",
          call. = FALSE)
  if (!"variants" %in% names(prep) || !"clinical" %in% names(prep))
    stop("prep must contain $variants and $clinical slots ",
          "(missing: ", paste(setdiff(c("variants","clinical"), names(prep)),
                                collapse = ", "), ").", call. = FALSE)
  if (!is.data.frame(tissue_df))
    stop("tissue_df must be a data.frame with Patient_ID, Gene, Variant_type.", call. = FALSE)
  need <- c("Patient_ID","Gene","Variant_type")
  miss <- setdiff(need, names(tissue_df))
  if (length(miss))
    stop("tissue_df is missing columns: ", paste(miss, collapse=", "), call. = FALSE)
  sort_by    <- match.arg(sort_by)
  sort_genes <- match.arg(sort_genes)

  # ---- 0. cohort restriction via `indications` -----------------------------
  if (!is.null(indications)) {
    if (!is.data.frame(prep$clinical))
      stop("indications: prep$clinical must be a data.frame.", call. = FALSE)
    ind_col <- if (!is.null(wrap) && wrap %in% names(prep$clinical)) wrap
               else if ("Indication" %in% names(prep$clinical)) "Indication"
               else if ("indication" %in% names(prep$clinical)) "indication"
               else stop("no Indication column in prep$clinical.", call. = FALSE)
    keep_pids <- as.character(prep$clinical$Patient_ID)[
      as.character(prep$clinical[[ind_col]]) %in% as.character(indications)]
    if (!length(keep_pids))
      stop(sprintf("no patients match %s in %s.",
                    paste(indications, collapse=", "), ind_col), call. = FALSE)
    prep$clinical <- prep$clinical[as.character(prep$clinical$Patient_ID) %in% keep_pids, ,
                                    drop = FALSE]
    if (is.data.frame(prep$variants))
      prep$variants <- prep$variants[as.character(prep$variants$Patient_ID) %in% keep_pids, ,
                                      drop = FALSE]
    tissue_df <- tissue_df[as.character(tissue_df$Patient_ID) %in% keep_pids, ,
                            drop = FALSE]
    message(sprintf("ctdna_concordance_oncoprint: restricted to %d indication(s) -> %d patient(s).",
                     length(indications), length(keep_pids)))
  }

  # ---- 1. classify + collapse ctDNA + tissue per patient -------------------
  cvar <- prep$variants
  if (!is.null(filter_scheme)) {
    prep_f <- ctdna_variant_filter(prep, filter_scheme,
                                     apply = TRUE, explain = FALSE,
                                     verbose = FALSE)
    cvar <- prep_f$variants
  }
  if (!is.null(visit) && !identical(tolower(as.character(visit)), "all") &&
      "Visit_name" %in% names(cvar)) {
    keep_visits <- .co_resolve_baseline(visit)
    normalise <- function(x) toupper(gsub("[[:space:]_\\-]+", "", trimws(x)))
    cvar <- cvar[normalise(cvar$Visit_name) %in% normalise(keep_visits), , drop = FALSE]
  }
  cvar$Alteration_class <- .co_classify_ctdna(cvar)
  ctdna_pg <- .co_collapse_per_patient(cvar[, c("Patient_ID","Gene","Alteration_class")])

  tissue_df$Alteration_class <- .co_classify_tissue(tissue_df)
  tissue_pg <- .co_collapse_per_patient(tissue_df[, c("Patient_ID","Gene","Alteration_class")])

  if (!is.null(alterations)) {
    unknown <- setdiff(alterations, names(.co_alt_palette))
    if (length(unknown))
      stop("alterations: unknown class(es): ", paste(unknown, collapse=", "),
            "\n  Valid: ", paste(names(.co_alt_palette), collapse=", "), call. = FALSE)
    ctdna_pg  <- ctdna_pg [ctdna_pg $Alteration_class %in% alterations, , drop = FALSE]
    tissue_pg <- tissue_pg[tissue_pg$Alteration_class %in% alterations, , drop = FALSE]
  }

  # ---- gene universe + sort ------------------------------------------------
  gene_group <- NULL
  if (!is.null(gene_sets)) {
    if (!is.list(gene_sets)) gene_sets <- list(genes = gene_sets)
    gene_universe <- unique(unlist(gene_sets, use.names = FALSE))
    gene_group <- unlist(mapply(function(nm, gs) setNames(rep(nm, length(gs)), gs),
                                 names(gene_sets), gene_sets, SIMPLIFY = FALSE),
                          use.names = TRUE)
  } else {
    gene_universe <- unique(c(ctdna_pg$Gene, tissue_pg$Gene))
  }
  ctdna_pg  <- ctdna_pg [ctdna_pg $Gene %in% gene_universe, , drop = FALSE]
  tissue_pg <- tissue_pg[tissue_pg$Gene %in% gene_universe, , drop = FALSE]

  patients_all <- unique(c(ctdna_pg$Patient_ID, tissue_pg$Patient_ID))
  wrap_lookup <- NULL
  if (!is.null(wrap)) {
    if (!is.data.frame(prep$clinical) || !wrap %in% names(prep$clinical))
      stop("wrap column '", wrap, "' not found in prep$clinical.", call. = FALSE)
    idx <- match(patients_all, as.character(prep$clinical$Patient_ID))
    wrap_vals <- as.character(prep$clinical[[wrap]])[idx]
    ok <- !is.na(wrap_vals) & nzchar(wrap_vals)
    if (sum(!ok))
      message(sprintf("ctdna_concordance_oncoprint: dropped %d patient(s) with NA %s.",
                       sum(!ok), wrap))
    patients_all <- patients_all[ok]
    wrap_lookup  <- setNames(wrap_vals[ok], patients_all)
  }

  cburden <- as.integer(table(factor(ctdna_pg$Patient_ID[ctdna_pg$Patient_ID %in% patients_all],
                                       levels = patients_all)))
  tburden <- as.integer(table(factor(tissue_pg$Patient_ID[tissue_pg$Patient_ID %in% patients_all],
                                       levels = patients_all)))
  names(cburden) <- names(tburden) <- patients_all
  ord_score <- switch(sort_by,
    ctdna_burden    = -cburden,
    tissue_burden   = -tburden,
    combined_burden = -(cburden + tburden),
    patient         = seq_along(patients_all))
  bucket <- if (!is.null(wrap_lookup)) wrap_lookup[patients_all]
            else rep("All", length(patients_all))
  keep_patients <- character(0); patient_bucket <- character(0)
  for (lv in unique(bucket)) {
    idx <- which(bucket == lv)
    o   <- order(ord_score[idx], patients_all[idx])
    idx <- idx[o]
    if (is.finite(max_patients_per_panel) && length(idx) > max_patients_per_panel)
      idx <- idx[seq_len(max_patients_per_panel)]
    keep_patients  <- c(keep_patients,  patients_all[idx])
    patient_bucket <- c(patient_bucket, bucket[idx])
  }
  n_dropped <- length(patients_all) - length(keep_patients)
  if (n_dropped > 0)
    message(sprintf("ctdna_concordance_oncoprint: kept top %d per panel; %d patient(s) not shown.",
                     max_patients_per_panel, n_dropped))
  patients_all <- keep_patients
  if (!is.null(wrap_lookup))
    wrap_lookup <- setNames(patient_bucket, patients_all)
  if (!length(patients_all))
    stop("no patients left to plot.", call. = FALSE)

  cgene_freq <- table(ctdna_pg$Gene[ctdna_pg$Patient_ID %in% patients_all])
  tgene_freq <- table(tissue_pg$Gene[tissue_pg$Patient_ID %in% patients_all])
  combined <- setNames(rep(0L, length(gene_universe)), gene_universe)
  combined[names(cgene_freq)] <- combined[names(cgene_freq)] + as.integer(cgene_freq)
  combined[names(tgene_freq)] <- combined[names(tgene_freq)] + as.integer(tgene_freq)
  if (sort_genes == "global") {
    gene_order <- names(sort(combined, decreasing = TRUE))
  } else if (sort_genes == "within_set" && !is.null(gene_group)) {
    gene_order <- unlist(lapply(names(gene_sets), function(nm) {
      gs <- gene_sets[[nm]]; gs <- gs[gs %in% gene_universe]
      gs[order(-combined[gs], gs)]
    }), use.names = FALSE)
  } else {
    gene_order <- gene_universe
  }

  # ---- long frame + main plot ---------------------------------------------
  row_labels <- as.vector(rbind(paste0(patients_all, " \u00B7 ctDNA"),
                                 paste0(patients_all, " \u00B7 Tissue")))
  ctdna_long  <- ctdna_pg [ctdna_pg $Patient_ID %in% patients_all &
                             ctdna_pg $Gene       %in% gene_order, , drop = FALSE]
  ctdna_long$Row <- if (nrow(ctdna_long))
                      paste0(ctdna_long$Patient_ID, " \u00B7 ctDNA") else character(0)
  tissue_long <- tissue_pg[tissue_pg$Patient_ID %in% patients_all &
                             tissue_pg$Gene      %in% gene_order, , drop = FALSE]
  tissue_long$Row <- if (nrow(tissue_long))
                       paste0(tissue_long$Patient_ID, " \u00B7 Tissue") else character(0)
  long <- rbind(ctdna_long, tissue_long)
  long$Row  <- factor(long$Row,  levels = rev(row_labels))
  long$Gene <- factor(long$Gene, levels = gene_order)
  present_classes <- intersect(names(.co_alt_palette),
                                unique(long$Alteration_class))
  long$Alteration_class <- factor(long$Alteration_class,
                                    levels = names(.co_alt_palette))

  bg <- expand.grid(Row = factor(rev(row_labels), levels = rev(row_labels)),
                     Gene = factor(gene_order, levels = gene_order),
                     stringsAsFactors = FALSE)
  if (!is.null(wrap_lookup)) {
    row_to_pid <- setNames(rep(patients_all, each = 2),
                           as.vector(rbind(paste0(patients_all, " \u00B7 ctDNA"),
                                            paste0(patients_all, " \u00B7 Tissue"))))
    bg$wrap_group   <- wrap_lookup[row_to_pid[as.character(bg$Row)]]
    long$wrap_group <- wrap_lookup[row_to_pid[as.character(long$Row)]]
  }
  full <- long[!long$Alteration_class %in% .co_thin_classes, , drop = FALSE]
  thin <- long[ long$Alteration_class %in% .co_thin_classes, , drop = FALSE]

  n_per_panel <- if (!is.null(wrap_lookup))
    max(table(wrap_lookup[patients_all])) else length(patients_all)
  if (is.na(show_patient_names))
    show_patient_names <- (n_per_panel <= 30)
  if (!length(present_classes)) present_classes <- names(.co_alt_palette)
  pal_used <- .co_alt_palette[present_classes]

  main <- ggplot2::ggplot(bg, ggplot2::aes(x = .data$Gene, y = .data$Row)) +
    ggplot2::geom_tile(fill = "#CCCCCC", colour = "white", linewidth = 0.3)
  if (nrow(full))
    main <- main + ggplot2::geom_tile(data = full,
                                        ggplot2::aes(fill = .data$Alteration_class),
                                        colour = "white", linewidth = 0.3)
  if (nrow(thin))
    main <- main + ggplot2::geom_tile(data = thin,
                                        ggplot2::aes(fill = .data$Alteration_class),
                                        colour = "white", linewidth = 0.3,
                                        height = 0.33)
  main <- main +
    ggplot2::scale_fill_manual(values = pal_used, drop = TRUE,
                                name = "Alteration",
                                limits = present_classes) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid    = ggplot2::element_blank(),
      panel.spacing = grid::unit(0.6, "lines"),
      axis.text     = ggplot2::element_text(size = 8),
      axis.text.x   = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y   = if (isTRUE(show_patient_names))
                        ggplot2::element_text() else ggplot2::element_blank(),
      axis.ticks    = ggplot2::element_blank(),
      strip.text.y  = ggplot2::element_text(face = "bold", angle = 0, hjust = 0),
      legend.position = legend_position,
      legend.key.size = grid::unit(0.4, "cm"))
  if (!is.null(wrap_lookup))
    main <- main + ggplot2::facet_grid(wrap_group ~ ., scales = "free_y",
                                         space = "free_y", switch = "y") +
              ggplot2::theme(strip.placement = "outside")

  # ---- annotation strip (left) --------------------------------------------
  ann_plot <- NULL; ann_cols <- character(0)
  if (is.null(patient_data) && !is.null(annotations) &&
      is.data.frame(prep$clinical)) {
    unknown <- setdiff(annotations, names(prep$clinical))
    if (length(unknown))
      warning("annotations: not in prep$clinical, dropping: ",
               paste(unknown, collapse=", "), call. = FALSE)
    annotations <- intersect(annotations, names(prep$clinical))
    if (length(annotations)) {
      pd_idx <- match(patients_all, as.character(prep$clinical$Patient_ID))
      patient_data <- data.frame(Patient_ID = patients_all,
                                   stringsAsFactors = FALSE)
      for (a in annotations)
        patient_data[[a]] <- prep$clinical[[a]][pd_idx]
    }
  }
  if (!is.null(patient_data) && is.data.frame(patient_data) &&
      nrow(patient_data)) {
    pd <- patient_data[match(patients_all, as.character(patient_data$Patient_ID)), , drop = FALSE]
    ann_cols <- setdiff(names(pd), "Patient_ID")
    if (length(ann_cols)) {
      ann_long <- do.call(rbind, lapply(ann_cols, function(a) {
        data.frame(
          Row   = as.vector(rbind(paste0(patients_all, " \u00B7 ctDNA"),
                                    paste0(patients_all, " \u00B7 Tissue"))),
          Ann   = a,
          Value = rep(as.character(pd[[a]]), each = 2),
          Patient_ID = rep(patients_all, each = 2),
          stringsAsFactors = FALSE)
      }))
      ann_long$Row <- factor(ann_long$Row, levels = rev(row_labels))
      ann_long$Ann <- factor(ann_long$Ann, levels = ann_cols)
      if (!is.null(wrap_lookup))
        ann_long$wrap_group <- wrap_lookup[ann_long$Patient_ID]

      .build_pal <- function(a, values) {
        v <- unique(values); v <- v[!is.na(v) & nzchar(v)]
        if (!length(v)) return(character(0))
        builtin <- if (a %in% names(.co_annot_palettes))
                     .co_annot_palettes[[a]] else character(0)
        pal <- builtin[intersect(names(builtin), v)]
        missing <- setdiff(v, names(pal))
        if (length(missing)) {
          extra <- if (requireNamespace("RColorBrewer", quietly = TRUE) &&
                        length(missing) <= 12) {
            cols <- RColorBrewer::brewer.pal(max(3, length(missing)), "Set3")
            setNames(cols[seq_along(missing)], missing)
          } else {
            setNames(grDevices::hcl.colors(length(missing), "Set 3"), missing)
          }
          pal <- c(pal, extra)
        }
        pal
      }

      if (!requireNamespace("ggnewscale", quietly = TRUE))
        stop("ctdna_concordance_oncoprint: annotations require `ggnewscale`. ",
              "install.packages(\"ggnewscale\")", call. = FALSE)

      ann_plot <- ggplot2::ggplot(mapping = ggplot2::aes(x = .data$Ann, y = .data$Row))
      for (i in seq_along(ann_cols)) {
        a   <- ann_cols[i]
        sub <- ann_long[ann_long$Ann == a, , drop = FALSE]
        if (!nrow(sub)) next
        pal_a <- .build_pal(a, sub$Value)
        ann_plot <- ann_plot +
          ggplot2::geom_tile(data = sub,
                               ggplot2::aes(fill = .data$Value),
                               colour = "white", linewidth = 0.3) +
          ggplot2::scale_fill_manual(values = pal_a,
                                       na.value = "#EEEEEE",
                                       name = a, drop = TRUE)
        if (i < length(ann_cols))
          ann_plot <- ann_plot + ggnewscale::new_scale_fill()
      }
      ann_plot <- ann_plot +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::theme(
          panel.grid    = ggplot2::element_blank(),
          panel.spacing = grid::unit(0.6, "lines"),
          axis.text.x   = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y   = if (isTRUE(show_patient_names))
                             ggplot2::element_text(size = 7) else ggplot2::element_blank(),
          axis.ticks    = ggplot2::element_blank(),
          strip.text    = ggplot2::element_blank(),
          legend.position = legend_position)
      if (!is.null(wrap_lookup))
        ann_plot <- ann_plot + ggplot2::facet_grid(wrap_group ~ .,
                                                     scales = "free_y",
                                                     space  = "free_y")
    }
  }

  # ---- per-gene frequency bar (top) ---------------------------------------
  freq_plot <- NULL
  if (isTRUE(show_freq_bar)) {
    n_pt <- length(patients_all)
    freq_df <- rbind(
      data.frame(Gene = gene_order, Source = "ctDNA",
                  Freq = as.integer(cgene_freq[gene_order]) / n_pt,
                  stringsAsFactors = FALSE),
      data.frame(Gene = gene_order, Source = "Tissue",
                  Freq = as.integer(tgene_freq[gene_order]) / n_pt,
                  stringsAsFactors = FALSE))
    freq_df$Freq[is.na(freq_df$Freq)] <- 0
    freq_df$Gene   <- factor(freq_df$Gene,   levels = gene_order)
    freq_df$Source <- factor(freq_df$Source, levels = c("ctDNA","Tissue"))
    freq_plot <- ggplot2::ggplot(freq_df,
                                    ggplot2::aes(x = .data$Gene,
                                                  y = .data$Freq,
                                                  fill = .data$Source)) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
      ggplot2::scale_fill_manual(values = c(ctDNA = "#2C7FB8", Tissue = "#D95F0E"),
                                   name = NULL) +
      ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                    expand = ggplot2::expansion(mult = c(0, 0.05))) +
      ggplot2::labs(x = NULL, y = "% altered") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),
        axis.text.x        = ggplot2::element_blank(),
        axis.ticks.x       = ggplot2::element_blank(),
        axis.text.y        = ggplot2::element_text(size = 8),
        legend.position    = legend_position)
  }

  # ---- compose with patchwork ---------------------------------------------
  n_ann <- length(ann_cols)
  ann_width <- if (n_ann > 0) max(0.10, 0.06 * n_ann) else 0
  main_composable <- if (!is.null(ann_plot))
    main + ggplot2::theme(axis.text.y = ggplot2::element_blank()) else main

  if (!is.null(freq_plot) && !is.null(ann_plot)) {
    p <- ann_plot + freq_plot + main_composable +
      patchwork::plot_layout(design  = "#B\nAC",
                              widths  = c(ann_width, 1),
                              heights = c(0.15, 1),
                              guides  = "collect") +
      patchwork::plot_annotation(title = title, subtitle = subtitle, caption = caption,
        theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold")))
  } else if (!is.null(freq_plot)) {
    p <- freq_plot / main +
      patchwork::plot_layout(heights = c(0.15, 1), guides = "collect") +
      patchwork::plot_annotation(title = title, subtitle = subtitle, caption = caption,
        theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold")))
  } else if (!is.null(ann_plot)) {
    p <- ann_plot + main_composable +
      patchwork::plot_layout(widths = c(ann_width, 1), guides = "collect") +
      patchwork::plot_annotation(title = title, subtitle = subtitle, caption = caption,
        theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold")))
  } else {
    p <- main +
      ggplot2::labs(title = title, subtitle = subtitle, caption = caption) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  }

  # ---- return -------------------------------------------------------------
  ctdna_mat <- {
    m <- matrix("", nrow = length(patients_all), ncol = length(gene_order),
                 dimnames = list(patients_all, gene_order))
    c_kept <- ctdna_pg[ctdna_pg$Patient_ID %in% patients_all &
                        ctdna_pg$Gene       %in% gene_order, , drop = FALSE]
    if (nrow(c_kept))
      m[cbind(match(c_kept$Patient_ID, patients_all),
              match(c_kept$Gene,       gene_order))] <-
        as.character(c_kept$Alteration_class)
    m
  }
  tissue_mat <- {
    m <- matrix("", nrow = length(patients_all), ncol = length(gene_order),
                 dimnames = list(patients_all, gene_order))
    t_kept <- tissue_pg[tissue_pg$Patient_ID %in% patients_all &
                         tissue_pg$Gene       %in% gene_order, , drop = FALSE]
    if (nrow(t_kept))
      m[cbind(match(t_kept$Patient_ID, patients_all),
              match(t_kept$Gene,       gene_order))] <-
        as.character(t_kept$Alteration_class)
    m
  }

  n_cols <- length(gene_order); n_rows_total <- 2 * length(patients_all)
  plot_dim <- c(width  = max(6, 0.16 * n_cols + 3),
                height = max(4, 0.16 * n_rows_total + 2 +
                                 ifelse(show_freq_bar, 1.2, 0)))

  attr(p, "ctdna_matrix")      <- ctdna_mat
  attr(p, "tissue_matrix")     <- tissue_mat
  attr(p, "plot_dim")          <- plot_dim
  attr(p, "patients")          <- patients_all
  attr(p, "genes")             <- gene_order
  attr(p, "long")              <- long
  attr(p, "n_dropped_by_topN") <- n_dropped
  class(p) <- unique(c("ctdna_concordance_oncoprint", class(p)))
  p
}


#' 50 pan-cancer genes (frequently mutated across TCGA/COSMIC)
#'
#' A character vector of 50 gene symbols commonly used as a pan-cancer landscape
#' set. Composition: solid-tumor drivers (TP53, PIK3CA, KRAS, APC, PTEN, BRAF,
#' etc.), DDR / HRR genes (BRCA1/2, PALB2, ATM, MSH2, MLH1), heme drivers
#' (DNMT3A, TET2, IDH1/2, NPM1, FLT3, JAK2), RTKs / fusion targets (EGFR,
#' ERBB2, MET, ALK, ROS1, RET), chromatin regulators (KMT2C, KMT2D, ARID1A,
#' SMARCA4), and indication-specific standouts (VHL, NOTCH1, STK11, KEAP1).
#'
#' @format Character vector of length 50.
#' @examples
#' head(cancer_50, 10)
#' @export
cancer_50 <- c(
  "TP53","PIK3CA","KRAS","APC","PTEN",
  "BRAF","ARID1A","KMT2D","EGFR","NRAS",
  "SETD2","FAT1","NF1","RB1","ATM",
  "CTNNB1","SMAD4","IDH1","FBXW7","CDKN2A",
  "NOTCH1","PIK3R1","STK11","KEAP1","BRCA2",
  "BRCA1","KMT2C","NFE2L2","CREBBP","SPOP",
  "VHL","HRAS","ERBB2","MET","GNAS",
  "FLT3","DNMT3A","TET2","IDH2","NPM1",
  "JAK2","PALB2","MYC","MYCN","ALK",
  "ROS1","RET","MSH2","MLH1","SMARCA4")
