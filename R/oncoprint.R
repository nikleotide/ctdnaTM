# =============================================================================
# v0.39.0 — ctdna_plot_oncoprint() internals
# =============================================================================
#
# Public entry point (ctdna_plot_oncoprint) is defined in ctdnaTM.R. Everything
# else lives here. The pipeline:
#
#   Phase 1: resolve_aliases     -> let label-based group_patients_by / wrap
#                                   refer to columns via annotation_labels.
#   Phase 2: resolve_gene_sets   -> name / vector / inline list / NULL.
#   Phase 3: classify_alterations -> 11-category Alteration_class (new palette).
#   Phase 4: apply_scheme        -> ctdna_stratify_recist() on any RECIST column
#                                   in use as a top annotation or sort key.
#   Phase 5: build_wrap          -> resolve wrap groups (optional).
#   Phase 6: build_global_matrix -> wide gene x patient matrix + gene order.
#   Phase 7: build_legend_state  -> cumulative annotation `at` + colours
#                                   across the WHOLE cohort, so all panels
#                                   share one legend covering every value.
#   Phase 8: render              -> dispatch to CH or ggplot engine; if wrap,
#                                   each engine builds one panel per level
#                                   sharing the legend state from phase 7.

# -----------------------------------------------------------------------------
# Constants — new 11-category alteration palette + RECIST palette/order
# -----------------------------------------------------------------------------

.oncoprint_alt_palette <- c(
  Focal_Amp      = "#FF0000",   # red          — focal_amplification
  Amp            = "#FFC0CB",   # pink         — amplification / aneuploid_amplification
  Homozygous_Del = "#0000FF",   # blue         — homozygous_deletion
  LOH            = "#ADD8E6",   # light blue   — loh_deletion
  Missense       = "#008000",   # green        — missense
  Truncating     = "#FFA500",   # orange       — splice_donor / splice_acceptor /
                                #                nonsense / frameshift /
                                #                stop_lost / start_lost
  InFrame        = "#98FB98",   # pale green   — inframe_insertion / _duplication /
                                #                _deletion / _indel
  Promoter       = "#A28520",   # ochre        — promoter
  LGR            = "#800080",   # purple       — Variant_type = "LGR"
  Fusion         = "#A52A2A",   # brown        — Variant_type = "Fusion"
  Other          = "#224A85"    # dark blue    — anything else (synonymous,
                                #                splice_region, NA, etc.)
)
# Cell-rendering style: thin centred bar for point mutations, full cell for
# structural / copy-number events.
.oncoprint_thin_classes <- c("Missense","Truncating","InFrame","Promoter")

#' Canonical alteration categories used by the oncoprint
#'
#' Returns the 11 canonical alteration-class names in display order
#' (the order the legend uses). Useful as the source of truth when you
#' want to pass an `alterations =` filter to
#' \code{\link{ctdna_plot_oncoprint}} without memorising the strings.
#'
#' @return Character vector of length 11.
#' @examples
#' ctdna_alteration_types()
#' # [1] "Focal_Amp" "Amp" "Homozygous_Del" "LOH" "Missense" "Truncating"
#' #     "InFrame" "Promoter" "LGR" "Fusion" "Other"
#'
#' # Use in ctdna_plot_oncoprint(alterations = ...):
#' sim <- ctdna_make_mock_study(n_patients = 20, seed = 1)
#' keep <- ctdna_alteration_types()
#' keep <- setdiff(keep, c("Other","Promoter"))  # drop two classes
#' ctdna_plot_oncoprint(sim$infinity_report, alterations = keep)
#' @export
ctdna_alteration_types <- function() {
  names(.oncoprint_alt_palette)
}

# v0.42.0: alteration-set resolver used by ctdna_plot_oncoprint().
# Centralised so the canonical list, error messages and
# whitelist/blacklist semantics live in one place.
#
# alterations  NULL or character vector
# flag_dots    list captured from ... by the caller; only names from
#              ctdna_alteration_types() that resolve to TRUE/FALSE are
#              treated as flags
#
# Returns: character vector of alteration types to keep (a subset of
# ctdna_alteration_types(), preserving the canonical order).
.resolve_oncoprint_alterations <- function(alterations, flag_dots) {
  all_types <- names(.oncoprint_alt_palette)

  # Detect the legitimate alteration flags vs unknown ... names
  flag_names <- intersect(names(flag_dots), all_types)
  unknown    <- setdiff(names(flag_dots), all_types)
  if (length(unknown) > 0L)
    stop("Unknown argument(s) passed to ctdna_plot_oncoprint(): ",
         paste(unknown, collapse = ", "),
         ". Valid alteration flag names: ",
         paste(all_types, collapse = ", "), call. = FALSE)

  flags <- flag_dots[flag_names]
  # Drop any flag whose value isn't a single TRUE/FALSE
  flags <- flags[vapply(flags, function(v)
    is.logical(v) && length(v) == 1L && !is.na(v), logical(1))]

  if (!is.null(alterations)) {
    if (!is.character(alterations) || !length(alterations))
      stop("`alterations` must be a non-empty character vector of names from ",
           "ctdna_alteration_types().", call. = FALSE)
    bad <- setdiff(alterations, all_types)
    if (length(bad) > 0L)
      stop("Unknown alteration(s): ", paste(bad, collapse = ", "),
           ". Available: ", paste(all_types, collapse = ", "),
           call. = FALSE)
    if (length(flags) > 0L)
      warning("Both `alterations` and individual flags were passed to ",
              "ctdna_plot_oncoprint(); using `alterations` and ignoring flags.",
              call. = FALSE)
    return(all_types[all_types %in% alterations])  # preserve canonical order
  }

  if (length(flags) == 0L) return(all_types)       # default: show all

  vals      <- vapply(flags, isTRUE, logical(1))
  positives <- names(flags)[vals]
  negatives <- names(flags)[!vals]

  if (length(positives) > 0L && length(negatives) > 0L) {
    warning("Mixed TRUE/FALSE alteration flags passed; using only positives: ",
            paste(positives, collapse = ", "), ".", call. = FALSE)
    return(all_types[all_types %in% positives])
  }
  if (length(positives) > 0L) return(all_types[all_types %in% positives])
  # only negatives -> show all except those
  all_types[!all_types %in% negatives]
}

# RECIST palette covers every notation the scheme collapse can emit AND
# every form ctdna_stratify_recist() can accept as input (raw RECIST,
# slash-form pre-collapsed labels from Response_Subcategory, plus-form
# labels emitted by the scheme).
.oncoprint_recist_pal <- c(
  "CR"="#0000FF","uCR"="#5F00FF","PR"="#0000FF","uPR"="#5F00FF",
  "SD"="#008000","PD"="#FF0000","NE"="#984EA3","NA"="#984EA3",
  "CR/PR"="#0000FF","uCR/PR"="#5F00FF",
  "CR+PR"="#0000FF",
  "PD/NE/NA"="#FF0000","PD+NE+NA"="#FF0000",
  "NE/NA"="#984EA3","NE+NA"="#984EA3",
  "R"="#0000FF","NR"="#FF0000"
)
.oncoprint_recist_order <- c(
  "CR","uCR","PR","uPR","CR/PR","uCR/PR","CR+PR",
  "SD",
  "PD","PD/NE/NA","PD+NE+NA",
  "NE","NE/NA","NE+NA","NA",
  "R","NR"
)

# Categorical fallback for non-RECIST annotation columns.
.oncoprint_set2_palette <- c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3",
                             "#A6D854","#FFD92F","#E5C494","#B3B3B3")

# -----------------------------------------------------------------------------
# Built-in gene-set registry — single source of truth lives in
# rules_library.R as .HRR14_GENES / .TSG_GENES / .RTK_GENES / etc.
# -----------------------------------------------------------------------------
.oncoprint_builtin_gene_sets <- function() {
  list(
    HRR14        = .HRR14_GENES,
    TSG          = .TSG_GENES,
    RTK          = .RTK_GENES,
    Cell_Cycle   = .CELL_CYCLE_GENES,
    TP53_pathway = .TP53_PATHWAY_GENES,
    MMR          = .MMR_GENES,
    PI3K         = .PI3K_GENES
  )
}

# -----------------------------------------------------------------------------
# Is this annotation column a RECIST-flavoured one?
# Triggered by either column name OR value content (>= 70 % of non-NA
# values are in the canonical RECIST palette).
# -----------------------------------------------------------------------------
.oncoprint_is_recist_col <- function(nm, vals) {
  if (toupper(nm) %in% c("RECIST","RECIST_GROUP","BOR","BCR",
                         "RESPONSE","RESPONSE_SUBCATEGORY",
                         "BEST_RESPONSE","AVALC"))
    return(TRUE)
  v <- as.character(vals); v <- v[!is.na(v) & nzchar(v)]
  if (!length(v)) return(FALSE)
  mean(v %in% names(.oncoprint_recist_pal)) >= 0.7
}

# -----------------------------------------------------------------------------
# Phase 1 — resolve an argument that may be either a column name (in
# patient_data or df) OR a display label from annotation_labels.
# -----------------------------------------------------------------------------
.oncoprint_resolve_alias <- function(x, arg_name, patient_data, df,
                                       annotation_labels) {
  if (is.null(x)) return(NULL)
  pd_cols <- if (!is.null(patient_data)) names(patient_data) else character(0)
  df_cols <- names(df)
  if (x %in% pd_cols) return(x)
  if (x %in% df_cols) return(x)
  if (!is.null(annotation_labels)) {
    hit <- names(annotation_labels)[annotation_labels == x]
    if (length(hit)) {
      if (hit[1] %in% pd_cols) return(hit[1])
      if (hit[1] %in% df_cols) return(hit[1])
    }
  }
  stop(sprintf(
    "ctdna_plot_oncoprint: `%s = '%s'` is not a column in patient_data or df, ",
    arg_name, x),
    "and not a known annotation label.\n",
    "  Available columns: ", paste(unique(c(pd_cols, df_cols)), collapse = ", "),
    if (!is.null(annotation_labels))
      paste0("\n  Known labels: ", paste(annotation_labels, collapse = ", "))
    else "",
    call. = FALSE)
}

# -----------------------------------------------------------------------------
# Phase 2 — resolve gene_sets:
#   NULL          -> single "All genes" set from df$Gene
#   character     -> one or more built-in gene-set names
#   named list    -> inline definition
# Returns: list(resolved = <named list of gene vectors>, gene_union)
# -----------------------------------------------------------------------------
.oncoprint_resolve_gene_sets <- function(df, gene_sets, gene_col) {
  builtin <- .oncoprint_builtin_gene_sets()
  resolved <- if (is.null(gene_sets)) {
    g <- df[[gene_col]]
    list(`All genes` = sort(unique(g[!is.na(g) & nzchar(g)])))
  } else if (is.character(gene_sets)) {
    unknown <- setdiff(gene_sets, names(builtin))
    if (length(unknown))
      stop("ctdna_plot_oncoprint: unknown gene-set name(s): ",
           paste(unknown, collapse = ", "),
           ". Built-ins: ", paste(names(builtin), collapse = ", "),
           ". Or pass a named list to define new ones.", call. = FALSE)
    builtin[gene_sets]
  } else if (is.list(gene_sets)) {
    if (is.null(names(gene_sets)) || any(!nzchar(names(gene_sets))))
      stop("ctdna_plot_oncoprint: inline `gene_sets` must be a NAMED list, ",
           "e.g. list(HRR = c('BRCA1','BRCA2')).",
           call. = FALSE)
    gene_sets
  } else {
    stop("ctdna_plot_oncoprint: `gene_sets` must be NULL, a character vector ",
         "of built-in names, or a named list.", call. = FALSE)
  }
  list(resolved = resolved,
       gene_union = unique(unlist(resolved, use.names = FALSE)))
}

# -----------------------------------------------------------------------------
# Phase 3 — classify each row to one of the 11 alteration categories. Uses
# Variant_type, Molecular_consequence and CNV_type per the v0.39.0 spec.
# Returns df with `Alteration_class` added.
# -----------------------------------------------------------------------------
.oncoprint_classify <- function(df) {
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
  df$Alteration_class <- cls
  df
}

# -----------------------------------------------------------------------------
# Phase 4 — apply RECIST scheme collapse via ctdna_stratify_recist() to any
# column in patient_data that is in use (top_annotations or group_patients_by)
# AND detected as RECIST. Returns the possibly-mutated patient_data along
# with the vector of column names that were collapsed.
# -----------------------------------------------------------------------------
.oncoprint_apply_scheme <- function(patient_data, top_annotations,
                                      group_patients_by, scheme) {
  recist_cols_used <- character(0)
  if (is.null(patient_data)) return(list(patient_data = patient_data,
                                          recist_cols_used = recist_cols_used))
  cand <- unique(c(top_annotations, group_patients_by))
  cand <- cand[cand %in% names(patient_data)]
  for (cn_ in cand) {
    if (.oncoprint_is_recist_col(cn_, patient_data[[cn_]])) {
      patient_data[[cn_]] <-
        as.character(ctdna_stratify_recist(patient_data[[cn_]],
                                             recist_grouping = scheme))
      recist_cols_used <- c(recist_cols_used, cn_)
    }
  }
  list(patient_data = patient_data, recist_cols_used = recist_cols_used)
}

# -----------------------------------------------------------------------------
# Phase 5 — resolve wrap groups (one panel per level). The wrap column may
# live in patient_data (preferred) or in df. Levels follow factor levels
# if available; otherwise alphabetical.
# Returns NULL if no wrap requested.
# -----------------------------------------------------------------------------
.oncoprint_resolve_wrap <- function(wrap, patient_data, pd_id_col, df,
                                      patient_col) {
  if (is.null(wrap)) return(NULL)
  if (!is.null(patient_data) && wrap %in% names(patient_data) &&
      !is.na(pd_id_col)) {
    w_pid <- as.character(patient_data[[pd_id_col]])
    w_raw <- patient_data[[wrap]]
    w_vals <- as.character(w_raw)
    lv_factor <- if (is.factor(w_raw)) levels(w_raw) else NULL
  } else if (wrap %in% names(df)) {
    w_pid <- as.character(df[[patient_col]])
    w_raw <- df[[wrap]]
    w_vals <- as.character(w_raw)
    ord <- order(is.na(w_vals)); w_pid <- w_pid[ord]; w_vals <- w_vals[ord]
    dup <- duplicated(w_pid); w_pid <- w_pid[!dup]; w_vals <- w_vals[!dup]
    lv_factor <- if (is.factor(w_raw)) levels(w_raw) else NULL
  } else {
    stop("ctdna_plot_oncoprint: `wrap` column '", wrap,
         "' not found in patient_data or df.", call. = FALSE)
  }
  keep <- !is.na(w_vals) & nzchar(w_vals)
  if (sum(!keep))
    message("ctdna_plot_oncoprint: ", sum(!keep),
            " patient(s) had NA in wrap column; dropped.")
  w_pid <- w_pid[keep]; w_vals <- w_vals[keep]
  levels_out <- if (!is.null(lv_factor)) lv_factor[lv_factor %in% w_vals]
                 else sort(unique(w_vals))
  list(levels = levels_out,
       pat_map = split(w_pid, w_vals)[levels_out])
}

# -----------------------------------------------------------------------------
# Phase 6 — global cohort, gene x patient matrix, gene order, optional
# row_split (when sort_genes = "within_set" and there are multiple sets).
# -----------------------------------------------------------------------------
.oncoprint_build_matrix <- function(df, patient_data, pd_id_col, patient_col,
                                      gene_col, variant_col, gene_union,
                                      resolved, sort_genes,
                                      explicit_cohort = NULL) {
  global_cohort <- if (!is.null(explicit_cohort)) {
    # v0.39.2: user passed show_all_patients = <char vector> — use as cohort
    explicit_cohort
  } else if (!is.null(patient_data) && !is.na(pd_id_col)) {
    unique(as.character(patient_data[[pd_id_col]]))
  } else {
    unique(as.character(df[[patient_col]]))
  }

  pid_df <- as.character(df[[patient_col]])
  all_genes <- intersect(gene_union, unique(df[[gene_col]]))
  global_mat <- matrix("", nrow = length(all_genes), ncol = length(global_cohort),
                       dimnames = list(all_genes, global_cohort))
  key <- paste(df[[gene_col]], pid_df, sep = "\u0001")
  spl <- split(df$Alteration_class, key)
  for (k in names(spl)) {
    parts <- strsplit(k, "\u0001", fixed = TRUE)[[1]]
    if (parts[1] %in% all_genes && parts[2] %in% global_cohort)
      global_mat[parts[1], parts[2]] <- paste(unique(spl[[k]]), collapse = ";")
  }
  global_freq <- rowSums(global_mat != "")
  global_mat  <- global_mat[global_freq > 0, , drop = FALSE]
  global_freq <- global_freq[global_freq > 0]
  if (nrow(global_mat) == 0L)
    stop("ctdna_plot_oncoprint: no alterations to plot after classification.",
         call. = FALSE)

  row_split <- NULL
  ordered_genes <- rownames(global_mat)
  if (identical(sort_genes, "global")) {
    ordered_genes <- ordered_genes[order(-global_freq)]
  } else if (identical(sort_genes, "within_set")) {
    set_of <- rep(NA_character_, length(ordered_genes))
    names(set_of) <- ordered_genes
    for (sn in names(resolved)) {
      hits <- ordered_genes %in% resolved[[sn]] & is.na(set_of)
      set_of[hits] <- sn
    }
    set_of[is.na(set_of)] <- "Unassigned"
    ord <- order(match(set_of, c(names(resolved), "Unassigned")), -global_freq)
    ordered_genes <- ordered_genes[ord]
    if (length(unique(set_of)) > 1L)
      row_split <- factor(set_of[ord],
                          levels = unique(c(names(resolved), "Unassigned")))
  }
  list(global_cohort = global_cohort,
       global_mat   = global_mat[ordered_genes, , drop = FALSE],
       ordered_genes = ordered_genes,
       row_split = row_split)
}

# -----------------------------------------------------------------------------
# Phase 7 — build the *cumulative* legend state for top annotations.
# Done from the whole cohort so wrap panels share one legend covering every
# value present anywhere, not just the values in the first panel.
# -----------------------------------------------------------------------------
.oncoprint_build_legend_state <- function(top_annotations, patient_data,
                                            recist_cols_used) {
  ann_col_global <- list(); ann_at_global <- list()
  if (is.null(top_annotations) || is.null(patient_data))
    return(list(col = ann_col_global, at = ann_at_global))

  for (ann in top_annotations) {
    if (!ann %in% names(patient_data)) {
      warning("ctdna_plot_oncoprint: annotation '", ann,
              "' not in patient_data; skipped.", call. = FALSE); next
    }
    raw <- patient_data[[ann]]
    uvals <- unique(as.character(raw))
    uvals <- uvals[!is.na(uvals) & nzchar(uvals)]
    if (!length(uvals)) next
    if (ann %in% recist_cols_used) {
      ord <- c(intersect(.oncoprint_recist_order, uvals),
               setdiff(uvals, .oncoprint_recist_order))
      cols <- vapply(ord, function(v)
        if (v %in% names(.oncoprint_recist_pal))
          .oncoprint_recist_pal[[v]] else "#9E9E9E",
        character(1))
      ann_col_global[[ann]] <- setNames(cols, ord)
    } else {
      ord <- if (is.factor(raw)) intersect(levels(raw), uvals) else sort(uvals)
      pal <- .oncoprint_set2_palette
      ann_col_global[[ann]] <- setNames(
        if (length(ord) <= length(pal)) pal[seq_along(ord)]
        else grDevices::colorRampPalette(pal)(length(ord)), ord)
    }
    ann_at_global[[ann]] <- ord
  }
  list(col = ann_col_global, at = ann_at_global)
}

# -----------------------------------------------------------------------------
# alter_fun factory for the CH engine: grey background; full cell for
# structural / CNV; thin centred bar for point mutations.
# -----------------------------------------------------------------------------
.oncoprint_alter_fun <- function(classes_present) {
  af <- list(background = function(x, y, w, h)
    grid::grid.rect(x, y,
                    w - grid::unit(0.5,"mm"), h - grid::unit(0.5,"mm"),
                    gp = grid::gpar(fill = "#CCCCCC", col = NA)))
  for (cl in classes_present) {
    local({
      cc <- cl; fill <- .oncoprint_alt_palette[[cc]]
      if (cc %in% .oncoprint_thin_classes) {
        af[[cc]] <<- function(x, y, w, h)
          grid::grid.rect(x, y, w - grid::unit(0.5,"mm"), h * 0.33,
                          gp = grid::gpar(fill = fill, col = NA))
      } else {
        af[[cc]] <<- function(x, y, w, h)
          grid::grid.rect(x, y,
                          w - grid::unit(0.5,"mm"), h - grid::unit(0.5,"mm"),
                          gp = grid::gpar(fill = fill, col = NA))
      }
    })
  }
  af
}

# -----------------------------------------------------------------------------
# Sort patients within a panel. Uses RECIST canonical order when the sort
# column was collapsed via the scheme; factor levels otherwise; alphabetical
# as final fallback. Burden is the tie-breaker.
# -----------------------------------------------------------------------------
.oncoprint_sort_patients <- function(sub, burden, group_patients_by,
                                       patient_data, pd_id_col,
                                       recist_cols_used) {
  if (!is.null(group_patients_by) && !is.null(patient_data) &&
      !is.na(pd_id_col) && group_patients_by %in% names(patient_data)) {
    raw_grp <- patient_data[[group_patients_by]]
    grp <- as.character(raw_grp)[
      match(colnames(sub), as.character(patient_data[[pd_id_col]]))]
    grp[is.na(grp) | !nzchar(grp)] <- "NA"
    level_order <- if (group_patients_by %in% recist_cols_used) {
      u <- unique(grp)
      c(intersect(.oncoprint_recist_order, u),
        setdiff(u, .oncoprint_recist_order))
    } else if (is.factor(raw_grp)) {
      intersect(levels(raw_grp), unique(grp))
    } else sort(unique(grp))
    sub[, order(match(grp, level_order), -burden), drop = FALSE]
  } else {
    sub[, order(-burden), drop = FALSE]
  }
}

# -----------------------------------------------------------------------------
# Build a HeatmapAnnotation for the CH engine — handles `annotation_labels`
# for both strip label and legend title, and uses the cumulative
# `at` / colour state so the legend lists every value present anywhere.
# `show_legend` flips off on all panels except the first so wrap renders
# get one shared legend.
# -----------------------------------------------------------------------------
.oncoprint_top_anno <- function(sub_cols, patient_data, pd_id_col,
                                  ann_col_global, ann_at_global,
                                  annotation_labels, show_legend) {
  if (!length(ann_col_global)) return(NULL)
  idx <- match(sub_cols, as.character(patient_data[[pd_id_col]]))
  ann_args <- list()
  for (ann in names(ann_col_global))
    ann_args[[ann]] <- as.character(patient_data[[ann]])[idx]

  ann_names  <- names(ann_args)
  ann_labels <- ann_names
  if (!is.null(annotation_labels)) {
    rn <- annotation_labels[ann_names]
    ann_labels <- ifelse(is.na(rn) | !nzchar(rn), ann_names, rn)
  }
  legend_param <- setNames(
    lapply(seq_along(ann_names), function(i) list(
      title  = ann_labels[i],
      at     = ann_at_global[[ann_names[i]]],
      labels = ann_at_global[[ann_names[i]]])),
    ann_names)

  do.call(ComplexHeatmap::HeatmapAnnotation,
    c(ann_args,
      list(col                       = ann_col_global[ann_names],
           na_col                    = "#E0E0E0",
           annotation_label          = ann_labels,
           annotation_name_side      = "left",
           annotation_name_gp        = grid::gpar(fontsize = 8),
           simple_anno_size          = grid::unit(4, "mm"),
           gap                       = grid::unit(0.5, "mm"),
           annotation_legend_param   = legend_param,
           show_legend               = show_legend)))
}

# -----------------------------------------------------------------------------
# Build a left-side gene-set membership annotation for the CH engine
# (multi-column display, one column per set, NA where the gene isn't a
# member). Only used when sort_genes is "global" or "none" and there is
# more than one set. In "within_set" mode we rely on row_split bands
# instead — showing both would be redundant.
#
# Returns a list with two slots:
#   $row_anno : a rowAnnotation (per-set legends suppressed)
#   $legend   : a single ComplexHeatmap::Legend titled "Gene sets" that
#               lists every set with its color. The caller stashes this
#               on the ctdna_oncoprint object so print.ctdna_oncoprint
#               can pass it to draw() via annotation_legend_list.
# Returns NULL when the annotation should not be drawn.
# -----------------------------------------------------------------------------
.oncoprint_left_anno <- function(ordered_genes, resolved, sort_genes) {
  set_names <- names(resolved)
  if (length(set_names) < 1L || identical(sort_genes, "within_set"))
    return(NULL)
  # Skip the auto "All genes" placeholder (gene_sets was NULL) — a
  # single-set legend titled "Gene sets" with one entry "All genes" is
  # noise, not information.
  if (length(set_names) == 1L && set_names[1] == "All genes")
    return(NULL)

  short_label <- function(s) {
    lookup <- c("HRR14" = "HRR14", "TSG" = "TSG", "RTK" = "RTK",
                "Cell_Cycle" = "Cell", "TP53_pathway" = "TP53",
                "MMR" = "MMR", "PI3K" = "PI3K")
    if (s %in% names(lookup)) lookup[[s]] else
      substr(gsub("_", " ", s), 1, 8)
  }
  set_pal <- c("#377EB8","#4DAF4A","#77EDDD","#E41A1C",
               "#984EA3","#FF7F00","#A65628","#F781BF",
               "#999999","#1F78B4")
  n <- length(set_names)
  fills <- if (n <= length(set_pal)) set_pal[seq_len(n)]
            else grDevices::colorRampPalette(set_pal)(n)
  names(fills) <- set_names

  ann_df       <- data.frame(row.names = ordered_genes, check.names = FALSE,
                             stringsAsFactors = FALSE)
  col_list_lt  <- list()
  for (s in set_names) {
    ann_df[[s]]        <- ifelse(ordered_genes %in% resolved[[s]], s, NA_character_)
    col_list_lt[[s]]   <- setNames(fills[[s]], s)
  }
  row_anno <- ComplexHeatmap::rowAnnotation(
    df                   = ann_df,
    col                  = col_list_lt,
    na_col               = "transparent",
    show_legend          = rep(FALSE, length(set_names)),   # v0.39.1: unified legend below
    annotation_label     = vapply(set_names, short_label, character(1)),
    annotation_name_side = "bottom",
    annotation_name_rot  = 45,
    annotation_name_gp   = grid::gpar(fontsize = 8),
    simple_anno_size     = grid::unit(3, "mm"))

  legend <- ComplexHeatmap::Legend(
    labels    = set_names,
    title     = "Gene sets",
    legend_gp = grid::gpar(fill = unname(fills)),
    title_gp  = grid::gpar(fontsize = 9, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 8))

  list(row_anno = row_anno, legend = legend)
}

# -----------------------------------------------------------------------------
# Phase 8a — ComplexHeatmap panel
# -----------------------------------------------------------------------------
.oncoprint_panel_ch <- function(global_mat, ordered_genes, row_split,
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

  sub <- .oncoprint_sort_patients(sub, burden, group_patients_by,
                                    patient_data, pd_id_col, recist_cols_used)

  top_anno  <- .oncoprint_top_anno(colnames(sub), patient_data, pd_id_col,
                                     ann_col_global, ann_at_global,
                                     annotation_labels,
                                     show_legend = is_first)
  # v0.39.1: .oncoprint_left_anno now returns list(row_anno, legend) or NULL.
  # The legend is collected at the top level and passed to draw() so all
  # panels share one combined "Gene sets" legend.
  la_info  <- if (is_first) .oncoprint_left_anno(rownames(sub),
                                                    resolved, sort_genes)
               else NULL
  left_anno <- if (!is.null(la_info)) la_info$row_anno else NULL

  right_anno <- if (isTRUE(show_freq_bar))
    ComplexHeatmap::rowAnnotation(
      row_barplot = ComplexHeatmap::anno_oncoprint_barplot(
        axis_param = list(side = "bottom", labels_rot = 0,
                          gp = grid::gpar(fontsize = 7))),
      show_annotation_name = FALSE,
      width = grid::unit(1.6, "cm")) else NULL

  col_width <- grid::unit(min(4, max(2, 200 / max(ncol(sub), 1))), "mm")

  ComplexHeatmap::oncoPrint(
    sub, alter_fun = alter_fun, col = col_vec,
    column_title         = panel_title,
    column_title_gp      = grid::gpar(fontsize = 11, fontface = "bold"),
    row_names_side       = "right",
    pct_side             = "left",
    row_names_gp         = grid::gpar(fontsize = 9),
    column_names_gp      = grid::gpar(fontsize = 6),
    pct_gp               = grid::gpar(fontsize = 8),
    left_annotation      = left_anno,
    right_annotation     = right_anno,
    top_annotation       = top_anno,
    row_order            = seq_len(nrow(sub)),
    column_order         = seq_len(ncol(sub)),
    row_split            = row_split,
    column_gap           = grid::unit(0.5,"mm"),
    row_gap              = if (!is.null(row_split)) grid::unit(3,"mm")
                            else grid::unit(0.5,"mm"),
    width                = col_width * ncol(sub),
    show_column_names    = isTRUE(show_patient_names),
    show_row_names       = is_last,
    show_pct             = TRUE,
    heatmap_legend_param = list(
      title    = "Alteration",
      at       = global_classes,
      labels   = global_classes,
      title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 8)),
    remove_empty_columns = FALSE,
    show_heatmap_legend  = is_first)
}

# -----------------------------------------------------------------------------
# Phase 8b — ggplot fallback panel. Same feature surface as the CH engine
# but rendered as a tile plot. `top_annotations` is rendered as extra
# horizontal strips above the main panel via patchwork-style row binding.
# Per-panel prevalence bar drawn on the right; alteration legend on the
# first panel only.
# -----------------------------------------------------------------------------
.oncoprint_panel_gg <- function(global_mat, ordered_genes, row_split,
                                  panel_pids, panel_title,
                                  is_first, is_last,
                                  patient_data, pd_id_col, patient_col,
                                  gene_col, variant_col,
                                  group_patients_by, recist_cols_used,
                                  ann_col_global, ann_at_global,
                                  annotation_labels,
                                  resolved, sort_genes,
                                  keep_zero_burden_cols, show_patient_names,
                                  show_freq_bar,
                                  global_classes) {
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
  sub <- .oncoprint_sort_patients(sub, burden, group_patients_by,
                                    patient_data, pd_id_col, recist_cols_used)

  # Long form for the main tile plot
  long <- do.call(rbind, lapply(seq_len(ncol(sub)), function(j) {
    pat <- colnames(sub)[j]
    do.call(rbind, lapply(seq_len(nrow(sub)), function(i) {
      cell <- sub[i, j]
      if (!nzchar(cell)) return(NULL)
      cls <- unlist(strsplit(cell, ";", fixed = TRUE))
      cls <- intersect(global_classes, cls)[1]   # take highest-priority
      data.frame(Patient = pat, Gene = rownames(sub)[i],
                  Class = cls, stringsAsFactors = FALSE)
    }))
  }))
  pat_lvls  <- colnames(sub)
  gene_lvls <- rev(rownames(sub))   # ggplot stacks top->bottom
  long$Patient <- factor(long$Patient, levels = pat_lvls)
  long$Gene    <- factor(long$Gene,    levels = gene_lvls)
  long$Class   <- factor(long$Class,   levels = global_classes)

  pal <- .oncoprint_alt_palette[global_classes]
  main <- ggplot2::ggplot(long,
                          ggplot2::aes(x = .data$Patient,
                                       y = .data$Gene,
                                       fill = .data$Class)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = pal, name = "Alteration",
                               drop = FALSE) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ctdna_theme() +
    ggplot2::theme(
      axis.text.x  = if (isTRUE(show_patient_names))
        ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7)
      else ggplot2::element_blank(),
      axis.ticks.x = if (isTRUE(show_patient_names))
        ggplot2::element_line() else ggplot2::element_blank(),
      panel.grid   = ggplot2::element_blank(),
      legend.position = if (is_first) "right" else "none") +
    ggplot2::labs(x = NULL, y = NULL, title = panel_title)

  # In within_set mode add a y-facet for the bands
  if (identical(sort_genes, "within_set") && !is.null(row_split)) {
    set_of <- as.character(row_split)
    names(set_of) <- rownames(sub)
    long$Gene_set <- factor(set_of[as.character(long$Gene)],
                            levels = unique(set_of))
    main <- main +
      ggplot2::facet_grid(rows = ggplot2::vars(.data$Gene_set),
                          scales = "free_y", space = "free_y",
                          switch = "y") +
      ggplot2::theme(strip.text.y.left = ggplot2::element_text(angle = 0),
                     strip.placement   = "outside",
                     panel.spacing.y   = grid::unit(2, "mm"))
  }

  # Annotation strips (above the panel) — one ggplot per annotation,
  # row-binded via patchwork-style facets isn't quite right; we use
  # individual ggplot rows assembled later. Keep it as a list so the
  # caller can decide layout.
  ann_plots <- list()
  if (length(ann_col_global)) {
    idx <- match(colnames(sub), as.character(patient_data[[pd_id_col]]))
    ann_names  <- names(ann_col_global)
    ann_labels <- ann_names
    if (!is.null(annotation_labels)) {
      rn <- annotation_labels[ann_names]
      ann_labels <- ifelse(is.na(rn) | !nzchar(rn), ann_names, rn)
    }
    for (i in seq_along(ann_names)) {
      ann  <- ann_names[i]
      vals <- as.character(patient_data[[ann]])[idx]
      d <- data.frame(Patient = factor(colnames(sub), levels = pat_lvls),
                       Value  = factor(vals, levels = ann_at_global[[ann]]),
                       Row    = ann_labels[i],
                       stringsAsFactors = FALSE)
      p_ann <- ggplot2::ggplot(d,
                                ggplot2::aes(x = .data$Patient, y = 1,
                                             fill = .data$Value)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.3) +
        ggplot2::scale_fill_manual(values = ann_col_global[[ann]],
                                     name  = ann_labels[i],
                                     drop  = FALSE,
                                     na.value = "#E0E0E0") +
        ggplot2::scale_x_discrete(expand = c(0, 0)) +
        ggplot2::scale_y_continuous(expand = c(0, 0),
                                      breaks = 1, labels = ann_labels[i]) +
        ctdna_theme() +
        ggplot2::theme(
          axis.text.x  = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          axis.title   = ggplot2::element_blank(),
          panel.grid   = ggplot2::element_blank(),
          legend.position = if (is_first) "right" else "none") +
        ggplot2::labs(y = NULL, x = NULL)
      ann_plots[[ann]] <- p_ann
    }
  }
  list(main = main, ann_plots = ann_plots,
       n_patients = ncol(sub), pat_lvls = pat_lvls)
}

# -----------------------------------------------------------------------------
# Assemble ggplot panels into one figure. If patchwork is installed, stack
# annotation strips above each panel and arrange panels side by side.
# Otherwise, return the first panel's main plot with a message — keeps the
# function from hard-failing when patchwork isn't there, since patchwork
# is in Suggests (matches v0.33 behaviour).
# -----------------------------------------------------------------------------
.oncoprint_assemble_gg <- function(panels, title, subtitle, caption,
                                     legend_position) {
  if (!length(panels)) stop("ctdna_plot_oncoprint: no ggplot panels.",
                              call. = FALSE)

  if (length(panels) == 1L && !length(panels[[1]]$ann_plots)) {
    p <- panels[[1]]$main +
      ggplot2::labs(title = title, subtitle = subtitle, caption = caption)
    return(p)
  }

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    message("ctdna_plot_oncoprint: patchwork not installed; the ggplot ",
            "engine is returning only the first panel's main tile plot. ",
            "Install patchwork for stacked annotations and wrap layouts.")
    p <- panels[[1]]$main +
      ggplot2::labs(title = title, subtitle = subtitle, caption = caption)
    return(p)
  }

  # Stack each panel's annotation strips above its main plot, then arrange
  # the panels horizontally. Each column's relative width = n_patients.
  cols <- lapply(panels, function(pp) {
    if (length(pp$ann_plots)) {
      strips <- Reduce(`/`, pp$ann_plots)
      strips / pp$main +
        patchwork::plot_layout(heights = c(rep(1, length(pp$ann_plots)),
                                            10))
    } else pp$main
  })
  rel_widths <- vapply(panels, function(pp) pp$n_patients, numeric(1))
  combined <- Reduce(`|`, cols) +
    patchwork::plot_layout(widths = rel_widths, guides = "collect") &
    ggplot2::theme(legend.position = legend_position)
  combined +
    patchwork::plot_annotation(title = title, subtitle = subtitle,
                                caption = caption)
}

# -----------------------------------------------------------------------------
# Auto-detect column names in df.
# -----------------------------------------------------------------------------
.oncoprint_autodetect_cols <- function(df, patient_col, gene_col, variant_col) {
  if (is.null(patient_col))
    patient_col <- intersect(c("Patient_ID", .o("subject"), "subject_id",
                                "SUBJID", "USUBJID"), names(df))[1]
  if (is.null(gene_col))
    gene_col    <- intersect(c("Gene", .o("gene"), "gene"), names(df))[1]
  if (is.null(variant_col))
    variant_col <- intersect(c("Variant_type", .o("alteration"),
                                "alteration_type"), names(df))[1]
  if (is.na(patient_col))
    stop("ctdna_plot_oncoprint: could not find a patient ID column in df. ",
         "Set `patient_col` explicitly.", call. = FALSE)
  if (is.na(gene_col))
    stop("ctdna_plot_oncoprint: could not find a gene column in df. ",
         "Set `gene_col` explicitly.", call. = FALSE)
  if (is.na(variant_col))
    stop("ctdna_plot_oncoprint: could not find a variant-type column in df. ",
         "Set `variant_col` explicitly.", call. = FALSE)
  list(patient_col = patient_col, gene_col = gene_col, variant_col = variant_col)
}
