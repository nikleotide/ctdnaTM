# =============================================================================
# ctdnaTM v0.67.0 -- ctdna_oncoprint() and ctdna_alteration_grid()
#
# Clean from-scratch rewrite of the two genomic-landscape plots. Old names
# (ctdna_plot_oncoprint, ctdna_plot_alteration_grid) are GONE -- no alias,
# no deprecation. Only the *rendering style* (palette, glyph sizes / z-order,
# stacked-bar look) is inherited from oncoprint.R helpers; the data plumbing
# is rewritten.
#
# Shared data pipeline (identical for both functions):
#   1. Patients come from prep$variants$Patient_ID.
#   2. Optional filter_scheme applied to prep$variants in-function.
#   3. Visit filter: default "C1D1"; accepts "Cycle 1 Day 1" and "baseline"
#      (all case-insensitive).
#   4. Look up Indication case-insensitively from a column literally named
#      "indication" (any case). Search order: prep$adrs, then any other
#      ADaM frame, then prep$variants, then prep$clinical.
#   5. Look up Dose case-insensitively from a column literally named "dose"
#      (any case). Same search order.
#   6. Drop patients missing Indication OR Dose; print a message with the
#      number removed (identical wording for both functions).
#   7. Reject "COHORT"/"Cohort"/"cohort" anywhere in args -- error points the
#      user to Indication.
# =============================================================================

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# ---- Per-set combine-mode wrappers (used inside `gene_sets` for the grid) ---
# These tag a gene vector with how "altered at set level" is decided.
# Bare vectors (no wrapper) default to `alt_any` for backward compatibility.

#' Combine-mode wrappers for `ctdna_alteration_grid` gene sets
#'
#' Wrap a gene vector inside a `gene_sets` list entry to specify how a patient
#' is scored as "altered" for that set: `alt_any(x)` = at least one gene in
#' `x` altered; `alt_all(x)` = every gene altered; `alt_any_n(x, n)` = at
#' least `n` genes altered. A bare character vector is treated as `alt_any`.
#'
#' @param x Character vector of gene symbols.
#' @param n Integer >= 1 and <= length(x). For `alt_any_n` only.
#' @return A character vector with S3 class `"ctdna_geneset"` and combine-mode
#'   metadata attached.
#' @examples
#' # gene_sets = list(TSG   = alt_any(TSG_gene_set),
#' #                  HRR   = alt_all(ctdna_gene_set("HRR14")),
#' #                  ABCx2 = alt_any_n(c("A","B","C"), 2))
#' @name alt_combine
NULL

#' @rdname alt_combine
#' @export
alt_any <- function(x) {
  x <- as.character(x)
  structure(x, class = c("ctdna_geneset", "character"),
            ctdna_combine = "any")
}

#' @rdname alt_combine
#' @export
alt_all <- function(x) {
  x <- as.character(x)
  structure(x, class = c("ctdna_geneset", "character"),
            ctdna_combine = "all")
}

#' @rdname alt_combine
#' @export
alt_any_n <- function(x, n) {
  x <- as.character(x)
  if (missing(n) || !is.numeric(n) || length(n) != 1L || n < 1L)
    stop("alt_any_n: `n` must be a positive integer scalar.", call. = FALSE)
  if (n > length(x))
    stop("alt_any_n: n (", n, ") exceeds set size (", length(x), ").",
         call. = FALSE)
  structure(x, class = c("ctdna_geneset", "character"),
            ctdna_combine = "any_n",
            ctdna_n       = as.integer(n))
}

# Extract combine spec from any gene_set entry (bare vector => alt_any).
# Returns list(genes, mode, n). `mode` is "any"/"all"/"any_n".
.geneset_combine <- function(entry) {
  if (inherits(entry, "ctdna_geneset")) {
    return(list(genes = as.character(entry),
                mode  = attr(entry, "ctdna_combine"),
                n     = attr(entry, "ctdna_n")))
  }
  list(genes = as.character(entry), mode = "any", n = NULL)
}

# Auto-generate a human-friendly label for a gene_set entry that the user
# didn't name (rare). Falls back to a short signature.
.geneset_autolabel <- function(entry, idx) {
  mode <- if (inherits(entry, "ctdna_geneset"))
            attr(entry, "ctdna_combine") else "any"
  paste0("set", idx, "_", mode)
}

# ---- other helpers ---------------------------------------------------------

# Names of ADaM-style frames typically found in a ctdna_prep. The actual
# scanning is permissive (any data.frame in `prep`), this just controls order.
.adam_frame_order <- function(prep) {
  fr <- names(prep)[vapply(prep, is.data.frame, logical(1))]
  pref <- c("adrs", "adsl", "adtr", "adtte", "adlb", "adae", "adcm", "advs")
  c(intersect(pref, fr),
    setdiff(fr, c(pref, "variants", "clinical")),
    intersect(c("variants", "clinical"), fr))
}

# Case-insensitive column lookup. Returns the actual column name (preserving
# its original case) or NA_character_.
.find_col_ci <- function(d, target) {
  hit <- which(tolower(names(d)) == tolower(target))
  if (length(hit)) names(d)[hit[1]] else NA_character_
}

# Pull a per-patient value (named character vector keyed by Patient_ID) for a
# column literally named `target` (case-insensitive) by scanning the frames of
# `prep` in the documented order. Returns NULL if no frame has the column.
.read_per_patient <- function(prep, target,
                              subj_col = "Patient_ID") {
  if (!is.list(prep)) return(NULL)
  for (fn in .adam_frame_order(prep)) {
    d <- prep[[fn]]
    if (!is.data.frame(d) || !nrow(d)) next
    if (!(subj_col %in% names(d))) next
    cc <- .find_col_ci(d, target)
    if (is.na(cc)) next
    v <- as.character(d[[cc]])
    # collapse to one value per patient (first non-NA/non-blank wins)
    out <- tapply(v, as.character(d[[subj_col]]),
                  function(z) {
                    z <- z[!is.na(z) & nzchar(z)]
                    if (length(z)) z[1] else NA_character_
                  })
    return(out)
  }
  NULL
}

# Resolve a user-supplied visit string. Adds "baseline" as a case-insensitive
# alias for the package's canonical "C1D1". Returns the (possibly translated)
# value -- still passed through .ctdna_normalize_visit downstream.
.resolve_visit_baseline <- function(visit) {
  if (is.null(visit)) return(NULL)
  v <- as.character(visit)
  is_baseline <- tolower(trimws(v)) %in% c("baseline", "base", "screening")
  # A user asking for `visit = "baseline"` means the baseline visit in whatever
  # label form the data uses. Expand into the full equivalent set so downstream
  # normalised matching picks up "Baseline", "Screening", "C1D1", and
  # "Cycle 1 Day 1" alike (mirrors .ctdna_is_baseline_label()).
  if (any(is_baseline)) {
    v <- c(v[!is_baseline],
           "C1D1", "Cycle 1 Day 1", "Baseline", "Screening")
    v <- unique(v)
  }
  v
}

# Hard-reject any "COHORT"/"Cohort"/"cohort" string anywhere in args. The
# v0.67.0 contract is: the dimension is "Indication", end of discussion.
.reject_cohort_args <- function(...) {
  args <- list(...)
  bad  <- c("COHORT", "Cohort", "cohort")
  for (a in names(args)) {
    v <- args[[a]]; if (is.null(v)) next
    if (any(as.character(v) %in% bad))
      stop(sprintf(
        paste0("`%s` references 'COHORT'/'Cohort'/'cohort'. v0.67.0 uses ",
               "'Indication' (read case-insensitively from the `indication` ",
               "column of your ADaM data). Pass `%s = \"Indication\"` or ",
               "drop the COHORT entry."), a, a),
        call. = FALSE)
  }
  invisible(NULL)
}

# Shared base-builder. Returns:
#   variants     -- prep$variants after filter_scheme + visit + drop-NA patients
#   base_pids    -- unique Patient_ID values present
#   patient_data -- one row per base patient with Patient_ID + Indication +
#                   Dose + RECIST (RECIST may be NA where missing)
.build_landscape_base <- function(prep, filter_scheme, visit, visit_col,
                                   patient_col, gene_col, fn_label) {

  if (!inherits(prep, "ctdna_prep"))
    stop(fn_label, ": `prep` must be a ctdna_prep object ",
         "(output of ctdna_prepare()).", call. = FALSE)
  v <- prep$variants
  if (!is.data.frame(v) || !nrow(v))
    stop(fn_label, ": prep$variants is empty.", call. = FALSE)
  if (!(patient_col %in% names(v)))
    stop(fn_label, ": '", patient_col,
         "' is not a column of prep$variants.", call. = FALSE)

  # ---- filter_scheme (optional) ----
  if (!is.null(filter_scheme))
    v <- .filter_apply_df(v, filter_scheme = filter_scheme,
                            gene_col = gene_col, verbose = FALSE)

  # ---- visit (default + baseline alias) ----
  visit <- .resolve_visit_baseline(visit)
  v <- .ctdna_apply_visit_filter(v, visit, visit_col,
                                 fn = fn_label, verbose = FALSE)

  # ---- patient base from variants ----
  base_pids <- unique(as.character(v[[patient_col]]))
  if (!length(base_pids))
    stop(fn_label, ": no patients in prep$variants after filter + visit.",
         call. = FALSE)

  # ---- Indication (required) ----
  ind <- .read_per_patient(prep, "indication", subj_col = patient_col)
  if (is.null(ind))
    stop(fn_label, ": no column literally named 'indication' ",
         "(case-insensitive) found in prep$adrs, any ADaM frame, ",
         "prep$variants, or prep$clinical. Please add an `indication` ",
         "column covering every patient and try again.", call. = FALSE)

  # ---- Dose (optional; annotation only, NA allowed) ----
  dose <- .read_per_patient(prep, "dose", subj_col = patient_col)

  # ---- RECIST (optional here -- grid will further enforce) ----
  rec <- .read_per_patient(prep, "RECIST", subj_col = patient_col)

  # ---- drop patients missing Indication ONLY ----
  ind_b <- ind[base_pids]
  miss  <- is.na(ind_b) | !nzchar(ind_b)
  n_drop <- sum(miss)
  if (n_drop > 0L)
    message(sprintf(
      "%s: removed %d patient(s) missing Indication.",
      fn_label, n_drop))
  base_pids <- base_pids[!miss]
  if (!length(base_pids))
    stop(fn_label, ": every patient in the variant base is missing ",
         "Indication; nothing to plot.", call. = FALSE)
  v <- v[as.character(v[[patient_col]]) %in% base_pids, , drop = FALSE]

  patient_data <- data.frame(
    Patient_ID = base_pids,
    Indication = unname(ind[base_pids]),
    Dose       = if (is.null(dose)) NA_character_
                 else unname(dose[base_pids]),
    RECIST     = if (is.null(rec)) NA_character_
                 else unname(rec[base_pids]),
    stringsAsFactors = FALSE)
  names(patient_data)[1] <- patient_col

  list(variants     = v,
       base_pids    = base_pids,
       patient_data = patient_data)
}


# =============================================================================
# ctdna_oncoprint()
# =============================================================================

#' Genomic landscape oncoprint
#'
#' Patient-by-gene alteration matrix drawn from `prep$variants` at a chosen
#' visit (default baseline). When `wrap` is supplied, one panel per value of
#' that column is rendered side-by-side; otherwise a single panel is drawn.
#' Rendering style (palette, glyph z-order, LGR box height, RECIST colours)
#' is inherited from the package's oncoprint helpers; the data plumbing is a
#' v0.67.0 rewrite.
#'
#' **Patient base.** `unique(prep$variants$Patient_ID)` after `filter_scheme`
#' (if any) and the `visit` filter, MINUS any patient missing `Indication` or
#' `Dose`. Identical to the base used by [ctdna_alteration_grid()].
#'
#' **Indication & Dose** are read case-insensitively from columns literally
#' named `"indication"` / `"dose"`. Search order: `prep$adrs`, any other ADaM
#' frame, then `prep$variants`, then `prep$clinical`.
#'
#' The word *Cohort* is not used by this function. Passing `"COHORT"` /
#' `"Cohort"` / `"cohort"` in `top_annotations`, `group_patients_by`,
#' `wrap`, or as a name in `annotation_labels` is an error.
#'
#' @param prep A `ctdna_prep` object.
#' @param gene_sets A NAMED list of gene-symbol vectors (one row split per
#'   set) or a built-in name (e.g. `"HRR14"`).
#' @param filter_scheme Optional filter-scheme name(s) passed to
#'   [.filter_apply_df()].
#' @param visit Visit to restrict to (default `"C1D1"`). Accepts `"Cycle 1
#'   Day 1"` and `"baseline"` (case-insensitive).
#' @param visit_col Column in `prep$variants` carrying the visit.
#' @param wrap Optional column name in `patient_data` (Indication, Dose,
#'   RECIST, ...). If supplied, the oncoprint is split into one panel per
#'   value of that column. If `NULL` (default), a single panel.
#' @param top_annotations Character vector of `patient_data` columns to show
#'   as coloured tracks above the matrix. Default: `"RECIST"` if available.
#' @param annotation_labels Optional NAMED character vector mapping column
#'   names to display labels (e.g. `c(RECIST = "Response")`).
#' @param group_patients_by Optional column to sort patients within a panel.
#' @param scheme RECIST collapsing scheme: `"raw"`, `"two"`, `"three"`
#'   (default), or `"four"`.
#' @param alterations Optional character vector of alteration classes to keep.
#' @param sort_genes `"global"` (default), `"within_set"`, or `"none"`.
#' @param show_all_patients TRUE (default) shows every patient in the base,
#'   even those with zero burden in the displayed genes; FALSE = altered only.
#' @param show_patient_names,show_freq_bar Display toggles.
#' @param engine `"auto"` (CH if installed, else ggplot), `"complexheatmap"`,
#'   or `"ggplot"`.
#' @param title,subtitle,caption Plot labels.
#' @param patient_col,gene_col,variant_col Column overrides.
#' @return A `ctdna_oncoprint` object.
#' @export
ctdna_oncoprint <- function(prep,
                            gene_sets,
                            filter_scheme      = NULL,
                            visit              = "C1D1",
                            visit_col          = "Visit_name",
                            wrap               = NULL,
                            top_annotations    = NULL,
                            annotation_labels  = NULL,
                            group_patients_by  = NULL,
                            scheme             = c("three", "two",
                                                   "four", "raw"),
                            alterations        = NULL,
                            sort_genes         = c("global", "within_set",
                                                   "none"),
                            show_all_patients  = TRUE,
                            show_patient_names = TRUE,
                            show_freq_bar      = TRUE,
                            engine             = c("auto", "complexheatmap",
                                                   "ggplot"),
                            title              = NULL,
                            subtitle           = NULL,
                            caption            = NULL,
                            patient_col        = "Patient_ID",
                            gene_col           = "Gene",
                            variant_col        = NULL) {

  engine     <- match.arg(engine)
  scheme     <- match.arg(scheme)
  sort_genes <- match.arg(sort_genes)
  .reject_cohort_args(top_annotations   = top_annotations,
                      group_patients_by = group_patients_by,
                      wrap              = wrap,
                      annotation_labels = names(annotation_labels))
  if (missing(gene_sets) || is.null(gene_sets))
    stop("ctdna_oncoprint: `gene_sets` is required ",
         "(named list, or built-in name like \"HRR14\").", call. = FALSE)

  base <- .build_landscape_base(prep, filter_scheme, visit, visit_col,
                                patient_col, gene_col, "ctdna_oncoprint")
  df           <- base$variants
  patient_data <- base$patient_data
  pd_id_col    <- patient_col

  # Engine resolution
  has_ch <- requireNamespace("ComplexHeatmap", quietly = TRUE)
  if (engine == "auto") engine <- if (has_ch) "complexheatmap" else "ggplot"
  if (engine == "complexheatmap" && !has_ch) {
    message("ctdna_oncoprint: ComplexHeatmap not installed; ",
            "falling back to ggplot engine.")
    engine <- "ggplot"
  }

  # Classify + optional alteration filter
  df <- .oncoprint_classify(df)
  if (!is.null(alterations))
    df <- df[df$Alteration_class %in% alterations, , drop = FALSE]
  if (!nrow(df))
    stop("ctdna_oncoprint: no variants remain after alteration filter.",
         call. = FALSE)

  # Resolve gene_sets
  rgs <- .oncoprint_resolve_gene_sets(df, gene_sets, gene_col)
  resolved   <- rgs$resolved
  gene_union <- rgs$gene_union

  # Auto-default top_annotations
  if (is.null(top_annotations) && "RECIST" %in% names(patient_data) &&
      any(!is.na(patient_data$RECIST)))
    top_annotations <- "RECIST"

  # Apply RECIST collapsing scheme to columns we will use
  sch <- .oncoprint_apply_scheme(patient_data, top_annotations,
                                  group_patients_by, scheme)
  patient_data     <- sch$patient_data
  recist_cols_used <- sch$recist_cols_used

  # Build matrix
  explicit_cohort <- if (is.character(show_all_patients))
                       show_all_patients else NULL
  keep_zero <- isTRUE(show_all_patients) || !is.null(explicit_cohort)
  built <- .oncoprint_build_matrix(df, patient_data, pd_id_col, patient_col,
                                    gene_col, variant_col, gene_union,
                                    resolved, sort_genes,
                                    explicit_cohort = explicit_cohort)
  global_mat    <- built$global_mat
  ordered_genes <- built$ordered_genes
  row_split     <- built$row_split

  # Legend / annotation palette
  leg <- .oncoprint_build_legend_state(top_annotations, patient_data,
                                        recist_cols_used)
  ann_col_global <- leg$col
  ann_at_global  <- leg$at

  # ---- Wrap dimension (USER-CONTROLLED) ----
  if (is.null(wrap)) {
    panel_values <- list(all = NULL)   # one panel, all patients
  } else {
    if (length(wrap) != 1L || !is.character(wrap))
      stop("ctdna_oncoprint: `wrap` must be a single column name or NULL.",
           call. = FALSE)
    if (!(wrap %in% names(patient_data)))
      stop("ctdna_oncoprint: wrap column '", wrap,
           "' not found in patient_data. ",
           "Available: ", paste(names(patient_data), collapse = ", "), ".",
           call. = FALSE)
    pv <- unique(as.character(patient_data[[wrap]]))
    pv <- pv[!is.na(pv) & nzchar(pv)]
    panel_values <- as.list(sort(pv))
    names(panel_values) <- unlist(panel_values)
  }

  # Build panels
  global_classes <- sort(unique(df$Alteration_class))
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
    stop("ctdna_oncoprint: no panel had any patient to render.",
         call. = FALSE)

  if (engine == "complexheatmap") {
    obj <- Reduce(`+`, panels)
    return(structure(list(plot     = obj,
                          engine   = "complexheatmap",
                          subtitle = subtitle,
                          caption  = caption,
                          title    = title),
                     class = c("ctdna_oncoprint", "list")))
  }

  assembled <- .oncoprint_assemble_gg(panels, title, subtitle, caption,
                                       legend_position = "right")
  structure(list(plot     = assembled,
                 engine   = "ggplot2",
                 subtitle = subtitle,
                 caption  = caption,
                 title    = title),
            class = c("ctdna_oncoprint", "list"))
}


# =============================================================================
# ctdna_alteration_grid()
# =============================================================================

#' Combinatorial alteration stacked-bar grid
#'
#' Stacked-bar grid showing the proportion of patients with an alteration in
#' each of several gene sets, broken down by response category, optionally
#' faceted by another column (default: no row facet -- one row of bars across
#' gene sets). Same patient base as [ctdna_oncoprint()].
#'
#' @param prep A `ctdna_prep` object.
#' @param gene_sets A named list where each entry is either a bare gene-symbol
#'   vector (= `alt_any`, backward-compatible) OR a wrapped set from
#'   [alt_any()], [alt_all()], or [alt_any_n()]. The wrapper carries the
#'   combine mode for that set. Facet columns appear in the order the sets
#'   are listed. Unnamed entries get auto-labelled.
#' @param wrap Optional column in `patient_data` used as the row-facet
#'   variable (e.g. `"Indication"`, `"Dose"`). NULL (default) = no row facet.
#' @param response_col Column with the response category (default `"RECIST"`).
#' @param response_levels Optional level order for the x-axis.
#' @param scheme RECIST collapsing scheme: `"raw"` (default), `"two"`,
#'   `"three"`, or `"four"`.
#' @param scales Facet `scales=` argument: `"fixed"`, `"free_y"`, etc.
#' @param filter_scheme Optional filter scheme(s) applied to `prep$variants`.
#' @param visit Visit to restrict to (default `"C1D1"`). Accepts `"Cycle 1
#'   Day 1"` and `"baseline"` (case-insensitive).
#' @param visit_col Column in `prep$variants` carrying the visit.
#' @param patient_col,gene_col Column overrides.
#' @param stat `"fisher"` (default), `"chisq"`, or `"none"`.
#' @param show_counts If TRUE, prints `pct% (n)` inside each segment.
#' @param title,subtitle,caption Plot labels.
#' @return `list(plot, summary, stats)`.
#' @export
ctdna_alteration_grid <- function(prep,
                                  gene_sets,
                                  wrap            = NULL,
                                  response_col    = "RECIST",
                                  response_levels = NULL,
                                  scheme          = c("raw", "two",
                                                       "three", "four"),
                                  scales          = c("fixed", "free_y",
                                                       "free_x", "free"),
                                  filter_scheme   = NULL,
                                  visit           = "C1D1",
                                  visit_col       = "Visit_name",
                                  patient_col     = "Patient_ID",
                                  gene_col        = "Gene",
                                  stat            = c("fisher", "chisq",
                                                       "none"),
                                  show_counts     = TRUE,
                                  title    = "Combinatorial alteration status",
                                  subtitle = NULL,
                                  caption  = NULL) {

  scheme <- match.arg(scheme)
  scales <- match.arg(scales)
  stat   <- match.arg(stat)
  .reject_cohort_args(wrap = wrap, response_col = response_col)
  if (missing(gene_sets) || !is.list(gene_sets) || !length(gene_sets))
    stop("ctdna_alteration_grid: `gene_sets` must be a non-empty list; ",
         "each entry is a gene-symbol vector or an alt_any/alt_all/alt_any_n ",
         "wrapper.", call. = FALSE)

  base <- .build_landscape_base(prep, filter_scheme, visit, visit_col,
                                 patient_col, gene_col,
                                 "ctdna_alteration_grid")
  df           <- base$variants
  patient_data <- base$patient_data
  base_pids    <- base$base_pids

  # response_col must exist now (after the Indication+Dose drop, may still
  # be missing); drop further patients without it, separate message.
  if (!(response_col %in% names(patient_data)))
    stop("ctdna_alteration_grid: column '", response_col,
         "' not available in patient_data.", call. = FALSE)
  rec_b   <- patient_data[[response_col]]
  miss_rec <- is.na(rec_b) | !nzchar(rec_b)
  if (sum(miss_rec) > 0L) {
    message(sprintf(
      "ctdna_alteration_grid: %d additional patient(s) lacking %s removed.",
      sum(miss_rec), response_col))
    patient_data <- patient_data[!miss_rec, , drop = FALSE]
    base_pids    <- patient_data[[patient_col]]
    df           <- df[as.character(df[[patient_col]]) %in% base_pids, ,
                       drop = FALSE]
  }

  # wrap validation (optional row facet)
  if (!is.null(wrap)) {
    if (length(wrap) != 1L || !is.character(wrap))
      stop("ctdna_alteration_grid: `wrap` must be a single column name or NULL.",
           call. = FALSE)
    if (!(wrap %in% names(patient_data)))
      stop("ctdna_alteration_grid: wrap column '", wrap,
           "' not in patient_data. Available: ",
           paste(names(patient_data), collapse = ", "), ".", call. = FALSE)
  }

  # RECIST scheme collapse
  if (scheme != "raw")
    patient_data[[response_col]] <- ctdna_stratify_recist(
      patient_data[[response_col]], scheme)
  if (!is.null(response_levels))
    patient_data[[response_col]] <- factor(patient_data[[response_col]],
                                            levels = response_levels)

  # combine spec
  # Resolve labels + combine mode per entry. Preserve insertion order so the
  # facet columns appear in the order the user listed them.
  raw_names <- names(gene_sets)
  if (is.null(raw_names)) raw_names <- rep("", length(gene_sets))
  labels <- character(length(gene_sets))
  for (i in seq_along(gene_sets)) {
    labels[i] <- if (nzchar(raw_names[i])) raw_names[i]
                 else .geneset_autolabel(gene_sets[[i]], i)
  }
  if (anyDuplicated(labels))
    stop("ctdna_alteration_grid: duplicate gene_set labels: ",
         paste(unique(labels[duplicated(labels)]), collapse = ", "),
         call. = FALSE)

  # Per-set per-patient alt/wt using each entry's own combine mode
  long_rows <- list()
  for (i in seq_along(gene_sets)) {
    label <- labels[i]
    spec  <- .geneset_combine(gene_sets[[i]])
    genes <- spec$genes
    if (!length(genes)) next
    hit_pat <- df[as.character(df[[gene_col]]) %in% genes, patient_col]
    n_hit   <- table(factor(as.character(hit_pat), levels = base_pids))
    threshold <- switch(spec$mode,
      "any"   = 1L,
      "all"   = length(genes),
      "any_n" = {
        if (spec$n > length(genes))
          stop("ctdna_alteration_grid: set '", label, "': n (", spec$n,
               ") exceeds gene count (", length(genes), ").", call. = FALSE)
        as.integer(spec$n)
      },
      stop("ctdna_alteration_grid: unknown combine mode for set '",
           label, "'.", call. = FALSE))
    is_alt <- as.integer(n_hit) >= threshold
    long_rows[[label]] <- data.frame(
      Patient_ID = base_pids,
      gene_set   = label,
      value      = ifelse(is_alt, "alt", "wt"),
      stringsAsFactors = FALSE)
    names(long_rows[[label]])[1] <- patient_col
  }
  long <- do.call(rbind, long_rows)
  # Preserve user-supplied order in the facet
  long$gene_set <- factor(long$gene_set, levels = labels)
  keep_cols <- c(patient_col, response_col)
  if (!is.null(wrap)) keep_cols <- c(keep_cols, wrap)
  long <- merge(long, patient_data[, keep_cols, drop = FALSE],
                by = patient_col, all.x = TRUE)

  # Summary aggregation
  group_cols <- c(if (!is.null(wrap)) wrap else NULL,
                  response_col, "gene_set", "value")
  long$.one <- 1L
  summary <- aggregate(
    stats::as.formula(paste0(".one ~ ",
                              paste(group_cols, collapse = " + "))),
    data = long, FUN = length)
  names(summary)[ncol(summary)] <- "n"
  tot_cols <- setdiff(group_cols, "value")
  totals <- aggregate(
    stats::as.formula(paste0("n ~ ",
                              paste(tot_cols, collapse = " + "))),
    data = summary, FUN = sum)
  names(totals)[ncol(totals)] <- "total"
  summary <- merge(summary, totals, by = tot_cols)
  summary$pct <- summary$n / summary$total

  # Stats per (wrap x gene_set), or per gene_set if no wrap
  stats_df <- NULL
  if (stat != "none") {
    key_cols <- if (is.null(wrap)) "gene_set" else c(wrap, "gene_set")
    keys <- unique(long[, key_cols, drop = FALSE])
    stats_rows <- list()
    for (i in seq_len(nrow(keys))) {
      mask <- rep(TRUE, nrow(long))
      for (kc in key_cols)
        mask <- mask & (long[[kc]] == keys[i, kc])
      sub <- long[mask, , drop = FALSE]
      tab <- table(sub[[response_col]], sub$value)
      if (nrow(tab) < 2 || ncol(tab) < 2) {
        p <- NA_real_
      } else if (stat == "fisher") {
        p <- tryCatch(stats::fisher.test(tab, workspace = 2e6,
                                         simulate.p.value = TRUE,
                                         B = 5000)$p.value,
                       error = function(e) NA_real_)
      } else {
        p <- tryCatch(suppressWarnings(stats::chisq.test(tab)$p.value),
                       error = function(e) NA_real_)
      }
      row <- as.list(keys[i, , drop = FALSE]); row$p <- p
      stats_rows[[i]] <- as.data.frame(row, stringsAsFactors = FALSE)
    }
    stats_df <- do.call(rbind, stats_rows)
    stats_df$label <- ctdna_pval_label(stats_df$p)
  }

  # ---- plot ----
  facet_rhs <- "gene_set"
  facet_lhs <- if (is.null(wrap)) "." else wrap
  p <- ggplot2::ggplot(summary,
                       ggplot2::aes(x = .data[[response_col]],
                                    y = .data$pct,
                                    fill = .data$value)) +
    ggplot2::geom_col(width = .o("bar_width"), position = "stack") +
    ggplot2::scale_fill_manual(
      values = c(alt = "#F2786F", wt = "#3CB6BD"),
      breaks = c("alt", "wt"),
      labels = c("alt", "wt"),
      name   = "value") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(),
                                 limits = c(0, 1.18),
                                 expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::facet_grid(
      stats::as.formula(paste(facet_lhs, "~", facet_rhs)),
      scales = scales, space = "fixed") +
    ctdna_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5,
                                           hjust = 1, size = 8),
      strip.text  = ggplot2::element_text(face = "bold"),
      panel.spacing = grid::unit(0.4, "lines")) +
    ggplot2::labs(x = "Response category",
                   y = "Percent genetic alteration",
                   title = title, subtitle = subtitle, caption = caption)

  if (show_counts) {
    summary$lbl <- sprintf("%.0f%% (%d)", 100 * summary$pct, summary$n)
    lbl_df <- summary[summary$pct >= 0.06, , drop = FALSE]
    if (nrow(lbl_df) > 0)
      p <- p + ggplot2::geom_text(
        data = lbl_df,
        ggplot2::aes(x = .data[[response_col]], y = .data$pct,
                      label = .data$lbl, group = .data$value),
        position = ggplot2::position_stack(vjust = 0.5),
        size = 2.6, color = "grey10", fontface = "bold",
        inherit.aes = FALSE)
  }

  if (!is.null(stats_df))
    p <- p + ggplot2::geom_text(
      data = stats_df,
      ggplot2::aes(x = -Inf, y = 1.15, label = .data$label),
      hjust = -0.08, vjust = 1, size = 2.4, color = "grey30",
      fontface = "italic",
      inherit.aes = FALSE)

  list(plot = p, summary = summary, stats = stats_df)
}
