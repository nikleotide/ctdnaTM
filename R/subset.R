# ============================================================================
# ctdnaTM -- ctdna_subset()
# ----------------------------------------------------------------------------
# v0.80.0: one place to carve a `ctdna_prep` down to a cohort of interest.
#
# Two kinds of restriction, deliberately distinguished:
#
#   * PATIENT-LEVEL (indication / dose / recist / mr / patient) -- resolves a
#     set of Patient_IDs and cascades to EVERY slot that carries a Patient_ID
#     column, so variants, assessments, IHC, DNA-seq and RNA-seq all stay
#     consistent with $samples and $clinical.
#
#   * ROW-LEVEL (visit / gene) -- filters rows inside the frames that carry the
#     relevant column. Patients are NOT dropped; a patient with no rows left
#     simply contributes nothing to that frame.
#
# All value matching goes through .ctdna_norm_ind(), so case and separators are
# irrelevant: dose = "8mg/kg" matches "8 mg/kg", indication = "nsq-nsclc"
# matches "Nsq_NSCLC".
# ============================================================================


# Resolve a column in a data frame by normalised name, trying each candidate in
# order. Returns the actual column name, or NA_character_.
.ctdna_subset_col <- function(df, candidates) {
  if (!is.data.frame(df)) return(NA_character_)
  nms <- names(df)
  key <- .ctdna_norm_ind(nms)
  for (cd in candidates) {
    if (is.null(cd) || !nzchar(cd)) next
    hit <- which(key == .ctdna_norm_ind(cd))
    if (length(hit)) return(nms[hit[1]])
  }
  NA_character_
}

# Patient_IDs whose value in `col` normalises into `values`. Searches
# $clinical first, then $samples, so a column present on only one still works.
.ctdna_subset_pids <- function(prep, candidates, values, what) {
  frames <- list(clinical = prep$clinical, samples = prep$samples)
  for (nm in names(frames)) {
    d <- frames[[nm]]
    if (!is.data.frame(d) || !"Patient_ID" %in% names(d)) next
    col <- .ctdna_subset_col(d, candidates)
    if (is.na(col)) next
    have <- as.character(d[[col]])
    sel  <- .ctdna_norm_ind(have) %in% .ctdna_norm_ind(values)
    if (!any(sel)) {
      avail <- sort(unique(have[!is.na(have)]))
      stop(sprintf(
        paste0("ctdna_subset: none of {%s} matched any value in `%s` (from $%s).\n",
               "  Available: %s\n",
               "  Matching ignores case and separators."),
        paste(values, collapse = ", "), col, nm,
        if (length(avail)) paste(avail, collapse = ", ") else "<none>"),
        call. = FALSE)
    }
    unmatched <- values[!(.ctdna_norm_ind(values) %in% .ctdna_norm_ind(have))]
    if (length(unmatched))
      warning(sprintf("ctdna_subset: %s matched nothing in `%s` and was ignored.",
                       paste(sprintf("`%s`", unmatched), collapse = ", "), col),
              call. = FALSE)
    return(unique(as.character(d$Patient_ID)[sel]))
  }
  stop(sprintf(
    paste0("ctdna_subset: no column for `%s` found in $clinical or $samples.\n",
           "  Looked for (case- and separator-insensitive): %s"),
    what, paste(candidates, collapse = ", ")), call. = FALSE)
}


#' Subset a ctdna_prep to a cohort of interest
#'
#' Carve a \code{ctdna_prep} object down to the patients, visits or genes you
#' want to analyse, keeping every slot mutually consistent. The result is a
#' \code{ctdna_prep} with the same structure and the same columns as the input
#' -- only the number of records changes -- so it drops straight into every
#' plotting and analysis function in the package.
#'
#' @details
#' Two kinds of restriction are available, and they behave differently on
#' purpose.
#'
#' \strong{Patient-level} (\code{indication}, \code{dose}, \code{recist},
#' \code{mr}, \code{patient}) resolves a set of \code{Patient_ID}s and cascades
#' the restriction to \emph{every} slot that carries a \code{Patient_ID}
#' column. Variants, assessments, IHC, tissue DNA-seq and RNA-seq all stay in
#' step with \code{$samples} and \code{$clinical}, so downstream statistics and
#' per-group \code{n} labels reflect the subset rather than the full study.
#'
#' \strong{Row-level} (\code{visit}, \code{gene}) filters rows inside the
#' frames that carry the relevant column. Patients are \emph{not} dropped: a
#' patient with no matching rows simply contributes nothing to that frame.
#'
#' Several arguments can be combined; they are AND-ed together. Every value is
#' matched case- and separator-insensitively (the same rule
#' \code{indication = } uses across the package), so \code{dose = "8mg/kg"}
#' matches a stored \code{"8 mg/kg"} and \code{indication = "nsq-nsclc"}
#' matches \code{"Nsq_NSCLC"}.
#'
#' A value that matches nothing raises an error listing what \emph{is}
#' available, rather than silently returning an empty object. If several values
#' are given and only some match, the call proceeds and warns about the rest.
#'
#' @param prep A \code{ctdna_prep} object, typically from
#'   \code{\link{ctdna_prepare}}.
#' @param indication Patient-level. Indication / cohort label(s) to keep, e.g.
#'   \code{"Nsq_NSCLC"} or \code{c("Nsq_NSCLC","ES_SCLC")}. Resolved against
#'   the canonical \code{Indication} column via \code{ctdna_opts("indication")}
#'   with fallbacks to \code{Indication}, \code{Cancertype} and
#'   \code{Cancer_Type}.
#' @param dose Patient-level. Dose level(s) to keep, e.g. \code{"8 mg/kg"}.
#'   Resolved via \code{ctdna_opts("dose")}.
#' @param recist Patient-level. Best-response category(ies) to keep, e.g.
#'   \code{"PD"} or \code{c("CR","PR")}. When \code{scheme} is supplied the
#'   raw RECIST values are collapsed first, so you can pass collapsed labels
#'   such as \code{"CR+PR"}.
#' @param scheme Optional RECIST-collapsing scheme applied before \code{recist}
#'   is matched. One of \code{"raw"} (default -- match the stored values),
#'   \code{"two"}, \code{"three"}, \code{"four"}, \code{"four_alt"},
#'   \code{"five"}, \code{"six"}.
#' @param mr Patient-level. Molecular-response call(s) to keep. Resolved via
#'   \code{ctdna_opts("mr")}, falling back to an \code{MR} column.
#' @param patient Patient-level. Explicit \code{Patient_ID} vector to keep.
#'   Matched exactly (no normalisation), since subject identifiers are
#'   case-significant.
#' @param visit Row-level. Visit label(s) to keep in the frames that carry a
#'   visit column, e.g. \code{"baseline"} or \code{c("C2D1","C3D1")}.
#' @param gene Row-level. Gene symbol(s) to keep in \code{$variants} and
#'   \code{$dnaseq}.
#' @param verbose Print a one-line summary of what was kept. Default
#'   \code{TRUE}.
#'
#' @return A \code{ctdna_prep} object. Same class, same slots, same columns as
#'   \code{prep}; only the number of records changes.
#'
#' @seealso \code{\link{ctdna_prepare}} to build a prep;
#'   \code{\link{ctdna_plot_baseline}} and friends, which all accept an
#'   \code{indication} argument for the same restriction inline.
#'
#' @examples
#' \dontrun{
#' sim  <- ctdna_make_mock_study(n_patients = 200, seed = 42)
#' prep <- ctdna_prepare(infinity_report = sim$infinity_report,
#'                       adam = list(adsl = sim$clinical, adtr = sim$adtr))
#'
#' # one cohort, then analyse it repeatedly
#' nsclc <- ctdna_subset(prep, indication = "NSCLC")
#' ctdna_plot_baseline(nsclc, group_by = "Dose", subgroup_by = "RECIST",
#'                     compare_by = "subgroup_within_group")
#' ctdna_plot_reduction(nsclc, scheme = "three")
#'
#' # combine restrictions (AND-ed)
#' ctdna_subset(prep, indication = c("NSCLC","BRCA"), dose = "High")
#'
#' # case and separators are irrelevant
#' ctdna_subset(prep, indication = "nsclc")
#'
#' # progressors only, collapsing RECIST first
#' ctdna_subset(prep, recist = "PD+NE+NA", scheme = "three")
#'
#' # row-level: baseline rows only, patients kept
#' ctdna_subset(prep, visit = "baseline")
#'
#' # row-level: restrict the variant tables to a gene panel
#' ctdna_subset(prep, gene = c("TP53", "PTEN", "RB1"))
#' }
#'
#' @export
ctdna_subset <- function(prep,
                          indication = NULL,
                          dose       = NULL,
                          recist     = NULL,
                          scheme     = c("raw","two","three","four",
                                          "four_alt","five","six"),
                          mr         = NULL,
                          patient    = NULL,
                          visit      = NULL,
                          gene       = NULL,
                          verbose    = TRUE) {
  if (!is.list(prep))
    stop("ctdna_subset: `prep` must be a ctdna_prep object (see ctdna_prepare()).",
         call. = FALSE)
  scheme <- match.arg(scheme)

  drop_empty <- function(x) x[!is.na(x) & nzchar(as.character(x))]
  indication <- drop_empty(indication); dose    <- drop_empty(dose)
  recist     <- drop_empty(recist);     mr      <- drop_empty(mr)
  patient    <- drop_empty(patient);    visit   <- drop_empty(visit)
  gene       <- drop_empty(gene)

  if (!length(indication) && !length(dose) && !length(recist) &&
      !length(mr) && !length(patient) && !length(visit) && !length(gene)) {
    if (verbose)
      message("ctdna_subset: no restriction supplied; returning prep unchanged.")
    return(prep)
  }

  n_before <- sum(vapply(prep, function(d)
    if (is.data.frame(d) && "Patient_ID" %in% names(d)) nrow(d) else 0L,
    integer(1)))
  pid_before <- if (is.data.frame(prep$clinical))
    length(unique(as.character(prep$clinical$Patient_ID))) else NA_integer_
  notes <- character(0)

  # ---- patient-level restrictions (intersected) ----------------------------
  keep_pids <- NULL
  add <- function(pids) {
    if (is.null(keep_pids)) pids else intersect(keep_pids, pids)
  }

  if (length(indication)) {
    keep_pids <- add(.ctdna_subset_pids(
      prep, c(.o("indication"), "Indication", "Cancertype", "Cancer_Type"),
      indication, "indication"))
    notes <- c(notes, sprintf("indication={%s}", paste(indication, collapse = ",")))
  }
  if (length(dose)) {
    keep_pids <- add(.ctdna_subset_pids(
      prep, c(.o("dose"), "Dose"), dose, "dose"))
    notes <- c(notes, sprintf("dose={%s}", paste(dose, collapse = ",")))
  }
  if (length(recist)) {
    rc <- .ctdna_subset_col(prep$clinical, c(.o("recist"), "RECIST"))
    if (is.na(rc)) rc <- .ctdna_subset_col(prep$samples, c(.o("recist"), "RECIST"))
    if (is.na(rc))
      stop("ctdna_subset: no RECIST column found in $clinical or $samples.",
           call. = FALSE)
    src <- if (is.data.frame(prep$clinical) && rc %in% names(prep$clinical))
      prep$clinical else prep$samples
    vals <- as.character(src[[rc]])
    if (!identical(scheme, "raw"))
      vals <- as.character(ctdna_stratify_recist(vals, scheme))
    sel <- .ctdna_norm_ind(vals) %in% .ctdna_norm_ind(recist)
    if (!any(sel)) {
      avail <- sort(unique(vals[!is.na(vals)]))
      stop(sprintf(
        paste0("ctdna_subset: none of {%s} matched any value in `%s`%s.\n",
               "  Available: %s"),
        paste(recist, collapse = ", "), rc,
        if (identical(scheme, "raw")) "" else sprintf(" (scheme = \"%s\")", scheme),
        paste(avail, collapse = ", ")), call. = FALSE)
    }
    keep_pids <- add(unique(as.character(src$Patient_ID)[sel]))
    notes <- c(notes, sprintf("recist={%s}%s", paste(recist, collapse = ","),
                               if (identical(scheme, "raw")) ""
                               else sprintf(" scheme=%s", scheme)))
  }
  if (length(mr)) {
    keep_pids <- add(.ctdna_subset_pids(
      prep, c(.o("mr"), "MR", "MR_call", "Molecular_response"), mr, "mr"))
    notes <- c(notes, sprintf("mr={%s}", paste(mr, collapse = ",")))
  }
  if (length(patient)) {
    # exact match -- subject identifiers are case-significant
    known <- unique(c(
      if (is.data.frame(prep$clinical)) as.character(prep$clinical$Patient_ID),
      if (is.data.frame(prep$samples))  as.character(prep$samples$Patient_ID)))
    miss <- setdiff(as.character(patient), known)
    if (length(miss) == length(patient))
      stop(sprintf("ctdna_subset: none of the %d requested Patient_ID(s) exist in this prep.",
                    length(patient)), call. = FALSE)
    if (length(miss))
      warning(sprintf("ctdna_subset: %d Patient_ID(s) not found and ignored: %s",
                       length(miss), paste(utils::head(miss, 5), collapse = ", ")),
              call. = FALSE)
    keep_pids <- add(intersect(as.character(patient), known))
    notes <- c(notes, sprintf("patient=%d id(s)", length(patient)))
  }

  if (!is.null(keep_pids)) {
    if (!length(keep_pids))
      stop("ctdna_subset: the combined restriction matched 0 patients. ",
           "Relax one of the arguments (they are AND-ed together).",
           call. = FALSE)
    for (slot in names(prep)) {
      d <- prep[[slot]]
      if (is.data.frame(d) && "Patient_ID" %in% names(d))
        prep[[slot]] <- d[as.character(d$Patient_ID) %in% keep_pids, , drop = FALSE]
    }
  }

  # ---- row-level restrictions ---------------------------------------------
  if (length(visit)) {
    # Only genuinely per-visit frames. $clinical is one row per patient and can
    # carry a constant enrolment label (e.g. Visit_name = "Screening"), so
    # filtering it to an on-treatment visit would drop every subject.
    per_visit <- setdiff(names(prep), c("clinical", "dictionary", "qc_removed",
                                         "filter_explanation"))
    want_baseline <- any(.ctdna_norm_ind(visit) == "baseline")
    matched <- character(0); skipped <- character(0)
    for (slot in per_visit) {
      d <- prep[[slot]]
      if (!is.data.frame(d)) next
      vc <- .ctdna_subset_col(d, c(.o("visit"), .o("time"), "Visit_name",
                                    "Visit", "Tumor_visit_name"))
      if (is.na(vc)) next
      vv <- as.character(d[[vc]])
      lv <- unique(vv[!is.na(vv)])
      # A single constant value is an enrolment/label column, not a visit axis.
      if (length(lv) <= 1L) next
      want <- .ctdna_norm_ind(visit)
      # Visit vocabularies differ between frames: $samples carries normalised
      # labels ("Baseline") while $variants keeps raw protocol labels
      # ("C1D1"). Resolve "baseline" to this frame's first visit level.
      if (want_baseline) {
        ord <- if (is.factor(d[[vc]])) levels(droplevels(d[[vc]])) else sort(lv)
        want <- unique(c(want, .ctdna_norm_ind(ord[1])))
      }
      keep <- .ctdna_norm_ind(vv) %in% want
      if (!any(keep)) { skipped <- c(skipped, slot); next }
      prep[[slot]] <- d[keep, , drop = FALSE]
      matched <- c(matched, slot)
    }
    if (!length(matched)) {
      av <- unique(unlist(lapply(per_visit, function(s) {
        d <- prep[[s]]; if (!is.data.frame(d)) return(NULL)
        vc <- .ctdna_subset_col(d, c(.o("visit"), .o("time"), "Visit_name",
                                      "Visit", "Tumor_visit_name"))
        if (is.na(vc)) NULL else unique(as.character(d[[vc]]))
      })))
      stop(sprintf(
        paste0("ctdna_subset: none of {%s} matched any visit in this prep.\n",
               "  Available: %s\n",
               "  Note: \"baseline\" resolves to each frame's first visit level."),
        paste(visit, collapse = ", "),
        paste(sort(av[!is.na(av)]), collapse = ", ")), call. = FALSE)
    }
    if (length(skipped))
      message("ctdna_subset: `visit` matched nothing in $",
              paste(skipped, collapse = ", $"),
              " (different visit vocabulary); left unfiltered.")
    notes <- c(notes, sprintf("visit={%s}", paste(visit, collapse = ",")))
  }

  if (length(gene)) {
    hit <- FALSE
    for (slot in c("variants", "dnaseq")) {
      d <- prep[[slot]]
      if (!is.data.frame(d)) next
      gc <- .ctdna_subset_col(d, c(.o("gene"), "Gene"))
      if (is.na(gc)) next
      keep <- .ctdna_norm_ind(as.character(d[[gc]])) %in% .ctdna_norm_ind(gene)
      prep[[slot]] <- d[keep, , drop = FALSE]
      hit <- TRUE
    }
    if (!hit)
      warning("ctdna_subset: no Gene column found in $variants or $dnaseq; ",
              "`gene` ignored.", call. = FALSE)
    else notes <- c(notes, sprintf("gene={%s}", paste(gene, collapse = ",")))
  }

  if (isTRUE(verbose)) {
    n_after <- sum(vapply(prep, function(d)
      if (is.data.frame(d) && "Patient_ID" %in% names(d)) nrow(d) else 0L,
      integer(1)))
    pid_after <- if (is.data.frame(prep$clinical))
      length(unique(as.character(prep$clinical$Patient_ID))) else NA_integer_
    message(sprintf(
      "ctdna_subset: %s -> %s patient(s), %s of %s row(s) across prep slots.",
      if (is.na(pid_before)) "?" else format(pid_before),
      if (is.na(pid_after))  "?" else format(pid_after),
      format(n_after), format(n_before)))
    if (length(notes)) message("  kept: ", paste(notes, collapse = "; "))
  }

  prep
}
