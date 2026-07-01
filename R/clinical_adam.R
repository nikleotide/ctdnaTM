# =============================================================================
# v0.35.0 — Clinical patient-data assembly (single function)
# =============================================================================
#
# ctdna_make_patient_data(...,
#                          Cancer_Type   = <REQUIRED 2-col data frame>,
#                          extra_cols    = NULL,
#                          verbose       = TRUE)
#
# Pulls clinical info from any number of data frames passed via `...`
# (typically adsl/adrs/adtr but custom frames are fine too) and assembles
# a per-patient frame with a fixed canonical schema:
#
#   Patient_ID, Dose, Cancer_Type, ARM, Response_Subcategory, Response
#   + any user-requested extra_cols
#
# Cancer_Type values come from the REQUIRED `Cancer_Type` argument
# (a 2-column data frame: ID + indication). This guarantees the user
# has explicitly curated their indication labels before plotting.
#
# Rows with empty Cancer_Type OR empty Dose are dropped (with messages).
#
# Response_Subcategory: 4-level factor (CR/PR / SD / PD / NE/NA) derived
# from BCR rows in adrs.
# Response: binary R / NR derived from Response_Subcategory
# (CR/PR -> R; SD/PD/NE/NA -> NR).
#
# Step 2 (harmonization) is GONE in v0.35.0. If the user has multiple
# labels for the same indication, they harmonize them BEFORE calling
# this function by editing the `Cancer_Type` data frame they pass in.


.clinical_env <- new.env(parent = emptyenv())

.clinical_default_opts <- function() list(
  # v0.38.0: ADaM-first ID priority (SUBJID before Patient_ID).
  # Order: ADaM standard (SUBJID, USUBJID) -> non-ADaM (Patient_ID).
  id_col_priority = c("SUBJID","USUBJID","Patient_ID","subject_id",
                       "PatientID","patient_id"),
  bcr_param       = "Best Confirmed Response by Investigator (RECIST 1.1)",
  bcr_value_col   = "AVALC",
  tumor_pchg_col       = "PCHG",
  tumor_parcat_col     = "PARCAT1",
  tumor_parcat_target  = c("Target Lesion","Target Lesions"))

.clinical_init <- function() {
  if (is.null(.clinical_env$cfg))
    .clinical_env$cfg <- .clinical_default_opts()
}

#' Get / set clinical data options
#'
#' Companion to \code{\link{ctdna_opts}} for the clinical-data
#' subsystem (\code{\link{ctdna_make_patient_data}}).
#'
#' @param ... Named option = value pairs to set, OR a single unnamed
#'   string to look up one option, OR nothing to return the full list.
#' @param .reset If \code{TRUE}, restore factory defaults.
#' @return Invisibly the (possibly updated) full options list.
#' @examples
#' # Inspect current clinical-ADaM defaults
#' clinical_opts()
#'
#' # Override one option for this session
#' clinical_opts(USUBJID_col = "USUBJID")
#' @export
clinical_opts <- function(..., .reset = FALSE) {
  .clinical_init()
  if (isTRUE(.reset)) {
    .clinical_env$cfg <- .clinical_default_opts()
    return(invisible(.clinical_env$cfg))
  }
  args <- list(...)
  if (length(args) == 0)
    return(invisible(.clinical_env$cfg))
  if (length(args) == 1L && (is.null(names(args)) || !nzchar(names(args)[1]))) {
    key <- args[[1]]
    if (!is.character(key) || length(key) != 1L)
      stop("To look up a single option, pass its name as a character string.",
            call. = FALSE)
    return(.clinical_env$cfg[[key]])
  }
  nm <- names(args)
  if (is.null(nm) || any(!nzchar(nm)))
    stop("All arguments to clinical_opts() must be named.", call. = FALSE)
  for (k in nm) {
    if (!k %in% names(.clinical_default_opts()))
      stop(sprintf("Unknown clinical_opts() key: '%s'. See ?clinical_opts.", k),
            call. = FALSE)
    .clinical_env$cfg[[k]] <- args[[k]]
  }
  invisible(.clinical_env$cfg)
}

.co <- function(key) { .clinical_init(); .clinical_env$cfg[[key]] }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Normalize a column name for fuzzy matching:
# strip non-alphanumerics, lowercase.
.norm_colname <- function(x) gsub("[^a-z0-9]", "", tolower(x))

# Patient ID column candidates (canonical forms after .norm_colname)
.id_canonical <- c("patientid","subjid","usubjid","subjectid",
                    "subject","subj","pid")

# Indication column candidates (canonical forms after .norm_colname)
.indication_canonical <- c("cancertype","indication","clinindication",
                            "diagnosis","tumortype")

# Find a column in a frame whose normalized name matches any in `targets`.
# Returns the actual column name, or NULL.
.find_col_by_norm <- function(df, targets) {
  nm  <- names(df)
  nrm <- .norm_colname(nm)
  hit <- which(nrm %in% targets)
  if (length(hit) == 0L) return(NULL)
  nm[hit[1]]
}

# RECIST -> 4-level factor (CR/PR / SD / PD / NE/NA)
.recist_to_subcategory <- function(x) {
  v <- as.character(x)
  out <- rep("NE/NA", length(v))
  for (i in seq_along(v)) {
    val <- v[i]
    if (is.na(val) || !nzchar(val) ||
        val %in% c("NE","NA","Not Evaluable","Not Available"))
      out[i] <- "NE/NA"
    else if (grepl("CR|PR|Complete|Partial", val, ignore.case = TRUE))
      out[i] <- "CR/PR"
    else if (grepl("^PD|Progressive", val, ignore.case = TRUE))
      out[i] <- "PD"
    else if (grepl("^SD|Stable", val, ignore.case = TRUE))
      out[i] <- "SD"
  }
  factor(out, levels = c("CR/PR","SD","PD","NE/NA"))
}

# Response_Subcategory -> binary R / NR
.subcategory_to_response <- function(f) {
  v <- as.character(f)
  out <- rep("NR", length(v))
  out[v == "CR/PR"] <- "R"
  factor(out, levels = c("R","NR"))
}

# Locate the patient ID column on a generic input frame
.find_id_col <- function(df, priority = NULL) {
  if (is.null(priority)) priority <- .co("id_col_priority")
  hit <- intersect(priority, names(df))
  if (length(hit) > 0L) return(hit[1])
  # Fall back to fuzzy match
  .find_col_by_norm(df, .id_canonical)
}

# Classify an input data frame by column signature
.classify_clinical_frame <- function(df) {
  nm <- names(df)
  has_param  <- "PARAM"  %in% nm
  has_avalc  <- "AVALC"  %in% nm
  has_parcat <- any(c("PARCAT1","PARCAT2") %in% nm)
  has_pchg   <- "PCHG"   %in% nm
  has_arm    <- any(c("ARM","ARMCD","ACTARM","ACTARMCD") %in% nm)
  if (has_param && has_parcat && has_pchg) return("adtr")
  if (has_param && has_avalc) return("adrs")
  if (any(c("STUDYID","SAFFL","ITTFL","TRT01P","TRT01A") %in% nm) && has_arm &&
      !has_param)
    return("adsl")
  "custom"
}

# Clean column values: unlist haven list-columns, trim whitespace, factor->char
.normalize_columns <- function(df) {
  for (col in names(df)) {
    v <- df[[col]]
    if (is.list(v) && !is.data.frame(v)) {
      lens <- vapply(v, length, integer(1))
      if (all(lens <= 1L)) {
        df[[col]] <- vapply(v, function(x)
          if (length(x) == 0L) NA else as.character(x[[1]]),
          character(1))
      } else {
        warning(sprintf(
          "Column '%s' is a list-column with multi-element rows; first element kept.",
          col), call. = FALSE)
        df[[col]] <- vapply(v, function(x)
          if (length(x) == 0L) NA_character_ else as.character(x[[1]]),
          character(1))
      }
    }
    if (is.factor(df[[col]])) df[[col]] <- as.character(df[[col]])
    if (is.character(df[[col]])) df[[col]] <- trimws(df[[col]])
  }
  df
}

# Extract dose from ARM/TRT01A free-text. Returns NA if no match.
.extract_dose_from_arm <- function(x) {
  v <- as.character(x)
  out <- rep(NA_character_, length(v))
  pos <- regexpr("[0-9]+(?:\\.[0-9]+)?\\s*mg(?:/kg)?",
                  v, perl = TRUE, ignore.case = TRUE)
  hit <- !is.na(pos) & pos != -1
  if (any(hit)) {
    m <- regmatches(v[hit], regexpr("[0-9]+(?:\\.[0-9]+)?\\s*mg(?:/kg)?",
                                       v[hit], perl = TRUE, ignore.case = TRUE))
    out[hit] <- m
  }
  out
}

# Source-column priority lists for the non-Cancer_Type canonical cols
.arm_priority    <- c("ARM", "ACTARM", "TRT01A", "treatment_arm")
.dose_priority   <- c("Dose","DOSE","dose")

# Extract a canonical value from frames_meta by walking source priority.
# First non-NA per patient wins. ADaM frames sorted first via caller.
.extract_canonical <- function(canonical_col, source_cols, frames_meta,
                                pat_ids, special_handlers = NULL) {
  out <- rep(NA_character_, length(pat_ids))
  for (fm in frames_meta) {
    df    <- fm$df
    kind  <- fm$kind
    id_col <- fm$id_col
    if (is.null(id_col)) next

    pulled <- NULL
    if (!is.null(special_handlers[[kind]])) {
      pulled <- special_handlers[[kind]](df, id_col, pat_ids)
    }
    if (is.null(pulled)) {
      sc <- intersect(source_cols, names(df))
      if (length(sc) == 0L) next
      use_col <- sc[1]
      idx <- match(pat_ids, .norm_pid(df[[id_col]]))
      pulled <- as.character(df[[use_col]])[idx]
      pulled <- trimws(pulled)
      pulled[!nzchar(pulled) | pulled == "NA"] <- NA_character_
    }
    to_fill <- is.na(out) & !is.na(pulled)
    out[to_fill] <- pulled[to_fill]
  }
  out
}


# ---------------------------------------------------------------------------
# Public function
# ---------------------------------------------------------------------------

#' Assemble a canonical patient data frame
#'
#' Takes any number of clinical data frames (ADaM datasets like
#' adsl/adrs/adtr, custom per-patient frames, anything with a patient
#' ID column) and returns a single per-patient frame with a fixed
#' canonical schema. Designed so the user never has to pre-merge or
#' pre-clean.
#'
#' \strong{The \code{Cancer_Type} argument is REQUIRED.} It carries
#' the canonical indication label for each patient. This forces the
#' user to curate their indication labels (e.g. collapse \dQuote{SCLC}
#' and \dQuote{Small cell lung cancer} to one canonical string)
#' \emph{before} calling this function. There is no fuzzy matching or
#' auto-collapse in v0.35.0+.
#'
#' The output always has these columns in this order:
#' \enumerate{
#'   \item \code{Patient_ID}
#'   \item \code{Dose} -- character; from a dedicated dose column if
#'     present, else regex-extracted from ARM / TRT01A text
#'     (e.g. \dQuote{8 mg/kg} from \dQuote{GSK227 8 mg/kg Q3W})
#'   \item \code{Cancer_Type} -- from the \code{Cancer_Type} argument
#'   \item \code{ARM} -- treatment arm
#'   \item \code{Response_Subcategory} -- 4-level factor:
#'     \code{CR/PR} / \code{SD} / \code{PD} / \code{NE/NA}, derived
#'     from raw BCR rows in adrs
#'   \item \code{Response} -- binary factor: \code{R} (responder) when
#'     \code{Response_Subcategory == "CR/PR"}, else \code{NR}
#' }
#' Plus any \code{extra_cols} you request, appended at the end.
#'
#' \strong{Rows with empty Cancer_Type OR empty Dose are dropped}
#' from the output, with messages naming how many rows were dropped
#' and showing up to 5 affected Patient_IDs each.
#'
#' @section The \code{Cancer_Type} argument:
#' Must be a data frame with at least 2 columns: patient ID and
#' indication label. Column names are auto-detected by fuzzy matching
#' (case + punctuation insensitive). Recognized ID column names
#' include \code{Patient_ID}, \code{SUBJID}, \code{USUBJID},
#' \code{Subject_ID}, \code{Patient-ID}, etc. Recognized indication
#' column names include \code{Cancer_Type}, \code{Cancertype},
#' \code{Indication}, \code{Clin_Indication}, \code{cancer-type}, etc.
#' Extra columns in this data frame are ignored -- only the ID and
#' indication are used.
#'
#' If your indication labels are messy (e.g. some patients have
#' \dQuote{SCLC} and others have \dQuote{Small cell lung cancer}),
#' harmonize them in the \code{Cancer_Type} data frame BEFORE calling
#' this function. Example:
#' \preformatted{
#' Master_Clinical$Clin_Indication <- ifelse(
#'   Master_Clinical$Clin_Indication \%in\% c("SCLC","Small cell lung cancer"),
#'   "SCLC", Master_Clinical$Clin_Indication)
#' }
#'
#' @param ... Any number of clinical data frames (ADaM or custom).
#'   Order doesn't matter except when frames disagree on a value --
#'   first frame in the call wins. ADaM frames are auto-detected and
#'   given priority on conflicts.
#' @param Cancer_Type \strong{REQUIRED}. A data frame with at least 2
#'   columns: patient ID + indication label. See description.
#' @param extra_cols Optional character vector of additional column
#'   names to pull through from any input frame. Appended at the end
#'   of the output.
#' @param verbose If \code{TRUE} (default), print \code{head()} and a
#'   Cancer_Type frequency table. Set \code{FALSE} for silent scripted
#'   use.
#' @return A \code{ctdna_patient_data} object (a data frame with the
#'   canonical schema described above).
#' @seealso \code{\link{clinical_opts}}, \code{\link{ctdna_oncoprint}}.
#' @examples
#' \dontrun{
#' # Curate indication labels first
#' Master_Clinical$Clin_Indication <- harmonize_my_labels(
#'   Master_Clinical$Clin_Indication)
#'
#' patient_df <- ctdna_make_patient_data(
#'   adsl, adrs, adtr,
#'   Cancer_Type = Master_Clinical[c("Patient_ID","Clin_Indication")],
#'   extra_cols  = c("AGE","SEX"))
#' }
#' @export
ctdna_make_patient_data <- function(..., Cancer_Type = NULL, extra_cols = NULL,
                                      verbose = TRUE) {
  # v0.38.0: Cancer_Type is now OPTIONAL.
  #   - If provided as a data frame: existing behavior (user-curated labels)
  #   - If NULL: auto-derive from ADSL using ctdna_opts("col_indication_adam")
  #     priority list (default PRTUMTY -> COHORT -> STUDYTRT)
  #   - When user data is supplied via `...`, it can override/supplement the
  #     ADaM-derived Cancer_Type
  use_adam_indication <- is.null(Cancer_Type)
  if (!use_adam_indication) {
    if (!is.data.frame(Cancer_Type))
      stop("ctdna_make_patient_data: `Cancer_Type` must be a data frame (or NULL ",
            "to derive from ADSL). Got: ", class(Cancer_Type)[1], call. = FALSE)
    if (ncol(Cancer_Type) < 2L)
      stop("ctdna_make_patient_data: `Cancer_Type` must have at least 2 columns ",
            "(patient ID + indication). Got ", ncol(Cancer_Type), " column(s).",
            call. = FALSE)
    if (nrow(Cancer_Type) == 0L)
      stop("ctdna_make_patient_data: `Cancer_Type` data frame is empty.",
            call. = FALSE)

    Cancer_Type <- .normalize_columns(Cancer_Type)

    ct_id_col   <- .find_col_by_norm(Cancer_Type, .id_canonical)
    ct_ind_col  <- .find_col_by_norm(Cancer_Type, .indication_canonical)
    if (is.null(ct_id_col))
      stop("ctdna_make_patient_data: `Cancer_Type` frame has no recognizable ",
            "patient ID column. Looked for (case/punct insensitive): ",
            "Patient_ID, SUBJID, USUBJID, Subject_ID. Got columns: ",
            paste(names(Cancer_Type), collapse = ", "), ".", call. = FALSE)
    if (is.null(ct_ind_col))
      stop("ctdna_make_patient_data: `Cancer_Type` frame has no recognizable ",
            "indication column. Looked for (case/punct insensitive): ",
            "Cancer_Type, Cancertype, Indication, Clin_Indication. Got columns: ",
            paste(names(Cancer_Type), collapse = ", "), ".", call. = FALSE)
    ct_ids  <- .norm_pid(Cancer_Type[[ct_id_col]])
    ct_vals <- as.character(Cancer_Type[[ct_ind_col]])
    ct_vals <- trimws(ct_vals)
    ct_vals[!nzchar(ct_vals) | ct_vals == "NA"] <- NA_character_

    # Drop duplicates in Cancer_Type frame (first occurrence wins)
    dup_mask <- duplicated(ct_ids)
    if (any(dup_mask)) {
      message(sprintf(
        "ctdna_make_patient_data: %d duplicate Patient_ID(s) in `Cancer_Type` frame; kept first occurrence.",
        sum(dup_mask)))
      ct_ids  <- ct_ids[!dup_mask]
      ct_vals <- ct_vals[!dup_mask]
    }
  } else {
    # Will derive ct_ids/ct_vals from ADSL below after we have frames classified
    ct_ids  <- character(0)
    ct_vals <- character(0)
    if (isTRUE(verbose))
      message("ctdna_make_patient_data: Cancer_Type not provided; will derive ",
              "from ADSL (priority: PRTUMTY -> COHORT -> STUDYTRT).")
  }

  # --- Validate ... frames --------------------------------------------------
  frames <- list(...)
  keep   <- vapply(frames, function(x) is.data.frame(x) && nrow(x) > 0L,
                    logical(1))
  frames <- frames[keep]
  if (length(frames) == 0L)
    stop("ctdna_make_patient_data: pass at least one non-empty data frame ",
          "in `...` (e.g. adsl/adrs/adtr). Cancer_Type alone is not enough.",
          call. = FALSE)
  frames <- lapply(frames, .normalize_columns)

  # Classify each frame
  frames_meta <- lapply(seq_along(frames), function(i) {
    df <- frames[[i]]
    list(idx = i, df = df,
         kind = .classify_clinical_frame(df),
         id_col = .find_id_col(df))
  })

  # Drop frames with no recognizable ID column
  no_id <- vapply(frames_meta, function(fm) is.null(fm$id_col), logical(1))
  if (any(no_id) && !all(no_id)) {
    # Some frames lack ID but others have it: warn and continue with the
    # ones that work.
    warning(sprintf(
      "ctdna_make_patient_data: %d input frame(s) had no recognizable patient ID column; ignored.",
      sum(no_id)), call. = FALSE)
    frames_meta <- frames_meta[!no_id]
  }
  if (length(frames_meta) == 0L || all(no_id)) {
    # v0.38.0: clear error path. No frame had an ID column we recognized,
    # so we tell the user exactly what to provide.
    id_priority <- .co("id_col_priority")
    cols_seen <- unique(unlist(lapply(list(...), names)))
    stop(
      "ctdna_make_patient_data: no input frame has a recognizable patient ID column.\n\n",
      "  Tried (in order, case-sensitive):\n    ",
      paste(id_priority, collapse = ", "), "\n\n",
      "  Columns seen across the input frame(s):\n    ",
      paste(utils::head(cols_seen, 30), collapse = ", "),
      if (length(cols_seen) > 30L) ", ..." else "", "\n\n",
      "To fix this, please provide:\n",
      "  1. The patient ID column name in your data (e.g. 'PT_NUM', 'subject').\n",
      "     Set it via:\n",
      "        clinical_opts(id_col_priority = c('YOUR_ID_COL', ",
      paste(shQuote(id_priority), collapse = ", "), "))\n",
      "  2. The cohort column name (e.g. 'COHORT', 'Cohort_Name').\n",
      "     Set it via:\n",
      "        ctdna_opts(cohort = 'YOUR_COHORT_COL')\n",
      "  3. The dose column name (e.g. 'DOSE', 'dose_level').\n",
      "     Set it via:\n",
      "        ctdna_opts(dose = 'YOUR_DOSE_COL')\n\n",
      "Then re-run ctdna_make_patient_data().",
      call. = FALSE)
  }

  # Sort ADaM frames first (priority on conflicts)
  kind_order <- c(adsl = 1L, adrs = 2L, adtr = 3L, custom = 4L)
  frames_meta <- frames_meta[order(vapply(frames_meta,
                                             function(fm) {
                                               v <- kind_order[fm$kind]
                                               if (is.na(v)) 4L else as.integer(v)
                                             },
                                             integer(1)))]

  if (isTRUE(verbose)) {
    counts <- table(vapply(frames_meta, `[[`, character(1), "kind"))
    message(sprintf(
      "ctdna_make_patient_data: %d input frame(s): %s",
      length(frames_meta),
      paste(sprintf("%s = %d", names(counts), counts), collapse = ", ")))
  }

  # --- Assemble the patient set --------------------------------------------
  # Patient universe: union of IDs in Cancer_Type + IDs in `...` frames
  ids_from_frames <- unique(unlist(lapply(frames_meta, function(fm)
    .norm_pid(fm$df[[fm$id_col]])), use.names = FALSE))
  all_ids <- unique(c(ct_ids, ids_from_frames))
  all_ids <- all_ids[!is.na(all_ids) & nzchar(all_ids)]

  # v0.38.0: when Cancer_Type wasn't supplied, derive it from ADSL by
  # walking a priority list of indication columns. Falls back to whichever
  # column exists first in this order: PRTUMTY, COHORT, STUDYTRT.
  if (use_adam_indication) {
    indication_priority <- .o("col_indication_adam") %||%
                            c("PRTUMTY","COHORT","STUDYTRT")
    derived <- .extract_canonical("_derived_indication", indication_priority,
                                    frames_meta, all_ids)
    derived <- trimws(derived)
    derived[!nzchar(derived) | derived == "NA"] <- NA_character_
    Cancer_Type_vec <- derived
    if (isTRUE(verbose)) {
      n_with <- sum(!is.na(Cancer_Type_vec))
      message(sprintf(
        "ctdna_make_patient_data: derived Cancer_Type from ADSL for %d/%d patient(s) using priority [%s].",
        n_with, length(all_ids), paste(indication_priority, collapse = " -> ")))
    }
  } else {
    # Pull Cancer_Type values aligned to all_ids
    Cancer_Type_vec <- ct_vals[match(all_ids, ct_ids)]
  }

  # --- BCR / Response (from adrs PARAM rows) -------------------------------
  bcr_param     <- .co("bcr_param")
  bcr_value_col <- .co("bcr_value_col")
  recist_from_adrs <- function(df, id_col, pat_ids) {
    if (!"PARAM" %in% names(df) || !bcr_value_col %in% names(df))
      return(NULL)
    sel <- df[!is.na(df$PARAM) & df$PARAM == bcr_param, , drop = FALSE]
    if (nrow(sel) == 0L) return(NULL)
    if ("ABLFL" %in% names(sel))
      sel <- sel[is.na(sel$ABLFL) | sel$ABLFL %in% c("Y",""), , drop = FALSE]
    sel <- sel[!duplicated(.norm_pid(sel[[id_col]])), , drop = FALSE]
    idx <- match(pat_ids, .norm_pid(sel[[id_col]]))
    as.character(sel[[bcr_value_col]])[idx]
  }

  raw_response <- .extract_canonical(
    "_raw_response", c("RECIST","BOR","Best_Response","Response_Subcategory"),
    frames_meta, all_ids,
    special_handlers = list(adrs = recist_from_adrs))

  # --- Dose -----------------------------------------------------------------
  Dose <- .extract_canonical("Dose", .dose_priority, frames_meta, all_ids)
  if (all(is.na(Dose))) {
    arm_vals <- .extract_canonical("_dose_via_arm", .arm_priority,
                                       frames_meta, all_ids)
    Dose <- .extract_dose_from_arm(arm_vals)
  }

  # --- ARM ------------------------------------------------------------------
  ARM <- .extract_canonical("ARM", .arm_priority, frames_meta, all_ids)

  # --- v0.38.0: expanded canonical columns from ADaM -----------------------
  # All derived from the ADSL columns the user listed (and similar in
  # ADRS/ADTR for backward compat). Each canonical column is filled from
  # the first source frame that has a non-NA value for that patient.
  arm_code     <- .extract_canonical("ARMCD",
                                       c("ARMCD","ACTARMCD"),
                                       frames_meta, all_ids)
  Sex          <- .extract_canonical("SEX",
                                       c("SEX","Sex"),
                                       frames_meta, all_ids)
  Age          <- .extract_canonical("AGE",
                                       c("AGE","AAGE1","AAGE2","Age"),
                                       frames_meta, all_ids)
  Race         <- .extract_canonical("RACE",
                                       c("RACE","ARACE","Race"),
                                       frames_meta, all_ids)
  Cohort       <- .extract_canonical("COHORT",
                                       c("COHORT","Cohort"),
                                       frames_meta, all_ids)
  TRTSDT       <- .extract_canonical("TRTSDT",
                                       c("TRTSDT","Treatment_start_date"),
                                       frames_meta, all_ids)
  TRTEDT       <- .extract_canonical("TRTEDT",
                                       c("TRTEDT","Treatment_end_date"),
                                       frames_meta, all_ids)
  DTHFL        <- .extract_canonical("DTHFL",
                                       c("DTHFL","Death_flag"),
                                       frames_meta, all_ids)
  DTHDT        <- .extract_canonical("DTHDT",
                                       c("DTHDT","Death_date"),
                                       frames_meta, all_ids)

  # --- Build canonical frame ------------------------------------------------
  out <- data.frame(
    Patient_ID  = all_ids,
    Dose        = Dose,
    Cancer_Type = Cancer_Type_vec,
    ARM         = ARM,
    stringsAsFactors = FALSE)
  out$Response_Subcategory <- .recist_to_subcategory(raw_response)
  out$Response             <- .subcategory_to_response(out$Response_Subcategory)

  # v0.38.0: append the broader canonical column set
  out$arm_code             <- arm_code
  out$Sex                  <- Sex
  out$Age                  <- suppressWarnings(as.numeric(Age))
  out$Race                 <- Race
  out$Cohort               <- Cohort
  out$Treatment_start_date <- TRTSDT
  out$Treatment_end_date   <- TRTEDT
  out$Death_flag           <- DTHFL
  out$Death_date           <- DTHDT

  # --- Extra cols -----------------------------------------------------------
  if (!is.null(extra_cols)) {
    if (!is.character(extra_cols))
      stop("`extra_cols` must be a character vector.", call. = FALSE)
    for (ec in unique(extra_cols)) {
      vals <- .extract_canonical(ec, ec, frames_meta, all_ids)
      out[[ec]] <- vals
      if (all(is.na(vals)) && isTRUE(verbose))
        message(sprintf(
          "ctdna_make_patient_data: extra_col '%s' has no source in any input frame; included as NA.",
          ec))
    }
  }

  # --- Drop rows with empty Cancer_Type OR empty Dose -----------------------
  .is_missing <- function(x) {
    s <- trimws(as.character(x))
    is.na(s) | !nzchar(s) |
      toupper(s) %in% c("NA","NOT AVAILABLE","UNKNOWN","N/A","NOT APPLICABLE")
  }
  ct_missing   <- .is_missing(out$Cancer_Type)
  dose_missing <- .is_missing(out$Dose)

  if (any(ct_missing)) {
    n_drop <- sum(ct_missing)
    ids <- out$Patient_ID[ct_missing]
    preview <- paste(utils::head(ids, 5), collapse = ", ")
    if (length(ids) > 5L) preview <- paste0(preview, ", ...")
    message(sprintf(
      "ctdna_make_patient_data: dropped %d patient(s) with empty Cancer_Type. IDs: %s",
      n_drop, preview))
  }
  if (any(dose_missing)) {
    n_drop <- sum(dose_missing)
    ids <- out$Patient_ID[dose_missing]
    preview <- paste(utils::head(ids, 5), collapse = ", ")
    if (length(ids) > 5L) preview <- paste0(preview, ", ...")
    message(sprintf(
      "ctdna_make_patient_data: dropped %d patient(s) with empty Dose. IDs: %s",
      n_drop, preview))
  }
  drop_mask <- ct_missing | dose_missing
  out <- out[!drop_mask, , drop = FALSE]
  rownames(out) <- NULL

  if (nrow(out) == 0L)
    stop("ctdna_make_patient_data: every patient was dropped (empty ",
          "Cancer_Type or empty Dose). Check your inputs.", call. = FALSE)

  out <- structure(out, class = c("ctdna_patient_data", class(out)))

  if (isTRUE(verbose)) {
    message(sprintf(
      "ctdna_make_patient_data: %d patient(s) in final output (%d columns).",
      nrow(out), ncol(out)))
    cat("\nFirst rows of the assembled patient_data:\n")
    print(utils::head(out))

    ct_table <- table(out$Cancer_Type, useNA = "ifany")
    if (length(ct_table) > 0L) {
      cat("\nCancer_Type values (with patient counts):\n")
      ct_df <- data.frame(
        Cancer_Type = names(ct_table),
        N           = as.integer(ct_table),
        stringsAsFactors = FALSE)
      ct_df <- ct_df[order(-ct_df$N), , drop = FALSE]
      rownames(ct_df) <- NULL
      print(ct_df)
    }
  }

  out
}


# ---------------------------------------------------------------------------
# Internal: resolve a patient_data arg passed to plot functions.
# Accepts NULL / data frame / ctdna_patient_data.
# A list of frames is REJECTED in v0.35.0+ -- users must call
# ctdna_make_patient_data() explicitly first because Cancer_Type is required.
# ---------------------------------------------------------------------------
.resolve_patient_data <- function(patient_data) {
  if (is.null(patient_data)) return(NULL)
  if (inherits(patient_data, "ctdna_patient_data")) return(patient_data)
  if (is.data.frame(patient_data)) return(patient_data)
  if (is.list(patient_data))
    stop("ctdna_oncoprint: passing a list of data frames as ",
          "`patient_data` is no longer supported (v0.35.0+). The patient ",
          "frame must be built explicitly first:\n\n",
          "  patient_df <- ctdna_make_patient_data(\n",
          "    adsl, adrs, adtr,\n",
          "    Cancer_Type = Master_Clinical[c(\"Patient_ID\",\"Cancer_Type\")])\n",
          "  ctdna_oncoprint(geno, patient_data = patient_df, ...)\n",
          call. = FALSE)
  stop("ctdna_oncoprint: `patient_data` must be NULL, a data frame, ",
        "or a ctdna_patient_data object.", call. = FALSE)
}
