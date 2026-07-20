# ============================================================================
# tree_reader.R — multi-file readers for dnaseq / rnaseq
# ----------------------------------------------------------------------------
# Two public entry points:
#
#   dnaseq_create()             — build a canonical dnaseq data.frame by
#                                  walking one or more directory trees, one
#                                  per alteration kind (SNV/Indel, CNV amp,
#                                  CNV del, LOH, Fusion, LGR, ...). Each kind
#                                  can have its own location and its own file
#                                  glob. Rows are stamped with the correct
#                                  Variant_type / CNV_type so the downstream
#                                  classifier + oncoprint pick them up
#                                  without any code changes.
#
#   .ctdna_read_delimited_tree()  — internal, RNAseq only. Kept so
#                                  rnaseq_loc / rnaseq_regex options-based
#                                  reading in ctdna_prepare() continues to
#                                  work. For dnaseq, this function is no
#                                  longer used — use dnaseq_create() instead.
#
# `ctdna_opts("dnaseq_loc")` and `ctdna_opts("dnaseq_regex")` have been
# removed in v0.76.0. All dnaseq tree reading now goes through
# dnaseq_create() (called directly, or via `ctdna_prepare(dnaseq = list(...))`
# which forwards the spec).
# ============================================================================

# ----- Column-canonicalisation defaults (unchanged behaviour) ---------------

.ctdna_dnaseq_defaults <- function() {
  list(
    dnaseq_gene_col             = "Gene Symbol",
    dnaseq_variant_type_col     = "Variant Type",
    dnaseq_vaf_col              = "Allelic Fraction",
    dnaseq_variant_effect_col   = "Variant Effect",
    dnaseq_protein_variant_col  = "Protein Variant",
    dnaseq_functional_class_col = "Functional Class",
    dnaseq_somatic_status_col   = "Somatic_status",
    dnaseq_qc_filter_col        = "Filter",
    dnaseq_cnv_type_col         = "CNV_type",
    dnaseq_min_vaf              = 0)
}
.ctdna_rnaseq_defaults <- function() {
  list(
    rnaseq_loc            = NULL,
    rnaseq_regex          = NULL,
    rnaseq_gene_col       = "Gene Symbol",
    rnaseq_value_col      = "TPM",
    rnaseq_qc_filter_col  = "Filter")
}

# ----- Low-level helpers (unchanged) ----------------------------------------

# Extract Sample_ID from filename given the filename portion of the glob.
# Example: pattern = "*.clinical.annotated.tsv",
# filename = "sampleA.clinical.annotated.tsv" -> "sampleA".
.ctdna_sample_id_from_filename <- function(filename, pattern) {
  esc <- gsub("([][{}().+^$|\\\\])", "\\\\\\1", pattern, perl = TRUE)
  esc <- sub("\\*", "(.*?)", esc)
  esc <- gsub("\\*", ".*", esc)
  m   <- regmatches(filename, regexec(paste0("^", esc, "$"), filename))
  if (length(m) && length(m[[1]]) >= 2) m[[1]][2] else NA_character_
}

# Robust delimited reader (TSV then CSV fallback).
.ctdna_read_delim <- function(path) {
  txt <- tryCatch(
    utils::read.table(path, sep = "\t", header = TRUE,
                       stringsAsFactors = FALSE, check.names = FALSE,
                       na.strings = c("", "NA", "NaN", "."),
                       comment.char = "", quote = "\""),
    error = function(e) NULL)
  if (!is.null(txt) && ncol(txt) > 1) return(txt)
  txt <- tryCatch(
    utils::read.table(path, sep = ",", header = TRUE,
                       stringsAsFactors = FALSE, check.names = FALSE,
                       na.strings = c("", "NA", "NaN", "."),
                       comment.char = "", quote = "\""),
    error = function(e) NULL)
  if (!is.null(txt) && ncol(txt) > 1) return(txt)
  stop(sprintf("could not parse %s as TSV or CSV", path), call. = FALSE)
}

# Union-rbind: pad each frame with missing columns as NA, then rbind.
.ctdna_rbind_union <- function(frames) {
  cols_all <- unique(unlist(lapply(frames, names)))
  frames <- lapply(frames, function(d) {
    miss <- setdiff(cols_all, names(d))
    for (m in miss) d[[m]] <- NA
    d[, cols_all, drop = FALSE]
  })
  do.call(rbind, frames)
}

# ----- Per-kind file walker -------------------------------------------------
# Reads every file under `loc` matching glob `pattern`, canonicalises the
# dnaseq column names via ctdna_opts(), applies optional Filter=="PASS" and
# min-VAF drops, and returns one data.frame. `Sample_ID` and `Patient_ID`
# are set from the wildcard capture in the filename.
.read_one_dnaseq_tree <- function(loc, pattern, min_vaf, verbose = TRUE) {
  if (is.null(loc) || !nzchar(loc))
    stop("dnaseq_create: `loc` is empty for pattern '", pattern, "'.",
         call. = FALSE)
  if (!dir.exists(loc))
    stop("dnaseq_create: loc does not exist: ", loc, call. = FALSE)

  full_glob <- file.path(loc, pattern)
  files <- Sys.glob(full_glob)
  if (!length(files))
    stop("dnaseq_create: no files match ", full_glob, call. = FALSE)

  if (isTRUE(verbose))
    message(sprintf("[dnaseq_create] %d file(s) match %s",
                     length(files), full_glob))

  o <- .ctdna_dnaseq_defaults()
  for (k in names(o)) o[[k]] <- tryCatch(
    ctdna_opts(k, default = o[[k]]),
    error = function(e) o[[k]])
  gene_col       <- o$dnaseq_gene_col
  var_type_col   <- o$dnaseq_variant_type_col
  vaf_col        <- o$dnaseq_vaf_col
  var_effect_col <- o$dnaseq_variant_effect_col
  prot_var_col   <- o$dnaseq_protein_variant_col
  func_class_col <- o$dnaseq_functional_class_col
  somatic_col    <- o$dnaseq_somatic_status_col
  qc_filter_col  <- o$dnaseq_qc_filter_col
  cnv_type_col   <- o$dnaseq_cnv_type_col

  fname_pat <- basename(pattern)
  pb <- if (isTRUE(verbose))
          utils::txtProgressBar(min = 0, max = length(files), style = 3)
        else NULL

  frames <- vector("list", length(files))
  errors <- character(0)
  n_dropped_qc  <- 0L
  n_dropped_vaf <- 0L
  for (i in seq_along(files)) {
    f <- files[i]
    if (!is.null(pb)) utils::setTxtProgressBar(pb, i)

    sid <- .ctdna_sample_id_from_filename(basename(f), fname_pat)
    if (is.na(sid) || !nzchar(sid))
      sid <- tools::file_path_sans_ext(basename(f))

    df <- tryCatch(.ctdna_read_delim(f), error = function(e) {
      errors <<- c(errors, sprintf("%s: %s", f, conditionMessage(e)))
      NULL
    })
    if (is.null(df) || !nrow(df)) next

    df$Sample_ID  <- sid
    df$Patient_ID <- sid

    rn <- function(from, to) {
      if (from %in% names(df) && !(to %in% names(df)))
        names(df)[names(df) == from] <<- to
    }
    rn(gene_col,       "Gene")
    rn(var_type_col,   "Variant_type")
    rn(vaf_col,        "Allelic_Fraction")
    rn(var_effect_col, "Variant_Effect")
    rn(prot_var_col,   "Protein_Variant")
    rn(func_class_col, "Functional_Class")
    rn(somatic_col,    "Somatic_status")
    rn(qc_filter_col,  "Filter")
    rn(cnv_type_col,   "CNV_type")

    if ("Filter" %in% names(df)) {
      pass <- toupper(as.character(df$Filter)) %in% c("PASS","")
      n_dropped_qc <- n_dropped_qc + sum(!pass)
      df <- df[pass, , drop = FALSE]
    }
    if ("Allelic_Fraction" %in% names(df) && min_vaf > 0) {
      vaf <- suppressWarnings(as.numeric(df$Allelic_Fraction))
      is_pct <- any(vaf > 1, na.rm = TRUE)
      if (is_pct) vaf <- vaf / 100
      keep <- !is.na(vaf) & vaf >= min_vaf
      n_dropped_vaf <- n_dropped_vaf + sum(!keep, na.rm = TRUE)
      df <- df[keep, , drop = FALSE]
    }
    if (!nrow(df)) next
    frames[[i]] <- df
  }
  if (!is.null(pb)) close(pb)

  frames <- frames[!vapply(frames, is.null, logical(1))]
  if (!length(frames))
    stop("dnaseq_create: every file was empty or failed under ", full_glob,
         ". First errors: ", paste(utils::head(errors, 3), collapse = "; "),
         call. = FALSE)

  out <- .ctdna_rbind_union(frames)
  if (isTRUE(verbose)) {
    message(sprintf("[dnaseq_create]   -> %d row(s) total", nrow(out)))
    if (n_dropped_qc)
      message(sprintf("[dnaseq_create]   -> dropped %d row(s) failing QC",
                       n_dropped_qc))
    if (n_dropped_vaf)
      message(sprintf("[dnaseq_create]   -> dropped %d row(s) below min VAF %.3f",
                       n_dropped_vaf, min_vaf))
    if (length(errors))
      warning(sprintf("[dnaseq_create] %d file(s) failed; e.g.: %s",
                       length(errors), errors[1]), call. = FALSE)
  }
  out
}

# Resolve a single regex_* argument to a normalised (loc, pattern) pair.
# Accepts: NULL, a bare string, or a list with `loc` and `pattern` keys.
.resolve_dnaseq_spec <- function(spec, default_loc, kind) {
  if (is.null(spec)) return(NULL)
  if (is.character(spec) && length(spec) == 1L && nzchar(spec)) {
    if (is.null(default_loc) || !nzchar(default_loc))
      stop("dnaseq_create: regex_", kind, " was given as a bare string but ",
           "top-level `loc` is missing. Either supply `loc = ...`, or pass ",
           "regex_", kind, " = list(loc = ..., pattern = ...).",
           call. = FALSE)
    return(list(loc = default_loc, pattern = spec))
  }
  if (is.list(spec) && all(c("loc","pattern") %in% names(spec))) {
    return(list(loc = as.character(spec$loc),
                pattern = as.character(spec$pattern)))
  }
  stop("dnaseq_create: regex_", kind, " must be NULL, a single character ",
       "string, or a list(loc = ..., pattern = ...). Got: ",
       paste(class(spec), collapse = "/"), ".", call. = FALSE)
}

# Registry of every kind this reader knows about, plus what to stamp on rows
# coming from that kind. Stamps are OVERWRITES (they win over any value the
# file itself may have carried), because if a user explicitly points regex_amp
# at an AMP file, they are asserting the kind.
.dnaseq_kind_registry <- function() {
  list(
    all         = list(stamp = NULL),
    snv         = list(stamp = NULL),
    cnv         = list(stamp = list(Variant_type = "CNV")),
    amp         = list(stamp = list(Variant_type = "CNV",
                                     CNV_type     = "amplification")),
    focal_amp   = list(stamp = list(Variant_type = "CNV",
                                     CNV_type     = "focal_amplification")),
    del         = list(stamp = list(Variant_type = "CNV",
                                     CNV_type     = "homozygous_deletion")),
    loh         = list(stamp = list(Variant_type = "CNV",
                                     CNV_type     = "loh_deletion")),
    fusion      = list(stamp = list(Variant_type = "Fusion")),
    lgr         = list(stamp = list(Variant_type = "LGR")))
}

#' Assemble a canonical dnaseq data.frame from one or more file trees
#'
#' @description
#' Reads variant-call files from disk and returns a single tidy data.frame
#' ready to pass to \code{\link{ctdna_prepare}} via its \code{dnaseq =}
#' argument. Handles the three common on-disk layouts:
#'
#' \enumerate{
#'   \item \strong{Combined} — everything (SNV, Indel, CNV, Fusion, LGR)
#'         lives in one file per sample. Use \code{regex_all} (or
#'         \code{regex_snv} if the file only carries small variants).
#'   \item \strong{Split by kind} — one file per sample for small variants,
#'         one file per sample for all CNV together. Use \code{regex_snv}
#'         and \code{regex_cnv}.
#'   \item \strong{Split by event} — one file per sample per event kind
#'         (e.g. \code{sampleA_AMP.tsv}, \code{sampleA_DEL.tsv},
#'         \code{sampleA_LOH.tsv}). Use \code{regex_amp}, \code{regex_del},
#'         \code{regex_loh}, \code{regex_focal_amp}, \code{regex_fusion},
#'         \code{regex_lgr} as needed.
#' }
#'
#' Each \code{regex_*} argument can be supplied in two shapes:
#' \itemize{
#'   \item A single glob pattern string — the file walk is anchored at the
#'         top-level \code{loc}.
#'   \item A \code{list(loc = "...", pattern = "...")} — the kind uses its
#'         own directory, ignoring the top-level \code{loc}. Different kinds
#'         may live under completely different trees.
#' }
#'
#' Rows coming from a per-event regex are stamped with the correct
#' \code{Variant_type} and \code{CNV_type} on the way in, so no downstream
#' code (classifier, oncoprint, filters) needs any changes — the values match
#' exactly what the internal \code{.oncoprint_classify()} helper expects.
#'
#' Glob rules: dots are literal, \code{*} matches any character except the
#' path separator. The FIRST \code{*} in the filename portion becomes
#' \code{Sample_ID}. If no wildcard is present, \code{Sample_ID} falls back
#' to the filename with the extension stripped.
#'
#' @param loc         Character path. Optional top-level default directory
#'   used by any \code{regex_*} passed as a bare string. If every
#'   \code{regex_*} carries its own \code{loc}, this can be left \code{NULL}.
#' @param regex_all   Glob or \code{list(loc, pattern)}. File(s) contain
#'   every kind of alteration — no stamping is done, rows pass through as-is.
#' @param regex_snv   Glob or \code{list(loc, pattern)}. File(s) contain
#'   SNV/Indel calls only. No stamping (Variant_type must be set in-file
#'   or by a preceding annotation step).
#' @param regex_cnv   Glob or \code{list(loc, pattern)}. File(s) contain all
#'   CNV events together. Rows are stamped with \code{Variant_type = "CNV"};
#'   \code{CNV_type} must come from a column in the file.
#' @param regex_amp   Glob or \code{list(loc, pattern)}. Amplification
#'   events. Rows stamped with \code{Variant_type = "CNV"}, \code{CNV_type =
#'   "amplification"}.
#' @param regex_focal_amp Glob or \code{list(loc, pattern)}. Focal
#'   amplification events. Rows stamped with \code{Variant_type = "CNV"},
#'   \code{CNV_type = "focal_amplification"}.
#' @param regex_del   Glob or \code{list(loc, pattern)}. Homozygous deletion
#'   events. Rows stamped with \code{Variant_type = "CNV"}, \code{CNV_type =
#'   "homozygous_deletion"}.
#' @param regex_loh   Glob or \code{list(loc, pattern)}. Loss-of-heterozygosity
#'   deletion events. Rows stamped with \code{Variant_type = "CNV"},
#'   \code{CNV_type = "loh_deletion"}.
#' @param regex_fusion Glob or \code{list(loc, pattern)}. Fusion events.
#'   Rows stamped with \code{Variant_type = "Fusion"}.
#' @param regex_lgr   Glob or \code{list(loc, pattern)}. Large genomic
#'   rearrangement events. Rows stamped with \code{Variant_type = "LGR"}.
#' @param verbose     Logical. Show per-file progress and summary counts.
#'   Default \code{TRUE}.
#'
#' @return A \code{data.frame} with at minimum \code{Sample_ID},
#'   \code{Patient_ID}, \code{Gene}, \code{Variant_type} — plus any of
#'   \code{CNV_type}, \code{Molecular_consequence}, \code{Allelic_Fraction},
#'   \code{Variant_Effect}, \code{Protein_Variant}, \code{Functional_Class},
#'   \code{Somatic_status} that were present in the source files. Column
#'   names come from \code{ctdna_opts()} (\code{dnaseq_gene_col},
#'   \code{dnaseq_variant_type_col}, etc.).
#'
#' @seealso \code{\link{ctdna_prepare}} for using the resulting frame,
#'   \code{\link{ctdna_opts}} for column-name customisation.
#'
#' @examples
#' \dontrun{
#'   # Mode 1: everything in one file per sample under a single tree
#'   d <- dnaseq_create(
#'     loc       = "/study/dnaseq",
#'     regex_all = "*.annotated.tsv")
#'
#'   # Mode 2: SNV/Indel + all-CNV, same tree
#'   d <- dnaseq_create(
#'     loc       = "/study/dnaseq",
#'     regex_snv = "*_variants.tsv",
#'     regex_cnv = "*_cnv.tsv")
#'
#'   # Mode 3: per-event files, mixed locations
#'   d <- dnaseq_create(
#'     loc          = "/study/dnaseq/snv",
#'     regex_snv    = "*_snv_indel.tsv",
#'     regex_amp    = list(loc = "/study/dnaseq/cnv",
#'                         pattern = "*_AMP.tsv"),
#'     regex_del    = list(loc = "/study/dnaseq/cnv",
#'                         pattern = "*_DEL.tsv"),
#'     regex_loh    = list(loc = "/study/dnaseq/cnv",
#'                         pattern = "*_LOH.tsv"),
#'     regex_fusion = list(loc = "/study/dnaseq/fusion",
#'                         pattern = "*.fusion.txt"),
#'     regex_lgr    = list(loc = "/study/dnaseq/lgr",
#'                         pattern = "*_LGR.tsv"))
#'
#'   # Then hand to ctdna_prepare():
#'   prep <- ctdna_prepare(infinity_report = inf, dnaseq = d, adam = adsl)
#' }
#'
#' @export
dnaseq_create <- function(loc              = NULL,
                           regex_all       = NULL,
                           regex_snv       = NULL,
                           regex_cnv       = NULL,
                           regex_amp       = NULL,
                           regex_focal_amp = NULL,
                           regex_del       = NULL,
                           regex_loh       = NULL,
                           regex_fusion    = NULL,
                           regex_lgr       = NULL,
                           verbose         = TRUE) {

  reg <- .dnaseq_kind_registry()
  # min_vaf pulled from opts once (shared across all kinds).
  min_vaf <- tryCatch(as.numeric(ctdna_opts("dnaseq_min_vaf", default = 0)),
                       error = function(e) 0)
  if (!is.finite(min_vaf) || min_vaf < 0 || min_vaf > 1) {
    warning("dnaseq_min_vaf out of range [0,1]; using 0.", call. = FALSE)
    min_vaf <- 0
  }

  supplied <- list(
    all         = regex_all,
    snv         = regex_snv,
    cnv         = regex_cnv,
    amp         = regex_amp,
    focal_amp   = regex_focal_amp,
    del         = regex_del,
    loh         = regex_loh,
    fusion      = regex_fusion,
    lgr         = regex_lgr)
  # Drop any kind that wasn't provided.
  supplied <- supplied[!vapply(supplied, is.null, logical(1))]
  if (!length(supplied))
    stop("dnaseq_create: no regex_* patterns provided. At least one of ",
         "regex_all / regex_snv / regex_cnv / regex_amp / regex_focal_amp / ",
         "regex_del / regex_loh / regex_fusion / regex_lgr must be set.",
         call. = FALSE)

  frames <- list()
  for (kind in names(supplied)) {
    s <- .resolve_dnaseq_spec(supplied[[kind]], loc, kind)
    if (isTRUE(verbose))
      message(sprintf("[dnaseq_create] kind='%s' -> loc=%s pattern=%s",
                       kind, s$loc, s$pattern))
    df <- .read_one_dnaseq_tree(loc = s$loc, pattern = s$pattern,
                                 min_vaf = min_vaf, verbose = verbose)
    stamps <- reg[[kind]]$stamp
    if (!is.null(stamps)) {
      for (col_ in names(stamps)) df[[col_]] <- stamps[[col_]]
    }
    df$Source_kind <- kind
    frames[[length(frames) + 1L]] <- df
  }

  out <- .ctdna_rbind_union(frames)
  if (isTRUE(verbose))
    message(sprintf("[dnaseq_create] combined: %d row(s) from %d kind(s).",
                     nrow(out), length(frames)))
  out
}

# ============================================================================
# INTERNAL: RNAseq tree reader — kept for ctdna_opts(rnaseq_loc/regex) and
# for the list-spec branch in ctdna_prepare(). Not used by dnaseq anymore.
# ============================================================================

#' @keywords internal
.ctdna_read_delimited_tree <- function(loc, regex,
                                        modality = c("rnaseq","dnaseq"),
                                        verbose = TRUE) {
  modality <- match.arg(modality)
  if (identical(modality, "dnaseq")) {
    # Back-compat shim: forward to dnaseq_create() with regex_all so
    # any legacy caller keeps working.
    return(dnaseq_create(loc = loc, regex_all = regex, verbose = verbose))
  }

  if (!dir.exists(loc))
    stop(".ctdna_read_delimited_tree: loc does not exist: ", loc, call. = FALSE)
  full_glob <- file.path(loc, regex)
  files <- Sys.glob(full_glob)
  if (!length(files))
    stop(".ctdna_read_delimited_tree: no files match ", full_glob, call. = FALSE)

  if (verbose)
    message(sprintf("[rnaseq reader] found %d file(s) under %s",
                     length(files), loc))

  fname_pat <- basename(regex)
  o <- .ctdna_rnaseq_defaults()
  for (k in names(o)) o[[k]] <- tryCatch(
    ctdna_opts(k, default = o[[k]]),
    error = function(e) o[[k]])
  gene_col       <- o$rnaseq_gene_col
  val_col        <- o$rnaseq_value_col
  qc_filter_col  <- o$rnaseq_qc_filter_col

  pb <- if (isTRUE(verbose))
          utils::txtProgressBar(min = 0, max = length(files), style = 3)
        else NULL

  frames <- vector("list", length(files))
  n_dropped_qc <- 0L
  errors <- character(0)

  for (i in seq_along(files)) {
    f <- files[i]
    if (!is.null(pb)) utils::setTxtProgressBar(pb, i)

    sid <- .ctdna_sample_id_from_filename(basename(f), fname_pat)
    if (is.na(sid) || !nzchar(sid))
      sid <- tools::file_path_sans_ext(basename(f))

    df <- tryCatch(.ctdna_read_delim(f), error = function(e) {
      errors <<- c(errors, sprintf("%s: %s", f, conditionMessage(e)))
      NULL
    })
    if (is.null(df) || !nrow(df)) next

    df$Sample_ID  <- sid
    df$Patient_ID <- sid

    rn <- function(from, to) {
      if (from %in% names(df) && !(to %in% names(df)))
        names(df)[names(df) == from] <<- to
    }
    rn(gene_col,      "Gene")
    rn(val_col,       "Value")
    rn(qc_filter_col, "Filter")

    if ("Filter" %in% names(df)) {
      pass <- toupper(as.character(df$Filter)) %in% c("PASS","")
      n_dropped_qc <- n_dropped_qc + sum(!pass)
      df <- df[pass, , drop = FALSE]
    }
    if (!nrow(df)) next
    frames[[i]] <- df
  }
  if (!is.null(pb)) close(pb)

  frames <- frames[!vapply(frames, is.null, logical(1))]
  if (!length(frames))
    stop(".ctdna_read_delimited_tree: every file was empty or failed. ",
          "First errors: ", paste(utils::head(errors, 3), collapse = "; "),
          call. = FALSE)

  out <- .ctdna_rbind_union(frames)
  if (verbose) {
    message(sprintf("[rnaseq reader] kept %d file(s), %d row(s) total",
                     length(frames), nrow(out)))
    if (n_dropped_qc)
      message(sprintf("[rnaseq reader] dropped %d row(s) failing QC filter",
                       n_dropped_qc))
    if (length(errors))
      warning(sprintf("[rnaseq reader] %d file(s) failed; e.g.: %s",
                       length(errors), errors[1]), call. = FALSE)
  }
  out
}
