# ============================================================================
# Multi-file readers for orthogonal modalities (dnaseq / rnaseq)  -- v0.75.10
# ----------------------------------------------------------------------------
# When each sample's variant calls live in a separate TSV under a per-sample
# directory tree, users can point ctdna_prepare() at the tree via a two-piece
# spec:
#
#   dnaseq = list(loc = "/mnt/etc1/etc2/",
#                 regex = "*/variant/annotated/*.clinical.annotated.tsv")
#
# Glob-style pattern. Dots are literal, `*` is any character except `/`.
# The FIRST `*` in the FILENAME becomes Sample_ID.
# ============================================================================

.ctdna_dnaseq_defaults <- function() {
  list(
    dnaseq_loc                  = NULL,
    dnaseq_regex                = NULL,
    dnaseq_gene_col             = "Gene Symbol",
    dnaseq_variant_type_col     = "Variant Type",
    dnaseq_vaf_col              = "Allelic Fraction",
    dnaseq_variant_effect_col   = "Variant Effect",
    dnaseq_protein_variant_col  = "Protein Variant",
    dnaseq_functional_class_col = "Functional Class",
    dnaseq_somatic_status_col   = "Somatic_status",
    dnaseq_qc_filter_col        = "Filter",
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

#' Read a directory tree of per-sample delimited files (dnaseq / rnaseq)
#'
#' Walks a directory root matching a glob pattern, reads every file, extracts
#' Sample_ID from the filename wildcard, canonicalises columns via
#' \code{ctdna_opts()} keys, applies optional QC filter and min-VAF, and
#' concatenates. Used internally by \code{ctdna_prepare()} when the user
#' supplies \code{dnaseq = list(loc = ..., regex = ...)} (or the equivalent
#' via \code{ctdna_opts(dnaseq_loc, dnaseq_regex)}).
#'
#' Glob rules: dots are literal, \code{*} matches any character except path
#' separator. The FIRST \code{*} in the filename (last path segment) becomes
#' the Sample_ID.
#'
#' @param loc      Directory root (character path).
#' @param regex    Glob pattern relative to \code{loc}.
#' @param modality \code{"dnaseq"} (default) or \code{"rnaseq"} -- selects
#'   which \code{ctdna_opts()} keys drive canonicalisation.
#' @param verbose  Show a progress bar and per-file status. Default
#'   \code{TRUE}.
#'
#' @return A \code{data.frame} with at minimum Sample_ID, Patient_ID, Gene,
#'   plus canonical modality columns.
#'
#' @keywords internal
.ctdna_read_delimited_tree <- function(loc, regex,
                                        modality = c("dnaseq","rnaseq"),
                                        verbose = TRUE) {
  modality <- match.arg(modality)
  if (!dir.exists(loc))
    stop(".ctdna_read_delimited_tree: loc does not exist: ", loc, call. = FALSE)

  full_glob <- file.path(loc, regex)
  files <- Sys.glob(full_glob)
  if (!length(files))
    stop(".ctdna_read_delimited_tree: no files match ", full_glob, call. = FALSE)

  if (verbose)
    message(sprintf("[%s reader] found %d file(s) under %s",
                     modality, length(files), loc))

  fname_pat <- basename(regex)
  if (modality == "dnaseq") {
    o <- .ctdna_dnaseq_defaults()
    for (k in names(o)) o[[k]] <- tryCatch(
      ctdna_opts(k, default = o[[k]]),
      error = function(e) o[[k]])
    gene_col        <- o$dnaseq_gene_col
    var_type_col    <- o$dnaseq_variant_type_col
    vaf_col         <- o$dnaseq_vaf_col
    var_effect_col  <- o$dnaseq_variant_effect_col
    prot_var_col    <- o$dnaseq_protein_variant_col
    func_class_col  <- o$dnaseq_functional_class_col
    somatic_col     <- o$dnaseq_somatic_status_col
    qc_filter_col   <- o$dnaseq_qc_filter_col
    min_vaf         <- as.numeric(o$dnaseq_min_vaf)
    if (!is.finite(min_vaf) || min_vaf < 0 || min_vaf > 1) {
      warning("dnaseq_min_vaf out of range [0,1]; using 0.", call. = FALSE)
      min_vaf <- 0
    }
  } else {
    o <- .ctdna_rnaseq_defaults()
    for (k in names(o)) o[[k]] <- tryCatch(
      ctdna_opts(k, default = o[[k]]),
      error = function(e) o[[k]])
    gene_col       <- o$rnaseq_gene_col
    val_col        <- o$rnaseq_value_col
    qc_filter_col  <- o$rnaseq_qc_filter_col
  }

  pb <- if (isTRUE(verbose))
          utils::txtProgressBar(min = 0, max = length(files), style = 3)
        else NULL

  frames <- vector("list", length(files))
  n_kept <- 0L; n_dropped_qc <- 0L; n_dropped_vaf <- 0L
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
    if (modality == "dnaseq") {
      rn(gene_col,       "Gene")
      rn(var_type_col,   "Variant_type")
      rn(vaf_col,        "Allelic_Fraction")
      rn(var_effect_col, "Variant_Effect")
      rn(prot_var_col,   "Protein_Variant")
      rn(func_class_col, "Functional_Class")
      rn(somatic_col,    "Somatic_status")
      rn(qc_filter_col,  "Filter")
    } else {
      rn(gene_col,      "Gene")
      rn(val_col,       "Value")
      rn(qc_filter_col, "Filter")
    }

    if ("Filter" %in% names(df)) {
      pass <- toupper(as.character(df$Filter)) %in% c("PASS","")
      n_dropped_qc <- n_dropped_qc + sum(!pass)
      df <- df[pass, , drop = FALSE]
    }
    if (modality == "dnaseq" && "Allelic_Fraction" %in% names(df) &&
        min_vaf > 0) {
      vaf <- suppressWarnings(as.numeric(df$Allelic_Fraction))
      is_pct <- any(vaf > 1, na.rm = TRUE)
      if (is_pct) vaf <- vaf / 100
      keep <- !is.na(vaf) & vaf >= min_vaf
      n_dropped_vaf <- n_dropped_vaf + sum(!keep, na.rm = TRUE)
      df <- df[keep, , drop = FALSE]
    }

    if (!nrow(df)) next
    frames[[i]] <- df
    n_kept <- n_kept + 1L
  }
  if (!is.null(pb)) close(pb)

  frames <- frames[!vapply(frames, is.null, logical(1))]
  if (!length(frames))
    stop(".ctdna_read_delimited_tree: every file was empty or failed. ",
          "First errors: ", paste(utils::head(errors, 3), collapse = "; "),
          call. = FALSE)

  out <- .ctdna_rbind_union(frames)
  if (verbose) {
    message(sprintf("[%s reader] kept %d file(s), %d row(s) total",
                     modality, n_kept, nrow(out)))
    if (n_dropped_qc)
      message(sprintf("[%s reader] dropped %d row(s) failing QC filter",
                       modality, n_dropped_qc))
    if (n_dropped_vaf)
      message(sprintf("[%s reader] dropped %d row(s) below min VAF %.3f",
                       modality, n_dropped_vaf, min_vaf))
    if (length(errors))
      warning(sprintf("[%s reader] %d file(s) failed; e.g.: %s",
                       modality, length(errors), errors[1]), call. = FALSE)
  }
  out
}
