# ============================================================================
# ctdnaTM -- package onAttach hook
# ============================================================================

.onAttach <- function(libname, pkgname) {
  v <- utils::packageVersion(pkgname)
  packageStartupMessage(
    sprintf("ctdnaTM v%s (ctDNA Deliverables) - developed by Hamid Nikbakht", v))
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    packageStartupMessage(
      "Tip: install ComplexHeatmap for publication-quality oncoprint:",
      "\n  install.packages('BiocManager')",
      "\n  BiocManager::install('ComplexHeatmap')"
    )
  }
  packageStartupMessage("Type `?ctdnaTM` for an overview, `ctdna_opts()` for config.")

  # v0.42.2: optional auto-check for updates. Opt-in via env var so
  # library() never makes surprise network calls. Set
  # Sys.setenv(CTDNATM_AUTO_UPDATE_CHECK = "TRUE") in ~/.Rprofile to
  # enable. Manual check via ctdna_check_for_updates() works regardless.
  if (identical(tolower(Sys.getenv("CTDNATM_AUTO_UPDATE_CHECK")), "true")) {
    tryCatch({
      res <- ctdna_check_for_updates(quiet = TRUE)
      if (!is.null(res) && isFALSE(res$up_to_date)) {
        packageStartupMessage(sprintf(
          "ctdnaTM: a newer version is available (%s; you have %s).\n  Download: %s",
          res$latest, res$installed, res$download_url))
      }
    }, error = function(e) NULL)
  }
}


# ============================================================================
# Indication label matching -- v0.78.1
# ----------------------------------------------------------------------------
# Cohort labels arrive with inconsistent capitalisation and separators across
# ADaM cuts, vendor files and analyst scripts: "Nsq_NSCLC", "Nsq-NSCLC",
# "nsq nsclc", "NSQ.NSCLC" all mean the same cohort. `.ctdna_norm_ind()`
# reduces a label to a comparison key so any of those forms match.
#
# Normalisation: lower-case, then drop every character that is not a letter
# or a digit. So "Nsq-NSCLC" -> "nsqnsclc" and "ES_SCLC" -> "essclc".
# ============================================================================
.ctdna_norm_ind <- function(x) {
  x <- as.character(x)
  gsub("[^a-z0-9]", "", tolower(x))
}

# Case-insensitive indication-column lookup. Returns the ACTUAL column name as
# it appears in `nms`, or NA_character_ if none of the candidates are present.
.ctdna_find_ind_col <- function(nms, col_hint = NULL) {
  cands <- c(col_hint, .o("indication"), "Indication", "indication",
             "Cancertype", "Cancer_Type", "cancer_type", "cancertype")
  cands <- unique(cands[!is.na(cands) & nzchar(cands)])
  key   <- .ctdna_norm_ind(nms)
  for (cd in cands) {
    hit <- which(key == .ctdna_norm_ind(cd))
    if (length(hit)) return(nms[hit[1]])
  }
  NA_character_
}

# ============================================================================
# .ctdna_filter_prep_by_indication(prep, indications, col_hint = NULL)
# ----------------------------------------------------------------------------
# Shared cohort-restriction helper for plot functions whose first argument
# is a `ctdna_prep` object.
#
# `indications = NULL` (or empty) -> return prep unchanged, i.e. ALL
# indications are used. This is the documented default everywhere.
#
# When set, keep only patients whose indication label matches one of the given
# values, and cascade the filter to every prep slot that has a Patient_ID
# column. Matching is case-insensitive and separator-insensitive.
# ============================================================================
.ctdna_filter_prep_by_indication <- function(prep, indications,
                                              col_hint = NULL) {
  # NULL / empty -> all indications (no filtering at all)
  if (is.null(indications) || !length(indications)) return(prep)
  indications <- indications[!is.na(indications) & nzchar(as.character(indications))]
  if (!length(indications)) return(prep)
  if (!is.list(prep) || !is.data.frame(prep$clinical)) return(prep)

  ind_col <- .ctdna_find_ind_col(names(prep$clinical), col_hint)
  if (is.na(ind_col))
    stop("indication: no indication column found in prep$clinical. Looked for ",
         "`", .o("indication"), "`, Indication, Cancertype, Cancer_Type ",
         "(case-insensitive). Add one before passing `indication`.",
         call. = FALSE)

  have <- as.character(prep$clinical[[ind_col]])
  sel  <- .ctdna_norm_ind(have) %in% .ctdna_norm_ind(indications)
  keep_pids <- as.character(prep$clinical$Patient_ID)[sel]

  # Fail loudly and immediately on a zero match. Previously this warned and
  # continued with an empty prep, which surfaced later as a confusing
  # "Auto-baseline detection found no candidate visits" error.
  if (!length(keep_pids)) {
    avail <- sort(unique(have[!is.na(have)]))
    stop(sprintf(
      paste0("indication: none of {%s} matched any value in `%s`.\n",
             "  Available: %s\n",
             "  Matching ignores case and separators, so \"%s\" would match \"%s\".\n",
             "  Pass indication = NULL (the default) to use every indication."),
      paste(indications, collapse = ", "), ind_col,
      if (length(avail)) paste(avail, collapse = ", ") else "<none>",
      if (length(avail)) gsub("_", "-", tolower(avail[1])) else "nsq-nsclc",
      if (length(avail)) avail[1] else "Nsq_NSCLC"),
      call. = FALSE)
  }

  # Report any requested value that matched nothing, even if others did.
  unmatched <- indications[!(.ctdna_norm_ind(indications) %in%
                               .ctdna_norm_ind(have))]
  if (length(unmatched))
    warning(sprintf("indication: %s matched no patients and was ignored.",
                    paste(sprintf("`%s`", unmatched), collapse = ", ")),
            call. = FALSE)

  n_before <- sum(vapply(prep, function(d)
    if (is.data.frame(d) && "Patient_ID" %in% names(d)) nrow(d) else 0L,
    integer(1)))
  for (slot in names(prep)) {
    d <- prep[[slot]]
    if (is.data.frame(d) && "Patient_ID" %in% names(d))
      prep[[slot]] <- d[as.character(d$Patient_ID) %in% keep_pids, , drop = FALSE]
  }
  n_after <- sum(vapply(prep, function(d)
    if (is.data.frame(d) && "Patient_ID" %in% names(d)) nrow(d) else 0L,
    integer(1)))
  n_dropped <- n_before - n_after
  if (n_dropped > 0)
    message(sprintf("indication: kept %d patient(s) in {%s}; dropped %d row(s) across prep slots.",
                     length(keep_pids), paste(indications, collapse = ", "),
                     n_dropped))
  prep
}


# ============================================================================
# .ctdna_filter_df_by_indication(df, indications, col_hint = NULL)
# ----------------------------------------------------------------------------
# Shared cohort-restriction helper for plot functions whose first argument
# is a plain data.frame (the integration `ctdna_plot_vs_*` family).
#
# `indications = NULL` -> return df unchanged (all indications). Matching is
# case- and separator-insensitive, as for the prep helper.
# ============================================================================
.ctdna_filter_df_by_indication <- function(df, indications,
                                             col_hint = NULL) {
  if (is.null(indications) || !length(indications)) return(df)
  indications <- indications[!is.na(indications) & nzchar(as.character(indications))]
  if (!length(indications)) return(df)
  if (!is.data.frame(df)) return(df)

  in_df <- .ctdna_find_ind_col(names(df), col_hint)
  if (is.na(in_df)) {
    warning("indication: no indication column found in the data.frame ",
            "(looked for `", .o("indication"), "`, Indication, Cancertype, ",
            "Cancer_Type; case-insensitive). Ignoring `indication`.",
            call. = FALSE)
    return(df)
  }

  have <- as.character(df[[in_df]])
  keep <- .ctdna_norm_ind(have) %in% .ctdna_norm_ind(indications)
  if (!any(keep)) {
    avail <- sort(unique(have[!is.na(have)]))
    stop(sprintf(
      paste0("indication: none of {%s} matched any value in `%s`.\n",
             "  Available: %s\n",
             "  Matching ignores case and separators."),
      paste(indications, collapse = ", "), in_df,
      if (length(avail)) paste(avail, collapse = ", ") else "<none>"),
      call. = FALSE)
  }
  out <- df[keep, , drop = FALSE]
  message(sprintf("indication: kept %d of %d row(s) in {%s}.",
                   nrow(out), nrow(df), paste(indications, collapse = ", ")))
  out
}


# ============================================================================
# .ctdna_resolve_indication(indication, indications)
# ----------------------------------------------------------------------------
# v0.78.0: `indication` is the documented argument name across the package
# (mirroring the historical `cancer_type=` argument). `indications` is kept as
# a working alias so existing scripts do not break.
#
# Accepts a single string or a character vector. If both are supplied and they
# disagree, `indication` wins and a message is emitted.
# ============================================================================
.ctdna_resolve_indication <- function(indication = NULL, indications = NULL) {
  a <- indication; b <- indications
  if (is.null(a) || !length(a)) return(b)
  if (is.null(b) || !length(b)) return(a)
  if (!setequal(as.character(a), as.character(b)))
    message("Both `indication` and `indications` were supplied and differ; ",
            "using `indication`. (`indications` is a deprecated alias.)")
  a
}


# ============================================================================
# .ctdna_resolve_facet(facet, wrap, fn = NULL)
# ----------------------------------------------------------------------------
# v0.79.0: `facet` and `wrap` are two names for the same concept -- "split the
# output into separate panels by this column". Historically the plot family
# used `facet` and the landscape family (oncoprint / alteration_grid /
# concordance_oncoprint) used `wrap`, so users had to remember which name went
# with which function. Both names are now accepted everywhere and resolve
# through this helper.
#
# If both are supplied and they differ, `facet` wins and a message is emitted.
# Returns NULL when neither is set (no panelling).
# ============================================================================
.ctdna_resolve_facet <- function(facet = NULL, wrap = NULL, fn = NULL) {
  a <- facet; b <- wrap
  if (is.null(a) || !length(a)) return(b)
  if (is.null(b) || !length(b)) return(a)
  if (!setequal(as.character(a), as.character(b)))
    message(sprintf(
      "%s`facet` and `wrap` were both supplied and differ; using `facet`. ",
      if (is.null(fn)) "" else paste0(fn, ": ")),
      "(They are aliases for the same panelling argument.)")
  a
}


# ============================================================================
# .ctdna_pooled_note(prep, dims, fn = NULL)
# ----------------------------------------------------------------------------
# v0.80.0: when a prep spans several indications but indication is not one of
# the plot's display dimensions, every cohort is silently pooled into the same
# boxes. That is a legitimate analysis, but it is easy to mistake a pooled plot
# for a single-cohort one -- so the fact is stated on the figure itself, in the
# caption, where it travels into reports.
#
# `dims` is the character vector of columns the caller is already displaying
# (group_by / subgroup_by / facet / x / colour). Returns a one-line string, or
# NULL when there is nothing to say: a single-indication prep, indication
# already on display, or the notice switched off via
# ctdna_opts(warn_pooled_dims = FALSE).
# ============================================================================
.ctdna_pooled_note <- function(prep, dims = character(0), fn = NULL) {
  on <- tryCatch(.o("warn_pooled_dims"), error = function(e) TRUE)
  if (!isTRUE(on)) return(NULL)

  df <- if (is.data.frame(prep)) prep else prep$samples
  if (!is.data.frame(df)) return(NULL)
  ind_col <- .ctdna_find_ind_col(names(df))
  if (is.na(ind_col)) return(NULL)

  lv <- unique(as.character(df[[ind_col]]))
  lv <- sort(lv[!is.na(lv) & nzchar(lv)])
  if (length(lv) < 2L) return(NULL)

  # already a display dimension? compare normalised, so "indication" matches
  # "Indication" and the opts-resolved name alike.
  dims <- dims[!vapply(dims, is.null, logical(1))]
  dims <- unlist(dims, use.names = FALSE)
  if (length(dims) && any(.ctdna_norm_ind(dims) == .ctdna_norm_ind(ind_col)))
    return(NULL)

  sprintf("Note: %d indications pooled (%s)", length(lv),
          paste(lv, collapse = ", "))
}
