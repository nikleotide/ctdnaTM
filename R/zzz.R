# ============================================================================
# ctdnaTM — package onAttach hook
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
# .ctdna_filter_prep_by_indication(prep, indications, col_hint = NULL)
# ----------------------------------------------------------------------------
# Shared cohort-restriction helper for plot functions whose first argument
# is a `ctdna_prep` object.
#
# `indications = NULL` -> return prep unchanged.
# When set, keep only patients whose Indication (or fallback column) is in
# the given values, and cascade the filter to every prep slot that has a
# Patient_ID column.
#
# Indication column lookup order in prep$clinical: `col_hint`, "Indication",
# "indication", "Cancertype", "Cancer_Type". First one that exists wins.
# ============================================================================
.ctdna_filter_prep_by_indication <- function(prep, indications,
                                              col_hint = NULL) {
  if (is.null(indications) || !length(indications)) return(prep)
  if (!is.list(prep) || !is.data.frame(prep$clinical)) return(prep)

  cands <- c(col_hint, "Indication", "indication", "Cancertype", "Cancer_Type")
  cands <- cands[!is.na(cands) & nzchar(cands)]
  ind_col <- cands[cands %in% names(prep$clinical)][1]
  if (is.na(ind_col))
    stop("indications: no indication column found in prep$clinical (looked ",
          "for: ", paste(cands, collapse = ", "),
          "). Add one before passing `indications`.", call. = FALSE)

  keep_pids <- as.character(prep$clinical$Patient_ID)[
    as.character(prep$clinical[[ind_col]]) %in% as.character(indications)]

  if (!length(keep_pids)) {
    warning(sprintf(
      "indications: no patients match {%s} in %s. Returning empty prep.",
      paste(indications, collapse = ", "), ind_col), call. = FALSE)
  }

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
    message(sprintf("indications: kept %d patient(s) in {%s}; dropped %d row(s) across prep slots.",
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
# Looks for an indication column in `df` itself. If not present, silently
# returns df unchanged (with a warning) since we have no prep to join against.
# ============================================================================
.ctdna_filter_df_by_indication <- function(df, indications,
                                             col_hint = NULL) {
  if (is.null(indications) || !length(indications)) return(df)
  if (!is.data.frame(df)) return(df)

  cands <- c(col_hint, "Indication", "indication", "Cancertype", "Cancer_Type")
  cands <- cands[!is.na(cands) & nzchar(cands)]
  in_df <- cands[cands %in% names(df)][1]
  if (is.na(in_df)) {
    warning("indications: no indication column found in the data.frame ",
             "(looked for: ", paste(cands, collapse = ", "),
             "). Filter is a no-op. Add an Indication column to the frame ",
             "before passing `indications`.", call. = FALSE)
    return(df)
  }

  n_before <- nrow(df)
  df <- df[as.character(df[[in_df]]) %in% as.character(indications), ,
            drop = FALSE]
  n_dropped <- n_before - nrow(df)
  if (n_dropped > 0)
    message(sprintf("indications: dropped %d row(s) outside {%s}.",
                     n_dropped, paste(indications, collapse = ", ")))
  df
}
