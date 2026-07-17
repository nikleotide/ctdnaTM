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
# .ctdna_filter_by_indication(df, indications, prep = NULL, col_hint = NULL)
# ----------------------------------------------------------------------------
# Shared cohort-restriction helper used by every user-facing plot function.
# Filters `df` (any prep-slot frame) to only patients whose indication is in
# `indications`.
#
# Lookup order for the indication column:
#   1. `col_hint` if provided and it exists in `df`.
#   2. "Indication" / "indication" / "Cancertype" / "Cancer_Type" in `df`.
#   3. If none of the above are in df but `prep$clinical` is available,
#      look them up there and translate Patient_ID -> Indication.
#   4. Otherwise -> stop() with an actionable message.
#
# `indications = NULL`  ->  no-op, returns df unchanged.
# ============================================================================
.ctdna_filter_by_indication <- function(df, indications,
                                          prep = NULL, col_hint = NULL) {
  if (is.null(indications) || !length(indications)) return(df)
  if (!is.data.frame(df))                            return(df)

  cands <- c(col_hint, "Indication", "indication", "Cancertype", "Cancer_Type")
  cands <- cands[!is.na(cands) & nzchar(cands)]

  # Fast path: indication column is already in df.
  in_df <- cands[cands %in% names(df)][1]
  if (!is.na(in_df)) {
    n_before  <- nrow(df)
    df        <- df[as.character(df[[in_df]]) %in% as.character(indications), ,
                     drop = FALSE]
    n_dropped <- n_before - nrow(df)
    if (n_dropped > 0)
      message(sprintf("indications: dropped %d row(s) outside {%s}.",
                       n_dropped, paste(indications, collapse = ", ")))
    return(df)
  }

  # Fallback: look up in prep$clinical, then filter by Patient_ID.
  if (!is.null(prep) && is.data.frame(prep$clinical) &&
      "Patient_ID" %in% names(df) &&
      "Patient_ID" %in% names(prep$clinical)) {
    in_cl <- cands[cands %in% names(prep$clinical)][1]
    if (!is.na(in_cl)) {
      keep_pids <- as.character(prep$clinical$Patient_ID)[
        as.character(prep$clinical[[in_cl]]) %in% as.character(indications)]
      n_before  <- nrow(df)
      df        <- df[as.character(df$Patient_ID) %in% keep_pids, ,
                       drop = FALSE]
      n_dropped <- n_before - nrow(df)
      if (n_dropped > 0)
        message(sprintf("indications: dropped %d row(s) outside {%s}.",
                         n_dropped, paste(indications, collapse = ", ")))
      return(df)
    }
  }

  stop("indications: no indication column found in df (looked for: ",
        paste(cands, collapse = ", "),
        "). Pass a prep object or add an Indication column to the frame.",
        call. = FALSE)
}
