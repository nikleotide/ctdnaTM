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
