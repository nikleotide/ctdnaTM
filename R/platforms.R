# ============================================================================
# ctdnaTM — platform preparation
# ============================================================================
# ctdna_prepare() is the entry point that takes RAW vendor data (the
# multiple files delivered per ctDNA platform) and produces the
# CANONICAL analysis-ready data structures (column names per
# ctdna_opts()) that the rest of the package consumes.
#
# For Guardant Health Infinity, raw inputs are:
#   - infinity_report      : per-sample per-alteration master report
#   - tf_change            : paired Visit_A vs Visit_B TF change (TF.methyl)
#   - panel_74_response    : paired 74-gene MR response
#   - panel_500_response   : paired 500-gene MR response
#   - methylation_response : paired methylation MR response (legacy.methyl)
#   - clinical             : per-patient efficacy / dose / RECIST / MR
#
# Canonical outputs (column names follow ctdna_opts()):
#   $ctdna       : long-format longitudinal TF
#                  (subject, time, dose, RECIST, methylTF, maxVAF,
#                   best_pct_change, MR)
#   $genomic_74  : baseline per-(subject × gene) calls
#                  (subject, gene, panel, alteration_type,
#                   mutation_status, vaf)
#   $genomic_500 : same shape, larger gene set
#
# Each platform registers a handler in `.platforms`. To add a new
# platform, append a list entry with $name and $prepare.
# ============================================================================


# ---- 74- and 500-gene panels ------------------------------------------------
# `.gh_genes_74` and `.gh_genes_500_only` are file-level constants
# defined in ctdnaTM.R (alongside ctdna_make_mock_study). Both files are
# loaded at package load so the constants are visible here.


# ---- Vendor variant_type → canonical alteration_type ------------------------
# Canonical values: {wt, snv, cnv_amp, cnv_loss, lgr, fusion}.
.gh_alt_type <- function(variant_type, copy_number) {
  out <- character(length(variant_type))
  out[variant_type == "SNV"]    <- "snv"
  out[variant_type == "Indel"]  <- "snv"           # collapse indel → snv class
  out[variant_type == "Fusion"] <- "fusion"
  out[variant_type == "LGR"]    <- "lgr"
  cnv_mask <- variant_type == "CNV"
  if (any(cnv_mask)) {
    cn <- copy_number[cnv_mask]
    out[cnv_mask] <- ifelse(!is.na(cn) & cn >= 4, "cnv_amp", "cnv_loss")
  }
  out
}


# ---- Build canonical ctDNA longitudinal frame from vendor data --------------
# Unstacks the paired Visit_A vs Visit_B representation into long format,
# joins per-sample maxVAF (from infinity_report) and per-patient clinical
# fields (RECIST, dose, best_pct_change, MR).
.gh_build_ctdna_long <- function(tf_change,
                                  infinity_report = NULL,
                                  clinical = NULL,
                                  qc_filter = TRUE,
                                  verbose = TRUE) {

  if (qc_filter) {
    n0 <- nrow(tf_change)
    tf_change <- tf_change[
      tf_change$MR_QC_status   %in% c("PASS","SUCCESS","Pass","pass") &
      tf_change$MR_quantifiable %in% c("Yes","TRUE","yes"), ]
    if (verbose)
      message(sprintf("  tf_change: kept %d / %d pairs after MR_QC filter",
                      nrow(tf_change), n0))
  }

  # Long-format: one row per (Patient × Visit). Both legs of the pair
  # contribute (A always baseline; B always later). Dedupe at the end.
  #
  # v0.41.6: the A-leg is, by construction, the baseline measurement.
  # Tag its `time_point` as `ctdna_opts("baseline")` (default
  # `"Baseline"`) so downstream code that filters
  # `df[time_point == ctdna_opts("baseline"), ]` (e.g.
  # `ctdna_metric_at()`, ratio computation in
  # `ctdna_plot_reduction()`) actually finds those rows. The original
  # cycle label of the A-leg is preserved in a sidecar column
  # `time_point_raw` for any code or downstream analysis that needs
  # the actual visit identifier (e.g., screening vs C1D1 baselines).
  bl_label <- .o("baseline")
  a_rows <- data.frame(
    subject_id      = tf_change$Patient_ID,
    time_point      = bl_label,                          # canonical "Baseline"
    time_point_raw  = tf_change$Visit_name_A,            # original cycle label
    methylTF        = tf_change$Methylation_tumor_fraction_percentage_A / 100,
    .date           = tf_change$Bloodcoll_date_A,
    stringsAsFactors = FALSE
  )
  b_rows <- data.frame(
    subject_id      = tf_change$Patient_ID,
    time_point      = tf_change$Visit_name_B,            # cycle label
    time_point_raw  = tf_change$Visit_name_B,
    methylTF        = tf_change$Methylation_tumor_fraction_percentage_B / 100,
    .date           = tf_change$Bloodcoll_date_B,
    stringsAsFactors = FALSE
  )
  long <- rbind(a_rows, b_rows)
  long <- long[!duplicated(long[, c("subject_id","time_point")]), ]
  long <- long[order(long$subject_id, long$.date), ]

  # ---- maxVAF + meanVAF per sample (from infinity_report) ----
  # v0.41.3: meanVAF added. The method registry has always offered
  # method = "meanvaf", but the platform handler never produced the
  # column, so method = "all" silently lost a panel and method =
  # "meanvaf" alone errored. Both are now computed in one aggregation
  # pass to keep the I/O cost flat.
  if (!is.null(infinity_report)) {
    ir <- infinity_report[infinity_report$Sample_status %in% c("PASS","SUCCESS","Pass","pass"), ]
    ir <- ir[!is.na(ir$VAF_percentage), ]
    if (nrow(ir) > 0) {
      vaf_agg <- stats::aggregate(VAF_percentage ~ Patient_ID + Visit_name,
                                   data = ir,
                                   FUN = function(v) c(max = max(v),
                                                       mean = mean(v)))
      # aggregate returns a matrix in the value column; split it
      vaf_df <- data.frame(
        subject_id     = vaf_agg$Patient_ID,
        time_point_raw = vaf_agg$Visit_name,            # raw cycle label
        maxVAF         = as.numeric(vaf_agg$VAF_percentage[, "max"])  / 100,
        meanVAF        = as.numeric(vaf_agg$VAF_percentage[, "mean"]) / 100,
        stringsAsFactors = FALSE)
      # v0.41.6: merge on `time_point_raw` (raw cycle label) so the
      # baseline row in `long` (whose `time_point` is now "Baseline" but
      # whose `time_point_raw` still holds the actual visit label e.g.
      # "C1D1") still picks up its VAF entry.
      long <- merge(long, vaf_df, by = c("subject_id","time_point_raw"),
                     all.x = TRUE)
    } else {
      long$maxVAF  <- NA_real_
      long$meanVAF <- NA_real_
    }
  } else {
    long$maxVAF  <- NA_real_
    long$meanVAF <- NA_real_
  }

  # ---- Clinical join (RECIST, dose, best_pct_change, MR) ----
  if (!is.null(clinical)) {
    cli <- data.frame(
      subject_id       = clinical$Patient_ID,
      dose             = if ("Dose" %in% names(clinical)) clinical$Dose
                         else if ("dose" %in% names(clinical)) clinical$dose
                         else NA_character_,
      RECIST           = clinical$RECIST,
      best_pct_change  = clinical$best_pct_change,
      MR               = clinical$MR,
      stringsAsFactors = FALSE
    )
    long <- merge(long, cli, by = "subject_id", all.x = TRUE)
    long <- long[order(long$subject_id, long$.date), ]
  } else {
    long$dose            <- NA_character_
    long$RECIST          <- NA_character_
    long$best_pct_change <- NA_real_
    long$MR              <- NA_character_
  }

  # Drop internal .date helper before returning
  long$.date <- NULL
  rownames(long) <- NULL

  # Reorder to a sensible canonical column order
  col_order <- c("subject_id","time_point","dose","RECIST","methylTF",
                  "maxVAF","meanVAF","best_pct_change","MR")
  col_order <- intersect(col_order, names(long))
  long <- long[, col_order, drop = FALSE]

  # Apply current ctdna_opts() renaming
  .apply_opts_names(long, c(
    subject_id      = "subject",
    time_point      = "time",
    dose            = "dose",
    RECIST          = "recist",
    methylTF        = "tf",
    maxVAF          = "maxvaf",
    meanVAF         = "meanvaf",
    best_pct_change = "best_change",
    MR              = "mr"
  ))
}


# ---- v0.43.3: per-variant column blacklist ----
# Columns that vary per-variant within a single sample (one row per variant).
# Used to split the consolidated InfinityReport into "per-sample" and
# "per-variant" slices.
.gh_variant_only_cols <- c(
  "Variant_type","Indel_type","Gene","Chromosome","Position","Exon",
  "Mut_aa","Mut_nt","Mut_cdna","Transcript","VAF_percentage",
  "Splice_effect","Somatic_status","Molecular_consequence",
  "Fusion_chrom_b","Fusion_gene_b","Fusion_position_a","Fusion_position_b",
  "Direction_a","Direction_b","Downstream_gene",
  "Copy_number","CNV_type",
  "COSMIC","dbSNP","ClinVar","ClinVarID",
  "Functional_impact","Mutant_allele_status","Mol_count","Alleletype",
  # canonicalized versions of the same things (in case alias map renamed them)
  "Alteration"
)

# v0.43.3: build $samples — every per-sample column from infinity_report
# (i.e. everything that ISN'T in the per-variant blacklist), deduplicated
# to one row per unique sample. Used by ctdna_sample_qc() and accessible
# directly for any per-sample biomarker (HRD_score, MSI_High, etc.).
.gh_samples_frame <- function(infinity_report) {
  if (!is.data.frame(infinity_report)) return(NULL)
  per_sample <- setdiff(names(infinity_report), .gh_variant_only_cols)
  if (length(per_sample) == 0L) return(NULL)
  sf <- unique(infinity_report[, per_sample, drop = FALSE])
  rownames(sf) <- NULL
  sf
}

# v0.43.3: build $variants — every per-variant row from the infinity_report
# (rows where Gene is non-empty), preserving ALL columns. The raw per-variant
# view; $genomic_74 / $genomic_500 are derived analysis views.
.gh_variants_frame <- function(infinity_report) {
  if (!is.data.frame(infinity_report)) return(NULL)
  if (!"Gene" %in% names(infinity_report)) return(NULL)
  g <- infinity_report$Gene
  mask <- !is.na(g) & nzchar(as.character(g))
  if (!any(mask)) return(NULL)
  vf <- infinity_report[mask, , drop = FALSE]
  rownames(vf) <- NULL
  vf
}


# ---- v0.43.2: build $ctdna directly from a consolidated InfinityReport ----
# Used when tf_change / methylation_response are not provided separately.
# Real-world consolidated InfinityReport CSVs include per-sample longitudinal
# columns (methylTF, maxVAF, TMB_score, ctDNA_detection_status) — dedupe to
# one row per (Patient_ID, Visit_name) and join clinical if available.
.gh_ctdna_from_infinity_report <- function(infinity_report, clinical = NULL,
                                            verbose = TRUE) {
  if (!is.data.frame(infinity_report)) return(NULL)
  if (!all(c("Patient_ID","Visit_name") %in% names(infinity_report))) {
    if (verbose)
      message("  cannot build $ctdna from infinity_report: missing Patient_ID or Visit_name")
    return(NULL)
  }

  # v0.43.3: use the per-variant blacklist to keep EVERY per-sample column
  # automatically. This avoids dropping HRD_score, MSI_High, TMB_category,
  # Cancertype, etc., that earlier versions silently lost.
  per_sample_cols <- setdiff(names(infinity_report), .gh_variant_only_cols)

  ct <- unique(infinity_report[, per_sample_cols, drop = FALSE])
  rownames(ct) <- NULL

  # Coerce known-numeric columns for downstream analyses
  for (nm in intersect(c("methylTF","maxVAF","meanVAF","TMB_score",
                          "HRD_score","cfDNA_ng","Plasma_ml_input",
                          "Plasma_ml_remaining"),
                        names(ct))) {
    ct[[nm]] <- suppressWarnings(as.numeric(ct[[nm]]))
  }

  # v0.43.4: derive the per-sample analysis columns downstream plot
  # functions expect. The raw InfinityReport doesn't contain these
  # directly; they have to be computed from the per-variant rows.
  #
  #   meanVAF    = mean(VAF_percentage) across that sample's variants
  #   n.somatic  = count of variants per sample (Gene non-empty)
  #   cfdna.conc = cfDNA_ng / Plasma_ml_input  (ng/mL plasma)
  ir <- infinity_report
  if (!"meanVAF" %in% names(ct) && "VAF_percentage" %in% names(ir)) {
    vaf_num <- suppressWarnings(as.numeric(ir$VAF_percentage))
    agg <- stats::aggregate(
      vaf_num,
      by = list(Patient_ID = ir$Patient_ID, Visit_name = ir$Visit_name),
      FUN = function(x) mean(x, na.rm = TRUE))
    names(agg)[3] <- "meanVAF"
    agg$meanVAF[is.nan(agg$meanVAF)] <- NA_real_
    ct <- merge(ct, agg, by = c("Patient_ID","Visit_name"),
                all.x = TRUE, sort = FALSE)
  }
  if (!"n.somatic" %in% names(ct) && "Gene" %in% names(ir)) {
    is_variant <- !is.na(ir$Gene) & nzchar(as.character(ir$Gene))
    agg <- stats::aggregate(
      is_variant,
      by = list(Patient_ID = ir$Patient_ID, Visit_name = ir$Visit_name),
      FUN = sum)
    names(agg)[3] <- "n.somatic"
    ct <- merge(ct, agg, by = c("Patient_ID","Visit_name"),
                all.x = TRUE, sort = FALSE)
  }
  if (!"cfdna.conc" %in% names(ct) &&
      all(c("cfDNA_ng","Plasma_ml_input") %in% names(ct))) {
    pm <- suppressWarnings(as.numeric(ct$Plasma_ml_input))
    ct$cfdna.conc <- ifelse(is.na(pm) | pm <= 0, NA_real_,
                              suppressWarnings(as.numeric(ct$cfDNA_ng)) / pm)
  }

  # v0.43.8: deterministic percentage→fraction conversion based on the
  # SOURCE COLUMN NAME, not on value ranges. If the canonical column was
  # canonicalized from a source name containing "percent" (e.g.
  # `Methylation_tumor_fraction_percentage` -> `methylTF`,
  # `Genomic_max_VAF_percentage` -> `maxVAF`), divide by 100 to convert
  # to fraction-scale. The library's plot functions apply percent_format()
  # which assumes fractions, so this is mandatory for correct display.
  #
  # v0.43.6's value-heuristic (max > 1 → convert) failed when ALL values
  # were < 1 but still represented percentages (e.g. methylTF = 0.5
  # meaning 0.5%, not 50%). Source-name detection has no such ambiguity.
  renames <- attr(infinity_report, "ctdna_renames")
  for (col in c("methylTF","maxVAF")) {
    if (col %in% names(ct) && !is.null(renames) && col %in% names(renames)) {
      src <- renames[[col]]
      if (grepl("percent", src, ignore.case = TRUE)) {
        ct[[col]] <- ct[[col]] / 100
        if (verbose)
          message(sprintf("  converted `%s` from percentage to fraction (source: `%s`)",
                          col, src))
      }
    }
  }
  # meanVAF is always derived from VAF_percentage above, so always /100
  if ("meanVAF" %in% names(ct) &&
      "VAF_percentage" %in% names(infinity_report)) {
    ct$meanVAF <- ct$meanVAF / 100
    if (verbose)
      message("  converted derived `meanVAF` from percentage to fraction (source: VAF_percentage)")
  }

  # Join clinical (if present and shares Patient_ID)
  if (!is.null(clinical) && is.data.frame(clinical) &&
      "Patient_ID" %in% names(clinical)) {
    clin_carry <- intersect(c(
      "Patient_ID","Dose","RECIST","Sex","Age","Indication",
      "ARM","ARMCD","TRT01A","TRT01P","ACTARM","ACTARMCD","COHORT",
      "indication_dose"
    ), names(clinical))
    if (length(clin_carry) >= 2L) {
      clin_pt <- unique(clinical[, clin_carry, drop = FALSE])
      clin_pt <- clin_pt[!duplicated(clin_pt$Patient_ID), , drop = FALSE]
      ct <- merge(ct, clin_pt, by = "Patient_ID",
                  all.x = TRUE, sort = FALSE)
    }
  }

  if (verbose)
    message(sprintf("  built $ctdna from infinity_report: %d samples × %d cols",
                    nrow(ct), ncol(ct)))
  ct
}


# ---- Build canonical baseline-only genomic frame ----------------------------
# One row per (subject × gene) for the panel's gene list. Detected
# alterations get mutation_status = "mut" and the appropriate
# alteration_type; non-detected get mutation_status = "wt".
.gh_build_genomic <- function(infinity_report,
                               panel_genes,
                               panel_label,
                               baseline_visit = NULL,
                               verbose = TRUE) {

  # v0.43.0: pure conversion. No Sample_status filter, no Visit_name
  # filter by default. Filtering moved to ctdna_sample_qc() so users
  # explicitly opt in. The legacy `baseline_visit` argument is kept for
  # backward compat: when set, restricts to those visits (case-insensitive,
  # exact match -- no fuzzy fallback). When NULL (default), uses every row.

  if (!is.null(baseline_visit) && length(baseline_visit) > 0L) {
    vn <- if ("Visit_name" %in% names(infinity_report))
            trimws(as.character(infinity_report$Visit_name))
          else character(nrow(infinity_report))
    visit_mask <- tolower(vn) %in% tolower(as.character(baseline_visit))
    ir <- infinity_report[visit_mask, , drop = FALSE]
    if (nrow(ir) == 0L && nrow(infinity_report) > 0L) {
      warning(sprintf(
        paste0(".gh_build_genomic[%s]: baseline_visit = c(%s) matched no rows. ",
               "Unique Visit_name values present: %s. Falling back to ALL rows."),
        panel_label,
        paste(shQuote(as.character(baseline_visit)), collapse=", "),
        paste(unique(vn), collapse=", ")), call. = FALSE)
      ir <- infinity_report
    }
  } else {
    ir <- infinity_report
  }

  all_patients <- sort(unique(ir$Patient_ID))
  if (length(all_patients) == 0)
    return(data.frame(subject_id=character(), gene=character(),
                       panel=character(), alteration_type=character(),
                       mutation_status=character(), vaf=numeric(),
                       stringsAsFactors = FALSE))

  detected <- ir[!is.na(ir$Gene) & ir$Gene %in% panel_genes, ]
  if (nrow(detected) > 0) {
    det_df <- data.frame(
      subject_id       = detected$Patient_ID,
      gene             = detected$Gene,
      panel            = panel_label,
      alteration_type  = .gh_alt_type(detected$Variant_type,
                                       detected$Copy_number),
      mutation_status  = "mut",
      vaf              = detected$VAF_percentage,
      stringsAsFactors = FALSE
    )

    # v0.43.3: carry through all per-variant annotation columns when
    # present. WT rows will get NA for these. Downstream code that uses
    # the minimal schema (subject_id/gene/alteration_type/mutation_status/vaf)
    # still works unchanged.
    extra_cols <- intersect(c(
      "Visit_name","Sample_ID","Sample_status",
      "Variant_type","Indel_type","Chromosome","Position","Exon",
      "Mut_aa","Mut_nt","Mut_cdna","Transcript",
      "Splice_effect","Somatic_status","Molecular_consequence",
      "Fusion_chrom_b","Fusion_gene_b","Fusion_position_a","Fusion_position_b",
      "Direction_a","Direction_b","Downstream_gene",
      "Copy_number","CNV_type",
      "COSMIC","dbSNP","ClinVar","ClinVarID",
      "Functional_impact","Mutant_allele_status","Mol_count","Alleletype"
    ), names(detected))
    for (col in extra_cols) {
      det_df[[col]] <- detected[[col]]
    }

    # When the same (subject, gene) has multiple variants, collapse to one
    # row keeping the highest-impact (snv > cnv_amp > fusion > others).
    impact_rank <- c(snv = 1, cnv_amp = 2, cnv_loss = 3, fusion = 4,
                      lgr = 5, wt = 6)
    det_df$.rank <- impact_rank[det_df$alteration_type]
    det_df <- det_df[order(det_df$subject_id, det_df$gene, det_df$.rank), ]
    det_df <- det_df[!duplicated(det_df[, c("subject_id","gene")]), ]
    det_df$.rank <- NULL
  } else {
    # Empty det_df with the FULL minimal schema (no extra cols)
    det_df <- data.frame(subject_id=character(), gene=character(),
                          panel=character(), alteration_type=character(),
                          mutation_status=character(), vaf=numeric(),
                          stringsAsFactors = FALSE)
  }

  # Enumerate wt rows for (Patient × Gene) combos that weren't detected.
  # WT rows get NA for any extension column the det_df has.
  grid <- expand.grid(subject_id = all_patients, gene = panel_genes,
                       stringsAsFactors = FALSE)
  det_keys <- paste(det_df$subject_id, det_df$gene)
  grid_keys <- paste(grid$subject_id, grid$gene)
  wt <- grid[!(grid_keys %in% det_keys), , drop = FALSE]
  wt$panel           <- panel_label
  wt$alteration_type <- "wt"
  wt$mutation_status <- "wt"
  wt$vaf             <- NA_real_
  # Fill in NA for any extra columns det_df has (so rbind works)
  extra_in_det <- setdiff(names(det_df), names(wt))
  for (col in extra_in_det) {
    template <- det_df[[col]]
    wt[[col]] <- if (is.numeric(template)) NA_real_
                  else if (is.integer(template)) NA_integer_
                  else NA_character_
  }
  wt <- wt[, names(det_df), drop = FALSE]

  out <- rbind(det_df, wt)
  out <- out[order(out$subject_id, out$gene), ]
  rownames(out) <- NULL
  if (verbose)
    message(sprintf("  %s: %d detected + %d wt = %d rows × %d cols",
                    panel_label, nrow(det_df), nrow(wt),
                    nrow(out), ncol(out)))

  .apply_opts_names(out, c(
    subject_id       = "subject",
    gene             = "gene",
    panel            = "panel",
    alteration_type  = "alteration",
    mutation_status  = "mutation"
  ))
}


# ---- Guardant Health Infinity prep ------------------------------------------
.prepare_gh_infinity <- function(infinity_report      = NULL,
                                  tf_change            = NULL,
                                  panel_74_response    = NULL,
                                  panel_500_response   = NULL,
                                  methylation_response = NULL,
                                  clinical             = NULL,
                                  qc_filter            = FALSE,
                                  baseline_visit       = NULL,
                                  verbose              = TRUE) {
  out <- list()

  # ctdna longitudinal frame: prefer tf_change (TF.methyl), fall back to
  # methylation_response (legacy.methyl). Warn if only legacy is present.
  src <- tf_change
  src_name <- "tf_change (TF.methyl)"
  if (is.null(src) && !is.null(methylation_response)) {
    warning("Only legacy.methyl is present (no TF.methyl tf_change). ",
            "Guardant Health is phasing out legacy.methyl in favor of ",
            "TF.methyl; treat results accordingly.")
    # Synthesize TF_A / TF_B from MR_score_percentage (legacy.methyl format
    # doesn't carry the explicit A/B TF columns).
    # Each row already has Visit_A baseline / Visit_B later; the MR score
    # is the Visit_B TF, so reconstruct A by joining baselines.
    mr <- methylation_response
    mr$Methylation_tumor_fraction_percentage_B <- mr$MR_score_percentage
    # For A, use the baseline pair (Visit_name_A=="C1D1") MR score from
    # the same patient — pull from the first row per patient.
    base <- mr[!duplicated(mr$Patient_ID), ]
    base_lookup <- stats::setNames(base$MR_score_percentage / 1, base$Patient_ID)
    # We don't actually have separate A TF; approximate as patient's max
    # MR_score across that patient's pairs (the most baseline-like value).
    pa <- stats::aggregate(MR_score_percentage ~ Patient_ID, data = mr, FUN = max)
    base_lookup <- stats::setNames(pa$MR_score_percentage, pa$Patient_ID)
    mr$Methylation_tumor_fraction_percentage_A <-
      base_lookup[as.character(mr$Patient_ID)]
    mr$Bloodcoll_date_A <- mr$Bloodcoll_date_B <- NA
    mr$MR_QC_status     <- mr$MR_QC_status
    mr$MR_quantifiable  <- mr$MR_quantifiable
    src <- mr
    src_name <- "methylation_response (legacy.methyl)"
  }

  if (!is.null(src)) {
    if (verbose) message("  building ctdna long-format from ", src_name)
    out$ctdna <- .gh_build_ctdna_long(
      tf_change       = src,
      infinity_report = infinity_report,
      clinical        = clinical,
      qc_filter       = qc_filter,
      verbose         = verbose
    )
  }

  # v0.43.2: fallback — if no tf_change / methylation_response was provided
  # but the infinity_report has the per-sample longitudinal columns
  # (methylTF, maxVAF, TMB_score, ctDNA_detection_status, etc.) we can
  # build $ctdna directly from it.
  if (is.null(out$ctdna) && !is.null(infinity_report)) {
    out$ctdna <- .gh_ctdna_from_infinity_report(infinity_report, clinical,
                                                  verbose = verbose)
  }

  # genomic_74: requires infinity_report (per-variant rows).
  if (!is.null(infinity_report)) {
    if (verbose) message("  building genomic_74 from infinity_report")
    out$genomic_74 <- .gh_build_genomic(
      infinity_report = infinity_report,
      panel_genes     = .gh_genes_74,
      panel_label     = "74genes",
      baseline_visit  = baseline_visit,
      verbose         = verbose
    )
    # v0.43.2: always also build genomic_500 from the infinity_report,
    # regardless of whether the separate panel_500_response / panel_74_response
    # summary files were provided. Real consolidated InfinityReport CSVs
    # contain variant rows for the full 500-gene panel.
    if (verbose) message("  building genomic_500 from infinity_report")
    out$genomic_500 <- .gh_build_genomic(
      infinity_report = infinity_report,
      panel_genes     = c(.gh_genes_74, .gh_genes_500_only),
      panel_label     = "500genes",
      baseline_visit  = baseline_visit,
      verbose         = verbose
    )
  }

  # v0.43.2: attach clinical as a standalone output frame (previously it
  # was only consumed by .gh_build_ctdna_long for joining; the standalone
  # frame is needed by downstream functions).
  if (!is.null(clinical) && is.data.frame(clinical)) {
    out$clinical <- clinical
  }

  # v0.43.3: per-sample frame — every per-sample column from the
  # infinity_report (not just QC). Used by ctdna_sample_qc() and any
  # downstream code that needs HRD_score, MSI_High, TMB_category,
  # Cancertype, etc.
  if (!is.null(infinity_report)) {
    out$samples <- .gh_samples_frame(infinity_report)
  }

  # v0.43.3: per-variant frame — every variant row from the
  # infinity_report, preserving all annotation columns (Chromosome,
  # Position, Mut_aa, Mut_nt, COSMIC, ClinVar, Functional_impact, fusion
  # fields, etc.). $genomic_74 / $genomic_500 are derived analysis views;
  # $variants is the raw access point.
  if (!is.null(infinity_report)) {
    out$variants <- .gh_variants_frame(infinity_report)
  }

  out
}


# ---- Platform registry ------------------------------------------------------
.platforms <- list(
  guardant_health_infinity = list(
    name    = "Guardant Health Infinity",
    prepare = .prepare_gh_infinity
  )
  # future platforms here, e.g.:
  # natera_signatera = list(name = "Natera Signatera", prepare = .prepare_natera_signatera),
)


# ---- Public API -------------------------------------------------------------

#' Prepare a ctDNA study bundle for downstream analysis
#'
#' Takes a study bundle in raw vendor format (e.g. the output of
#' \code{\link{ctdna_make_mock_study}} or a hand-built list whose names
#' match those keys) and returns the analysis-ready canonical frames.
#' The function auto-detects which modalities are present in
#' \code{sim} and processes only those — missing modalities are
#' skipped with an informational message; nothing errors.
#'
#' \strong{Modalities handled (v0.40.0)}
#' \itemize{
#'   \item \strong{ctDNA bundle} (any of \code{infinity_report},
#'     \code{tf_change}, \code{panel_74_response},
#'     \code{panel_500_response}, \code{methylation_response}) →
#'     dispatched to the platform handler. Produces \code{ctdna},
#'     \code{genomic_74}, \code{genomic_500} in canonical column
#'     names per \code{\link{ctdna_opts}}.
#'   \item \strong{ihc} → cleaned + canonicalized.
#'   \item \strong{dnaseq} (tissue WES/WGS alteration calls) →
#'     cleaned + canonicalized. Also accepts the legacy key
#'     \code{wes}.
#'   \item \strong{rnaseq} (list with \code{counts}, \code{tpm},
#'     \code{sample_metadata}) → reshaped to canonical long-form
#'     \code{expression} (one row per subject × gene, in
#'     \code{log2(TPM + 1)}).
#'   \item \strong{tmb} → derived from \code{infinity_report}'s
#'     \code{TMB_score} column (per-sample). Skipped if no
#'     \code{infinity_report} is supplied.
#'   \item \strong{clinical} → if present in \code{sim}, joined into
#'     the ctdna frame by the platform handler.
#' }
#'
#' \strong{Column naming.} Output frames use the canonical column
#' names defined by \code{ctdna_opts()} (e.g. \code{subject_id} for
#' the patient identifier, \code{gene}, \code{expression}). Pass
#' \code{new_col_names} to override any of these for downstream
#' compatibility (e.g. when the rest of your pipeline expects
#' \code{Patient_ID}).
#'
#' @param sim Named list of vendor-format frames. Recognized keys:
#'   \code{infinity_report}, \code{tf_change},
#'   \code{panel_74_response}, \code{panel_500_response},
#'   \code{methylation_response}, \code{clinical}, \code{ihc},
#'   \code{dnaseq} (or \code{wes}), \code{rnaseq}. Any subset
#'   accepted; nothing required.
#' @param vendor Vendor identifier. Default \code{"guardant_health"}.
#' @param technology Technology identifier. Default \code{"infinity"}.
#' @param qc_filter If \code{TRUE} (default), drop rows that fail
#'   vendor QC during the platform transformation.
#' @param cleanup If \code{TRUE} (default), run light janitorial steps
#'   (trim whitespace in ID columns, coerce factor IDs to character,
#'   drop empty rows) on each output frame.
#' @param filter_scheme Filter scheme name passed to
#'   \code{\link{ctdna_variant_filter}} for alteration-level filtering
#'   on the produced genomic frames. Default \code{"default"}; pass
#'   \code{NULL} to skip.
#' @param new_col_names Optional NAMED list mapping
#'   \code{old_column_name = "new_column_name"}. After all other
#'   processing, every output data frame is scanned and any column
#'   whose name matches a key in this list is renamed to the
#'   corresponding value. Use it to bridge to a downstream pipeline
#'   that expects different names, e.g.
#'   \code{new_col_names = list(subject_id = "Patient_ID",
#'   gene = "Gene_symbol")}. Default \code{NULL}.
#' @param verbose If \code{TRUE} (default), print informational
#'   messages about which modalities were prepared or skipped.
#' @return A named list. The \code{platform} attribute records the
#'   resolved platform key. Possible elements (each present only when
#'   the corresponding input was supplied in \code{sim}):
#'   \itemize{
#'     \item \code{ctdna} — long-format ctDNA frame
#'     \item \code{genomic_74}, \code{genomic_500} — panel alteration
#'       frames
#'     \item \code{ihc} — long-format IHC frame
#'     \item \code{dnaseq} — tissue alteration calls
#'     \item \code{expression} — long-format expression frame
#'     \item \code{tmb} — per-sample TMB frame
#'   }
#' @examples
#' sim  <- ctdna_make_mock_study(n_patients = 15, seed = 1)
#' prep <- ctdna_prepare(sim, verbose = FALSE)
#' names(prep)        # ctdna, mutations, tmb, ihc, expression, ...
#' head(prep$samples)
#'
#' # With pre-filtering applied to mutation/CNV data
#' prep <- ctdna_prepare(sim, filter_scheme = "default", verbose = FALSE)
#' @export


# -----------------------------------------------------------------------------
# v0.42.8 / v0.42.9 / v0.42.10 diagnostic + canonicalization helpers
# -----------------------------------------------------------------------------

# v0.42.10: column-name aliases. For each canonical column name the platform
# handler reads, the value is a character vector of accepted aliases.
# Matching is CASE-INSENSITIVE (the first alias is the canonical form a
# matched column is renamed TO).
.ctdna_column_aliases <- function() {
  list(
    Patient_ID = c("Patient_ID","patient_id","PatientID","patientid",
                   "Subject_ID","subject_id","SubjectID","subjectid",
                   "SUBJID","subjid","SubjID","SUBJ_ID","Subj_ID",
                   "USUBJID","usubjid","UniqueSubjID","Unique_Subject_ID",
                   "PT_ID","pt_id","PtID","ptid"),
    Sample_ID  = c("Sample_ID","sample_id","SampleID","sampleid",
                   "SAMPLEID","SAMPLE_ID",
                   "GHSampleID","GHsampleid","GH_Sample_ID"),
    Visit_name = c("Visit_name","visit_name","VisitName","visitname",
                   "Visit","visit","VISIT","VisitID","Visit_ID",
                   "VISITN","AVISIT","avisit","AVISITN","avisitn",
                   "Timepoint","timepoint","Time_point","TIME_POINT",
                   "Cycle_Day","cycle_day","TIMEPOINT","TP"),
    Sample_status = c("Sample_status","sample_status","SampleStatus",
                      "samplestatus","QC_status","qc_status","QCStatus",
                      "Status","status","STATUS","Sample_QC"),
    Gene = c("Gene","gene","GENE","gene_symbol","Gene_symbol","GeneSymbol",
             "Hugo_Symbol","HUGO_SYMBOL","HugoSymbol","Symbol","SYMBOL"),
    Variant_type = c("Variant_type","variant_type","VariantType",
                     "Variant_Type","VARTYPE","VarType","variant",
                     "Mutation_type","mutation_type","MutationType",
                     "Alt_type","alt_type"),
    Alteration = c("Alteration","alteration","ALTERATION","Variant",
                   "VARIANT","Mutation","mutation","Alt"),
    Copy_number = c("Copy_number","copy_number","CopyNumber","copynumber",
                    "CN","cn","COPY_NUMBER","Copies"),
    VAF_percentage = c("VAF_percentage","vaf_percentage","VAF",
                       "vaf","VAFPercentage","VAF_pct","VAF_pct100",
                       "Variant_allele_frequency"),
    TF_change = c("TF_change","tf_change","TFChange","TumorFraction_change",
                  "TF_delta"),
    Methylation_tumor_fraction_percentage = c(
                 "Methylation_tumor_fraction_percentage",
                 "methylation_tumor_fraction_percentage",
                 "Methylation_tumor_fraction",
                 "methylTF","methyl_tf","methylationTF","methylation_TF",
                 "methyl_TF","mTF"),
    Genomic_max_VAF_percentage = c(
               "Genomic_max_VAF_percentage","genomic_max_vaf_percentage",
               "Genomic_max_VAF",
               "Max_VAF_percentage","max_vaf_percentage",
               "maxVAF","max_vaf","maxvaf","Max_VAF","MAX_VAF"),
    Mean_VAF_percentage = c(
               "Mean_VAF_percentage","mean_vaf_percentage",
               "meanVAF","mean_vaf","meanvaf","Mean_VAF","MEAN_VAF"),
    TMB_score = c("TMB_score","tmb_score","TMB","tmb","TMB_Score","TMBScore",
                  "Tumor_Mutational_Burden"),
    HRD_score = c("HRD_score","hrd_score","HRD","hrd","HRD_Score","HRDScore",
                  "Genomic_HRD_score","genomic_hrd_score","HRD_status","hrd_status",
                  "Homologous_Recombination_Deficiency"),
    Dose = c("Dose","dose","DOSE","Dose_level","dose_level","DoseLevel",
             "DOSE_LEVEL","Dose_Level",
             # v0.43.7: ARM (descriptive dose label) wins over
             # ARMCD (short code: "B", "C"). Per-modality rename can
             # override (e.g. clinical = list(data=df, dose_col="ARMCD")).
             "ARM","arm",
             "TRT01P","trt01p","TRT01A","trt01a",
             "ACTARM","actarm",
             "Treatment","treatment",
             "indication_dose","Indication_Dose",
             "ARMCD","armcd","ACTARMCD","actarmcd"),
    RECIST = c("RECIST","recist","Recist","RECIST_response","recist_response",
               "BestResponse","best_response","BOR","bor","BORS",
               "Best_Overall_Response","best_overall_response",
               "OBJECTIVE_RESPONSE","Objective_response","Response_RECIST",
               "AVALC"),  # AVALC sometimes used in ADRS ADaM
    Sex = c("Sex","sex","SEX","Gender","gender","GENDER","sx"),
    Age = c("Age","age","AGE","AgeAtBaseline","age_at_baseline","AGEDY"),
    Indication = c("Indication","indication","INDICATION","disease",
                   "Disease","Cancer_type","cancer_type","TUMTYPE","tumtype",
                   "Tumor_type","CANCER_TYPE","TumorType")
  )
}

# v0.42.12: resolve a user-provided rename key (e.g. "patient_id_col",
# "patient_id", "PatientID", "subject_id", "SUBJID") to a canonical
# column name. Returns NA_character_ if no canonical can be matched.
.ctdna_resolve_rename_key <- function(key, aliases = NULL) {
  if (is.null(aliases)) aliases <- .ctdna_column_aliases()
  if (is.null(key) || !nzchar(key)) return(NA_character_)

  # Strip optional "_col" suffix
  bare_key <- sub("_col$", "", key, ignore.case = TRUE)

  # 1) Exact match against canonical names
  if (key %in% names(aliases)) return(key)
  if (bare_key %in% names(aliases)) return(bare_key)

  # 2) Case-insensitive match against canonical names
  ci_canon <- which(tolower(names(aliases)) == tolower(bare_key))
  if (length(ci_canon) == 1L) return(names(aliases)[ci_canon])
  # Try the key as-is too (without _col strip)
  ci_canon2 <- which(tolower(names(aliases)) == tolower(key))
  if (length(ci_canon2) == 1L) return(names(aliases)[ci_canon2])

  # 3) Case-insensitive match against any alias for any canonical
  for (canonical in names(aliases)) {
    alts_lc <- tolower(aliases[[canonical]])
    if (tolower(bare_key) %in% alts_lc) return(canonical)
    if (tolower(key)       %in% alts_lc) return(canonical)
  }

  NA_character_
}

# v0.42.12: split a modality slot input into a (data, renames) spec.
#   - data frame / coercible -> list(data = df, renames = char(0))
#   - list with `data` key + other named entries -> list(data, renames)
.ctdna_extract_modality_spec <- function(x, slot_name = "<slot>") {
  if (is.null(x)) return(NULL)

  # Pure data frame or coercible (matrix, tibble, etc.)
  if (is.data.frame(x))
    return(list(data = x, renames = character(0)))
  if (is.matrix(x))
    return(list(data = .ctdna_to_data_frame(x, slot_name),
                renames = character(0)))

  # List-form spec: must have `data` key
  if (is.list(x)) {
    if (!"data" %in% names(x)) {
      # Could be a tibble or arrow Table (list-like under the hood) - coerce
      tryCatch({
        df <- .ctdna_to_data_frame(x, slot_name)
        return(list(data = df, renames = character(0)))
      }, error = function(e) NULL)
      stop(sprintf(
        paste0("ctdna_prepare: `%s` is a list without a `data` key. ",
               "To pass a frame + per-modality renames, use the spec form: ",
               "%s = list(data = your_frame, patient_id_col = \"SUBJID\", ...). ",
               "To pass a bundle of multiple modalities, use the top-level ",
               "`bundle = list(...)` argument."),
        slot_name, slot_name), call. = FALSE)
    }
    df <- .ctdna_to_data_frame(x$data, paste0(slot_name, "$data"))

    rename_keys <- setdiff(names(x), "data")
    renames <- character(0)
    for (k in rename_keys) {
      v <- x[[k]]
      if (is.null(v)) next
      if (!is.character(v) || length(v) != 1L || !nzchar(v))
        stop(sprintf(
          "ctdna_prepare: `%s$%s` must be a single non-empty string (the source column name to rename FROM). Got: %s",
          slot_name, k, paste(class(v), collapse="/")), call. = FALSE)
      renames[k] <- v
    }
    return(list(data = df, renames = renames))
  }

  stop(sprintf(
    "ctdna_prepare: `%s` must be a data frame, matrix, tibble, or a list ",
    "with a `data` key. Got class %s.",
    slot_name, paste(class(x), collapse="/")), call. = FALSE)
}


# v0.42.11: priority-based canonicalization. The alias list for each
# canonical column is ALSO a priority order (first = highest priority).
# When multiple columns match the same canonical, the highest-priority
# alias wins. A warning() is emitted so the user knows which was chosen.
#
# Optional `patient_id_override` lets the caller force a specific column
# (by literal name) to be used as Patient_ID — useful when the auto-pick
# would choose the wrong one or when the column has a custom name not
# in the alias list.
.ctdna_canonicalize_columns <- function(df, frame_name = "<frame>",
                                          patient_id_override = NULL,
                                          renames = character(0),
                                          verbose = TRUE) {
  if (!is.data.frame(df) || ncol(df) == 0L) return(df)
  aliases <- .ctdna_column_aliases()
  nms     <- names(df)
  nms_lc  <- tolower(nms)
  renames_done <- character(0)

  # --- Step 0 (v0.42.12): apply per-modality renames first ---
  # These come from the list-form spec, e.g.
  #   infinity_report = list(data = df, patient_id_col = "SUBJID", maxvaf = "MaxAF")
  # Each (target_key = source_col) maps the source column in df to the
  # canonical column resolved from target_key. Per-modality renames have
  # the HIGHEST precedence (above patient_id_override and auto-detection).
  if (length(renames) > 0L) {
    for (target_key in names(renames)) {
      source_col <- renames[[target_key]]
      canonical  <- .ctdna_resolve_rename_key(target_key, aliases)
      if (is.na(canonical)) {
        warning(sprintf(
          paste0("ctdna_prepare: in `%s`, rename key `%s` doesn't match any ",
                 "known canonical column. Known canonicals: %s. Skipping."),
          frame_name, target_key,
          paste(names(aliases), collapse=", ")),
          call. = FALSE)
        next
      }
      idx <- which(nms_lc == tolower(source_col))
      if (length(idx) == 0L) {
        warning(sprintf(
          paste0("ctdna_prepare: in `%s`, rename source column `%s` ",
                 "(for canonical `%s`) not found. Skipping."),
          frame_name, source_col, canonical),
          call. = FALSE)
        next
      }
      old_name <- nms[idx[1]]
      if (old_name == canonical) next  # already canonical, no-op

      # If a different column already has the canonical name, shadow it
      if (canonical %in% nms && nms[which(nms == canonical)[1]] != old_name) {
        shadow_idx  <- which(nms == canonical)[1]
        shadow_name <- paste0(canonical, "__shadowed_by_", old_name)
        names(df)[shadow_idx] <- shadow_name
        if (verbose)
          message(sprintf(
            " - note: in `%s`, existing column `%s` shadowed by user rename (`%s` -> `%s`); kept as `%s`",
            frame_name, canonical, old_name, canonical, shadow_name))
        nms     <- names(df)
        nms_lc  <- tolower(nms)
        idx     <- which(nms == old_name)  # reacquire index of source col
      }

      names(df)[idx[1]] <- canonical
      renames_done <- c(renames_done, setNames(old_name, canonical))
      if (verbose)
        message(sprintf(" - applied per-modality rename in `%s`: `%s` -> `%s`",
                        frame_name, old_name, canonical))
      nms     <- names(df)
      nms_lc  <- tolower(nms)
    }
  }

  # --- Step 1: honor explicit (global) patient_id_override ---
  if (!is.null(patient_id_override) && nzchar(patient_id_override) &&
      !"Patient_ID" %in% renames_done) {   # skip if per-modality already set Patient_ID
    hit_idx <- which(nms_lc == tolower(patient_id_override))
    if (length(hit_idx) == 0L) {
      if (verbose)
        message(sprintf(" - note: patient_id_col = `%s` was specified but not found in `%s`. Falling back to alias-based detection.",
                        patient_id_override, frame_name))
    } else {
      old_name <- nms[hit_idx[1]]
      if (old_name != "Patient_ID") {
        if ("Patient_ID" %in% nms) {
          existing_idx <- which(nms == "Patient_ID")[1]
          shadow_name <- paste0("Patient_ID__shadowed_by_", old_name)
          names(df)[existing_idx] <- shadow_name
          if (verbose)
            message(sprintf(" - note: in `%s`, existing `Patient_ID` shadowed by patient_id_col override `%s` (kept as `%s`)",
                            frame_name, old_name, shadow_name))
          nms <- names(df)
          hit_idx <- which(nms == old_name)
        }
        names(df)[hit_idx[1]] <- "Patient_ID"
        renames_done <- c(renames_done, setNames(old_name, "Patient_ID"))
        if (verbose)
          message(sprintf(" - canonicalized `Patient_ID` in `%s`: forced `%s` -> `Patient_ID` via patient_id_col override",
                          frame_name, old_name))
        nms     <- names(df)
        nms_lc  <- tolower(nms)
      }
    }
  }

  # --- Step 2: priority-based alias matching for everything else ---
  for (canonical in names(aliases)) {
    if (canonical %in% nms) next
    alts    <- aliases[[canonical]]
    alts_lc <- tolower(alts)

    matches_in_priority <- integer(0)
    for (j in seq_along(alts_lc)) {
      hit <- which(nms_lc == alts_lc[j])
      if (length(hit) > 0L && !(hit[1] %in% matches_in_priority)) {
        matches_in_priority <- c(matches_in_priority, hit[1])
      }
    }
    if (length(matches_in_priority) == 0L) next

    chosen_idx <- matches_in_priority[1]
    other_idx  <- matches_in_priority[-1]
    old_name   <- nms[chosen_idx]

    if (length(other_idx) > 0L) {
      other_nms <- nms[other_idx]
      warning(sprintf(
        paste0("ctdna_prepare: in frame `%s`, multiple columns match canonical ",
               "`%s`: %s. Using `%s` (highest priority); the others are kept ",
               "under their original names. Pass a per-modality rename like ",
               "`%s = list(data = ..., %s_col = '<name>')` to override."),
        frame_name, canonical,
        paste(c(old_name, other_nms), collapse = ", "),
        old_name, frame_name, tolower(canonical)),
        call. = FALSE)
    }

    names(df)[chosen_idx] <- canonical
    renames_done <- c(renames_done, setNames(old_name, canonical))
    if (verbose && length(other_idx) == 0L)
      message(sprintf(" - canonicalized `%s` in `%s`: `%s` -> `%s`",
                      canonical, frame_name, old_name, canonical))
    nms     <- names(df)
    nms_lc  <- tolower(nms)
  }

  if (length(renames_done) > 0L)
    attr(df, "ctdna_renames") <- renames_done
  df
}


# v0.42.11: aggressive coercion to plain data.frame. Accepts tibbles,
# data.tables, matrices, list-of-vectors — anything as.data.frame() handles.
# Strips tibble/data.table subclasses (downstream relies on base subsetting).
.ctdna_to_data_frame <- function(x, name = "<unknown>") {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) {
    if (!identical(class(x), "data.frame")) {
      # Strip tibble/data.table/etc. subclass
      x <- as.data.frame(x, stringsAsFactors = FALSE)
    }
    return(x)
  }
  if (is.matrix(x) || is.list(x) || inherits(x, c("tbl","tbl_df","data.table","arrow_table","Table"))) {
    out <- tryCatch(
      as.data.frame(x, stringsAsFactors = FALSE),
      error = function(e) {
        stop(sprintf(
          "ctdna_prepare: could not convert `%s` (class %s) to a data.frame: %s",
          name, paste(class(x), collapse="/"), conditionMessage(e)),
          call. = FALSE)
      })
    return(out)
  }
  stop(sprintf(
    "ctdna_prepare: `%s` must be a data frame or convertible to one (got class %s).",
    name, paste(class(x), collapse="/")), call. = FALSE)
}


# Modality expectations table. For each modality slot, we encode:
#   `required`  - canonical columns that MUST be present (after canonicalization)
#   `essential` - canonical columns the platform handler actually reads
#   `typical`   - other expected columns
#   `desc`      - human-readable description
.ctdna_modality_specs <- function() {
  list(
    infinity_report = list(
      # Without these, .gh_build_genomic / .gh_build_ctdna can't run at all:
      required  = c("Patient_ID"),
      # Without these, the platform handler produces empty output:
      essential = c("Visit_name", "Sample_status", "Variant_type", "Gene"),
      typical   = c("Alteration", "TF_change", "TMB_score", "VAF_percentage",
                    "Copy_number"),
      desc      = "GH Infinity variant/sample report (per-row per-variant)"
    ),
    tf_change = list(
      required  = c("Patient_ID"),
      essential = c("Visit_name"),
      typical   = c("TF_change","methylTF"),
      desc      = "Tumor-fraction time series"
    ),
    panel_74_response = list(
      required  = c("Patient_ID"),
      essential = c("Visit_name"),
      typical   = c("panel_74_response","Response"),
      desc      = "Panel-74 longitudinal response table"
    ),
    panel_500_response = list(
      required  = c("Patient_ID"),
      essential = c("Visit_name"),
      typical   = c("panel_500_response","Response"),
      desc      = "Panel-500 longitudinal response table"
    ),
    methylation_response = list(
      required  = c("Patient_ID"),
      essential = c("Visit_name"),
      typical   = c("methylation_response","methylTF"),
      desc      = "Methylation longitudinal response table"
    ),
    clinical = list(
      required  = c("Patient_ID"),
      essential = character(0),
      typical   = c("Dose","RECIST","Sex","Age","Indication"),
      desc      = "Per-patient clinical / demographics frame"
    ),
    ihc = list(
      required  = c("Patient_ID"),
      essential = character(0),
      typical   = c("Marker","IHC_score","Visit_name"),
      desc      = "IHC marker frame"
    ),
    dnaseq = list(
      required  = c("Patient_ID"),
      essential = c("Gene"),
      typical   = c("Variant_type","Alteration"),
      desc      = "DNA sequencing variant frame (e.g. WES)"
    ),
    wes = list(
      required  = c("Patient_ID"),
      essential = c("Gene"),
      typical   = c("Variant_type","Alteration"),
      desc      = "Whole-exome sequencing variant frame (alias of dnaseq)"
    ),
    rnaseq = list(
      required  = c("Patient_ID"),
      essential = c("Gene"),
      typical   = c("expression","log2_tpm","TPM"),
      desc      = "RNA expression frame (long format)"
    )
  )
}

# Walk each modality frame and verify column shape. Hard-errors when `strict`,
# otherwise warns. Returns invisibly the list of diagnostics for testing.
.ctdna_check_modality_shapes <- function(bundle, strict = TRUE, verbose = TRUE) {
  specs <- .ctdna_modality_specs()
  diag  <- list()
  problems <- character(0)

  for (k in names(bundle)) {
    df  <- bundle[[k]]
    sp  <- specs[[k]]
    if (is.null(sp)) {
      if (verbose)
        message(" - note: '", k, "' is not a known modality slot; passing ",
                "through. Known slots: ", paste(names(specs), collapse = ", "), ".")
      next
    }
    if (!is.data.frame(df)) {
      msg <- sprintf("frame passed as `%s` is not a data frame (got: %s).",
                     k, class(df)[1])
      problems <- c(problems, msg)
      next
    }

    nms       <- names(df)
    miss_req  <- setdiff(sp$required, nms)
    miss_ess  <- setdiff(sp$essential, nms)
    has_typ   <- intersect(sp$typical, nms)

    if (length(miss_req) > 0L) {
      problems <- c(problems, sprintf(
        "frame passed as `%s` is missing REQUIRED column(s): %s. (Modality: %s.)",
        k, paste(miss_req, collapse = ", "), sp$desc))
      next
    }

    # Essential columns: needed for the platform handler to produce anything.
    # If MORE THAN HALF are missing, this frame can't really be processed as
    # this modality. (One missing essential = warning; many missing = fatal.)
    if (length(miss_ess) > 0L) {
      missing_ratio <- length(miss_ess) / length(sp$essential)
      if (missing_ratio > 0.5) {
        problems <- c(problems, sprintf(
          paste0("frame passed as `%s` is missing %d of %d ESSENTIAL ",
                 "column(s) [%s]. The platform handler will produce empty ",
                 "output. Is this the right frame for the `%s` slot? ",
                 "Expected modality: %s."),
          k, length(miss_ess), length(sp$essential),
          paste(miss_ess, collapse = ", "), k, sp$desc))
        next
      } else if (verbose) {
        message(" - note: frame for `", k, "` is missing some essential ",
                "columns (", paste(miss_ess, collapse = ", "),
                ") but may still be partially usable.")
      }
    }

    if (length(has_typ) == 0L && length(sp$typical) > 0L) {
      problems <- c(problems, sprintf(
        paste0("frame passed as `%s` has NONE of the typical columns for ",
               "this modality. Is this the right frame? Expected typical ",
               "columns: %s. (Modality: %s.)"),
        k, paste(sp$typical, collapse = ", "), sp$desc))
      next
    }

    diag[[k]] <- list(status = "ok",
                       found_typical   = has_typ,
                       missing_typical = setdiff(sp$typical, has_typ),
                       missing_essential = miss_ess)
  }

  if (length(problems) > 0L) {
    msg <- paste0("ctdna_prepare: input shape validation failed.\n  - ",
                  paste(problems, collapse = "\n  - "),
                  "\n  See ?ctdna_prepare for the column shape each modality ",
                  "expects. To downgrade these errors to warnings, pass ",
                  "`strict = FALSE`.")
    if (strict) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }

  invisible(diag)
}

# Detect when the same R object was passed as more than one modality
# (e.g. ctdna_prepare(infinity_report = df, clinical = df) with `df` identical).
# This is almost always a user mistake: the same wide table cannot serve as
# two different modality frames.
.ctdna_check_duplicate_frames <- function(bundle, strict = TRUE, verbose = TRUE) {
  if (length(bundle) < 2L) return(invisible(NULL))
  keys <- names(bundle)
  seen <- list()
  for (i in seq_along(bundle)) {
    matched <- FALSE
    for (g in seq_along(seen)) {
      rep_key <- seen[[g]][1]
      if (identical(bundle[[keys[i]]], bundle[[rep_key]])) {
        seen[[g]] <- c(seen[[g]], keys[i])
        matched <- TRUE
        break
      }
    }
    if (!matched) seen[[length(seen) + 1L]] <- keys[i]
  }
  dups <- Filter(function(g) length(g) >= 2L, seen)
  if (length(dups) == 0L) return(invisible(NULL))

  msg_lines <- vapply(dups, function(grp) sprintf(
    "the same data frame was passed under slots: %s",
    paste(grp, collapse = ", ")), character(1))
  msg <- paste0(
    "ctdna_prepare: duplicate-frame detection failed.\n  - ",
    paste(msg_lines, collapse = "\n  - "),
    "\n  Each modality slot expects a DIFFERENT, narrow per-modality frame. ",
    "If you have a wide table containing multiple modalities mixed together ",
    "(e.g. Master_Clinical with clinical + infinity + genomic columns in one), ",
    "split it into per-modality subsets BEFORE calling ctdna_prepare(). ",
    "See ?ctdna_prepare for the column shape each modality expects. ",
    "To downgrade this to a warning, pass `strict = FALSE`.")

  if (strict) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  invisible(NULL)
}




# In v0.40.0 these renamed Patient_ID -> ctdna_opts("subject"). In v0.41.0
# the default subject is "Patient_ID" so no rename is needed — these helpers
# now just hand the vendor frame back unchanged. The rename surface, if a
# user really needs it, is `new_col_names` in ctdna_prepare().

.prep_ihc <- function(df, c_subj) {
  if (!is.data.frame(df)) return(df)
  df
}

.prep_dnaseq <- function(df, c_subj) {
  if (!is.data.frame(df)) return(df)
  df
}

# Reshape the rnaseq list (counts, tpm, sample_metadata) into a long-form
# expression data frame with columns named per ctdna_opts():
#   .o("subject"), .o("gene"), .o("expression"), and a sample_id column.
.prep_expression <- function(rnaseq, c_subj, c_gene, c_expr) {
  if (!is.list(rnaseq) || is.null(rnaseq$tpm))
    return(NULL)
  tpm <- rnaseq$tpm
  meta <- rnaseq$sample_metadata
  n_genes <- nrow(tpm)
  n_samples <- ncol(tpm)
  long <- data.frame(
    Sample_ID = rep(colnames(tpm), each = n_genes),
    Gene      = rep(rownames(tpm), times = n_samples),
    expr      = log2(as.vector(tpm) + 1),
    stringsAsFactors = FALSE
  )
  if (is.data.frame(meta) && "Sample_ID" %in% names(meta) &&
      "Patient_ID" %in% names(meta)) {
    long$Patient_ID <- meta$Patient_ID[match(long$Sample_ID, meta$Sample_ID)]
  } else {
    long$Patient_ID <- long$Sample_ID
  }
  # Rename to the configured opts names (subject / gene / expression)
  names(long)[names(long) == "Patient_ID"] <- c_subj
  names(long)[names(long) == "Gene"]       <- c_gene
  names(long)[names(long) == "expr"]       <- c_expr
  attr(long, "units") <- "log2_tpm_plus_one"
  long
}

# Per-sample TMB pulled straight from the raw infinity_report. Pass through
# vendor names (Patient_ID, Sample_ID, Visit_name, TMB_score) by default;
# rename only when ctdna_opts() has been pointed elsewhere.
.prep_tmb <- function(inf_report, c_subj) {
  cols <- intersect(c("Patient_ID","Sample_ID","Visit_name","TMB_score"),
                     names(inf_report))
  if (!"TMB_score" %in% cols) return(NULL)
  df <- unique(inf_report[, cols, drop = FALSE])
  if ("Patient_ID" %in% names(df) && c_subj != "Patient_ID")
    names(df)[names(df) == "Patient_ID"] <- c_subj
  rownames(df) <- NULL
  df
}

#' List ctDNA platforms supported by `ctdna_prepare()`
#'
#' Returns the registry of platform handlers known to
#' \code{\link{ctdna_prepare}}. Each row is a `(platform, name)` pair
#' the user can pass as `platform =` to dispatch to a vendor-specific
#' preparation routine (currently `"gh_infinity"` for Guardant Infinity
#' Vantage). Mostly diagnostic — call when adapting the package to a
#' new study and you want to confirm the registry is loaded.
#'
#' @return A data frame with `platform` (the vendor_technology key) and
#'   `name` (human-readable label) columns.
#' @examples
#' # See which platforms ctdna_prepare() supports
#' ctdna_platforms()
#' @export
ctdna_platforms <- function() {
  if (length(.platforms) == 0)
    return(data.frame(platform = character(), name = character(),
                      stringsAsFactors = FALSE))
  data.frame(
    platform = names(.platforms),
    name     = vapply(.platforms, `[[`, character(1), "name"),
    stringsAsFactors = FALSE
  )
}


#' Get the path to the bundled ctdnaTM logo
#'
#' Returns the on-disk path to the logo PNG that ships with the package
#' (in `inst/extdata`). Use it in Rmd headers / chunks without having
#' to remember the full \code{system.file()} call.
#'
#' Two common idioms:
#' \preformatted{
#' # In an Rmd chunk:
#' knitr::include_graphics(ctdna_logo())
#'
#' # As Markdown inline image (paste the path the function returns):
#' cat(sprintf("![](\%s)", ctdna_logo()))
#' }
#'
#' To replace the logo, overwrite the file at the returned path with
#' your own PNG. The change takes effect immediately (no rebuild
#' needed).
#'
#' @return Character path to the bundled logo PNG, or NA if the file
#'   is missing from the installation.
#' @examples
#' ctdna_logo()
#' @export
ctdna_logo <- function() {
  p <- system.file("extdata", "logo.png", package = "ctdnaTM")
  if (!nzchar(p)) return(NA_character_)
  p
}
