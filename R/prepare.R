# =============================================================================
# ctdnaTM v0.44.0 -- unified preparation
#
# ctdna_prepare() gathers whatever modalities the user has, canonicalizes
# column names to ctdna_opts() values, and returns a classed `ctdna_prep`
# list with TWO ctDNA frames keyed by their row grain:
#
#   prep$samples   one row per Patient_ID x Visit_name x Sample_ID
#                  (TF%, VAF%, TMB, QC, Dose, RECIST, ...)
#   prep$variants  one row per Patient_ID x Visit_name x Sample_ID x variant
#                  (Gene, Variant_type, VAF_percentage, annotation, ...)
#   prep$clinical  the clinical frame (per patient / per visit)
#   prep$dictionary  variable / column / frame / source / units / purpose
#   prep$<modality>  rnaseq / dnaseq / ihc / ... (added later)
#
# Per-sample QC used to run inside prepare (Sample_QC arg, removed in v0.75.0) and CASCADES: a sample
# that fails QC is removed from $samples AND all its variants from $variants.
# Variant *filtering* (the rule/scheme system) is separate and only ever
# touches $variants -- it never drops a sample.
# =============================================================================


# ---- Registry: single source of truth --------------------------------------
# Drives the dictionary. One row per canonical variable that may appear in the
# ctDNA frames. `column` is resolved live from ctdna_opts() so renaming an opt
# updates the dictionary automatically.
.ctdna_registry <- function() {
  o <- function(k) tryCatch(.o(k), error = function(e) NA_character_)
  r <- function(variable, column, frame, grain, sources, units, purpose)
    data.frame(variable = variable, column = column, frame = frame,
               grain = grain, source = sources, units = units,
               purpose = purpose, stringsAsFactors = FALSE)
  rbind(
    r("Patient",                o("subject"), "samples,variants", "key",
      "all", "id", "join key; every analysis"),
    r("Visit",                  o("time"), "samples,variants", "key",
      "infinity>tf_change>panels", "label", "longitudinal axis; baseline detection"),
    r("Sample",                 o("sample"), "samples,variants", "key",
      "infinity>tf_change>panels", "id", "per-sample key; QC cascade key"),
    r("Baseline flag",          "Baseline", "samples,variants", "per-sample",
      "derived", "logical", "marks baseline draw; genomic baseline views"),
    r("Methylation TF",         o("tf"), "samples", "per-sample",
      "tf_change>infinity", "percent 0-100", "TF plots, mr_recist, mut_methyl"),
    r("Genomic max VAF",        o("maxvaf"), "samples", "per-sample",
      "infinity>panel500>panel74", "percent 0-100", "reduction, longitudinal"),
    r("Mean VAF",               o("meanvaf"), "samples", "per-sample",
      "infinity(agg)>panels", "percent 0-100", "mut_methyl"),
    r("Somatic variant count",  o("n_somatic"), "samples", "per-sample",
      "infinity(agg)", "count", "burden summaries"),
    r("MR score",               o("mr"), "samples", "per-sample",
      "tf_change>methylation>panels", "percent 0-100", "mr_recist"),
    r("ctDNA % change",         "Ctdna_percentage_change", "samples", "per-sample",
      "tf_change>panels", "percent", "pct_change_by_dose_time"),
    r("TMB",                    o("tmb"), "samples", "per-sample",
      "infinity", "mut/Mb", "TMB plots, mut_methyl"),
    r("cfDNA input",            o("cfdna_conc"), "samples", "per-sample",
      "infinity", "ng", "sample QC"),
    r("Plasma input",           "Plasma_ml_input", "samples", "per-sample",
      "infinity", "mL", "sample QC"),
    r("Sample status",          o("col_sample_status"), "samples", "per-sample",
      "infinity", "label", "sample QC"),
    r("Dose",                   o("dose"), "samples,clinical", "per-patient",
      "clinical(adsl)>ctdna", "label", "faceting / grouping"),
    r("RECIST",                 o("recist"), "samples,clinical", "per-patient/visit",
      "clinical(adrs)>ctdna", "category", "response plots"),
    r("Cancer type",            o("cancertype"), "samples,clinical", "per-patient",
      "clinical>ctdna", "label", "indication faceting"),
    r("Gene",                   o("gene"), "variants", "per-variant",
      "infinity", "symbol", "oncoprint, alteration_grid"),
    r("Variant type",           o("alteration"), "variants", "per-variant",
      "infinity", "category", "oncoprint, alteration_grid"),
    r("Variant VAF",            o("vaf"), "variants", "per-variant",
      "infinity", "percent 0-100", "oncoprint VAF, alteration_grid"),
    r("Copy number",            o("copy_number"), "variants", "per-variant",
      "infinity", "integer", "CNV calls")
  )
}

# Canonical per-sample columns that prepare assembles into $samples (in order).
.ctdna_sample_cols <- function() {
  c(.o("subject"), .o("sample"), .o("time"), "Visit_name_raw", "Baseline",
    .o("tf"), .o("maxvaf"), .o("meanvaf"), .o("n_somatic"),
    .o("mr"), "Ctdna_percentage_change",
    .o("tmb"), "TMB_category", "MSI_High", "HRD_score", "ctDNA_detection_status",
    .o("col_sample_status"), .o("cfdna_conc"), "Plasma_ml_input",
    .o("cancertype"),
    .o("dose"), .o("recist"), .o("best_change"), .o("lesion"))
}

# ---- small helpers ----------------------------------------------------------
.coalesce2 <- function(primary, secondary) {
  if (is.null(primary))   return(secondary)
  if (is.null(secondary)) return(primary)
  out <- primary
  na <- is.na(out)
  out[na] <- secondary[na]
  out
}

# robust baseline detector: "Cycle 1 Day 1", "C1D1", "Baseline", "Screening"
.ctdna_is_baseline_label <- function(x) {
  if (is.null(x)) return(logical(0))
  xx <- tolower(trimws(as.character(x)))
  grepl("^cycle\\s*0*1\\s*day\\s*0*1$", xx) |
    grepl("^c0*1d0*1$", xx) |
    xx %in% c("baseline", "screening", "c1d1", "cycle 1 day 1")
}

# melt a paired (A/B) GH response frame into per-sample long form
.melt_paired <- function(df, verbose = FALSE) {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(NULL)
  nms    <- names(df)
  a_cols <- grep("_A$", nms, value = TRUE)
  stems  <- sub("_A$", "", a_cols)
  stems  <- stems[paste0(stems, "_B") %in% nms]
  if (!length(stems)) return(NULL)
  pair_level <- setdiff(nms, c(paste0(stems, "_A"), paste0(stems, "_B")))
  carry_keys <- intersect(c("Patient_ID", "Study_ID"), pair_level)
  attach_b   <- setdiff(pair_level, carry_keys)   # pair-level -> follow-up (B)

  legs <- lapply(c("A", "B"), function(L) {
    sub <- df[, paste0(stems, "_", L), drop = FALSE]
    names(sub) <- stems
    pl <- df[, pair_level, drop = FALSE]
    if (L == "A" && length(attach_b))
      pl[attach_b] <- lapply(pl[attach_b], function(z) z[NA_integer_ + seq_along(z)])
    out <- cbind(pl, sub, .leg = L, stringsAsFactors = FALSE)
    out
  })
  long <- do.call(rbind, legs)
  rownames(long) <- NULL
  long
}

# Build a Sample_ID-keyed lookup of one column from a per-sample contribution.
.by_sample <- function(df, col, key = "Sample_ID") {
  if (is.null(df) || !is.data.frame(df) || !col %in% names(df) ||
      !key %in% names(df)) return(NULL)
  d <- df[!is.na(df[[key]]), c(key, col)]
  d <- d[!duplicated(d[[key]]) | !is.na(d[[col]]), , drop = FALSE]
  # prefer non-NA value per key
  d <- d[order(d[[key]], is.na(d[[col]])), ]
  d <- d[!duplicated(d[[key]]), ]
  stats::setNames(d[[col]], d[[key]])
}

.coalesce_from <- function(samples, lookup, target, key = "Sample_ID",
                           primary_is_samples = TRUE) {
  if (is.null(lookup) || !length(lookup)) return(samples)
  inc <- unname(lookup[as.character(samples[[key]])])
  if (!target %in% names(samples)) samples[[target]] <- NA
  samples[[target]] <- if (primary_is_samples)
    .coalesce2(samples[[target]], inc) else .coalesce2(inc, samples[[target]])
  samples
}


# ---- per-variant frame ------------------------------------------------------
.ctdna_build_variants <- function(inf) {
  v <- inf
  v$Baseline <- .ctdna_is_baseline_label(v[[.o("time")]])
  # keep variant rows (a detected alteration); samples with no variant are
  # represented in $samples, not here.
  g <- .o("gene")
  if (g %in% names(v)) v <- v[!is.na(v[[g]]) & nzchar(as.character(v[[g]])), , drop = FALSE]
  rownames(v) <- NULL
  v
}

# ---- per-sample frame -------------------------------------------------------
.ctdna_build_samples <- function(inf, tf, m500, m74, methyl, verbose = TRUE) {
  c_subj <- .o("subject"); c_samp <- .o("sample"); c_time <- .o("time")

  # 1) spine: one row per Sample_ID from the infinity report (per-sample cols
  #    are repeated across that sample's variant rows -> distinct collapses).
  keep <- intersect(.ctdna_sample_cols(), names(inf))
  spine <- inf[, unique(c(c_subj, c_samp, c_time, keep)), drop = FALSE]
  spine <- spine[!duplicated(spine[[c_samp]]), , drop = FALSE]

  # 2) derive per-sample aggregates from inf variant rows
  if (.o("vaf") %in% names(inf)) {
    vv <- suppressWarnings(as.numeric(inf[[.o("vaf")]]))
    agg_mean <- tapply(vv, inf[[c_samp]], function(z) mean(z, na.rm = TRUE))
    agg_n    <- tapply(inf[[.o("gene")]], inf[[c_samp]],
                       function(z) sum(!is.na(z) & nzchar(as.character(z))))
    spine <- .coalesce_from(spine, agg_mean, .o("meanvaf"), c_samp, TRUE)
    spine[[.o("n_somatic")]] <- unname(agg_n[as.character(spine[[c_samp]])])
  }

  # 3) melt response frames -> per-sample contributions, canonicalize names
  melt_one <- function(df, nm) {
    m <- .melt_paired(df)
    if (is.null(m)) return(NULL)
    m <- .ctdna_canonicalize_columns(m, nm, verbose = FALSE)
    if ("Visit_name" %in% names(m)) m[[c_time]] <- m[["Visit_name"]]
    m$Baseline <- m$.leg == "A"
    m
  }
  tf_m     <- melt_one(tf,     "tf_change")
  m500_m   <- melt_one(m500,   "panel_500_response")
  m74_m    <- melt_one(m74,    "panel_74_response")
  methyl_m <- melt_one(methyl, "methylation_response")

  # 4) coalesce per source priority (clinical handled later; among ctDNA:
  #    tf_change > infinity > genomic(500>74); infinity already in spine)
  # methylTF: tf_change WINS over infinity (spine)
  spine <- .coalesce_from(spine, .by_sample(tf_m, .o("tf"), c_samp),
                          .o("tf"), c_samp, primary_is_samples = FALSE)
  # maxVAF: infinity (spine) WINS over panels(500>74)
  spine <- .coalesce_from(spine, .by_sample(m500_m, .o("maxvaf"), c_samp),
                          .o("maxvaf"), c_samp, TRUE)
  spine <- .coalesce_from(spine, .by_sample(m74_m, .o("maxvaf"), c_samp),
                          .o("maxvaf"), c_samp, TRUE)
  # meanVAF: infinity agg (spine) WINS over panels
  spine <- .coalesce_from(spine, .by_sample(m500_m, .o("meanvaf"), c_samp),
                          .o("meanvaf"), c_samp, TRUE)
  spine <- .coalesce_from(spine, .by_sample(m74_m, .o("meanvaf"), c_samp),
                          .o("meanvaf"), c_samp, TRUE)
  # MR score: tf_change > methylation > panels
  spine <- .coalesce_from(spine, .by_sample(tf_m, .o("mr"), c_samp),
                          .o("mr"), c_samp, FALSE)
  spine <- .coalesce_from(spine, .by_sample(methyl_m, .o("mr"), c_samp),
                          .o("mr"), c_samp, TRUE)
  spine <- .coalesce_from(spine, .by_sample(m500_m, .o("mr"), c_samp),
                          .o("mr"), c_samp, TRUE)
  # ctDNA % change: tf_change > panels
  spine <- .coalesce_from(spine, .by_sample(tf_m, "Ctdna_percentage_change", c_samp),
                          "Ctdna_percentage_change", c_samp, FALSE)

  # 5) variant-free / response-only sample retention: any Sample_ID seen in a
  #    melt but absent from the infinity spine becomes its own per-sample row.
  add_missing <- function(spine, m) {
    if (is.null(m) || !c_samp %in% names(m)) return(spine)
    new_ids <- setdiff(stats::na.omit(unique(m[[c_samp]])), spine[[c_samp]])
    if (!length(new_ids)) return(spine)
    add <- m[match(new_ids, m[[c_samp]]), , drop = FALSE]
    blank <- spine[rep(1, length(new_ids)), , drop = FALSE]; blank[] <- NA
    for (col in intersect(names(spine), names(add))) blank[[col]] <- add[[col]]
    blank[[c_samp]] <- new_ids
    rbind(spine, blank)
  }
  for (m in list(tf_m, m500_m, m74_m, methyl_m)) spine <- add_missing(spine, m)

  # 6) baseline relabel (single frame now; no cross-frame mismatch possible)
  if (c_time %in% names(spine)) {
    spine[["Visit_name_raw"]] <- spine[[c_time]]
    bl <- .ctdna_is_baseline_label(spine[[c_time]])
    spine[[c_time]][bl] <- .o("baseline")
    spine[["Baseline"]]  <- spine[[c_time]] == .o("baseline")
  }
  rownames(spine) <- NULL
  spine
}

# ---- clinical ---------------------------------------------------------------
# Accepts a single frame OR a named list (names give priority, e.g.
# list(adsl=, adrs=, adtr=)). Columns coalesce in list order (first wins).
.ctdna_collapse_subject <- function(f, key) {
  if (!key %in% names(f) || !anyDuplicated(f[[key]])) return(f)
  # one row per subject; first NON-NA value per column (handles multi-row adrs)
  parts <- split(f, f[[key]])
  out <- do.call(rbind, lapply(parts, function(g) {
    g[1, , drop = FALSE] -> r
    for (col in names(g)) {
      v <- g[[col]][!is.na(g[[col]])]
      if (length(v)) r[[col]] <- v[1]
    }
    r
  }))
  rownames(out) <- NULL
  out
}

.ctdna_build_clinical <- function(clinical, verbose = TRUE) {
  if (is.null(clinical)) return(NULL)
  if (is.data.frame(clinical))
    return(.ctdna_canonicalize_columns(clinical, "clinical", verbose = verbose))
  if (is.list(clinical)) {
    key <- .o("subject")
    frames <- lapply(seq_along(clinical), function(i)
      .ctdna_collapse_subject(
        .ctdna_canonicalize_columns(.ctdna_to_data_frame(clinical[[i]],
          names(clinical)[i] %||% paste0("clinical[[", i, "]]")),
          names(clinical)[i] %||% "clinical", verbose = verbose), key))
    base <- frames[[1]]
    if (length(frames) > 1) for (k in 2:length(frames)) {
      f <- frames[[k]]
      if (!key %in% names(base) || !key %in% names(f)) next
      idx <- match(base[[key]], f[[key]])
      # new columns from f, plus fill NAs in shared columns (clinical priority
      # already set by list order: earlier frame wins where non-NA)
      newc <- setdiff(names(f), names(base))
      for (col in newc) base[[col]] <- f[[col]][idx]
      shared <- intersect(setdiff(names(f), key), names(base))
      for (col in shared) base[[col]] <- .coalesce2(base[[col]], f[[col]][idx])
    }
    return(base)
  }
  .ctdna_canonicalize_columns(.ctdna_to_data_frame(clinical, "clinical"),
                              "clinical", verbose = verbose)
}

# ---- v0.51.0: ADTR (per-visit tumour size) ingestion -----------------------
# Canonicalise an ADTR frame (SUBJID->Patient_ID, AVISIT->Visit_name), filter to
# the target-lesion SUM rows, and keep Patient_ID/Visit_name/AVAL(SLD)/PCHG.
.ctdna_build_adtr <- function(adtr, verbose = TRUE) {
  if (is.null(adtr)) return(NULL)
  df <- .ctdna_canonicalize_columns(.ctdna_to_data_frame(adtr, "adtr"), "adtr",
                                    verbose = verbose)
  subj <- .o("subject"); tt <- .o("time")
  # honour the raw AVISIT for the tumour-assessment visit: if canonicalization did
  # not land it on Visit_name, take it verbatim from a raw visit column.
  if (!(tt %in% names(df))) {
    for (vc in c("AVISIT", "Tumor_visit_name", "VISIT", "AVISITN"))
      if (vc %in% names(df)) { df[[tt]] <- df[[vc]]; break }
  }
  if (!(subj %in% names(df))) {
    if (verbose) message("ctdna_prepare: adtr lacks a subject column after canonicalization; skipping.")
    return(NULL)
  }
  if (!(tt %in% names(df))) {
    if (verbose) message("ctdna_prepare: adtr lacks a visit (AVISIT) column after canonicalization; skipping.")
    return(NULL)
  }
  if (!("PCHG" %in% names(df))) {
    if (verbose) message("ctdna_prepare: adtr has no PCHG column; skipping.")
    return(NULL)
  }
  # v0.56.0 (user-mandated): keep EVERY adtr record that has a non-NA PCHG and a
  # non-empty AVISIT. NO PARCAT1 / sum-of-diameters / dedup filtering -- nothing
  # else is removed. PARAM / PARAMCD / PARCAT1 (and the raw AVISIT) are carried
  # through verbatim so the caller can subset later if desired.
  out <- data.frame(.subj = as.character(df[[subj]]),
                    .visit = as.character(df[[tt]]),   # raw AVISIT, verbatim
                    check.names = FALSE, stringsAsFactors = FALSE)
  names(out) <- c(subj, tt)
  out[[.o("lesion")]]     <- if ("AVAL" %in% names(df)) suppressWarnings(as.numeric(df$AVAL)) else NA_real_
  out[[.o("lesion_pct")]] <- suppressWarnings(as.numeric(df$PCHG))
  out$Tumor_param    <- if ("PARAM"   %in% names(df)) as.character(df$PARAM)   else NA_character_
  out$Tumor_paramcd  <- if ("PARAMCD" %in% names(df)) as.character(df$PARAMCD) else NA_character_
  out$Tumor_category <- if ("PARCAT1" %in% names(df)) as.character(df$PARCAT1) else NA_character_
  keep <- !is.na(out[[.o("lesion_pct")]]) &
          !is.na(out[[tt]]) & nzchar(trimws(out[[tt]]))
  out <- out[keep, , drop = FALSE]
  if (!nrow(out)) {
    if (verbose) message("ctdna_prepare: adtr has no rows with both a non-NA PCHG and an AVISIT; skipping.")
    return(NULL)
  }
  rownames(out) <- NULL
  out
}

# ---- v0.56.0: route an ADaM list by grain ------------------------------------
# Per-patient frames (adsl, collapsed adrs, any 1-row-per-subject frame) build
# $clinical; per-visit frames (adtr tumour size, labs, ...) build $assessments.
# Routing is by domain name first (adtr -> assessments; adsl/adrs -> clinical),
# then by grain for unknown domains (a visit column + >1 row per subject ->
# assessments). Returns a named clinical list (for .ctdna_build_clinical) and a
# list of tagged parts (for .ctdna_build_assessments).
.ctdna_route_adam <- function(adam, verbose = TRUE) {
  if (is.null(adam) || !length(adam)) return(list(clinical = NULL, assessments = NULL))
  if (is.data.frame(adam)) adam <- list(adsl = adam)
  nm <- names(adam); if (is.null(nm)) nm <- rep("", length(adam))
  nm[!nzchar(nm)] <- paste0("adam", which(!nzchar(nm)))
  subj <- .o("subject"); tt <- .o("time")
  clin <- list(); ass <- list(); log <- character()
  for (i in seq_along(adam)) {
    fr <- adam[[i]]; nmi <- tolower(nm[i]); if (is.null(fr)) next
    if (grepl("adtr", nmi)) {
      ass[[length(ass) + 1]] <- list(frame = fr, type = "tumor", name = nm[i])
      log <- c(log, sprintf("%s -> $assessments (tumour size)", nm[i])); next
    }
    if (grepl("adsl|adrs", nmi)) {
      clin[[nm[i]]] <- fr
      log <- c(log, sprintf("%s -> $clinical (per-patient)", nm[i])); next
    }
    f <- tryCatch(.ctdna_canonicalize_columns(.ctdna_to_data_frame(fr, nm[i]), nm[i],
                  verbose = FALSE), error = function(e) NULL)
    multivisit <- !is.null(f) && subj %in% names(f) && tt %in% names(f) &&
                  anyDuplicated(as.character(f[[subj]])) > 0
    if (multivisit) {
      ass[[length(ass) + 1]] <- list(frame = fr, type = "generic", name = nm[i])
      log <- c(log, sprintf("%s -> $assessments (per-visit)", nm[i]))
    } else {
      clin[[nm[i]]] <- fr
      log <- c(log, sprintf("%s -> $clinical (per-patient)", nm[i]))
    }
  }
  if (verbose && length(log))
    message("ctdna_prepare: ADaM routing:\n  ", paste(log, collapse = "\n  "))
  list(clinical = if (length(clin)) clin else NULL,
       assessments = if (length(ass)) ass else NULL)
}

# ---- v0.56.0: build the long $assessments frame ------------------------------
# Patient_ID x <visit> x parameter. The tumour-assessment visit is kept in its OWN
# column `Tumor_visit_name` (raw AVISIT, verbatim -- never normalised and never
# merged into the ctDNA `Visit_name`), because tumour imaging and ctDNA draws are
# on different schedules. Tumour rows carry the opts-named Lesion_sum_mm +
# Tumor_size_pct_change so the per-visit tumour token resolves with no special-
# casing; generic per-visit frames melt numeric measures into Value under their
# own `Visit_name`. ALWAYS returns a frame (empty, with the full schema, if there
# are no usable rows) so $assessments is always present.
.ctdna_empty_assessments <- function() {
  subj <- .o("subject"); les <- .o("lesion"); lp <- .o("lesion_pct")
  out <- data.frame(a = character(0), Tumor_visit_name = character(0),
    Visit_name = character(0), Param_code = character(0), Param = character(0),
    Category = character(0), Unit = character(0), Source = character(0),
    Value = numeric(0), b = numeric(0), c = numeric(0),
    check.names = FALSE, stringsAsFactors = FALSE)
  names(out) <- c(subj, "Tumor_visit_name", "Visit_name", "Param_code", "Param",
                  "Category", "Unit", "Source", "Value", les, lp)
  out
}

.ctdna_build_assessments <- function(parts, verbose = TRUE) {
  if (is.null(parts) || !length(parts)) return(.ctdna_empty_assessments())
  subj <- .o("subject"); tt <- .o("time"); les <- .o("lesion"); lp <- .o("lesion_pct")
  rows <- list()
  for (p in parts) {
    if (identical(p$type, "tumor")) {
      af <- .ctdna_build_adtr(p$frame, verbose = verbose)
      if (is.null(af) || !nrow(af)) {
        warning("ctdna_prepare: adtr produced no usable tumour rows; $assessments ",
                "will be created but empty for the tumour part. Check PARCAT1 (Target), ",
                "a sum-of-diameters PARAM, and a non-NA PCHG.", call. = FALSE); next
      }
      r <- data.frame(check.names = FALSE, stringsAsFactors = FALSE,
        a = as.character(af[[subj]]),
        Tumor_visit_name = as.character(af[[tt]]),   # raw AVISIT, verbatim
        Param_code = if ("Tumor_paramcd" %in% names(af) && !all(is.na(af$Tumor_paramcd)))
                       as.character(af$Tumor_paramcd) else "SUMDIAM",
        Param = if ("Tumor_param" %in% names(af)) as.character(af$Tumor_param)
                else "Sum of target lesion diameters",
        Category = if ("Tumor_category" %in% names(af)) as.character(af$Tumor_category) else NA_character_,
        Unit = "mm", Source = "adtr", Value = NA_real_)
      names(r)[1] <- subj
      r[[les]] <- if (les %in% names(af)) af[[les]] else NA_real_
      r[[lp]]  <- if (lp  %in% names(af)) af[[lp]]  else NA_real_
      rows[[length(rows) + 1]] <- r
    } else {
      f <- tryCatch(.ctdna_canonicalize_columns(.ctdna_to_data_frame(p$frame, p$name),
                    p$name, verbose = FALSE), error = function(e) NULL)
      if (is.null(f) || !all(c(subj, tt) %in% names(f))) next
      num <- setdiff(names(f)[vapply(f, is.numeric, logical(1))], c(subj, tt))
      for (cc in num) {
        r <- data.frame(check.names = FALSE, stringsAsFactors = FALSE,
          a = as.character(f[[subj]]), Visit_name = as.character(f[[tt]]),
          Param_code = cc, Param = cc, Category = NA_character_,
          Unit = NA_character_, Source = p$name,
          Value = suppressWarnings(as.numeric(f[[cc]])))
        names(r)[1] <- subj
        r[[les]] <- NA_real_; r[[lp]] <- NA_real_
        rows[[length(rows) + 1]] <- r
      }
    }
  }
  if (!length(rows)) return(.ctdna_empty_assessments())
  allc <- Reduce(union, lapply(rows, names))
  rows <- lapply(rows, function(d) {
    for (m in setdiff(allc, names(d)))
      d[[m]] <- if (m %in% c(les, lp, "Value")) NA_real_ else NA_character_
    d[, allc, drop = FALSE] })
  out <- do.call(rbind, rows); rownames(out) <- NULL
  out
}

# Normalise a visit label for matching: "Cycle 2 Day 1"/"C02D01"/"cycle2day1" -> "C2D1".
.ctdna_normalize_visit <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("CYCLE\\s*0*([0-9]+)\\s*DAY\\s*0*([0-9]+)", "C\\1D\\2", x)
  x <- gsub("C0+([0-9]+)D0+([0-9]+)", "C\\1D\\2", x)
  x <- gsub("[[:space:]_\\-]+", "", x)
  x
}

# Join per-visit ADTR columns onto samples by Patient_ID x Visit_name, matching
# visit labels tolerantly (format-insensitive) so "Cycle 2 Day 1" == "C2D1".
.ctdna_join_adtr <- function(samples, adtr) {
  if (is.null(adtr) || is.null(samples)) return(samples)
  subj <- .o("subject"); tt <- .o("time")
  if (!all(c(subj, tt) %in% names(samples))) return(samples)
  idx <- match(paste(as.character(samples[[subj]]), .ctdna_normalize_visit(samples[[tt]])),
               paste(as.character(adtr[[subj]]),    .ctdna_normalize_visit(adtr[[tt]])))
  for (col in setdiff(names(adtr), c(subj, tt)))
    samples[[col]] <- adtr[[col]][idx]
  samples
}

# join clinical Dose/RECIST/Cancertype onto samples (clinical WINS)
.ctdna_join_clinical <- function(samples, clinical) {
  if (is.null(clinical) || is.null(samples)) return(samples)
  key <- .o("subject")
  if (!key %in% names(clinical) || !key %in% names(samples)) return(samples)
  idx <- match(samples[[key]], clinical[[key]])
  for (col in c(.o("dose"), .o("recist"), .o("cancertype"), .o("best_change"))) {
    if (col %in% names(clinical)) {
      inc <- clinical[[col]][idx]
      samples[[col]] <- .coalesce2(inc, samples[[col]])   # clinical wins
    }
  }
  samples
}

# ---- sample QC (per-sample; cascades to variants in the caller) -------------
.ctdna_sample_pass <- function(samples,
                               sample_status = .o("sample_qc_status"),
                               plasma_min    = .o("sample_qc_plasma_min"),
                               cfdna_min     = .o("sample_qc_cfdna_min")) {
  n <- nrow(samples); pass <- rep(TRUE, n)
  st <- .o("col_sample_status")
  if (!is.null(sample_status) && length(sample_status) && st %in% names(samples))
    pass <- pass & (tolower(as.character(samples[[st]])) %in% tolower(sample_status) |
                      is.na(samples[[st]]))
  pcol <- "Plasma_ml_input"
  if (!is.null(plasma_min) && plasma_min > 0 && pcol %in% names(samples))
    pass <- pass & (suppressWarnings(as.numeric(samples[[pcol]])) >= plasma_min |
                      is.na(samples[[pcol]]))
  ccol <- .o("cfdna_conc")
  if (!is.null(cfdna_min) && cfdna_min > 0 && ccol %in% names(samples))
    pass <- pass & (suppressWarnings(as.numeric(samples[[ccol]])) >= cfdna_min |
                      is.na(samples[[ccol]]))
  pass
}

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a


# ---- dictionary -------------------------------------------------------------
.ctdna_build_dictionary <- function(prep) {
  reg <- .ctdna_registry()
  present <- function(col, frames) {
    fr <- strsplit(frames, ",")[[1]]
    any(vapply(fr, function(f)
      is.data.frame(prep[[f]]) && col %in% names(prep[[f]]), logical(1)))
  }
  keep <- vapply(seq_len(nrow(reg)), function(i)
    is.na(reg$column[i]) || present(reg$column[i], reg$frame[i]), logical(1))
  reg <- reg[keep, , drop = FALSE]
  rownames(reg) <- NULL
  reg
}

# =============================================================================
# ctdna_prepare()
# =============================================================================

#' Prepare a study bundle into canonical ctDNA frames
#'
#' Gathers GH Infinity ctDNA modalities (plus optional clinical / IHC /
#' DNA-seq / RNA-seq), canonicalizes every column to its \code{ctdna_opts()}
#' name, and returns a classed \code{ctdna_prep} list with a per-sample frame
#' (\code{$samples}), a per-variant frame (\code{$variants}), the clinical
#' frame (\code{$clinical}), a per-visit non-variant \code{$assessments} frame,
#' and a data \code{$dictionary}. \code{infinity_report} is required.
#'
#' All inputs are routed to the right frame automatically (\dQuote{under the
#' hood}); the caller never names the internal frames. ADaM datasets go in one
#' \code{adam = list(adsl = , adrs = , adtr = , ...)} argument and are routed by
#' grain: one-row-per-patient frames (adsl, best-response adrs) build
#' \code{$clinical}; per-visit frames (adtr tumour size, labs, ...) build
#' \code{$assessments}.
#'
#' @param infinity_report Required. The GH Infinity consolidated report.
#' @param `74_genomic_panel`,`500_genomic_panel` Optional GH 74-/500-gene panel
#'   response files (paired A/B wide format). Names start with a digit, so in a
#'   call they need backticks, e.g. \code{`74_genomic_panel` = df}.
#' @param tf_change Optional GH tumour-fraction change file.
#' @param Methylation_GuardantResponseRUOReport Optional GH methylation response
#'   (MR) RUO report.
#' @param adam Optional named list of ADaM datasets, e.g.
#'   \code{list(adsl = , adrs = , adtr = )}. Any ADaM domain may be supplied;
#'   each is routed by grain (per-patient -> \code{$clinical}; per-visit ->
#'   \code{$assessments}). ADTR drives \code{tumor_size_change} (per visit, in
#'   \code{$assessments}) and \code{best_tumor_size_change} (nadir, in
#'   \code{$clinical}).
#' @param verbose Logical. When \code{TRUE}, prints how each \code{adam} frame
#'   was routed.
#' @param ... Reserved (back-compatible aliases: \code{panel_74_response},
#'   \code{panel_500_response}, \code{methylation_response}, \code{clinical},
#'   \code{adtr}, \code{ihc}, \code{dnaseq}, \code{rnaseq}, \code{wes},
#'   \code{bundle}).
#' @return An object of class \code{ctdna_prep}.
#' @export
ctdna_prepare <- function(infinity_report                       = NULL,
                          `74_genomic_panel`                    = NULL,
                          `500_genomic_panel`                   = NULL,
                          tf_change                             = NULL,
                          Methylation_GuardantResponseRUOReport = NULL,
                          adam                                  = NULL,
                          verbose                               = TRUE,
                          ...) {

  dots <- list(...)
  `%or%` <- function(a, b) if (!is.null(a)) a else b
  # ---- ctDNA frames (primary names + back-compatible aliases via ...) ----
  args <- list(
    infinity_report      = infinity_report      %or% dots$infinity_report,
    tf_change            = tf_change            %or% dots$tf_change,
    panel_74_response    = `74_genomic_panel`   %or% dots$panel_74_response  %or% dots$panel_74,
    panel_500_response   = `500_genomic_panel`  %or% dots$panel_500_response %or% dots$panel_500,
    methylation_response = Methylation_GuardantResponseRUOReport %or% dots$methylation_response,
    ihc = dots$ihc, dnaseq = dots$dnaseq, rnaseq = dots$rnaseq, wes = dots$wes)
  args <- args[!vapply(args, is.null, logical(1))]
  bundle <- dots$bundle
  if (!is.null(bundle)) {
    if (!is.list(bundle) || is.null(names(bundle)) || any(!nzchar(names(bundle))))
      stop("ctdna_prepare: `bundle` must be a NAMED list of frames.", call. = FALSE)
    args[setdiff(names(bundle), names(args))] <- bundle[setdiff(names(bundle), names(args))]
  }

  # ---- assemble the ADaM input list (primary `adam=` + back-compat clinical/adtr) ----
  adam_in <- if (is.null(adam)) list()
             else if (is.data.frame(adam)) list(adsl = adam)
             else as.list(adam)
  if (!is.null(dots$clinical)) {
    cl <- dots$clinical
    if (is.list(cl) && !is.data.frame(cl)) adam_in <- c(adam_in, cl)
    else adam_in[["adsl"]] <- adam_in[["adsl"]] %or% cl
  }
  if (!is.null(dots$adtr)) adam_in[["adtr"]] <- adam_in[["adtr"]] %or% dots$adtr
  routed <- .ctdna_route_adam(adam_in, verbose = verbose)

  if (is.null(args$infinity_report))
    stop("ctdna_prepare: `infinity_report` is required (the GH Infinity report). ",
         "Pass at least infinity_report = <df>.", call. = FALSE)


  # ---- coerce + canonicalize each ctDNA frame ----
  canon <- function(x, nm) if (is.null(x)) NULL else
    .ctdna_canonicalize_columns(.ctdna_to_data_frame(x, nm), nm, verbose = verbose)
  inf    <- canon(args$infinity_report,      "infinity_report")
  tf     <- if (!is.null(args$tf_change))            .ctdna_to_data_frame(args$tf_change, "tf_change") else NULL
  m500   <- if (!is.null(args$panel_500_response))   .ctdna_to_data_frame(args$panel_500_response, "panel_500_response") else NULL
  m74    <- if (!is.null(args$panel_74_response))    .ctdna_to_data_frame(args$panel_74_response, "panel_74_response") else NULL
  methyl <- if (!is.null(args$methylation_response)) .ctdna_to_data_frame(args$methylation_response, "methylation_response") else NULL

  if (verbose) message("ctdna_prepare: building per-sample and per-variant frames")
  samples  <- .ctdna_build_samples(inf, tf, m500, m74, methyl, verbose = verbose)
  variants <- .ctdna_build_variants(inf)
  clinical_out <- .ctdna_build_clinical(routed$clinical, verbose = verbose)
  samples  <- .ctdna_join_clinical(samples, clinical_out)

  # ---- v0.56.0: per-visit, non-variant measurements (ADTR tumour size, labs, ...)
  #      live in their OWN long frame ($assessments) at Patient_ID x Visit x param
  #      grain. $samples stays PURE ctDNA (one row per ctDNA sample). The per-patient
  #      tumour nadir/best is derived here and stored in $clinical.
  subj <- .o("subject"); tt <- .o("time"); lp <- .o("lesion_pct"); les <- .o("lesion"); bc <- .o("best_change")
  assessments <- .ctdna_build_assessments(routed$assessments, verbose = verbose)
  if (!is.null(assessments) && lp %in% names(assessments)) {
    # per-patient nadir of tumour PCHG (independent of any ctDNA visit alignment)
    tum <- assessments[is.finite(assessments[[lp]]), , drop = FALSE]
    if (nrow(tum)) {
      nad <- tapply(tum[[lp]], as.character(tum[[subj]]),
                    function(z) { z <- z[is.finite(z)]; if (length(z)) min(z) else NA_real_ })
      if (!is.null(clinical_out) && subj %in% names(clinical_out)) {
        clinical_out[["best_lesion_change"]] <- as.numeric(nad[as.character(clinical_out[[subj]])])
        if (!(bc %in% names(clinical_out)) || all(is.na(clinical_out[[bc]])))
          clinical_out[[bc]] <- as.numeric(nad[as.character(clinical_out[[subj]])])
        # refresh the clinical columns broadcast onto $samples
        samples <- .ctdna_join_clinical(samples, clinical_out)
      }
    }
  }

  out <- list(samples = samples, variants = variants)
  if (!is.null(clinical_out))  out$clinical    <- clinical_out
  if (!is.null(assessments))   out$assessments <- assessments

  # v0.63.0: cohort / indication grouping across the package (oncoprint,
  # alteration grid, faceting) keys off the canonical `Indication` column.
  # Tell the user to create it if it is absent. It is derived from an adrs
  # `indication` column (or any alias the canonicaliser folds into `Indication`,
  # e.g. Cancer_type / disease / TumorType).
  if (isTRUE(verbose) &&
      (is.null(clinical_out) || !"Indication" %in% names(clinical_out)))
    message("ctdna_prepare: no `Indication` column in the clinical data. ",
            "Oncoprint/grid cohort and indication grouping rely on it -- add an ",
            "`indication` column to your adrs (e.g. Nsq_NSCLC / HNSCC / ES_SCLC / ",
            "mCRPC) so it is carried through as canonical `Indication`.")

  # ---- extra modalities (reuse v0.43 modality preppers) ----
  c_subj <- .o("subject"); c_gene <- .o("gene"); c_expr <- .o("expression")
  qc_state <- c(samples = FALSE, variants = FALSE)
  if (!is.null(args$ihc))    { out$ihc    <- .prep_ihc(.ctdna_to_data_frame(args$ihc, "ihc"), c_subj); qc_state["ihc"] <- FALSE }
  dnaseq_in <- args$dnaseq %||% args$wes
  if (!is.null(dnaseq_in))   { out$dnaseq <- .prep_dnaseq(.ctdna_to_data_frame(dnaseq_in, "dnaseq"), c_subj); qc_state["dnaseq"] <- FALSE }
  if (!is.null(args$rnaseq)) { out$rnaseq <- .prep_expression(args$rnaseq, c_subj, c_gene, c_expr); qc_state["rnaseq"] <- FALSE }

  out <- structure(out, class = "ctdna_prep",
                   platform = "guardant_health_infinity",
                   qc_state = qc_state)

  # record pre-(sample QC) counts so print can show before -> after
  attr(out, "n_pre_qc") <- c(
    samples  = nrow(out$samples),
    variants = if (is.data.frame(out$variants)) nrow(out$variants) else NA_integer_)

  # ---- sample QC removed from prepare in v0.75.0 -- Sample_QC arg gone ----
  # Users call ctdna_sample_qc() explicitly:
  #     prep <- ctdna_prepare(...)
  #     prep <- ctdna_sample_qc(prep)
  if (isTRUE(verbose))
    message("ctdna_prepare: prep$variants is raw (sample QC not applied). ",
            "Call ctdna_sample_qc(prep) before analysis.")

  out$dictionary <- .ctdna_build_dictionary(out)
  out
}

# =============================================================================
# ctdna_sample_qc() -- re-runnable per-sample QC over the whole prep, with
# cascade to variants. Call with just the prep object.
# =============================================================================

#' Apply per-sample QC across a prepared bundle (cascades to variants)
#'
#' Runs per-sample QC on \code{prep$samples} using \code{ctdna_opts()}
#' thresholds (override per call), drops failing samples, and removes every
#' variant of a dropped sample from \code{prep$variants} (cascade keyed on
#' \code{Sample_ID}). Re-runnable after adding data.
#'
#' @param prep A \code{ctdna_prep} object.
#' @param sample_status,plasma_min,cfdna_min QC thresholds; default to
#'   \code{ctdna_opts()}.
#' @param verbose Logical.
#' @return The filtered \code{ctdna_prep}.
#' @export
ctdna_sample_qc <- function(prep,
                            sample_status = .o("sample_qc_status"),
                            plasma_min    = .o("sample_qc_plasma_min"),
                            cfdna_min     = .o("sample_qc_cfdna_min"),
                            verbose       = TRUE) {
  if (!inherits(prep, "ctdna_prep"))
    stop("ctdna_sample_qc: `prep` must be a ctdna_prep object from ctdna_prepare().",
         call. = FALSE)
  c_samp <- .o("sample"); c_subj <- .o("subject")
  s <- prep$samples
  pass <- .ctdna_sample_pass(s, sample_status, plasma_min, cfdna_min)
  kept_ids  <- s[[c_samp]][pass]
  removed_samples  <- s[!pass, , drop = FALSE]
  prep$samples <- s[pass, , drop = FALSE]

  removed_variants <- prep$variants[0, , drop = FALSE]
  if (is.data.frame(prep$variants) && c_samp %in% names(prep$variants)) {
    drop_v <- !(prep$variants[[c_samp]] %in% kept_ids)
    removed_variants <- prep$variants[drop_v, , drop = FALSE]
    prep$variants    <- prep$variants[!drop_v, , drop = FALSE]
  }

  # clinical cascade: a patient with no surviving sample is dropped from
  # $clinical and recorded.
  removed_clinical <- NULL
  if (is.data.frame(prep$clinical) && c_subj %in% names(prep$clinical)) {
    surviving_pts <- unique(prep$samples[[c_subj]])
    drop_c <- !(prep$clinical[[c_subj]] %in% surviving_pts)
    removed_clinical <- prep$clinical[drop_c, , drop = FALSE]
    prep$clinical    <- prep$clinical[!drop_c, , drop = FALSE]
    rownames(prep$clinical) <- NULL
  }

  rownames(prep$samples) <- NULL
  if (is.data.frame(prep$variants)) rownames(prep$variants) <- NULL
  rownames(removed_samples) <- NULL; rownames(removed_variants) <- NULL

  # accumulate across repeated QC runs
  prior <- prep$qc_removed
  prep$qc_removed <- list(
    samples  = if (is.null(prior)) removed_samples  else rbind(prior$samples,  removed_samples),
    variants = if (is.null(prior)) removed_variants else rbind(prior$variants, removed_variants),
    clinical = if (is.null(prior)) removed_clinical else rbind(prior$clinical, removed_clinical))

  qs <- attr(prep, "qc_state"); qs[c("samples","variants")] <- TRUE
  attr(prep, "qc_state") <- qs
  if (verbose) {
    message(sprintf("ctdna_sample_qc: %d of %d samples pass (%.1f%%); %d variants and %d clinical record(s) cascaded out.",
                    sum(pass), length(pass), 100 * mean(pass),
                    nrow(removed_variants), if (is.null(removed_clinical)) 0L else nrow(removed_clinical)))
    if (any(!pass))
      message("  removed sample IDs: ",
              paste(utils::head(removed_samples[[c_samp]], 20), collapse = ", "),
              if (sum(!pass) > 20) sprintf(" ... (+%d more)", sum(!pass) - 20) else "")
  }
  prep$dictionary <- .ctdna_build_dictionary(prep)
  prep
}

# =============================================================================
# ctdna_prep_add() -- add a modality to an existing prep; flags it un-QC'd.
# =============================================================================

#' Add a modality to an existing prep (flags it as needing QC)
#'
#' @param prep A \code{ctdna_prep} object.
#' @param ... One named modality, e.g. \code{rnaseq = mat} or
#'   \code{infinity_report = df} (new ctDNA re-merges into samples/variants).
#' @param verbose Logical.
#' @return The updated \code{ctdna_prep}; prints a notice to run
#'   \code{ctdna_sample_qc(prep)}.
#' @export
ctdna_prep_add <- function(prep, ..., verbose = TRUE) {
  if (!inherits(prep, "ctdna_prep"))
    stop("ctdna_prep_add: `prep` must be a ctdna_prep object.", call. = FALSE)
  new <- list(...)
  if (!length(new) || is.null(names(new)) || any(!nzchar(names(new))))
    stop("ctdna_prep_add: pass exactly one named modality, e.g. rnaseq = mat.",
         call. = FALSE)
  c_subj <- .o("subject"); c_gene <- .o("gene"); c_expr <- .o("expression")
  qs <- attr(prep, "qc_state")
  for (nm in names(new)) {
    x <- new[[nm]]
    if (nm == "ihc")          prep$ihc    <- .prep_ihc(.ctdna_to_data_frame(x, nm), c_subj)
    else if (nm == "dnaseq" || nm == "wes") prep$dnaseq <- .prep_dnaseq(.ctdna_to_data_frame(x, nm), c_subj)
    else if (nm == "rnaseq")  prep$rnaseq <- .prep_expression(x, c_subj, c_gene, c_expr)
    else prep[[nm]] <- .ctdna_to_data_frame(x, nm)
    qs[nm] <- FALSE
    message(sprintf("ctdna_prep_add: '%s' added but NOT yet QC'd. Run ctdna_sample_qc(prep) before analysis.", nm))
  }
  attr(prep, "qc_state") <- qs
  prep$dictionary <- .ctdna_build_dictionary(prep)
  prep
}

# ---- grain accessor: functions accept whole `prep` and pull their grain -----
.ctdna_grain <- function(x, grain = c("samples", "variants")) {
  grain <- match.arg(grain)
  if (inherits(x, "ctdna_prep")) return(x[[grain]])
  if (is.data.frame(x)) return(x)   # already a frame (back-compat path)
  stop("expected a ctdna_prep object or a data frame.", call. = FALSE)
}

# ---- print method -----------------------------------------------------------
#' @export
print.ctdna_prep <- function(x, ...) {
  qs  <- attr(x, "qc_state")
  pre <- attr(x, "n_pre_qc")
  cat(sprintf("<ctdna_prep>  platform: %s\n", attr(x, "platform") %||% "?"))
  frames <- names(x)[vapply(x, is.data.frame, logical(1))]
  for (f in frames) {
    d <- x[[f]]
    if (f %in% c("samples", "variants") && isTRUE(qs[[f]]) &&
        !is.null(pre) && f %in% names(pre) && !is.na(pre[[f]])) {
      # before -> after sample-based QC
      cat(sprintf("  $%-10s %6d -> %-6d x %-3d  Sample Based QC'd\n",
                  f, pre[[f]], nrow(d), ncol(d)))
    } else {
      flag <- if (f %in% names(qs))
                (if (isTRUE(qs[[f]])) "Sample Based QC'd" else "needs Sample QC") else ""
      cat(sprintf("  $%-10s %6d x %-3d  %s\n", f, nrow(d), ncol(d), flag))
    }
  }
  if (!is.null(x$dictionary))
    cat(sprintf("  %d variables in $dictionary (variable / column / frame / source / units / purpose)\n",
                nrow(x$dictionary)))
  invisible(x)
}
