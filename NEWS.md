# ctdnaTM 0.75.9

* **Removed the legacy `ctdna_make_patient_data()` API entirely.** `R/clinical_adam.R` (679 lines — the function + its private helpers, options infrastructure, and public `clinical_opts()`) is gone. Nothing in the current pipeline called it; users on the modern `ctdna_prepare(adam = list(...))` API are unaffected. If you were still using it, switch to `ctdna_prepare()` — `prep$clinical` now provides the same per-patient frame.
* **Removed all hardcoded proprietary references** (GSK, trial identifiers, protocol numbers, cohort strings, treatment-arm dose labels) from package sources, docs, NEWS, README, and the shipped deliverable/starter templates. Doc-comment examples that used real compound names have been rewritten to be generic.
* **Bundled logo restored.** `inst/extdata/logo.png` now ships with the package again, so `ctdna_logo()` returns a valid path after installation. A second copy lives at repo-root `images/logo.png` for the GitHub page (excluded from the tarball via `.Rbuildignore`).
* **Logo at top of README and core deliverables Rmd.** The core deliverable now uses `knitr::include_graphics(ctdnaTM::ctdna_logo())` in its first chunk.


* **`prep$clinical` is now the UNION of every ADaM frame you pass, not just the first one.** Previously `.ctdna_build_clinical()` left-joined subsequent frames into `frames[[1]]` — any patient present in `adrs` or `adtr` but not in `adsl` was silently dropped. Fixed to compute the union of `Patient_ID` across every non-empty frame, then coalesce columns in list order (first non-NA wins).
* **`adtr` frames now contribute their patient IDs to `prep$clinical`** in addition to `$assessments`. Per the user directive that the oncoprint should include every patient the user provides, `adtr`-only patients now appear in the landscape cohort.
* Empty frames (0 rows after canonicalization) are skipped silently.


* **Landscape plots — Dose is now optional.** `.build_landscape_base()` previously dropped every patient missing Indication OR Dose. Now only Indication is required (needed for the wrap-facet key); Dose becomes an optional annotation column that can be `NA` without dropping the patient. RECIST was already optional.


* **Landscape plots — `visit = "baseline"` now matches every equivalent label.** `.resolve_visit_baseline()` used to collapse `"baseline"` down to just `"C1D1"`; the downstream normaliser then dropped any variant whose `Visit_name` was `"Baseline"`, `"Screening"`, or another baseline-equivalent form even though `.ctdna_is_baseline_label()` elsewhere in the package recognises them. Now the resolver expands `"baseline"` into the full equivalent set (`"C1D1"`, `"Cycle 1 Day 1"`, `"Baseline"`, `"Screening"`) so downstream matching picks up all four forms. Fixes the "N patients missing from oncoprint" bug seen when the Guardant report labels baseline samples as `"Baseline"` or `"Screening"` rather than `"C1D1"`.


* **`ctdna_prepare()` no longer filters variants.** `prep$variants` is now the `infinity_report` verbatim (with only the `Baseline` flag column added). The previous "drop rows where Gene is empty" filter is gone — nothing is removed at prep time.
* **`ctdna_sample_qc()` cascade tightened.**
  * Variants are removed only when their `Sample_ID` matches a sample that failed QC. Variants with `NA` or unmatched `Sample_ID` are kept.
  * The clinical cascade is removed. `prep$clinical` is no longer touched by sample QC — patients whose samples all failed QC remain in `$clinical` (they show as all-wildtype columns in the oncoprint) instead of being moved to `qc_removed$clinical`.
* **Landscape base reverted** to its v0.75.3 behavior (v0.75.4's `qc_removed$clinical` reunite was replaced by fixing the source).

# ctdnaTM 0.75.4

* Landscape plots widened the cohort by reuniting `prep$qc_removed$clinical` (superseded in 0.75.5 by fixing the source: `ctdna_prepare` no longer filters and sample QC no longer cascades to clinical).


* **Landscape plots (`ctdna_oncoprint`, `ctdna_alteration_grid`) — widen the cohort.** `.build_landscape_base()` now:
  * Reunites `prep$qc_removed$clinical` rows into the local clinical frame for patients whose variants survived sample QC. Patients archived for non-sample-QC reasons (cohort filter, ADSL mismatch, etc.) no longer disappear from the landscape plots.
  * Requires only `Indication` (the wrap-facet key); `Dose` and `RECIST` are now optional annotations and can be `NA` without dropping the patient.
* Reports how many patients were restored: `"restored N patient(s) from prep$qc_removed$clinical (variants survived sample QC)."`
* Trajectory plots (`baseline`, `reduction`, `longitudinal`, `boxplot`, `scatter`) are unchanged — they still work off the strict `$samples`/`$variants` cohort where sample+variant pairing matters.


* **Pipeline registry corrections** — cross-modality steps re-pointed at the frames that actually hold the data:
  * `vs_tmb` needs `$samples` only (TMB values are a column on `prep$samples`).
  * `vs_mutation` needs `$samples` **and** `$variants` (mutation status comes from variants, not from an external `$genomic_74` frame).
  * `vs_ihc` and `vs_expression` still need `$samples` plus the user-attached `$ihc` / `$expression` frames.
* Every version now ships the full Rmd/md set on every iteration: `README.md`, `NEWS.md`, three vignettes (`walkthrough`, `function_tour`, `filter_internals`), the core deliverables Rmd, and the analysis starter.

# ctdnaTM 0.75.2

* **Pipeline registry rewrite** — `prep$ctdna` is gone; every step now checks the actual frame it consumes.
* 5 ghost registry entries removed (pointed at deleted functions: `baseline_dose_recist`, `mr_recist`, `longitudinal_recist`, `reduction_vs_tumor`, `mut_methyl`).
* 11 remaining entries rewritten. Sample-level plots use `needs = "samples"`. Cross-modality entries name the actual frame they need.
* Skip messages now name the actual frame: `"needs prep$ihc (not found in prep)"` etc.
* `ctdna_pipeline()` `@param prep` doc rewritten to list current prep frames.
* All 13 `prep$ctdna` references in doc-string examples swept to `prep$samples`.

# ctdnaTM 0.75.1

* Documentation refresh across all Rmd/md files (except NEWS).
* README updated with all four biomarker schemes + branch-label note.
* Function-tour vignette refreshed: `TSG_Tier3` added, branch labels explained with example output.
* NEW vignette `ctdnaTM_walkthrough.Rmd` — one-page end-to-end.
* Fresh analysis starter for real data.

# ctdnaTM 0.75.0

* **Dead-code sweep**: 269 lines of orphan roxygen removed; dead helpers deleted; sentinel removed; deprecated aliases removed; `Sample_QC` arg removed from `ctdna_prepare`.
* **`TSG_Tier3()`** added — most restrictive TSG scheme (TP53/RB1/PTEN). No LGR, no fusion, no rare-syn carve-in.
* **Branch labels** — `create_filtering_scheme()` accepts `branches = list(name = expr, ...)`. `filtering_criteria` column now shows `PASS::TSG_Tier1(somatic_nonsyn)`.
* **New vignette** `ctdnaTM_filter_internals.Rmd` — object shapes, evaluator, engine internals, extension recipes.
* **Runnable `@examples`** on all four biomarker schemes.


## Documentation release

* Comprehensive **function-tour vignette** (`vignette("ctdnaTM_function_tour")`) covering every user-facing function with runnable examples on the mock study.
* New **README.md** with a five-minute tour, filter surface summary, and grammar quick-reference.
* This **NEWS.md**.
* `@param` coverage completed on the remaining user-facing functions (`ctdna_pipeline`).
* All 151 exports now have help pages (autocomplete + `?fn` tooltips populated).
* No API changes vs 0.73.1.

# ctdnaTM 0.73.1

* **CHIP mutations excluded from every built-in scheme** (`TSG_Tier1`, `TSG_Tier2`, `HRR14` all now `not(rule_ch())` in their gates).
* `TSG_Tier2` body reverted to the old `scheme_TSG` structure with the LGR path swapped to the new intra-chromosomal rule.
* `HRR14` body reverted to the old `scheme_HRR` structure with the LGR path swapped to the new intra-genic rule.

# ctdnaTM 0.73.0

* Biomarker schemes redone.
* Old `scheme_TSG` and `scheme_HRR` **removed**.
* Three new built-ins:
  * `TSG_Tier1()` — most permissive (TP53/RB1/PTEN).
  * `TSG_Tier2()` — intermediate (TP53/RB1/PTEN).
  * `HRR14()` — HRR biomarker on 14 HRR genes.

# ctdnaTM 0.72.0

* Filtering API consolidated to `ctdna_variant_filter(prep, filter_scheme, apply, explain)` — both `apply` and `explain` mandatory; passing both `FALSE` errors.
* `ctdna_qc_filter` renamed to `ctdna_sample_qc` (affects `prep$samples` and cascades to `prep$variants`).
* Removed `ctdna_apply_filter` and `ctdna_explain_filter`.
* Sample QC removed from `ctdna_variant_filter`; call `ctdna_sample_qc` explicitly beforehand.

# ctdnaTM 0.71.0

* Filtering makeover: prep-in / prep-out for all filter functions.
* `filtering_criteria` column stored on the annotated snapshot only, not on live `prep$variants`.
* Sample QC removed from `ctdna_prepare()`; `Sample_QC` default flipped to `FALSE`.
* Old `ctdna_filter_apply(df, ...)` removed from exports.

# ctdnaTM 0.70.0

* `ctdna_pipeline()` redesign: results accessible as `res$<key>$plot`; every top-level argument flows to steps that accept it via an `opts` bag.

# ctdnaTM 0.69.0 – 0.66.0

* Clean rewrite of the two genomic-landscape plots.
* `ctdna_oncoprint(prep, gene_sets, ...)` and `ctdna_alteration_grid(prep, gene_sets, ...)`.
* `alt_any()`, `alt_all()`, `alt_any_n(n)` combine-mode wrappers for gene-set entries.
