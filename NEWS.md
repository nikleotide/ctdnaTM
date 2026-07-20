# ctdnaTM 0.76.1

## Fix — gene names now render on `ctdna_concordance_oncoprint_core()`

* On real data (longer gene names like `SMARCA4`, `RAD51B`, `BARD1`,
  `TERT_promoter`), the gene-column labels at the bottom of the last
  panel were being visually covered by the gene-set (TSG / HRR) annotation
  strip. Root cause: ComplexHeatmap's built-in `show_column_names = TRUE`
  on a Heatmap that is followed by another `HeatmapAnnotation` in a
  `%v%` stack has a layout race that even `auto_adjust = FALSE` and
  extra `heatmap_height` slack don't resolve.
* v0.76.1 draws gene names as an explicit `HeatmapAnnotation` stack
  element (using `anno_text`) inserted between the last heatmap panel
  and the gene-set strip. Both now coexist with a small visual gap and
  neither can cover the other regardless of gene-name length.
* New parameter: **`show_gene_set_annotation`** (default `TRUE`) — set
  to `FALSE` to suppress the bottom TSG/HRR strip entirely without
  losing gene names.

# ctdnaTM 0.76.0

## New — file-tree reader for multi-kind dnaseq layouts

* **`dnaseq_create()`** (exported) — assembles a canonical dnaseq
  `data.frame` from one or more directory trees. Handles the three
  common on-disk layouts in one call: (1) everything combined per
  sample, (2) split by kind (SNV/Indel vs. all-CNV), (3) split by event
  (one file each for AMP / DEL / LOH / FUSION / LGR).
* Each `regex_*` argument accepts either a bare glob or a
  `list(loc, pattern)` — different alteration kinds can live under
  completely different directory trees.
* Rows read via kind-specific regexes (`regex_amp`, `regex_del`,
  `regex_loh`, `regex_focal_amp`, `regex_fusion`, `regex_lgr`) are
  stamped with the correct canonical `Variant_type` / `CNV_type` so
  the downstream classifier + oncoprint pick them up without any
  further code changes.
* `Source_kind` column added to output for auditability (tells you
  which regex each row came from).
* `ctdna_prepare()` and `ctdna_prep_add()` accept the extended dnaseq
  spec inline, e.g.:
  ```r
  ctdna_prepare(
    infinity_report = inf,
    dnaseq = list(loc       = "/study/dnaseq/snv",
                  regex_snv = "*_snv.tsv",
                  regex_amp = list(loc = "/study/dnaseq/cnv",
                                   pattern = "*_AMP.tsv"),
                  regex_del = list(loc = "/study/dnaseq/cnv",
                                   pattern = "*_DEL.tsv")),
    adam = adsl)
  ```
* Back-compat: the older `list(loc, regex)` shape maps to `regex_all`
  and continues to work.

## New — dnaseq alteration-type dictionary

* **`dnaseq_dict()`** (exported) — canonical mapping from raw vendor
  strings to ctdnaTM's alteration-type universe used by every oncoprint
  / concordance / stats function. Three sub-dictionaries: `variant_type`,
  `molecular_consequence`, `cnv_type`. API mirrors `ctdna_opts()`:
  ```r
  dnaseq_dict()                                        # list all
  dnaseq_dict("variant_type")                          # sub-dict
  dnaseq_dict("molecular_consequence.missense")        # one entry
  dnaseq_dict(molecular_consequence.missense = "MyVendor_MSs")  # append
  ```
* **`dnaseq_dict_add()`**, **`dnaseq_dict_reset()`**,
  **`dnaseq_dict_lookup()`** — companion helpers.
* Matching is **case-insensitive** and normalises
  whitespace/underscore/dash, so `"stop_gained"`, `"Stop Gained"`,
  `"STOP-GAINED"` all match. Combined consequences with `+` separator
  (e.g. `"frameshift_variant+splice_donor_variant+intron_variant"`)
  are split, each part looked up, and the MOST-DAMAGING component
  wins via a fixed priority order.
* Shipped defaults cover all 52 real-world `Variant_Effect` strings
  observed in production dnaseq data, plus `AMP` / `DEL` in CNA Type
  and `NONE` / `SILENT` / `MISSENSE` / `NONSENSE` in Functional_Class.
* **Automatic normalisation in `ctdna_prepare()` and
  `ctdna_prep_add()`** — the internal `.prep_dnaseq()` now runs
  `.dnaseq_dict_apply()` on every dnaseq frame before it lands in
  `prep$dnaseq`. Raw vendor values are preserved as `Variant_type_raw`,
  `Molecular_consequence_raw`, `CNV_type_raw` columns for auditability.
  Canonical values land in the standard-named columns that
  `.oncoprint_classify()` reads.

## New — ctDNA vs Tissue concordance oncoprint

* **`ctdna_concordance_oncoprint_core()`** (exported) — the transposed
  variant of `ctdna_oncoprint()` for paired ctDNA + Tissue analysis.
  Rows are patient pairs (each patient's ctDNA row directly above its
  Tissue row). Top annotation per panel: two side-by-side stacked
  alteration bars per gene, one for ctDNA and one for Tissue, each with
  its own % label using its own sample count as denominator. Right
  annotation per panel: per-patient Jaccard concordance bar with a
  perpendicular n-altered actual-count scale line.
* Panel title, gene-set annotation strip, left patient annotations
  (Dose, RECIST), and gene column names at the bottom all match the
  original `ctdna_oncoprint()` layout, just transposed.
* Alteration classes inherit LGR at 75% cell height, thin classes
  (Missense/Truncating/InFrame/Promoter) at 33%, CNV/Fusion at full
  cell, z-order CNV → LGR → Truncating → SNV top from the shared
  `.oncoprint_alter_fun`.

## Removed

* `ctdna_opts("dnaseq_loc")` and `ctdna_opts("dnaseq_regex")` — the
  options-driven fallback for dnaseq tree reading. Pass tree specs
  inline via `ctdna_prepare(dnaseq = list(...))` instead. RNAseq
  equivalents (`rnaseq_loc`, `rnaseq_regex`) are unchanged.

# ctdnaTM 0.75.10

## New

* `ctdna_concordance_oncoprint()` — transposed oncoprint pairing ctDNA calls
  with orthogonal tissue-genomic calls per patient. Genes in columns, patients
  in rows with paired ctDNA / tissue tracks. Full annotation surface matching
  `ctdna_oncoprint()`: patient annotations, per-gene frequency bar, alteration
  filter, cohort restriction via `indications`. Requires `patchwork` and
  `ggnewscale`.
* `cancer_50` — exported character vector of 50 pan-cancer genes.
* Multi-file dnaseq / rnaseq reader wired into `ctdna_prepare()`. Pass
  `dnaseq = list(loc = "dir", regex = "*.tsv")` and prep will glob the tree,
  read every TSV/CSV, canonicalize columns via `ctdna_opts()`, apply optional
  QC filter and min-VAF, and concatenate. Progress bar via
  `utils::txtProgressBar`.
* New `ctdna_opts()` keys for tissue-DNA / RNAseq column mapping:
  `dnaseq_loc`, `dnaseq_regex`, `dnaseq_gene_col` (default `"Gene Symbol"`),
  `dnaseq_variant_type_col` (`"Variant Type"`), `dnaseq_vaf_col`
  (`"Allelic Fraction"`), `dnaseq_variant_effect_col`,
  `dnaseq_protein_variant_col`, `dnaseq_functional_class_col`,
  `dnaseq_somatic_status_col`, `dnaseq_qc_filter_col` (`"Filter"`),
  `dnaseq_min_vaf` (default `0`). Same set with `rnaseq_` prefix.

## Changed

* **`indications` argument added to every user-facing plot function** — 12
  functions total. One-line cohort restriction (e.g.
  `ctdna_boxplot(prep, y = "tf", x = "Visit_name", indications = "NSCLC")`).
  Every function's signature is otherwise byte-identical to v0.75.9 — no
  arguments reordered, no defaults changed, no other behavior touched.

  Functions patched: `ctdna_plot_baseline`, `ctdna_plot_reduction`,
  `ctdna_plot_longitudinal`, `ctdna_plot_scatter`,
  `ctdna_plot_pct_change_by_dose_time`, `ctdna_boxplot`, `ctdna_oncoprint`,
  `ctdna_alteration_grid`, `ctdna_plot_vs_expression`,
  `ctdna_plot_vs_mutation`, `ctdna_plot_vs_tmb`, `ctdna_plot_vs_ihc`.

  Two shared internal helpers back all of them:
  `.ctdna_filter_prep_by_indication()` for prep-taking functions and
  `.ctdna_filter_df_by_indication()` for the integration `vs_*` functions.

## Depends

* `ggnewscale` added to `Imports:` (required by
  `ctdna_concordance_oncoprint()` for per-annotation legend splitting).

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
