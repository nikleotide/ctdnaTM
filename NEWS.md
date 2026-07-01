# ctdnaTM 0.75.3

* **Pipeline registry corrections** — cross-modality steps re-pointed at the frames that actually hold the data:
  * `vs_tmb` needs `$samples` only (TMB values are a column on `prep$samples`).
  * `vs_mutation` needs `$samples` **and** `$variants` (mutation status comes from variants, not from an external `$genomic_74` frame).
  * `vs_ihc` and `vs_expression` still need `$samples` plus the user-attached `$ihc` / `$expression` frames.
* Every version now ships the full Rmd/md set on every iteration: `README.md`, `NEWS.md`, three vignettes (`walkthrough`, `function_tour`, `filter_internals`), the core deliverables Rmd, and the SOW36 analysis starter.

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
* Fresh SOW36 analysis starter for real data.

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

* Biomarker schemes redone per the GSK ctDNA doc.
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
