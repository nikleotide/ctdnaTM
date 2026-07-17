# ctdnaTM 0.75.10

## New

* **`ctdna_concordance_oncoprint()`** — transposed oncoprint pairing ctDNA
  calls with orthogonal tissue-genomic calls (WES/WGS/panel) per patient.
  Genes in columns, patients in rows with two stacked tracks per patient
  (ctDNA on top, tissue underneath). Same glyphs and palette as
  `ctdna_plot_oncoprint()`, and the same annotation surface —
  `patient_data` / `annotations` (clinical strips at left), `show_freq_bar`
  (per-gene prevalence bar at top), `alterations` (alteration-type filter),
  `sort_genes`, `indications`, `subtitle` / `caption` / `legend_position`.
  Returns a ggplot (via patchwork) with class
  `ctdna_concordance_oncoprint / patchwork / gg / ggplot`; auxiliary data
  attached via attributes (`ctdna_matrix`, `tissue_matrix`, `plot_dim`,
  `patients`, `genes`, `long`, `n_dropped_by_topN`). Requires the
  `patchwork` and `ggnewscale` packages.

* **`cancer_50`** — exported character vector of 50 pan-cancer genes
  commonly used as a landscape set (solid-tumor drivers, DDR/HRR, heme
  drivers, RTK/fusion targets, chromatin regulators). Use in
  `ctdna_concordance_oncoprint()` or any function that takes a `gene_sets`
  argument.

* **Multi-file readers in `ctdna_prepare()`.** The `dnaseq` and `rnaseq`
  slots now accept a two-piece spec `list(loc = "dir", regex = "*.tsv")`.
  The reader walks the tree, matches by glob (dots literal, `*` = any
  non-separator char), extracts Sample_ID from the filename wildcard,
  canonicalizes columns via `ctdna_opts()` keys, and concatenates. Progress
  bar via `utils::txtProgressBar`. QC filter (default column `"Filter"` =
  `"PASS"`) and minimum-VAF (`dnaseq_min_vaf`, default `0`) applied at
  read-time if configured.

* **New `ctdna_opts()` keys** for orthogonal-modality column mapping:
  `dnaseq_loc`, `dnaseq_regex`, `dnaseq_gene_col` (default
  `"Gene Symbol"`), `dnaseq_variant_type_col` (`"Variant Type"`),
  `dnaseq_vaf_col` (`"Allelic Fraction"`), `dnaseq_variant_effect_col`
  (`"Variant Effect"`), `dnaseq_protein_variant_col` (`"Protein Variant"`),
  `dnaseq_functional_class_col` (`"Functional Class"`),
  `dnaseq_somatic_status_col` (`"Somatic_status"`, optional),
  `dnaseq_qc_filter_col` (`"Filter"`, optional), `dnaseq_min_vaf`
  (default `0`). Same set with `rnaseq_` prefix
  (`rnaseq_gene_col`, `rnaseq_value_col`, `rnaseq_qc_filter_col`, plus
  `rnaseq_loc` and `rnaseq_regex`).

## Changed

* **`indications` argument added to every user-facing plot function**
  (16 total) — one-line cohort restriction without needing to pre-filter
  the input data.frame or the `prep` object.

  Applied to: `ctdna_plot_oncoprint`, `ctdna_plot_alteration_grid`,
  `ctdna_plot_baseline`, `ctdna_plot_baseline_dose_recist`,
  `ctdna_plot_reduction`, `ctdna_plot_mr_recist`, `ctdna_plot_longitudinal`,
  `ctdna_plot_longitudinal_recist`, `ctdna_plot_reduction_vs_tumor`,
  `ctdna_plot_mut_methyl`, `ctdna_plot_pct_change_by_dose_time`,
  `ctdna_boxplot`, `ctdna_plot_vs_tmb`, `ctdna_plot_vs_mutation`,
  `ctdna_plot_vs_ihc`, `ctdna_plot_vs_expression`.

  A single shared internal helper `.ctdna_filter_by_indication()` backs all
  of them so behavior is identical everywhere: same column-lookup fallback
  (`Indication` -> `indication` -> `Cancertype` -> `Cancer_Type`, first in
  the frame; if not present, join through `prep$clinical`), same message
  when patients are dropped. **No other behavior of these functions
  changed.**

* **Roxygen documentation sweep.** Every exported function's `@param`
  block now covers every argument, including the new `indications` arg,
  so `?function_name` and RStudio tab-autocomplete surface the full
  argument set.

## Depends

* `ggnewscale` added to `Imports:` (required by
  `ctdna_concordance_oncoprint()` for per-annotation legend splitting).
  `patchwork` was already in `Imports:` and is now used more broadly.
