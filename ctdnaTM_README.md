# ctdnaTM

**ctDNA analysis toolkit**

Standardised plotting, statistics, and multi-modal integration helpers
for the eight core ctDNA deliverables. Built around the
**Guardant Health Infinity** platform but works with any vendor that
produces the same data shape.


---

## Installation

```r
# install.packages("remotes")
remotes::install_github("hnikbakht/ctdnaTM")

# Strongly recommended for publication-quality oncoprints:
install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

library(ctdnaTM)
# > ctdnaTM v0.19.0 - ctDNA Deliverables

```

## Quick start

```r
sim <- ctdna_make_mock_study(n_patients = 40, seed = 1)
mm  <- ctdna_make_example_multimodal(n_subjects = 40, seed = 1)

out <- ctdna_run_pipeline(
  ctdna_df = mm$ctdna,
  expr_df  = mm$expression,  tmb_df   = mm$tmb,
  ihc_df   = mm$ihc,         g74_df   = mm$genomic_74,
  sim      = sim,
  legend_position = "bottom"
)
out$baseline           # D1
out$mr_recist          # D4
out$oncoprint          # ComplexHeatmap-backed genomic landscape
out$alteration_grid    # Combinatorial alteration grid
```

---

# What's new in v0.19.0

- **`ctdna_plot_oncoprint()` now uses `ComplexHeatmap::oncoPrint()`** when
  installed: real publication-quality output with proper legends, right-
  margin per-gene barplot, top per-patient burden bar, stacked annotation
  tracks, and gene-set row splits. Falls back to a clean ggplot tile
  plot when ComplexHeatmap is not installed.
- **No-overlap audit** across every plot: n-labels moved under x-axis ticks
  (no longer overlap data), long stats subtitles soft-wrap on `;`/`|`/` - `
  separators, MR-threshold line label moved left of data, the "ctDNA
  clearance" annotation no longer repeats in every facet of D6, the
  alteration-grid bar labels no longer stack on top of each other.

# What's new in v0.12.0

- **Universal `scheme` argument** across every RECIST-aware plot.
  Statistics follow the chosen scheme: `"two"` (R vs NR), `"three"`
  (CR/PR vs SD vs PD/NE/NA), `"four"` (CR/PR vs uCR/PR vs SD vs PD/NE/NA).
- **`ctdna_metric_at()`** + new `at` / `metric` / `mr_threshold` arguments
  on `ctdna_plot_reduction()` — compute per-subject ctDNA value as
  ratio / raw TF / percent change at "best" / "last" / a named visit.

---

# Function reference

48 exported functions, all prefixed with `ctdna_`.

## 1. Configuration — `ctdna_opts()`

Single entry point for all package configuration. Five modes:

```r
ctdna_opts()                                 # print full registry
ctdna_opts("tf")                             # get a single key
ctdna_opts(c("tf","maxaf"))                  # get multiple keys
ctdna_opts(tf = "TF_v3", display_floor = 1e-3)   # set
ctdna_opts(.reset = TRUE)                    # reset to defaults
ctdna_opts(load = cfg)                       # bulk-load from a list
```

**Save & restore a study config:**

```r
ctdna_opts(tf = "TF_methyl_v3", maxaf = "MAX_AF")
saveRDS(ctdna_opts(), "study_config.rds")
# In a fresh session:
cfg <- readRDS("study_config.rds")
ctdna_opts(load = cfg)
```

**Notable keys.** Column names (`subject`, `time`, `dose`, `recist`,
`tf`, `maxaf`, `tmb`, `mr`, ...); numeric thresholds (`display_floor`,
`loq`); plot defaults (`stat_position`, `legend_position`, `show_stats`,
`show_n`); units (`expression_unit` = `"log2_tpm_plus_one"`); multi-modal
column names.

---

## 2. Mock data — `ctdna_make_mock_study()` / `ctdna_make_example_multimodal()`

```r
# Full vendor bundle (5 Guardant Infinity files + clinical + IHC + WES + RNA-seq)
sim <- ctdna_make_mock_study(
  n_patients   = 40,
  visits       = c("C1D1","C2D1","C3D1","C5D1","EOT"),
  cancer_types = c("NSCLC","CRC","BRCA","mCRPC","GBM"),
  doses        = c("Low","Mid","High"),
  seed         = 42)

# Canonical multi-modal frames (faster, already in canonical shape)
mm <- ctdna_make_example_multimodal(n_subjects = 30, seed = 42)
```

---

## 3. Vendor → canonical — `ctdna_prepare()`

```r
prepped <- ctdna_prepare(
  infinity_report    = sim$infinity_report,
  tf_change          = sim$tf_change,
  panel_74_response  = sim$panel_74_response,
  panel_500_response = sim$panel_500_response,
  clinical           = sim$clinical,
  qc_filter          = TRUE,
  verbose            = FALSE)
# names(prepped): "ctdna", "genomic_74", "genomic_500"
```

`ctdna_assemble_modalities()` validates required columns and assembles
a multi-modal database (auto-canonicalizes expression to log2(TPM+1)).

---

## 4. Filtering — `ctdna_filter_alterations()`

Unified filter for any alteration table. Every parameter accepts `FALSE`
to disable that filter:

```r
ctdna_filter_alterations(
  sim$infinity_report,
  variant_type      = c("SNV","Indel","CNV","Fusion","LGR"),
  somatic_status    = c("Somatic","Likely somatic"),
  functional_impact = c("High","Moderate"),
  max_gnomad_af     = 0.01,
  min_vaf           = 0.5,
  min_cn_amp        = 4,
  max_cn_loss       = 1,
  verbose           = FALSE)
```

---

## 5. Helpers

| Function | Purpose |
|---|---|
| `ctdna_floor_tf(x)` | Replace zeros with `display_floor` for log plots |
| `ctdna_stratify_recist(x, scheme)` | Collapse to "two" / "three" / "four" categories |
| `ctdna_ratio(df)` | Per-subject per-timepoint TF / baseline TF ratios |
| `ctdna_summary(df)` | Per-subject baseline / best / last values |
| `ctdna_compute_mr(df, at, threshold)` | Molecular response classifier |
| `ctdna_metric_at(df, metric, at)` | Per-subject ratio / TF / pct_change at best / last / a named visit |
| `ctdna_merge_modalities(ctdna_df, modality_df, by)` | Merge per-subject |
| `ctdna_to_log2tpm(df)` | Convert TPM ↔ log2(TPM+1) |

---

## 6. Statistics

| Function | Purpose |
|---|---|
| `ctdna_pval_label(p)` | Vectorized `"p=0.012"` / `"p<0.001"` |
| `ctdna_run_group_test(df, value, group, method)` | Auto / Wilcoxon / Kruskal / ANOVA / t |
| `ctdna_overall_group_test(df, value, group)` | Overall test |
| `ctdna_pairwise_ranksum(df, value, group, method)` | Pairwise Wilcoxon |
| `ctdna_cor_test(x, y, method)` | Spearman / Pearson / Kendall |
| `ctdna_paired_test(df, value, time, method)` | Paired Wilcoxon / t |
| `ctdna_cohen_kappa(tab)` | Cohen's kappa for a 2x2 |
| `ctdna_concordance_test(tab, method)` | Agreement % + kappa + (2x2) McNemar |

---

## 7. Display options (every plot accepts these)

| Option | Default | Effect |
|---|---|---|
| `scheme` | `"two"` / `"three"` / `"four"` | RECIST stratification (drives both colours AND statistics) |
| `stat_position` | `"subtitle"` | Stats text placement: `"subtitle"` / `"caption"` / `"on_plot"` / `"none"` |
| `legend_position` | `"right"` | `"right"`/`"left"`/`"top"`/`"bottom"`/`"none"` |
| `show_stats` | `TRUE` | Compute / display the stats label |
| `show_n` | `TRUE` | Annotate group sizes under x-axis ticks (no overlap with data) |
| `title`, `subtitle`, `caption`, `xlab`, `ylab` | — | Standard ggplot labels |

Set as global default via `ctdna_opts(stat_position = "caption", legend_position = "bottom")`.

Long stats subtitles automatically soft-wrap onto multiple lines on
`;`, `|`, or ` - ` separators (width ≈ 70 chars).

---

## 8. The 8 core deliverable plots

| Function | Purpose |
|---|---|
| **D1** `ctdna_plot_baseline(df, group, vars, scheme, stat, ...)` | Boxplot of baseline values grouped by any clinical variable |
| **D2** `ctdna_plot_baseline_dose_recist(df, scheme, stat, ...)` | Baseline methylTF dodged by dose × RECIST |
| **D3** `ctdna_plot_reduction(df, scheme, at, metric, mr_threshold, ...)` | Per-patient ctDNA reduction with MR threshold |
| **D4** `ctdna_plot_mr_recist(df, mr_timepoint, mr_threshold, stat, ...)` | MR↔RECIST 2×2 concordance heatmap |
| **D5** `ctdna_plot_longitudinal(df, stat, ...)` | Per-subject spaghetti + median per dose |
| **D6** `ctdna_plot_longitudinal_recist(df, scheme, stat, ...)` | Longitudinal trajectory faceted by dose × RECIST |
| **D7** `ctdna_plot_reduction_vs_tumor(df, scheme, stat, ...)` | Scatter of ctDNA best-ratio vs % tumor change (RECIST=colour, dose=shape) |
| **D8** `ctdna_plot_mut_methyl(df, stat, ...)` | Per-baseline maxAF vs methylTF |

D4 returns `$plot`, `$table_plot` (publication-style contingency table),
`$agreement`, `$kappa`, `$mcnemar_p`.

### D3 — three metric scales × three timepoint modes

```r
# ratio (default): TF_at / TF_baseline. 0.1 = TF dropped to 10% of baseline.
#   Pairs with mr_threshold (default 0.5 = "TF dropped by at least 50%").
ctdna_plot_reduction(ctd, scheme = "three", metric = "ratio",
                       at = "best", mr_threshold = 0.5)

# Absolute TF at a specific visit (log y)
ctdna_plot_reduction(ctd, scheme = "three", metric = "tf", at = "C5D1")

# Percent change at the last on-treatment visit (linear y)
ctdna_plot_reduction(ctd, scheme = "three", metric = "pct_change", at = "last")
```

---

## 9. Cross-modality plots

All four accept `ctdna_df`, the modality data frame, and standard
display options. They use `ctdna_metric` (default `"best_ratio"`) as
the per-subject ctDNA summary.

**What is `best_ratio`?** For each subject across their on-treatment
visits, `TF(visit) / TF(baseline)` and pick the **lowest** ratio — the
patient's deepest molecular response. `0.1` means TF dropped to 10% of
baseline.

```r
ctdna_plot_vs_expression(mm$ctdna, mm$expression, gene = "PDL1")
ctdna_plot_vs_ihc(mm$ctdna, mm$ihc, marker = "PDL1")
ctdna_plot_vs_mutation(mm$ctdna, mm$genomic_74, gene = "TP53", panel = "74genes")
ctdna_plot_vs_tmb(mm$ctdna, mm$tmb)
```

### `ctdna_mutation_concordance(mut_df, panel_a, panel_b, ...)`

Per-(subject × gene) call agreement between two assays on the same
patients (e.g. 74-gene vs 500-gene panels). Plot shows per-gene
agreement % (sorted, vertical labels with n underneath each gene);
subtitle reports overall agreement / Cohen's kappa / N.

```r
mut_long <- rbind(mm$genomic_74, mm$genomic_500)
mc <- ctdna_mutation_concordance(mut_long,
  panel_a = "74genes", panel_b = "500genes", stat = "both")
mc$plot         # bar chart
mc$table_plot   # publication-style 2x2
```

---

## 10. Genomic landscape

### `ctdna_plot_oncoprint(df, gene_sets, patient_data, top_annotations, ...)`

When **ComplexHeatmap** is installed (recommended), uses
`ComplexHeatmap::oncoPrint()` under the hood: per-gene right-margin
barplot, top per-patient burden bar, stacked annotation tracks, row
splits per gene set, and a clean merged legend at the bottom. Falls
back to a ggplot tile plot otherwise.

```r
ctdna_plot_oncoprint(sim$infinity_report,
  gene_sets = list(
    RTK_RAS = c("KRAS","EGFR","BRAF","MET"),
    PI3K    = c("PIK3CA","PTEN"),
    TP53    = "TP53"),
  patient_data    = sim$clinical,
  top_annotations = c("RECIST","Dose","Cancertype"),
  filter          = list(clinvar = c("Pathogenic","Likely pathogenic"),
                          functional_impact = FALSE,
                          verbose = FALSE))
```

**Engine selection.** `engine = "auto"` (default) picks ComplexHeatmap
if available; `"complexheatmap"` forces it; `"ggplot"` forces the
fallback. Returns a `ctdna_oncoprint` object; `print()` it to draw.

**Gene-set API.** Pass `gene_sets` as a *named list*:
- `NULL` (default) → all genes under one heading.
- Single-entry list → all those genes under that name.
- Multi-entry list → one row split per set with a side label.

### `ctdna_plot_alteration_grid(df, gene_sets, patient_data, ...)`

Stacked-bar grid: rows = indication, columns = gene set, x within each
panel = response category, fill = alt vs wt, with Fisher p-values per
panel and `%alt (n_alt / n_total)` labels above each bar.

```r
ctdna_plot_alteration_grid(sim$infinity_report,
  gene_sets = list(KRAS = "KRAS",
                   HRR  = c("BRCA1","BRCA2","ATM","CHEK2","PALB2"),
                   TP53 = "TP53"),
  patient_data    = sim$clinical,
  response_col    = "RECIST",
  indication_col  = "Cancertype",
  filter          = list(functional_impact = c("High","Moderate"),
                          verbose = FALSE))$plot
```

---

## 11. Theme, scales, palettes

```r
ctdna_colors()              # named character: CR/PR, uCR/PR, SD, PD/NE/NA, R, NR
ctdna_dose_palette(n = 3)

ctdna_theme()               # base theme
ctdna_scale_recist()        # colour scale for RECIST
ctdna_scale_recist_fill()
ctdna_scale_dose()          # colour scale for dose (alone)
ctdna_scale_dose_fill()
ctdna_scale_dose_shape()    # shape scale for dose (when RECIST is on colour)
```

**Visual encoding convention** (every plot follows this):

| Feature | Aesthetic |
|---|---|
| RECIST | colour (always when shown) |
| Dose + RECIST shown together | shape (dose) |
| Dose alone | colour |

---

## 12. Publication-ready tables — `ctdna_table()`

Renders any data frame / matrix / `table` as a styled ggplot table
panel with title, subtitle, caption, alternating row stripes, bold
header with underline, top/bottom rules, and optional diagonal-cell
tint:

```r
tab <- as.data.frame(matrix(c(12, 3, 1, 9), 2, 2,
  dimnames = list(MR = c("Response","NonResponse"),
                  RECIST = c("CR/PR","SD/PD/NE"))))
ctdna_table(tab,
  title          = "MR × RECIST contingency",
  subtitle       = "n = 25; Cohen's kappa = 0.62",
  caption        = "MR computed at C5D1; cutoff = 50% TF drop",
  highlight_diag = TRUE)
```

---

## 13. End-to-end pipeline — `ctdna_run_pipeline()`

```r
out <- ctdna_run_pipeline(
  ctdna_df = mm$ctdna,
  expr_df  = mm$expression,  tmb_df = mm$tmb,
  ihc_df   = mm$ihc,         g74_df = mm$genomic_74,
  sim      = sim,
  stat_position   = "subtitle",
  legend_position = "bottom",
  scheme   = "three",
  compose  = FALSE,
  verbose  = TRUE)

out$baseline; out$mr_recist; out$oncoprint; out$alteration_grid
```

Per-plot overrides via `plot_overrides`:

```r
out <- ctdna_run_pipeline(
  ctdna_df = mm$ctdna,
  plot_overrides = list(
    plot_baseline = list(vars = c("methylTF","maxAF"), stat = "wilcox"),
    plot_vs_tmb   = list(panel = "500genes")
  ))
```

---

## 14. Expression units

The canonical unit is **log2(TPM+1)** throughout. If your data is raw
TPM, tag the data frame with `attr(df, "units") <- "tpm"` and
`ctdna_to_log2tpm()` will convert it. `ctdna_assemble_modalities()`
calls this automatically.

---

## Citation

Nikbakht H. (2026). _ctdnaTM: ctDNA Deliverables_.
R package version 0.19.0. https://github.com/hnikbakht/ctdnaTM

## License

MIT.

---

