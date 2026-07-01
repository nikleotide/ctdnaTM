# ctdnaTM

Standardised plotting, statistics, and multi-modal integration helpers for the eight core ctDNA deliverables. Built around the Guardant Health Infinity platform, extensible to other vendors.

Suggested — install `ComplexHeatmap` for publication-quality oncoprints:

```r
install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
```

## Install

```r
install.packages("devtools")
devtools::install_github("nikleotide/ctdnaTM")
```

## Five-minute tour

```r
library(ctdnaTM)

# 1. Mock study (real usage: replace with your infinity_report + ADaM)
sim  <- ctdna_make_mock_study(n_patients = 60, seed = 1)

# 2. Assemble the prep object
prep <- ctdna_prepare(
  infinity_report = sim$infinity_report,
  adam            = list(adsl = sim$clinical))

# 3. Sample QC (idempotent; drops failing samples + their variants)
prep <- ctdna_sample_qc(prep)

# 4. Variant filtering with the two mandatory flags
prep <- ctdna_variant_filter(prep,
          filter_scheme = c("TSG_Tier1", "HRR14"),
          apply         = TRUE,   # drop failing rows
          explain       = TRUE)   # store annotated snapshot

# 5. Landscape plots
prep$adrs <- data.frame(Patient_ID = prep$clinical$Patient_ID,
                        indication = "NSCLC")
op <- ctdna_oncoprint(prep,
        gene_sets = list(TSG   = c("TP53","RB1","PTEN"),
                         HRR14 = ctdna_gene_set("HRR14")),
        engine    = "ggplot")
print(op$plot)

# 6. Or run the whole deliverable pipeline
res <- ctdna_pipeline(prep)
res$oncoprint$plot
res$alteration_grid$plot
```

## Filter surface

Four public functions:

| Function | Purpose |
|:--|:--|
| `ctdna_sample_qc(prep)` | Sample QC + cascade to variants. Idempotent. |
| `ctdna_variant_filter(prep, filter_scheme, apply, explain)` | Variant filter. Both flags mandatory. |
| `ctdna_create_scheme(name, genes, criteria, combine)` | User-friendly scheme builder. |
| `create_filtering_scheme(...)` | Full grammar builder (used by all built-ins). |

## Built-in schemes

| Scheme | Scope | Notes |
|:--|:--|:--|
| `TSG_Tier1()` | TP53 / RB1 / PTEN | Most permissive. Any impactful somatic short variant, fusion, or deleterious LGR. |
| `TSG_Tier2()` | TP53 / RB1 / PTEN | Intermediate. Somatic LoF, ClinVar P/LP + missense/inframe, LOH or homozygous del, intra-chromosomal LGR. |
| `TSG_Tier3()` | TP53 / RB1 / PTEN | Most restrictive. Germline P/LP, somatic LoF, ClinVar P/LP + missense/inframe, LOH or homozygous del. No LGR, no fusion, no rare-syn carve-in. |
| `HRR14()` | HRR14 gene list | ClinVar P/LP, truncating + rare + not BRCA2 K3326X, homozygous del, somatic deleterious intra-genic LGR. |
| `scheme_basic()` | all genes | GH Infinity "impactful alterations" rules. |

All four biomarker schemes exclude ClinVar benign and CHIP unconditionally. Each scheme's body is broken into **named branches** (`germline_PLP`, `somatic_LoF`, `intra_gene_LGR`, ...) that show up in the `filtering_criteria` column as `PASS::TSG_Tier1(somatic_nonsyn)`.

## Grammar (build your own schemes)

```r
sch <- create_filtering_scheme(
  # Scope
  rule(Gene = c("TP53","RB1","PTEN")),
  # Global exclusions
  not(rule_clinvar_benign()),
  not(rule_ch()),
  # Body: any-of branches
  anyOf(
    allOf(rule_germline(), rule_clinvar_path()),
    allOf(rule_somatic(),  rule_truncating()),
    allOf(rule_somatic(),  rule_is_LGR_deleterious())),
  name = "my_tsg", overwrite = TRUE)

prep <- ctdna_variant_filter(prep, "my_tsg", apply = TRUE, explain = TRUE)
```

Combinators: `rule()`, `allOf()`, `anyOf()`, `not()`. About 60 pre-built predicates (`rule_somatic`, `rule_truncating`, `rule_clinvar_path`, `rule_is_LGR_deleterious`, etc.); see `?ctdna_rules_library`.

## Documentation

- Function tour vignette: `vignette("ctdnaTM_function_tour", package = "ctdnaTM")`
- Every exported function has a help page: `?ctdna_variant_filter`, `?ctdna_oncoprint`, `?TSG_Tier1`, ...
- Changelog: `NEWS.md`
