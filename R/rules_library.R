# =============================================================================
# v1.0.0 -- Pre-shipped rule library + built-in schemes
# =============================================================================
#
# This file defines the pre-built ctdna_rule atoms users can plug into
# their schemes without writing them from scratch, plus the two
# canonical built-in schemes (HRR and TSG) defined per the published
# specifications.
#
# Naming convention:
#   rule_<column>_<value>  -- single-column constraint on a canonical value
#   rule_genes_<name>      -- gene-set membership rule
#   rule_ind_<name>        -- indication rule (uses ctdna_opts cancertype_dictionary)
#   rule_<composite>       -- multi-rule composite (e.g. rule_truncating)
#   scheme_<name>          -- full filtering scheme
#
# All rule objects are exported so users can use them directly.
# Rules are constructed lazily inside .ctdna_rule_library() so that
# changes to ctdna_opts (e.g. col_clinvar) are picked up at use time
# rather than load time.

# -----------------------------------------------------------------------------
# Variant_type rules (5 values: SNV, Indel, CNV, Fusion, LGR)
# -----------------------------------------------------------------------------

#' Pre-built rule predicates for ctdna_filter
#'
#' A family of zero-argument constructors that return
#' \code{ctdna_rule} objects matching common single-column constraints
#' on the Guardant Infinity variant report. Combine them with
#' \code{\link{allOf}}, \code{\link{anyOf}}, \code{\link{not}}
#' to build filter schemes.
#'
#' Naming convention:
#' \itemize{
#'   \item \code{rule_is_<TYPE>()} -- \code{Variant_type} predicate
#'     (SNV, Indel, CNV, Fusion, LGR).
#'   \item \code{rule_lgr_<value>()} -- \code{LGR_subtype} predicate
#'     (deletion, tandem duplication, inversion).
#'   \item \code{rule_cnv_<value>()} -- \code{CNV_type} predicate
#'     (focal_amp, aneuploid_amp, amp, homozyg_del, loh_del).
#'   \item \code{rule_<molecular_consequence>()} -- \code{Molecular_consequence}
#'     predicate (missense, nonsense, frameshift, splice_donor,
#'     splice_acceptor, splice_region, start_lost, stop_lost,
#'     synonymous, promoter, inframe_indel, ...).
#'   \item \code{rule_sample_<status>()} -- \code{Sample_status} predicate
#'     (success, fail).
#'   \item \code{rule_<origin>()} -- \code{Somatic_status} predicate
#'     (germline, somatic, ch).
#'   \item \code{rule_clinvar_<class>()} -- ClinVar predicate
#'     (benign, path, vus, conflict).
#'   \item \code{rule_genes_<set>()} -- gene-set membership rule using
#'     the built-in constants (HRR14, TSG, RTK, Cell_Cycle,
#'     TP53_pathway, MMR, PI3K).
#'   \item \code{rule_ind_<indication>()} -- indication predicate
#'     (NSCLC, BRCA, SCLC, HNSCC, CRC, mCRPC, GBM) resolved via
#'     \code{ctdna_opts("cancertype_dictionary")}.
#'   \item \code{rule_rare_<AF>()} / \code{rule_common_<AF>()} --
#'     gnomAD allele-frequency thresholds (\code{rule_rare_001},
#'     \code{rule_rare_0001}, \code{rule_common_001}).
#'   \item Functional-impact atoms: \code{rule_deleterious},
#'     \code{rule_reversion}, \code{rule_biallelic}.
#'   \item Parameterised: \code{rule_amp_thresh(cn_threshold = 4)} --
#'     CNV amp with a minimum copy-number cutoff.
#'   \item Composite predicates (multiple atoms combined for
#'     convenience): \code{rule_truncating()},
#'     \code{rule_is_LGR_deleterious()},
#'     \code{rule_splice_event_mc()}.
#'   \item Curated one-offs (data-specific carve-outs from the
#'     biomarker specs): \code{rule_BRCA2_K3326()}.
#' }
#'
#' All rules are lazy: each constructor returns the rule fresh, so
#' changes to \code{ctdna_opts()} (e.g. \code{col_clinvar}) are picked
#' up at use time rather than binding at library-load time.
#'
#' @return A \code{ctdna_rule} object suitable for
#'   \code{\link{create_filtering_scheme}} or (via
#'   \code{\link{ctdna_variant_filter}}) direct application to a
#'   \code{ctdna_prep} object.
#'
#' @seealso \code{\link{allOf}}, \code{\link{anyOf}}, \code{\link{not}}
#'   (combinators); \code{\link{create_filtering_scheme}} (build a
#'   named scheme from rules); \code{\link{ctdna_create_scheme}}
#'   (a keyword-driven wrapper that hides the DSL);
#'   \code{\link{TSG_Tier1}}, \code{\link{HRR14}},
#'   \code{\link{scheme_basic}} (built-in schemes composed from these
#'   atoms).
#'
#' @examples
#' # Two ways to express "SNV missense not in ClinVar benign":
#' r1 <- allOf(rule_is_SNV(), rule_missense(),
#'             not(rule_clinvar_benign()))
#' print(r1)
#'
#' # Composite predicates (built from atoms):
#' print(rule_truncating())         # nonsense + frameshift + stop/start_lost + splice
#' print(rule_splice_event_mc())    # splice_donor / splice_acceptor / splice_region
#' print(rule_is_LGR_deleterious()) # LGR with a clinically meaningful subtype
#'
#' # Gene-set + AF-threshold combination:
#' r2 <- allOf(rule_genes_HRR14(), rule_somatic(), rule_rare_001())
#' print(r2)
#'
#' # Amp threshold - parameterised atom:
#' print(rule_amp_thresh(cn_threshold = 6))
#'
#' # Assemble into a named scheme; register in the session catalog.
#' sch <- create_filtering_scheme(
#'   rule_is_SNV(),
#'   anyOf(rule_missense(), rule_truncating()),
#'   not(rule_germline()),
#'   name = "snv_coding_somatic")
#' print(sch)
#' @name ctdna_rules_library
NULL


#' @rdname ctdna_rules_library
#' @export
rule_is_SNV     <- function() rule(Variant_type = "SNV")
#' @rdname ctdna_rules_library
#' @export
rule_is_Indel   <- function() rule(Variant_type = "Indel")
#' @rdname ctdna_rules_library
#' @export
rule_is_CNV     <- function() rule(Variant_type = "CNV")
#' @rdname ctdna_rules_library
#' @export
rule_is_Fusion  <- function() rule(Variant_type = "Fusion")
#' @rdname ctdna_rules_library
#' @export
rule_is_LGR     <- function() rule(Variant_type = "LGR")

# ---- LGR subtype rules ----------------------------------------------------
# LGR subtype is INFERRED from Direction_a / Direction_b columns
# (Guardant Health Infinity has no pre-computed LGR-subtype column).
# Pairing convention:
#   Direction_a == -1, Direction_b ==  1  =>  deletion
#   Direction_a ==  1, Direction_b == -1  =>  tandem_duplication
#   Direction_a == +-1, Direction_b same  =>  inversion
# Direction values are coerced to numeric at eval time so the rules
# work whether the source columns are stored as integer or character.
# Rows with NA / missing Direction columns return FALSE (drop).

#' @rdname ctdna_rules_library
#' @export
rule_lgr_deletion <- function()
  allOf(rule_is_LGR(),
          rule(Direction_a = -1),
          rule(Direction_b =  1))

#' @rdname ctdna_rules_library
#' @export
rule_lgr_tandem_duplication <- function()
  allOf(rule_is_LGR(),
          rule(Direction_a =  1),
          rule(Direction_b = -1))

#' @rdname ctdna_rules_library
#' @export
rule_lgr_inversion <- function()
  allOf(rule_is_LGR(),
          anyOf(
            allOf(rule(Direction_a =  1), rule(Direction_b =  1)),
            allOf(rule(Direction_a = -1), rule(Direction_b = -1))))

# CNV_type rules (5 canonical values)
#' @rdname ctdna_rules_library
#' @export
rule_cnv_focal_amp     <- function() rule(CNV_type = "focal_amplification")
#' @rdname ctdna_rules_library
#' @export
rule_cnv_aneuploid_amp <- function() rule(CNV_type = "aneuploid_amplification")
#' @rdname ctdna_rules_library
#' @export
rule_cnv_amp           <- function() rule(CNV_type = "amplification")
#' @rdname ctdna_rules_library
#' @export
rule_cnv_homozyg_del   <- function() rule(CNV_type = "homozygous_deletion")
#' @rdname ctdna_rules_library
#' @export
rule_cnv_loh_del       <- function() rule(CNV_type = "loh_deletion")

# Splice_effect rules
#' @rdname ctdna_rules_library
#' @export
rule_splice_acceptor <- function() rule(Splice_effect = "splice_acceptor_variant")
#' @rdname ctdna_rules_library
#' @export
rule_splice_donor    <- function() rule(Splice_effect = "splice_donor_variant")
#' @rdname ctdna_rules_library
#' @export
rule_splice_region   <- function() rule(Splice_effect = "splice_region_variant")

# Sample_status rules
#' @rdname ctdna_rules_library
#' @export
rule_sample_success  <- function() rule(Sample_status = "SUCCESS")
#' @rdname ctdna_rules_library
#' @export
rule_sample_fail     <- function() rule(Sample_status = "FAIL")

# Somatic_status rules (germline, somatic, somatic_putative_ch)
#' @rdname ctdna_rules_library
#' @export
rule_germline        <- function() rule(Somatic_status = "germline")
#' @rdname ctdna_rules_library
#' @export
rule_somatic         <- function() rule(Somatic_status = "somatic")
#' @rdname ctdna_rules_library
#' @export
rule_ch              <- function() rule(Somatic_status = "somatic_putative_ch")

# Molecular_consequence rules (single-value atoms)
#' @rdname ctdna_rules_library
#' @export
rule_missense        <- function() rule(Molecular_consequence = "missense")
#' @rdname ctdna_rules_library
#' @export
rule_synonymous      <- function() rule(Molecular_consequence = "synonymous")
#' @rdname ctdna_rules_library
#' @export
rule_nonsense        <- function() rule(Molecular_consequence = "nonsense")
#' @rdname ctdna_rules_library
#' @export
rule_frameshift      <- function() rule(Molecular_consequence = "frameshift")
#' @rdname ctdna_rules_library
#' @export
rule_start_lost      <- function() rule(Molecular_consequence = "start_lost")
#' @rdname ctdna_rules_library
#' @export
rule_stop_lost       <- function() rule(Molecular_consequence = "stop_lost")
#' @rdname ctdna_rules_library
#' @export
rule_promoter        <- function() rule(Molecular_consequence = "promoter")
#' @rdname ctdna_rules_library
#' @export
rule_inframe_indel   <- function()
  rule(Molecular_consequence = c("inframe_indel","inframe_insertion",
                                  "inframe_deletion","inframe_duplication"))
#' @rdname ctdna_rules_library
#' @export
rule_splice_event_mc <- function()
  rule(Molecular_consequence = c("splice_acceptor","splice_donor",
                                  "splice_region","splice_event"))

# Composite: LGR + deleterious Functional_impact
# (Per v1.0.2 spec: LGRs are kept only when annotated as deleterious.)
#' @rdname ctdna_rules_library
#' @export
rule_is_LGR_deleterious <- function()
  allOf(rule_is_LGR(), rule_deleterious())

# Amplification with a copy-number threshold (default 4).
# Use with `cn_threshold` to override; >=4 by default.

#' Amplification rule with a copy-number threshold
#'
#' Composite rule that matches focal, aneuploid, OR generic
#' amplification calls AND requires \code{Copy_number} to be at or
#' above \code{cn_threshold}.
#'
#' @param cn_threshold Minimum integer copy number required. Default
#'   \code{4}.
#' @return A \code{ctdna_rule} object.
#' @rdname ctdna_rules_library
#' @export
rule_amp_thresh <- function(cn_threshold = 4)
  allOf(anyOf(rule_cnv_focal_amp(),
                  rule_cnv_aneuploid_amp(),
                  rule_cnv_amp()),
          rule(Copy_number = list(op = ">=", value = cn_threshold)))

# Truncating / null variants -- the canonical LoF set.
# v0.37.0: added stop_lost (nonstop) to align with HRR + TSG specs.
#' @rdname ctdna_rules_library
#' @export
rule_truncating <- function() {
  rule(Molecular_consequence = c("frameshift","nonsense","start_lost",
                                  "stop_lost","splice_acceptor","splice_donor"))
}

# ClinVar rules (canonical pathogenicity buckets)
#' @rdname ctdna_rules_library
#' @export
rule_clinvar_benign <- function()
  rule(ClinVar = c("Benign","Likely_benign","Benign/Likely_benign"))
#' @rdname ctdna_rules_library
#' @export
rule_clinvar_path <- function()
  rule(ClinVar = c("Pathogenic","Likely_pathogenic",
                    "Pathogenic/Likely_pathogenic",
                    "Pathogenic/Likely_pathogenic,_other",
                    "Likely_pathogenic,_other",
                    "Likely_pathogenic,_risk_factor"))
#' @rdname ctdna_rules_library
#' @export
rule_clinvar_vus <- function() rule(ClinVar = "Uncertain_significance")
#' @rdname ctdna_rules_library
#' @export
rule_clinvar_conflict <- function()
  rule(ClinVar = "Conflicting_interpretations_of_pathogenicity")

# Functional_impact rules
#' @rdname ctdna_rules_library
#' @export
rule_deleterious <- function() rule(Functional_impact = "deleterious")
#' @rdname ctdna_rules_library
#' @export
rule_reversion   <- function() rule(Functional_impact = c("reversion","reversion_cis"))

# Mutant_allele_status rules
#' @rdname ctdna_rules_library
#' @export
rule_biallelic   <- function() rule(Mutant_allele_status = "biallelic")

# Population-frequency rules (gnomAD-based; require ctdna_annotate_population_freq)
# v0.38.0: rule_rare_001() now honors ctdna_opts("rare_variant_filter") and
# ctdna_opts("rare_variant_threshold"). When rare_variant_filter = FALSE,
# the rule effectively always passes (matches everything regardless of
# gnomAD_AF), letting users disable population-frequency filtering
# without rebuilding schemes.
#' @rdname ctdna_rules_library
#' @export
rule_rare_001 <- function() {
  if (isTRUE(.o("rare_variant_filter"))) {
    thresh <- .o("rare_variant_threshold") %||% 0.01
    rule(gnomAD_AF = list(op = "<", value = thresh))
  } else {
    # rare_variant_filter disabled -> always-pass node
    structure(list(type = "true"),
              class = c("ctdna_rule_combinator","ctdna_rule"))
  }
}
#' @rdname ctdna_rules_library
#' @export
rule_rare_0001   <- function() rule(gnomAD_AF = list(op = "<",  value = 0.001))
#' @rdname ctdna_rules_library
#' @export
rule_common_001  <- function() rule(gnomAD_AF = list(op = ">=", value = 0.01))


# -----------------------------------------------------------------------------
# Gene-set rules
# -----------------------------------------------------------------------------

# HRR14 (per the HRR scheme spec)
.HRR14_GENES <- c("BRCA1","BRCA2","ATM","BARD1","BRIP1","CDK12","CHEK1","CHEK2",
                   "FANCL","PALB2","RAD51B","RAD51C","RAD51D","RAD54L")
#' @rdname ctdna_rules_library
#' @export
rule_genes_HRR14  <- function() rule(Gene = .HRR14_GENES)

# TSG (per the TSG scheme spec -- extended literature list)
# v0.35.0: BRCA1 and BRCA2 removed (they live in HRR14 / HRR14 instead).
.TSG_GENES <- c("TP53","RB1","PTEN","APC","NF1","NF2","VHL","CDKN2A","CDKN2B",
                 "STK11","SMAD4","CDH1","WT1","BAP1","ARID1A",
                 "ARID1B","ARID2","SETD2","KMT2D","KMT2C","TSC1","TSC2","FBXW7")
#' @rdname ctdna_rules_library
#' @export
rule_genes_TSG   <- function() rule(Gene = .TSG_GENES)

# RTK / MAPK pathway
.RTK_GENES <- c("EGFR","ERBB2","ERBB3","ERBB4","MET","ALK","ROS1","RET","KIT",
                 "PDGFRA","PDGFRB","FGFR1","FGFR2","FGFR3","FGFR4","NTRK1",
                 "NTRK2","NTRK3","KRAS","NRAS","HRAS","BRAF","RAF1","MAP2K1",
                 "MAP2K2","NF1","PTPN11")
#' @rdname ctdna_rules_library
#' @export
rule_genes_RTK <- function() rule(Gene = .RTK_GENES)

# Cell cycle
.CELL_CYCLE_GENES <- c("CDKN2A","CDKN2B","CDKN1A","CDKN1B","CDK4","CDK6","CCND1",
                        "CCND2","CCND3","CCNE1","RB1","E2F1","E2F3","MDM2","MDM4",
                        "TP53","TP73")
#' @rdname ctdna_rules_library
#' @export
rule_genes_Cell_Cycle <- function() rule(Gene = .CELL_CYCLE_GENES)

# TP53 pathway
.TP53_PATHWAY_GENES <- c("TP53","TP63","TP73","MDM2","MDM4","ATM","ATR","CHEK1","CHEK2")
#' @rdname ctdna_rules_library
#' @export
rule_genes_TP53_pathway <- function() rule(Gene = .TP53_PATHWAY_GENES)

# MMR
.MMR_GENES <- c("MLH1","MSH2","MSH6","PMS2","MSH3","MLH3","POLE","POLD1",
                 "EXO1","RFC1","PCNA")
#' @rdname ctdna_rules_library
#' @export
rule_genes_MMR <- function() rule(Gene = .MMR_GENES)

# PI3K
.PI3K_GENES <- c("PIK3CA","PIK3CB","PIK3R1","PIK3R2","PTEN","AKT1","AKT2","AKT3",
                  "MTOR","TSC1","TSC2","RHEB","RPTOR","RICTOR","STK11","INPP4B",
                  "PHLPP1","PHLPP2")
#' @rdname ctdna_rules_library
#' @export
rule_genes_PI3K <- function() rule(Gene = .PI3K_GENES)

#' Look up one or more built-in gene sets
#'
#' Returns the gene-symbol vector(s) for a built-in gene set. Same
#' registry the oncoprint's character-vector \code{gene_sets} argument
#' resolves against, exposed as a small helper so you can mix built-in
#' sets with custom ones without reaching into the package internals.
#'
#' @param name Character scalar or vector of built-in gene-set names.
#'   Built-ins (case-sensitive): \code{"HRR14"}, \code{"TSG"},
#'   \code{"RTK"}, \code{"Cell_Cycle"}, \code{"TP53_pathway"},
#'   \code{"MMR"}, \code{"PI3K"}.
#' @return A character vector of gene symbols when \code{name} has
#'   length 1; a NAMED list with one element per requested set when
#'   \code{name} has length > 1.
#' @examples
#' ctdna_gene_set("HRR14")
#' ctdna_gene_set(c("HRR14","TSG"))
#'
#' # Mix a built-in with a custom set inside ctdna_oncoprint:
#' ## ctdna_oncoprint(df,
#' ##   gene_sets = list(TSG  = ctdna_gene_set("TSG"),
#' ##                    Mine = c("TP53","RB1","PTEN")))
#' @rdname ctdna_rules_library
#' @export
ctdna_gene_set <- function(name) {
  if (!is.character(name) || !length(name))
    stop("ctdna_gene_set: `name` must be a non-empty character vector.",
         call. = FALSE)
  registry <- list(
    HRR14        = .HRR14_GENES,
    TSG          = .TSG_GENES,
    RTK          = .RTK_GENES,
    Cell_Cycle   = .CELL_CYCLE_GENES,
    TP53_pathway = .TP53_PATHWAY_GENES,
    MMR          = .MMR_GENES,
    PI3K         = .PI3K_GENES
  )
  unknown <- setdiff(name, names(registry))
  if (length(unknown))
    stop("ctdna_gene_set: unknown gene-set name(s): ",
         paste(unknown, collapse = ", "),
         ". Built-ins: ", paste(names(registry), collapse = ", "), ".",
         call. = FALSE)
  if (length(name) == 1L) registry[[name]] else registry[name]
}


# -----------------------------------------------------------------------------
# Indication rules (use the cancertype_dictionary from ctdna_opts)
# -----------------------------------------------------------------------------

# Returns the canonical-name vector for a given indication, expanded
# via the user-configurable dictionary in ctdna_opts.
.indication_values <- function(canonical) {
  dict <- tryCatch(.o("cancertype_dictionary"), error = function(e) NULL)
  if (is.null(dict) || !canonical %in% names(dict))
    return(canonical)
  vals <- dict[[canonical]]
  if (is.null(vals) || (length(vals) == 1 && is.na(vals)))
    return(canonical)
  unique(c(canonical, as.character(vals)))
}

#' @rdname ctdna_rules_library
#' @export
rule_ind_NSCLC <- function() rule(Cancertype = .indication_values("NSCLC"))
#' @rdname ctdna_rules_library
#' @export
rule_ind_BRCA  <- function() rule(Cancertype = .indication_values("BRCA"))
#' @rdname ctdna_rules_library
#' @export
rule_ind_SCLC  <- function() rule(Cancertype = .indication_values("SCLC"))
#' @rdname ctdna_rules_library
#' @export
rule_ind_HNSCC <- function() rule(Cancertype = .indication_values("HNSCC"))
#' @rdname ctdna_rules_library
#' @export
rule_ind_CRC   <- function() rule(Cancertype = .indication_values("CRC"))
#' @rdname ctdna_rules_library
#' @export
rule_ind_mCRPC <- function() rule(Cancertype = .indication_values("mCRPC"))
#' @rdname ctdna_rules_library
#' @export
rule_ind_GBM   <- function() rule(Cancertype = .indication_values("GBM"))


# -----------------------------------------------------------------------------
# Curated exclusion rules
# -----------------------------------------------------------------------------

# BRCA2 p.K3326* -- knowledge-based benign exclusion per the HRR spec
#' @rdname ctdna_rules_library
#' @export
rule_BRCA2_K3326 <- function()
  rule(Gene = "BRCA2", Mut_aa = c("p.K3326*","K3326*","p.Lys3326*"))


# =============================================================================
# Built-in schemes
# =============================================================================

#' HRR biomarker scheme (14 HRR genes)
#'
#' Built-in scheme implementing the HRR biomarker rules on the HRR14
#' gene list (BRCA1, BRCA2, ATM, BARD1, BRIP1, CDK12, CHEK1, CHEK2,
#' FANCL, PALB2, RAD51B, RAD51C, RAD51D, RAD54L). Drops ClinVar
#' benign / likely-benign. Passes: germline or somatic ClinVar P/LP,
#' OR somatic truncating LoF, OR somatic homozygous deletion, OR
#' somatic deleterious intra-genic LGR (deletion, tandem-duplication,
#' inversion).
#'
#' \strong{Side effect.} Calling \code{HRR14()} both returns the
#' scheme object AND registers it in the session's scheme catalog
#' under the name \code{"HRR14"}, so downstream code can refer to it
#' by name via \code{ctdna_variant_filter(prep, "HRR14", ...)}.
#'
#' @return A \code{ctdna_filtering_scheme} object registered as \code{HRR14}.
#' @seealso \code{\link{TSG_Tier1}}, \code{\link{TSG_Tier2}},
#'   \code{\link{TSG_Tier3}}, \code{\link{scheme_basic}},
#'   \code{\link{ctdna_variant_filter}}.
#' @examples
#' sch <- HRR14()
#' print(sch)
#'
#' sim  <- ctdna_make_mock_study(n_patients = 30, seed = 1)
#' prep <- ctdna_prepare(infinity_report = sim$infinity_report,
#'                       adam = list(adsl = sim$clinical), verbose = FALSE)
#' prep <- ctdna_sample_qc(prep, verbose = FALSE)
#' prep <- ctdna_variant_filter(prep, "HRR14",
#'                              apply = TRUE, explain = TRUE, verbose = FALSE)
#' print(head(prep$filter_explanation$HRR14$filtering_criteria))
#' @export
HRR14 <- function() {
  # Deleterious intra-GENIC LGR: same gene on both breakpoints (NEW LGR path).
  lgr_intra_gene <- allOf(
    rule_is_LGR(),
    rule_deleterious(),
    anyOf(rule_lgr_deletion(),
          rule_lgr_tandem_duplication(),
          rule_lgr_inversion()),
    rule(Gene = list(op = "==col", value = "Fusion_gene_b")))

  create_filtering_scheme(
    rule_genes_HRR14(),
    not(rule_clinvar_benign()),
    not(rule_ch()),                     # CHIP always excluded
    branches = list(
      clinvar_PLP                = rule_clinvar_path(),
      truncating_rare_notK3326   = allOf(rule_truncating(), rule_rare_001(),
                                         not(rule_BRCA2_K3326())),
      intra_gene_LGR             = allOf(rule_somatic(), lgr_intra_gene),
      homozyg_del                = allOf(rule_is_CNV(),
                                         rule_cnv_homozyg_del())),
    name        = "HRR14",
    description = paste0(
      "HRR biomarker on the 14 HRR genes. Drops ClinVar benign and CHIP. ",
      "Passes: germline or somatic ClinVar P/LP, OR truncating + rare(<1%) ",
      "not-BRCA2-K3326X, OR somatic deleterious intra-genic LGR ",
      "(deletion, tandem-duplication, inversion; Gene == Fusion_gene_b), ",
      "OR homozygous deletion."),
    category    = "gene_set",
    overwrite   = TRUE)
}

#' TSG biomarker scheme -- Tier 1 (most permissive)
#'
#' Built-in scheme implementing Tier 1 TSG rules on TP53 / RB1 / PTEN.
#' Drops ClinVar benign / likely-benign. Passes: germline P/LP, OR
#' somatic non-synonymous (missense / inframe / LoF, plus rare P/LP
#' synonymous carve-in), OR somatic LOH or homozygous deletion, OR
#' somatic deleterious LGR (any subtype), OR somatic fusion.
#'
#' \strong{Side effect.} Calling \code{TSG_Tier1()} both returns the
#' scheme object AND registers it in the session's scheme catalog
#' under the name \code{"TSG_Tier1"}, so downstream code can refer to
#' it by name via \code{ctdna_variant_filter(prep, "TSG_Tier1", ...)}
#' without keeping the returned object in scope. Re-invoking the
#' constructor is idempotent.
#'
#' Compared with \code{\link{TSG_Tier2}}: Tier 1 keeps somatic
#' missense/inframe without requiring ClinVar-P/LP support and keeps
#' fusions.
#'
#' @return A \code{ctdna_filtering_scheme} registered as \code{TSG_Tier1}.
#' @seealso \code{\link{TSG_Tier2}}, \code{\link{TSG_Tier3}},
#'   \code{\link{HRR14}}, \code{\link{scheme_basic}},
#'   \code{\link{ctdna_variant_filter}}.
#' @examples
#' # Constructor returns and registers the scheme
#' sch <- TSG_Tier1()
#' print(sch)
#'
#' # Typical use: pass the name to ctdna_variant_filter().
#' sim  <- ctdna_make_mock_study(n_patients = 30, seed = 1)
#' prep <- ctdna_prepare(infinity_report = sim$infinity_report,
#'                       adam = list(adsl = sim$clinical), verbose = FALSE)
#' prep <- ctdna_sample_qc(prep, verbose = FALSE)
#' prep <- ctdna_variant_filter(prep, "TSG_Tier1",
#'                              apply = TRUE, explain = TRUE, verbose = FALSE)
#' print(head(prep$filter_explanation$TSG_Tier1$filtering_criteria))
#' @export
TSG_Tier1 <- function() {
  # Non-synonymous impactful short variants. Rare (<1%) ClinVar P/LP
  # synonymous is carved back in per the biomarker doc.
  nonsyn_impactful <- anyOf(
    rule(Molecular_consequence = c(
        "missense",
        "inframe_indel", "inframe_insertion",
        "inframe_deletion", "inframe_duplication",
        "nonsense", "frameshift", "stop_lost", "start_lost",
        "splice_acceptor", "splice_donor")),
    allOf(rule(Molecular_consequence = "synonymous"),
          rule_clinvar_path(),
          rule_rare_001()))

  create_filtering_scheme(
    rule(Gene = c("TP53","RB1","PTEN")),
    not(rule_clinvar_benign()),
    not(rule_ch()),                     # CHIP always excluded
    branches = list(
      germline_PLP        = allOf(rule_germline(), rule_clinvar_path()),
      somatic_nonsyn      = allOf(rule_somatic(), nonsyn_impactful),
      LOH_or_homozyg_del  = allOf(rule_somatic(), anyOf(rule_cnv_loh_del(),
                                                        rule_cnv_homozyg_del())),
      deleterious_LGR     = allOf(rule_somatic(), rule_is_LGR_deleterious()),
      fusion              = allOf(rule_somatic(), rule_is_Fusion())),
    name        = "TSG_Tier1",
    description = paste0(
      "TSG Tier 1 (most permissive) on TP53/RB1/PTEN. Drops ClinVar ",
      "benign and CHIP. Passes: germline P/LP, OR somatic non-synonymous, ",
      "OR somatic LOH/homozygous deletion, OR somatic deleterious LGR ",
      "(any subtype), OR somatic fusion."),
    category    = "gene_set",
    overwrite   = TRUE)
}


#' TSG biomarker scheme -- Tier 2 (intermediate)
#'
#' Built-in scheme implementing Tier 2 TSG rules on TP53 / RB1 / PTEN.
#' Stricter than Tier 1: somatic short variants must be LoF or
#' ClinVar-P/LP missense/inframe; deleterious LGR must be intra-chromosomal
#' (deletion, tandem-duplication, inversion).
#'
#' \strong{Side effect.} Calling \code{TSG_Tier2()} both returns the
#' scheme object AND registers it in the session's scheme catalog
#' under the name \code{"TSG_Tier2"}, so downstream code can refer to
#' it by name via \code{ctdna_variant_filter(prep, "TSG_Tier2", ...)}.
#'
#' @return A \code{ctdna_filtering_scheme} registered as \code{TSG_Tier2}.
#' @seealso \code{\link{TSG_Tier1}}, \code{\link{TSG_Tier3}},
#'   \code{\link{HRR14}}, \code{\link{scheme_basic}},
#'   \code{\link{ctdna_variant_filter}}.
#' @examples
#' sch <- TSG_Tier2()
#' print(sch)
#'
#' sim  <- ctdna_make_mock_study(n_patients = 30, seed = 1)
#' prep <- ctdna_prepare(infinity_report = sim$infinity_report,
#'                       adam = list(adsl = sim$clinical), verbose = FALSE)
#' prep <- ctdna_sample_qc(prep, verbose = FALSE)
#' prep <- ctdna_variant_filter(prep, "TSG_Tier2",
#'                              apply = TRUE, explain = TRUE, verbose = FALSE)
#' print(head(prep$filter_explanation$TSG_Tier2$filtering_criteria))
#' @export
TSG_Tier2 <- function() {
  # Non-synonymous impactful (missense + inframe forms). Rare P/LP synonymous
  # carve-in -- carried over from the old TSG_Tier2.
  nonsynonymous <- rule(Molecular_consequence =
    c("missense","inframe_indel","inframe_insertion",
      "inframe_deletion","inframe_duplication"))
  rare_path_synonymous <- allOf(
    rule(Molecular_consequence = "synonymous"),
    rule_clinvar_path(),
    rule_rare_001())

  # Deleterious intra-CHROMOSOMAL LGR (Tier 2 doc spec): all three subtypes
  # require same chromosome.
  lgr_intra_chrom <- allOf(
    rule_is_LGR(),
    rule_deleterious(),
    anyOf(rule_lgr_deletion(),
          rule_lgr_tandem_duplication(),
          rule_lgr_inversion()),
    rule(Chromosome = list(op = "==col", value = "Fusion_chrom_b")))

  create_filtering_scheme(
    rule(Gene = c("TP53","RB1","PTEN")),
    not(rule_clinvar_benign()),
    not(rule_ch()),                     # CHIP always excluded
    branches = list(
      germline_PLP           = allOf(rule_germline(), rule_clinvar_path()),
      somatic_LoF            = allOf(rule_somatic(), rule_truncating()),
      somatic_PLP_missense   = allOf(rule_somatic(), rule_clinvar_path(),
                                     anyOf(nonsynonymous, rare_path_synonymous)),
      LOH_or_homozyg_del     = allOf(rule_somatic(), rule_is_CNV(),
                                     anyOf(rule_cnv_homozyg_del(),
                                           rule_cnv_loh_del())),
      intra_chrom_LGR        = allOf(rule_somatic(), lgr_intra_chrom)),
    name        = "TSG_Tier2",
    description = paste0(
      "TSG Tier 2 (intermediate) on TP53/RB1/PTEN. Drops ClinVar benign ",
      "and CHIP. Passes: germline P/LP, OR somatic LoF, OR somatic ",
      "ClinVar-P/LP + (missense|inframe|rare-P/LP-synonymous), OR somatic ",
      "LOH/homozygous deletion, OR somatic deleterious intra-chromosomal ",
      "LGR (deletion, tandem-duplication, inversion; all Chromosome == ",
      "Fusion_chrom_b)."),
    category    = "gene_set",
    overwrite   = TRUE)
}


#' TSG biomarker scheme -- Tier 3 (most restrictive)
#'
#' Most-restrictive TSG scheme on TP53 / RB1 / PTEN. Drops ClinVar
#' benign and CHIP. Passes: germline P/LP, OR somatic LoF, OR somatic
#' ClinVar-P/LP + (missense/inframe), OR somatic LOH/homozygous
#' deletion. No LGR path, no fusion path, no rare-P/LP-synonymous
#' carve-in.
#'
#' Compared with \code{\link{TSG_Tier2}}: drops the LGR path and the
#' rare-P/LP-synonymous carve-in; ClinVar-P/LP missense/inframe path
#' is kept as-is.
#'
#' \strong{Side effect.} Calling \code{TSG_Tier3()} both returns the
#' scheme object AND registers it in the session's scheme catalog
#' under the name \code{"TSG_Tier3"}, so downstream code can refer to
#' it by name via \code{ctdna_variant_filter(prep, "TSG_Tier3", ...)}.
#'
#' @return A \code{ctdna_filtering_scheme} registered as \code{TSG_Tier3}.
#' @seealso \code{\link{TSG_Tier1}}, \code{\link{TSG_Tier2}},
#'   \code{\link{HRR14}}, \code{\link{scheme_basic}},
#'   \code{\link{ctdna_variant_filter}}.
#' @examples
#' sch <- TSG_Tier3()
#' print(sch)
#'
#' sim  <- ctdna_make_mock_study(n_patients = 30, seed = 1)
#' prep <- ctdna_prepare(infinity_report = sim$infinity_report,
#'                       adam = list(adsl = sim$clinical), verbose = FALSE)
#' prep <- ctdna_sample_qc(prep, verbose = FALSE)
#' prep <- ctdna_variant_filter(prep, "TSG_Tier3",
#'                              apply = TRUE, explain = TRUE, verbose = FALSE)
#' print(head(prep$filter_explanation$TSG_Tier3$filtering_criteria))
#' @export
TSG_Tier3 <- function() {
  # ClinVar P/LP + missense-or-inframe only (no rare-syn carve-in).
  missense_or_inframe_path <- allOf(
    rule_clinvar_path(),
    rule(Molecular_consequence = c(
        "missense",
        "inframe_indel", "inframe_insertion",
        "inframe_deletion", "inframe_duplication")))

  create_filtering_scheme(
    rule(Gene = c("TP53","RB1","PTEN")),
    not(rule_clinvar_benign()),
    not(rule_ch()),                     # CHIP always excluded
    branches = list(
      germline_PLP           = allOf(rule_germline(), rule_clinvar_path()),
      somatic_LoF            = allOf(rule_somatic(), rule_truncating()),
      somatic_PLP_missense   = allOf(rule_somatic(), missense_or_inframe_path),
      LOH_or_homozyg_del     = allOf(rule_somatic(), rule_is_CNV(),
                                     anyOf(rule_cnv_homozyg_del(),
                                           rule_cnv_loh_del()))),
    name        = "TSG_Tier3",
    description = paste0(
      "TSG Tier 3 (most restrictive) on TP53/RB1/PTEN. Drops ClinVar ",
      "benign and CHIP. Passes: germline P/LP, OR somatic LoF, OR somatic ",
      "ClinVar-P/LP + missense/inframe, OR somatic LOH/homozygous deletion. ",
      "No LGR, no fusion, no rare-P/LP-synonymous carve-in."),
    category    = "gene_set",
    overwrite   = TRUE)
}

# Returns a named list mapping rule-name -> function that produces the rule.
# Used by inspection tools.
.ctdna_rule_library <- function() {
  list(
    # Variant_type
    rule_is_SNV         = rule_is_SNV,
    rule_is_Indel       = rule_is_Indel,
    rule_is_CNV         = rule_is_CNV,
    rule_is_Fusion      = rule_is_Fusion,
    rule_is_LGR         = rule_is_LGR,
    rule_is_LGR_deleterious = rule_is_LGR_deleterious,
    rule_amp_thresh     = rule_amp_thresh,
    # CNV_type
    rule_cnv_focal_amp     = rule_cnv_focal_amp,
    rule_cnv_aneuploid_amp = rule_cnv_aneuploid_amp,
    rule_cnv_amp           = rule_cnv_amp,
    rule_cnv_homozyg_del   = rule_cnv_homozyg_del,
    rule_cnv_loh_del       = rule_cnv_loh_del,
    # Splice_effect
    rule_splice_acceptor = rule_splice_acceptor,
    rule_splice_donor    = rule_splice_donor,
    rule_splice_region   = rule_splice_region,
    # Sample_status
    rule_sample_success  = rule_sample_success,
    rule_sample_fail     = rule_sample_fail,
    # Somatic_status
    rule_germline        = rule_germline,
    rule_somatic         = rule_somatic,
    rule_ch              = rule_ch,
    # Molecular_consequence
    rule_missense        = rule_missense,
    rule_synonymous      = rule_synonymous,
    rule_nonsense        = rule_nonsense,
    rule_frameshift      = rule_frameshift,
    rule_start_lost      = rule_start_lost,
    rule_promoter        = rule_promoter,
    rule_inframe_indel   = rule_inframe_indel,
    rule_splice_event_mc = rule_splice_event_mc,
    rule_truncating      = rule_truncating,
    # ClinVar
    rule_clinvar_benign  = rule_clinvar_benign,
    rule_clinvar_path    = rule_clinvar_path,
    rule_clinvar_vus     = rule_clinvar_vus,
    rule_clinvar_conflict= rule_clinvar_conflict,
    # Functional_impact
    rule_deleterious     = rule_deleterious,
    rule_reversion       = rule_reversion,
    # Mutant_allele_status
    rule_biallelic       = rule_biallelic,
    # gnomAD
    rule_rare_001        = rule_rare_001,
    rule_rare_0001       = rule_rare_0001,
    rule_common_001      = rule_common_001,
    # Gene sets
    rule_genes_HRR14        = rule_genes_HRR14,
    rule_genes_TSG          = rule_genes_TSG,
    rule_genes_RTK          = rule_genes_RTK,
    rule_genes_Cell_Cycle   = rule_genes_Cell_Cycle,
    rule_genes_TP53_pathway = rule_genes_TP53_pathway,
    rule_genes_MMR          = rule_genes_MMR,
    rule_genes_PI3K         = rule_genes_PI3K,
    # Indications
    rule_ind_NSCLC = rule_ind_NSCLC,
    rule_ind_BRCA  = rule_ind_BRCA,
    rule_ind_SCLC  = rule_ind_SCLC,
    rule_ind_HNSCC = rule_ind_HNSCC,
    rule_ind_CRC   = rule_ind_CRC,
    rule_ind_mCRPC = rule_ind_mCRPC,
    rule_ind_GBM   = rule_ind_GBM,
    # Curated exclusions
    rule_BRCA2_K3326 = rule_BRCA2_K3326)
}

# Returns the built-in schemes as a named list.
.ctdna_builtin_schemes <- function() {
  list(
    TSG_Tier1    = TSG_Tier1(),
    TSG_Tier2    = TSG_Tier2(),
    TSG_Tier3    = TSG_Tier3(),
    HRR14        = HRR14(),
    scheme_basic = scheme_basic())
}


#' Basic ruleset for retaining impactful Guardant Health Infinity alterations
#'
#' Built-in v1.0.0 scheme translating the canonical GH Infinity
#' \code{func_retainImpactfulAlterations} logic into the composable
#' rule system. Used as the default filter when no scheme is named, and
#' as the catch-all for genes that fall outside any gene-set scheme
#' (e.g. when \code{filter_scheme = "HRR14"}, non-HRR genes get
#' \code{scheme_basic} applied to them).
#'
#' \strong{Per-variant-type logic (matches the GH function):}
#' \itemize{
#'   \item \strong{SNV / Indel:} somatic, NOT synonymous, ClinVar NOT in
#'     \{Likely_benign, Benign/Likely_benign, Benign,
#'     Uncertain_significance\}, VAF >= \code{somatic_vaf_threshold}.
#'   \item \strong{SNV / Indel germline (optional):} for genes in
#'     \code{germline_genes}, also keep germline variants that are NOT
#'     synonymous AND ClinVar in \{Pathogenic, Pathogenic/Likely_pathogenic,
#'     Likely_pathogenic\}.
#'   \item \strong{CNV:} somatic; amplification / aneuploid_amplification /
#'     focal_amplification require \code{Copy_number >=
#'     cna_copynumber_threshold}; all other CNV types (homozyg_del,
#'     loh_del, etc.) pass without copy-number threshold.
#'   \item \strong{Fusion:} somatic, VAF >= \code{somatic_vaf_threshold}.
#'   \item \strong{LGR:} somatic, VAF >= \code{somatic_vaf_threshold},
#'     \code{Functional_impact == "deleterious"}.
#' }
#'
#' @param somatic_vaf_threshold Minimum VAF percentage for somatic
#'   SNV/Indel/Fusion/LGR variants. Default 0 (no threshold).
#' @param cna_copynumber_threshold Minimum copy number for amplification-
#'   class CNVs. Default 4.
#' @param germline_genes Character vector of gene symbols for which germline
#'   pathogenic/likely-pathogenic variants are kept. Default BRCA1/BRCA2.
#'   Pass NULL to disable germline handling entirely (somatic-only).
#' @return A \code{ctdna_filtering_scheme} object.
#' @seealso \code{\link{HRR14}}, \code{\link{TSG_Tier1}},
#'   \code{\link{TSG_Tier2}}, \code{\link{TSG_Tier3}},
#'   \code{\link{ctdna_create_scheme}},
#'   \code{\link{ctdna_variant_filter}}.
#' @examples
#' # Default: permissive somatic + BRCA1/BRCA2 germline P/LP.
#' sch <- scheme_basic()
#' print(sch)
#'
#' # Tune: require somatic VAF >= 1%, only call amplifications with CN >= 6,
#' # and expand germline handling to a wider HRR list.
#' sch_strict <- scheme_basic(
#'   somatic_vaf_threshold    = 1.0,
#'   cna_copynumber_threshold = 6,
#'   germline_genes           = c("BRCA1","BRCA2","PALB2","ATM","CHEK2"))
#' print(sch_strict)
#'
#' # Somatic-only variant: pass germline_genes = NULL.
#' sch_som <- scheme_basic(germline_genes = NULL)
#' print(sch_som)
#' @rdname ctdna_rules_library
#' @export
scheme_basic <- function(somatic_vaf_threshold   = 0,
                          cna_copynumber_threshold = 4,
                          germline_genes           = c("BRCA1","BRCA2")) {

  # ---- SNV / Indel: somatic, impactful ----------------------------------
  snv_indel_somatic <- allOf(
    rule(Variant_type = c("SNV","Indel")),
    rule_somatic(),
    not(rule_synonymous()),
    not(rule(ClinVar = c("Likely_benign","Benign/Likely_benign","Benign",
                          "Uncertain_significance"))),
    rule(VAF_percentage = list(op = ">=", value = somatic_vaf_threshold)))

  # ---- SNV / Indel germline P/LP (optional, on configured genes) -------
  snv_indel_germline <- if (!is.null(germline_genes) && length(germline_genes) > 0) {
    allOf(
      rule(Variant_type = c("SNV","Indel")),
      rule_germline(),
      rule(Gene = germline_genes),
      not(rule_synonymous()),
      rule(ClinVar = c("Pathogenic","Pathogenic/Likely_pathogenic",
                        "Likely_pathogenic")))
  } else NULL

  # ---- CNV: somatic; amp-class needs copy-number threshold -------------
  amp_classes <- c("amplification","aneuploid_amplification","focal_amplification")
  cnv_branch <- allOf(
    rule_is_CNV(),
    rule_somatic(),
    anyOf(
      # Non-amp CNV (homozyg_del, loh_del, etc.) passes without threshold
      not(rule(CNV_type = amp_classes)),
      # Amp-class needs Copy_number >= threshold
      allOf(rule(CNV_type = amp_classes),
              rule(Copy_number = list(op = ">=",
                                        value = cna_copynumber_threshold)))))

  # ---- Fusion: somatic, VAF threshold ----------------------------------
  fusion_branch <- allOf(
    rule_is_Fusion(),
    rule_somatic(),
    rule(VAF_percentage = list(op = ">=", value = somatic_vaf_threshold)))

  # ---- LGR: somatic, VAF threshold, deleterious ------------------------
  lgr_branch <- allOf(
    rule_is_LGR(),
    rule_somatic(),
    rule(VAF_percentage = list(op = ">=", value = somatic_vaf_threshold)),
    rule_deleterious())

  # ---- Top-level: any of the variant-type branches ---------------------
  branches <- list(snv_indel_somatic)
  if (!is.null(snv_indel_germline)) branches <- c(branches, list(snv_indel_germline))
  branches <- c(branches, list(cnv_branch, fusion_branch, lgr_branch))

  root <- do.call(anyOf, branches)

  create_filtering_scheme(
    root,
    name        = "scheme_basic",
    description = paste0("Basic ruleset for impactful GH Infinity alterations ",
                          "(translated from func_retainImpactfulAlterations_GH_Infinity). ",
                          "Per-variant-type logic: SNV/Indel somatic+impactful (+optional ",
                          "germline P/LP on configured genes); CNV somatic with copy-number ",
                          "threshold for amplification classes; Fusion + LGR somatic with ",
                          "VAF threshold (LGR also requires deleterious)."),
    category    = "user",
    scope       = NA,  # scheme_basic applies to ALL genes (no scope)
    overwrite   = TRUE)
}

