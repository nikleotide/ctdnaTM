# =============================================================================
# v1.0.0 — Pre-shipped rule library + built-in schemes
# =============================================================================
#
# This file defines the pre-built ctdna_rule atoms users can plug into
# their schemes without writing them from scratch, plus the two
# canonical built-in schemes (HRR and TSG) defined per the published
# specifications.
#
# Naming convention:
#   rule_<column>_<value>  — single-column constraint on a canonical value
#   rule_genes_<name>      — gene-set membership rule
#   rule_ind_<name>        — indication rule (uses ctdna_opts cancertype_dictionary)
#   rule_<composite>       — multi-rule composite (e.g. rule_truncating)
#   scheme_<name>          — full filtering scheme
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
#' \code{\link{rule_and}}, \code{\link{rule_or}}, \code{\link{rule_not}}
#' or the convenience \code{\link{anyOf}} helper to build filter
#' schemes.
#'
#' Naming convention:
#' \itemize{
#'   \item \code{rule_is_<TYPE>()} — \code{Variant_type} predicate
#'     (SNV, Indel, CNV, Fusion, LGR).
#'   \item \code{rule_lgr_<value>()} — \code{LGR_subtype} predicate
#'     (deletion, tandem duplication, inversion).
#'   \item \code{rule_cnv_<value>()} — \code{CNV_type} predicate
#'     (focal_amp, aneuploid_amp, amp, homozyg_del, loh_del).
#'   \item \code{rule_<molecular_consequence>()} — \code{Molecular_consequence}
#'     predicate (missense, nonsense, frameshift, splice_donor,
#'     splice_acceptor, splice_region, start_lost, stop_lost,
#'     synonymous, promoter, inframe_indel, ...).
#'   \item \code{rule_sample_<status>()} — \code{Sample_status} predicate
#'     (success, fail).
#'   \item \code{rule_<origin>()} — \code{Variant_origin} predicate
#'     (germline, somatic, ch).
#'   \item \code{rule_clinvar_<class>()} — \code{ClinVar_class}
#'     predicate (benign, path).
#'   \item Composites: \code{rule_truncating()},
#'     \code{rule_is_LGR_deleterious()},
#'     \code{rule_splice_event_mc()}.
#' }
#'
#' All rules are lazy: each constructor returns the rule fresh, so
#' changes to \code{ctdna_opts()} (e.g. \code{col_clinvar}) are picked
#' up at use time.
#'
#' @return A \code{ctdna_rule} object suitable for
#'   \code{\link{create_filtering_scheme}}.
#' @examples
#' # Two ways to express "SNV missense not in ClinVar benign":
#' rule_and(rule_is_SNV(), rule_missense(),
#'          rule_not(rule_clinvar_benign()))
#'
#' # Composite predicates (built from atoms):
#' rule_truncating()         # nonsense + frameshift + stop_lost + start_lost
#' rule_splice_event_mc()    # splice_donor / splice_acceptor / splice_region
#' rule_is_LGR_deleterious() # LGR with a clinically meaningful subtype
#'
#' # Use in a scheme:
#' create_filtering_scheme(
#'   rule_is_SNV(),
#'   rule_or(rule_missense(), rule_truncating()),
#'   rule_not(rule_germline()),
#'   name = "snv_coding_somatic")
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

# Truncating / null variants — the canonical LoF set.
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

# TSG (per the TSG scheme spec — extended literature list)
# v0.35.0: BRCA1 and BRCA2 removed (they live in scheme_HRR / HRR14 instead).
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
#' # Mix a built-in with a custom set inside ctdna_plot_oncoprint:
#' ## ctdna_plot_oncoprint(df,
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

# BRCA2 p.K3326* — knowledge-based benign exclusion per the HRR spec
#' @rdname ctdna_rules_library
#' @export
rule_BRCA2_K3326 <- function()
  rule(Gene = "BRCA2", Mut_aa = c("p.K3326*","K3326*","p.Lys3326*"))


# =============================================================================
# Built-in schemes
# =============================================================================

#' HRR (Homologous Recombination Repair deficiency) scheme
#'
#' Built-in v1.0.0 scheme for identifying HRRm / BRCAm classifications
#' from short variants, structural variants, and copy-number calls.
#' Implements the published 8-step decision logic via the composable
#' rule system:
#' \itemize{
#'   \item Restrict to the HRR14 gene list.
#'   \item Drop benign-annotated calls.
#'   \item Drop CHIP variants.
#'   \item Pass if ClinVar pathogenic OR (truncating AND rare AND
#'     not a known knowledge-based exclusion) OR LGR (deleterious by
#'     definition in this gene set) OR homozygous deletion.
#' }
#' Population-frequency rule requires \code{ctdna_annotate_population_freq()}
#' to have been run on the input frame; otherwise the rule silently
#' evaluates to FALSE for the affected path.
#'
#' @return A \code{ctdna_filtering_scheme} object.
#' @seealso \code{\link{scheme_TSG}}, \code{\link{ctdna_filter_apply}}.
#' @examples
#' \dontrun{
#' df_filt <- ctdna_filter_apply(df, filter_scheme = "scheme_HRR")
#' # or directly:
#' mask <- ctdnaTM:::.eval_rule(scheme_HRR(), df)
#' }
#' @rdname ctdna_rules_library
#' @export
scheme_HRR <- function() {
  create_filtering_scheme(
    rule_genes_HRR14(),
    not(rule_clinvar_benign()),
    not(rule_ch()),
    anyOf(
      rule_clinvar_path(),
      allOf(rule_truncating(), rule_rare_001(), not(rule_BRCA2_K3326())),
      rule_is_LGR_deleterious(),                                 # v1.0.2: deleterious LGR only
      allOf(rule_is_CNV(), rule_cnv_homozyg_del())),
    name        = "scheme_HRR",
    description = paste0("HRR / BRCAm classification (HRR14 gene set; ClinVar ",
                          "pathogenic OR truncating/rare OR deleterious LGR OR ",
                          "homozygous deletion). LGR kept only when annotated deleterious."),
    category    = "gene_set",
    overwrite   = TRUE)  # built-in constructor; safe to re-register on repeat calls
}

#' TSG (Tumor Suppressor Gene) scheme
#'
#' Built-in v1.0.0 scheme for identifying TSG-deleterious classifications
#' from short variants and structural alterations on a curated
#' tumor-suppressor gene list. Implements the published TSG decision logic:
#' \itemize{
#'   \item Restrict to the TSG gene list (BRCA1/BRCA2 deliberately
#'     excluded -- they belong to \code{scheme_HRR}).
#'   \item Drop benign-annotated calls.
#'   \item Drop CHIP variants.
#'   \item Pass if (germline AND ClinVar pathogenic) OR
#'     (somatic AND \{truncating OR missense/nonsynonymous/inframe_indel\})
#'     OR (somatic AND deleterious LGR) OR (somatic AND
#'     \{homozygous deletion OR LoH deletion\}).
#' }
#' \strong{v0.35.0 changes:} BRCA1/BRCA2 removed from gene list;
#' LGR path now requires somatic; CNV path now requires somatic;
#' somatic short-variant path B widened beyond strict truncating to
#' include missense, nonsynonymous, and inframe_indel calls.
#'
#' @return A \code{ctdna_filtering_scheme} object.
#' @seealso \code{\link{scheme_HRR}}, \code{\link{ctdna_filter_apply}}.
#' @examples
#' \dontrun{
#' df_filt <- ctdna_filter_apply(df, filter_scheme = "scheme_TSG")
#' }
#' @rdname ctdna_rules_library
#' @export
scheme_TSG <- function() {
  # Non-synonymous (coding-impactful, non-LoF): missense + all inframe_*
  # forms. Per spec, we EXCLUDE promoter / splice_region / splice_event /
  # NA-annotated rows from this bucket (those are filtered out by the
  # absence of their annotation in this list).
  nonsynonymous <- rule(Molecular_consequence =
    c("missense","inframe_indel","inframe_insertion",
      "inframe_deletion","inframe_duplication"))

  # Rare-pathogenic synonymous carve-out (TSG point 3): a synonymous
  # variant that is ClinVar pathogenic/likely-pathogenic AND population
  # rare (<1%) is pulled into the non-synonymous bucket.
  rare_path_synonymous <- allOf(
    rule(Molecular_consequence = "synonymous"),
    rule_clinvar_path(),
    rule_rare_001())

  # Intragenic LGR (TSG point 9): deleterious-annotated structural
  # variants of three subtypes inferred from Direction_a/Direction_b.
  # Deletion + inversion: require intrachromosomal (Chromosome ==
  # Fusion_chrom_b). Tandem duplication: require intragenic (Gene ==
  # Fusion_gene_b). Direction columns absent or NA -> row drops.
  intragenic_LGR <- allOf(
    rule_is_LGR(),
    rule_deleterious(),
    anyOf(
      allOf(rule_lgr_deletion(),
              rule(Chromosome = list(op = "==col",
                                       value = "Fusion_chrom_b"))),
      allOf(rule_lgr_inversion(),
              rule(Chromosome = list(op = "==col",
                                       value = "Fusion_chrom_b"))),
      allOf(rule_lgr_tandem_duplication(),
              rule(Gene = list(op = "==col",
                                value = "Fusion_gene_b")))))

  create_filtering_scheme(
    rule_genes_TSG(),
    not(rule_clinvar_benign()),
    not(rule_ch()),
    anyOf(
      # Point 5: germline + ClinVar pathogenic
      allOf(rule_germline(), rule_clinvar_path()),
      # Point 6: somatic LoF (LoF set incl. stop_lost via rule_truncating)
      allOf(rule_somatic(), rule_truncating()),
      # Point 7: somatic + ClinVar pathogenic + missense/inframe (or
      # rare-pathogenic synonymous carved into the non-synonymous bucket)
      allOf(rule_somatic(),
              rule_clinvar_path(),
              anyOf(nonsynonymous, rare_path_synonymous)),
      # Point 8: LOH-deletion or homozygous deletion (somatic gating)
      allOf(rule_somatic(), rule_is_CNV(),
              anyOf(rule_cnv_homozyg_del(), rule_cnv_loh_del())),
      # Point 9: deleterious intragenic LGR (somatic gating)
      allOf(rule_somatic(), intragenic_LGR)),
    name        = "scheme_TSG",
    description = paste0("TSG-deleterious classification (curated TSG list, ",
                          "BRCA1/BRCA2 excluded). Pass if: germline+P/LP OR ",
                          "somatic+LoF OR somatic+P/LP+(missense|inframe|rare-P/LP-synonymous) ",
                          "OR somatic+(homozyg_del|loh_del) OR somatic+deleterious ",
                          "intragenic-LGR (deletion/inversion intrachromosomal, ",
                          "tandem-dup intragenic)."),
    category    = "gene_set",
    overwrite   = TRUE)
}


# -----------------------------------------------------------------------------
# Rule library catalog (for ctdna_filter_schemes_list)
# -----------------------------------------------------------------------------

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
    scheme_HRR   = scheme_HRR(),
    scheme_TSG   = scheme_TSG(),
    scheme_basic = scheme_basic())
}


#' Basic ruleset for retaining impactful Guardant Health Infinity alterations
#'
#' Built-in v1.0.0 scheme translating the canonical GH Infinity
#' \code{func_retainImpactfulAlterations} logic into the composable
#' rule system. Used as the default filter when no scheme is named, and
#' as the catch-all for genes that fall outside any gene-set scheme
#' (e.g. when \code{filter_scheme = "scheme_HRR"}, non-HRR genes get
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
#' @seealso \code{\link{scheme_HRR}}, \code{\link{scheme_TSG}},
#'   \code{\link{ctdna_filter_apply}}.
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

