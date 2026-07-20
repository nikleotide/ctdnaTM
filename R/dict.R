# ============================================================================
# dnaseq alteration-type dictionary  --  v0.76.0
# ----------------------------------------------------------------------------
# Single source of truth for mapping raw vendor strings coming out of dnaseq
# (and ctdna) files into ctdnaTM's canonical alteration-type universe. Runs
# automatically inside ctdna_prepare() / ctdna_prep_add() (via .prep_dnaseq)
# so every downstream function -- .oncoprint_classify() in particular -- sees
# consistent canonical values regardless of which vendor / pipeline / column
# the data came from.
#
# Public API mirrors ctdna_opts():
#   dnaseq_dict()                                       # list all
#   dnaseq_dict("variant_type")                         # get sub-dict
#   dnaseq_dict("molecular_consequence.missense")       # get one entry
#   dnaseq_dict(molecular_consequence.missense = c("MyVendor_MSs"))  # append
#   dnaseq_dict_add("cnv_type", "focal_amplification", "HighGain")
#   dnaseq_dict_reset()
#   dnaseq_dict_lookup("molecular_consequence", raw_values_vector)
# ============================================================================

.dnaseq_dict_env <- new.env(parent = emptyenv())

# Normalize a raw value into a match-key: lowercase, whitespace/underscore/
# dash all collapsed to a single "_". Everything else is preserved so we
# still distinguish e.g. "5_prime_UTR_variant" from something else.
.dnaseq_dict_norm_key <- function(x) {
  x <- as.character(x)
  x <- tolower(trimws(x))
  x <- gsub("[[:space:]_\\-]+", "_", x, perl = TRUE)
  x
}

# Shipped defaults. Covers the 52 Variant_Effect strings + the CNA Type +
# Variant_type + Functional_Class value sets provided from real data.
.dnaseq_dict_defaults <- function() {
  list(
    variant_type = list(
      SNV    = c("SNV","point_mutation","substitution","MNV"),
      Indel  = c("Indel","INDEL","insertion","deletion"),
      CNV    = c("CNV","CNA","copy_number","copy_number_variant"),
      Fusion = c("Fusion","gene_fusion","bidirectional_gene_fusion","SV"),
      LGR    = c("LGR","large_rearrangement","large_genomic_rearrangement")
    ),

    molecular_consequence = list(
      missense            = c("missense","missense_variant","MISSENSE"),
      nonsense            = c("nonsense","stop_gained","NONSENSE"),
      frameshift          = c("frameshift","frameshift_variant"),
      splice_donor        = c("splice_donor","splice_donor_variant"),
      splice_acceptor     = c("splice_acceptor","splice_acceptor_variant"),
      stop_lost           = c("stop_lost"),
      start_lost          = c("start_lost"),
      inframe_insertion   = c("inframe_insertion","in_frame_insertion"),
      inframe_deletion    = c("inframe_deletion","in_frame_deletion"),
      inframe_duplication = c("inframe_duplication","in_frame_duplication"),
      inframe_indel       = c("inframe_indel","in_frame_indel"),
      promoter            = c("promoter","promoter_variant","tert_promoter",
                              "5_prime_UTR_premature_start_codon_gain_variant"),
      silent              = c("silent","synonymous","synonymous_variant",
                              "SILENT","NONE",
                              "intron_variant","upstream_gene_variant",
                              "downstream_gene_variant",
                              "5_prime_UTR_variant","3_prime_UTR_variant",
                              "intergenic_region","non_coding_exon_variant",
                              "exon_variant","splice_region_variant",
                              "intragenic_variant")
    ),

    cnv_type = list(
      focal_amplification = c("focal_amplification","AMP","focal_amp"),
      amplification       = c("amplification","aneuploid_amplification","GAIN"),
      homozygous_deletion = c("homozygous_deletion","DEL","hom_del"),
      loh_deletion        = c("loh_deletion","LOH","loss_of_heterozygosity")
    )
  )
}

# Priority for combined consequences. Highest damage first. Used to pick
# ONE canonical class from combined strings like
# "frameshift_variant+splice_donor_variant+intron_variant".
.dnaseq_dict_priority <- function() {
  c("frameshift","nonsense","splice_donor","splice_acceptor",
    "stop_lost","start_lost",
    "inframe_indel","inframe_deletion","inframe_insertion","inframe_duplication",
    "missense","promoter","silent")
}

.dnaseq_dict_get_state <- function() {
  if (!exists("dict", envir = .dnaseq_dict_env))
    assign("dict", .dnaseq_dict_defaults(), envir = .dnaseq_dict_env)
  get("dict", envir = .dnaseq_dict_env)
}

# Build a flat normalized-alias -> canonical map for one sub-dict.
.dnaseq_dict_flat <- function(sub) {
  d <- .dnaseq_dict_get_state()[[sub]]
  out <- character(0)
  for (canonical in names(d)) {
    for (alias in d[[canonical]]) {
      out[.dnaseq_dict_norm_key(alias)] <- canonical
    }
  }
  out
}

#' Vendor -> canonical alteration-type dictionary
#'
#' @description
#' The dictionary that maps raw vendor strings (`"missense_variant"`,
#' `"stop_gained"`, `"AMP"`, `"CNA"`, `"gene_fusion"`, `"frameshift_variant+
#' splice_donor_variant+intron_variant"`, ...) into ctdnaTM's canonical
#' alteration-type universe used by every oncoprint / concordance / stats
#' function.
#'
#' Three sub-dictionaries:
#' \describe{
#'   \item{variant_type}{Normalises `Variant_type` / `Variant Type` columns.
#'     Canonical: `SNV`, `Indel`, `CNV`, `Fusion`, `LGR`.}
#'   \item{molecular_consequence}{Normalises `Variant_Effect` /
#'     `Molecular_consequence` / `Functional_Class` columns. Canonical:
#'     `missense`, `nonsense`, `frameshift`, `splice_donor`,
#'     `splice_acceptor`, `stop_lost`, `start_lost`, `inframe_*`,
#'     `promoter`, `silent` (silent covers intronic/UTR/synonymous/
#'     intergenic/non-coding effects that shouldn't appear in oncoprint).}
#'   \item{cnv_type}{Normalises `CNV_type` / `CNA_type` / `CNA Type`
#'     columns. Canonical: `focal_amplification`, `amplification`,
#'     `homozygous_deletion`, `loh_deletion`.}
#' }
#'
#' Matching is CASE-INSENSITIVE and normalises whitespace/underscore/dash,
#' so `"stop_gained"`, `"Stop Gained"`, `"STOP-GAINED"` all match. Combined
#' consequences with `+` separator (e.g. `"frameshift_variant+
#' splice_donor_variant+intron_variant"`) are split, each part looked up,
#' and the MOST-DAMAGING component wins via a fixed priority order
#' (frameshift > nonsense > splice > stop/start_lost > inframe_* > missense
#' > promoter > silent).
#'
#' Runs automatically inside \code{\link{ctdna_prepare}} and
#' \code{\link{ctdna_prep_add}} before the dnaseq frame lands in
#' `prep$dnaseq`, so downstream classifier + oncoprint work with clean
#' canonical values with zero user action needed. Direct manipulation via
#' this function is only needed when your vendor uses strings not covered
#' by the shipped defaults.
#'
#' @param ... Zero args to list the whole dict; one character key
#'   (`"sub"` or `"sub.entry"`) to fetch; OR named args of shape
#'   `sub.entry = c(alias, ...)` to APPEND vendor-specific aliases to a
#'   canonical entry (duplicates removed case-insensitively).
#'
#' @return The whole dictionary, the requested slice, or (when setting)
#'   the updated dict invisibly.
#'
#' @seealso \code{\link{dnaseq_dict_add}}, \code{\link{dnaseq_dict_reset}},
#'   \code{\link{dnaseq_dict_lookup}}, \code{\link{ctdna_opts}}.
#'
#' @examples
#' \dontrun{
#'   dnaseq_dict()                                       # whole dict
#'   dnaseq_dict("variant_type")                         # one sub-dict
#'   dnaseq_dict("molecular_consequence.missense")       # one entry
#'
#'   # Append vendor-specific aliases
#'   dnaseq_dict(molecular_consequence.missense = c("MyVendor_Missense"))
#'   dnaseq_dict(cnv_type.focal_amplification   = c("HighGain","HighAmp"))
#'
#'   dnaseq_dict_reset()                                 # back to defaults
#' }
#'
#' @export
dnaseq_dict <- function(...) {
  args <- list(...)
  dict <- .dnaseq_dict_get_state()

  # Zero args -> whole dict
  if (!length(args)) return(dict)

  # Named args -> set/append
  if (length(names(args)) && all(nzchar(names(args)))) {
    for (k in names(args)) {
      parts <- strsplit(k, ".", fixed = TRUE)[[1]]
      if (length(parts) != 2L)
        stop("dnaseq_dict: setter keys must be 'sub.entry' (e.g. ",
             "'molecular_consequence.missense'). Got: '", k, "'.",
             call. = FALSE)
      sub <- parts[1]; entry <- parts[2]
      if (!sub %in% names(dict))
        stop("dnaseq_dict: unknown sub-dict '", sub,
             "'. Valid: ", paste(names(dict), collapse = ", "), ".",
             call. = FALSE)
      if (!entry %in% names(dict[[sub]]))
        stop("dnaseq_dict: unknown canonical entry '", entry,
             "' in sub-dict '", sub, "'. Valid: ",
             paste(names(dict[[sub]]), collapse = ", "), ".",
             call. = FALSE)
      new_v <- as.character(args[[k]])
      old_v <- dict[[sub]][[entry]]
      combined <- c(old_v, new_v)
      keep <- !duplicated(.dnaseq_dict_norm_key(combined))
      dict[[sub]][[entry]] <- combined[keep]
    }
    assign("dict", dict, envir = .dnaseq_dict_env)
    return(invisible(dict))
  }

  # Positional single character key -> get
  if (length(args) == 1L && is.character(args[[1]])) {
    key <- args[[1]]
    parts <- strsplit(key, ".", fixed = TRUE)[[1]]
    if (length(parts) == 1L) {
      if (!parts %in% names(dict))
        stop("dnaseq_dict: unknown sub-dict '", parts,
             "'. Valid: ", paste(names(dict), collapse = ", "), ".",
             call. = FALSE)
      return(dict[[parts]])
    }
    if (length(parts) == 2L) {
      sub <- parts[1]; entry <- parts[2]
      if (!sub %in% names(dict))
        stop("dnaseq_dict: unknown sub-dict '", sub, "'.", call. = FALSE)
      if (!entry %in% names(dict[[sub]]))
        stop("dnaseq_dict: unknown entry '", entry, "' in '", sub, "'.",
             call. = FALSE)
      return(dict[[sub]][[entry]])
    }
    stop("dnaseq_dict: key must be 'sub' or 'sub.entry'.", call. = FALSE)
  }

  stop("dnaseq_dict: pass no args (list), one character key (get), ",
       "or named args (set/append).", call. = FALSE)
}

#' Add a vendor-specific alias to the dnaseq dictionary
#'
#' Convenience wrapper for
#' \code{dnaseq_dict("sub.entry" = aliases)}.
#'
#' @param sub One of `"variant_type"`, `"molecular_consequence"`,
#'   `"cnv_type"`.
#' @param canonical The canonical entry name (must already exist).
#' @param aliases Character vector of vendor strings to add.
#'
#' @return The updated dictionary, invisibly.
#' @export
dnaseq_dict_add <- function(sub, canonical, aliases) {
  key <- paste(sub, canonical, sep = ".")
  do.call(dnaseq_dict, setNames(list(as.character(aliases)), key))
}

#' Reset the dnaseq dictionary to shipped defaults
#' @return The (now-default) dictionary, invisibly.
#' @export
dnaseq_dict_reset <- function() {
  assign("dict", .dnaseq_dict_defaults(), envir = .dnaseq_dict_env)
  invisible(.dnaseq_dict_get_state())
}

#' Look up raw vendor value(s) and return the canonical name(s)
#'
#' Handles combined `+`-separated values for `sub == "molecular_consequence"`
#' by splitting and returning the highest-priority component.
#'
#' @param sub One of `"variant_type"`, `"molecular_consequence"`,
#'   `"cnv_type"`.
#' @param value Character vector of raw vendor strings. `NA` and empty
#'   strings return `NA_character_`.
#'
#' @return Character vector of canonical names, same length as `value`.
#'   Unmapped inputs return `NA_character_`.
#'
#' @examples
#' \dontrun{
#'   dnaseq_dict_lookup("molecular_consequence",
#'     c("missense_variant",
#'       "frameshift_variant+splice_donor_variant+intron_variant",
#'       "STOP GAINED",
#'       "banana"))
#'   # [1] "missense" "frameshift" "nonsense" NA
#' }
#'
#' @export
dnaseq_dict_lookup <- function(sub, value) {
  dict <- .dnaseq_dict_get_state()
  if (!sub %in% names(dict))
    stop("dnaseq_dict_lookup: unknown sub-dict '", sub, "'.", call. = FALSE)
  flat <- .dnaseq_dict_flat(sub)

  if (identical(sub, "molecular_consequence")) {
    priority <- .dnaseq_dict_priority()
    vapply(as.character(value), function(v) {
      if (is.na(v) || !nzchar(v)) return(NA_character_)
      parts <- strsplit(v, "+", fixed = TRUE)[[1]]
      keys  <- .dnaseq_dict_norm_key(parts)
      canons <- unname(flat[keys])
      canons <- canons[!is.na(canons)]
      if (!length(canons)) return(NA_character_)
      idx <- which(priority %in% canons)
      if (!length(idx)) return(NA_character_)
      priority[min(idx)]
    }, character(1))
  } else {
    vapply(as.character(value), function(v) {
      if (is.na(v) || !nzchar(v)) return(NA_character_)
      k <- .dnaseq_dict_norm_key(v)
      out <- flat[[k]]
      if (is.null(out) || is.na(out)) NA_character_ else unname(out)
    }, character(1))
  }
}

# ============================================================================
# INTERNAL: apply the dictionary to a dnaseq data.frame.
# Wired into .prep_dnaseq() so ctdna_prepare() / ctdna_prep_add() run this
# automatically. Users typically don't call this directly.
#
# Behaviour:
#   * Source columns are searched in this priority (first hit wins):
#       Variant_type          <- Variant_type | Variant Type | variant_type
#       Molecular_consequence <- Molecular_consequence | Variant_Effect |
#                                Functional_Class | ... (space variants)
#       CNV_type              <- CNV_type | CNA_type | CNA Type | ...
#   * Original values are preserved as *_raw columns for auditability.
#   * Canonical values are written to standard-named columns
#     (Variant_type, Molecular_consequence, CNV_type) that
#     .oncoprint_classify() reads.
# ============================================================================

.dnaseq_dict_apply <- function(df, verbose = TRUE) {
  if (!is.data.frame(df) || !nrow(df)) return(df)

  vt_col <- head(intersect(c("Variant_type","Variant Type","variant_type"),
                            names(df)), 1L)
  mc_col <- head(intersect(c("Molecular_consequence","Variant_Effect",
                               "Functional_Class",
                               "Molecular Consequence","Variant Effect",
                               "Functional Class"),
                            names(df)), 1L)
  cn_col <- head(intersect(c("CNV_type","CNA_type","CNA Type","CNV Type",
                               "cna_type","cnv_type"),
                            names(df)), 1L)

  if (length(vt_col)) {
    df$Variant_type_raw <- df[[vt_col]]
    df$Variant_type     <- unname(dnaseq_dict_lookup("variant_type",
                                                       df[[vt_col]]))
  }
  if (length(mc_col)) {
    df$Molecular_consequence_raw <- df[[mc_col]]
    df$Molecular_consequence     <- unname(dnaseq_dict_lookup(
      "molecular_consequence", df[[mc_col]]))
  }
  if (length(cn_col)) {
    df$CNV_type_raw <- df[[cn_col]]
    df$CNV_type     <- unname(dnaseq_dict_lookup("cnv_type", df[[cn_col]]))
  }

  if (isTRUE(verbose)) {
    tot  <- nrow(df)
    n_vt <- if ("Variant_type"          %in% names(df)) sum(!is.na(df$Variant_type)) else 0L
    n_mc <- if ("Molecular_consequence" %in% names(df)) sum(!is.na(df$Molecular_consequence)) else 0L
    n_cn <- if ("CNV_type"              %in% names(df)) sum(!is.na(df$CNV_type)) else 0L
    message(sprintf(
      "[dnaseq_dict] normalised %d row(s): Variant_type=%d, Molecular_consequence=%d, CNV_type=%d",
      tot, n_vt, n_mc, n_cn))
  }
  df
}
