# ============================================================================
# ctdnaTM — population-frequency annotation (v0.25.0)
# ============================================================================
# Annotates an alteration data frame with gnomAD population frequencies.
#
# Architecture: cache-then-annotate.
#   Stage A — build per-gene caches by querying gnomAD's GraphQL API.
#   Stage B — annotate rows by looking up values in the caches (no API).
#
# Four variant-type paths, three caches:
#   1. Short variants (SNV / indel) → gnomad_r2_1 (GRCh37). Cache:
#      .ctdna_cache_shortvar.  Lookup by chrom-pos-ref-alt.
#   2. CNVs (focal_amp, aneuploid_amp, homozygous_del, loh_del, amp) →
#      gnomAD-CNV (GRCh38). Cache: .ctdna_cache_cnv. Lookup by gene + type.
#   3. SV — Fusion (Variant_type == "Fusion") → gnomad_sv_r2_1.
#      Cache: .ctdna_cache_sv. Lookup by gene + breakpoint proximity.
#   4. SV — LGR (Variant_type == "LGR") → gnomad_sv_r2_1.
#      Cache: .ctdna_cache_sv. Lookup by gene with type filter.
#
# Anything else (NA Mut_nt that isn't a CNV/Fusion/LGR) gets NA with
# reason "no gnomAD source for this variant type" in the summary.

# ---- Caches (session-scoped) -----------------------------------------------
# Each cache key includes dataset + assembly to avoid cross-version collisions.
.ctdna_cache_shortvar <- new.env(parent = emptyenv())
.ctdna_cache_cnv      <- new.env(parent = emptyenv())
.ctdna_cache_sv       <- new.env(parent = emptyenv())
.ctdna_cache_meta     <- new.env(parent = emptyenv())  # save/load metadata

.ctdna_popfreq_cache_clear <- function() {
  for (e in list(.ctdna_cache_shortvar, .ctdna_cache_cnv,
                  .ctdna_cache_sv, .ctdna_cache_meta))
    rm(list = ls(envir = e, all.names = TRUE), envir = e)
}

# ---- Utilities --------------------------------------------------------------

.ctdna_normalize_chrom <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^chromosome[\\s_]*", "", x, ignore.case = TRUE, perl = TRUE)
  x <- sub("^chr[\\s_]*",        "", x, ignore.case = TRUE, perl = TRUE)
  toupper(x)
}

.ctdna_normalize_pos <- function(x) {
  # Strip commas and whitespace; cast to integer; warn on coercion loss
  y <- as.character(x)
  y <- gsub("[,\\s]", "", y, perl = TRUE)
  out <- suppressWarnings(as.integer(y))
  out
}

# Format a seconds value as "Xm Ys" (or "X.Xs" if <60s). Used at the end
# of each progress bar to report time spent on that path.
.ctdna_format_elapsed <- function(secs) {
  secs <- as.numeric(secs)
  if (is.na(secs) || secs < 0) return("?")
  if (secs < 60) return(sprintf("%.1fs", secs))
  m <- floor(secs / 60)
  s <- secs - 60 * m
  sprintf("%dm %.0fs", as.integer(m), s)
}

# Universal Mut_nt parser (carries over from v0.24.0). Returns (ref, alt)
# with NA for rows that don't have explicit ref>alt nucleotides.
.ctdna_parse_mut_nt <- function(mut_nt) {
  s <- as.character(mut_nt)
  s <- sub("\\s*\\([^)]*\\)\\s*$", "", s)                         # parenthetical
  s <- sub("^chr[0-9XYM]+:", "", s, ignore.case = TRUE)           # chr:
  s <- sub("^[gcnmr]\\.", "", s, ignore.case = TRUE)              # HGVS prefix
  s <- sub("^[0-9]+", "", s)                                       # leading pos
  m <- regmatches(s, regexec("^([ACGTN-]+)>([ACGTN-]+)$", s,
                              ignore.case = TRUE))
  ref <- vapply(m, function(x) if (length(x) >= 3) toupper(x[2]) else NA_character_,
                 character(1))
  alt <- vapply(m, function(x) if (length(x) >= 3) toupper(x[3]) else NA_character_,
                 character(1))
  data.frame(ref = ref, alt = alt, stringsAsFactors = FALSE)
}

# Guardant CNV-subtype (from CNV_type column) → gnomAD-CNV type
.ctdna_cnv_type_map <- function(vtype) {
  v <- tolower(as.character(vtype))
  ifelse(grepl("amp|amplif", v), "DUP",
         ifelse(grepl("del|loh", v), "DEL", NA_character_))
}

# Classify rows into routing paths. Returns a character vector of:
#   "short" | "cnv" | "fusion" | "lgr" | "none"
#
# CNV routing (Guardant Infinity convention):
#   - `Variant_type == "CNV"` flags the row as a copy-number call
#   - the subtype (focal_amplification, homozygous_deletion, etc.) lives
#     in the separate `CNV_type` column (default cnv_type_col = "CNV_type")
#   - .ctdna_cnv_type_map() is then called on CNV_type at lookup time
# For backward compatibility, the 5 Guardant CNV subtype strings appearing
# DIRECTLY in Variant_type are also treated as CNV.
.ctdna_classify_rows <- function(df, mut_nt_col, vtype_col) {
  vtype <- if (vtype_col %in% names(df)) as.character(df[[vtype_col]])
            else rep(NA_character_, nrow(df))
  parsed <- .ctdna_parse_mut_nt(df[[mut_nt_col]])
  has_explicit_refalt <- !is.na(parsed$ref) & !is.na(parsed$alt)
  vlow <- tolower(vtype)
  # Primary CNV flag: Variant_type == "CNV" (Guardant Infinity convention)
  is_cnv_primary <- vlow %in% c("cnv")
  # Backward-compat: Variant_type carrying the subtype directly
  is_cnv_legacy  <- vlow %in% c("aneuploid_amplification","loh_deletion",
                                 "focal_amplification","homozygous_deletion",
                                 "amplification")
  is_cnv    <- is_cnv_primary | is_cnv_legacy
  is_fusion <- vlow %in% c("fusion")
  is_lgr    <- vlow %in% c("lgr")
  cls <- rep("none", nrow(df))
  cls[has_explicit_refalt] <- "short"
  cls[is_cnv]    <- "cnv"
  cls[is_fusion] <- "fusion"
  cls[is_lgr]    <- "lgr"
  cls
}

# ---- API failure handler ----
# ---- API: short-variant gene query (gnomad_r2_1 / GRCh37) ------------------
.ctdna_gnomad_query_short_gene <- function(gene_symbol,
                                             dataset = "gnomad_r2_1",
                                             reference_genome = "GRCh37",
                                             timeout_s = 10) {
  query <- sprintf('
{
  gene(gene_symbol: "%s", reference_genome: %s) {
    variants(dataset: %s) {
      variant_id
      pos
      exome  { ac an }
      genome { ac an }
    }
  }
}', gene_symbol, reference_genome, dataset)

  resp <- tryCatch(
    httr::POST("https://gnomad.broadinstitute.org/api",
                body = list(query = query),
                encode = "json",
                httr::timeout(timeout_s),
                httr::user_agent("ctdnaTM R package")),
    error = function(e) NULL)
  if (is.null(resp)) return(NULL)
  if (httr::status_code(resp) == 429) {
    Sys.sleep(2); resp <- tryCatch(
      httr::POST("https://gnomad.broadinstitute.org/api",
                  body = list(query = query), encode = "json",
                  httr::timeout(timeout_s),
                  httr::user_agent("ctdnaTM R package")),
      error = function(e) NULL)
    if (is.null(resp) || httr::status_code(resp) != 200) return(NULL)
  } else if (httr::status_code(resp) != 200) return(NULL)

  body <- tryCatch(httr::content(resp, as = "parsed", encoding = "UTF-8"),
                    error = function(e) NULL)
  if (is.null(body) || !is.null(body$errors)) return(NULL)
  variants <- body$data$gene$variants
  if (is.null(variants))
    return(data.frame(variant_id = character(0), gnomAD_AF = numeric(0),
                       stringsAsFactors = FALSE))
  ex_ac <- vapply(variants, function(v) v$exome$ac  %||% 0L, integer(1))
  ex_an <- vapply(variants, function(v) v$exome$an  %||% 0L, integer(1))
  gn_ac <- vapply(variants, function(v) v$genome$ac %||% 0L, integer(1))
  gn_an <- vapply(variants, function(v) v$genome$an %||% 0L, integer(1))
  tan   <- ex_an + gn_an
  af    <- ifelse(tan > 0, (ex_ac + gn_ac) / tan, NA_real_)
  data.frame(
    variant_id = vapply(variants, function(v) v$variant_id, character(1)),
    gnomAD_AF  = af,
    stringsAsFactors = FALSE)
}

# ---- API: short-variant BATCHED query via GraphQL aliasing -----------------
# Sends one HTTP POST containing N aliased gene{...} sub-queries
# (g1: gene(...), g2: gene(...), ...). Returns a NAMED list keyed by
# gene symbol, each entry either a data frame (success) or NULL (failure).
# A whole-batch failure (network, non-200) returns a list of NULLs.
.ctdna_gnomad_query_short_batch <- function(gene_symbols,
                                              dataset = "gnomad_r2_1",
                                              reference_genome = "GRCh37",
                                              timeout_s = 30) {
  if (length(gene_symbols) == 0) return(list())
  # Build aliased query
  blocks <- vapply(seq_along(gene_symbols), function(i) {
    sprintf('  g%d: gene(gene_symbol: "%s", reference_genome: %s) {
    variants(dataset: %s) { variant_id pos exome{ac an} genome{ac an} }
  }', i, gene_symbols[i], reference_genome, dataset)
  }, character(1))
  query <- paste0("{\n", paste(blocks, collapse = "\n"), "\n}")

  resp <- tryCatch(
    httr::POST("https://gnomad.broadinstitute.org/api",
                body = list(query = query),
                encode = "json",
                httr::timeout(timeout_s),
                httr::user_agent("ctdnaTM R package")),
    error = function(e) NULL)
  if (is.null(resp)) return(setNames(vector("list", length(gene_symbols)),
                                       gene_symbols))
  if (httr::status_code(resp) == 429) {
    Sys.sleep(2); resp <- tryCatch(
      httr::POST("https://gnomad.broadinstitute.org/api",
                  body = list(query = query), encode = "json",
                  httr::timeout(timeout_s),
                  httr::user_agent("ctdnaTM R package")),
      error = function(e) NULL)
    if (is.null(resp) || httr::status_code(resp) != 200)
      return(setNames(vector("list", length(gene_symbols)), gene_symbols))
  } else if (httr::status_code(resp) != 200) {
    return(setNames(vector("list", length(gene_symbols)), gene_symbols))
  }

  body <- tryCatch(httr::content(resp, as = "parsed", encoding = "UTF-8"),
                    error = function(e) NULL)
  if (is.null(body))
    return(setNames(vector("list", length(gene_symbols)), gene_symbols))
  # gnomAD returns partial data + errors when some aliases fail. We treat
  # per-alias NULL as a single-gene failure but keep successful aliases.
  out <- setNames(vector("list", length(gene_symbols)), gene_symbols)
  for (i in seq_along(gene_symbols)) {
    g <- gene_symbols[i]
    node <- body$data[[paste0("g", i)]]
    if (is.null(node)) { out[[g]] <- NULL; next }
    variants <- node$variants
    if (is.null(variants)) {
      out[[g]] <- data.frame(variant_id = character(0),
                              gnomAD_AF = numeric(0),
                              stringsAsFactors = FALSE)
      next
    }
    ex_ac <- vapply(variants, function(v) v$exome$ac  %||% 0L, integer(1))
    ex_an <- vapply(variants, function(v) v$exome$an  %||% 0L, integer(1))
    gn_ac <- vapply(variants, function(v) v$genome$ac %||% 0L, integer(1))
    gn_an <- vapply(variants, function(v) v$genome$an %||% 0L, integer(1))
    tan   <- ex_an + gn_an
    af    <- ifelse(tan > 0, (ex_ac + gn_ac) / tan, NA_real_)
    out[[g]] <- data.frame(
      variant_id = vapply(variants, function(v) v$variant_id, character(1)),
      gnomAD_AF  = af,
      stringsAsFactors = FALSE)
  }
  out
}

# ---- API: CNV gene query (gnomad_cnv_v4 / GRCh38) --------------------------
# Schema fields: copy_number_variants on gene with sc, sn, sf, type.
# gnomAD-CNV is GRCh38-only; we query by gene symbol so coordinates don't
# directly matter. The returned 'sf' is the site frequency (AF analogue).
.ctdna_gnomad_query_cnv_gene <- function(gene_symbol,
                                           dataset = "gnomad_cnv_r4",
                                           timeout_s = 10) {
  query <- sprintf('
{
  gene(gene_symbol: "%s", reference_genome: GRCh38) {
    copy_number_variants(dataset: %s) {
      variant_id
      type
      sc
      sn
      sf
    }
  }
}', gene_symbol, dataset)
  resp <- tryCatch(
    httr::POST("https://gnomad.broadinstitute.org/api",
                body = list(query = query),
                encode = "json",
                httr::timeout(timeout_s),
                httr::user_agent("ctdnaTM R package")),
    error = function(e) NULL)
  if (is.null(resp)) return(NULL)
  if (httr::status_code(resp) != 200) {
    # Possible: schema field name changed across gnomAD releases. Return
    # NULL; caller treats as API failure for this gene.
    return(NULL)
  }
  body <- tryCatch(httr::content(resp, as = "parsed", encoding = "UTF-8"),
                    error = function(e) NULL)
  if (is.null(body) || !is.null(body$errors)) return(NULL)
  cnvs <- body$data$gene$copy_number_variants
  if (is.null(cnvs))
    return(data.frame(variant_id = character(0), type = character(0),
                       gnomAD_AF = numeric(0), stringsAsFactors = FALSE))
  data.frame(
    variant_id = vapply(cnvs, function(v) v$variant_id %||% NA_character_,
                         character(1)),
    type       = vapply(cnvs, function(v) v$type %||% NA_character_,
                         character(1)),
    gnomAD_AF  = vapply(cnvs, function(v) {
      sf <- v$sf
      if (is.null(sf)) NA_real_ else as.numeric(sf)
    }, numeric(1)),
    stringsAsFactors = FALSE)
}

# ---- API: CNV BATCHED query via GraphQL aliasing ---------------------------
.ctdna_gnomad_query_cnv_batch <- function(gene_symbols,
                                            dataset = "gnomad_cnv_r4",
                                            timeout_s = 30) {
  if (length(gene_symbols) == 0) return(list())
  blocks <- vapply(seq_along(gene_symbols), function(i) {
    sprintf('  g%d: gene(gene_symbol: "%s", reference_genome: GRCh38) {
    copy_number_variants(dataset: %s) { variant_id type sc sn sf }
  }', i, gene_symbols[i], dataset)
  }, character(1))
  query <- paste0("{\n", paste(blocks, collapse = "\n"), "\n}")
  resp <- tryCatch(
    httr::POST("https://gnomad.broadinstitute.org/api",
                body = list(query = query), encode = "json",
                httr::timeout(timeout_s),
                httr::user_agent("ctdnaTM R package")),
    error = function(e) NULL)
  if (is.null(resp) || httr::status_code(resp) != 200)
    return(setNames(vector("list", length(gene_symbols)), gene_symbols))
  body <- tryCatch(httr::content(resp, as = "parsed", encoding = "UTF-8"),
                    error = function(e) NULL)
  if (is.null(body))
    return(setNames(vector("list", length(gene_symbols)), gene_symbols))
  out <- setNames(vector("list", length(gene_symbols)), gene_symbols)
  for (i in seq_along(gene_symbols)) {
    g <- gene_symbols[i]
    node <- body$data[[paste0("g", i)]]
    if (is.null(node)) { out[[g]] <- NULL; next }
    cnvs <- node$copy_number_variants
    if (is.null(cnvs)) {
      out[[g]] <- data.frame(variant_id = character(0), type = character(0),
                              gnomAD_AF = numeric(0), stringsAsFactors = FALSE)
      next
    }
    out[[g]] <- data.frame(
      variant_id = vapply(cnvs, function(v) v$variant_id %||% NA_character_,
                           character(1)),
      type       = vapply(cnvs, function(v) v$type %||% NA_character_,
                           character(1)),
      gnomAD_AF  = vapply(cnvs, function(v) {
        sf <- v$sf
        if (is.null(sf)) NA_real_ else as.numeric(sf)
      }, numeric(1)),
      stringsAsFactors = FALSE)
  }
  out
}

# ---- API: SV gene query (gnomad_sv_r2_1 / GRCh37) --------------------------
.ctdna_gnomad_query_sv_gene <- function(gene_symbol,
                                          dataset = "gnomad_sv_r2_1",
                                          reference_genome = "GRCh37",
                                          timeout_s = 10) {
  query <- sprintf('
{
  gene(gene_symbol: "%s", reference_genome: %s) {
    structural_variants(dataset: %s) {
      variant_id
      type
      chrom
      pos
      end
      ac
      an
      af
    }
  }
}', gene_symbol, reference_genome, dataset)
  resp <- tryCatch(
    httr::POST("https://gnomad.broadinstitute.org/api",
                body = list(query = query),
                encode = "json",
                httr::timeout(timeout_s),
                httr::user_agent("ctdnaTM R package")),
    error = function(e) NULL)
  if (is.null(resp)) return(NULL)
  if (httr::status_code(resp) != 200) return(NULL)
  body <- tryCatch(httr::content(resp, as = "parsed", encoding = "UTF-8"),
                    error = function(e) NULL)
  if (is.null(body) || !is.null(body$errors)) return(NULL)
  svs <- body$data$gene$structural_variants
  if (is.null(svs))
    return(data.frame(variant_id = character(0), type = character(0),
                       chrom = character(0), pos = integer(0),
                       end = integer(0), gnomAD_AF = numeric(0),
                       stringsAsFactors = FALSE))
  data.frame(
    variant_id = vapply(svs, function(v) v$variant_id %||% NA_character_,
                         character(1)),
    type       = vapply(svs, function(v) v$type %||% NA_character_,
                         character(1)),
    chrom      = vapply(svs, function(v) as.character(v$chrom %||% NA),
                         character(1)),
    pos        = vapply(svs, function(v) {
                          p <- v$pos %||% NA_integer_
                          if (is.null(p)) NA_integer_ else as.integer(p)
                        }, integer(1)),
    end        = vapply(svs, function(v) {
                          p <- v$end %||% NA_integer_
                          if (is.null(p)) NA_integer_ else as.integer(p)
                        }, integer(1)),
    gnomAD_AF  = vapply(svs, function(v) {
                          a <- v$af %||% NA_real_
                          if (is.null(a)) NA_real_ else as.numeric(a)
                        }, numeric(1)),
    stringsAsFactors = FALSE)
}

# ---- API: SV BATCHED query via GraphQL aliasing ----------------------------
.ctdna_gnomad_query_sv_batch <- function(gene_symbols,
                                           dataset = "gnomad_sv_r2_1",
                                           reference_genome = "GRCh37",
                                           timeout_s = 30) {
  if (length(gene_symbols) == 0) return(list())
  blocks <- vapply(seq_along(gene_symbols), function(i) {
    sprintf('  g%d: gene(gene_symbol: "%s", reference_genome: %s) {
    structural_variants(dataset: %s) { variant_id type chrom pos end ac an af }
  }', i, gene_symbols[i], reference_genome, dataset)
  }, character(1))
  query <- paste0("{\n", paste(blocks, collapse = "\n"), "\n}")
  resp <- tryCatch(
    httr::POST("https://gnomad.broadinstitute.org/api",
                body = list(query = query), encode = "json",
                httr::timeout(timeout_s),
                httr::user_agent("ctdnaTM R package")),
    error = function(e) NULL)
  if (is.null(resp) || httr::status_code(resp) != 200)
    return(setNames(vector("list", length(gene_symbols)), gene_symbols))
  body <- tryCatch(httr::content(resp, as = "parsed", encoding = "UTF-8"),
                    error = function(e) NULL)
  if (is.null(body))
    return(setNames(vector("list", length(gene_symbols)), gene_symbols))
  out <- setNames(vector("list", length(gene_symbols)), gene_symbols)
  for (i in seq_along(gene_symbols)) {
    g <- gene_symbols[i]
    node <- body$data[[paste0("g", i)]]
    if (is.null(node)) { out[[g]] <- NULL; next }
    svs <- node$structural_variants
    if (is.null(svs)) {
      out[[g]] <- data.frame(variant_id = character(0), type = character(0),
                              chrom = character(0), pos = integer(0),
                              end = integer(0), gnomAD_AF = numeric(0),
                              stringsAsFactors = FALSE)
      next
    }
    out[[g]] <- data.frame(
      variant_id = vapply(svs, function(v) v$variant_id %||% NA_character_,
                           character(1)),
      type       = vapply(svs, function(v) v$type %||% NA_character_,
                           character(1)),
      chrom      = vapply(svs, function(v) as.character(v$chrom %||% NA),
                           character(1)),
      pos        = vapply(svs, function(v) {
                            p <- v$pos %||% NA_integer_
                            if (is.null(p)) NA_integer_ else as.integer(p)
                          }, integer(1)),
      end        = vapply(svs, function(v) {
                            p <- v$end %||% NA_integer_
                            if (is.null(p)) NA_integer_ else as.integer(p)
                          }, integer(1)),
      gnomAD_AF  = vapply(svs, function(v) {
                            a <- v$af %||% NA_real_
                            if (is.null(a)) NA_real_ else as.numeric(a)
                          }, numeric(1)),
      stringsAsFactors = FALSE)
  }
  out
}

# %||% helper if not defined
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Save the gnomAD population-frequency cache to an RDS file
#'
#' Writes all session caches (short-variant, CNV, SV) and metadata to a
#' single RDS file. Useful for sharing the cache across colleagues or
#' persisting it across R sessions in environments where the home cache
#' is wiped.
#'
#' @param path File path. Convention: \code{*.rds}.
#' @return Invisibly, the path.
#' @examples
#' # Round-trip cache to disk
#' tmp <- tempfile(fileext = ".rds")
#' ctdna_save_popfreq_cache(tmp)
#' ctdna_load_popfreq_cache(tmp)
#' @export
ctdna_save_popfreq_cache <- function(path) {
  obj <- list(
    shortvar = as.list(.ctdna_cache_shortvar),
    cnv      = as.list(.ctdna_cache_cnv),
    sv       = as.list(.ctdna_cache_sv),
    meta     = list(
      saved_at = Sys.time(),
      package_version = as.character(utils::packageVersion("ctdnaTM"))
    )
  )
  saveRDS(obj, file = path)
  message(sprintf(
    "Saved popfreq cache to %s (shortvar: %d genes, cnv: %d genes, sv: %d genes).",
    path, length(obj$shortvar), length(obj$cnv), length(obj$sv)))
  invisible(path)
}

#' Load a gnomAD population-frequency cache from an RDS file
#'
#' Restores a cache previously written by
#' \code{\link{ctdna_save_popfreq_cache}}. Use this to avoid hitting
#' the gnomAD API repeatedly across sessions, or to share a cache
#' between collaborators on the same study.
#'
#' @param path File path written by \code{ctdna_save_popfreq_cache()}.
#' @param merge If \code{TRUE} (default), merge into the current session
#'   cache (existing keys are preserved). If \code{FALSE}, the loaded
#'   cache REPLACES the current session cache.
#' @return Invisibly, the metadata list from the file.
#' @examples
#' # Round-trip
#' tmp <- tempfile(fileext = ".rds")
#' ctdna_save_popfreq_cache(tmp)
#' meta <- ctdna_load_popfreq_cache(tmp)
#' meta$saved_at
#' @export
ctdna_load_popfreq_cache <- function(path, merge = TRUE) {
  if (!file.exists(path)) stop("File not found: ", path)
  obj <- readRDS(path)
  if (!isTRUE(merge)) .ctdna_popfreq_cache_clear()
  load_into <- function(env, lst) {
    for (k in names(lst)) {
      if (!merge || !exists(k, envir = env, inherits = FALSE))
        assign(k, lst[[k]], envir = env)
    }
  }
  load_into(.ctdna_cache_shortvar, obj$shortvar %||% list())
  load_into(.ctdna_cache_cnv,      obj$cnv      %||% list())
  load_into(.ctdna_cache_sv,       obj$sv       %||% list())
  message(sprintf(
    "Loaded popfreq cache from %s (shortvar: %d, cnv: %d, sv: %d). Save timestamp: %s",
    path, length(obj$shortvar), length(obj$cnv), length(obj$sv),
    as.character(obj$meta$saved_at %||% NA)))
  invisible(obj$meta)
}

#' Annotate alterations with gnomAD population allele frequencies
#'
#' Routes rows to one of four gnomAD query paths based on Variant_type:
#' SNV/indel (short variants), CNV, SV-Fusion, SV-LGR. Each path uses
#' gene-batched queries with separate per-gene caches.
#'
#' @section Variant routing:
#' \itemize{
#'   \item Rows with explicit \code{ref>alt} in \code{Mut_nt} → short
#'     variant cache (\code{gnomad_r2_1}, GRCh37).
#'   \item \code{Variant_type == "CNV"} → CNV cache
#'     (\code{gnomad_cnv_r4}, GRCh38 by necessity — gnomAD has no GRCh37
#'     CNV dataset). The CNV subtype is read from the separate
#'     \code{CNV_type} column (\code{cnv_type_col} arg) and mapped to
#'     gnomAD-CNV DUP/DEL. For backward compatibility, the 5 Guardant
#'     CNV subtype strings appearing directly in \code{Variant_type} are
#'     also treated as CNV.
#'   \item \code{Variant_type == "Fusion"} → SV cache (\code{gnomad_sv_r2_1}),
#'     matched against gene + Fusion_position_b ± \code{fusion_bp_tolerance}.
#'   \item \code{Variant_type == "LGR"} → SV cache, gene-level lookup.
#'   \item Anything else → \code{NA} with reason "no gnomAD source".
#' }
#'
#' @section Filter NA-safety:
#' Rows that end up with \code{gnomAD_AF = NA} (any reason) are kept by
#' \code{ctdna_filter_apply()} (NA-safe by default).
#'
#' @param df Alteration data frame.
#' @param chrom_col,pos_col,mut_nt_col,gene_col,vtype_col Column names.
#'   Defaults match the Guardant Infinity report.
#' @param cnv_type_col Column carrying the CNV subtype for rows where
#'   \code{Variant_type == "CNV"}. Default \code{"CNV_type"} (Guardant
#'   Infinity convention). Expected values: \code{aneuploid_amplification},
#'   \code{loh_deletion}, \code{focal_amplification},
#'   \code{homozygous_deletion}, \code{amplification}. The function maps
#'   these to gnomAD-CNV \code{DUP} / \code{DEL} for lookup.
#' @param fusion_gene_b_col,fusion_pos_a_col,fusion_pos_b_col Fusion-specific
#'   column names. Defaults: \code{"Fusion_gene_b"}, \code{"Fusion_position_a"},
#'   \code{"Fusion_position_b"}.
#' @param gnomad_af_col Output column name for the AF value.
#' @param dataset Short-variant dataset. Default \code{"gnomad_r2_1"}.
#' @param reference_genome Short-variant assembly. Default \code{"GRCh37"}.
#' @param cnv_dataset CNV dataset. Default \code{"gnomad_cnv_r4"} (GRCh38).
#' @param sv_dataset SV dataset. Default \code{"gnomad_sv_r2_1"}.
#' @param fusion_bp_tolerance Breakpoint tolerance for matching a Guardant
#'   fusion against gnomAD-SV BND records, in bp. Default 500.
#' @param variant_fallback If \code{TRUE}, when the short-variant gene-batched
#'   query has no hit for a row, fall back to per-variant
#'   \code{variant(variantId: ...)} queries. Off by default (v0.25.0
#'   architectural decision to keep API workload predictable).
#' @param sleep_time Seconds to pause between gnomAD API calls (rate
#'   limiting). Default 0.3. Increase if you see HTTP 429 responses.
#' @param batch_size Number of genes to query per HTTP request via
#'   GraphQL aliasing. Default 10. Larger values reduce round-trips;
#'   too-large values risk hitting gnomAD's per-request size limit or
#'   timing out. Each path (short variants, CNV, SV) uses the same
#'   batch size.
#' @param retry Controls automatic retry of failed-gene API calls.
#'   Accepts:
#'   \itemize{
#'     \item \code{TRUE} (default) - always retry failed genes once.
#'     \item \code{FALSE} - never retry; leave failed rows as NA.
#'     \item Numeric in (0, 1] - retry only if the success rate fell
#'       BELOW this threshold. For example, \code{retry = 0.95} means
#'       "skip retry only if at least 95\% of genes succeeded". Useful
#'       for accepting a partial result on small transient hiccups but
#'       still retrying on big failures.
#'   }
#'   The y/n prompt from earlier versions is replaced by this argument.
#' @param overwrite If \code{TRUE}, refetch rows that already have a value.
#' @param verbose Show progress messages. Default \code{TRUE}.
#' @return The input frame with \code{gnomAD_AF}, \code{gnomAD_version},
#'   and \code{gnomAD_id} columns populated. \code{gnomAD_id} is the
#'   gnomAD variant identifier of the matched record (e.g.
#'   \code{"17-7577538-G-A"} for a short variant) and is \code{NA}
#'   when no record was matched.
#' @export
ctdna_annotate_population_freq <- function(df,
        chrom_col          = "Chromosome",
        pos_col            = "Position",
        mut_nt_col         = "Mut_nt",
        gene_col           = "Gene",
        vtype_col          = "Variant_type",
        cnv_type_col       = "CNV_type",
        fusion_gene_b_col  = "Fusion_gene_b",
        fusion_pos_a_col   = "Fusion_position_a",
        fusion_pos_b_col   = "Fusion_position_b",
        gnomad_af_col      = "gnomAD_AF",
        dataset            = "gnomad_r2_1",
        reference_genome   = "GRCh37",
        cnv_dataset        = "gnomad_cnv_r4",
        sv_dataset         = "gnomad_sv_r2_1",
        fusion_bp_tolerance = 500L,
        variant_fallback   = FALSE,
        sleep_time         = 0.3,
        batch_size         = 10L,
        retry              = TRUE,
        overwrite          = FALSE,
        verbose            = TRUE) {
  if (!is.numeric(sleep_time) || sleep_time < 0)
    stop("`sleep_time` must be a non-negative number (seconds).")
  if (!is.numeric(batch_size) || batch_size < 1)
    stop("`batch_size` must be a positive integer.")
  batch_size <- as.integer(batch_size)
  # Validate retry:
  #   TRUE  -> always retry failed genes
  #   FALSE -> never retry
  #   numeric in (0, 1] -> retry only if success rate fell BELOW this threshold
  #     e.g. retry = 0.95 means "skip retry only if at least 95% succeeded"
  if (is.logical(retry)) {
    if (length(retry) != 1L || is.na(retry))
      stop("`retry` (logical) must be a single TRUE or FALSE.")
  } else if (is.numeric(retry)) {
    if (length(retry) != 1L || is.na(retry) || retry <= 0 || retry > 1)
      stop("`retry` (numeric) must be a single value in (0, 1] ",
           "indicating the success-rate threshold below which to retry.")
  } else {
    stop("`retry` must be TRUE, FALSE, or a numeric threshold in (0, 1].")
  }
  if (!is.data.frame(df)) stop("`df` must be a data frame.")
  needed <- c(chrom_col, pos_col, mut_nt_col, gene_col)
  missing_cols <- setdiff(needed, names(df))
  if (length(missing_cols) > 0)
    stop(sprintf("Required column(s) missing from df: %s.",
                  paste(missing_cols, collapse = ", ")))

  if (!gnomad_af_col   %in% names(df)) df[[gnomad_af_col]] <- NA_real_
  if (!"gnomAD_version" %in% names(df)) df$gnomAD_version  <- NA_character_
  if (!"gnomAD_id"      %in% names(df)) df$gnomAD_id       <- NA_character_
  # Track per-row classification reason — useful for summary
  df$.gnomAD_class <- rep(NA_character_, nrow(df))

  # Normalize chrom/pos once
  df$.chrom_norm <- .ctdna_normalize_chrom(df[[chrom_col]])
  df$.pos_norm  <- .ctdna_normalize_pos(df[[pos_col]])

  parsed <- .ctdna_parse_mut_nt(df[[mut_nt_col]])
  df$.ref <- parsed$ref
  df$.alt <- parsed$alt

  classes <- .ctdna_classify_rows(df, mut_nt_col, vtype_col)

  # Identify rows that need fetching (skip already-populated unless overwrite)
  needs_af <- if (overwrite) rep(TRUE, nrow(df))
              else is.na(df[[gnomad_af_col]])

  failed_genes_short <- character()
  failed_genes_cnv   <- character()
  failed_genes_sv    <- character()

  pause <- function() Sys.sleep(sleep_time)

  # ---- Path 1: short variants -----------------------------------------
  short_idx <- which(classes == "short" & needs_af)
  if (length(short_idx) > 0) {
    genes <- unique(df[[gene_col]][short_idx])
    genes <- genes[!is.na(genes) & nzchar(as.character(genes))]
    if (verbose && length(genes) > 0)
      message(sprintf(
        "[1/3] Short variants: fetching %d gene(s) from %s (batch_size=%d) ...",
        length(genes), dataset, batch_size))
    t_start <- Sys.time()
    if (verbose && length(genes) > 0)
      pb <- utils::txtProgressBar(min = 0, max = length(genes), style = 3)
    # Build worklist (skip already-cached)
    cached_keys <- ls(envir = .ctdna_cache_shortvar)
    to_fetch <- genes[!paste(genes, dataset, reference_genome, sep = "|") %in%
                        cached_keys]
    n_cached_initial <- length(genes) - length(to_fetch)
    if (verbose && n_cached_initial > 0)
      utils::setTxtProgressBar(pb, n_cached_initial)
    consecutive_fails <- 0L
    circuit_broken <- FALSE
    seen <- n_cached_initial
    if (length(to_fetch) > 0) {
      starts <- seq.int(1L, length(to_fetch), by = batch_size)
      for (st in starts) {
        chunk <- to_fetch[st:min(st + batch_size - 1L, length(to_fetch))]
        if (circuit_broken) {
          failed_genes_short <- c(failed_genes_short, chunk)
          seen <- seen + length(chunk)
          if (verbose) utils::setTxtProgressBar(pb, seen)
          next
        }
        res <- .ctdna_gnomad_query_short_batch(
          chunk, dataset = dataset, reference_genome = reference_genome)
        all_null <- all(vapply(res, is.null, logical(1)))
        for (g in chunk) {
          if (is.null(res[[g]])) {
            failed_genes_short <- c(failed_genes_short, g)
          } else {
            assign(paste(g, dataset, reference_genome, sep = "|"),
                    res[[g]], envir = .ctdna_cache_shortvar)
          }
        }
        if (all_null) {
          consecutive_fails <- consecutive_fails + 1L
          if (consecutive_fails >= 3L) {
            circuit_broken <- TRUE
            if (verbose) message(sprintf(
              "  [circuit-break] 3 consecutive batch failures; skipping remaining %d gene(s).",
              length(to_fetch) - (st + length(chunk) - 1L)))
          }
        } else {
          consecutive_fails <- 0L
        }
        seen <- seen + length(chunk)
        if (verbose) utils::setTxtProgressBar(pb, seen)
        pause()
      }
    }
    if (verbose && length(genes) > 0) {
      close(pb)
      message(sprintf("  [1/3] elapsed: %s",
                       .ctdna_format_elapsed(
                         as.numeric(difftime(Sys.time(), t_start, units = "secs")))))
    }

    # Annotate rows from cache
    for (r in short_idx) {
      g <- df[[gene_col]][r]
      if (is.na(g) || !nzchar(as.character(g))) {
        df$.gnomAD_class[r] <- "no_gene"; next
      }
      key <- paste(g, dataset, reference_genome, sep = "|")
      if (!exists(key, envir = .ctdna_cache_shortvar)) {
        df$.gnomAD_class[r] <- "api_failed_short"; next
      }
      tab <- get(key, envir = .ctdna_cache_shortvar)
      row_key <- sprintf("%s-%s-%s-%s",
                          df$.chrom_norm[r], df$.pos_norm[r],
                          df$.ref[r], df$.alt[r])
      hit <- match(row_key, tab$variant_id)
      if (is.na(hit)) {
        df[[gnomad_af_col]][r] <- 0  # absent from gnomAD's gene records = rare
        df$gnomAD_version[r]   <- dataset
        df$gnomAD_id[r]        <- NA_character_   # nothing matched in gnomAD
        df$.gnomAD_class[r]    <- "absent"
      } else {
        af <- tab$gnomAD_AF[hit]
        df[[gnomad_af_col]][r] <- if (is.na(af)) 0 else af
        df$gnomAD_version[r]   <- dataset
        df$gnomAD_id[r]        <- tab$variant_id[hit]
        df$.gnomAD_class[r]    <- "annotated"
      }
    }
  }

  # ---- Path 2: CNV ----------------------------------------------------
  cnv_idx <- which(classes == "cnv" & needs_af)
  if (length(cnv_idx) > 0) {
    genes <- unique(df[[gene_col]][cnv_idx])
    genes <- genes[!is.na(genes) & nzchar(as.character(genes))]
    if (verbose && length(genes) > 0)
      message(sprintf(
        "[2/3] CNVs: fetching %d gene(s) from %s (batch_size=%d) ...",
        length(genes), cnv_dataset, batch_size))
    t_start <- Sys.time()
    if (verbose && length(genes) > 0)
      pb <- utils::txtProgressBar(min = 0, max = length(genes), style = 3)
    cached_keys <- ls(envir = .ctdna_cache_cnv)
    to_fetch <- genes[!paste(genes, cnv_dataset, sep = "|") %in% cached_keys]
    n_cached_initial <- length(genes) - length(to_fetch)
    if (verbose && n_cached_initial > 0)
      utils::setTxtProgressBar(pb, n_cached_initial)
    consecutive_fails <- 0L
    circuit_broken <- FALSE
    seen <- n_cached_initial
    if (length(to_fetch) > 0) {
      starts <- seq.int(1L, length(to_fetch), by = batch_size)
      for (st in starts) {
        chunk <- to_fetch[st:min(st + batch_size - 1L, length(to_fetch))]
        if (circuit_broken) {
          failed_genes_cnv <- c(failed_genes_cnv, chunk)
          seen <- seen + length(chunk)
          if (verbose) utils::setTxtProgressBar(pb, seen)
          next
        }
        res <- .ctdna_gnomad_query_cnv_batch(chunk, dataset = cnv_dataset)
        all_null <- all(vapply(res, is.null, logical(1)))
        for (g in chunk) {
          if (is.null(res[[g]])) {
            failed_genes_cnv <- c(failed_genes_cnv, g)
          } else {
            assign(paste(g, cnv_dataset, sep = "|"),
                    res[[g]], envir = .ctdna_cache_cnv)
          }
        }
        if (all_null) {
          consecutive_fails <- consecutive_fails + 1L
          if (consecutive_fails >= 3L) {
            circuit_broken <- TRUE
            if (verbose) message(sprintf(
              "  [circuit-break CNV] 3 consecutive batch failures; skipping remaining %d gene(s).",
              length(to_fetch) - (st + length(chunk) - 1L)))
          }
        } else {
          consecutive_fails <- 0L
        }
        seen <- seen + length(chunk)
        if (verbose) utils::setTxtProgressBar(pb, seen)
        pause()
      }
    }
    if (verbose && length(genes) > 0) {
      close(pb)
      message(sprintf("  [2/3] elapsed: %s",
                       .ctdna_format_elapsed(
                         as.numeric(difftime(Sys.time(), t_start, units = "secs")))))
    }

    # Annotate from cache. Aggregation: highest AF among matching type.
    # CNV subtype is read from CNV_type column (cnv_type_col arg).
    # If CNV_type column is absent OR the row's CNV_type is NA, fall
    # back to reading Variant_type (handles legacy data where the
    # subtype string was carried in Variant_type itself).
    has_cnv_type_col <- cnv_type_col %in% names(df)
    for (r in cnv_idx) {
      g <- df[[gene_col]][r]
      if (is.na(g) || !nzchar(as.character(g))) {
        df$.gnomAD_class[r] <- "no_gene"; next
      }
      key <- paste(g, cnv_dataset, sep = "|")
      if (!exists(key, envir = .ctdna_cache_cnv)) {
        df$.gnomAD_class[r] <- "api_failed_cnv"; next
      }
      tab <- get(key, envir = .ctdna_cache_cnv)
      # Pull subtype from CNV_type, else Variant_type as fallback
      sub <- if (has_cnv_type_col) df[[cnv_type_col]][r] else NA_character_
      if (is.na(sub) || !nzchar(as.character(sub)))
        sub <- df[[vtype_col]][r]
      want_type <- .ctdna_cnv_type_map(sub)
      if (nrow(tab) == 0 || is.na(want_type)) {
        df[[gnomad_af_col]][r] <- 0
        df$gnomAD_version[r]   <- cnv_dataset
        df$gnomAD_id[r]        <- NA_character_
        df$.gnomAD_class[r]    <- "absent"
        next
      }
      hits <- tab[tab$type == want_type, , drop = FALSE]
      if (nrow(hits) == 0) {
        df[[gnomad_af_col]][r] <- 0
        df$gnomAD_version[r]   <- cnv_dataset
        df$gnomAD_id[r]        <- NA_character_
        df$.gnomAD_class[r]    <- "absent"
      } else {
        # which.max picks the first NA-skipping max; the variant_id of
        # that record is what we report alongside the AF
        best <- which.max(hits$gnomAD_AF)
        df[[gnomad_af_col]][r] <- hits$gnomAD_AF[best]
        df$gnomAD_version[r]   <- cnv_dataset
        df$gnomAD_id[r]        <- hits$variant_id[best]
        df$.gnomAD_class[r]    <- "annotated"
      }
    }
  }

  # ---- Path 3: SV (fusion + LGR) --------------------------------------
  sv_idx <- which(classes %in% c("fusion","lgr") & needs_af)
  if (length(sv_idx) > 0) {
    genes <- unique(df[[gene_col]][sv_idx])
    genes <- genes[!is.na(genes) & nzchar(as.character(genes))]
    if (verbose && length(genes) > 0)
      message(sprintf(
        "[3/3] SVs (fusion + LGR): fetching %d gene(s) from %s (batch_size=%d) ...",
        length(genes), sv_dataset, batch_size))
    t_start <- Sys.time()
    if (verbose && length(genes) > 0)
      pb <- utils::txtProgressBar(min = 0, max = length(genes), style = 3)
    cached_keys <- ls(envir = .ctdna_cache_sv)
    to_fetch <- genes[!paste(genes, sv_dataset, reference_genome, sep = "|")
                        %in% cached_keys]
    n_cached_initial <- length(genes) - length(to_fetch)
    if (verbose && n_cached_initial > 0)
      utils::setTxtProgressBar(pb, n_cached_initial)
    consecutive_fails <- 0L
    circuit_broken <- FALSE
    seen <- n_cached_initial
    if (length(to_fetch) > 0) {
      starts <- seq.int(1L, length(to_fetch), by = batch_size)
      for (st in starts) {
        chunk <- to_fetch[st:min(st + batch_size - 1L, length(to_fetch))]
        if (circuit_broken) {
          failed_genes_sv <- c(failed_genes_sv, chunk)
          seen <- seen + length(chunk)
          if (verbose) utils::setTxtProgressBar(pb, seen)
          next
        }
        res <- .ctdna_gnomad_query_sv_batch(
          chunk, dataset = sv_dataset, reference_genome = reference_genome)
        all_null <- all(vapply(res, is.null, logical(1)))
        for (g in chunk) {
          if (is.null(res[[g]])) {
            failed_genes_sv <- c(failed_genes_sv, g)
          } else {
            assign(paste(g, sv_dataset, reference_genome, sep = "|"),
                    res[[g]], envir = .ctdna_cache_sv)
          }
        }
        if (all_null) {
          consecutive_fails <- consecutive_fails + 1L
          if (consecutive_fails >= 3L) {
            circuit_broken <- TRUE
            if (verbose) message(sprintf(
              "  [circuit-break SV] 3 consecutive batch failures; skipping remaining %d gene(s).",
              length(to_fetch) - (st + length(chunk) - 1L)))
          }
        } else {
          consecutive_fails <- 0L
        }
        seen <- seen + length(chunk)
        if (verbose) utils::setTxtProgressBar(pb, seen)
        pause()
      }
    }
    if (verbose && length(genes) > 0) {
      close(pb)
      message(sprintf("  [3/3] elapsed: %s",
                       .ctdna_format_elapsed(
                         as.numeric(difftime(Sys.time(), t_start, units = "secs")))))
    }

    # Annotate rows from cache
    for (r in sv_idx) {
      g <- df[[gene_col]][r]
      if (is.na(g) || !nzchar(as.character(g))) {
        df$.gnomAD_class[r] <- "no_gene"; next
      }
      key <- paste(g, sv_dataset, reference_genome, sep = "|")
      if (!exists(key, envir = .ctdna_cache_sv)) {
        df$.gnomAD_class[r] <- "api_failed_sv"; next
      }
      tab <- get(key, envir = .ctdna_cache_sv)
      if (nrow(tab) == 0) {
        df[[gnomad_af_col]][r] <- 0
        df$gnomAD_version[r]   <- sv_dataset
        df$gnomAD_id[r]        <- NA_character_
        df$.gnomAD_class[r]    <- "absent"
        next
      }
      if (classes[r] == "fusion") {
        # Match gene + BND type + Fusion_position_b within tolerance
        if (!(fusion_pos_b_col %in% names(df)) ||
            is.na(df[[fusion_pos_b_col]][r])) {
          df[[gnomad_af_col]][r] <- 0
          df$gnomAD_version[r]   <- sv_dataset
          df$gnomAD_id[r]        <- NA_character_
          df$.gnomAD_class[r]    <- "absent"
          next
        }
        target <- as.integer(df[[fusion_pos_b_col]][r])
        cand <- tab[tab$type == "BND", , drop = FALSE]
        if (nrow(cand) == 0) {
          df[[gnomad_af_col]][r] <- 0
          df$gnomAD_version[r]   <- sv_dataset
          df$gnomAD_id[r]        <- NA_character_
          df$.gnomAD_class[r]    <- "absent"
          next
        }
        # SVs are gene-A scoped; the BND's `end` field carries the mate
        # breakpoint (gnomAD-SV convention).
        d_pos <- abs(cand$pos - target)
        d_end <- abs(cand$end - target)
        d_min <- pmin(d_pos, d_end, na.rm = FALSE)
        in_tol <- !is.na(d_min) & d_min <= fusion_bp_tolerance
        if (!any(in_tol)) {
          df[[gnomad_af_col]][r] <- 0
          df$gnomAD_version[r]   <- sv_dataset
          df$gnomAD_id[r]        <- NA_character_
          df$.gnomAD_class[r]    <- "absent"
        } else {
          cand_in <- cand[in_tol, , drop = FALSE]
          best <- which.max(cand_in$gnomAD_AF)
          df[[gnomad_af_col]][r] <- cand_in$gnomAD_AF[best]
          df$gnomAD_version[r]   <- sv_dataset
          df$gnomAD_id[r]        <- cand_in$variant_id[best]
          df$.gnomAD_class[r]    <- "annotated"
        }
      } else {
        # LGR — gene-level lookup, take highest AF among DEL/DUP records
        cand <- tab[tab$type %in% c("DEL","DUP","INV","INS"), , drop = FALSE]
        if (nrow(cand) == 0) {
          df[[gnomad_af_col]][r] <- 0
          df$gnomAD_version[r]   <- sv_dataset
          df$gnomAD_id[r]        <- NA_character_
          df$.gnomAD_class[r]    <- "absent"
        } else {
          best <- which.max(cand$gnomAD_AF)
          df[[gnomad_af_col]][r] <- cand$gnomAD_AF[best]
          df$gnomAD_version[r]   <- sv_dataset
          df$gnomAD_id[r]        <- cand$variant_id[best]
          df$.gnomAD_class[r]    <- "annotated"
        }
      }
    }
  }

  # ---- Rows that didn't fit any path ----------------------------------
  none_idx <- which(classes == "none" & needs_af)
  if (length(none_idx) > 0) {
    df$.gnomAD_class[none_idx] <- "no_source"
    # gnomAD_AF stays NA (NA-safe at the filter)
  }

  # ---- Retry decision (driven by `retry` argument, no prompt) ---------
  total_failed_genes <- length(unique(c(failed_genes_short,
                                          failed_genes_cnv,
                                          failed_genes_sv)))
  total_genes_attempted <-
    length(unique(c(unique(df[[gene_col]][short_idx]),
                     unique(df[[gene_col]][cnv_idx]),
                     unique(df[[gene_col]][sv_idx])))) %||% 0
  # Decide whether to retry:
  #   retry == TRUE  -> retry whenever there are failed genes
  #   retry == FALSE -> never
  #   retry numeric  -> only retry if success rate fell BELOW threshold
  do_retry <- FALSE
  if (total_failed_genes > 0 && total_genes_attempted > 0) {
    if (isTRUE(retry)) {
      do_retry <- TRUE
    } else if (is.numeric(retry)) {
      success_rate <- 1 - (total_failed_genes / total_genes_attempted)
      do_retry <- success_rate < retry
      if (verbose)
        message(sprintf(
          "  retry decision: %d/%d failed (success rate %.1f%%); threshold %.1f%% -> %s",
          total_failed_genes, total_genes_attempted,
          100 * success_rate, 100 * retry,
          if (do_retry) "RETRYING" else "skipping retry"))
    }  # retry == FALSE -> do_retry stays FALSE
  }

  # Wrap retry + classification in tryCatch so a crash here doesn't
  # discard the work done in Stage A/B above. Whatever's been computed
  # so far in `df` is returned with a warning.
  tryCatch({
  if (do_retry) {
      if (verbose)
        message(sprintf("Retrying %d failed gene(s)...", total_failed_genes))
      retry_short <- failed_genes_short
      retry_cnv   <- failed_genes_cnv
      retry_sv    <- failed_genes_sv

      # Helper: run a batched retry loop with progress bar + elapsed time
      run_retry <- function(label, genes_v, query_fn, cache_env, key_fmt) {
        if (length(genes_v) == 0) return(invisible(NULL))
        if (verbose) message(sprintf("  retry %s: %d gene(s) ...",
                                       label, length(genes_v)))
        t_r <- Sys.time()
        if (verbose)
          pb_r <- utils::txtProgressBar(min = 0, max = length(genes_v),
                                          style = 3)
        seen <- 0L
        starts <- seq.int(1L, length(genes_v), by = batch_size)
        for (st in starts) {
          chunk <- genes_v[st:min(st + batch_size - 1L, length(genes_v))]
          res <- query_fn(chunk)
          for (g in chunk) {
            if (!is.null(res[[g]]))
              assign(key_fmt(g), res[[g]], envir = cache_env)
          }
          seen <- seen + length(chunk)
          if (verbose) utils::setTxtProgressBar(pb_r, seen)
          pause()
        }
        if (verbose) {
          close(pb_r)
          message(sprintf("    %s elapsed: %s", label,
                           .ctdna_format_elapsed(
                             as.numeric(difftime(Sys.time(), t_r, units = "secs")))))
        }
      }

      run_retry("short-variant", retry_short,
                 function(ch) .ctdna_gnomad_query_short_batch(
                   ch, dataset = dataset, reference_genome = reference_genome),
                 .ctdna_cache_shortvar,
                 function(g) paste(g, dataset, reference_genome, sep = "|"))
      run_retry("CNV", retry_cnv,
                 function(ch) .ctdna_gnomad_query_cnv_batch(
                   ch, dataset = cnv_dataset),
                 .ctdna_cache_cnv,
                 function(g) paste(g, cnv_dataset, sep = "|"))
      run_retry("SV", retry_sv,
                 function(ch) .ctdna_gnomad_query_sv_batch(
                   ch, dataset = sv_dataset, reference_genome = reference_genome),
                 .ctdna_cache_sv,
                 function(g) paste(g, sv_dataset, reference_genome, sep = "|"))

      # Re-annotate the previously-failed rows after retry refreshes caches.
      # FIX (v0.26.0): the previous version called this function recursively
      # and then did `df[mark, ] <- df_remain`. That assignment fails with
      # "Can't recycle input of size X to size Y" because the recursive call
      # STRIPS internal helper columns (.gnomAD_class, .chrom_norm, .pos_norm,
      # .ref, .alt = 5 columns) before returning, so the returned frame has
      # 5 fewer columns than the parent's `df`. Fix: write back column-by-
      # column for ONLY the output columns we care about, never whole rows.
      mark <- which(df$.gnomAD_class %in%
                       c("api_failed_short","api_failed_cnv","api_failed_sv"))
      if (length(mark) > 0) {
        df_remain <- df[mark, , drop = FALSE]
        df_remain[[gnomad_af_col]] <- NA_real_
        df_remain$gnomAD_version   <- NA_character_
        df_remain$gnomAD_id        <- NA_character_
        df_remain$.gnomAD_class    <- NA_character_
        df_remain <- ctdna_annotate_population_freq(df_remain,
            chrom_col = chrom_col, pos_col = pos_col, mut_nt_col = mut_nt_col,
            gene_col = gene_col, vtype_col = vtype_col,
            cnv_type_col = cnv_type_col,
            fusion_gene_b_col = fusion_gene_b_col,
            fusion_pos_a_col = fusion_pos_a_col,
            fusion_pos_b_col = fusion_pos_b_col,
            gnomad_af_col = gnomad_af_col,
            dataset = dataset, reference_genome = reference_genome,
            cnv_dataset = cnv_dataset, sv_dataset = sv_dataset,
            fusion_bp_tolerance = fusion_bp_tolerance,
            variant_fallback = variant_fallback,
            sleep_time = sleep_time,
            batch_size = batch_size,
            retry = FALSE,  # never recurse infinitely
            overwrite = TRUE, verbose = FALSE)
        # Column-by-column write — output columns only, no internal cols.
        df[[gnomad_af_col]][mark] <- df_remain[[gnomad_af_col]]
        df$gnomAD_version[mark]   <- df_remain$gnomAD_version
        df$gnomAD_id[mark]        <- df_remain$gnomAD_id
        # Re-derive .gnomAD_class for the marked rows from the new AF values:
        # NA AF means still-failed (keep api_failed_*); non-NA AF means it
        # was annotated (AF > 0) or absent (AF == 0). We can't fully recover
        # the api_failed_short vs api_failed_cnv vs api_failed_sv split from
        # AF alone, but the user only sees the AGGREGATE failure count in
        # the summary, so collapse all still-failing rows to "api_failed".
        new_af <- df[[gnomad_af_col]][mark]
        new_cls <- ifelse(is.na(new_af),
                           "api_failed",
                           ifelse(new_af > 0, "annotated", "absent"))
        df$.gnomAD_class[mark] <- new_cls
      }
  }
  },
  error = function(e) {
    warning(sprintf(
      "ctdna_annotate_population_freq: retry/classification step failed (%s). ",
      conditionMessage(e)),
      "Returning the partially-annotated frame from Stage A/B. ",
      "The caches are intact; you can save them with ctdna_save_popfreq_cache().",
      call. = FALSE)
  })

  # ---- Summary --------------------------------------------------------
  if (verbose) {
    cls <- df$.gnomAD_class
    n_ann  <- sum(cls == "annotated", na.rm = TRUE)
    n_abs  <- sum(cls == "absent",    na.rm = TRUE)
    n_fail <- sum(cls %in% c("api_failed",
                              "api_failed_short","api_failed_cnv","api_failed_sv"),
                   na.rm = TRUE)
    n_none <- sum(cls == "no_source", na.rm = TRUE)
    n_nog  <- sum(cls == "no_gene",   na.rm = TRUE)
    n_skip <- nrow(df) - n_ann - n_abs - n_fail - n_none - n_nog
    message(sprintf(
"ctdna_annotate_population_freq summary:
  Annotated with AF from gnomAD     : %d
  Treated as AF=0 (absent from gnomAD): %d
  Kept NA (API failure, retry-eligible): %d
  Kept NA (no gnomAD source for type)  : %d
  Kept NA (missing gene symbol)        : %d
  Already-populated (skipped)          : %d",
      n_ann, n_abs, n_fail, n_none, n_nog, n_skip))
  }

  # Remove internal columns (defensive: tolerate any missing column,
  # since a tryCatch above may have aborted before all were created).
  for (col in c(".gnomAD_class",".chrom_norm",".pos_norm",".ref",".alt"))
    if (col %in% names(df)) df[[col]] <- NULL
  df
}
