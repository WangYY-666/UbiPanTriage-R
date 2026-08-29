# ============================================================================
# annotation.R - gene feature annotation (gene ubiquitin feature annotation)
# Ubiquitin functional tier scores + gene type classification (E1/E2/E3/DUB/UBD/ULD/Other).
# ============================================================================

#' Ubiquitin functional tier scores (ubiquitin functional hierarchy scores)
#'
#' Coarse biological hierarchy of the ubiquitination cascade. E3 ligases
#' confer substrate specificity and are the most experimentally tractable
#' targets; DUBs reverse ubiquitination; E2/E1 are the catalytic chain;
#' UBD/ULD are readers/modifiers rather than core enzymes. Values are the
#' default ubiquitin dimension score and are constant across cancers.
#' @format named numeric vector
#' @source IUCCD2.0 classification, tier scores re-designed per biological hierarchy.
"ubi_function_scores" <- c(
  "E3"    = 1.00,   # ubiquitin ligase - substrate specificity
  "DUB"   = 0.90,   # deubiquitinating enzyme - removes ubiquitin
  "E2"    = 0.80,   # ubiquitin-conjugating enzyme - catalytic carrier
  "E1"    = 0.75,   # ubiquitin-activating enzyme - cascade initiator
  "ULD"   = 0.70,   # ubiquitin-like domain - modifier/regulator
  "UBD"   = 0.65,   # ubiquitin-binding domain - reader
  "Other" = 0.60    # other ubiquitin-related gene
)

#' Map a primary classification to a coarse ubiquitin type (E1/E2/E3/DUB/UBD/ULD/Other)
#' @noRd
.ubi_type_from_primary <- function(primary) {
  out <- rep(NA_character_, length(primary))
  p <- as.character(primary)
  out[grepl("^E1", p)] <- "E1"
  out[grepl("^E2", p)] <- "E2"
  out[grepl("^E3", p)] <- "E3"
  out[grepl("^DUB", p)] <- "DUB"
  out[grepl("^UBD", p)] <- "UBD"
  out[grepl("^ULD", p)] <- "ULD"
  out[is.na(out) & !is.na(p) & nzchar(p) & p != "Unknown"] <- "Other"
  out[is.na(out) & !is.na(p) & p == "Unknown"] <- "Other"
  out
}

#' Full display name of a ubiquitin type (full name for figures/tables)
#' @noRd
.ubi_type_full <- function(type) {
  tbl <- c(
    E1    = "Ubiquitin-activating enzyme (E1)",
    E2    = "Ubiquitin-conjugating enzyme (E2)",
    E3    = "Ubiquitin ligase (E3)",
    DUB   = "Deubiquitinating enzyme (DUB)",
    UBD   = "Ubiquitin-binding domain (UBD)",
    ULD   = "Ubiquitin-like domain (ULD)",
    Other = "Other ubiquitin-related gene")
  unname(tbl[type])
}

#' Tier score of a ubiquitin type (tier score of one type)
#' @noRd
.ubi_tier_score <- function(type) unname(ubi_function_scores[type])

#' Extract primary class from a family string (extract primary class from Family)
#' @noRd
.extract_primary <- function(fam) {
  if (is.na(fam) || fam == "") return(NA_character_)
  cls <- strsplit(fam, "; ", fixed = TRUE)[[1]]
  prio <- c("E1", "E2", "E3 activity", "E3 adaptor", "DUB", "UBD", "ULD")
  for (p in prio) {
    hit <- cls[startsWith(cls, p)]
    if (length(hit) > 0) return(hit[1])
  }
  cls[1]
}

#' Legacy fine-grained family score (kept for cache compatibility)
#' @noRd
.family_score <- function(cls) {
  if (is.na(cls) || cls == "") return(0.3)
  ty <- .ubi_type_from_primary(cls)
  if (is.na(ty)) 0.3 else unname(ubi_function_scores[ty])
}

#' Annotate genes with ubiquitin feature scores (annotate genes with ubiquitin feature)
#'
#' Joins user genes with the bundled IUCCD2.0-based annotation: coarse
#' ubiquitin type (E1/E2/E3/DUB/UBD/ULD/Other), full display name, family
#' classification and the tier-based feature score (0-1). Genes without
#' annotation get \code{NA} scores (custom genes are not penalized).
#' @param cache cancer cache (cancer cache, must contain an annotation table)
#' @param genes gene symbols (gene symbols)
#' @return data.frame with columns gene, family, primary, ubi_type,
#'   ubi_type_full, ubi_score, feature_score
#' @examples
#' \donttest{
#' luad <- load_cancer_data("LUAD")
#' head(annotate_genes(luad, c("MDM2", "TRIM44", "EGFR")))
#' }
#' @export
annotate_genes <- function(cache, genes) {
  genes <- unique(trimws(as.character(genes)))
  ann  <- cache$annotation
  out  <- data.frame(gene = genes, stringsAsFactors = FALSE)
  if (is.null(ann)) {
    out$family <- NA_character_
    out$primary <- NA_character_
    out$ubi_type <- NA_character_
    out$ubi_type_full <- NA_character_
    out$ubi_score <- NA_real_
    out$feature_score <- NA_real_
    return(out)
  }
  m <- match(out$gene, ann$gene)
  out$family         <- ann$family[m]
  out$primary        <- ann$primary[m]
  if ("ubi_type" %in% colnames(ann)) {
    out$ubi_type      <- ann$ubi_type[m]
    out$ubi_type_full <- ann$ubi_type_full[m]
    out$ubi_score     <- as.numeric(ann$ubi_score[m])
  } else {
    # fallback: recompute from primary (older caches)
    ty <- .ubi_type_from_primary(out$primary)
    out$ubi_type      <- ty
    out$ubi_type_full <- .ubi_type_full(ty)
    out$ubi_score     <- .ubi_tier_score(ty)
  }
  out$feature_score <- out$ubi_score
  out
}

#' Feature score for a set of genes (feature score, NA kept for custom genes)
#' @noRd
.feature_scores <- function(annot) annot$feature_score
