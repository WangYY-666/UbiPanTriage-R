# ============================================================================
# scoring_combine.R - combine four dimensions & rank (four-dimension scoring: ubiquitin/basic/immune/metabolic)
# ============================================================================

#' Default parameters for the basic dimension (default basic parameters)
#' @export
default_basic_params <- function() list(
  weights = c(expression = 1, diff = 1, survive = 1, roc = 1),
  expression = list(lo = 5, peak_lo = 10, peak_hi = 100, hi = 1000),
  diff = list(methods = c("limma", "edgeR", "deseq2", "pair"),
              method_weights = c(limma = 1, edgeR = 1, deseq2 = 1, pair = 0.5),
              padj_cutoff = NULL),
  survive = list(p_cutoff = 0.05),
  roc = list(timepoints = c(1, 3, 5))
)

#' Default parameters for the immune dimension (default immune parameters)
#' @export
default_immune_params <- function() list(method = "timer", weights = NULL)

#' Default parameters for the metabolic dimension (default metabolic parameters)
#' @export
default_metabolic_params <- function() list(weights = NULL)

#' Default parameters for score combination (default combine parameters)
#'
#' Four dimensions: ubiquitin / basic / immune / metabolic, combined with a
#' weighted mean (equal weights 0.25 x 4 by default) or optionally a harmonic
#' mean. The ubiquitin dimension uses the functional tier score and can be
#' merged with an optional activity-evidence bonus.
#' @export
default_combine_params <- function() list(
  dimensions = c("ubiquitin", "basic", "immune", "metabolic"),
  method = "weighted",                          # weighted or harmonic
  dimension_weights = c(ubiquitin = 0.25, basic = 0.25,
                        immune = 0.25, metabolic = 0.25),
  ubi_activity = FALSE,                         # add set2/set3 activity bonus
  dual_bonus = FALSE,                           # optional dual-association bonus
  rank_in = "universe"                          # universe or query
)

#' Validate dimension weights (validate dimension weights)
#'
#' Dimension weights must be numeric, within 0.1-1.0, with one decimal place,
#' and sum to 1 (as required by the interactive interface).
#' @param w named numeric weights for the four dimensions
#' @param dims dimension names
#' @return list(ok, message, weights)
#' @export
check_dim_weights <- function(w, dims = c("ubiquitin", "basic", "immune", "metabolic")) {
  if (is.null(w)) return(list(ok = FALSE, message = "weights are required"))
  if (length(dims) == 0) return(list(ok = FALSE, message = "no dimensions selected"))
  w <- w[dims]
  if (any(is.na(w)) || length(w) != length(dims))
    return(list(ok = FALSE, message = "one weight per dimension is required"))
  if (any(!is.numeric(w)))
    return(list(ok = FALSE, message = "weights must be numeric"))
  if (any(w < 0.1 | w > 1.0))
    return(list(ok = FALSE, message = "weights must be within 0.1 - 1.0"))
  if (any(abs(w - round(w, 1)) > 1e-6))
    return(list(ok = FALSE, message = "weights must have one decimal place"))
  if (abs(sum(w) - 1) > 1e-6)
    return(list(ok = FALSE, message = "weights must sum to 1"))
  list(ok = TRUE, message = "", weights = w)
}

#' Percentile of query scores against a background (percentile against background)
#' @noRd
.pct_vs_bg <- function(q, bg) {
  if (is.null(bg)) return(pct_rank(q))
  n <- sum(!is.na(c(bg, q)))
  rank(c(bg, q), ties.method = "min", na.last = "keep")[seq_along(q)] / n
}

#' Merge params with defaults (merge params with defaults)
#' @noRd
.merge_params <- function(p, defaults) {
  if (is.null(p)) return(defaults)
  for (nm in names(defaults)) if (is.null(p[[nm]])) p[[nm]] <- defaults[[nm]]
  p
}

#' Ubiquitin dimension score (ubiquitin dimension score)
#'
#' Default: the functional tier score of the gene (E3=1.00, DUB=0.90, E2=0.80,
#' E1=0.75, ULD=0.70, UBD=0.65, Other=0.60), constant across cancers.
#' With \code{use_activity = TRUE}, genes also annotated in the IUCCD2.0
#' activity evidence sets get a small bonus. Non-ubiquitin genes get NA
#' (dimension excluded for those genes).
#' @param cache cancer cache
#' @param genes gene symbols
#' @param use_activity add activity-evidence bonus (default FALSE)
#' @return data.frame(gene, Ubi_Score, Ubi_Type, Ubi_Type_Full)
#' @export
calc_ubi_score <- function(cache, genes, use_activity = FALSE) {
  genes <- unique(as.character(genes))
  ann <- annotate_genes(cache, genes)
  score <- ann$ubi_score
  if (use_activity && all(c("in_set2", "in_set3") %in% colnames(cache$annotation))) {
    a <- cache$annotation
    b <- ifelse(a$gene %in% ann$gene & a$in_set2 & a$in_set3, 0.15,
         ifelse(a$gene %in% ann$gene & (a$in_set2 | a$in_set3), 0.08, 0))
    m <- match(genes, a$gene)
    score <- score + b[m]
    score <- pmin(score, 1)
  }
  data.frame(gene = genes, Ubi_Score = as.numeric(score),
             Ubi_Type = ann$ubi_type, Ubi_Type_Full = ann$ubi_type_full,
             stringsAsFactors = FALSE)
}

#' Four-dimensional gene scoring for one cancer (four-dimension gene scoring)
#'
#' The main scoring entry point. Computes Ubiquitin / Basic / Immune /
#' Metabolic dimension scores for the queried genes and combines them into a
#' single 0-1 score with ranking. Weights are re-normalized per gene over the
#' dimensions that are available (e.g. ubiquitin is excluded for custom genes,
#' basic sub-items are re-weighted when a cancer lacks diff/normal tissue).
#' @param cache cancer cache (see \code{\link{load_cancer_data}})
#' @param genes gene symbols
#' @param basic_params parameters for the basic dimension
#' @param immune_params parameters for the immune dimension
#' @param metabolic_params parameters for the metabolic dimension
#' @param combine_params parameters for combining dimensions
#' @return list with \code{scores} (main table), \code{immune},
#'   \code{metabolic}, \code{annotation}, \code{params}
#' @examples
#' \donttest{
#' luad <- load_cancer_data("LUAD")
#' r <- score_genes(luad, c("MDM2", "TRIM44", "UBE2C"))
#' print(r$scores[, c("Gene", "Ubi_Score", "Basic_Score", "Immune_Score", "Metabolic_Score")])
#' }
#' @export
score_genes <- function(cache, genes,
                        basic_params = default_basic_params(),
                        immune_params = default_immune_params(),
                        metabolic_params = default_metabolic_params(),
                        combine_params = default_combine_params()) {
  basic_params   <- .merge_params(basic_params,   default_basic_params())
  immune_params  <- .merge_params(immune_params,  default_immune_params())
  metabolic_params <- .merge_params(metabolic_params, default_metabolic_params())
  combine_params <- .merge_params(combine_params, default_combine_params())

  genes <- unique(trimws(as.character(genes)))
  genes <- genes[!is.na(genes) & genes != ""]
  avail <- rownames(cache$expr_tpm)
  ck <- check_genes(genes, avail)
  if (length(ck$genes) == 0) stop("none of the input genes were found in the cache")
  genes <- ck$genes

  dims <- intersect(combine_params$dimensions, c("ubiquitin", "basic", "immune", "metabolic"))
  if (length(dims) == 0) stop("at least one dimension must be selected")

  # ---- dimension scores for query genes ----
  ubi <- if ("ubiquitin" %in% dims)
    calc_ubi_score(cache, genes, use_activity = isTRUE(combine_params$ubi_activity)) else NULL
  basic <- if ("basic" %in% dims)
    calc_basic_score(cache, genes,
                     weights = basic_params$weights,
                     expression = basic_params$expression,
                     diff = basic_params$diff,
                     survive = basic_params$survive,
                     roc = basic_params$roc) else NULL
  imm <- if ("immune" %in% dims)
    calc_immune_score(cache, genes, immune_params$method, immune_params$weights) else NULL
  meta <- if ("metabolic" %in% dims)
    calc_metabolic_score(cache, genes, metabolic_params$weights) else NULL

  # ---- background scores for percentile denominators ----
  bg <- if (combine_params$rank_in == "universe") cache$universe else NULL
  bg_ubi <- bg_basic <- bg_imm <- bg_meta <- NULL
  if (!is.null(bg)) {
    if ("ubiquitin" %in% dims) bg_ubi <- calc_ubi_score(cache, bg,
        use_activity = isTRUE(combine_params$ubi_activity))$Ubi_Score
    if ("basic" %in% dims) bg_basic <- calc_basic_score(cache, bg,
        weights = basic_params$weights, expression = basic_params$expression,
        diff = basic_params$diff, survive = basic_params$survive,
        roc = basic_params$roc)$Basic_Score
    if ("immune" %in% dims) bg_imm <- calc_immune_score(cache, bg,
        immune_params$method, immune_params$weights)$Immune_Score
    if ("metabolic" %in% dims) bg_meta <- calc_metabolic_score(cache, bg,
        metabolic_params$weights)$Metabolic_Score
  }

  out <- data.frame(Gene = genes, stringsAsFactors = FALSE)
  pct_cols <- character(0)
  if ("ubiquitin" %in% dims) {
    out$Ubi_Score <- ubi$Ubi_Score[match(genes, ubi$gene)]
    out$Ubi_Type <- ubi$Ubi_Type[match(genes, ubi$gene)]
    out$Ubi_Type_Full <- ubi$Ubi_Type_Full[match(genes, ubi$gene)]
    out$Ubi_Pct <- .pct_vs_bg(out$Ubi_Score, bg_ubi)
    pct_cols <- c(pct_cols, "Ubi_Pct")
  }
  if ("basic" %in% dims) {
    out$Basic_Score <- basic$Basic_Score[match(genes, basic$gene)]
    out$Missing <- basic$Missing[match(genes, basic$gene)]
    out$Basic_Pct <- .pct_vs_bg(out$Basic_Score, bg_basic)
    pct_cols <- c(pct_cols, "Basic_Pct")
  }
  if ("immune" %in% dims) {
    out$Immune_Score <- imm$Immune_Score[match(genes, imm$gene)]
    out$Top_immune_cell <- imm$Top_immune_cell[match(genes, imm$gene)]
    out$Cor_value1 <- imm$Cor_value1[match(genes, imm$gene)]
    out$Immune_Pct <- .pct_vs_bg(out$Immune_Score, bg_imm)
    pct_cols <- c(pct_cols, "Immune_Pct")
  }
  if ("metabolic" %in% dims) {
    out$Metabolic_Score <- meta$Metabolic_Score[match(genes, meta$gene)]
    out$Top_metabolic_trait <- meta$Top_metabolic_trait[match(genes, meta$gene)]
    out$Cor_value2 <- meta$Cor_value2[match(genes, meta$gene)]
    out$Metabolic_Pct <- .pct_vs_bg(out$Metabolic_Score, bg_meta)
    pct_cols <- c(pct_cols, "Metabolic_Pct")
  }
  # ubiquitin missing flag (non-ubiquitin gene -> "ubi" in Missing)
  if ("ubiquitin" %in% dims && !is.null(out$Missing)) {
    i <- is.na(out$Ubi_Score)
    out$Missing[i] <- ifelse(is.na(out$Missing[i]), "ubi", paste(out$Missing[i], "ubi", sep = ","))
  }

  # ---- combined score (per-gene weight re-normalization over available dims) ----
  pc <- as.matrix(out[, pct_cols, drop = FALSE])
  dw <- combine_params$dimension_weights[tolower(sub("_Pct$", "", pct_cols))]
  dw[is.na(dw)] <- 1
  if (combine_params$method == "weighted") {
    wm <- matrix(rep(dw, each = nrow(pc)), nrow = nrow(pc))
    wm[is.na(pc)] <- 0
    wt <- rowSums(wm)
    combined <- rowSums(pc * wm, na.rm = TRUE) / ifelse(wt > 0, wt, 1)
  } else {
    combined <- apply(pc, 1, harm_mean, eps = 0.001)
  }
  if (combine_params$dual_bonus && ncol(pc) >= 2) {
    allmin <- apply(pc, 1, min, na.rm = TRUE)
    combined <- combined + ifelse(allmin > 0.7, 0.2, ifelse(allmin > 0.5, 0.1, 0))
  }
  out$Combined_Score <- combined
  out$Combined_Rank <- rank(-combined, ties.method = "min")

  # ---- attach annotation & detail ----
  ann <- annotate_genes(cache, genes)
  colnames(ann)[colnames(ann) == "gene"] <- "Gene"
  ann <- ann[, intersect(c("Gene", "family", "primary", "ubi_type",
                           "ubi_type_full", "ubi_score"), colnames(ann)), drop = FALSE]
  out <- merge(ann, out, by = "Gene", all.y = TRUE, sort = FALSE)
  if (!is.null(basic)) {
    detail_basic <- basic[match(genes, basic$gene),
                          c("gene", "Expression_Score", "Diff_Score", "Survive_Score",
                            "ROC_Score", "HR", "HR_p", "direction", "AUC_surv", "AUC_diag")]
    colnames(detail_basic)[1] <- "Gene"
    out <- merge(out, detail_basic, by = "Gene", all.x = TRUE, sort = FALSE)
  }
  list(scores = out,
       immune = imm,
       metabolic = meta,
       annotation = ann,
       params = list(basic = basic_params, immune = immune_params,
                     metabolic = metabolic_params, combine = combine_params))
}

#' Multi-cancer scoring (multi-cancer scoring)
#'
#' Runs \code{\link{score_genes}} for several cancer caches and returns a long
#' gene x cancer table plus a per-gene pan-cancer summary.
#' @param caches named list of cancer caches
#' @param genes gene symbols
#' @inheritParams score_genes
#' @return list with \code{scores_long} (gene x cancer) and \code{summary}
#' @export
score_genes_multi <- function(caches, genes,
                              basic_params = default_basic_params(),
                              immune_params = default_immune_params(),
                              metabolic_params = default_metabolic_params(),
                              combine_params = default_combine_params()) {
  per <- lapply(names(caches), function(ca) {
    s <- score_genes(caches[[ca]], genes, basic_params, immune_params,
                     metabolic_params, combine_params)
    s$scores$Cancer <- ca
    s
  })
  names(per) <- names(caches)
  keep <- c("Gene", "Cancer", "Ubi_Score", "Basic_Score", "Immune_Score",
            "Metabolic_Score", "Ubi_Pct", "Basic_Pct", "Immune_Pct",
            "Metabolic_Pct", "Combined_Score", "Missing")
  long <- do.call(rbind, lapply(per, function(s) s$scores[, intersect(keep, colnames(s$scores))]))
  rownames(long) <- NULL
  score_cols <- intersect(c("Ubi_Score", "Basic_Score", "Immune_Score",
                            "Metabolic_Score", "Combined_Score"), colnames(long))
  pct_cols <- intersect(c("Ubi_Pct", "Basic_Pct", "Immune_Pct", "Metabolic_Pct"),
                        colnames(long))
  summ <- aggregate(long[, score_cols, drop = FALSE],
                    by = list(Gene = long$Gene), mean, na.rm = TRUE)
  if (length(pct_cols) > 0) {
    sp <- aggregate(long[, pct_cols, drop = FALSE],
                    by = list(Gene = long$Gene), mean, na.rm = TRUE)
    summ <- merge(summ, sp, by = "Gene")
  }
  summ$Pan_Rank <- rank(-summ$Combined_Score, ties.method = "min")
  summ <- summ[order(summ$Pan_Rank), , drop = FALSE]
  list(scores_long = long, summary = summ, per_cancer = per)
}
