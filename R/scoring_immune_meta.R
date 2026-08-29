# ============================================================================
# scoring_immune_meta.R - immune & metabolic dimension scoring
# Immune & metabolic dimension scoring (免疫与代谢维度评分)
# ============================================================================

#' Tumor samples (patient-deduplicated) shared by expression and trait tables
#' 取表达矩阵与特征矩阵共有的肿瘤样本（同一患者去重），样本量不足时返回 NULL。
#' @param cache cancer cache (癌种缓存)
#' @param trait_df trait data.frame with an ID column (含 ID 列的特征表)
#' @noRd
.tumor_common <- function(cache, trait_df) {
  si <- cache$sample_info
  tum <- si$sample[si$type == "Tumor"]
  pid <- substr(tum, 1, 12)
  tum_keep <- tum[!duplicated(pid)]
  pids <- substr(tum_keep, 1, 12)
  tdf <- trait_df
  tdf$pid <- substr(tdf$ID, 1, 12)
  tdf <- tdf[!duplicated(tdf$pid), , drop = FALSE]
  common <- intersect(pids, tdf$pid)
  if (length(common) < 10) return(NULL)
  x <- cache$expr_tpm[, tum_keep[match(common, pids)], drop = FALSE]
  y <- tdf[match(common, tdf$pid), , drop = FALSE]
  list(x = x, y = y)
}

#' Gene-trait spearman correlation matrix (基因-特征 Spearman 相关矩阵)
#' Uses pre-computed values for background genes and computes on the fly
#' for user-supplied genes. 背景基因直接读取预计算矩阵，新基因即时计算。
#' @param cache cancer cache (癌种缓存)
#' @param genes gene symbols (基因符号)
#' @param trait_df trait data.frame with an ID column (含 ID 列的特征表)
#' @param pre_mat pre-computed correlation matrix (预计算相关矩阵)
#' @noRd
.gene_trait_cor <- function(cache, genes, trait_df, pre_mat) {
  genes <- unique(as.character(genes))
  if (is.null(pre_mat) && is.null(trait_df)) return(NULL)
  out <- NULL
  if (!is.null(pre_mat)) {
    pre <- pre_mat[intersect(genes, rownames(pre_mat)), , drop = FALSE]
    if (nrow(pre) > 0) out <- pre
  }
  new_genes <- setdiff(genes, if (is.null(out)) character(0) else rownames(out))
  if (length(new_genes) > 0 && !is.null(trait_df)) {
    tc <- .tumor_common(cache, trait_df)
    if (!is.null(tc)) {
      present <- intersect(new_genes, rownames(tc$x))
      if (length(present) > 0) {
        x <- log2(tc$x[present, , drop = FALSE] + 1)
        y <- as.matrix(tc$y[, setdiff(colnames(tc$y), c("ID", "pid")), drop = FALSE])
        cm <- spearman_matrix(x, y)
        if (!is.null(cm)) {
          rownames(cm) <- rownames(x)
          out <- if (is.null(out)) cm else rbind(out, cm)
        }
      }
    }
  }
  if (is.null(out)) {
    out <- matrix(NA_real_, nrow = length(genes), ncol = 0,
                  dimnames = list(genes, character(0)))
    return(out)
  }
  miss <- setdiff(genes, rownames(out))
  if (length(miss) > 0) {
    m <- matrix(NA_real_, nrow = length(miss), ncol = ncol(out),
                dimnames = list(miss, colnames(out)))
    out <- rbind(out, m)
  }
  out[genes, , drop = FALSE]
}

#' Immune infiltration correlation (免疫浸润相关矩阵)
#' @param cache cancer cache (癌种缓存)
#' @param genes gene symbols (基因符号)
#' @param method immune method: timer / cibersort / mcp / ssgsea / xcell
#' @return matrix genes x immune traits (基因 x 免疫特征的相关系数矩阵)
#' @export
immune_cor <- function(cache, genes, method = "timer") {
  trait_df <- cache$immune[[method]]
  pre_mat  <- cache$pre$cor_immune[[method]]
  .gene_trait_cor(cache, genes, trait_df, pre_mat)
}

#' Metabolic pathway correlation (代谢途径相关矩阵)
#' @param cache cancer cache (癌种缓存)
#' @param genes gene symbols (基因符号)
#' @return matrix genes x pathways (基因 x 代谢途径的相关系数矩阵)
#' @export
metabolic_cor <- function(cache, genes) {
  .gene_trait_cor(cache, genes, cache$metabolism, cache$pre$cor_meta)
}

#' Immune dimension score (免疫维度评分)
#'
#' Weighted mean of |spearman r| between the gene and immune infiltration
#' traits; default weights follow the LUAD example (CD8/Macrophage/DC
#' emphasized). Top trait and signed correlation are also reported.
#' 免疫分 = 基因与免疫浸润特征加权 |Spearman r| 的加权均值；权重默认按 LUAD 示例侧重 CD8 / Macrophage / DC，并报告主要相关细胞与相关方向。
#' @param cache cancer cache (癌种缓存)
#' @param genes gene symbols (基因符号)
#' @param method immune method (免疫浸润方法)
#' @param weights named weights per trait (各特征的命名权重；NULL 表示组内等权)
#' @return data.frame with Immune_Score, Top_immune_cell, Cor_value1 and
#'   per-trait correlations (免疫分、主要相关细胞、相关系数及各特征相关值)
#' @export
calc_immune_score <- function(cache, genes, method = "timer", weights = NULL) {
  genes <- unique(as.character(genes))
  cm <- immune_cor(cache, genes, method)
  if (is.null(cm) || ncol(cm) == 0) {
    return(data.frame(gene = genes, Immune_Score = 0,
                      Top_immune_cell = NA_character_, Cor_value1 = NA_real_,
                      stringsAsFactors = FALSE))
  }
  w <- weights %||% default_immune_weights(method)
  if (is.null(w)) w <- stats::setNames(rep(1, ncol(cm)), colnames(cm))
  wfull <- w[colnames(cm)]; wfull[is.na(wfull)] <- 1
  score <- row_wmean(cm, wfull)
  tops <- apply(cm, 1, top_trait)
  out <- data.frame(
    gene = genes,
    Immune_Score = as.numeric(score),
    Top_immune_cell = unlist(lapply(tops, `[[`, "feature")),
    Cor_value1 = as.numeric(unlist(lapply(tops, `[[`, "cor"))),
    stringsAsFactors = FALSE)
  colnames(cm) <- paste0("Imm_", colnames(cm))
  cbind(out, as.data.frame(cm, check.names = FALSE))
}

#' Metabolic dimension score (代谢维度评分)
#'
#' Weighted mean of |spearman r| between the gene and 7 GSVA metabolic
#' pathways. 代谢分 = 基因与 7 条 GSVA 代谢途径加权 |Spearman r| 的加权均值。
#' @param cache cancer cache (癌种缓存)
#' @param genes gene symbols (基因符号)
#' @param weights named weights per pathway (各途径的命名权重；NULL 表示等权)
#' @return data.frame with Metabolic_Score, Top_metabolic_trait, Cor_value2
#'   and per-pathway correlations (代谢分、主要相关途径、相关系数及各途径相关值)
#' @export
calc_metabolic_score <- function(cache, genes, weights = NULL) {
  genes <- unique(as.character(genes))
  cm <- metabolic_cor(cache, genes)
  if (is.null(cm) || ncol(cm) == 0) {
    return(data.frame(gene = genes, Metabolic_Score = 0,
                      Top_metabolic_trait = NA_character_, Cor_value2 = NA_real_,
                      stringsAsFactors = FALSE))
  }
  w <- weights %||% default_metabolic_weights()
  if (is.null(w)) w <- stats::setNames(rep(1, ncol(cm)), colnames(cm))
  wfull <- w[colnames(cm)]; wfull[is.na(wfull)] <- 1
  score <- row_wmean(cm, wfull)
  tops <- apply(cm, 1, top_trait)
  out <- data.frame(
    gene = genes,
    Metabolic_Score = as.numeric(score),
    Top_metabolic_trait = unlist(lapply(tops, `[[`, "feature")),
    Cor_value2 = as.numeric(unlist(lapply(tops, `[[`, "cor"))),
    stringsAsFactors = FALSE)
  colnames(cm) <- paste0("Meta_", colnames(cm))
  cbind(out, as.data.frame(cm, check.names = FALSE))
}
