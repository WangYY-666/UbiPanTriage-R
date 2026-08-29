# ============================================================================
# utils.R - shared helper functions (通用工具函数)
# ============================================================================

#' Default value for NULL
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Linear mapping from [x_min, x_max] to [0, 1] (线性映射到 0-1)
#' @noRd
linear_map <- function(x, x_min, x_max) {
  (x - x_min) / (x_max - x_min)
}

#' Map ranks to scores in [0.5, 1] (排名映射为 0.5-1 分)
#' @param ranks integer ranks, 1 = best (整数排名，1 为最优)
#' @noRd
rank_to_score <- function(ranks) {
  if (length(ranks) == 0 || all(is.na(ranks))) return(rep(NA_real_, length(ranks)))
  (ranks - 1) / (max(ranks, na.rm = TRUE) - 1) * 0.5 + 0.5
}

#' Percentile rank, 1 = best (百分位排名，1 为最优)
#' @param x numeric scores (数值向量)
#' @param na.rm drop NAs (是否忽略缺失值)
#' @noRd
pct_rank <- function(x, na.rm = TRUE) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  r <- rank(x, ties.method = "min", na.last = "keep")
  n <- sum(!is.na(x))
  r / n
}

#' Harmonic mean with stability epsilon (带稳定项的调和平均数)
#' @noRd
harm_mean <- function(x, eps = 0.001) {
  x <- x[!is.na(x) & x > 0]
  if (length(x) == 0) return(NA_real_)
  length(x) / sum(1 / (x + eps))
}

#' Normalize vector to [0,1] (向量归一化到 0-1)
#' @noRd
normalize01 <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (diff(r) == 0) return(rep(0, length(x)))
  (x - r[1]) / diff(r)
}

#' Check and clean user-provided gene symbols (校验并清洗基因名)
#' @param genes character vector of gene symbols (基因符号)
#' @param universe available gene symbols in the cache (缓存中可用的基因)
#' @return list with found / missing genes (返回找到与缺失的基因列表)
#' @export
check_genes <- function(genes, universe) {
  genes <- unique(trimws(as.character(genes)))
  genes <- genes[!is.na(genes) & genes != ""]
  found   <- intersect(genes, universe)
  missing <- setdiff(genes, universe)
  list(genes = found, missing = missing)
}

#' Fast spearman correlation matrix between two matrices (快速 Spearman 相关)
#' @param x matrix, rows = genes (基因矩阵)
#' @param y matrix, rows = samples (特征矩阵，样本行为主)
#' @return correlation matrix genes x traits (基因 x 特征 相关矩阵)
#' @noRd
spearman_matrix <- function(x, y) {
  x <- x[apply(x, 1, function(v) stats::sd(v) > 0), , drop = FALSE]
  y <- y[, apply(y, 2, function(v) stats::sd(v) > 0), drop = FALSE]
  if (nrow(x) == 0 || ncol(y) == 0) return(NULL)
  stats::cor(t(x), y, method = "spearman")
}

#' Asymptotic p-value of spearman correlation (Spearman 相关渐近 p 值)
#' @param r correlation value (相关系数)
#' @param n sample size (样本量)
#' @noRd
cor_pvalue <- function(r, n) {
  2 * stats::pt(-abs(r) * sqrt((n - 2) / (1 - r^2 + 1e-12)), df = n - 2)
}

#' Robust row-weighted mean (稳健的行加权均值)
#' @noRd
row_wmean <- function(mat, w) {
  w <- w[names(w) %in% colnames(mat)]
  m <- mat[, names(w), drop = FALSE]
  if (ncol(m) == 0) return(rep(0, nrow(mat)))
  rowSums(abs(m) * rep(w, each = nrow(m)), na.rm = TRUE) / sum(w)
}

#' Pick the trait with the strongest absolute correlation (取相关性最强的特征)
#' @noRd
top_trait <- function(v) {
  if (all(is.na(v)) || length(v) == 0) return(list(feature = NA_character_, cor = NA_real_))
  i <- which.max(abs(v))
  list(feature = names(v)[i], cor = unname(v[i]))
}
#' Map scores to percentile ranks in [0, 1] (分数映射为 0-1 百分位，越高越接近 1)
#'
#' Best score -> 1, worst -> 1/N, NAs are kept. Used to put the basic /
#' immune / metabolic dimensions on a comparable 0-1 scale before
#' cross-cancer aggregation.
#' 将分数转换为 0-1 百分位（最优为 1，最差为 1/N），用于跨癌种汇总前
#' 将基础 / 免疫 / 代谢维度统一到可比的 0-1 尺度。
#' @param x numeric vector of scores (分数向量)
#' @return numeric vector of percentile ranks (百分位向量)
#' @export
rank_to_pct <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  r <- rank(-x, ties.method = "min", na.last = "keep")
  n <- sum(!is.na(x))
  (n - r + 1) / n
}