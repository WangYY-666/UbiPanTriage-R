# ============================================================================
# cancer_data.R - local cancer cache access (本地癌种缓存访问)
# ============================================================================
#' Method-specific suffix used in immune table column names (方法列名后缀)
#' @noRd
.method_suffix <- function(method) {
  switch(method,
    timer     = "_TIMER",
    cibersort = "_CIBERSORT",
    mcp       = "_MCPcounter",
    xcell     = "_xCell",
    ssgsea    = "")
}
#' Strip method suffix from a trait name (去除方法后缀)
#' @noRd
.strip_suffix <- function(x, method) {
  sfx <- .method_suffix(method)
  if (sfx == "") x else sub(paste0(sfx, "$"), "", x)
}
#' Post-process a loaded cache: clean trait names (载入后清洗特征名)
#' @noRd
.tidy_cache <- function(cache) {
  # recompute ubiquitin annotation with the current tier scheme (recompute ubi annotation)
  if (!is.null(cache$annotation) && "primary" %in% colnames(cache$annotation)) {
    a <- cache$annotation
    if (!"ubi_type" %in% colnames(a)) {
      ty <- .ubi_type_from_primary(a$primary)
      a$ubi_type <- ty
      a$ubi_type_full <- .ubi_type_full(ty)
      a$ubi_score <- .ubi_tier_score(ty)
    }
    cache$annotation <- a
  }
  for (m in names(cache$immune)) {
    sfx <- .method_suffix(m)
    if (sfx != "") {
      colnames(cache$immune[[m]]) <- sub(paste0(sfx, "$"), "", colnames(cache$immune[[m]]))
      cm <- cache$pre$cor_immune[[m]]
      if (!is.null(cm)) {
        colnames(cm) <- sub(paste0(sfx, "$"), "", colnames(cm))
        cache$pre$cor_immune[[m]] <- cm
      }
    }
    # drop meta columns of CIBERSORT output (CIBERSORT 输出的元数据列)
    drop <- grepl("^(P.value|Correlation|RMSE)$", colnames(cache$immune[[m]]))
    if (any(drop)) cache$immune[[m]] <- cache$immune[[m]][, !drop, drop = FALSE]
    cm <- cache$pre$cor_immune[[m]]
    if (!is.null(cm)) {
      drop <- grepl("^(P.value|Correlation|RMSE)$", colnames(cm))
      if (any(drop)) cache$pre$cor_immune[[m]] <- cm[, !drop, drop = FALSE]
    }
    # xCell overall scores (immune/stroma/microenvironment) are not cell-type
    # traits and are excluded from cell scoring (xCell 整体评分不参与免疫细胞评分)
    if (m == "xcell") {
      xdrop <- grepl("^(ImmuneScore|StromaScore|MicroenvironmentScore)$",
                     colnames(cache$immune[[m]]))
      if (any(xdrop)) cache$immune[[m]] <- cache$immune[[m]][, !xdrop, drop = FALSE]
      cmx <- cache$pre$cor_immune[[m]]
      if (!is.null(cmx)) {
        xdrop <- grepl("^(ImmuneScore|StromaScore|MicroenvironmentScore)$", colnames(cmx))
        if (any(xdrop)) cache$pre$cor_immune[[m]] <- cmx[, !xdrop, drop = FALSE]
      }
    }
  }
  cache
}
#' List available cancer types (列出可用癌种)
#'
#' Scans \code{data_dir} for pre-computed cancer cache files
#' (\code{<CANCER>_cache.rds}). Each cache is built by
#' \code{data-raw/build_luad_cache.R} (or its pan-cancer equivalents).
#' 扫描本地缓存目录，列出已构建好数据的癌种。
#'
#' @param data_dir directory containing the cache files; defaults to the
#'   bundled \code{inst/extdata} folder of the package.
#'   缓存目录，默认使用包内置的 extdata 目录。
#' @return character vector of TCGA cancer codes (TCGA 癌种缩写).
#' @examples
#' available_cancers()
#' @export
available_cancers <- function(data_dir = system.file("extdata", package = "UbiPanTriage")) {
  fs <- if (dir.exists(data_dir)) {
    list.files(data_dir, pattern = "_cache\\.rds$", full.names = TRUE)
  } else character(0)
  if (length(fs) == 0)
    message("No pre-computed cancer caches found under ", data_dir,
            ". Run UbiPanTriage::download_extdata() once to install the ",
            "pan-cancer data (~4.2 GB).")
  sub("_cache\\.rds$", "", basename(fs))
}
#' Load a cancer cache (载入癌种缓存)
#'
#' Loads the pre-computed local data (expression, clinical, differential
#' expression, immune infiltration, metabolism GSVA, per-gene statistics) for
#' one cancer type. 载入某个癌种的预计算本地数据，用于即时评分。
#'
#' @param cancer TCGA cancer code, e.g. "LUAD" (TCGA 癌种缩写，如 "LUAD").
#' @param data_dir directory containing the cache files (缓存目录).
#' @return a list (cache object) with all data needed for scoring.
#' @examples
#' luad <- load_cancer_data("LUAD")
#' @export
load_cancer_data <- function(cancer = "LUAD",
                             data_dir = system.file("extdata", package = "UbiPanTriage")) {
  f <- file.path(data_dir, paste0(cancer, "_cache.rds"))
  if (!file.exists(f)) stop("cache not found for cancer: ", cancer,
                            " (available: ", paste(available_cancers(data_dir), collapse = ", "), "). ",
                            "If you installed from GitHub, run UbiPanTriage::download_extdata() first.")
  .tidy_cache(readRDS(f))
}
#' Available gene symbols in a cache (缓存中的全部基因)
#' @param cache cancer cache (癌种缓存)
#' @return character vector (基因符号向量)
#' @export
cancer_genes <- function(cache) {
  rownames(cache$expr_tpm)
}
#' Gene expression (log2 TPM+1) for one gene across samples (单基因样本表达)
#' @param cache cancer cache (癌种缓存)
#' @param gene gene symbol (基因符号)
#' @return named numeric vector, names = sample IDs (命名数值向量)
#' @export
gene_expression <- function(cache, gene) {
  if (!gene %in% rownames(cache$expr_tpm)) {
    return(stats::setNames(rep(NA_real_, ncol(cache$expr_tpm)), colnames(cache$expr_tpm)))
  }
  log2(cache$expr_tpm[gene, ] + 1)
}
#' Expression data.frame for plotting: sample, type, expression (绘图用表达数据)
#' @noRd
expr_for_plot <- function(cache, gene) {
  v <- gene_expression(cache, gene)
  d <- data.frame(sample = names(v), expression = as.numeric(v),
                  stringsAsFactors = FALSE)
  d <- merge(d, cache$sample_info[, c("sample", "type")], by = "sample")
  d[d$type %in% c("Tumor", "Normal"), , drop = FALSE]
}
#' Immune trait names for a method (某免疫浸润方法下的细胞特征名)
#' @export
immune_traits <- function(cache, method = "timer") {
  d <- cache$immune[[method]]
  if (is.null(d)) return(character(0))
  .strip_suffix(setdiff(colnames(d), "ID"), method)
}
#' Metabolic pathway names (代谢途径名称)
#' @export
metabolic_traits <- function(cache) {
  setdiff(colnames(cache$metabolism), "ID")
}
#' Default immune weights per method (各免疫方法的默认权重)
#' @export
default_immune_weights <- function(method = "timer") {
  # NULL -> all traits weighted equally within every method (NULL = equal weight inside each method)
  NULL
}
#' Default metabolic weights (代谢途径默认权重)
#' @export
default_metabolic_weights <- function() {
  # NULL -> all 7 pathways weighted equally (NULL = equal weight across pathways)
  NULL
}

