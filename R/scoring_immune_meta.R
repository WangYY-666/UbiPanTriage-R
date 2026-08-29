# ============================================================================
# scoring_immune_meta.R - immune & metabolic dimension scoring
# 閸忓秶鏌呮稉搴濆敩鐠嬨垻娣惔锕佺槑閸掑棴绱欓崺鍝勬礈鐞涖劏鎻稉搴″帳閻ゎ偅韫堝☉?娴狅綀闃块柅鏂跨窞閻ㄥ嫬濮為弶鍐祲閸忚櫕鈧嶇礆
# ============================================================================

#' Tumor samples (patient-deduplicated) shared by expression and trait tables
#' 閸欐牞銆冩潏鍙ョ瑢閻楃懓绶涚悰銊ュ彙閺堝娈戦懖璺ㄦЙ閺嶉攱婀伴敍鍫熷瘻閹綀鈧懎骞撻柌宥忕礉娑撳酣顣╃拋锛勭暬娑撯偓閼疯揪绱?#' @noRd
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

#' Gene-trait spearman correlation matrix (閸╁搫娲?閻楃懓绶?Spearman 閻╃鍙ч惌鈺呮█)
#' Uses pre-computed values for background genes and computes on the fly
#' for user-supplied genes. 閼冲本娅欓崺鍝勬礈閻劑顣╃拋锛勭暬閸婄》绱濋悽銊﹀煕閸╁搫娲滈崡铏鐠侊紕鐣婚妴?#' @noRd
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

#' Immune infiltration correlation (閸忓秶鏌呭ù鍛婇紟閻╃鍙ч幀褏鐓╅梼?
#' @param cache cancer cache (閻у瞼顫掔紓鎾崇摠)
#' @param genes gene symbols (閸╁搫娲滅粭锕€褰?
#' @param method immune method: timer / cibersort / mcp / ssgsea / xcell
#' @return matrix genes x immune traits (閸╁搫娲?x 閸忓秶鏌呯紒鍡氬劒閻╃鍙ч惌鈺呮█)
#' @export
immune_cor <- function(cache, genes, method = "timer") {
  trait_df <- cache$immune[[method]]
  pre_mat  <- cache$pre$cor_immune[[method]]
  .gene_trait_cor(cache, genes, trait_df, pre_mat)
}

#' Metabolic pathway correlation (娴狅綀闃块柅鏂跨窞閻╃鍙ч幀褏鐓╅梼?
#' @param cache cancer cache (閻у瞼顫掔紓鎾崇摠)
#' @param genes gene symbols (閸╁搫娲滅粭锕€褰?
#' @return matrix genes x pathways (閸╁搫娲?x 娴狅綀闃块柅鏂跨窞閻╃鍙ч惌鈺呮█)
#' @export
metabolic_cor <- function(cache, genes) {
  .gene_trait_cor(cache, genes, cache$metabolism, cache$pre$cor_meta)
}

#' Immune dimension score (閸忓秶鏌呯紒鏉戝鐠囧嫬鍨?
#'
#' Weighted mean of |spearman r| between the gene and immune infiltration
#' traits; default weights follow the LUAD example (CD8/Macrophage/DC
#' emphasized). Top trait and signed correlation are also reported.
#' 閸╁搫娲滄稉搴″帳閻ゎ偅韫堝☉锔惧瀵颁胶娈戦崝鐘虫綀 |閻╃鍙х化缁樻殶| 閸у洤鈧》绱辨妯款吇閺夊啴鍣搁崥?LUAD 缁€杞扮伐閵?#'
#' @param cache cancer cache (閻у瞼顫掔紓鎾崇摠)
#' @param genes gene symbols (閸╁搫娲滅粭锕€褰?
#' @param method immune method (閸忓秶鏌呭ù鍛婇紟閺傝纭?
#' @param weights named weights per trait (閸氬嫬鍘ら悿顐ょ矎閼崇偞娼堥柌宥忕幢NULL 娴ｈ法鏁ゆ妯款吇)
#' @return data.frame with Immune_Score, Top_immune_cell, Cor_value1 and
#'   per-trait correlations (閸忓秶鏌呯拠鍕瀻閵嗕焦娓跺铏规祲閸忓磭绮忛懗鐐偓浣烘祲閸忓磭閮撮弫鏉垮挤閸氬嫮绮忛懗鐐垫祲閸?
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

#' Metabolic dimension score (娴狅綀闃跨紒鏉戝鐠囧嫬鍨?
#'
#' Weighted mean of |spearman r| between the gene and 7 GSVA metabolic
#' pathways. 閸╁搫娲滄稉?7 閺夆€插敩鐠嬨垽鈧柨绶為惃鍕閺?|閻╃鍙х化缁樻殶| 閸у洤鈧鈧?#'
#' @param cache cancer cache (閻у瞼顫掔紓鎾崇摠)
#' @param genes gene symbols (閸╁搫娲滅粭锕€褰?
#' @param weights named weights per pathway (閸氬嫪鍞拫銏も偓鏂跨窞閺夊啴鍣搁敍姹礥LL 娴ｈ法鏁ゆ妯款吇)
#' @return data.frame with Metabolic_Score, Top_metabolic_trait, Cor_value2
#'   and per-pathway correlations (娴狅綀闃跨拠鍕瀻閵嗕焦娓跺铏规祲閸忔娊鈧柨绶為崣濠傛倗閼奉亞娴夐崗?
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