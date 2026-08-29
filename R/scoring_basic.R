# ============================================================================
# scoring_basic.R - basic dimension scoring (basic dimension: expression/diff/survival/ROC)
# ============================================================================

#' Expression score (expression score)
#'
#' Bell-shaped preference for moderately expressed genes, based on mean FPKM
#' in tumor tissues: no score below \code{lo} or above \code{hi}; full score
#' between \code{peak_lo} and \code{peak_hi}.
#' @param fpkm_mean mean FPKM per gene
#' @param lo lower detection threshold
#' @param peak_lo lower bound of the full-score plateau
#' @param peak_hi upper bound of the full-score plateau
#' @param hi upper detection threshold
#' @return numeric vector in [0,1]
#' @export
calc_expression_score <- function(fpkm_mean, lo = 5, peak_lo = 10,
                                  peak_hi = 100, hi = 1000) {
  x <- as.numeric(fpkm_mean)
  s <- rep(0, length(x))
  ok <- !is.na(x)
  i1 <- ok & x >= lo & x < peak_lo
  i2 <- ok & x >= peak_lo & x <= peak_hi
  i3 <- ok & x > peak_hi & x <= hi
  s[i1] <- linear_map(x[i1], lo, peak_lo)
  s[i2] <- 1
  s[i3] <- 1 - linear_map(x[i3], peak_hi, hi)
  s[is.na(s)] <- 0
  s
}

#' Differential expression score (differential expression score)
#'
#' For each chosen DEG method, the score is the rank-based percentile of
#' |log2FC| among all genes measured by that method (0-1). Methods are
#' combined with a weighted mean; \code{padj_cutoff} optionally zeroes
#' non-significant genes.
#' @param cache cancer cache
#' @param genes gene symbols
#' @param methods DEG methods to use: "limma", "edgeR", "deseq2", "pair"
#' @param method_weights named numeric weights for each method
#' @param padj_cutoff if set, genes with adjusted p >= cutoff score 0
#' @return data.frame(gene, Diff_Score)
#' @export
calc_diff_score <- function(cache, genes, methods = c("limma", "edgeR", "deseq2", "pair"),
                            method_weights = c(limma = 1, edgeR = 1, deseq2 = 1, pair = 0.5),
                            padj_cutoff = NULL) {
  genes <- unique(as.character(genes))
  methods <- intersect(methods, names(cache$pre$diff_pct))
  if (length(methods) == 0) {
    return(data.frame(gene = genes, Diff_Score = 0, stringsAsFactors = FALSE))
  }
  mw <- method_weights[methods]
  mw[is.na(mw)] <- 1
  scores <- matrix(NA_real_, nrow = length(genes), ncol = length(methods),
                   dimnames = list(genes, methods))
  for (j in seq_along(methods)) {
    m <- methods[j]
    d <- cache$pre$diff_pct[[m]]
    idx <- match(genes, d$gene)
    s <- d$pct[idx]
    if (!is.null(padj_cutoff)) {
      padj <- d$padj[idx]
      s[is.na(padj) | padj >= padj_cutoff] <- 0
    }
    scores[, j] <- s
  }
  wsum <- rowSums(scores * rep(mw, each = nrow(scores)), na.rm = TRUE)
  wtot <- rowSums(!is.na(scores) * rep(mw, each = nrow(scores)))
  out <- ifelse(wtot > 0, wsum / wtot, 0)
  data.frame(gene = genes, Diff_Score = out, stringsAsFactors = FALSE)
}

#' Survival score from univariate Cox regression (survival score from Cox)
#'
#' Significant genes (default p < 0.05) are ranked by the magnitude of the
#' log HR within each direction (protective HR<1 / risk HR>1) and mapped to
#' 0.5-1; non-significant genes score 0.
#' @param cache cancer cache
#' @param genes gene symbols
#' @param p_cutoff significance threshold
#' @return data.frame(gene, Survive_Score, HR, HR_p)
#' @export
calc_survive_score <- function(cache, genes, p_cutoff = 0.05) {
  genes <- unique(as.character(genes))
  cox <- cache$pre$cox
  idx <- match(genes, cox$gene)
  p  <- cox$p[idx]
  hr <- exp(cox$HR_log[idx])
  sig <- !is.na(p) & p < p_cutoff & !is.na(hr) & hr > 0
  out <- rep(0, length(genes))
  dir <- rep(NA_character_, length(genes))
  dir[sig & hr < 1]  <- "prot"
  dir[sig & hr > 1]  <- "risk"
  mag <- rep(NA_real_, length(genes))
  mag[sig] <- abs(cox$HR_log[idx][sig])
  for (d in c("prot", "risk")) {
    i <- which(dir == d)
    if (length(i) > 0) {
      r <- rank(-mag[i], ties.method = "min")
      out[i] <- rank_to_score(r)
    }
  }
  data.frame(gene = genes, Survive_Score = out,
             HR = hr, HR_p = p, direction = dir, stringsAsFactors = FALSE)
}

#' ROC score (survival time-dependent AUC + tumor/normal diagnostic AUC)
#'
#' Survival AUCs are time-dependent ROC AUCs (timeROC) at the chosen years;
#' the diagnostic AUC comes from tumor vs normal ROC. Both are averaged when
#' available; otherwise the available one is used.
#' @param cache cancer cache
#' @param genes gene symbols
#' @param timepoints follow-up years for survival AUC
#' @return data.frame(gene, ROC_Score, AUC_surv, AUC_diag)
#' @export
calc_roc_score <- function(cache, genes, timepoints = c(1, 3, 5)) {
  genes <- unique(as.character(genes))
  sa <- cache$pre$surv_auc
  da <- cache$pre$diag_auc
  idx <- match(genes, sa$gene)
  cols <- paste0("auc", timepoints)
  cols <- intersect(cols, colnames(sa))
  if (length(cols) == 0) stop("no survival AUC columns in cache")
  m <- as.matrix(sa[idx, cols, drop = FALSE])
  surv_auc <- rowMeans(m, na.rm = TRUE)
  surv_auc[is.nan(surv_auc)] <- NA_real_
  diag_auc <- da$auc[match(genes, da$gene)]
  both <- !is.na(surv_auc) & !is.na(diag_auc)
  score <- rep(0, length(genes))
  score[both] <- (surv_auc[both] + diag_auc[both]) / 2
  score[!both & !is.na(surv_auc)] <- surv_auc[!both & !is.na(surv_auc)]
  score[!both & !is.na(diag_auc)] <- diag_auc[!both & !is.na(diag_auc)]
  data.frame(gene = genes, ROC_Score = score,
             AUC_surv = surv_auc, AUC_diag = diag_auc, stringsAsFactors = FALSE)
}

#' Basic dimension score (basic dimension score)
#'
#' Weighted mean of Expression / Diff / Survival / ROC sub-scores.
#' Sub-items that are structurally unavailable for a cancer (e.g. no normal
#' tissue for diff/diagnostic ROC, no survival follow-up) are dropped and
#' the remaining weights are re-normalized, so cancers with incomplete data
#' are NOT systematically under-scored. Unavailable items are reported in
#' the \code{Missing} column (e.g. "diff,diag").
#' @param cache cancer cache
#' @param genes gene symbols
#' @param weights named weights: c(expression, diff, survive, roc)
#' @param expression list(lo, peak_lo, peak_hi, hi) thresholds
#' @param diff list(methods, method_weights, padj_cutoff)
#' @param survive list(p_cutoff)
#' @param roc list(timepoints)
#' @return data.frame with sub-scores, Basic_Score and Missing detail
#' @export
calc_basic_score <- function(cache, genes,
                             weights = c(expression = 1, diff = 1, survive = 1, roc = 1),
                             expression = list(lo = 5, peak_lo = 10, peak_hi = 100, hi = 1000),
                             diff = list(methods = c("limma", "edgeR", "deseq2", "pair"),
                                         method_weights = c(limma = 1, edgeR = 1, deseq2 = 1, pair = 0.5),
                                         padj_cutoff = NULL),
                             survive = list(p_cutoff = 0.05),
                             roc = list(timepoints = c(1, 3, 5))) {
  genes <- unique(as.character(genes))
  fpkm <- cache$pre$fpkm_tumor_mean[genes]
  es  <- calc_expression_score(fpkm, expression$lo %||% 5, expression$peak_lo %||% 10,
                               expression$peak_hi %||% 100, expression$hi %||% 1000)
  ds  <- calc_diff_score(cache, genes, diff$methods %||% c("limma", "edgeR", "deseq2", "pair"),
                         diff$method_weights %||% c(limma = 1, edgeR = 1, deseq2 = 1, pair = 0.5),
                         diff$padj_cutoff)
  ss  <- calc_survive_score(cache, genes, survive$p_cutoff %||% 0.05)
  rs  <- calc_roc_score(cache, genes, roc$timepoints %||% c(1, 3, 5))
  sub <- data.frame(gene = genes,
                    Expression_Score = es,
                    Diff_Score = ds$Diff_Score[match(genes, ds$gene)],
                    Survive_Score = ss$Survive_Score[match(genes, ss$gene)],
                    ROC_Score = rs$ROC_Score[match(genes, rs$gene)],
                    stringsAsFactors = FALSE)
  # ---- structural availability per sub-item ----
  has_os   <- !is.null(cache$clinical) && nrow(cache$clinical) > 0 &&
              any(!is.na(cache$clinical$OS.time) & cache$clinical$OS.time > 0)
  has_norm <- sum(cache$sample_info$type == "Normal") > 1
  avail_item <- c(expression = TRUE,
                  diff      = isTRUE(cache$diff_available),
                  survive   = has_os,
                  roc       = has_os || has_norm)
  m <- as.matrix(sub[, c("Expression_Score", "Diff_Score", "Survive_Score", "ROC_Score")])
  item_cols <- c(expression = "Expression_Score", diff = "Diff_Score",
                  survive = "Survive_Score", roc = "ROC_Score")
  for (nm in names(avail_item)) if (!avail_item[nm]) m[, item_cols[[nm]]] <- NA_real_
  # ---- per-gene weighted mean with re-normalization over available items ----
  w <- weights[c("expression", "diff", "survive", "roc")]
  w[is.na(w)] <- 0
  w[!avail_item] <- 0
  wm <- matrix(rep(w, each = nrow(m)), nrow = nrow(m))
  wm[is.na(m)] <- 0
  wt <- rowSums(wm)
  score <- rowSums(m * wm, na.rm = TRUE) / ifelse(wt > 0, wt, 1)
  score[wt == 0] <- 0
  sub$Basic_Score <- score
  miss_lbl <- c(expression = "expr", diff = "diff", survive = "surv", roc = "diag")
  miss <- apply(m, 1, function(r) paste(unname(miss_lbl[names(avail_item)[is.na(r)]]), collapse = ","))
  sub$Missing <- ifelse(miss == "", NA_character_, miss)
  # mark genes without a diagnostic AUC (no normal tissue etc.)
  i_diag <- is.na(rs$AUC_diag[match(genes, rs$gene)])
  sub$Missing[i_diag] <- ifelse(is.na(sub$Missing[i_diag]), "diag",
                                paste(sub$Missing[i_diag], "diag", sep = ","))
  # ---- detail columns ----
  sub$HR <- ss$HR[match(genes, ss$gene)]
  sub$HR_p <- ss$HR_p[match(genes, ss$gene)]
  sub$direction <- ss$direction[match(genes, ss$gene)]
  sub$AUC_surv <- rs$AUC_surv[match(genes, rs$gene)]
  sub$AUC_diag <- rs$AUC_diag[match(genes, rs$gene)]
  sub
}