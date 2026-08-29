options(warn=-1)
.libPaths(c("D:/PackagesR/R_LIBS", .libPaths()))
suppressMessages({library(survival); library(pROC); library(timeROC); library(parallel)})
`%||%` <- function(x, y) if (is.null(x)) y else x
# ============================================================================
# build_web_lookup.R
# Build the compact web lookup library for pages 2-4 (arbitrary genes):
#   per-gene RDS files: scores + cor vectors + expression slices
#   plus genes_index / sample_info / clinical shared files
# Usage (from project root):
#   Rscript UbiPanTriage/data-raw/build_web_lookup.R                 (all cancers)
#   Rscript UbiPanTriage/data-raw/build_web_lookup.R --only=LUAD,BRCA
#   Rscript UbiPanTriage/data-raw/build_web_lookup.R --max-genes=2000  (test mode)
#   Rscript UbiPanTriage/data-raw/build_web_lookup.R --rebuild
# ============================================================================

# ---- locate project root ----
find_root <- function() {
  sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  cands <- c(getwd(), if (length(sf)) dirname(normalizePath(sf)) else NULL)
  cands <- c(cands, file.path(cands, ".."), file.path(cands, "../.."))
  cands <- unique(normalizePath(cands[nzchar(cands)]))
  hit <- cands[vapply(cands, function(p) dir.exists(file.path(p, "pancaner_data")), logical(1))]
  if (!length(hit)) stop("cannot locate pancaner_data/ ; run Rscript from project root")
  hit[1]
}
args <- commandArgs(TRUE)
get_arg <- function(pat, def) {
  v <- sub(pat, "", args[grepl(pat, args)])
  if (length(v) == 0) def else v
}
only <- unlist(strsplit(get_arg("^--only=", ""), ","))
only <- only[nzchar(only)]
max_genes <- suppressWarnings(as.integer(get_arg("^--max-genes=", ""))); if (is.na(max_genes)) max_genes <- Inf
rebuild <- "--rebuild" %in% args
n_cores <- as.integer(Sys.getenv("UBI_CORES", get_arg("^--cores=", "8")))
min_frac <- as.numeric(get_arg("^--min-frac=", "0.2"))
min_tpm <- as.numeric(get_arg("^--min-tpm=", "1"))

root <- find_root()
pd <- file.path(root, "pancaner_data")
outdir <- file.path(root, "UbiPanTriage", "inst", "extdata", "lookup")
partdir <- file.path(root, "UbiPanTriage", "inst", "extdata", "lookup", ".partial")
dir.create(file.path(outdir, "gene"), recursive = TRUE, showWarnings = FALSE)
dir.create(partdir, recursive = TRUE, showWarnings = FALSE)
cat("root:", root, "\nlookup outdir:", outdir, "\n")

# ---- helpers (same logic as build_pancancer_cache.R) ----
encode_expr <- function(df) {                       # integer-encoded log2(x+1)*100
  m <- as.matrix(df); storage.mode(m) <- "double"
  mx <- apply(m, 1, max, na.rm = TRUE)
  m <- m[is.finite(mx) & mx >= 0.5, , drop = FALSE]
  v <- round(log2(m + 1) * 100); storage.mode(v) <- "integer"; v
}
load_matrix_processed <- function(d, ca, unit, obj) {
  f <- file.path(d, paste0("TCGA_", ca, "_matrix_", unit, "_processed.Rdata"))
  if (!file.exists(f)) return(NULL)
  e <- new.env(); load(f, envir = e)
  if (is.null(e[[obj]])) NULL else encode_expr(e[[obj]])
}
std_diff <- function(df, gene_col, logfc_col, p_col, padj_col) {
  gene <- if (gene_col == "rownames") rownames(df) else df[[gene_col]]
  out <- data.frame(gene = as.character(gene), logFC = as.numeric(df[[logfc_col]]),
                    p = as.numeric(df[[p_col]]), padj = as.numeric(df[[padj_col]]),
                    stringsAsFactors = FALSE)
  out[!is.na(out$gene) & is.finite(out$logFC), , drop = FALSE]
}
cor_one <- function(xgenes, ydf, xmat) {
  if (is.null(ydf)) return(NULL)            # cancer without this method's data
  ydf$pid <- substr(ydf$ID, 1, 12)
  ydf <- ydf[!duplicated(ydf$pid), , drop = FALSE]
  common <- intersect(colnames(xmat), ydf$pid)
  if (length(common) < 10) return(NULL)
  x <- xmat[xgenes, common, drop = FALSE]
  y <- as.matrix(ydf[match(common, ydf$pid), setdiff(colnames(ydf), c("ID", "pid")), drop = FALSE])
  x <- x[apply(x, 1, function(v) sd(v) > 0), , drop = FALSE]
  y <- y[, apply(y, 2, function(v) sd(v) > 0), drop = FALSE]
  if (nrow(x) == 0 || ncol(y) == 0) return(NULL)
  cor(t(x), y, method = "spearman")
}
calc_expression_score <- function(x, lo = 5, peak_lo = 10, peak_hi = 100, hi = 1000) {
  s <- rep(0, length(x)); ok <- !is.na(x)
  i1 <- ok & x >= lo & x < peak_lo; i2 <- ok & x >= peak_lo & x <= peak_hi; i3 <- ok & x > peak_hi & x <= hi
  s[i1] <- (x[i1] - lo) / (peak_lo - lo); s[i2] <- 1
  s[i3] <- 1 - (x[i3] - peak_hi) / (hi - peak_hi); s[is.na(s)] <- 0; s
}
row_wmean <- function(mat, w) {
  w <- w[names(w) %in% colnames(mat)]; m <- mat[, names(w), drop = FALSE]
  if (ncol(m) == 0) return(rep(0, nrow(mat)))
  rowSums(abs(m) * rep(w, each = nrow(m)), na.rm = TRUE) / sum(w)
}
top_trait <- function(v) {
  if (all(is.na(v)) || length(v) == 0) return(list(feature = NA_character_, cor = NA_real_))
  i <- which.max(abs(v)); list(feature = names(v)[i], cor = unname(v[i]))
}
rank_to_score <- function(ranks) {
  if (length(ranks) == 0 || all(is.na(ranks))) return(rep(NA_real_, length(ranks)))
  (ranks - 1) / (max(ranks, na.rm = TRUE) - 1) * 0.5 + 0.5
}
# int16 pack (log2(x+1)*100, -32768 = NA)
pack16 <- function(v) { v[is.na(v)] <- -32768; writeBin(as.integer(round(v)), raw(), size = 2) }

CANCER_NAMES <- c(
  ACC="Adrenocortical carcinoma",BLCA="Bladder urothelial carcinoma",BRCA="Breast invasive carcinoma",
  CESC="Cervical squamous cell carcinoma",CHOL="Cholangiocarcinoma",COAD="Colon adenocarcinoma",
  DLBC="Diffuse large B-cell lymphoma",ESCA="Esophageal carcinoma",GBM="Glioblastoma multiforme",
  HNSC="Head and neck squamous cell carcinoma",KICH="Kidney chromophobe",KIRC="Kidney renal clear cell carcinoma",
  KIRP="Kidney renal papillary cell carcinoma",LAML="Acute myeloid leukemia",LGG="Brain lower grade glioma",
  LIHC="Liver hepatocellular carcinoma",LUAD="Lung adenocarcinoma",LUSC="Lung squamous cell carcinoma",
  MESO="Mesothelioma",OV="Ovarian serous cystadenocarcinoma",PAAD="Pancreatic adenocarcinoma",
  PCPG="Pheochromocytoma and paraganglioma",PRAD="Prostate adenocarcinoma",READ="Rectum adenocarcinoma",
  SARC="Sarcoma",SKCM="Skin cutaneous melanoma",STAD="Stomach adenocarcinoma",TGCT="Testicular germ cell tumors",
  THCA="Thyroid carcinoma",THYM="Thymoma",UCEC="Uterine corpus endometrial carcinoma",UCS="Uterine carcinosarcoma",
  UVM="Uveal melanoma")

# ubiquitin annotation from the existing precomputed scores
pcs <- readRDS(file.path(root, "UbiPanTriage", "inst", "extdata", "pancancer_ubi_scores.rds"))
ubi_ann <- unique(pcs$scores[, c("Gene", "Ubi_Score", "Ubi_Type", "Ubi_Type_Full"), drop = FALSE])

# per-gene workers (globals staged per cancer, same as build_pancancer_cache.R)
one_cox <- function(g) {
  f <- as.formula(paste("Surv(OS.time, OS) ~", safe[g]))
  fit <- tryCatch(coxph(f, data = cox_df), error = function(e) NULL)
  if (is.null(fit)) return(c(NA_real_, NA_real_))
  s <- summary(fit); c(s$coefficients[1, "coef"], s$coefficients[1, "Pr(>|z|)"])
}
one_roc <- function(g) {
  vapply(c(1, 3, 5), function(tp) {
    r <- tryCatch(timeROC(T = surv_df$OS.time, delta = surv_df$OS, marker = surv_df[[safe[g]]],
                          cause = 1, weighting = "aalen", time = tp, ROC = FALSE), error = function(e) NULL)
    if (is.null(r)) NA_real_ else as.numeric(r$AUC[2])
  }, numeric(1))
}
one_diag <- function(g) {
  vals <- c(as.numeric(norm_log2[g, ]), as.numeric(expr_tumor_log2[g, ]))
  rr <- tryCatch(roc(diag_type ~ vals, levels = c("Normal", "Tumor"), direction = ">", quiet = TRUE),
                 error = function(e) NULL)
  if (is.null(rr)) NA_real_ else as.numeric(rr$auc)
}

build_one <- function(ca) {
  d <- file.path(pd, paste0("TCGA-", ca))
  pf <- file.path(partdir, paste0(ca, ".rds"))
  if (!rebuild && file.exists(pf)) { cat("skip", ca, "(partial exists)\n"); return(TRUE) }
  cat("\n[", ca, "] loading ...\n", sep = "")
  e <- new.env(); load(file.path(d, paste0("TCGA_", ca, "_matrix_tpm_processed.Rdata")), envir = e)
  tpm_raw <- as.matrix(e$expr_tpm)                       # full doubles
  st <- substr(colnames(tpm_raw), 14, 16)
  sample_info <- data.frame(sample = colnames(tpm_raw),
    type = ifelse(st %in% c("01A","01B","01C","02A","02B","03A","03B","05A","06A","06B","07A","07B","07C","01R"), "Tumor",
           ifelse(st %in% c("11A","11B","11C"), "Normal", "Other")), stringsAsFactors = FALSE)
  mx <- apply(tpm_raw, 1, max, na.rm = TRUE)
  tpm_raw <- tpm_raw[is.finite(mx) & mx >= 0.5, , drop = FALSE]
  cat("  tpm genes:", nrow(tpm_raw), " samples:", ncol(tpm_raw), "\n")

  # expressed universe: TPM >= min_tpm in >= min_frac of tumor samples
  tum <- sample_info$sample[sample_info$type == "Tumor"]
  tum <- tum[!duplicated(substr(tum, 1, 12))]
  frac <- apply(tpm_raw[, tum, drop = FALSE], 1, function(v) mean(v >= min_tpm, na.rm = TRUE))
  expr_genes <- sort(rownames(tpm_raw)[frac >= min_frac])
  if (length(expr_genes) == 0) stop("no genes pass the expression filter")
  if (is.finite(max_genes)) expr_genes <- head(expr_genes, max_genes)
  cat("  expressed genes (TPM>=", min_tpm, " in >=", min_frac, "):", length(expr_genes), "\n", sep = "")

  # int-encoded full matrices (aligned to tpm sample order)
  align_cols <- function(m) {
    if (is.null(m)) return(NULL)
    m <- m[, intersect(sample_info$sample, colnames(m)), drop = FALSE]
    m[, sample_info$sample[sample_info$sample %in% colnames(m)], drop = FALSE]
  }
  tpm_int <- encode_expr(tpm_raw); tpm_int <- align_cols(tpm_int)
  cnt_int <- align_cols(load_matrix_processed(d, ca, "count", "expr_count"))
  fpkm_int <- align_cols(load_matrix_processed(d, ca, "fpkm", "expr_fpkm"))

  # FPKM tumor mean (expression score source)
  e <- new.env(); load(file.path(d, paste0("TCGA_", ca, "_matrix_fpkm.Rdata")), envir = e)
  mat_fpkm <- e$mat_fpkm
  ts <- intersect(sample_info$sample[sample_info$type == "Tumor"], colnames(mat_fpkm))
  fpkm_tumor_mean <- setNames(numeric(0), character(0))
  if (length(ts) > 0) {
    ft <- mat_fpkm[, c("symbol", ts), drop = FALSE]
    ft[, -1] <- lapply(ft[, -1, drop = FALSE], function(x) as.numeric(x))
    ft <- ft[!is.na(ft$symbol) & ft$symbol != "", , drop = FALSE]
    fm <- rowMeans(as.matrix(ft[, -1]), na.rm = TRUE); names(fm) <- ft$symbol
    fpkm_tumor_mean <- tapply(fm, names(fm), mean)
  }

  # differential results
  diff <- list(); diffFile <- file.path(d, paste0(ca, "_diff.Rdata"))
  if (file.exists(diffFile)) {
    e <- new.env(); load(diffFile, envir = e)
    add <- function(tag, obj, gc, lfc, pc, pj) if (!is.null(e[[obj]]))
      diff[[tag]] <<- std_diff(as.data.frame(e[[obj]]), gc, lfc, pc, pj)
    add("limma",  "res_limma",     "rownames", "logFC",          "P.Value",    "adj.P.Val")
    add("edgeR",  "res_edgeR",     "rownames", "logFC",          "PValue",     "FDR")
    add("deseq2", "res_deseq2_df", "gene_id",  "log2FoldChange", "pvalue",     "padj")
    add("pair",   "res_pair_df",   "gene_id",  "log2FoldChange", "pvalue",     "padj")
  }
  diff <- diff[lengths(diff) > 0]
  diff_available <- length(diff) > 0

  # clinical
  e <- new.env(); load(file.path(d, paste0("clin_", ca, "_pro.Rdata")), envir = e)
  clin <- as.data.frame(e$clin_pro)
  clin <- clin[, intersect(c("sample", "OS", "OS.time"), colnames(clin)), drop = FALSE]
  clin <- clin[!is.na(clin$sample) & clin$OS %in% c(0, 1), , drop = FALSE]
  clin <- clin[!duplicated(clin$sample), , drop = FALSE]

  # immune + metabolism
  e <- new.env(); load(file.path(d, paste0(ca, "_immune_cells.Rdata")), envir = e)
  immune <- list(timer = if (is.null(e$timer_res)) NULL else as.data.frame(e$timer_res),
                 cibersort = if (is.null(e$cibersort_res)) NULL else as.data.frame(e$cibersort_res),
                 mcp = if (is.null(e$mcp_res)) NULL else as.data.frame(e$mcp_res),
                 ssgsea = if (is.null(e$ssgsea_res)) NULL else as.data.frame(e$ssgsea_res),
                 xcell = if (is.null(e$xcell_res)) NULL else as.data.frame(e$xcell_res))
  method_sfx <- function(m) switch(m, timer = "_TIMER", cibersort = "_CIBERSORT",
                                   mcp = "_MCPcounter", xcell = "_xCell", ssgsea = "")
  tidy_immune <- function(df, m) {
    if (is.null(df)) return(NULL)
    df <- df[!is.na(df$ID) & df$ID != "", , drop = FALSE]
    keep <- vapply(df, is.numeric, logical(1)); keep["ID"] <- TRUE
    df <- df[, keep, drop = FALSE]; df <- df[!duplicated(df$ID), , drop = FALSE]
    sfx <- method_sfx(m)
    if (sfx != "") colnames(df) <- sub(paste0(sfx, "$"), "", colnames(df))
    drop <- grepl("^(P.value|Correlation|RMSE)$", colnames(df))
    if (any(drop)) df <- df[, !drop, drop = FALSE]
    if (m == "xcell") {
      xdrop <- grepl("^(ImmuneScore|StromaScore|MicroenvironmentScore)$", colnames(df))
      if (any(xdrop)) df <- df[, !xdrop, drop = FALSE]
    }
    df
  }
  immune <- lapply(names(immune), function(m) tidy_immune(immune[[m]], m))
  names(immune) <- c("timer", "cibersort", "mcp", "ssgsea", "xcell")
  e <- new.env(); load(file.path(d, paste0(ca, "_metabolism7_GSVA.Rdata")), envir = e)
  metabolism <- as.data.frame(e$GSVA_Metabolism)
  metabolism <- metabolism[!is.na(metabolism$ID) & metabolism$ID != "", , drop = FALSE]
  metabolism <- metabolism[!duplicated(metabolism$ID), , drop = FALSE]

  # tumor log2 TPM (patient IDs)
  tum12 <- substr(tum, 1, 12)
  expr_tumor_log2 <- log2(tpm_raw[expr_genes, tum, drop = FALSE] + 1)
  colnames(expr_tumor_log2) <- tum12
  clin_tumor <- clin[clin$sample %in% tum12, , drop = FALSE]
  clin_tumor <- clin_tumor[!is.na(clin_tumor$OS.time) & is.finite(clin_tumor$OS.time), , drop = FALSE]
  gene_df <- as.data.frame(t(expr_tumor_log2), check.names = FALSE)
  gene_df$sample <- rownames(gene_df)
  cox_df <- merge(clin_tumor, gene_df, by = "sample")
  cox_df <- cox_df[!is.na(cox_df$OS) & !is.na(cox_df$OS.time) & cox_df$OS.time > 0, , drop = FALSE]
  safe <- make.names(expr_genes); names(safe) <- expr_genes
  colnames(cox_df)[4:ncol(cox_df)] <- safe[colnames(cox_df)[4:ncol(cox_df)]]
  surv_df <- cox_df; surv_df$OS.time <- surv_df$OS.time / 365
  W <- list(safe = safe, cox_df = cox_df, surv_df = surv_df, expr_tumor_log2 = expr_tumor_log2,
            norm_log2 = NULL, diag_type = NULL)
  list2env(W, envir = globalenv())
  if (!is.null(cl)) clusterExport(cl, names(W), envir = globalenv())

  cx <- if (is.null(cl)) vapply(expr_genes, one_cox, numeric(2)) else parSapply(cl, expr_genes, one_cox)
  cat("  cox done, sig:", sum(cx[2, ] < 0.05, na.rm = TRUE), "\n")
  rc <- if (is.null(cl)) sapply(expr_genes, one_roc) else parSapply(cl, expr_genes, one_roc)
  cat("  timeROC done\n")
  nor <- sample_info$sample[sample_info$type == "Normal"]
  diag_auc <- setNames(rep(NA_real_, length(expr_genes)), expr_genes)
  if (length(nor) > 1) {
    norm_log2 <- log2(tpm_raw[expr_genes, nor, drop = FALSE] + 1)
    diag_type <- c(rep("Normal", ncol(norm_log2)), rep("Tumor", ncol(expr_tumor_log2)))
    assign("norm_log2", norm_log2, envir = globalenv()); assign("diag_type", diag_type, envir = globalenv())
    if (!is.null(cl)) clusterExport(cl, c("norm_log2", "diag_type"), envir = globalenv())
    dd <- if (is.null(cl)) vapply(expr_genes, one_diag, numeric(1)) else parSapply(cl, expr_genes, one_diag)
    diag_auc[expr_genes] <- dd
  }
  cat("  diag ROC done, non-NA:", sum(!is.na(diag_auc)), "\n")

  cor_immune <- lapply(immune, function(dd) cor_one(expr_genes, dd, expr_tumor_log2))
  cor_meta   <- cor_one(expr_genes, metabolism, expr_tumor_log2)
  diff_pct <- lapply(diff, function(dd) {
    dd <- dd[dd$gene %in% expr_genes, , drop = FALSE]
    if (nrow(dd) == 0) return(NULL)
    r <- rank(-abs(dd$logFC), ties.method = "min"); dd$pct <- r / nrow(dd)
    dd[, c("gene", "logFC", "p", "padj", "pct"), drop = FALSE]
  })
  fpkm_tumor_mean <- fpkm_tumor_mean[names(fpkm_tumor_mean) %in% expr_genes]
  cat("  stats computed\n")

  # ---- assemble per-gene score table ----
  has_os   <- nrow(clin) > 0 && any(!is.na(clin$OS.time) & clin$OS.time > 0)
  has_norm <- sum(sample_info$type == "Normal") > 1
  avail_item <- c(expression = TRUE, diff = diff_available, survive = has_os, roc = (has_os || has_norm))
  fpk <- fpkm_tumor_mean[expr_genes]
  es <- calc_expression_score(fpk)
  diff_eff <- lapply(diff_pct, function(dd) dd$pct[match(expr_genes, dd$gene)])
  diff_w <- c(limma = 1, edgeR = 1, deseq2 = 1, pair = 0.5)
  dm <- do.call(cbind, diff_eff)
  if (!is.null(dm) && nrow(dm) > 0) dm[is.na(dm)] <- 0
  dw <- diff_w[names(diff_eff)]; dw[is.na(dw)] <- 1
  ds <- if (!is.null(dm) && ncol(dm) > 0)
          rowSums(dm * rep(dw, each = nrow(dm))) / sum(dw)
        else rep(NA_real_, length(expr_genes))
  pv <- do.call(cbind, lapply(diff_pct, function(dd) dd$p[match(expr_genes, dd$gene)]))
  hr <- exp(cx[1, ]); hp <- cx[2, ]
  sig <- !is.na(hp) & hp < 0.05 & !is.na(hr) & hr > 0
  dirn <- rep(NA_character_, length(expr_genes)); dirn[sig & hr < 1] <- "prot"; dirn[sig & hr > 1] <- "risk"
  ss <- rep(0, length(expr_genes))
  for (di in c("prot", "risk")) {
    i <- which(dirn == di); if (length(i) > 0) { r <- rank(-abs(cx[1, i]), ties.method = "min"); ss[i] <- rank_to_score(r) }
  }
  sa <- do.call(rbind, lapply(1:3, function(k) rc[k, ])); surv_auc <- colMeans(sa, na.rm = TRUE)
  surv_auc[is.nan(surv_auc)] <- NA_real_
  both <- !is.na(surv_auc) & !is.na(diag_auc)
  rocs <- rep(0, length(expr_genes)); rocs[both] <- (surv_auc[both] + diag_auc[both]) / 2
  rocs[!both & !is.na(surv_auc)] <- surv_auc[!both & !is.na(surv_auc)]
  rocs[!both & !is.na(diag_auc)] <- diag_auc[!both & !is.na(diag_auc)]
  m4 <- cbind(es, ds, ss, rocs); colnames(m4) <- c("Expression_Score", "Diff_Score", "Survive_Score", "ROC_Score")
  m4cols <- c(expression = "Expression_Score", diff = "Diff_Score",
              survive = "Survive_Score", roc = "ROC_Score")
  for (nm in names(avail_item)) if (!avail_item[nm]) m4[, m4cols[nm]] <- NA_real_
  w4 <- c(expression = 1, diff = 1, survive = 1, roc = 1); w4[!avail_item] <- 0
  wm <- matrix(rep(w4, each = nrow(m4)), nrow = nrow(m4)); wm[is.na(m4)] <- 0
  wt <- rowSums(wm); basic <- rowSums(m4 * wm, na.rm = TRUE) / ifelse(wt > 0, wt, 1); basic[wt == 0] <- 0
  miss_lbl <- c(expression = "expr", diff = "diff", survive = "surv", roc = "diag")
  miss <- apply(m4, 1, function(r) paste(unname(miss_lbl[names(avail_item)[is.na(r)]]), collapse = ","))
  miss <- ifelse(miss == "", NA_character_, miss)
  miss[is.na(diag_auc)] <- ifelse(is.na(miss[is.na(diag_auc)]), "diag",
                                  paste(miss[is.na(diag_auc)], "diag", sep = ","))

  # immune / metabolic columns (all methods, all traits)
  imm_df <- data.frame(Gene = expr_genes, stringsAsFactors = FALSE)
  imm_default <- rep(NA_real_, length(expr_genes)); meta_df <- NULL
  imm_top <- rep(NA_character_, length(expr_genes)); imm_topc <- rep(NA_real_, length(expr_genes))
  for (m in names(cor_immune)) {
    cm <- cor_immune[[m]]
    if (is.null(cm)) next
    tr <- colnames(cm)
    cols <- cm[expr_genes, , drop = FALSE]
    colnames(cols) <- paste0("Imm_", m, "_", tr)
    imm_df <- cbind(imm_df, cols)
    if (m == "timer" || is.na(imm_default[1])) {
      sc <- row_wmean(cm[expr_genes, , drop = FALSE], stats::setNames(rep(1, length(tr)), tr))
      imm_default <- sc
      tt <- apply(cols, 1, top_trait)
      imm_top <- vapply(tt, `[[`, character(1), "feature"); imm_topc <- vapply(tt, `[[`, numeric(1), "cor")
    }
  }
  if (!is.null(cor_meta)) {
    cm <- cor_meta[expr_genes, , drop = FALSE]
    meta_df <- cm; colnames(meta_df) <- paste0("Meta_", colnames(meta_df))
    meta_score <- row_wmean(cm, stats::setNames(rep(1, ncol(cm)), colnames(cm)))
    mt <- apply(meta_df, 1, top_trait)
    meta_top <- vapply(mt, `[[`, character(1), "feature"); meta_topc <- vapply(mt, `[[`, numeric(1), "cor")
  } else {
    meta_score <- rep(NA_real_, length(expr_genes)); meta_top <- rep(NA_character_, length(expr_genes))
    meta_topc <- rep(NA_real_, length(expr_genes))
  }

  scores <- data.frame(Gene = expr_genes,
                       Ubi_Score = ubi_ann$Ubi_Score[match(expr_genes, ubi_ann$Gene)],
                       Ubi_Type_Full = ubi_ann$Ubi_Type_Full[match(expr_genes, ubi_ann$Gene)],
                       Basic_Score = basic, Missing = miss,
                       Immune_Score = imm_default, Top_immune_cell = imm_top, Cor_value1 = imm_topc,
                       Metabolic_Score = meta_score, Top_metabolic_trait = meta_top, Cor_value2 = meta_topc,
                       Expression_Score = m4[, "Expression_Score"], Diff_Score = m4[, "Diff_Score"],
                       Survive_Score = m4[, "Survive_Score"], ROC_Score = m4[, "ROC_Score"],
                       HR = hr, HR_p = hp, direction = dirn,
                       AUC_surv = surv_auc, AUC_diag = diag_auc,
                       stringsAsFactors = FALSE)
  scores <- cbind(scores, imm_df[, -1, drop = FALSE])
  if (!is.null(meta_df)) scores <- cbind(scores, meta_df)
  rownames(scores) <- NULL
  cat("  scores dim:", dim(scores), "\n")

  # ---- expression slices (int16, aligned to sample order) ----
  slices <- list()
  get_row <- function(m, g) {
    if (is.null(m) || !g %in% rownames(m)) return(NULL)
    pack16(m[g, , drop = TRUE])
  }
  slices$tpm   <- lapply(expr_genes, get_row, m = tpm_int)
  slices$count <- lapply(expr_genes, get_row, m = cnt_int)
  slices$fpkm  <- lapply(expr_genes, get_row, m = fpkm_int)
  names(slices$tpm) <- names(slices$count) <- names(slices$fpkm) <- expr_genes

  partial <- list(cancer = ca, cancer_name = CANCER_NAMES[[ca]] %||% ca,
                  sample_order = sample_info$sample, sample_type = sample_info$type,
                  clinical = clin, genes = expr_genes, scores = scores, slices = slices)
  saveRDS(partial, pf, compress = "gzip")
  cat("  partial saved:", pf, round(file.info(pf)$size / 1e6, 1), "MB\n")
  rm(list = c("expr_tumor_log2", "cox_df", "surv_df", "norm_log2", "gene_df"), envir = globalenv())
  invisible(TRUE)
}

# ---- merge partials into per-gene files (chunked; low memory) ----
merge_all <- function() {
  pfs <- list.files(partdir, pattern = "\\.rds$", full.names = TRUE)
  if (length(pfs) == 0) stop("no partial files; run the per-cancer phase first")
  cas <- sub("\\.rds$", "", basename(pfs))
  cat("\nmerging", length(pfs), "cancers into per-gene files ...\n")
  read_genes <- function(f) readRDS(f)$genes
  all_genes <- sort(unique(unlist(lapply(pfs, read_genes))))
  cat("union genes:", length(all_genes), "\n")

  # genes whose symbol contains filename-invalid chars need a sanitized file name
  file_map <- character()
  bad <- grepl("[/\\:*?\"<>|]", all_genes)
  if (any(bad)) {
    used <- all_genes[!bad]
    for (g in all_genes[bad]) {
      safe <- gsub("[/\\:*?\"<>|]", "_", g)
      while (safe %in% used) safe <- paste0(safe, "_")
      used <- c(used, safe)
      file_map[[g]] <- safe
    }
    cat("sanitized filenames for", length(file_map), "genes\n")
  }

  si_all <- clin_all <- vector("list", length(cas)); names(si_all) <- names(clin_all) <- cas
  for (ca in cas) {
    p <- readRDS(file.path(partdir, paste0(ca, ".rds")))
    si_all[[ca]] <- data.frame(sample = p$sample_order, type = p$sample_type, stringsAsFactors = FALSE)
    clin_all[[ca]] <- p$clinical
    rm(p); gc()
  }
  saveRDS(si_all, file.path(outdir, "sample_info.rds"), compress = "gzip")
  saveRDS(clin_all, file.path(outdir, "clinical.rds"), compress = "gzip")
  cat("sample_info / clinical written\n")

  chunk <- 1500
  nchunk <- ceiling(length(all_genes) / chunk)
  cancer_list <- vector("list", length(all_genes)); names(cancer_list) <- all_genes
  failed <- character()
  for (k in seq_len(nchunk)) {
    gs <- all_genes[((k - 1) * chunk + 1):min(k * chunk, length(all_genes))]
    pk <- vector("list", length(cas)); names(pk) <- cas
    for (ca in cas) {
      p <- readRDS(file.path(partdir, paste0(ca, ".rds")))
      j <- match(gs, p$genes); j <- j[!is.na(j)]
      pk[[ca]] <- list(genes = p$genes[j],
                       scores = p$scores[j, , drop = FALSE],
                       slices = lapply(p$slices, function(s) s[j]))
      rm(p)
    }
    for (g in gs) {
      per <- list()
      for (ca in cas) {
        q <- pk[[ca]]; j <- match(g, q$genes)
        if (is.na(j)) next
        row <- q$scores[j, , drop = FALSE]; row$Cancer <- ca
        per[[ca]] <- list(scores = row,
                          expr = list(tpm = q$slices$tpm[[j]],
                                      count = q$slices$count[[j]],
                                      fpkm = q$slices$fpkm[[j]]))
      }
      if (length(per) == 0) next
      fn <- if (g %in% names(file_map)) file_map[[g]] else g
      ok <- tryCatch({
        saveRDS(list(gene = g, per_cancer = per),
                file.path(outdir, "gene", paste0(fn, ".rds")), compress = "gzip")
        TRUE
      }, error = function(e) {
        cat("SAVEFAIL", g, ":", conditionMessage(e), " nchar:", nchar(g), " chunk:", k, "\n")
        FALSE
      })
      if (ok) cancer_list[[g]] <- names(per) else failed <- c(failed, g)
    }
    cat(sprintf("  chunk %d/%d done (%d genes)\n", k, nchunk, length(gs)))
    rm(pk); gc()
  }
  keep <- setdiff(all_genes, failed)
  if (length(file_map)) saveRDS(file_map, file.path(outdir, "file_map.rds"), compress = "gzip")
  idx <- data.frame(Gene = keep,
                    Cancers = vapply(cancer_list[keep], paste, character(1),
                                     collapse = ","),
                    stringsAsFactors = FALSE)
  saveRDS(idx, file.path(outdir, "genes_index.rds"), compress = "gzip")
  cat("merge done. per-gene files:", length(keep), " failed:", length(failed), "\n")
}
# ---- main ----
cancers <- sort(sub("^TCGA-", "", list.dirs(pd, recursive = FALSE, full.names = FALSE)))
if (length(only)) cancers <- intersect(cancers, only)
n_cores <- min(n_cores, max(1, detectCores()))
cl <- NULL
if (n_cores > 1 && length(cancers) > 1) {
  cl <- makeCluster(min(n_cores, detectCores()), outfile = "")
  clusterEvalQ(cl, { .libPaths(c("D:/PackagesR/R_LIBS", .libPaths())); library(survival); library(timeROC); library(pROC) })
}
cat("cancers:", length(cancers), " cores:", if (is.null(cl)) 1 else length(cl), "\n")
for (ca in cancers) {
  tryCatch(build_one(ca), error = function(e) cat("!! FAILED", ca, ":", conditionMessage(e), "\n"))
}
if (!is.null(cl)) stopCluster(cl)
merge_all()
cat("\nALL DONE\n")





