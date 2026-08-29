# ============================================================================
# build_pancancer_cache.R
# Build per-cancer caches for all TCGA cancers under pancaner_data/.
# Each cache is saved to UbiPanTriage/inst/extdata/<CANCER>_cache.rds
# Run from project root:
#   Rscript UbiPanTriage/data-raw/build_pancancer_cache.R             (skip existing)
#   Rscript UbiPanTriage/data-raw/build_pancancer_cache.R --rebuild   (force all)
#   Rscript UbiPanTriage/data-raw/build_pancancer_cache.R --only=BRCA,OV
# Resumable: finished cancers are skipped unless --rebuild.
# Parallel: set env var UBI_CORES (default 8) to tune worker count.
# ============================================================================

.libPaths(c("D:/PackagesR/R_LIBS", .libPaths()))
suppressMessages({
  library(dplyr); library(tidyr); library(tibble); library(stringr)
  library(survival); library(pROC); library(timeROC); library(readxl)
  library(parallel)
})
`%||%` <- function(x, y) if (is.null(x)) y else x
# integer-encoded log2(x+1)*100 matrix (genes x samples), ~1% expression precision
encode_expr <- function(df) {
  m <- as.matrix(df)
  storage.mode(m) <- "double"
  mx <- apply(m, 1, max, na.rm = TRUE)
  m <- m[is.finite(mx) & mx >= 0.5, , drop = FALSE]
  v <- round(log2(m + 1) * 100)
  storage.mode(v) <- "integer"
  v
}
load_matrix_processed <- function(d, ca, unit, obj) {
  f <- file.path(d, paste0("TCGA_", ca, "_matrix_", unit, "_processed.Rdata"))
  if (!file.exists(f)) return(NULL)
  e <- new.env(); load(f, envir = e)
  if (is.null(e[[obj]])) NULL else encode_expr(e[[obj]])
}

# ---- locate project root (folder containing pancaner_data/) ----
find_root <- function() {
  sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  cands <- c(getwd(), if (length(sf)) dirname(normalizePath(sf)) else NULL)
  cands <- c(cands, file.path(cands, ".."), file.path(cands, "../.."))
  cands <- unique(normalizePath(cands[nzchar(cands)]))
  hit <- cands[vapply(cands, function(p) dir.exists(file.path(p, "pancaner_data")), logical(1))]
  if (!length(hit)) stop("cannot locate pancaner_data/ ; run Rscript from project root")
  hit[1]
}
root  <- find_root()
pd    <- file.path(root, "pancaner_data")
outdir <- file.path(root, "UbiPanTriage", "inst", "extdata")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
cat("root:", root, "\n")# ---- shared helper: standardize a differential result table ----
std_diff <- function(df, gene_col, logfc_col, p_col, padj_col) {
  gene <- if (gene_col == "rownames") rownames(df) else df[[gene_col]]
  out <- data.frame(gene = as.character(gene),
                    logFC = as.numeric(df[[logfc_col]]),
                    p = as.numeric(df[[p_col]]),
                    padj = as.numeric(df[[padj_col]]), stringsAsFactors = FALSE)
  out[!is.na(out$gene) & is.finite(out$logFC), , drop = FALSE]
}

# ---- shared helper: gene x trait Spearman correlation ----
cor_one <- function(xgenes, ydf, xmat) {
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

# ---- per-gene worker functions (globals exported per cancer) ----
one_cox  <- function(g) {
  f <- as.formula(paste("Surv(OS.time, OS) ~", safe[g]))
  fit <- tryCatch(coxph(f, data = cox_df), error = function(e) NULL)
  if (is.null(fit)) return(c(NA_real_, NA_real_))
  s <- summary(fit)
  c(s$coefficients[1, "coef"], s$coefficients[1, "Pr(>|z|)"])
}
one_roc <- function(g) {
  vapply(c(1, 3, 5), function(tp) {
    r <- tryCatch(timeROC(T = surv_df$OS.time, delta = surv_df$OS, marker = surv_df[[safe[g]]],
                          cause = 1, weighting = "aalen", time = tp, ROC = FALSE),
                  error = function(e) NULL)
    if (is.null(r)) NA_real_ else as.numeric(r$AUC[2])
  }, numeric(1))
}
one_diag <- function(g) {
  vals <- c(as.numeric(norm_log2[g, ]), as.numeric(expr_tumor_log2[g, ]))
  rr <- tryCatch(roc(diag_type ~ vals, levels = c("Normal", "Tumor"), direction = ">", quiet = TRUE),
                 error = function(e) NULL)
  if (is.null(rr)) NA_real_ else as.numeric(rr$auc)
}# ---- ubiquitin annotation (shared across cancers; from IUCCD2.0) ----
e <- new.env(); load(file.path(root, "ubiquitin_genes.RData"), envir = e)
ubi_genes <- e$ubiquitin_genes; set2 <- e$ubiquitin_genes2; set3 <- e$ubiquitin_genes3
xlsx <- file.path(root, "IUCCD2.0 ubiquitin genes.xlsx")
raw <- as.data.frame(read_excel(xlsx, sheet = 1, col_names = FALSE), stringsAsFactors = FALSE)
dat <- raw[-(1:2), , drop = FALSE]            # row1 = title, row2 = header
annot_raw <- data.frame(gene = as.character(dat[[3]]), uniprot = as.character(dat[[2]]),
                        ensembl = as.character(dat[[4]]), family = as.character(dat[[5]]),
                        species = as.character(dat[[6]]), stringsAsFactors = FALSE)
annot_raw <- annot_raw[!is.na(annot_raw$gene) & annot_raw$species == "Homo sapiens", , drop = FALSE]
annot_raw <- annot_raw[!duplicated(annot_raw$gene), , drop = FALSE]

prio <- c("E1", "E2", "E3 activity", "E3 adaptor", "DUB", "UBD", "ULD")
extract_primary <- function(fam) {
  if (is.na(fam) || fam == "") return(NA_character_)
  cls <- strsplit(fam, "; ", fixed = TRUE)[[1]]
  for (p in prio) { hit <- cls[startsWith(cls, p)]; if (length(hit) > 0) return(hit[1]) }
  cls[1]
}
annot_raw$primary <- vapply(annot_raw$family, extract_primary, character(1))

fn_scores <- c(
  "E1" = 1.0, "E2/UBC" = 0.9, "E2/UEV" = 0.85, "E2/Other" = 0.8,
  "E3 activity/RING/RING" = 0.95, "E3 activity/HECT" = 0.92, "E3 activity/RBR" = 0.90,
  "E3 activity/RING/U-box" = 0.88, "E3 activity/RING/PHD" = 0.85, "E3 activity/Other" = 0.75,
  "E3 adaptor/Cullin RING/SCF/F-box" = 0.82, "E3 adaptor/Cullin RING/DCX/DWD" = 0.80,
  "E3 adaptor/Cullin RING/ECS/SOCS_VHL_BC-box" = 0.78, "E3 adaptor/Cullin RING/BCR/BTB_3-box" = 0.76,
  "E3 adaptor/Cullin RING/BCR/BTB_Other" = 0.74, "E3 adaptor/Cullin RING/APC_C/APC_C" = 0.80,
  "E3 adaptor/Cullin RING/Other" = 0.70,
  "DUB/USP" = 0.75, "DUB/OTU" = 0.72, "DUB/UCH" = 0.70, "DUB/JAMM" = 0.68,
  "DUB/ULP" = 0.65, "DUB/Josephin" = 0.65, "DUB/Other" = 0.60,
  "UBD/Alpha-Helix/UBA" = 0.55, "UBD/Alpha-Helix/UIM" = 0.52, "UBD/Alpha-Helix/CUE" = 0.50,
  "UBD/ZnF/UBZ" = 0.53, "UBD/ZnF/NZF" = 0.51, "UBD/Other/Other" = 0.45,
  "ULD/UDP/UFD/PB1" = 0.40, "ULD/UDP/UFD/UBX" = 0.42, "ULD/UBL/SUMO" = 0.38, "ULD/UBL/ATG8" = 0.36)
score_family <- function(cls) {
  if (is.na(cls)) return(0.3)
  for (k in names(fn_scores)) if (grepl(k, cls, fixed = TRUE)) return(unname(fn_scores[k]))
  if (grepl("E1", cls)) return(0.95)
  if (grepl("E2", cls)) return(0.85)
  if (grepl("E3 activity", cls)) return(0.80)
  if (grepl("E3 adaptor", cls)) return(0.70)
  if (grepl("DUB", cls)) return(0.65)
  if (grepl("UBD", cls)) return(0.50)
  if (grepl("ULD", cls)) return(0.40)
  0.3
}
annot_raw$fn_score <- vapply(annot_raw$primary, score_family, numeric(1))
annot_raw$in_set2 <- annot_raw$gene %in% set2
annot_raw$in_set3 <- annot_raw$gene %in% set3
annot_raw$set_bonus <- ifelse(annot_raw$in_set2 & annot_raw$in_set3, 0.15,
                       ifelse(annot_raw$in_set2 | annot_raw$in_set3, 0.08, 0))
annot_raw$ubi_score <- pmin(annot_raw$fn_score + annot_raw$set_bonus, 1)

annot <- data.frame(gene = ubi_genes, stringsAsFactors = FALSE) %>%
  left_join(annot_raw[, c("gene", "family", "primary", "fn_score",
                          "in_set2", "in_set3", "set_bonus", "ubi_score")], by = "gene")
for (cc in c("family", "primary", "fn_score", "in_set2", "in_set3", "set_bonus", "ubi_score"))
  annot[[cc]][is.na(annot[[cc]])] <- if (cc %in% c("family", "primary")) "Unknown" else
    if (cc %in% c("in_set2", "in_set3")) FALSE else 0.3
cat("ubiquitin annotation ready, genes:", nrow(annot), "\n")# ---- full TCGA cancer names ----
CANCER_NAMES <- c(
  ACC = "Adrenocortical carcinoma", BLCA = "Bladder urothelial carcinoma",
  BRCA = "Breast invasive carcinoma", CESC = "Cervical squamous cell carcinoma",
  CHOL = "Cholangiocarcinoma", COAD = "Colon adenocarcinoma",
  DLBC = "Diffuse large B-cell lymphoma", ESCA = "Esophageal carcinoma",
  GBM = "Glioblastoma multiforme", HNSC = "Head and neck squamous cell carcinoma",
  KICH = "Kidney chromophobe", KIRC = "Kidney renal clear cell carcinoma",
  KIRP = "Kidney renal papillary cell carcinoma", LAML = "Acute myeloid leukemia",
  LGG = "Brain lower grade glioma", LIHC = "Liver hepatocellular carcinoma",
  LUAD = "Lung adenocarcinoma", LUSC = "Lung squamous cell carcinoma",
  MESO = "Mesothelioma", OV = "Ovarian serous cystadenocarcinoma",
  PAAD = "Pancreatic adenocarcinoma", PCPG = "Pheochromocytoma and paraganglioma",
  PRAD = "Prostate adenocarcinoma", READ = "Rectum adenocarcinoma",
  SARC = "Sarcoma", SKCM = "Skin cutaneous melanoma",
  STAD = "Stomach adenocarcinoma", TGCT = "Testicular germ cell tumors",
  THCA = "Thyroid carcinoma", THYM = "Thymoma",
  UCEC = "Uterine corpus endometrial carcinoma", UCS = "Uterine carcinosarcoma",
  UVM = "Uveal melanoma")

# ---- build one cancer cache ----
build_one <- function(ca) {
  d <- file.path(pd, paste0("TCGA-", ca))
  cat("\n[", ca, "] loading ...\n", sep = "")
  e <- new.env(); load(file.path(d, paste0("TCGA_", ca, "_matrix_tpm_processed.Rdata")), envir = e)
  expr_tpm_full <- encode_expr(e$expr_tpm)
  expr_tpm <- as.matrix(e$expr_tpm)
  expr_count <- load_matrix_processed(d, ca, "count", "expr_count")
  expr_fpkm  <- load_matrix_processed(d, ca, "fpkm", "expr_fpkm")
  st <- substr(colnames(expr_tpm), 14, 16)
  sample_info <- data.frame(sample = colnames(expr_tpm),
    type = ifelse(st %in% c("01A", "01B", "01C", "02A", "02B", "03A", "03B",
                           "05A", "06A", "06B", "07A", "07B", "07C", "01R"), "Tumor",
           ifelse(st %in% c("11A", "11B", "11C"), "Normal", "Other")), stringsAsFactors = FALSE)
  mx <- apply(expr_tpm, 1, max, na.rm = TRUE)
  expr_tpm <- expr_tpm[mx >= 0.5, , drop = FALSE]
  all_genes <- rownames(expr_tpm)
  cat("  tpm genes:", nrow(expr_tpm), " samples:", ncol(expr_tpm), "\n")

  universe <- sort(intersect(ubi_genes, rownames(expr_tpm)))
  expr_tpm <- expr_tpm[universe, , drop = FALSE]
  cat("  universe:", length(universe), "\n")

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
  cat("  fpkm genes:", length(fpkm_tumor_mean), "\n")

  # differential results (may be absent for some cancers)
  diff <- list()
  diffFile <- file.path(d, paste0(ca, "_diff.Rdata"))
  if (file.exists(diffFile)) {
    e <- new.env(); load(diffFile, envir = e)
    add <- function(tag, obj, gc, lfc, pc, pj) {
      if (!is.null(e[[obj]])) diff[[tag]] <<- std_diff(as.data.frame(e[[obj]]), gc, lfc, pc, pj)
    }
    add("limma",  "res_limma",     "rownames", "logFC",          "P.Value",    "adj.P.Val")
    add("edgeR",  "res_edgeR",     "rownames", "logFC",          "PValue",     "FDR")
    add("deseq2", "res_deseq2_df", "gene_id",  "log2FoldChange", "pvalue",     "padj")
    add("pair",   "res_pair_df",   "gene_id",  "log2FoldChange", "pvalue",     "padj")
  }
  diff <- diff[lengths(diff) > 0]
  diff_available <- length(diff) > 0
  cat("  diff:", if (diff_available) paste(names(diff), collapse = ",") else "NONE", "\n")  # clinical
  e <- new.env(); load(file.path(d, paste0("clin_", ca, "_pro.Rdata")), envir = e)
  clin <- as.data.frame(e$clin_pro) %>%
    dplyr::select(any_of(c("sample", "OS", "OS.time"))) %>%
    filter(OS %in% c(0, 1)) %>% distinct(sample, .keep_all = TRUE)

  # immune infiltration (5 methods) + metabolism GSVA (7 pathways)
  e <- new.env(); load(file.path(d, paste0(ca, "_immune_cells.Rdata")), envir = e)
  immune <- list(timer = if (is.null(e$timer_res)) NULL else as.data.frame(e$timer_res),
                 cibersort = if (is.null(e$cibersort_res)) NULL else as.data.frame(e$cibersort_res),
                 mcp = if (is.null(e$mcp_res)) NULL else as.data.frame(e$mcp_res),
                 ssgsea = if (is.null(e$ssgsea_res)) NULL else as.data.frame(e$ssgsea_res),
                 xcell = if (is.null(e$xcell_res)) NULL else as.data.frame(e$xcell_res))
  immune <- immune[!vapply(immune, is.null, logical(1))]
  tidy_immune <- function(df) {
    df <- df[!is.na(df$ID) & df$ID != "", , drop = FALSE]
    keep <- vapply(df, is.numeric, logical(1)); keep["ID"] <- TRUE
    df <- df[, keep, drop = FALSE]
    df[!duplicated(df$ID), , drop = FALSE]
  }
  immune <- lapply(immune, tidy_immune)
  e <- new.env(); load(file.path(d, paste0(ca, "_metabolism7_GSVA.Rdata")), envir = e)
  metabolism <- as.data.frame(e$GSVA_Metabolism)
  metabolism <- metabolism[!is.na(metabolism$ID) & metabolism$ID != "", , drop = FALSE]
  metabolism <- metabolism[!duplicated(metabolism$ID), , drop = FALSE]

  # tumor log2 TPM (patient IDs)
  tum <- sample_info$sample[sample_info$type == "Tumor"]
  tum <- tum[!duplicated(substr(tum, 1, 12))]
  expr_tumor_log2 <- log2(expr_tpm[, tum, drop = FALSE] + 1)
  colnames(expr_tumor_log2) <- substr(colnames(expr_tumor_log2), 1, 12)
  clin_tumor <- clin[clin$sample %in% colnames(expr_tumor_log2), , drop = FALSE]
  clin_tumor <- clin_tumor[!is.na(clin_tumor$OS.time) & is.finite(clin_tumor$OS.time), , drop = FALSE]
  expr_genes <- intersect(universe, rownames(expr_tumor_log2))
  cat("  tumor samples:", ncol(expr_tumor_log2), " with OS:", nrow(clin_tumor), "\n")  # Cox (per gene) -> logHR + p
  gene_df <- as.data.frame(t(expr_tumor_log2[expr_genes, , drop = FALSE]), check.names = FALSE)
  gene_df$sample <- rownames(gene_df)
  cox_df <- merge(clin_tumor, gene_df, by = "sample")
  cox_df <- cox_df[!is.na(cox_df$OS) & !is.na(cox_df$OS.time) & cox_df$OS.time > 0, , drop = FALSE]
  safe <- make.names(expr_genes); names(safe) <- expr_genes
  colnames(cox_df)[4:ncol(cox_df)] <- safe[colnames(cox_df)[4:ncol(cox_df)]]
  surv_df <- cox_df; surv_df$OS.time <- surv_df$OS.time / 365
  W <- list(safe = safe, cox_df = cox_df, surv_df = surv_df,
            expr_tumor_log2 = expr_tumor_log2, norm_log2 = NULL, diag_type = NULL)
  list2env(W, envir = globalenv())      # stage data in master globalenv (per-task closures pick it up)
  if (!is.null(cl)) clusterExport(cl, names(W), envir = globalenv())
  cx <- if (is.null(cl)) vapply(expr_genes, one_cox, numeric(2)) else parSapply(cl, expr_genes, one_cox)
  cat("  cox done, sig(p<0.05):", sum(cx[2, ] < 0.05, na.rm = TRUE), "\n")
  rc <- if (is.null(cl)) sapply(expr_genes, one_roc) else parSapply(cl, expr_genes, one_roc)
  cat("  timeROC done\n")

  # diagnostic ROC (tumor vs normal), only if normal samples exist
  nor <- sample_info$sample[sample_info$type == "Normal"]
  diag_auc <- setNames(rep(NA_real_, length(expr_genes)), expr_genes)
  if (length(nor) > 1 && ncol(expr_tumor_log2) > 1) {
    norm_log2 <- log2(expr_tpm[, nor, drop = FALSE] + 1)
    norm_log2 <- norm_log2[expr_genes, , drop = FALSE]
    diag_type <- c(rep("Normal", ncol(norm_log2)), rep("Tumor", ncol(expr_tumor_log2)))
    assign("norm_log2", norm_log2, envir = globalenv())
    assign("diag_type", diag_type, envir = globalenv())
    if (!is.null(cl)) clusterExport(cl, c("norm_log2", "diag_type"), envir = globalenv())
    dd <- if (is.null(cl)) vapply(expr_genes, one_diag, numeric(1)) else parSapply(cl, expr_genes, one_diag)
    diag_auc[expr_genes] <- dd
  }
  n_diag <- sum(!is.na(diag_auc))
  if (length(nor) > 1 && n_diag == 0) warning("diag ROC all NA for ", ca, " (normals=", length(nor), ")", call. = FALSE)
  cat("  diag ROC done, non-NA:", n_diag, "\n")  # correlations: immune (per method) + metabolism (7 pathways)
  cor_immune <- lapply(immune, function(dd) cor_one(expr_genes, dd, expr_tumor_log2))
  cor_meta   <- cor_one(expr_genes, metabolism, expr_tumor_log2)

  # differential percentile within expressed universe
  diff_pct <- lapply(diff, function(dd) {
    dd <- dd[dd$gene %in% expr_genes, , drop = FALSE]
    if (nrow(dd) == 0) return(NULL)
    r <- rank(-abs(dd$logFC), ties.method = "min")
    dd$pct <- r / nrow(dd)
    dd[, c("gene", "logFC", "p", "padj", "pct"), drop = FALSE]
  })

  fpkm_tumor_mean <- fpkm_tumor_mean[names(fpkm_tumor_mean) %in% expr_genes]

  # assemble + save
  cache <- list(
    cancer = ca, cancer_name = CANCER_NAMES[[ca]] %||% ca,
    sample_info = sample_info, expr_tpm = expr_tpm, all_genes = all_genes,
    expr_count = expr_count, expr_tpm_full = expr_tpm_full, expr_fpkm = expr_fpkm,
    clinical = clin, diff = diff, diff_available = diff_available,
    immune = immune, metabolism = metabolism, annotation = annot, universe = expr_genes,
    pre = list(
      fpkm_tumor_mean = fpkm_tumor_mean,
      cox = data.frame(gene = expr_genes, HR_log = unname(cx[1, ]), p = unname(cx[2, ]),
                       stringsAsFactors = FALSE),
      surv_auc = data.frame(gene = expr_genes, auc1 = unname(rc[1, ]), auc3 = unname(rc[2, ]),
                            auc5 = unname(rc[3, ]), stringsAsFactors = FALSE),
      diag_auc = data.frame(gene = expr_genes, auc = unname(diag_auc[expr_genes]),
                            stringsAsFactors = FALSE),
      diff_pct = diff_pct, cor_immune = cor_immune, cor_meta = cor_meta
    )
  )
  out <- file.path(outdir, paste0(ca, "_cache.rds"))
  saveRDS(cache, out, compress = "gzip")
  cat("  saved:", out, " size:", round(file.info(out)$size / 1e6, 1), "MB\n")
  invisible(TRUE)
}

# ---- main loop ----
args <- commandArgs(TRUE)
rebuild <- "--rebuild" %in% args
only <- unlist(strsplit(sub("^--only=", "", args[grepl("^--only=", args)]), ","))
cancers <- sort(sub("^TCGA-", "", list.dirs(pd, recursive = FALSE, full.names = FALSE)))
if (length(only)) cancers <- intersect(cancers, only)
n_cores <- as.integer(Sys.getenv("UBI_CORES", "8"))
cl <- NULL
if (n_cores > 1 && length(cancers) > 1) {
  cl <- makeCluster(min(n_cores, detectCores()), outfile = "")
  clusterEvalQ(cl, { .libPaths(c("D:/PackagesR/R_LIBS", .libPaths())); library(survival); library(timeROC); library(pROC) })
}
cat("cancers to build:", length(cancers), " cores:", if (is.null(cl)) 1 else length(cl), "\n")
for (ca in cancers) {
  out <- file.path(outdir, paste0(ca, "_cache.rds"))
  if (!rebuild && file.exists(out)) { cat("skip", ca, "(cache exists)\n"); next }
  tryCatch(build_one(ca), error = function(e) cat("!! FAILED", ca, ":", conditionMessage(e), "\n"))
}
if (!is.null(cl)) stopCluster(cl)
cat("\nALL DONE\n")
