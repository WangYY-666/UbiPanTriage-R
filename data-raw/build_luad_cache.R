# ============================================================================
# build_luad_cache.R
# Build the LUAD local cache used by the UbiPanTriage package.
# ?? LUAD ??????? inst/extdata/LUAD_cache.rds?
# Run: Rscript data-raw/build_luad_cache.R  (from project root)
# ============================================================================
# local R library (adjust to your machine)
.libPaths(c("D:/PackagesR/R_LIBS", .libPaths()))
suppressMessages({
  library(dplyr); library(tidyr); library(tibble); library(stringr)
  library(survival); library(pROC); library(timeROC); library(readxl)
})

cands <- c(normalizePath("../../"), normalizePath(".."), normalizePath("."))
root  <- cands[vapply(cands, function(x) dir.exists(file.path(x, "example")), logical(1))][1]
example <- file.path(root, "example")
outfile <- file.path(root, "UbiPanTriage", "inst", "extdata", "LUAD_cache.rds")
stopifnot(dir.exists(example))

cat("[1/8] loading gene sets ...\n")
e <- new.env(); load(file.path(example, "ubiquitin_genes.RData"), envir = e)
ubi_genes   <- e$ubiquitin_genes
set2 <- e$ubiquitin_genes2; set3 <- e$ubiquitin_genes3

cat("[2/8] loading expression (TPM) ...\n")
e <- new.env(); load(file.path(example, "TCGA_LUAD_matrix_tpm_processed.Rdata"), envir = e)
expr_tpm <- as.matrix(e$expr_tpm)                    # symbol x sample
samp_type <- substr(colnames(expr_tpm), 14, 16)
sample_info <- data.frame(
  sample = colnames(expr_tpm),
  type   = ifelse(samp_type %in% c("01A","01B","01C","02A"), "Tumor",
            ifelse(samp_type %in% c("11A","11B"), "Normal", "Other")),
  stringsAsFactors = FALSE
)
max_tpm <- apply(expr_tpm, 1, max, na.rm = TRUE)     # keep reasonably expressed genes
expr_tpm <- expr_tpm[max_tpm >= 0.5, , drop = FALSE]
cat("  genes kept:", nrow(expr_tpm), " samples:", ncol(expr_tpm), "\n")

cat("[3/8] loading FPKM (expression score) ...\n")
e <- new.env(); load(file.path(example, "TCGA_LUAD_matrix_fpkm.Rdata"), envir = e)
mat_fpkm <- e$mat_fpkm                              # Ensembl x sample, col1 = symbol
tumor_samples <- sample_info$sample[sample_info$type == "Tumor"]
fpkm_t <- mat_fpkm[, c("symbol", tumor_samples), drop = FALSE]
fpkm_t[, -1] <- lapply(fpkm_t[, -1, drop = FALSE], function(x) as.numeric(x))
fpkm_t <- fpkm_t[!is.na(fpkm_t$symbol) & fpkm_t$symbol != "", , drop = FALSE]
fpkm_mean <- rowMeans(as.matrix(fpkm_t[, -1]), na.rm = TRUE)
names(fpkm_mean) <- fpkm_t$symbol
fpkm_tumor_mean <- tapply(fpkm_mean, names(fpkm_mean), mean)   # dedupe by symbol

cat("[4/8] standardizing diff results ...\n")
e <- new.env(); load(file.path(example, "LUAD_diff.Rdata"), envir = e)
std_diff <- function(df, gene_col, logfc_col, p_col, padj_col) {
  if (gene_col == "rownames") df$gene <- rownames(df) else df$gene <- df[[gene_col]]
  out <- data.frame(gene = as.character(df$gene),
                    logFC = as.numeric(df[[logfc_col]]),
                    p = as.numeric(df[[p_col]]),
                    padj = as.numeric(df[[padj_col]]), stringsAsFactors = FALSE)
  out[!is.na(out$gene) & is.finite(out$logFC), , drop = FALSE]
}
diff <- list(
  limma  = std_diff(as.data.frame(e$res_limma),     "rownames", "logFC", "P.Value",   "adj.P.Val"),
  edgeR  = std_diff(as.data.frame(e$res_edgeR),     "rownames", "logFC", "PValue",    "FDR"),
  deseq2 = std_diff(as.data.frame(e$res_deseq2_df), "gene_id",  "log2FoldChange", "pvalue", "padj"),
  pair   = std_diff(as.data.frame(e$res_pair_df),   "gene_id",  "log2FoldChange", "pvalue", "padj")
)
for (m in names(diff)) cat("  ", m, ":", nrow(diff[[m]]), "genes\n")

cat("[5/8] clinical + immune + metabolism ...\n")
e <- new.env(); load(file.path(example, "clin_LUAD_pro.Rdata"), envir = e)
clin <- as.data.frame(e$clin_pro) %>%
  dplyr::select(sample, OS, OS.time) %>%
  filter(OS %in% c(0, 1)) %>%
  distinct(sample, .keep_all = TRUE)

e <- new.env(); load(file.path(example, "LUAD_immune_cells.Rdata"), envir = e)
immune <- list(
  timer     = as.data.frame(e$timer_res),
  cibersort = as.data.frame(e$cibersort_res),
  mcp       = as.data.frame(e$mcp_res),
  ssgsea    = as.data.frame(e$ssgsea_res),
  xcell     = as.data.frame(e$xcell_res)
)
e <- new.env(); load(file.path(example, "LUAD_metabolism7_GSVA.Rdata"), envir = e)
metabolism <- as.data.frame(e$GSVA_Metabolism)

tidy_immune <- function(df) {
  df <- df[!is.na(df$ID) & df$ID != "", , drop = FALSE]
  keep <- vapply(df, is.numeric, logical(1)); keep["ID"] <- TRUE
  df <- df[, keep, drop = FALSE]
  df[!duplicated(df$ID), , drop = FALSE]
}
immune <- lapply(immune, tidy_immune)
metabolism <- metabolism[!is.na(metabolism$ID) & metabolism$ID != "", , drop = FALSE]
metabolism <- metabolism[!duplicated(metabolism$ID), , drop = FALSE]

cat("[6/8] ubiquitin annotation (IUCCD2.0) ...\n")
xlsx <- file.path(example, "IUCCD2.0 ubiquitin genes.xlsx")
raw  <- as.data.frame(read_excel(xlsx, sheet = 1, col_names = FALSE), stringsAsFactors = FALSE)
dat  <- raw[-(1:2), , drop = FALSE]                 # row1 = title, row2 = header
annot_raw <- data.frame(
  gene     = as.character(dat[[3]]),                # Gene name
  uniprot  = as.character(dat[[2]]),
  ensembl  = as.character(dat[[4]]),
  family   = as.character(dat[[5]]),
  species  = as.character(dat[[6]]),
  stringsAsFactors = FALSE
)
annot_raw <- annot_raw[!is.na(annot_raw$gene) & annot_raw$species == "Homo sapiens", , drop = FALSE]
annot_raw <- annot_raw[!duplicated(annot_raw$gene), , drop = FALSE]

prio <- c("E1", "E2", "E3 activity", "E3 adaptor", "DUB", "UBD", "ULD")
extract_primary <- function(fam) {
  if (is.na(fam) || fam == "") return(NA_character_)
  cls <- strsplit(fam, "; ", fixed = TRUE)[[1]]
  for (p in prio) {
    hit <- cls[startsWith(cls, p)]
    if (length(hit) > 0) return(hit[1])
  }
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
  "ULD/UDP/UFD/PB1" = 0.40, "ULD/UDP/UFD/UBX" = 0.42, "ULD/UBL/SUMO" = 0.38, "ULD/UBL/ATG8" = 0.36
)
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
  left_join(annot_raw[, c("gene","family","primary","fn_score","in_set2","in_set3","set_bonus","ubi_score")], by = "gene")
annot$family[is.na(annot$family)] <- "Unknown"
annot$primary[is.na(annot$primary)] <- "Unknown"
annot$fn_score[is.na(annot$fn_score)] <- 0.3
annot$in_set2[is.na(annot$in_set2)] <- FALSE
annot$in_set3[is.na(annot$in_set3)] <- FALSE
annot$set_bonus[is.na(annot$set_bonus)] <- 0
annot$ubi_score[is.na(annot$ubi_score)] <- 0.3

# ---------------------------------------------------------------------------
cat("[7/8] precomputing per-gene stats ...\n")
universe <- sort(intersect(ubi_genes, rownames(expr_tpm)))
cat("  universe genes:", length(universe), "\n")

tum <- sample_info$sample[sample_info$type == "Tumor"]
tum_pid <- substr(tum, 1, 12)
tum_samp <- tum[!duplicated(tum_pid)]
expr_tumor_log2 <- log2(expr_tpm[, tum_samp, drop = FALSE] + 1)
colnames(expr_tumor_log2) <- substr(colnames(expr_tumor_log2), 1, 12)

clin_tumor <- clin[clin$sample %in% colnames(expr_tumor_log2), , drop = FALSE]
clin_tumor <- clin_tumor[!is.na(clin_tumor$OS.time) & is.finite(clin_tumor$OS.time), , drop = FALSE]
cat("  tumor samples:", ncol(expr_tumor_log2), " with OS data:", nrow(clin_tumor), "\n")

expr_genes <- intersect(universe, rownames(expr_tumor_log2))
expr_gene_df <- as.data.frame(t(expr_tumor_log2[expr_genes, , drop = FALSE]), check.names = FALSE)
expr_gene_df$sample <- rownames(expr_gene_df)
cox_df <- merge(clin_tumor, expr_gene_df, by = "sample")
cox_df <- cox_df[!is.na(cox_df$OS) & !is.na(cox_df$OS.time) & cox_df$OS.time > 0, , drop = FALSE]

safe <- make.names(expr_genes); names(safe) <- expr_genes
colnames(cox_df)[4:ncol(cox_df)] <- safe[colnames(cox_df)[4:ncol(cox_df)]]

cox_p <- setNames(rep(NA_real_, length(expr_genes)), expr_genes)
cox_hr <- setNames(rep(NA_real_, length(expr_genes)), expr_genes)
for (g in expr_genes) {
  f <- as.formula(paste("Surv(OS.time, OS) ~", safe[g]))
  fit <- tryCatch(coxph(f, data = cox_df), error = function(e) NULL)
  if (!is.null(fit)) {
    s <- summary(fit)
    cox_p[g]  <- s$coefficients[1, "Pr(>|z|)"]
    cox_hr[g] <- s$coefficients[1, "coef"]           # log HR
  }
}
cat("  cox done, significant(p<0.05):", sum(cox_p < 0.05, na.rm = TRUE), "\n")

surv_df <- cox_df; surv_df$OS.time <- surv_df$OS.time / 365
auc1 <- auc3 <- auc5 <- setNames(rep(NA_real_, length(expr_genes)), expr_genes)
for (g in expr_genes) {
  for (tp in c(1, 3, 5)) {
    r <- tryCatch(timeROC(T = surv_df$OS.time, delta = surv_df$OS, marker = surv_df[[safe[g]]],
                          cause = 1, weighting = "aalen", time = tp, ROC = FALSE),
                  error = function(e) NULL)
    if (!is.null(r)) {
      a <- as.numeric(r$AUC[2])
      if (tp == 1) auc1[g] <- a else if (tp == 3) auc3[g] <- a else auc5[g] <- a
    }
  }
}
cat("  timeROC done\n")

nor <- sample_info$sample[sample_info$type == "Normal"]
diag_auc <- setNames(rep(NA_real_, length(expr_genes)), expr_genes)
if (length(nor) > 0) {
  norm_log2 <- log2(expr_tpm[, nor, drop = FALSE] + 1)
  type <- c(rep("Normal", ncol(norm_log2)), rep("Tumor", ncol(expr_tumor_log2)))
  for (g in expr_genes) {
    vals <- c(as.numeric(norm_log2[g, ]), as.numeric(expr_tumor_log2[g, ]))
    rr <- tryCatch(roc(type ~ vals, levels = c("Normal", "Tumor"), direction = ">", quiet = TRUE),
                   error = function(e) NULL)
    diag_auc[g] <- if (!is.null(rr)) as.numeric(rr$auc) else NA_real_
  }
}
cat("  diagnostic ROC done\n")

cor_one <- function(xgenes, ydf) {
  ydf$pid <- substr(ydf$ID, 1, 12)
  ydf <- ydf[!duplicated(ydf$pid), , drop = FALSE]
  common <- intersect(colnames(expr_tumor_log2), ydf$pid)
  if (length(common) < 10) return(NULL)
  x <- expr_tumor_log2[xgenes, common, drop = FALSE]
  y <- as.matrix(ydf[match(common, ydf$pid), setdiff(colnames(ydf), c("ID", "pid")), drop = FALSE])
  x <- x[apply(x, 1, function(v) sd(v) > 0), , drop = FALSE]
  y <- y[, apply(y, 2, function(v) sd(v) > 0), drop = FALSE]
  if (nrow(x) == 0 || ncol(y) == 0) return(NULL)
  cor(t(x), y, method = "spearman")
}
cor_immune <- lapply(immune, function(d) cor_one(expr_genes, d))
cor_meta <- cor_one(expr_genes, metabolism)

fpkm_tumor_mean <- fpkm_tumor_mean[names(fpkm_tumor_mean) %in% expr_genes]

diff_pct <- lapply(diff, function(d) {
  d <- d[d$gene %in% expr_genes, , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  r <- rank(-abs(d$logFC), ties.method = "min")
  d$pct <- r / nrow(d)
  d[, c("gene", "logFC", "p", "padj", "pct"), drop = FALSE]
})

cat("[8/8] assembling cache ...\n")
cache <- list(
  cancer      = "LUAD",
  cancer_name = "Lung adenocarcinoma",
  sample_info = sample_info,
  expr_tpm    = expr_tpm,
  clinical    = clin,
  diff        = diff,
  immune      = immune,
  metabolism  = metabolism,
  annotation  = annot,
  universe    = expr_genes,
  pre = list(
    fpkm_tumor_mean = fpkm_tumor_mean,
    cox = data.frame(gene = expr_genes, HR_log = unname(cox_hr[expr_genes]),
                     p = unname(cox_p[expr_genes]), stringsAsFactors = FALSE),
    surv_auc = data.frame(gene = expr_genes, auc1 = unname(auc1[expr_genes]),
                          auc3 = unname(auc3[expr_genes]), auc5 = unname(auc5[expr_genes]),
                          stringsAsFactors = FALSE),
    diag_auc = data.frame(gene = expr_genes, auc = unname(diag_auc[expr_genes]),
                          stringsAsFactors = FALSE),
    diff_pct = diff_pct,
    cor_immune = cor_immune,
    cor_meta   = cor_meta
  )
)
dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
saveRDS(cache, outfile, compress = "gzip")
cat("cache written:", outfile, " size:", round(file.info(outfile)$size / 1e6, 1), "MB\n")

