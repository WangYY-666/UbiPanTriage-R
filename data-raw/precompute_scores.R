# ============================================================================
# precompute_scores.R
# Pre-compute the four-dimension scores for every ubiquitin gene x 33 cancers.
# Output: UbiPanTriage/inst/extdata/pancancer_ubi_scores.rds
# v0.3.1: keeps ALL sub-scores (Expression/Diff/Survive/ROC, per-method
#   immune traits, 7 metabolic pathways), per-dimension ranks within each
#   cancer, and a gene-level pan-cancer summary (rank ranges + top cancers).
# The Shiny app reads this table for instant gene-set ranking queries.
# Run from project root:
#   Rscript UbiPanTriage/data-raw/precompute_scores.R
# ============================================================================
.libPaths(c("E:/LIANXI_R/Bulk RNAseq/UBI Cancer/Immunometabolism_ubiquitin_LUAD/Ubi_panTriage R/.Rlibs",
            "D:/PackagesR/R_LIBS", .libPaths()))
suppressMessages({ library(UbiPanTriage); library(parallel) })

root <- "E:/LIANXI_R/Bulk RNAseq/UBI Cancer/Immunometabolism_ubiquitin_LUAD/Ubi_panTriage R"
outdir <- file.path(root, "UbiPanTriage", "inst", "extdata")
data_dir <- outdir
cancers <- sort(available_cancers(data_dir))
cat("cancers:", length(cancers), "\n")

n_cores <- as.integer(Sys.getenv("UBI_CORES", "8"))
cl <- NULL
if (n_cores > 1) {
  cl <- makeCluster(min(n_cores, detectCores()), outfile = "")
  clusterEvalQ(cl, { .libPaths(c("E:/LIANXI_R/Bulk RNAseq/UBI Cancer/Immunometabolism_ubiquitin_LUAD/Ubi_panTriage R/.Rlibs",
                                 "D:/PackagesR/R_LIBS", .libPaths())); library(UbiPanTriage) })
}

one <- function(ca, dd) {
  cache <- load_cancer_data(ca, dd)
  genes <- cache$universe
  pref <- c("timer", "cibersort", "mcp", "ssgsea", "xcell")
  hit <- intersect(pref, names(cache$immune))
  imm_method <- if (length(hit) == 0) "timer" else hit[1]
  ip <- default_immune_params()
  ip$method <- imm_method
  r <- score_genes(cache, genes, immune_params = ip)
  sc <- r$scores
  # merge per-trait immune columns (default method) with a method prefix
  if (!is.null(r$immune)) {
    imm_df <- r$immune
    ic <- imm_df[, grep("^Imm_", colnames(imm_df)), drop = FALSE]
    if (ncol(ic) > 0) {
      colnames(ic) <- paste0("Imm_", imm_method, "_", sub("^Imm_", "", colnames(ic)))
      ic$gene <- imm_df$gene
      sc <- merge(sc, ic, by.x = "Gene", by.y = "gene", all.x = TRUE, sort = FALSE)
    }
  }
  # merge per-trait metabolic columns (7 GSVA pathways)
  if (!is.null(r$metabolic)) {
    mt_df <- r$metabolic
    mc <- mt_df[, grep("^Meta_", colnames(mt_df)), drop = FALSE]
    if (ncol(mc) > 0) {
      mc$gene <- mt_df$gene
      sc <- merge(sc, mc, by.x = "Gene", by.y = "gene", all.x = TRUE, sort = FALSE)
    }
  }
  # add immune trait columns of the other available methods (instant method switch)
  for (m in setdiff(hit, imm_method)) {
    cm <- tryCatch(immune_cor(cache, genes, m), error = function(e) NULL)
    if (is.null(cm) || ncol(cm) == 0) next
    colnames(cm) <- paste0("Imm_", m, "_", colnames(cm))
    cm <- as.data.frame(cm[genes, , drop = FALSE])
    sc <- cbind(sc, cm)
  }
  # basic/immune/metabolic are mapped to 0-1 percentile ranks within the
  # cancer (best = 1, worst = 1/N) so the three dimensions are comparable
  # across cancers and with the ubiquitin tier score (0-1)
  sc$Basic_Score     <- rank_to_pct(sc$Basic_Score)
  sc$Immune_Score    <- rank_to_pct(sc$Immune_Score)
  sc$Metabolic_Score <- rank_to_pct(sc$Metabolic_Score)
  # per-dimension ranks within this cancer (higher score = rank 1)
  sc$Ubi_Rank      <- rank(-sc$Ubi_Score, ties.method = "min")
  sc$Basic_Rank    <- rank(-sc$Basic_Score, ties.method = "min")
  sc$Immune_Rank   <- rank(-sc$Immune_Score, ties.method = "min")
  sc$Metabolic_Rank <- rank(-sc$Metabolic_Score, ties.method = "min")
  sc$Cancer <- ca
  sc$n_tumor <- sum(cache$sample_info$type == "Tumor")
  sc$n_normal <- sum(cache$sample_info$type == "Normal")
  sc$has_diff <- isTRUE(cache$diff_available)
  sc$has_normal <- sc$n_normal > 0
  sc$Immune_Method <- imm_method
  sc
}
res <- if (is.null(cl)) lapply(cancers, one, data_dir) else parLapply(cl, cancers, one, dd = data_dir)
if (!is.null(cl)) stopCluster(cl)

# unify columns (cancers use different immune methods -> different trait columns)
all_cols <- unique(unlist(lapply(res, colnames)))
res2 <- lapply(res, function(x) {
  miss <- setdiff(all_cols, colnames(x))
  if (length(miss) > 0) for (cc in miss) x[[cc]] <- NA
  x[, all_cols, drop = FALSE]
})
scores <- do.call(rbind, res2)
rownames(scores) <- NULL

# ---- annotation table (gene-level, cancer-invariant) ----
ann <- unique(scores[, c("Gene", "Ubi_Type", "Ubi_Type_Full", "Ubi_Score")])

# ---- pan-cancer summary: means, rank ranges and top-cancer lists ----
dim_cols <- c("Ubi_Score", "Basic_Score", "Immune_Score", "Metabolic_Score")
summ <- aggregate(scores[, dim_cols, drop = FALSE],
                  by = list(Gene = scores$Gene), mean, na.rm = TRUE)
rank_cols <- c(Ubi = "Ubi_Rank", Basic = "Basic_Rank", Immune = "Immune_Rank",
               Metabolic = "Metabolic_Rank")
for (nm in names(rank_cols)) {
  rc <- rank_cols[[nm]]
  lo <- aggregate(list(v = scores[[rc]]), by = list(Gene = scores$Gene), min, na.rm = TRUE)
  hi <- aggregate(list(v = scores[[rc]]), by = list(Gene = scores$Gene), max, na.rm = TRUE)
  summ[[paste0(nm, "_Rank_Min")]] <- lo$v[match(summ$Gene, lo$Gene)]
  summ[[paste0(nm, "_Rank_Max")]] <- hi$v[match(summ$Gene, hi$Gene)]
}
# per-gene top-3 cancers by each dimension
pick_top <- function(col, n = 3) {
  g <- split(scores, scores$Gene)
  vapply(g, function(d) {
    d <- d[!is.na(d[[col]]), , drop = FALSE]
    if (nrow(d) == 0) return("")
    d <- d[order(-d[[col]]), , drop = FALSE]
    paste(head(d$Cancer, n), collapse = ",")
  }, character(1))
}
for (nm in c("Basic", "Immune", "Metabolic")) {
  summ[[paste0(nm, "_Top")]] <- pick_top(paste0(nm, "_Score"))[summ$Gene]
}
summ$Combined_Score <- rowMeans(summ[, c("Basic_Score", "Immune_Score", "Metabolic_Score")], na.rm = TRUE)
summ$Pan_Rank <- rank(-summ$Combined_Score, ties.method = "min")
summ <- summ[order(summ$Pan_Rank), , drop = FALSE]
summ <- merge(summ, ann, by = "Gene", all.x = TRUE, sort = FALSE)

cancer_meta <- unique(scores[, c("Cancer", "n_tumor", "n_normal", "has_diff", "has_normal", "Immune_Method")])

out <- list(version = "0.3.1",
            built = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            annotation = ann,
            scores = scores,
            summary = summ,
            cancer_meta = cancer_meta)
saveRDS(out, file.path(outdir, "pancancer_ubi_scores.rds"), compress = "gzip")
cat("saved pancancer_ubi_scores.rds:", round(file.info(file.path(outdir, "pancancer_ubi_scores.rds"))$size / 1e6, 1), "MB\n")
cat("rows:", nrow(scores), " genes:", length(unique(scores$Gene)),
    " cols:", ncol(scores), "\n")
cat("ALL DONE\n")