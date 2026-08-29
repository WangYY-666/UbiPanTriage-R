# Build lookup/ubi_univ.rds (long ubiquitin-gene score table, new scoring) and
# lookup/pcs_meta.rds (cancer_meta + annotation only). Run after the per-cancer
# partial files exist (build_web_lookup.R merge step). Read-only on inputs.
.libPaths(c("D:/PackagesR/R_LIBS", .libPaths()))
root <- "E:/LIANXI_R/Bulk RNAseq/UBI Cancer/Immunometabolism_ubiquitin_LUAD/Ubi_panTriage R"
lookup <- file.path(root, "UbiPanTriage", "inst", "extdata", "lookup")
partdir <- file.path(lookup, ".partial")
pfs <- list.files(partdir, pattern = "\\.rds$", full.names = TRUE)
if (length(pfs) == 0) stop("no partial files")
cas <- sort(sub("\\.rds$", "", basename(pfs)))
# universe = genes with a ubiquitin score in at least one cancer (new build)
uni <- lapply(pfs, function(f) {
  p <- readRDS(f); s <- p$scores
  s$Gene[!is.na(s$Ubi_Score)]
})
univ <- sort(unique(unlist(uni)))
cat("ubiquitin universe genes:", length(univ), " cancers:", length(cas), "\n")
blocks <- lapply(cas, function(ca) {
  p <- readRDS(file.path(partdir, paste0(ca, ".rds")))
  s <- p$scores[p$scores$Gene %in% univ, , drop = FALSE]
  s$Cancer <- ca
  s
})
allc <- unique(unlist(lapply(blocks, colnames)))
blocks <- lapply(blocks, function(b) {
  miss <- setdiff(allc, colnames(b))
  for (nm in miss) b[[nm]] <- NA
  b[, allc, drop = FALSE]
})
long <- do.call(rbind, blocks)
rownames(long) <- NULL
cat("ubi_univ dim:", dim(long), "\n")
saveRDS(long, file.path(lookup, "ubi_univ.rds"), compress = "gzip")
old <- readRDS(file.path(root, "UbiPanTriage", "inst", "extdata", "pancancer_ubi_scores.rds"))
saveRDS(list(version = old$version, built = old$built,
             cancer_meta = old$cancer_meta, annotation = old$annotation),
        file.path(lookup, "pcs_meta.rds"), compress = "gzip")
cat("pcs_meta saved. DONE\n")
