# ============================================================================
# update_expr_caches.R
# Add full-gene count / tpm / fpkm expression matrices to existing caches.
# Values are stored as integer round(log2(x + 1) * 100) to halve the file
# size while keeping ~1% expression precision (enough for boxplots).
# Output: UbiPanTriage/inst/extdata/<CANCER>_cache.rds (fields expr_count,
#         expr_tpm_full, expr_fpkm). Resumable: finished cancers are skipped.
# Run from project root:
#   Rscript UbiPanTriage/data-raw/update_expr_caches.R
#   Rscript UbiPanTriage/data-raw/update_expr_caches.R --only=LUAD,BRCA
# ============================================================================
.libPaths(c("D:/PackagesR/R_LIBS", .libPaths()))
suppressMessages(library(parallel))

`%||%` <- function(x, y) if (is.null(x)) y else x
find_root <- function() {
  sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  cands <- c(getwd(), if (length(sf)) dirname(normalizePath(sf)) else NULL)
  cands <- c(cands, file.path(cands, ".."), file.path(cands, "../.."))
  cands <- unique(normalizePath(cands[nzchar(cands)]))
  hit <- cands[vapply(cands, function(p) dir.exists(file.path(p, "pancaner_data")), logical(1))]
  if (!length(hit)) stop("cannot locate pancaner_data/ ; run Rscript from project root")
  hit[1]
}
root   <- find_root()
pd     <- file.path(root, "pancaner_data")
outdir <- file.path(root, "UbiPanTriage", "inst", "extdata")

# integer-encoded log2(x+1)*100 matrix (genes x samples)
encode_expr <- function(df) {
  m <- as.matrix(df)
  storage.mode(m) <- "double"
  mx <- apply(m, 1, max, na.rm = TRUE)
  m <- m[is.finite(mx) & mx >= 0.5, , drop = FALSE]
  v <- round(log2(m + 1) * 100)
  storage.mode(v) <- "integer"
  v
}

load_matrix <- function(d, ca, unit) {
  f <- file.path(d, paste0("TCGA_", ca, "_matrix_", unit, "_processed.Rdata"))
  if (!file.exists(f)) return(NULL)
  e <- new.env(); load(f, envir = e)
  obj <- if (unit == "tpm") "expr_tpm" else if (unit == "count") "expr_count" else "expr_fpkm"
  if (is.null(e[[obj]])) return(NULL)
  encode_expr(e[[obj]])
}

update_one <- function(ca) {
  d <- file.path(pd, paste0("TCGA-", ca))
  cf <- file.path(outdir, paste0(ca, "_cache.rds"))
  if (!file.exists(cf)) { cat("skip", ca, "(no cache)\n"); return(invisible(FALSE)) }
  cache <- readRDS(cf)
  if (!is.null(cache$expr_tpm_full)) { cat("skip", ca, "(expr already added)\n"); return(invisible(FALSE)) }
  cat("\n[", ca, "] adding expression matrices ...\n", sep = "")
  for (unit in c("count", "tpm", "fpkm")) {
    m <- load_matrix(d, ca, unit)
    field <- switch(unit, count = "expr_count", tpm = "expr_tpm_full", fpkm = "expr_fpkm")
    if (is.null(m)) { cat("  ", unit, "MISSING\n"); next }
    cache[[field]] <- m
    cat("  ", unit, ":", nrow(m), "x", ncol(m), "\n")
  }
  cache$version <- paste(c(cache$version, "expr+count/tpm/fpkm"), collapse = ";")
  saveRDS(cache, cf, compress = "gzip")
  cat("  saved:", cf, " size:", round(file.info(cf)$size / 1e6, 1), "MB\n")
  invisible(TRUE)
}

args <- commandArgs(TRUE)
only <- unlist(strsplit(sub("^--only=", "", args[grepl("^--only=", args)]), ","))
cancers <- sort(sub("^TCGA-", "", list.dirs(pd, recursive = FALSE, full.names = FALSE)))
if (length(only)) cancers <- intersect(cancers, only)
cat("cancers to update:", length(cancers), "\n")
n_cores <- as.integer(Sys.getenv("UBI_CORES", "6"))
if (n_cores > 1 && length(cancers) > 1) {
  cl <- makeCluster(min(n_cores, detectCores()), outfile = "")
  clusterEvalQ(cl, { .libPaths(c("D:/PackagesR/R_LIBS", .libPaths())) })
  clusterExport(cl, c("pd", "outdir", "update_one", "encode_expr", "load_matrix"), envir = environment())
  invisible(parLapply(cl, cancers, update_one))
  stopCluster(cl)
} else {
  invisible(lapply(cancers, update_one))
}
cat("\nALL DONE\n")