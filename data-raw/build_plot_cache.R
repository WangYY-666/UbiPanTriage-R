# ============================================================================
# build_plot_cache.R - pre-render per-cancer plots (box / roc / km) as small
# PNG files into inst/extdata/plotcache, so the web app can compose them into
# grids without recomputing survminer / pROC / boxplots on the (2-core 4GB)
# server. PNG (~0.1-0.2 MB each) is memory-light; caching ggplot objects as
# RDS was tested and rejected (~16 MB RAM per object, too heavy for 4GB).
#
# Usage (from the project root):
#   Rscript UbiPanTriage/data-raw/build_plot_cache.R                  # demo list
#   Rscript UbiPanTriage/data-raw/build_plot_cache.R --all-ubi        # all ubiquitin genes
#   Rscript UbiPanTriage/data-raw/build_plot_cache.R --genes TP53,EGFR --cancers LUAD,BRCA
#   Rscript UbiPanTriage/data-raw/build_plot_cache.R --out /path/to/cache
#   Rscript UbiPanTriage/data-raw/build_plot_cache.R --kinds box,roc   # subset
#
# Existing files are skipped, so re-running resumes where it stopped.
# Disk estimate per gene x cancer: ~0.2 MB (box + roc + km). A full ubiquitin
# pre-render (1503 genes x 33 cancers) is ~30 GB and is NOT recommended; the
# app's runtime write-through cache covers the rest on demand.
# ============================================================================

# --- locate package root (data-raw/..) ---
.args  <- commandArgs(FALSE)
.file  <- sub("^--file=", "", .args[grepl("^--file=", .args)])
root <- if (length(.file)) dirname(normalizePath(.file)) else getwd()
root <- dirname(root)
data_dir <- file.path(root, "inst", "extdata")
if (!dir.exists(data_dir)) stop("package inst/extdata not found under ", root)

# --- R package library: reuse the project local lib when present ---
proj <- dirname(root)
for (lp in c(file.path(proj, ".Rlibs"), "D:/PackagesR/R_LIBS"))
  if (dir.exists(lp)) .libPaths(c(lp, .libPaths()))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  i <- grep(paste0("^--", name, "="), args)
  if (length(i)) sub(paste0("^--", name, "="), "", args[i[1]]) else default
}
has_flag <- function(name) paste0("--", name) %in% args

out_dir  <- get_arg("out", file.path(data_dir, "plotcache"))
cancers  <- get_arg("cancers", NULL)
kinds    <- strsplit(get_arg("kinds", "box,roc,km"), ",")[[1]]
kinds    <- intersect(kinds, c("box", "roc", "km"))
if (length(kinds) == 0) stop("no valid kinds (box/roc/km)")

# --- dependencies ---
suppressPackageStartupMessages(library(UbiPanTriage))

# --- lookup data ---
lookup <- file.path(data_dir, "lookup")
si       <- readRDS(file.path(lookup, "sample_info.rds"))
clinical <- readRDS(file.path(lookup, "clinical.rds"))
fm       <- readRDS(file.path(lookup, "file_map.rds"))
if (is.null(cancers)) cancers <- sort(names(si))
cancers <- intersect(cancers, names(si))
if (length(cancers) == 0) stop("no valid cancers selected")

# --- gene list ---
demo_genes <- c("TRIM44", "USP7", "MDM2", "TP53", "EGFR", "MDM4", "UBE2C",
                "UBE2D1", "BRCA1", "BRCA2", "CDH1", "KRAS", "PTEN", "VHL",
                "PIK3CA", "MYC", "STAT3", "B2M")
if (has_flag("all-ubi")) {
  ubi <- readRDS(file.path(lookup, "ubi_univ.rds"))
  genes <- sort(unique(ubi$Gene))
} else if (!is.null(get_arg("file", NULL))) {
  genes <- trimws(readLines(get_arg("file", "")))
} else if (!is.null(get_arg("genes", NULL))) {
  genes <- trimws(unlist(strsplit(get_arg("genes", ""), "[,; \t]+")))
} else {
  genes <- demo_genes
}
genes <- unique(genes[nzchar(genes)])
genes <- intersect(genes, sort(unique(readRDS(file.path(lookup, "genes_index.rds"))$Gene)))
if (length(genes) == 0) stop("no genes found in the lookup library")
cat("genes:", length(genes), "| cancers:", length(cancers),
    "| kinds:", paste(kinds, collapse = ","), "\n")
cat("expected disk (all kinds):", sprintf("%.2f GB",
    length(genes) * length(cancers) * 0.2 / 1024), "\n\n")

# --- per-gene cache constructor (same as the app) ---
get_gene_data <- function(g) {
  fn <- if (!is.null(fm) && g %in% names(fm)) fm[[g]] else g
  readRDS(file.path(lookup, "gene", paste0(fn, ".rds")))
}
.read_i16 <- function(r) {
  if (is.null(r)) return(NULL)
  v <- readBin(r, integer(), n = length(r) %/% 2L, size = 2L, signed = TRUE)
  v[v == -32768L] <- NA_integer_
  v
}
.gene_cache <- function(g, ca) {
  gd <- get_gene_data(g)
  pc <- gd$per_cancer[[ca]]
  if (is.null(pc)) return(NULL)
  s <- si[[ca]]
  one <- function(v, nm) if (is.null(v)) NULL else matrix(v, nrow = 1, dimnames = list(nm, s$sample))
  tpm_i <- .read_i16(pc$expr$tpm)
  tpm <- tpm_i / 100
  pre_imm <- list()
  for (m in c("timer", "cibersort", "mcp", "ssgsea", "xcell")) {
    cols <- grep(paste0("^Imm_", m, "_"), colnames(pc$scores), value = TRUE)
    if (length(cols))
      pre_imm[[m]] <- matrix(as.numeric(pc$scores[1, cols]), nrow = 1,
                             dimnames = list(g, sub(paste0("^Imm_", m, "_"), "", cols)))
  }
  pre_meta <- NULL
  mcols <- grep("^Meta_", colnames(pc$scores), value = TRUE)
  if (length(mcols))
    pre_meta <- matrix(as.numeric(pc$scores[1, mcols]), nrow = 1,
                       dimnames = list(g, sub("^Meta_", "", mcols)))
  list(cancer = ca, sample_info = s, clinical = clinical[[ca]],
       expr_tpm = one(2^tpm - 1, g), expr_tpm_full = one(tpm_i, g),
       expr_count = one(.read_i16(pc$expr$count), g),
       expr_fpkm  = one(.read_i16(pc$expr$fpkm), g),
       pre = list(cor_immune = pre_imm, cor_meta = pre_meta))
}

# --- render one plot kind to a PNG file (returns the path or NULL) ---
render_one <- function(kind, cache, gene, f) {
  obj <- switch(kind,
                box = tryCatch(plot_expr_box(cache, gene), error = function(e) NULL),
                roc = tryCatch(plot_roc_diag(cache, gene), error = function(e) NULL),
                km  = tryCatch({ p <- plot_km(cache, gene); p$plot }, error = function(e) NULL))
  if (is.null(obj)) return(NULL)
  grDevices::png(f, width = 3.9, height = 3.3, units = "in", res = 130)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(obj)
  f
}

# --- main loop (skip existing) ---
dir.create(out_dir, recursive = TRUE)
t_start <- Sys.time()
done <- 0L; skipped <- 0L; failed <- 0L
for (gi in seq_along(genes)) {
  g <- genes[gi]
  for (ca in cancers) {
    cache <- tryCatch(.gene_cache(g, ca), error = function(e) NULL)
    if (is.null(cache)) next
    for (k in kinds) {
      f <- file.path(out_dir, k, paste0(g, "__", ca, ".png"))
      if (file.exists(f)) { skipped <- skipped + 1L; next }
      ok <- tryCatch(!is.null(render_one(k, cache, g, f)), error = function(e) FALSE)
      if (ok) done <- done + 1L else failed <- failed + 1L
    }
  }
  if (gi %% 5 == 0 || gi == length(genes))
    cat(sprintf("  %4d/%d genes done | cached %d skipped %d failed %d | elapsed %.1f min\n",
                gi, length(genes), done, skipped, failed,
                as.numeric(Sys.time() - t_start, units = "mins")))
}
el <- as.numeric(Sys.time() - t_start, units = "mins")
mb <- sum(file.info(list.files(out_dir, recursive = TRUE, full.names = TRUE))$size) / 1024^2
cat(sprintf("\nDONE: %d cached, %d skipped, %d failed in %.1f min | total %.1f MB -> %s\n",
            done, skipped, failed, el, mb, out_dir))
