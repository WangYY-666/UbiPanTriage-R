# ============================================================================
# make_readme_figures.R
# Regenerate the README example figures with the CURRENT package code so the
# gallery matches the Shiny app outputs (page 1: gene-set ranking; page 2:
# single-gene query). Run from the package root:
#   Rscript data-raw/make_readme_figures.R
# Requires the installed package (remotes::install_local or R CMD INSTALL)
# and the 33-cancer caches (inst/extdata, e.g. via download_extdata()).
# ============================================================================
suppressPackageStartupMessages(library(UbiPanTriage))
suppressPackageStartupMessages(library(ggplot2))

root <- normalizePath(file.path("..", ".."), winslash = "/")   # workspace root (run from data-raw/)
data_dir <- file.path(root, "UbiPanTriage", "inst", "extdata")
if (!dir.exists(data_dir)) stop("caches not found: ", data_dir)
fig_dir  <- file.path("..", "man", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
if (length(list.files(fig_dir, pattern = "\\.png$", full.names = TRUE)) > 0)
  file.remove(list.files(fig_dir, pattern = "\\.png$", full.names = TRUE))

## ---------- data ----------
cancers <- c("LUAD", "BRCA", "LIHC", "COAD", "KIRC")
caches  <- setNames(lapply(cancers, function(ca) load_cancer_data(ca, data_dir)),
                    cancers)
luad <- caches[["LUAD"]]

## representative ubiquitination-related gene set (E3 / E2 / DUB / UBD)
genes <- c("MUL1", "MDM2", "TRIM44", "NEDD4", "UBE2C", "FBXW7", "USP7", "USP14", "UBB")
ann   <- annotate_genes(luad, genes)
stopifnot(all(ann$ubi_type %in% c("E1", "E2", "E3", "DUB", "UBD", "ULD", "Other")))

res   <- score_genes_multi(caches, genes)
summ  <- merge(res$summary, ann[, c("gene", "ubi_type_full")], by.x = "Gene", by.y = "gene")
long  <- res$scores_long
dims3 <- c(Basic = "Basic_Score", Immune = "Immune_Score", Metabolic = "Metabolic_Score")
facet <- do.call(rbind, lapply(names(dims3), function(nm)
  data.frame(Gene = long$Gene, Cancer = long$Cancer, Dim = nm,
             Score = long[[dims3[[nm]]]], stringsAsFactors = FALSE)))
d4 <- do.call(rbind, lapply(names(dims3), function(nm)
  data.frame(Gene = long$Gene, Dim = nm, Score = long[[dims3[[nm]]]],
             Combined = long$Combined_Score, stringsAsFactors = FALSE)))

## ============ Page 1: gene-set ranking (Shiny tab 1) ============
ggsave(file.path(fig_dir, "fig_geneset_bubble.png"),
       plot_3d_bubble(summ, type_col = "ubi_type_full"), width = 7, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "fig_geneset_heatmap.png"),
       plot_facet_heatmap(facet), width = 9, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, "fig_geneset_ranking.png"),
       plot_ranking_grid(summ, top_n = nrow(summ)), width = 10, height = 7, dpi = 300)
ggsave(file.path(fig_dir, "fig_geneset_alluvial.png"),
       plot_alluvial(long, top_n = nrow(summ)), width = 9, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "fig_geneset_4d.png"),
       plot_4d_dot(d4), width = 7, height = 5.5, dpi = 300)

## ============ Page 2: single-gene query - MUL1 ============
gene <- "MUL1"
gsub <- long[long$Gene == gene, , drop = FALSE]
vals <- c(Ubiquitin = mean(gsub$Ubi_Score, na.rm = TRUE),
          Basic     = mean(gsub$Basic_Score, na.rm = TRUE),
          Immune    = mean(gsub$Immune_Score, na.rm = TRUE),
          Metabolic = mean(gsub$Metabolic_Score, na.rm = TRUE))
png(file.path(fig_dir, "fig_mul1_radar.png"), width = 2400, height = 1800, res = 300)
plot_radar_fmsb(vals, gene = gene)
dev.off()

n_tumor <- function(ca) length(unique(substr(caches[[ca]]$sample_info$sample[
  caches[[ca]]$sample_info$type == "Tumor"], 1, 12)))
imm <- do.call(rbind, lapply(cancers, function(ca) {
  cm <- immune_cor(caches[[ca]], gene, "timer")
  data.frame(Cancer = ca, Trait = colnames(cm), Cor = as.numeric(cm[1, ]),
             P = UbiPanTriage:::cor_pvalue(as.numeric(cm[1, ]), n_tumor(ca)),
             stringsAsFactors = FALSE)
}))
meta <- do.call(rbind, lapply(cancers, function(ca) {
  cm <- metabolic_cor(caches[[ca]], gene)
  data.frame(Cancer = ca, Trait = colnames(cm), Cor = as.numeric(cm[1, ]),
             P = UbiPanTriage:::cor_pvalue(as.numeric(cm[1, ]), n_tumor(ca)),
             stringsAsFactors = FALSE)
}))
ggsave(file.path(fig_dir, "fig_mul1_immune.png"),
       plot_cor_heatmap(imm, title = paste0(gene, " - immune infiltration correlation (TIMER)")),
       width = 7.5, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "fig_mul1_metabolic.png"),
       plot_cor_heatmap(meta, title = paste0(gene, " - metabolic pathway correlation")),
       width = 7.5, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "fig_mul1_expr.png"),
       plot_expr_pancancer(caches, gene, "tpm", show_normal = TRUE),
       width = 9, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "fig_mul1_box.png"),
       plot_expr_box(luad, gene), width = 5, height = 4.5, dpi = 300)
ggsave(file.path(fig_dir, "fig_mul1_roc.png"),
       plot_roc_diag(luad, gene), width = 5, height = 4.5, dpi = 300)
km <- plot_km(luad, gene)
ggsave(file.path(fig_dir, "fig_mul1_km.png"), km$plot, width = 5.5, height = 4.5, dpi = 300)

dup <- gsub[, c("Cancer", "Basic_Score", "Immune_Score", "Metabolic_Score")]
png(file.path(fig_dir, "fig_mul1_upset.png"), width = 2400, height = 1800, res = 300)
print(plot_upset_3dim(dup, threshold = 0.5))
dev.off()

## ---------- summary table for README ----------
tab <- summ[order(-summ$Combined_Score), c("Gene", "ubi_type_full", "Ubi_Score",
                                           "Basic_Score", "Immune_Score",
                                           "Metabolic_Score", "Combined_Score", "Pan_Rank")]
write.csv(tab, file.path(fig_dir, "example_summary.csv"), row.names = FALSE)
print(tab)
cat("\nDONE -", length(list.files(fig_dir, pattern = "\\.png$")), "figures written to", fig_dir, "\n")


