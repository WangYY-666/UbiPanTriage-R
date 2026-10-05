# ============================================================================
# UbiPanTriage Shiny app - pan-cancer ubiquitin / immunometabolic scoring
#   Page 1: ubiquitin gene-set ranking (4 dimensions, adjustable weights)
#   Page 2: ubiquitin single-gene query (radar / correlations / expression /
#           boxplot / ROC / KM / UpSet / inference)
#   Pages 2-4 read the pre-computed per-gene lookup library (inst/extdata/lookup)
#   so any TCGA-expressed gene loads in seconds; fine-grained weight control and
#   on-the-fly re-scoring are provided by the UbiPanTriage R package instead.
#   Page 3: custom gene-set scoring (optional user feature scores)
#   Page 4: custom single-gene query (same plots as page 2)
#   Page 5: usage guide
# Default language: English; switch to Chinese in the top-right corner.
# ============================================================================
options(sass.cache = FALSE)
options(sass.cache_dir = tempdir())

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(ggplot2)
  library(dplyr)
  library(DT)
  library(patchwork)
  library(pROC)
  library(survminer)
  library(shinybusy)
  library(UbiPanTriage)
})

# Long-running servers: bslib compiles the theme CSS into a sub-directory of
# tempdir() and serves it from there, so an OS /tmp cleanup (systemd-tmpfiles /
# tmpwatch / cloud security agents, ~10 days) would break every page with
# "The output directory '/tmp/RtmpXXXX/bslib-XXXX' does not exist".
# The keep-alive timer that prevents this is armed by run_shiny_app(); it must
# NOT be started here: a pending `later` timer that is hours away makes
# shiny::testServer(), which loads this very file, wait for the whole delay.
# The durable fix is TMPDIR pointing outside /tmp (see deploy/Dockerfile).

# ----------------------------------------------------------------------------
# data (pre-computed pan-cancer scores v0.3+)
# ----------------------------------------------------------------------------
data_dir  <- getOption("ubi.data_dir", system.file("extdata", package = "UbiPanTriage"))
message("[plotcache] package pre-rendered PNGs: ",
        length(list.files(file.path(data_dir, "plotcache"), recursive = TRUE, pattern = "[.]png$")))
if (!requireNamespace("png", quietly = TRUE))
  message("[plotcache] WARNING: R package 'png' is not installed; plot grids fall back to live computation")

# compact web lookup library (per-gene files, loaded on demand).
# pages 1-4 read precomputed per-gene data here instead of full cancer caches.
lookup_dir <- file.path(data_dir, "lookup")
if (!dir.exists(lookup_dir))
  stop("web lookup library not found: ", lookup_dir,
       ". Build it with UbiPanTriage/data-raw/build_web_lookup.R first.")

pcs        <- readRDS(file.path(lookup_dir, "pcs_meta.rds"))   # cancer_meta + annotation
scores_tbl <- readRDS(file.path(lookup_dir, "ubi_univ.rds"))   # ubiquitin universe (new scoring)
cancer_meta <- pcs$cancer_meta
ubi_genes  <- sort(unique(scores_tbl$Gene))
gene_index <- readRDS(file.path(lookup_dir, "genes_index.rds"))
lookup_sample_info <- readRDS(file.path(lookup_dir, "sample_info.rds"))
lookup_clinical   <- readRDS(file.path(lookup_dir, "clinical.rds"))
gene_file_map <- tryCatch(readRDS(file.path(lookup_dir, "file_map.rds")),
                          error = function(e) NULL)
web_genes <- sort(unique(gene_index$Gene))
cancers    <- sort(names(lookup_sample_info))
cancer_names <- stats::setNames(cancer_meta$Cancer, cancer_meta$Cancer)   # fallback
if ("cancer_name" %in% colnames(cancer_meta))
  cancer_names <- stats::setNames(cancer_meta$cancer_name, cancer_meta$Cancer)

.dim_cols <- c(ubiquitin = "Ubi_Score", basic = "Basic_Score",
               immune = "Immune_Score", metabolic = "Metabolic_Score")

# default cancer selection for the single-gene pages (fast first load)
.default_cancers <- c("BRCA", "LUAD", "LUSC")

lookup_env <- new.env(parent = emptyenv())
# cap in-memory gene files: 3000 files can exceed 4GB RAM on small servers;
# 800 keeps a single-gene-set query comfortably inside a 4GB container
lookup_env_max <- 800L
get_gene_data <- function(g) {
  if (!exists(g, envir = lookup_env, inherits = FALSE)) {
    if (length(ls(lookup_env, all.names = TRUE)) >= lookup_env_max)
      rm(list = ls(lookup_env, all.names = TRUE), envir = lookup_env)
    fn <- if (!is.null(gene_file_map) && g %in% names(gene_file_map)) gene_file_map[[g]] else g
    assign(g, readRDS(file.path(lookup_dir, "gene", paste0(fn, ".rds"))), envir = lookup_env)
  }
  get(g, envir = lookup_env, inherits = FALSE)
}
# read int16-packed log2(x+1)*100 vector (NA marker -32768) -> integer vector
.read_i16 <- function(r) {
  if (is.null(r)) return(NULL)
  v <- readBin(r, integer(), n = length(r) %/% 2L, size = 2L, signed = TRUE)
  v[v == -32768L] <- NA_integer_
  v
}
# int16 log2(x+1)*100 -> log2(x+1)
.decode16 <- function(r) .read_i16(r) / 100
# tiny per-gene pseudo-cache consumed by the package plot functions
.gene_cache <- function(g, ca) {
  gd <- get_gene_data(g)
  pc <- gd$per_cancer[[ca]]
  if (is.null(pc)) return(NULL)
  si <- lookup_sample_info[[ca]]
  one <- function(v, nm) {
    if (is.null(v)) return(NULL)
    matrix(v, nrow = 1, dimnames = list(nm, si$sample))
  }
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
  list(cancer = ca, sample_info = si, clinical = lookup_clinical[[ca]],
       expr_tpm = one(2^tpm - 1, g),
       expr_tpm_full = one(tpm_i, g),
       expr_count = one(.read_i16(pc$expr$count), g),
       expr_fpkm  = one(.read_i16(pc$expr$fpkm), g),
       pre = list(cor_immune = pre_imm, cor_meta = pre_meta))
}

# ----------------------------------------------------------------------------
# persistent per-cancer plot cache (box / roc / km PNG files)
#   - the package may ship pre-rendered PNGs under inst/extdata/plotcache
#   - runtime misses are computed once and written through to a writable user
#     cache (tools::R_user_dir), so repeat queries just read small PNG files
#     and compose them into a grid - no survminer / pROC / boxplot compute
#   - PNG per (gene x cancer) is ~0.1-0.2 MB, memory-light on 4GB servers
#     (caching ggplot objects as RDS was tested: ~16 MB RAM each, too heavy)
# ----------------------------------------------------------------------------
plot_cache_read <- function(gene, cancer, kind) {
  for (d in .plot_cache_dirs()) {
    f <- file.path(d, kind, paste0(gene, "__", cancer, ".png"))
    if (file.exists(f)) return(f)
  }
  NULL
}
.plot_cache_dirs <- function() {
  c(file.path(data_dir, "plotcache"),
    file.path(tools::R_user_dir("UbiPanTriage", "cache"), "plotcache"))
}
.plot_cache_write <- function(gene, cancer, kind, p) {
  d <- file.path(tools::R_user_dir("UbiPanTriage", "cache"), "plotcache")
  dir.create(file.path(d, kind), recursive = TRUE, showWarnings = FALSE)
  f <- file.path(d, kind, paste0(gene, "__", cancer, ".png"))
  grDevices::png(f, width = 3.9, height = 3.3, units = "in", res = 130)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(p)
  invisible(f)
}
# per-cancer PNG path: from cache, or compute + write-through on a miss
.plot_cached_png <- function(gene, cache, ca, kind, fun) {
  f <- plot_cache_read(gene, ca, kind)
  if (is.null(f)) {
    p <- tryCatch(fun(cache, gene), error = function(e) NULL)
    if (!is.null(p)) f <- .plot_cache_write(gene, ca, kind, p)
  }
  f
}
# first immune method available for a cancer's score row
.imm_method_avail <- function(sc, method) {
  if (any(grepl(paste0("^Imm_", method, "_"), colnames(sc)))) return(method)
  for (m in c("timer", "cibersort", "mcp", "ssgsea", "xcell"))
    if (any(grepl(paste0("^Imm_", m, "_"), colnames(sc)))) return(m)
  method
}
cancers_no_normal <- sort(unique(cancer_meta$Cancer[!cancer_meta$has_normal]))
cancers_no_diff   <- sort(unique(cancer_meta$Cancer[!cancer_meta$has_diff]))
cancers_incomplete <- sort(unique(c(cancers_no_normal, cancers_no_diff)))

# immune trait names per method (from the pre-computed table columns)
imm_traits_for <- function(method) {
  pre <- paste0("Imm_", method, "_")
  cc <- grep(paste0("^", pre), colnames(scores_tbl), value = TRUE)
  sub(pre, "", cc)
}
meta_traits <- function() sub("^Meta_", "", grep("^Meta_", colnames(scores_tbl), value = TRUE))

# ----------------------------------------------------------------------------
# language dictionary (default English)
# ----------------------------------------------------------------------------
L <- list(
  en = list(
    app_title = "UbiPanTriage - Pan-cancer ubiquitin immunometabolic scoring",
    tab1 = "1. Ubiquitin gene-set ranking", tab2 = "2. Ubiquitin single gene",
    tab3 = "3. Custom gene set", tab4 = "4. Custom single gene", tab5 = "5. Usage guide",
    more = "More settings", run = "Run scoring", download = "Download summary (CSV)",
    download_detail = "Download details (Excel)", download_detail_csv = "Download details (CSV)",
    example_dl = "Download example gene file (CSV)",
    reset = "Reset to equal", rel_note = "Sub-item weights are relative (0 - 1, one decimal); they are normalized internally. Set 0 to exclude an item.",
    w_note = "Each dimension weight: 0.1 - 1.0, one decimal place, sum = 1.",
    p1_intro = "Score and rank ubiquitin genes across selected cancers. Four dimensions (ubiquitin / basic / immune / metabolic) are combined with user-adjustable weights; a usage explanation is on page 5.",
    p1_genes = "Ubiquitin genes (multiple)", p1_example = "Example: select all ubiquitin genes",
    p1_example_note = "The example fills in all ubiquitin genes without expanding the box; download the list below.",
    p1_cancers = "Cancers (multiple)", p1_dims = "Score dimensions",
    p1_agg = "Cross-cancer aggregation", agg_equal = "Equal mean (default)",
    agg_n = "Weighted by sample size", p1_wmode = "Dimension weights",
    w_equal = "Equal (1:1:1:1)", w_manual = "Manual",
    p1_imm = "Immune method (TIMER default)",
    p1_basic_w = "Basic sub-item weights (expression / diff / survival / ROC)",
    p1_cell_w = "Immune cell weights (within the selected method)",
    p1_meta_w = "Metabolic pathway weights",
    p1_topn = "Top N genes", p1_thr = "High-score threshold",
    thr_median = "Median (per cancer)", thr_fixed = "Fixed value",
    p1_table = "Gene summary (aggregated across selected cancers)",
    p1_contrib_note = "High: number / list of selected cancers where the gene ranks above the cancer-specific median in that dimension.",
    p1_rank = "Ranking of input genes by dimension (bar charts)",
    p1_4d = "Score bubble plot (x: immune, y: metabolic, size: basic)", p1_facet = "Score heatmap by dimension",
    p1_foot = "Cancers without normal tissue / DEG data are re-normalized over the available basic sub-items (*).",
    p1_missing_cancers = "Cancers lacking normal tissue or DEG data",
    p1_weight_err = "Invalid dimension weights: each must be 0.1 - 1.0 with one decimal place, and the sum must equal 1. Please re-enter.",
    p1_rel_err = "Invalid sub-item weight: values must be between 0 and 1 with one decimal place.",
    p1_no_sel = "Please select at least one gene, one cancer and one dimension.",
    p1_no_data = "No matching data; please check your selection.",
    p1_remind = "Selected cancers lacking normal tissue / DEG data:",
    p1_topnote = "You selected %d genes; only the top %d are shown in the plots.",
    c_gene = "Gene", c_type = "Ubiquitin type", c_type_full = "Ubiquitin class (full name)",
    c_ubi = "Ubi", c_basic = "Basic", c_immune = "Immune", c_meta = "Metabolic",
    c_top_imm = "Top immune cell", c_top_meta = "Top metabolic pathway",
    c_comb = "Combined", c_rank = "Rank", c_pct = "Pct (vs all)", c_nmiss = "n*",
    c_contrib = "High-score cancers", c_missing = "Missing",
    dim_ubiquitin = "Ubiquitin", dim_basic = "Basic", dim_immune = "Immune", dim_metabolic = "Metabolic",
    p2_intro = "Query one ubiquitin gene across cancers: dimension radar, immune / metabolic correlations, expression, tumor-vs-normal boxplot, ROC, survival and an automatic conservative inference.",
    p2_gene = "Gene (single ubiquitin gene)", p2_cancers = "Cancers (multiple)",
    p2_imm = "Immune method (TIMER default)", p2_dims = "Radar dimensions",
    p2_plots = "Plots to show",
    p2_show_radar = "Radar chart", p2_show_imm = "Immune correlation heatmap",
    p2_show_meta = "Metabolic correlation heatmap", p2_show_expr = "Expression comparison",
    p2_show_box = "Tumor vs normal boxplot", p2_show_roc = "Diagnostic ROC",
    p2_show_km = "KM survival curves", p2_show_upset = "UpSet (high basic/immune/metabolic)",
    p2_expr_unit = "Expression unit (TPM default)", p2_expr_normal = "Expression samples",
    expr_tumor = "Tumor only", expr_both = "Tumor + Normal",
    p2_scores = "Scores across selected cancers",
    p2_infer = "Pan-cancer inference",
    p2_non_ubi = "This gene is not in the ubiquitin gene list. Single-gene query only supports ubiquitin-related genes; please use the Custom gene set / Custom single gene pages.",
    p2_progress = "Loading cancer data and computing ...",
    p2_not_found = "The gene was not found in the selected cancers.",
    p2_no_input = "Please type a ubiquitin-related gene and select at least one cancer.",
    p4_no_input = "Please type a gene and select at least one cancer.",
    p2_no_plot = "Not enough data to draw this plot for the selected cancers.",
    dl_no_plot = "No plot available to download yet - please run the query first.",
    p2_imm_note = "Immune scores of cancers without TIMER data fall back to the first available method (see the notes on page 5).",
    p2_upset_note = "UpSet: cancers where the gene passes the high-score rule in basic / immune / metabolic dimensions.",
    p3_intro = "Score any gene set (e.g. phospho-genes, your own list) with the basic / immune / metabolic dimensions; optionally provide your own feature scores (0 - 1) as a fourth dimension.",
    p3_mode_paste = "Paste gene list", p3_mode_file = "Upload file",
    p3_paste = "Genes (one per line), or gene,score pairs when feature scores are used",
    p3_file = "Upload (.txt / .csv / .tsv / .xlsx): first column = genes, optional second column = feature scores (0 - 1)",
    p3_has_score = "My input includes feature scores (second column, 0 - 1)",
    p3_cancers = "Cancers (multiple)", p3_dims = "Score dimensions (feature = your scores)",
    p3_imm = "Immune method", p3_table = "Gene summary (aggregated across selected cancers)",
    p3_missing = "Genes not found in any selected cancer", p3_clip = "Feature scores were clipped to 0 - 1.",
    p3_no_gene = "Please provide at least one gene (paste or upload).",
    p3_renamed = "Auto-matched gene symbols", p3_not_in_lib = "Genes not in the database (check the downloadable gene list)",
    p3_na_note = "\u2014 means no data for that dimension in the selected cancers.",
    dl_genes = "Download all genes (CSV)",
    dl_genes_note = "Full gene symbol list available in the web library (39278 genes).",
    c_cov = "Coverage",
    c_feature = "Feature",
    p4_intro = "Query one custom (non-ubiquitin) gene: optional feature score (0 - 1), dimension radar, immune / metabolic correlations, expression, boxplot, ROC, KM and inference.",
    p4_gene = "Gene (single gene)", p4_feature_opt = "Optional feature score (0 - 1)",
    p4_feature_note = "Leave empty if you have no prior feature score.",
    p4_cancers = "Cancers (multiple)", p4_imm = "Immune method",
    p4_dims = "Radar dimensions", p4_plots = "Plots to show",
    p4_progress = "Loading cancer data and computing ...",
    p4_not_found = "The gene was not found in the selected cancers.",
    p4_gene_renamed = "Gene symbol auto-matched: %s -> %s.",
    p4_gene_unknown = "This gene is not in the web library. Please check the spelling (TCGA may use a different symbol), or download the full gene list to confirm; the R package UbiPanTriage can compute any gene on the fly.",
    p5_title = "Usage guide",
    p5_what = "What does this tool do?",
    p5_what_txt = "UbiPanTriage scores a gene or a gene set in each selected TCGA cancer on four dimensions: (1) Ubiquitin - the functional tier of the gene in the ubiquitination cascade (E3 ligase = 1.00, DUB = 0.90, E2 = 0.80, E1 = 0.75, ULD = 0.70, UBD = 0.65, other = 0.60). This tier ranking is only a rough reference based on the cascade hierarchy; if you prefer your own ranking, supply feature scores in the custom gene-set pages; (2) Basic - evidence from expression level, differential expression (tumor vs normal), univariate Cox survival and diagnostic ROC; (3) Immune - weighted correlations between the gene and immune infiltration traits (TIMER / CIBERSORT / MCPcounter / ssGSEA / xCell); (4) Metabolic - weighted correlations with 7 GSVA metabolic pathways. Within each cancer the basic / immune / metabolic scores are mapped to 0-1 percentile ranks (best = 1) so the three dimensions are comparable across cancers; the ubiquitin tier score keeps its own 0-1 scale.",
    p5_data = "Data source and disclaimer",
    p5_data_txt = "All expression, clinical and annotation data come from The Cancer Genome Atlas (TCGA; https://portal.gdc.cancer.gov/). The interpretation of the raw data belongs to TCGA; this tool only provides a computational analysis method based on those public data and makes no claim over the data themselves.",
    p5_how = "How is the score interpreted?",
    p5_how_txt = "A high basic score suggests the gene is expressed, deregulated and/or prognosis-related; it does NOT by itself guarantee that cell-behavior or tumor-model experiments will succeed - always combine it with the literature, your own experience and clinical data. A high immune score suggests association with the immune microenvironment, pointing the project toward immunology. A high metabolic score suggests metabolic involvement and supports metabolism-oriented studies. The overall score is a weighted combination; all weights are adjustable and the default is equal weighting. Scores and the automatic inference are hypothesis-generating suggestions only - always confirm with your own experiments.",
    p5_weight = "Weight rules",
    p5_weight_txt = "Dimension weights: 0.1 - 1.0, one decimal place, and the sum must equal 1. Sub-item weights (basic sub-items, immune cells, metabolic pathways): 0 - 1 with one decimal place; they are relative and are normalized internally, so set 0 to exclude an item. Within each immune method the cells are equally weighted by default (TIMER is shown by default; other methods can be chosen in More settings). Cross-cancer aggregation is the equal mean by default and can be switched to sample-size weighting.",
    p5_missing = "Missing data",
    p5_missing_txt = "Some cancers have no normal tissue or no DEG results (ACC, DLBC, LAML, LGG, MESO, OV, TGCT, UCS, UVM). For these cancers the basic score is re-normalized over the available sub-items and marked with * in the tables, so the basic score is not systematically lower; still, when comparing across cancers, treat starred values with care. LAML has no TIMER infiltration data; immune scores fall back to ssGSEA.",
    p5_web_r = "Web app vs R package",
    p5_web_r_txt = "The web app is a fast, pre-computed pan-cancer lookup: every gene expressed in TCGA (about 20,000 per cancer) is scored once offline, and queries read the cached results so pages 1-4 load in seconds. All pages share the same pre-computed scores, so numbers are consistent across pages. If a gene symbol is not recognized, download the full gene list (39,278 genes) from pages 3-4 to check the exact symbol (case and version suffixes like .1 are matched automatically). To keep the server responsive, fine-grained controls (sub-item weights, per-cell / per-pathway weights, custom thresholds, re-computation on new data, offline batch analysis) are provided by the companion R package UbiPanTriage. Use the web app for quick exploration; use the R package when you need full control or want to re-score with your own parameters.",
    p5_faq = "FAQ",
    p5_faq_txt = "Q: Why is the immune method default TIMER? A: TIMER is the most widely cited deconvolution method in TCGA studies; the five methods are internally equally weighted, and you can switch the method in More settings. Q: What does a score of 1.000 / 0.123 mean? A: For basic / immune / metabolic, the score is the gene's percentile rank within the scored gene universe of that cancer (best = 1); the ubiquitin score is the tier score. Q: Why is my gene's basic score starred (*)? A: That cancer lacks normal tissue / DEG / survival data, so the basic score was re-normalized over the available sub-items. Q: Can I use non-ubiquitin genes? A: Yes - use pages 3 and 4, which run the same scoring without the ubiquitin dimension. Q: Where are the sub-scores (expression / diff / survival / ROC, immune cells, pathways) and ranks? A: They are all included in the Download details (CSV / Excel) files. The single-gene score tables additionally list the most correlated immune cell type and metabolic pathway as columns.",
    p5_ack = "Acknowledgements",
    p5_ack_txt = "We thank The Cancer Genome Atlas (TCGA) for the public multi-omics data, and the IUCCD2.0 database for the curated ubiquitin gene annotations used in this project. We also thank the open-source R community and the developers of the R packages that this project depends on.",
    p5_note = "Everything shown here is a computational hypothesis; experimental validation is required before drawing conclusions."
  ),
  zh = list(
    app_title = "UbiPanTriage - 泛癌泛素免疫代谢评分系统",
    tab1 = "1. 泛素基因集评分排序", tab2 = "2. 泛素单基因查询",
    tab3 = "3. 自定义基因集", tab4 = "4. 自定义单基因", tab5 = "5. 使用说明",
    more = "更多设置", run = "开始评分", download = "下载汇总 (CSV)",
    download_detail = "下载细分详情 (Excel)", download_detail_csv = "下载细分详情 (CSV)",
    example_dl = "下载示例基因文件 (CSV)",
    reset = "恢复等权", rel_note = "细分权重为相对权重（0-1，保留一位小数），内部自动归一化；设为 0 表示剔除该项目。",
    w_note = "各维度权重需在 0.1-1.0 之间、保留一位小数，且总和为 1。",
    p1_intro = "对所选泛素基因在所选癌种中进行评分与排序。四个维度（泛素 / 基础 / 免疫 / 代谢）权重均可调整；评分标准说明见第 5 页。",
    p1_genes = "泛素基因（可多选）", p1_example = "示例：全选全部泛素基因",
    p1_example_note = "点击示例会在输入框中填入全部泛素基因但不会展开列表；完整基因文件可在下方下载。",
    p1_cancers = "癌种（可多选）", p1_dims = "评分维度",
    p1_agg = "跨癌种汇总方式", agg_equal = "等权平均（默认）",
    agg_n = "按样本量加权", p1_wmode = "维度权重",
    w_equal = "等权（1:1:1:1）", w_manual = "手动设置",
    p1_imm = "免疫浸润方法（默认 TIMER）",
    p1_basic_w = "基础维度细分权重（表达 / 差异 / 生存 / ROC）",
    p1_cell_w = "免疫细胞权重（所选方法内）",
    p1_meta_w = "代谢途径权重",
    p1_topn = "前 N 个基因", p1_thr = "高分阈值",
    thr_median = "中位数（按癌种）", thr_fixed = "固定数值",
    p1_table = "基因汇总（所选癌种汇总）",
    p1_contrib_note = "高分癌种：所选癌种中，该基因在某维度得分高于该癌种全部泛素基因中位数的癌种数量与列表。",
    p1_rank = "输入基因各维度评分排序（条形图）",
    p1_4d = "评分气泡图（横轴：免疫、纵轴：代谢、气泡大小：基础）", p1_facet = "各维度评分热图",
    p1_foot = "缺少正常组织 / 差异数据的癌种，其基础分在可用细分项目上重归一化（*）。",
    p1_missing_cancers = "缺少正常组织或差异数据的癌种",
    p1_weight_err = "维度权重无效：每个需在 0.1-1.0 之间、保留一位小数，且总和必须为 1，请重新输入。",
    p1_rel_err = "细分权重无效：数值需在 0-1 之间并保留一位小数。",
    p1_no_sel = "请至少选择一个基因、一个癌种和一个评分维度。",
    p1_no_data = "未找到匹配数据，请检查您的选择。",
    p1_remind = "所选癌种中缺少正常组织 / 差异数据的癌种：",
    p1_topnote = "共选择了 %d 个基因，绘图仅展示评分前 %d 的基因。",
    c_gene = "基因", c_type = "泛素类型", c_type_full = "泛素类别（全称）",
    c_ubi = "泛素", c_basic = "基础", c_immune = "免疫", c_meta = "代谢",
    c_top_imm = "相关性最高免疫细胞", c_top_meta = "相关性最高代谢通路",
    c_comb = "综合", c_rank = "排名", c_pct = "百分位（vs 全部）", c_nmiss = "n*",
    c_contrib = "高分癌种", c_missing = "缺失",
    dim_ubiquitin = "泛素", dim_basic = "基础", dim_immune = "免疫", dim_metabolic = "代谢",
    p2_intro = "查询单个泛素基因在泛癌中的表现：维度雷达图、免疫 / 代谢相关性、表达比较、肿瘤 vs 正常箱线图、ROC、生存曲线与保守推论。",
    p2_gene = "基因（单个泛素基因）", p2_cancers = "癌种（可多选）",
    p2_imm = "免疫浸润方法（默认 TIMER）", p2_dims = "雷达图维度",
    p2_plots = "显示的图表",
    p2_show_radar = "雷达图", p2_show_imm = "免疫细胞相关性热图",
    p2_show_meta = "代谢途径相关性热图", p2_show_expr = "表达比较",
    p2_show_box = "肿瘤 vs 正常箱线图", p2_show_roc = "诊断 ROC",
    p2_show_km = "KM 生存曲线", p2_show_upset = "UpSet（基础/免疫/代谢高分癌种）",
    p2_expr_unit = "表达数据单位（默认 TPM）", p2_expr_normal = "表达样本",
    expr_tumor = "仅肿瘤", expr_both = "肿瘤 + 正常",
    p2_scores = "所选癌种评分",
    p2_infer = "泛癌推论",
    p2_non_ubi = "该基因不在泛素基因列表中。单基因查询仅支持泛素相关基因，其他基因请使用自定义基因集或自定义单基因页面。",
    p2_progress = "正在载入癌种数据并计算……",
    p2_not_found = "所选癌种中未检出该基因。",
    p2_no_input = "请输入要查询的泛素基因，并至少选择一个癌种。",
    p4_no_input = "请输入基因，并至少选择一个癌种。",
    p2_no_plot = "所选癌种数据不足，无法绘制该图。",
    dl_no_plot = "暂无可下载的图表，请先完成查询。",
    p2_imm_note = "缺少 TIMER 数据的癌种，免疫分自动回退到第一种可用方法（详见第 5 页说明）。",
    p2_upset_note = "UpSet：展示该基因在哪些癌种中同时满足基础 / 免疫 / 代谢的高分规则。",
    p3_intro = "可对任意基因集（如磷酸化相关基因、您关注的基因列表）进行基础 / 免疫 / 代谢三维评分；可选提供您自己的特征打分（0-1）作为第四维。",
    p3_mode_paste = "粘贴基因列表", p3_mode_file = "上传文件",
    p3_paste = "基因列表（每行一个）；若包含特征分，格式为 gene,score",
    p3_file = "上传（.txt / .csv / .tsv / .xlsx）：第一列为基因，可选的第二列为特征分（0-1）",
    p3_has_score = "我的输入包含特征分（第二列，0-1）",
    p3_cancers = "癌种（可多选）", p3_dims = "评分维度（特征 = 您提供的打分）",
    p3_imm = "免疫浸润方法", p3_table = "基因汇总（所选癌种汇总）",
    p3_missing = "在所选癌种中均未检出的基因", p3_clip = "特征分已裁剪至 0-1 范围。",
    p3_no_gene = "请至少提供一个基因（粘贴或上传）。",
    p3_renamed = "已自动匹配的基因名", p3_not_in_lib = "以下基因不在数据库中（请对照下载的基因表核对符号）",
    p3_na_note = "\u2014 表示该维度在所选癌种中无数据。",
    dl_genes = "下载全部基因表（CSV）",
    dl_genes_note = "网页文库收录的全部基因符号（39278 个），查询前可先对照。",
    c_cov = "覆盖度",
    c_feature = "特征",
    p4_intro = "查询任意自定义（非泛素）单基因：可选特征分（0-1）、维度雷达图、免疫 / 代谢相关性、表达比较、箱线图、ROC、KM 与推论。",
    p4_gene = "基因（单个基因）", p4_feature_opt = "可选特征分（0-1）",
    p4_feature_note = "若无先验特征分可留空。",
    p4_cancers = "癌种（可多选）", p4_imm = "免疫浸润方法",
    p4_dims = "雷达图维度", p4_plots = "显示的图表",
    p4_progress = "正在载入癌种数据并计算……",
    p4_not_found = "所选癌种中未检出该基因。",
    p4_gene_renamed = "已自动匹配基因名：%s -> %s。",
    p4_gene_unknown = "该基因不在网页文库中。请检查拼写（TCGA 注释可能使用不同符号），或下载全部基因表核对；R 包 UbiPanTriage 可对任意基因进行实时计算。",
    p5_title = "使用说明",
    p5_what = "这个工具做什么？",
    p5_what_txt = "UbiPanTriage 对所选基因（或基因集）在所选 TCGA 癌种中计算四个维度评分：(1) 泛素维度——基因在泛素化级联中的功能层级（E3 连接酶 = 1.00，DUB = 0.90，E2 = 0.80，E1 = 0.75，ULD = 0.70，UBD = 0.65，其他 = 0.60）。该层级赋分仅是基于级联功能层级的粗略参考；如有自己的排序思路，可在自定义基因集页面自行提供特征分；(2) 基础维度——表达水平、差异表达（肿瘤 vs 正常）、单因素 Cox 生存与诊断 ROC 的证据；(3) 免疫维度——基因与免疫浸润特征（TIMER / CIBERSORT / MCPcounter / ssGSEA / xCell）的加权相关性；(4) 代谢维度——与 7 条 GSVA 代谢途径的加权相关性。在每个癌种内部，基础 / 免疫 / 代谢得分映射为 0-1 百分位（最优为 1），使三个维度在跨癌种时可比；泛素层级分保持自身 0-1 尺度。",
    p5_data = "数据来源与声明",
    p5_data_txt = "所有表达、临床与注释数据均来自癌症基因组图谱（TCGA；https://portal.gdc.cancer.gov/）。原始数据的解释权归 TCGA 所有；本工具仅基于这些公开数据提供一种计算方法，不对数据本身做任何主张。",
    p5_how = "如何解读评分？",
    p5_how_txt = "基础分高提示该基因表达充分、存在差异和/或与预后相关，但并不代表细胞行为学或成瘤实验必然成功，仍需结合文献、个人经验与临床数据综合判断；免疫分高提示与免疫微环境相关，课题可向免疫方向倾斜；代谢分高提示与代谢相关，可开展代谢类科研命题。综合分为各维度的加权组合，默认等权，权重均可调整。评分与自动推论仅为计算假设，具体结论务必以实验为准。",
    p5_weight = "权重规则",
    p5_weight_txt = "维度权重：0.1-1.0、保留一位小数、总和必须为 1。细分权重（基础四子项、免疫细胞、代谢途径）：0-1、保留一位小数，为相对权重，内部自动归一化，设为 0 即剔除该项目。每种免疫方法内部默认等权（默认展示 TIMER，可在更多设置中选择其他方法）。跨癌种汇总默认为等权平均，可切换为按样本量加权。",
    p5_missing = "数据缺失说明",
    p5_missing_txt = "部分癌种没有正常组织或差异分析结果（ACC、DLBC、LAML、LGG、MESO、OV、TGCT、UCS、UVM）。这些癌种的基础分会在可用细分项目上重归一化并以 * 标注，因此基础分不会系统性偏低；但跨癌种比较时对带 * 的值需谨慎。LAML 无 TIMER 浸润数据，免疫分自动回退到 ssGSEA。",
    p5_web_r = "网页版与 R 包的分工",
    p5_web_r_txt = "网页版提供快速、预计算的泛癌查询：TCGA 中每个表达的基因（每癌种约 2 万个）都会离线预先完成全部评分，查询时直接读取缓存结果，因此页面 1-4 数秒内即可加载。页面 1-4 使用同一套预计算评分，数值口径完全一致。若输入的基因提示查不到，可在页面 3/4 下载全部基因表（39278 个）核对符号（大小写、.1 等版本后缀会自动匹配）。为保证服务器响应速度，精细控制（细分权重、细胞 / 通路权重、自定义阈值、用新数据重新计算、离线批量分析）由配套 R 包 UbiPanTriage 提供（install.packages 或 remotes::install_local 安装）。日常快速探索请使用网页版；需要完整控制或按自己的参数重新评分时，请使用 R 包。",
    p5_faq = "常见问题",
    p5_faq_txt = "问：为什么默认用 TIMER？答：TIMER 是 TCGA 研究中最常引用的反卷积方法；五种方法内部等权，可在更多设置中切换。问：评分 1.000 / 0.123 代表什么？答：基础 / 免疫 / 代谢维度是该基因在该癌种内、所在基因集中的百分位（最优为 1）；泛素维度为功能层级分。问：为什么我的基因基础分带 * 号？答：该癌种缺少正常组织 / 差异 / 生存数据，基础分已在可用细分项目上重归一化。问：可以用非泛素基因吗？答：可以，使用第 3、4 页即可，评分流程相同，只是不含泛素维度。问：细分评分在哪里看？答：基础各子项、免疫细胞、代谢途径的细分及各维度排名均在“下载细分详情 (CSV / Excel)”中给出；单基因评分表还会直接列出相关性最高的免疫细胞与代谢通路两列。",
    p5_ack = "致谢",
    p5_ack_txt = "感谢癌症基因组图谱（TCGA）提供的公开多组学数据，以及 IUCCD2.0 泛素基因数据库提供的泛素基因注释。同时感谢开源 R 社区及本项目所依赖的 R 包开发者们。",
    p5_note = "页面所有内容均为计算假设，结论需经实验验证后方可成立。"
  )
)
tr <- function(key, lng) L[[lng]][[key]] %||% key
# ----------------------------------------------------------------------------
# small helpers
# ----------------------------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

.fmt3 <- function(x) ifelse(is.na(x), "\u2014", sprintf("%.3f", as.numeric(x)))

# rbind per-cancer blocks by column-name union (missing columns -> NA)
.bind_rows_fill <- function(blocks) {
  blocks <- Filter(Negate(is.null), blocks)
  if (length(blocks) == 0) return(NULL)
  allc <- unique(unlist(lapply(blocks, colnames)))
  blocks <- lapply(blocks, function(b) {
    miss <- setdiff(allc, colnames(b))
    for (nm in miss) b[[nm]] <- NA
    b[, allc, drop = FALSE]
  })
  do.call(rbind, blocks)
}

# resolve user-entered symbols: exact -> uppercase -> strip version suffix
.normalize_gene <- function(g) {
  if (is.na(g) || !nzchar(g)) return(NA_character_)
  if (g %in% web_genes) return(g)
  u <- toupper(g)
  if (u %in% web_genes) return(u)
  g2 <- sub("\\..*$", "", g)
  if (g2 %in% web_genes) return(g2)
  g3 <- sub("\\..*$", "", u)
  if (g3 %in% web_genes) return(g3)
  NA_character_
}

# survival statistics for a gene across cancers (ubiquitin table or lookup fallback)
.surv_stats <- function(gene, cas) {
  st <- scores_tbl[scores_tbl$Gene == gene & scores_tbl$Cancer %in% cas,
                   c("Cancer", "HR", "HR_p", "direction"), drop = FALSE]
  if (nrow(st) > 0) return(st)
  rows <- lapply(cas, function(ca) {
    pc <- tryCatch(get_gene_data(gene)$per_cancer[[ca]], error = function(e) NULL)
    if (is.null(pc)) return(NULL)
    s <- pc$scores
    data.frame(Cancer = ca, HR = s$HR, HR_p = s$HR_p, direction = s$direction,
               stringsAsFactors = FALSE)
  })
  st2 <- .bind_rows_fill(rows)
  if (is.null(st2))
    st2 <- data.frame(Cancer = character(), HR = numeric(), HR_p = numeric(),
                      direction = character(), stringsAsFactors = FALSE)
  st2
}

.eq_weights <- function(dims) stats::setNames(rep(1 / length(dims), length(dims)), dims)

.manual_default <- function(dims) {
  n <- length(dims)
  if (n == 1) return(stats::setNames(1, dims))
  if (n == 2) return(stats::setNames(c(0.5, 0.5), dims))
  if (n == 3) return(stats::setNames(c(0.4, 0.3, 0.3), dims))
  stats::setNames(rep(0.25, n), dims)
}

.pct_better <- function(x) {
  n <- sum(!is.na(x))
  if (n == 0) return(rep(NA_real_, length(x)))
  (n - rank(x, na.last = "keep", ties.method = "min") + 1) / n
}

.app_na <- function(msg) {
  ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_text(label = msg, size = 4.5, color = "grey40") +
    ggplot2::theme_void() + ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1)
}

# PNG / PDF download buttons shown beside a plot (可视化旁下载按钮)
.plot_dl <- function(id) {
  tags$span(style = "float:right;font-weight:normal;",
            downloadButton(paste0(id, "_png"), "PNG", class = "btn-sm"),
            " ",
            downloadButton(paste0(id, "_pdf"), "PDF", class = "btn-sm"))
}

# weighted mean over available values
.wmean_na <- function(x, w) {
  i <- !is.na(x) & w > 0
  if (!any(i)) return(NA_real_)
  sum(x[i] * w[i]) / sum(w[i])
}

# sub-item weights: 0-1, one decimal (relative, no sum constraint)
.check_rel_weights <- function(w) {
  if (length(w) == 0) return(list(ok = TRUE))
  if (any(is.na(w)) || any(!is.numeric(w))) return(list(ok = FALSE))
  if (any(w < 0 | w > 1)) return(list(ok = FALSE))
  if (any(abs(w - round(w, 1)) > 1e-6)) return(list(ok = FALSE))
  list(ok = TRUE)
}

.parse_genes <- function(txt, has_score = FALSE) {
  if (has_score) {
    lines <- trimws(unlist(strsplit(txt, "[\r\n;]+")))
    lines <- lines[nzchar(lines)]
    if (length(lines) == 0) return(data.frame(Gene = character(0), Score = numeric(0)))
    parts <- strsplit(lines, "[,\t ]+")
    gene <- vapply(parts, `[`, "", 1)
    score <- suppressWarnings(as.numeric(vapply(parts, `[`, "", 2)))
    return(data.frame(Gene = gene, Score = score, stringsAsFactors = FALSE))
  }
  # gene lists may be separated by newline, comma, tab, space or semicolon
  g <- trimws(unlist(strsplit(txt, "[,\t ;\r\n]+")))
  g <- g[nzchar(g)]
  data.frame(Gene = g, Score = NA_real_, stringsAsFactors = FALSE)
}

.read_genes_file <- function(path, ext) {
  if (ext == "xlsx") {
    d <- openxlsx::read.xlsx(path, colNames = FALSE)
  } else {
    sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
    d <- utils::read.table(path, header = FALSE, sep = sep, fill = TRUE,
                           stringsAsFactors = FALSE, quote = "\"", comment.char = "")
  }
  if (ncol(d) == 0) return(data.frame(Gene = character(0), Score = numeric(0)))
  d[] <- lapply(d, as.character)
  gene <- trimws(d[[1]])
  score <- rep(NA_real_, length(gene))
  if (ncol(d) >= 2) {
    s <- suppressWarnings(as.numeric(trimws(d[[2]])))
    if (any(!is.na(s))) score <- s
    hdr <- toupper(gene[1]) %in% c("GENE", "SYMBOL", "FEATURE") && is.na(s[1])
    if (hdr) { gene <- gene[-1]; score <- score[-1] }
  }
  data.frame(Gene = gene[nzchar(gene)], Score = score[nzchar(gene)], stringsAsFactors = FALSE)
}

# ---- recompute dimension scores from the stored sub-scores ----
.recompute_dims <- function(d, dims, basic_w = NULL, imm_method = "timer",
                            cell_w = NULL, meta_w = NULL) {
  if ("basic" %in% dims) {
    bcols <- c("Expression_Score", "Diff_Score", "Survive_Score", "ROC_Score")
    bw <- if (is.null(basic_w)) stats::setNames(rep(1, 4), c("expression", "diff", "survive", "roc")) else basic_w
    wv <- bw[c("expression", "diff", "survive", "roc")]; wv[is.na(wv)] <- 0
    m <- as.matrix(d[, bcols, drop = FALSE])
    wm <- matrix(rep(wv, each = nrow(m)), nrow = nrow(m)); wm[is.na(m)] <- 0
    wt <- rowSums(wm)
    d$Basic_Score <- ifelse(wt > 0, rowSums(m * wm, na.rm = TRUE) / wt, 0)
  } else d$Basic_Score <- NA_real_
  if ("immune" %in% dims) {
    pre <- paste0("Imm_", imm_method, "_")
    icols <- grep(paste0("^", pre), colnames(d), value = TRUE)
    if (length(icols) > 0) {
      m <- abs(as.matrix(d[, icols, drop = FALSE]))
      traits <- sub(pre, "", icols)
      wv <- if (is.null(cell_w)) rep(1, length(traits)) else cell_w[traits]
      wv[is.na(wv)] <- 1
      s <- rowSums(m * rep(wv, each = nrow(m)), na.rm = TRUE) / sum(wv)
      # cancers without the requested method keep their stored (fallback) score
      has <- rowSums(!is.na(m)) > 0
      old <- d$Immune_Score
      d$Immune_Score <- ifelse(has, s, old)
    }
  } else d$Immune_Score <- NA_real_
  if ("metabolic" %in% dims) {
    mcols <- grep("^Meta_", colnames(d), value = TRUE)
    if (length(mcols) > 0) {
      m <- abs(as.matrix(d[, mcols, drop = FALSE]))
      traits <- sub("^Meta_", "", mcols)
      wv <- if (is.null(meta_w)) rep(1, length(traits)) else meta_w[traits]
      wv[is.na(wv)] <- 1
      d$Metabolic_Score <- rowSums(m * rep(wv, each = nrow(m)), na.rm = TRUE) / sum(wv)
    } else d$Metabolic_Score <- NA_real_
  } else d$Metabolic_Score <- NA_real_
  if (!"ubiquitin" %in% dims) d$Ubi_Score <- NA_real_
  d
}


# per-cancer ranks within the universe (higher score = rank 1)
.add_ranks <- function(du) {
  for (col in c("Ubi_Score", "Basic_Score", "Immune_Score", "Metabolic_Score", "Combined_Score")) {
    nm <- sub("_Score$", "", col)
    du[[paste0(nm, "_Rank")]] <- ave(-du[[col]], du$Cancer,
                                     FUN = function(x) rank(x, ties.method = "min", na.last = "keep"))
  }
  du
}

# map basic/immune/metabolic scores to 0-1 percentile ranks within each
# cancer (best = 1, worst = 1/N); the ubiquitin tier score keeps its scale
.pct_map_dims <- function(du) {
  if (!"Cancer" %in% colnames(du) || nrow(du) == 0) return(du)
  for (col in c("Basic_Score", "Immune_Score", "Metabolic_Score"))
    if (col %in% colnames(du))
      du[[col]] <- ave(du[[col]], du$Cancer, FUN = rank_to_pct)
  du
}
# query rows take the universe percentile values (按基因+癌种匹配百分位)
.pct_match <- function(d, du) {
  key <- paste(du$Gene, du$Cancer)
  for (col in c("Basic_Score", "Immune_Score", "Metabolic_Score"))
    d[[col]] <- du[[col]][match(paste(d$Gene, d$Cancer), key)]
  d
}

# gene-level summary from a recomputed long table (query genes d, universe du)
.gene_summary <- function(d, du, agg, cas, dims, genes) {
  wts <- rep(1, length(cas)); names(wts) <- cas
  if (agg == "n") {
    wts <- cancer_meta$n_tumor[match(cas, cancer_meta$Cancer)]
    names(wts) <- cas; wts[is.na(wts)] <- 1
  }
  dim_scores <- c("Ubi_Score", "Basic_Score", "Immune_Score", "Metabolic_Score")
  g <- split(d, d$Gene)
  sm <- do.call(rbind, lapply(g, function(x) {
    data.frame(Gene = x$Gene[1],
               n_cancers = length(unique(x$Cancer)),
               Ubi_Score = .wmean_na(x$Ubi_Score, wts[x$Cancer]),
               Basic_Score = .wmean_na(x$Basic_Score, wts[x$Cancer]),
               Immune_Score = .wmean_na(x$Immune_Score, wts[x$Cancer]),
               Metabolic_Score = .wmean_na(x$Metabolic_Score, wts[x$Cancer]),
               Combined_Score = .wmean_na(x$Combined_Score, wts[x$Cancer]),
               n_missing = sum(!is.na(x$Missing) & nzchar(x$Missing)),
               stringsAsFactors = FALSE)
  }))
  rownames(sm) <- NULL
  # universe combined means for the percentile
  ug <- split(du, du$Gene)
  ucomb <- vapply(ug, function(x) .wmean_na(x$Combined_Score, wts[x$Cancer]), numeric(1))
  sm$Pct <- .pct_better(ucomb[sm$Gene])
  # high-score cancer contribution per dimension
  for (nm in c("Basic", "Immune", "Metabolic", "Ubi")) {
    col <- paste0(nm, "_Score")
    med <- tapply(du[[col]], du$Cancer, median, na.rm = TRUE)
    contrib <- vapply(seq_len(nrow(sm)), function(i) {
      x <- d[d$Gene == sm$Gene[i], , drop = FALSE]
      hi <- cas[!is.na(x[[col]][match(cas, x$Cancer)]) &
                  x[[col]][match(cas, x$Cancer)] > med[cas]]
      if (length(hi) == 0) return("")
      paste0(length(hi), "/", length(cas), ": ", paste(head(hi, 3), collapse = ", "),
             if (length(hi) > 3) ", ..." else "")
    }, character(1))
    sm[[paste0(nm, "_high_cancers")]] <- contrib
  }
  # annotation
  ann <- pcs$annotation
  sm <- merge(sm, ann[, c("Gene", "Ubi_Type", "Ubi_Type_Full")], by = "Gene", all.x = TRUE)
  sm$Rank <- rank(-sm$Combined_Score, ties.method = "min")
  sm[order(sm$Rank), , drop = FALSE]
}

# sticky notice bar html
.notice_html <- function(lng, cas) {
  miss <- intersect(cas, cancers_incomplete)
  if (length(miss) == 0) return("")
  txt <- tr("p1_remind", lng)
  div(class = "notice-bar",
      tags$b(txt), " ", paste(miss, collapse = ", "),
      tags$span(class = "notice-close", onclick = "this.parentNode.style.display='none'", "x"))
}

# format score columns to 3 decimals for exports (导出表格分数保留三位小数)
.fmt_tbl <- function(d) {
  for (cl in grep("Score$", colnames(d), value = TRUE))
    if (is.numeric(d[[cl]])) d[[cl]] <- round(d[[cl]], 3)
  d
}
# detail table for export: drop sample-count columns, add high-cancer lists
.clean_detail <- function(long, summ) {
  drop <- c("n_tumor", "n_normal", "has_diff", "has_normal", "Immune_Method")
  long <- long[, setdiff(colnames(long), drop), drop = FALSE]
  long <- .fmt_tbl(long)
  for (nm in c("Ubi_high_cancers", "Basic_high_cancers",
               "Immune_high_cancers", "Metabolic_high_cancers")) {
    if (nm %in% colnames(summ)) {
      m <- match(long$Gene, summ$Gene)
      long[[nm]] <- summ[[nm]][m]
    }
  }
  long
}
# xlsx detail export (long scores + gene summary)
.write_detail_xlsx <- function(file, long, summ) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Scores by cancer")
  openxlsx::writeData(wb, 1, long)
  openxlsx::addWorksheet(wb, "Gene summary")
  openxlsx::writeData(wb, 2, summ)
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
}

# ---- top correlated immune cell / metabolic pathway per row ----------------
# highest |correlation| within the selected immune method / the 7 pathways
.top_traits <- function(d, method) {
  n <- nrow(d)
  out <- data.frame(Top_immune_cell = rep(NA_character_, n),
                    Top_immune_cor  = rep(NA_real_, n),
                    Top_metabolic_trait = rep(NA_character_, n),
                    Top_metabolic_cor   = rep(NA_real_, n),
                    stringsAsFactors = FALSE)
  pick <- function(m, prefix) {
    if (is.null(m) || ncol(m) == 0)
      return(list(name = rep(NA_character_, n), cor = rep(NA_real_, n)))
    a <- abs(m); a[is.na(a)] <- -Inf
    idx <- max.col(a, ties.method = "first")
    valid <- rowSums(is.finite(a)) > 0
    nm <- sub(prefix, "", colnames(m)[idx])
    cr <- m[cbind(seq_len(n), idx)]
    nm[!valid] <- NA_character_; cr[!valid] <- NA_real_
    list(name = nm, cor = cr)
  }
  ipre <- paste0("Imm_", method, "_")
  icols <- grep(paste0("^", ipre), colnames(d), value = TRUE)
  if (length(icols)) {
    z <- pick(as.matrix(d[, icols, drop = FALSE]), ipre)
    out$Top_immune_cell <- z$name; out$Top_immune_cor <- z$cor
  }
  mcols <- grep("^Meta_", colnames(d), value = TRUE)
  if (length(mcols)) {
    z <- pick(as.matrix(d[, mcols, drop = FALSE]), "^Meta_")
    out$Top_metabolic_trait <- z$name; out$Top_metabolic_cor <- z$cor
  }
  out
}

# ---- detailed export table: every sub-score of every dimension -------------
# one row per gene x cancer with the basic sub-items (+ HR / AUC), every immune
# cell correlation of the selected method, the 7 metabolic pathways, the top
# correlated cell / pathway and the ranks (a gene-set table also carries the
# high-score-cancer lists). Ordered for readability, numeric columns rounded.
.detail_table <- function(long, method, dims) {
  if (is.null(long) || nrow(long) == 0) return(long)
  tt <- .top_traits(long, method)
  out <- data.frame(Gene = long$Gene, Cancer = long$Cancer, stringsAsFactors = FALSE)
  if ("Ubi_Type_Full" %in% colnames(long)) out$Ubi_Type_Full <- long$Ubi_Type_Full
  add <- function(col, name = col)
    if (col %in% colnames(long)) out[[name]] <<- long[[col]]
  if ("ubiquitin" %in% dims) add("Ubi_Score")
  if ("basic" %in% dims)
    for (cl in c("Basic_Score", "Expression_Score", "Diff_Score", "Survive_Score",
                 "ROC_Score", "HR", "HR_p", "AUC_surv", "AUC_diag")) add(cl)
  if ("immune" %in% dims) {
    add("Immune_Score")
    out$Top_immune_cell <- tt$Top_immune_cell
    out$Top_immune_cor  <- tt$Top_immune_cor
    ipre <- paste0("Imm_", method, "_")
    for (cl in grep(paste0("^", ipre), colnames(long), value = TRUE))
      out[[paste0("ImmCor_", sub(ipre, "", cl))]] <- long[[cl]]
  }
  if ("metabolic" %in% dims) {
    add("Metabolic_Score")
    out$Top_metabolic_trait <- tt$Top_metabolic_trait
    out$Top_metabolic_cor   <- tt$Top_metabolic_cor
    for (cl in grep("^Meta_", colnames(long), value = TRUE))
      out[[paste0("MetaCor_", sub("^Meta_", "", cl))]] <- long[[cl]]
  }
  if ("feature" %in% dims) add("Feature_Score")
  for (cl in c("Combined_Score", "Ubi_Rank", "Basic_Rank", "Immune_Rank",
               "Metabolic_Rank", "Combined_Rank", "Missing")) add(cl)
  for (nm in c("Ubi_high_cancers", "Basic_high_cancers",
               "Immune_high_cancers", "Metabolic_high_cancers"))
    if (nm %in% colnames(long)) out[[nm]] <- long[[nm]]
  for (cl in colnames(out)) if (is.numeric(out[[cl]])) out[[cl]] <- round(out[[cl]], 3)
  out
}

# ---- single-gene score summary (display + "download summary" share it) -----
.scores_summary_gene <- function(d, method, include_ubi = TRUE, feature = FALSE) {
  tt <- .top_traits(d, method)
  out <- data.frame(Cancer = d$Cancer, stringsAsFactors = FALSE)
  if (include_ubi) {
    if ("Ubi_Type_Full" %in% colnames(d)) out$Ubi_Type_Full <- d$Ubi_Type_Full
    if ("Ubi_Score" %in% colnames(d)) out$Ubi_Score <- d$Ubi_Score
  }
  out$Basic_Score <- d$Basic_Score
  if ("Missing" %in% colnames(d)) out$Missing <- d$Missing
  out$Immune_Score <- d$Immune_Score
  out$Top_immune_cell <- tt$Top_immune_cell
  out$Metabolic_Score <- d$Metabolic_Score
  out$Top_metabolic_trait <- tt$Top_metabolic_trait
  if (feature && "Feature_Score" %in% colnames(d)) out$Feature_Score <- d$Feature_Score
  out$Combined_Score <- d$Combined_Score
  if ("Combined_Rank" %in% colnames(d)) out$Combined_Rank <- d$Combined_Rank
  for (cl in colnames(out)) if (is.numeric(out[[cl]])) out[[cl]] <- round(out[[cl]], 3)
  out
}

# ----------------------------------------------------------------------------
# UI builders (rebuilt on language switch; values restored from rv$sel)
# ----------------------------------------------------------------------------
p1_ui <- function(lng, sel) {
  fluidRow(
    column(3,
      card(
        card_header(tr("tab1", lng)),
        p(tr("p1_intro", lng), class = "text-muted small"),
        selectizeInput("p1_genes", tr("p1_genes", lng), choices = NULL, multiple = TRUE,
                       selected = sel$genes %||% character(0),
                       options = list(placeholder = "MDM2, TRIM44, ...")),
        actionButton("p1_example_all", tr("p1_example", lng), width = "100%"),
        p(tr("p1_example_note", lng), class = "text-muted small"),
        downloadButton("p1_example_dl", tr("example_dl", lng), width = "100%"),
        br(), br(),
        selectizeInput("p1_cancers", tr("p1_cancers", lng), choices = cancers, multiple = TRUE,
                       selected = sel$cancers %||% cancers),
        accordion(open = FALSE,
          accordion_panel(tr("more", lng),
            checkboxGroupInput("p1_dims", tr("p1_dims", lng), inline = TRUE,
              choices = stats::setNames(
                c("ubiquitin", "basic", "immune", "metabolic"),
                c(tr("dim_ubiquitin", lng), tr("dim_basic", lng),
                  tr("dim_immune", lng), tr("dim_metabolic", lng))),
              selected = sel$dims %||% c("ubiquitin", "basic", "immune", "metabolic")),
            selectInput("p1_imm", tr("p1_imm", lng),
              choices = c("TIMER" = "timer", "CIBERSORT" = "cibersort", "MCPcounter" = "mcp",
                          "ssGSEA" = "ssgsea", "xCell" = "xcell"),
              selected = sel$imm %||% "timer"),
            selectInput("p1_agg", tr("p1_agg", lng),
              choices = stats::setNames(c("equal", "n"), c(tr("agg_equal", lng), tr("agg_n", lng))),
              selected = sel$agg %||% "equal"),
            radioButtons("p1_wmode", tr("p1_wmode", lng), inline = TRUE,
              choices = stats::setNames(c("equal", "manual"), c(tr("w_equal", lng), tr("w_manual", lng))),
              selected = sel$wmode %||% "equal"),
            uiOutput("p1_weights_ui"),
            p(tr("w_note", lng), class = "text-muted small"),
            hr(),
            tags$b(tr("p1_basic_w", lng)),
            fluidRow(
              column(6, numericInput("p1_bw_expression", "Expression", 1, 0, 1, 0.1)),
              column(6, numericInput("p1_bw_diff", "Diff", 1, 0, 1, 0.1)),
              column(6, numericInput("p1_bw_survive", "Survival", 1, 0, 1, 0.1)),
              column(6, numericInput("p1_bw_roc", "ROC", 1, 0, 1, 0.1))
            ),
            hr(),
            tags$b(tr("p1_cell_w", lng)),
            div(style = "margin-bottom:12px;", uiOutput("p1_cell_w_ui")),
            tags$b(tr("p1_meta_w", lng)),
            div(style = "margin-bottom:12px;", uiOutput("p1_meta_w_ui")),
            p(tr("rel_note", lng), class = "text-muted small"),
            hr(),
            numericInput("p1_topn", tr("p1_topn", lng), value = min(sel$topn %||% 20, 20),
                         min = 5, max = 20, step = 5),
            selectInput("p1_thr_mode", tr("p1_thr", lng),
              choices = stats::setNames(c("median", "fixed"),
                c(tr("thr_median", lng), tr("thr_fixed", lng))),
              selected = sel$thr_mode %||% "median"),
            numericInput("p1_thr", NULL, value = sel$thr %||% 0.5, min = 0, max = 1, step = 0.05)
          )
        ),
        br(),
        actionButton("p1_run", tr("run", lng), class = "btn-primary", width = "100%")
      )
    ),
    column(9,
      htmlOutput("p1_notice"),
      card(card_header(tr("p1_table", lng)),
           DTOutput("p1_table"),
           br(),
           downloadButton("p1_download", tr("download", lng)),
           downloadButton("p1_detail_csv", tr("download_detail_csv", lng)),
           downloadButton("p1_detail_dl", tr("download_detail", lng)),
           p(paste0(tr("p1_missing_cancers", lng), ": ",
                    paste(cancers_incomplete, collapse = ", "), ". ", tr("p1_foot", lng)),
             class = "text-muted small"),
           htmlOutput("p1_topnote")),
      card(card_header(tr("p1_4d", lng), .plot_dl("p1_4d")),
           div(class = "plot-wrap", plotOutput("p1_4d", height = "640px"))),
      card(card_header(tr("p1_facet", lng), .plot_dl("p1_facet")),
           div(class = "plot-wrap", plotOutput("p1_facet", height = "720px"))),
      card(card_header(tr("p1_rank", lng), .plot_dl("p1_ranking")),
           div(class = "plot-wrap", plotOutput("p1_ranking", height = "760px")))
    )
  )
}

p2_ui <- function(lng, sel) {
  fluidRow(
    column(3,
      card(
        card_header(tr("tab2", lng)),
        p(tr("p2_intro", lng), class = "text-muted small"),
        selectizeInput("p2_gene", tr("p2_gene", lng), choices = NULL, selected = sel$gene %||% "MDM2"),
        selectizeInput("p2_cancers", tr("p2_cancers", lng), choices = cancers, multiple = TRUE,
                       selected = sel$cancers %||% .default_cancers),
        selectInput("p2_imm", tr("p2_imm", lng),
          choices = c("TIMER" = "timer", "CIBERSORT" = "cibersort", "MCPcounter" = "mcp",
                      "ssGSEA" = "ssgsea", "xCell" = "xcell"),
          selected = sel$imm %||% "timer"),
        accordion(open = FALSE,
          accordion_panel(tr("more", lng),
            checkboxGroupInput("p2_dims", tr("p2_dims", lng), inline = TRUE,
              choices = stats::setNames(
                c("ubiquitin", "basic", "immune", "metabolic"),
                c(tr("dim_ubiquitin", lng), tr("dim_basic", lng),
                  tr("dim_immune", lng), tr("dim_metabolic", lng))),
              selected = sel$dims %||% c("ubiquitin", "basic", "immune", "metabolic")),
            checkboxGroupInput("p2_plots", tr("p2_plots", lng),
              choices = stats::setNames(
                c("radar", "imm", "meta", "expr", "box", "roc", "km", "upset"),
                c(tr("p2_show_radar", lng), tr("p2_show_imm", lng),
                  tr("p2_show_meta", lng), tr("p2_show_expr", lng),
                  tr("p2_show_box", lng), tr("p2_show_roc", lng),
                  tr("p2_show_km", lng), tr("p2_show_upset", lng))),
              selected = sel$plots %||% c("radar", "imm", "meta", "expr", "box", "roc", "km", "upset")),
            selectInput("p2_expr_unit", tr("p2_expr_unit", lng),
              choices = c("TPM" = "tpm", "Count" = "count", "FPKM" = "fpkm"),
              selected = sel$expr_unit %||% "tpm"),
            radioButtons("p2_expr_normal", tr("p2_expr_normal", lng), inline = TRUE,
              choices = stats::setNames(c("both", "tumor"),
                c(tr("expr_both", lng), tr("expr_tumor", lng))),
              selected = sel$expr_normal %||% "both"),
            selectInput("p2_thr_mode", tr("p1_thr", lng),
              choices = stats::setNames(c("median", "fixed"),
                c(tr("thr_median", lng), tr("thr_fixed", lng))),
              selected = sel$thr_mode %||% "median"),
            numericInput("p2_thr", NULL, value = sel$thr %||% 0.5, min = 0, max = 1, step = 0.05)
          )
        ),
        br(),
        actionButton("p2_run", tr("run", lng), class = "btn-primary", width = "100%")
      )
    ),
    column(9,
      htmlOutput("p2_notice"),
      card(card_header(tr("p2_scores", lng)),
           DTOutput("p2_scores"),
           br(),
           downloadButton("p2_download", tr("download", lng)),
           downloadButton("p2_detail_csv", tr("download_detail_csv", lng)),
           downloadButton("p2_detail_dl", tr("download_detail", lng)),
           p(tr("p2_imm_note", lng), class = "text-muted small")),
      card(card_header(tr("p2_show_radar", lng), .plot_dl("p2_radar")),
           div(class = "plot-wrap", plotOutput("p2_radar", height = "540px"))),
      fluidRow(
        column(6, card(card_header(tr("p2_show_imm", lng), .plot_dl("p2_imm_heat")),
                       div(class = "plot-wrap", plotOutput("p2_imm_heat", height = "640px")))),
        column(6, card(card_header(tr("p2_show_meta", lng), .plot_dl("p2_meta_heat")),
                       div(class = "plot-wrap", plotOutput("p2_meta_heat", height = "640px"))))
      ),
      card(card_header(tr("p2_show_expr", lng), .plot_dl("p2_expr")),
           div(class = "plot-wrap", plotOutput("p2_expr", height = "720px"))),
      card(card_header(tr("p2_show_box", lng), .plot_dl("p2_box")),
           div(class = "plot-wrap", imageOutput("p2_box", height = "720px"))),
      card(card_header(tr("p2_show_roc", lng), .plot_dl("p2_roc")),
           div(class = "plot-wrap", imageOutput("p2_roc", height = "720px"))),
      card(card_header(tr("p2_show_km", lng), .plot_dl("p2_km")),
           div(class = "plot-wrap", imageOutput("p2_km", height = "720px"))),
      card(card_header(tr("p2_show_upset", lng), .plot_dl("p2_upset")),
           div(class = "plot-wrap", plotOutput("p2_upset", height = "520px")),
           p(tr("p2_upset_note", lng), class = "text-muted small")),
      card(card_header(tr("p2_infer", lng)),
           pre(id = "p2_infer", class = "inference-box", textOutput("p2_infer_txt")))
    )
  )
}

p3_ui <- function(lng, sel) {
  fluidRow(
    column(3,
      card(
        card_header(tr("tab3", lng)),
        p(tr("p3_intro", lng), class = "text-muted small"),
        radioButtons("p3_mode", NULL, inline = TRUE,
          choices = stats::setNames(c("paste", "file"),
            c(tr("p3_mode_paste", lng), tr("p3_mode_file", lng))),
          selected = sel$mode %||% "paste"),
        conditionalPanel(condition = "input.p3_mode == 'paste'",
          textAreaInput("p3_paste", tr("p3_paste", lng), rows = 8,
                        value = sel$paste %||% "MDM2\nTRIM44\nUSP7")),
        conditionalPanel(condition = "input.p3_mode == 'file'",
          fileInput("p3_file", tr("p3_file", lng), accept = c(".txt", ".csv", ".tsv", ".xlsx"))),
        checkboxInput("p3_has_score", tr("p3_has_score", lng), value = sel$has_score %||% FALSE),
        downloadButton("dl_genes3", tr("dl_genes", lng), class = "btn-sm"),
        p(tr("dl_genes_note", lng), class = "text-muted small"),
        selectizeInput("p3_cancers", tr("p3_cancers", lng), choices = cancers, multiple = TRUE,
                       selected = sel$cancers %||% cancers),
        accordion(open = FALSE,
          accordion_panel(tr("more", lng),
            checkboxGroupInput("p3_dims", tr("p3_dims", lng), inline = TRUE,
              choices = stats::setNames(
                c("basic", "immune", "metabolic", "feature"),
                c(tr("dim_basic", lng), tr("dim_immune", lng),
                  tr("dim_metabolic", lng), tr("c_feature", lng))),
              selected = sel$dims %||% c("basic", "immune", "metabolic")),
            selectInput("p3_imm", tr("p3_imm", lng),
              choices = c("TIMER" = "timer", "CIBERSORT" = "cibersort", "MCPcounter" = "mcp",
                          "ssGSEA" = "ssgsea", "xCell" = "xcell"),
              selected = sel$imm %||% "timer"),
            selectInput("p3_agg", tr("p1_agg", lng),
              choices = stats::setNames(c("equal", "n"), c(tr("agg_equal", lng), tr("agg_n", lng))),
              selected = sel$agg %||% "equal"),
            radioButtons("p3_wmode", tr("p1_wmode", lng), inline = TRUE,
              choices = stats::setNames(c("equal", "manual"), c(tr("w_equal", lng), tr("w_manual", lng))),
              selected = sel$wmode %||% "equal"),
            uiOutput("p3_weights_ui"),
            p(tr("w_note", lng), class = "text-muted small"),
            hr(),
            tags$b(tr("p1_basic_w", lng)),
            fluidRow(
              column(6, numericInput("p3_bw_expression", "Expression", 1, 0, 1, 0.1)),
              column(6, numericInput("p3_bw_diff", "Diff", 1, 0, 1, 0.1)),
              column(6, numericInput("p3_bw_survive", "Survival", 1, 0, 1, 0.1)),
              column(6, numericInput("p3_bw_roc", "ROC", 1, 0, 1, 0.1))
            ),
            hr(),
            tags$b(tr("p1_cell_w", lng)),
            div(style = "margin-bottom:12px;", uiOutput("p3_cell_w_ui")),
            tags$b(tr("p1_meta_w", lng)),
            div(style = "margin-bottom:12px;", uiOutput("p3_meta_w_ui")),
            p(tr("rel_note", lng), class = "text-muted small"),
            hr(),
            numericInput("p3_topn", tr("p1_topn", lng), value = min(sel$topn %||% 20, 20),
                         min = 5, max = 20, step = 5),
            selectInput("p3_thr_mode", tr("p1_thr", lng),
              choices = stats::setNames(c("median", "fixed"),
                c(tr("thr_median", lng), tr("thr_fixed", lng))),
              selected = sel$thr_mode %||% "median"),
            numericInput("p3_thr", NULL, value = sel$thr %||% 0.5, min = 0, max = 1, step = 0.05)
          )
        ),
        br(),
        actionButton("p3_run", tr("run", lng), class = "btn-primary", width = "100%")
      )
    ),
    column(9,
      htmlOutput("p3_notice"),
      card(card_header(tr("p3_table", lng)),
           DTOutput("p3_table"),
           p(tr("p3_na_note", lng), class = "text-muted small"),
           br(),
           downloadButton("p3_download", tr("download", lng)),
           downloadButton("p3_detail_csv", tr("download_detail_csv", lng)),
           downloadButton("p3_detail_dl", tr("download_detail", lng)),
           htmlOutput("p3_topnote")),
      card(card_header(tr("p1_4d", lng), .plot_dl("p3_4d")),
           div(class = "plot-wrap", plotOutput("p3_4d", height = "640px"))),
      card(card_header(tr("p1_facet", lng), .plot_dl("p3_facet")),
           div(class = "plot-wrap", plotOutput("p3_facet", height = "720px"))),
      card(card_header(tr("p1_rank", lng), .plot_dl("p3_ranking")),
           div(class = "plot-wrap", plotOutput("p3_ranking", height = "760px"))),
      htmlOutput("p3_missing_txt")
    )
  )
}

p4_ui <- function(lng, sel) {
  fluidRow(
    column(3,
      card(
        card_header(tr("tab4", lng)),
        p(tr("p4_intro", lng), class = "text-muted small"),
        textInput("p4_gene", tr("p4_gene", lng), value = sel$gene %||% "CDK1"),
        downloadButton("dl_genes4", tr("dl_genes", lng), class = "btn-sm"),
        p(tr("dl_genes_note", lng), class = "text-muted small"),
        numericInput("p4_feature", tr("p4_feature_opt", lng),
                     value = if (is.null(sel$feature)) NA_real_ else sel$feature,
                     min = 0, max = 1, step = 0.05),
        p(tr("p4_feature_note", lng), class = "text-muted small"),
        selectizeInput("p4_cancers", tr("p4_cancers", lng), choices = cancers, multiple = TRUE,
                       selected = sel$cancers %||% .default_cancers),
        selectInput("p4_imm", tr("p4_imm", lng),
          choices = c("TIMER" = "timer", "CIBERSORT" = "cibersort", "MCPcounter" = "mcp",
                      "ssGSEA" = "ssgsea", "xCell" = "xcell"),
          selected = sel$imm %||% "timer"),
        accordion(open = FALSE,
          accordion_panel(tr("more", lng),
            checkboxGroupInput("p4_dims", tr("p4_dims", lng), inline = TRUE,
              choices = stats::setNames(
                c("basic", "immune", "metabolic", "feature"),
                c(tr("dim_basic", lng), tr("dim_immune", lng),
                  tr("dim_metabolic", lng), tr("c_feature", lng))),
              selected = sel$dims %||% c("basic", "immune", "metabolic")),
            checkboxGroupInput("p4_plots", tr("p4_plots", lng),
              choices = stats::setNames(
                c("radar", "imm", "meta", "expr", "box", "roc", "km", "upset"),
                c(tr("p2_show_radar", lng), tr("p2_show_imm", lng),
                  tr("p2_show_meta", lng), tr("p2_show_expr", lng),
                  tr("p2_show_box", lng), tr("p2_show_roc", lng),
                  tr("p2_show_km", lng), tr("p2_show_upset", lng))),
              selected = sel$plots %||% c("radar", "imm", "meta", "expr", "box", "roc", "km", "upset")),
            selectInput("p4_expr_unit", tr("p2_expr_unit", lng),
              choices = c("TPM" = "tpm", "Count" = "count", "FPKM" = "fpkm"),
              selected = sel$expr_unit %||% "tpm"),
            radioButtons("p4_expr_normal", tr("p2_expr_normal", lng), inline = TRUE,
              choices = stats::setNames(c("both", "tumor"),
                c(tr("expr_both", lng), tr("expr_tumor", lng))),
              selected = sel$expr_normal %||% "both"),
            selectInput("p4_thr_mode", tr("p1_thr", lng),
              choices = stats::setNames(c("median", "fixed"),
                c(tr("thr_median", lng), tr("thr_fixed", lng))),
              selected = sel$thr_mode %||% "median"),
            numericInput("p4_thr", NULL, value = sel$thr %||% 0.5, min = 0, max = 1, step = 0.05)
          )
        ),
        br(),
        actionButton("p4_run", tr("run", lng), class = "btn-primary", width = "100%")
      )
    ),
    column(9,
      htmlOutput("p4_notice"),
      card(card_header(tr("p2_scores", lng)),
           DTOutput("p4_scores"),
           br(),
           downloadButton("p4_download", tr("download", lng)),
           downloadButton("p4_detail_csv", tr("download_detail_csv", lng)),
           downloadButton("p4_detail_dl", tr("download_detail", lng)),
           p(tr("p2_imm_note", lng), class = "text-muted small")),
      card(card_header(tr("p2_show_radar", lng), .plot_dl("p4_radar")),
           div(class = "plot-wrap", plotOutput("p4_radar", height = "540px"))),
      fluidRow(
        column(6, card(card_header(tr("p2_show_imm", lng), .plot_dl("p4_imm_heat")),
                       div(class = "plot-wrap", plotOutput("p4_imm_heat", height = "640px")))),
        column(6, card(card_header(tr("p2_show_meta", lng), .plot_dl("p4_meta_heat")),
                       div(class = "plot-wrap", plotOutput("p4_meta_heat", height = "640px"))))
      ),
      card(card_header(tr("p2_show_expr", lng), .plot_dl("p4_expr")),
           div(class = "plot-wrap", plotOutput("p4_expr", height = "720px"))),
      card(card_header(tr("p2_show_box", lng), .plot_dl("p4_box")),
           div(class = "plot-wrap", imageOutput("p4_box", height = "720px"))),
      card(card_header(tr("p2_show_roc", lng), .plot_dl("p4_roc")),
           div(class = "plot-wrap", imageOutput("p4_roc", height = "720px"))),
      card(card_header(tr("p2_show_km", lng), .plot_dl("p4_km")),
           div(class = "plot-wrap", imageOutput("p4_km", height = "720px"))),
      card(card_header(tr("p2_show_upset", lng), .plot_dl("p4_upset")),
           div(class = "plot-wrap", plotOutput("p4_upset", height = "520px")),
           p(tr("p2_upset_note", lng), class = "text-muted small")),
      card(card_header(tr("p2_infer", lng)),
           pre(id = "p4_infer", class = "inference-box", textOutput("p4_infer_txt")))
    )
  )
}

p5_ui <- function(lng) {
  fluidRow(
    column(8, offset = 2,
      card(card_header(tr("p5_title", lng)),
        tags$b(tr("p5_what", lng)), br(),
        p(tr("p5_what_txt", lng)),
        hr(),
        tags$b(tr("p5_data", lng)), br(),
        p(tr("p5_data_txt", lng)),
        hr(),
        tags$b(tr("p5_how", lng)), br(),
        p(tr("p5_how_txt", lng)),
        hr(),
        tags$b(tr("p5_weight", lng)), br(),
        p(tr("p5_weight_txt", lng)),
        hr(),
        tags$b(tr("p5_missing", lng)), br(),
        p(tr("p5_missing_txt", lng)),
        hr(),
        tags$b(tr("p5_web_r", lng)), br(),
        p(tr("p5_web_r_txt", lng)),
        hr(),
        tags$b(tr("p5_faq", lng)), br(),
        p(tr("p5_faq_txt", lng)),
        hr(),
        tags$b(tr("p5_ack", lng)), br(),
        p(tr("p5_ack_txt", lng)),
        hr(),
        p(tr("p5_note", lng), class = "text-muted small"))
    )
  )
}

build_ui <- function(lng, sel) {
  bslib::navset_pill(id = "main_tabs", selected = sel$tab %||% "p1",
    bslib::nav_panel(tr("tab1", lng), value = "p1", p1_ui(lng, sel$p1 %||% list())),
    bslib::nav_panel(tr("tab2", lng), value = "p2", p2_ui(lng, sel$p2 %||% list())),
    bslib::nav_panel(tr("tab3", lng), value = "p3", p3_ui(lng, sel$p3 %||% list())),
    bslib::nav_panel(tr("tab4", lng), value = "p4", p4_ui(lng, sel$p4 %||% list())),
    bslib::nav_panel(tr("tab5", lng), value = "p5", p5_ui(lng))
  )
}

app_theme <- tryCatch({
  bslib::bs_theme_dependencies(bslib::bs_theme(version = 5, bootswatch = "flatly",
                                              primary = "#2C5FA8"))
  bslib::bs_theme(version = 5, bootswatch = "flatly", primary = "#2C5FA8")
}, error = function(e) NULL)

ui <- fluidPage(
  theme = app_theme,
  tags$head(tags$title("UbiPanTriage.com"), tags$style(HTML("
    .topbar { display:flex; align-items:center; justify-content:space-between;
              background:#1F3B57; color:#fff; padding:10px 16px; border-radius:8px; margin-bottom:14px; }
    .topbar .title { font-size:18px; font-weight:600; }
    .topbar .lang { width:150px; margin:0; }
    .inference-box { background:#f6f8fa; border:1px solid #e1e4e8; border-radius:6px;
                     padding:12px; white-space:pre-wrap; font-size:13px; line-height:1.6; }
    .card { border:1px solid #dfe3e8; border-radius:8px; margin-bottom:14px; background:#fff; }
    .card-header { background:#f4f6f8; border-bottom:1px solid #dfe3e8; padding:10px 14px;
                   font-weight:600; border-radius:8px 8px 0 0; }
    .card-body { padding:12px 14px; }
    .notice-bar { position:sticky; top:0; z-index:10; background:#FFF8E1; border:1px solid #F0C36D;
                  border-radius:6px; padding:8px 14px; margin-bottom:12px; font-size:13px;
                  box-shadow:0 2px 6px rgba(0,0,0,.06); }
    .notice-close { float:right; cursor:pointer; color:#B08A3E; font-weight:700; padding-left:8px; }
    .plot-wrap { position:relative; min-height:80px; }
    .plot-wrap .shiny-plot-output { transition:opacity .25s; }
    .plot-wrap .shiny-plot-output[data-shiny-output-busy='true'] { opacity:.25; }
    .plot-wrap::before { content:''; position:absolute; top:50%; left:50%; width:36px; height:36px;
                         margin:-18px 0 0 -18px; border:3px solid #dfe3e8; border-top-color:#2C5FA8;
                         border-radius:50%; animation:ubi-spin .9s linear infinite; display:none; z-index:5; }
    .plot-wrap:has(.shiny-plot-output[data-shiny-output-busy='true'])::before { display:block; }
    @keyframes ubi-spin { to { transform:rotate(360deg); } }
    .nav-pills .nav-link { font-weight:600; }
  "))),
  div(class = "topbar",
      div(class = "title", "UbiPanTriage.com"),
      div(class = "lang", selectInput("lang", NULL, c("English" = "en", "简体中文" = "zh"),
                                      selected = "en", width = "150px"))),
  uiOutput("app_ui"),
  add_busy_bar(color = "#2C5FA8")
)

# ----------------------------------------------------------------------------
# server
# ----------------------------------------------------------------------------
.dim_colmap <- c(ubiquitin = "Ubi_Score", basic = "Basic_Score",
                 immune = "Immune_Score", metabolic = "Metabolic_Score",
                 feature = "Feature_Score")

.combined_scores <- function(d, dims, dim_w) {
  cols <- .dim_colmap[dims]
  cols <- cols[!is.na(cols)]
  if (length(cols) == 0) return(rep(NA_real_, nrow(d)))
  m <- as.matrix(d[, cols, drop = FALSE])
  w <- dim_w[dims]
  wm <- matrix(rep(w, each = nrow(m)), nrow = nrow(m)); wm[is.na(m)] <- 0
  wt <- rowSums(wm)
  ifelse(wt > 0, rowSums(m * wm, na.rm = TRUE) / wt, NA_real_)
}

# aggregated per-gene x dimension long table for the 4D dot plot
.agg_dot <- function(d, dims, agg, cas) {
  wts <- rep(1, length(cas)); names(wts) <- cas
  if (agg == "n") { wts <- cancer_meta$n_tumor[match(cas, cancer_meta$Cancer)]; names(wts) <- cas; wts[is.na(wts)] <- 1 }
  cols <- .dim_colmap[dims]; cols <- cols[!is.na(cols)]
  gs <- split(d, d$Gene)
  rows <- lapply(gs, function(x) {
    data.frame(Gene = x$Gene[1],
               Dim = names(cols),
               Score = vapply(cols, function(cl) .wmean_na(x[[cl]], wts[x$Cancer]), numeric(1)),
               Combined = .wmean_na(x$Combined_Score, wts[x$Cancer]),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# per gene x cancer x dimension long table for the facet heatmap
.facet_long <- function(d, dims) {
  cols <- .dim_colmap[dims]; cols <- cols[!is.na(cols)]
  rows <- lapply(names(cols), function(nm)
    data.frame(Gene = d$Gene, Cancer = d$Cancer, Dim = nm, Score = d[[cols[[nm]]]],
               stringsAsFactors = FALSE))
  do.call(rbind, rows)
}

server <- function(input, output, session) {
  lang <- reactive(input$lang %||% "en")
  rv <- reactiveValues(sel = list())

  # ---- PNG / PDF download support: each renderPlot stores its last output ----
  plot_env <- new.env(parent = emptyenv())
  .set_plot <- function(id, p) assign(id, p, envir = plot_env)
  .get_plot <- function(id) {
    if (exists(id, envir = plot_env, inherits = FALSE)) get(id, envir = plot_env, inherits = FALSE) else NULL
  }
  .reg_plot_dl <- function(id, w = 8, h = 6) {
    reg_one <- function(fmt) {
      output[[paste0(id, "_", fmt)]] <- downloadHandler(
        filename = function() paste0(id, "_", Sys.Date(), ".", fmt),
        content = function(con) {
          p <- tryCatch(.get_plot(id), error = function(e) NULL)
          if (is.null(p)) {
            showNotification(tr("dl_no_plot", lang()), type = "warning", duration = 6); return(NULL)
          }
          # box / roc / km grids are stored as composed PNG file paths
          if (is.character(p)) {
            if (fmt == "png") { file.copy(p, con, overwrite = TRUE); return(NULL) }
            im <- png::readPNG(p)
            asp <- dim(im)[1L] / dim(im)[2L]
            grDevices::pdf(con, width = w, height = h)
            on.exit(grDevices::dev.off(), add = TRUE)
            grid::grid.newpage()
            if (asp > 1) grid::grid.raster(im, width = 0.98 / asp, height = 0.98)
            else grid::grid.raster(im, width = 0.98, height = 0.98 * asp)
            return(NULL)
          }
          if (fmt == "png") grDevices::png(con, width = w, height = h, units = "in", res = 300)
          else grDevices::pdf(con, width = w, height = h)
          on.exit(grDevices::dev.off(), add = TRUE)
          if (is.function(p)) p()
          else if (inherits(p, c("ggplot", "patchwork", "upset"))) print(p)
        })
    }
    for (fmt in c("png", "pdf")) reg_one(fmt)
  }

  # register PNG/PDF download handlers for every plot (buttons active on load;
  # they warn "no plot yet" until the corresponding renderPlot has run)
  .reg_plot_dl("p1_4d", 8, 6)
  .reg_plot_dl("p1_facet", 9, 6)
  .reg_plot_dl("p1_ranking", 9, 7)
  .reg_plot_dl("p2_radar", 7, 6)
  .reg_plot_dl("p2_imm_heat", 8, 7)
  .reg_plot_dl("p2_meta_heat", 8, 6)
  .reg_plot_dl("p2_expr", 10, 6)
  .reg_plot_dl("p2_box", 8, 6)
  .reg_plot_dl("p2_roc", 8, 6)
  .reg_plot_dl("p2_km", 9, 7)
  .reg_plot_dl("p2_upset", 8, 6)
  .reg_plot_dl("p3_4d", 8, 6)
  .reg_plot_dl("p3_facet", 9, 6)
  .reg_plot_dl("p3_ranking", 9, 7)
  .reg_plot_dl("p4_radar", 7, 6)
  .reg_plot_dl("p4_imm_heat", 8, 7)
  .reg_plot_dl("p4_meta_heat", 8, 6)
  .reg_plot_dl("p4_expr", 10, 6)
  .reg_plot_dl("p4_box", 8, 6)
  .reg_plot_dl("p4_roc", 8, 6)
  .reg_plot_dl("p4_km", 9, 7)
  .reg_plot_dl("p4_upset", 8, 6)

  # box / roc / km web grids: compose the per-cancer PNG panels into one image
  # with pure grid (see compose_png_grid) and serve via renderImage. The same
  # composed file powers the PNG/PDF download buttons.
  .grid_img <- function(r, kind, title, fun, id) {
    paths <- if (is.null(r)) character(0) else .grid_png_paths(r, kind, fun)
    out_dir <- file.path(tools::R_user_dir("UbiPanTriage", "cache"), "webplots")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out <- file.path(out_dir, paste0(id, ".png"))
    if (!(length(paths) && compose_png_grid(paths, out)))
      .grid_na_write(out, title)
    .set_plot(id, out)
    list(src = out, contentType = "image/png", deleteFile = FALSE)
  }

  # ---- state persistence (values survive the language-switch UI rebuild) ----
  observe({ rv$sel$tab <- input$main_tabs })
  observe({ if (!is.null(input$p1_genes)) rv$sel$p1$genes <- input$p1_genes })
  observe({ rv$sel$p1$cancers <- input$p1_cancers })
  observe({ rv$sel$p1$dims <- input$p1_dims })
  observe({ rv$sel$p1$imm <- input$p1_imm })
  observe({ rv$sel$p1$agg <- input$p1_agg })
  observe({ rv$sel$p1$wmode <- input$p1_wmode })
  observe({ rv$sel$p1$topn <- input$p1_topn })
  observe({ rv$sel$p1$thr_mode <- input$p1_thr_mode })
  observe({ rv$sel$p1$thr <- input$p1_thr })
  observe({ rv$sel$p2$gene <- input$p2_gene })
  observe({ rv$sel$p2$cancers <- input$p2_cancers })
  observe({ rv$sel$p2$imm <- input$p2_imm })
  observe({ rv$sel$p2$dims <- input$p2_dims })
  observe({ rv$sel$p2$plots <- input$p2_plots })
  observe({ rv$sel$p2$expr_unit <- input$p2_expr_unit })
  observe({ rv$sel$p2$expr_normal <- input$p2_expr_normal })
  observe({ rv$sel$p2$thr_mode <- input$p2_thr_mode })
  observe({ rv$sel$p2$thr <- input$p2_thr })
  observe({ rv$sel$p3$mode <- input$p3_mode })
  observe({ rv$sel$p3$paste <- input$p3_paste })
  observe({ rv$sel$p3$has_score <- input$p3_has_score })
  observe({ rv$sel$p3$cancers <- input$p3_cancers })
  observe({ rv$sel$p3$dims <- input$p3_dims })
  observe({ rv$sel$p3$imm <- input$p3_imm })
  observe({ rv$sel$p3$agg <- input$p3_agg })
  observe({ rv$sel$p3$wmode <- input$p3_wmode })
  observe({ rv$sel$p3$topn <- input$p3_topn })
  observe({ rv$sel$p3$thr_mode <- input$p3_thr_mode })
  observe({ rv$sel$p3$thr <- input$p3_thr })
  observe({ rv$sel$p4$gene <- input$p4_gene })
  observe({ rv$sel$p4$feature <- input$p4_feature })
  observe({ rv$sel$p4$cancers <- input$p4_cancers })
  observe({ rv$sel$p4$imm <- input$p4_imm })
  observe({ rv$sel$p4$dims <- input$p4_dims })
  observe({ rv$sel$p4$plots <- input$p4_plots })
  observe({ rv$sel$p4$expr_unit <- input$p4_expr_unit })
  observe({ rv$sel$p4$expr_normal <- input$p4_expr_normal })
  observe({ rv$sel$p4$thr_mode <- input$p4_thr_mode })
  observe({ rv$sel$p4$thr <- input$p4_thr })

  # ---- weight input renderers ----
  .weight_inputs <- function(prefix, traits, ids) {
    lapply(seq_along(traits), function(i) {
      div(style = "display:inline-block;width:48%;margin-bottom:6px;",
          numericInput(paste0(prefix, "_", ids[i]), traits[i],
                       value = 1, min = 0, max = 1, step = 0.1))
    })
  }
  .trait_id <- function(x) gsub("[^A-Za-z0-9_]", "_", x)

  output$p1_weights_ui <- renderUI({
    dims <- input$p1_dims
    if (length(dims) == 0 || !identical(input$p1_wmode, "manual")) return(NULL)
    w <- .manual_default(dims)
    lapply(dims, function(d) numericInput(paste0("p1_w_", d), tr(paste0("dim_", d), lang()),
                                          value = w[d], min = 0.1, max = 1, step = 0.1))
  })
  output$p3_weights_ui <- renderUI({
    dims <- input$p3_dims
    if (length(dims) == 0 || !identical(input$p3_wmode, "manual")) return(NULL)
    w <- .manual_default(dims)
    lapply(dims, function(d) numericInput(paste0("p3_w_", d), tr(paste0("dim_", d), lang()),
                                          value = w[d], min = 0.1, max = 1, step = 0.1))
  })
  output$p1_cell_w_ui <- renderUI({
    method <- input$p1_imm %||% "timer"
    trs <- imm_traits_for(method)
    if (length(trs) == 0) return(NULL)
    .weight_inputs("p1_cell", trs, vapply(trs, .trait_id, character(1)))
  })
  output$p1_meta_w_ui <- renderUI({
    trs <- meta_traits()
    if (length(trs) == 0) return(NULL)
    .weight_inputs("p1_meta", trs, vapply(trs, .trait_id, character(1)))
  })
  output$p3_cell_w_ui <- renderUI({
    method <- input$p3_imm %||% "timer"
    trs <- imm_traits_for(method)
    if (length(trs) == 0) return(NULL)
    .weight_inputs("p3_cell", trs, vapply(trs, .trait_id, character(1)))
  })
  output$p3_meta_w_ui <- renderUI({
    trs <- meta_traits()
    if (length(trs) == 0) return(NULL)
    .weight_inputs("p3_meta", trs, vapply(trs, .trait_id, character(1)))
  })

  # server-side selectize for the large gene lists
  observe({ updateSelectizeInput(session, "p1_genes", choices = ubi_genes, server = TRUE) })
  observe({ updateSelectizeInput(session, "p2_gene", choices = ubi_genes, server = TRUE) })

  # main UI (rebuilt only when the language changes)
  output$app_ui <- renderUI({
    lng <- lang()
    sel <- isolate(rv$sel)
    build_ui(lng, sel)
  })

  # ---- shared: read and validate weights ----
  get_weights <- function(prefix, dims) {
    if (!identical(input[[paste0(prefix, "_wmode")]], "manual"))
      return(list(ok = TRUE, equal = TRUE, weights = NULL))
    w <- sapply(dims, function(d) input[[paste0(prefix, "_w_", d)]], USE.NAMES = TRUE)
    chk <- check_dim_weights(w, dims)
    if (!chk$ok) {
      showNotification(tr("p1_weight_err", lang()), type = "error", duration = 12)
      return(list(ok = FALSE))
    }
    list(ok = TRUE, equal = FALSE, weights = chk$weights)
  }
  get_basic_weights <- function(prefix) {
    w <- c(expression = input[[paste0(prefix, "_bw_expression")]],
           diff = input[[paste0(prefix, "_bw_diff")]],
           survive = input[[paste0(prefix, "_bw_survive")]],
           roc = input[[paste0(prefix, "_bw_roc")]])
    if (any(is.na(w))) return(stats::setNames(rep(1, 4), c("expression", "diff", "survive", "roc")))
    chk <- .check_rel_weights(w)
    if (!chk$ok) { showNotification(tr("p1_rel_err", lang()), type = "error", duration = 10); return(NULL) }
    w
  }
  get_cell_weights <- function(prefix, method) {
    trs <- imm_traits_for(method)
    if (length(trs) == 0) return(NULL)
    ids <- vapply(trs, .trait_id, character(1))
    vals <- lapply(seq_along(trs), function(i) input[[paste0(prefix, "_cell_", ids[i])]])
    w <- unlist(vals, use.names = FALSE)
    # weight inputs not yet rendered -> default to equal weights (1) (权重输入未就绪时默认等权 1)
    if (length(w) != length(trs) || anyNA(w) || !is.numeric(w))
      return(stats::setNames(rep(1, length(trs)), trs))
    names(w) <- trs
    chk <- .check_rel_weights(w)
    if (!chk$ok) { showNotification(tr("p1_rel_err", lang()), type = "error", duration = 10); return(NULL) }
    w
  }
  get_meta_weights <- function(prefix) {
    trs <- meta_traits()
    if (length(trs) == 0) return(NULL)
    ids <- vapply(trs, .trait_id, character(1))
    vals <- lapply(seq_along(trs), function(i) input[[paste0(prefix, "_meta_", ids[i])]])
    w <- unlist(vals, use.names = FALSE)
    # weight inputs not yet rendered -> default to equal weights (1) (权重输入未就绪时默认等权 1)
    if (length(w) != length(trs) || anyNA(w) || !is.numeric(w))
      return(stats::setNames(rep(1, length(trs)), trs))
    names(w) <- trs
    chk <- .check_rel_weights(w)
    if (!chk$ok) { showNotification(tr("p1_rel_err", lang()), type = "error", duration = 10); return(NULL) }
    w
  }
  # non-blocking reminder when cancers with missing data are selected
  remind_missing <- function(cas) {
    miss <- intersect(cas, cancers_incomplete)
    if (length(miss) > 0)
      showNotification(paste0(tr("p1_remind", lang()), " ", paste(miss, collapse = ", "), " (*)"),
                       type = "message", duration = 6)
  }
  observeEvent(input$p1_cancers, { remind_missing(input$p1_cancers) })
  observeEvent(input$p2_cancers, { remind_missing(input$p2_cancers) })
  observeEvent(input$p3_cancers, { remind_missing(input$p3_cancers) })
  observeEvent(input$p4_cancers, { remind_missing(input$p4_cancers) })

  # ==========================================================================
  # Page 1: ubiquitin gene-set ranking
  # ==========================================================================
  p1_res <- reactiveVal(NULL)

  observeEvent(input$p1_example_all, {
    updateSelectizeInput(session, "p1_genes", selected = ubi_genes)
  })

  output$p1_example_dl <- downloadHandler(
    filename = function() "ubi_example_genes.csv",
    content = function(con) {
      ann <- pcs$annotation
      write.csv(ann[, c("Gene", "Ubi_Type", "Ubi_Type_Full")], con, row.names = FALSE)
    })

  observeEvent(input$p1_run, {
    genes <- input$p1_genes %||% character(0)
    cas   <- input$p1_cancers %||% character(0)
    dims  <- input$p1_dims %||% character(0)
    if (length(genes) == 0 || length(cas) == 0 || length(dims) == 0) {
      showNotification(tr("p1_no_sel", lang()), type = "warning", duration = 10); return()
    }
    if (length(setdiff(genes, ubi_genes)) > 0) {
      showNotification(tr("p2_non_ubi", lang()), type = "warning", duration = 12); return()
    }
    method <- input$p1_imm %||% "timer"
    agg    <- input$p1_agg %||% "equal"
    gw <- get_weights("p1", dims); if (!gw$ok) return()
    weights <- if (gw$equal) .eq_weights(dims) else gw$weights
    bw <- get_basic_weights("p1"); if (is.null(bw)) return()
    cw <- get_cell_weights("p1", method)
    mw <- get_meta_weights("p1")
    if (is.null(cw) && length(imm_traits_for(method)) > 0) return()

    withProgress(message = tr("progress", lang()), value = 0.05, {
      d  <- scores_tbl[scores_tbl$Gene %in% genes & scores_tbl$Cancer %in% cas, , drop = FALSE]
      du <- scores_tbl[scores_tbl$Cancer %in% cas, , drop = FALSE]
      if (nrow(d) == 0) { showNotification(tr("p1_no_data", lang()), type = "warning"); return() }
      incProgress(0.2)
      d  <- .recompute_dims(d, dims, bw, method, cw, mw)
      du <- .recompute_dims(du, dims, bw, method, cw, mw)
      incProgress(0.45)
      du <- .pct_map_dims(du)
      d  <- .pct_match(d, du)
      d$Combined_Score  <- .combined_scores(d, dims, weights)
      du$Combined_Score <- .combined_scores(du, dims, weights)
      incProgress(0.65)
      sm <- .gene_summary(d, du, agg, cas, dims, genes)
      du2 <- .add_ranks(du)
      detail <- du2[du2$Gene %in% genes, , drop = FALSE]
      incProgress(0.95)
      p1_res(list(summ = sm, detail = detail, d = d, du = du,
                  dims = dims, method = method, agg = agg, weights = weights,
                  genes = genes, cas = cas))
    })
  })

  output$p1_topnote <- renderText({
    r <- p1_res(); if (is.null(r)) return("")
    n <- length(r$genes); topn <- min(input$p1_topn %||% 20, 20)
    if (n <= topn) return("")
    sprintf(tr("p1_topnote", lang()), n, topn)
  })
  output$p1_notice <- renderUI({
    r <- p1_res(); if (is.null(r)) return(NULL)
    .notice_html(lang(), r$cas)
  })

  output$p1_table <- renderDT({
    r <- p1_res(); if (is.null(r)) return(NULL)
    lng <- lang(); d <- r$summ; sel <- r$dims
    tab <- data.frame(Gene = d$Gene, stringsAsFactors = FALSE)
    if ("ubiquitin" %in% sel) tab$Type <- d$Ubi_Type_Full
    if ("ubiquitin" %in% sel) tab$Ubi <- .fmt3(d$Ubi_Score)
    if ("basic" %in% sel)     tab$Basic <- .fmt3(d$Basic_Score)
    if ("immune" %in% sel)    tab$Immune <- .fmt3(d$Immune_Score)
    if ("metabolic" %in% sel) tab$Metabolic <- .fmt3(d$Metabolic_Score)
    tab$Combined <- .fmt3(d$Combined_Score)
    tab$Rank <- d$Rank
    cn <- c(tr("c_gene", lng))
    if ("ubiquitin" %in% sel) cn <- c(cn, tr("c_type_full", lng), tr("c_ubi", lng))
    if ("basic" %in% sel)     cn <- c(cn, tr("c_basic", lng))
    if ("immune" %in% sel)    cn <- c(cn, tr("c_immune", lng))
    if ("metabolic" %in% sel) cn <- c(cn, tr("c_meta", lng))
    cn <- c(cn, tr("c_comb", lng), tr("c_rank", lng))
    colnames(tab) <- cn
    datatable(tab, rownames = FALSE, extensions = "Buttons",
              options = list(pageLength = 15, scrollX = TRUE, dom = "Bfrtip",
                             buttons = c("copy", "csv", "excel")))
  })

  output$p1_4d <- renderPlot({ p <- {

    r <- p1_res(); if (is.null(r)) return(NULL)
    plot_3d_bubble(r$summ,
                   type_col = if ("Ubi_Type_Full" %in% colnames(r$summ)) "Ubi_Type_Full" else NULL)
  
}
  .set_plot("p1_4d", p)
  p})

  output$p1_facet <- renderPlot({ p <- {

    r <- p1_res(); if (is.null(r)) return(NULL)
    dims3 <- setdiff(r$dims, "ubiquitin")
    if (length(dims3) == 0) return(.app_na(tr("p2_no_plot", lang())))
    topn <- min(input$p1_topn %||% 20, nrow(r$summ))
    top <- r$summ$Gene[seq_len(topn)]
    fl <- .facet_long(r$d[r$d$Gene %in% top, , drop = FALSE], dims3)
    plot_facet_heatmap(fl)
  
}
  .set_plot("p1_facet", p)
  p})

  output$p1_ranking <- renderPlot({ p <- {

    r <- p1_res(); if (is.null(r)) return(NULL)
    topn <- min(input$p1_topn %||% 20, nrow(r$summ))
    plot_ranking_grid(r$summ, top_n = topn)
  
}
  .set_plot("p1_ranking", p)
  p})

  output$p1_download <- downloadHandler(
    filename = function() paste0("ubi_gene_scores_", Sys.Date(), ".csv"),
    content = function(con) {
      r <- p1_res(); if (is.null(r)) return(NULL)
      write.csv(.fmt_tbl(r$summ), con, row.names = FALSE)
    })
  output$p1_detail_csv <- downloadHandler(
    filename = function() paste0("ubi_gene_detail_", Sys.Date(), ".csv"),
    content = function(con) {
      r <- p1_res(); if (is.null(r)) return(NULL)
      write.csv(.detail_table(r$detail, r$method, r$dims), con, row.names = FALSE)
    })
  output$p1_detail_dl <- downloadHandler(
    filename = function() paste0("ubi_gene_detail_", Sys.Date(), ".xlsx"),
    content = function(con) {
      r <- p1_res(); if (is.null(r)) return(NULL)
      .write_detail_xlsx(con, .detail_table(r$detail, r$method, r$dims), .fmt_tbl(r$summ))
    })

  # ==========================================================================
  # Page 2: ubiquitin single-gene query
  # ==========================================================================
  p2_res <- reactiveVal(NULL)

  # shared per-cancer plot collection (used by pages 2 and 4)
  .collect_plots <- function(cas, gene, method, unit, show_normal, prog = NULL) {
    imm <- meta <- list(); box <- roc <- km <- expr <- caches <- list()
    for (i in seq_along(cas)) {
      if (!is.null(prog)) prog(i / length(cas))
      ca <- cas[i]; cache <- .gene_cache(gene, ca)
      if (is.null(cache)) next
      cm <- tryCatch(immune_cor(cache, gene, method), error = function(e) NULL)
      if (!is.null(cm) && ncol(cm) > 0) {
        n <- length(unique(substr(cache$sample_info$sample[cache$sample_info$type == "Tumor"], 1, 12)))
        imm[[ca]] <- data.frame(Cancer = ca, Trait = colnames(cm), Cor = as.numeric(cm[1, ]),
                                P = UbiPanTriage:::cor_pvalue(as.numeric(cm[1, ]), n),
                                stringsAsFactors = FALSE)
      }
      mt <- tryCatch(metabolic_cor(cache, gene), error = function(e) NULL)
      if (!is.null(mt) && ncol(mt) > 0) {
        n <- length(unique(substr(cache$sample_info$sample[cache$sample_info$type == "Tumor"], 1, 12)))
        meta[[ca]] <- data.frame(Cancer = ca, Trait = colnames(mt), Cor = as.numeric(mt[1, ]),
                                 P = UbiPanTriage:::cor_pvalue(as.numeric(mt[1, ]), n),
                                 stringsAsFactors = FALSE)
      }
      b <- .plot_cached_png(gene, cache, ca, "box", function(cch, g) plot_expr_box(cch, g))
      if (!is.null(b)) box[[ca]] <- b
      rr <- .plot_cached_png(gene, cache, ca, "roc", function(cch, g) plot_roc_diag(cch, g))
      if (!is.null(rr)) roc[[ca]] <- rr
      kk <- .plot_cached_png(gene, cache, ca, "km", function(cch, g) { p <- plot_km(cch, g); p$plot })
      if (!is.null(kk)) km[[ca]] <- kk
      caches[[ca]] <- cache
    }
    list(imm = if (length(imm)) do.call(rbind, imm) else NULL,
         meta = if (length(meta)) do.call(rbind, meta) else NULL,
         box = box, roc = roc, km = km, caches = caches)
  }

  observeEvent(input$p2_run, {
    gene <- input$p2_gene %||% ""
    cas  <- input$p2_cancers %||% character(0)
    if (!nzchar(gene) || length(cas) == 0) {
      showNotification(tr("p2_no_input", lang()), type = "warning", duration = 10); return()
    }
    if (!gene %in% ubi_genes) {
      showNotification(tr("p2_non_ubi", lang()), type = "warning", duration = 12); return()
    }
    method <- input$p2_imm %||% "timer"
    dims   <- input$p2_dims %||% c("ubiquitin", "basic", "immune", "metabolic")
    unit   <- input$p2_expr_unit %||% "tpm"
    show_normal <- identical(input$p2_expr_normal, "both")
    withProgress(message = tr("p2_progress", lang()), value = 0.05, {
      d  <- scores_tbl[scores_tbl$Gene == gene & scores_tbl$Cancer %in% cas, , drop = FALSE]
      du <- scores_tbl[scores_tbl$Cancer %in% cas, , drop = FALSE]
      if (nrow(d) == 0) { showNotification(tr("p2_not_found", lang()), type = "warning"); return() }
      d  <- .recompute_dims(d, dims, imm_method = method)
      du <- .recompute_dims(du, dims, imm_method = method)
      du <- .pct_map_dims(du)
      d  <- .pct_match(d, du)
      d$Combined_Score  <- .combined_scores(d, dims, .eq_weights(dims))
      du$Combined_Score <- .combined_scores(du, dims, .eq_weights(dims))
      incProgress(0.15)
      du2 <- .add_ranks(du)
      rk <- du2[du2$Gene == gene, c("Cancer", "Ubi_Rank", "Basic_Rank", "Immune_Rank",
                                    "Metabolic_Rank", "Combined_Rank"), drop = FALSE]
      d$Ubi_Rank       <- rk$Ubi_Rank[match(d$Cancer, rk$Cancer)]
      d$Basic_Rank     <- rk$Basic_Rank[match(d$Cancer, rk$Cancer)]
      d$Immune_Rank    <- rk$Immune_Rank[match(d$Cancer, rk$Cancer)]
      d$Metabolic_Rank <- rk$Metabolic_Rank[match(d$Cancer, rk$Cancer)]
      d$Combined_Rank  <- rk$Combined_Rank[match(d$Cancer, rk$Cancer)]
      med <- c(
        tapply(du$Basic_Score, du$Cancer, median, na.rm = TRUE),
        tapply(du$Immune_Score, du$Cancer, median, na.rm = TRUE),
        tapply(du$Metabolic_Score, du$Cancer, median, na.rm = TRUE))
      med_nm <- sub("\\..*$", "", names(med))
      med_dim <- c(rep("Basic", length(cas)), rep("Immune", length(cas)), rep("Metabolic", length(cas)))
      med <- stats::setNames(as.numeric(med), paste0(med_nm, ".", med_dim))
      plots <- .collect_plots(cas, gene, method, unit, show_normal,
                              prog = function(f) incProgress(f * 0.8 + 0.2))
      incProgress(1)
      p2_res(list(d = d, du = du, med = med, gene = gene, cas = cas,
                  method = method, dims = dims, unit = unit, show_normal = show_normal,
                  caches = stats::setNames(lapply(cas, function(ca) .gene_cache(gene, ca)), cas),
                  imm = plots$imm, meta = plots$meta, box = plots$box,
                  roc = plots$roc, km = plots$km))
    })
  })

  output$p2_notice <- renderUI({
    r <- p2_res(); if (is.null(r)) return(NULL)
    .notice_html(lang(), r$cas)
  })

  output$p2_scores <- renderDT({
    r <- p2_res(); if (is.null(r)) return(NULL)
    lng <- lang(); s <- .scores_summary_gene(r$d, r$method, include_ubi = TRUE)
    tab <- data.frame(Cancer = s$Cancer, Type = s$Ubi_Type_Full,
                      Ubi = .fmt3(s$Ubi_Score),
                      Basic = ifelse(!is.na(s$Missing) & nzchar(s$Missing),
                                     paste0(.fmt3(s$Basic_Score), "*"), .fmt3(s$Basic_Score)),
                      Immune = .fmt3(s$Immune_Score),
                      TopImm = s$Top_immune_cell,
                      Metabolic = .fmt3(s$Metabolic_Score),
                      TopMeta = s$Top_metabolic_trait,
                      Combined = .fmt3(s$Combined_Score),
                      Rank = s$Combined_Rank,
                      stringsAsFactors = FALSE)
    colnames(tab) <- c("Cancer", tr("c_type_full", lng), tr("c_ubi", lng), tr("c_basic", lng),
                       tr("c_immune", lng), tr("c_top_imm", lng), tr("c_meta", lng),
                       tr("c_top_meta", lng), tr("c_comb", lng), tr("c_rank", lng))
    datatable(tab, rownames = FALSE, extensions = "Buttons",
              options = list(pageLength = 15, scrollX = TRUE, dom = "Bfrtip",
                             buttons = c("copy", "csv", "excel")))
  })

  output$p2_radar <- renderPlot({
    r <- p2_res(); if (is.null(r)) return(NULL)
    if (!"radar" %in% (input$p2_plots %||% character(0))) return(NULL)
    d <- r$d
    vals <- sapply(r$dims, function(x) mean(d[[.dim_cols[x]]], na.rm = TRUE))
    labs <- c(ubiquitin = "Ubiquitin", basic = "Basic", immune = "Immune", metabolic = "Metabolic")
    names(vals) <- labs[names(vals)]
    .set_plot("p2_radar", function() plot_radar_fmsb(vals, gene = r$gene))
    plot_radar_fmsb(vals, gene = r$gene)
  })

  output$p2_imm_heat <- renderPlot({ p <- {

    r <- p2_res(); if (is.null(r)) return(NULL)
    if (!"imm" %in% (input$p2_plots %||% character(0))) return(NULL)
    if (is.null(r$imm)) return(.app_na(tr("p2_no_plot", lang())))
    plot_cor_heatmap(r$imm, title = paste0(r$gene, " - immune infiltration correlation (", r$method, ")"))
  
}
  .set_plot("p2_imm_heat", p)
  p})

  output$p2_meta_heat <- renderPlot({ p <- {

    r <- p2_res(); if (is.null(r)) return(NULL)
    if (!"meta" %in% (input$p2_plots %||% character(0))) return(NULL)
    if (is.null(r$meta)) return(.app_na(tr("p2_no_plot", lang())))
    plot_cor_heatmap(r$meta, title = paste0(r$gene, " - metabolic pathway correlation"))
  
}
  .set_plot("p2_meta_heat", p)
  p})

  output$p2_expr <- renderPlot({ p <- {

    r <- p2_res(); if (is.null(r)) return(NULL)
    if (!"expr" %in% (input$p2_plots %||% character(0))) return(NULL)
    plot_expr_pancancer(r$caches, r$gene, r$unit, r$show_normal)
  
}
  .set_plot("p2_expr", p)
  p}, height = 640)

  output$p2_box <- renderImage({
    r <- p2_res()
    if (is.null(r) || !"box" %in% (input$p2_plots %||% character(0)))
      return(.grid_img(NULL, "box", tr("p2_show_box", lang()), NULL, "p2_box"))
    .grid_img(r, "box", paste0(r$gene, " - tumor vs normal expression"),
              function(cch, g) plot_expr_box(cch, g), "p2_box")
  }, deleteFile = FALSE)

  output$p2_roc <- renderImage({
    r <- p2_res()
    if (is.null(r) || !"roc" %in% (input$p2_plots %||% character(0)))
      return(.grid_img(NULL, "roc", tr("p2_show_roc", lang()), NULL, "p2_roc"))
    .grid_img(r, "roc", paste0(r$gene, " - diagnostic ROC (tumor vs normal)"),
              function(cch, g) plot_roc_diag(cch, g), "p2_roc")
  }, deleteFile = FALSE)

  output$p2_km <- renderImage({
    r <- p2_res()
    if (is.null(r) || !"km" %in% (input$p2_plots %||% character(0)))
      return(.grid_img(NULL, "km", tr("p2_show_km", lang()), NULL, "p2_km"))
    .grid_img(r, "km", paste0(r$gene, " - KM survival curves (median split)"),
              function(cch, g) { p <- plot_km(cch, g); p$plot }, "p2_km")
  }, deleteFile = FALSE)

  output$p2_upset <- renderPlot({ p <- {

    r <- p2_res(); if (is.null(r)) return(NULL)
    if (!"upset" %in% (input$p2_plots %||% character(0))) return(NULL)
    if (identical(input$p2_thr_mode, "fixed"))
      plot_upset_3dim(r$d, threshold = input$p2_thr %||% 0.5)
    else
      plot_upset_3dim(r$d, medians = r$med)
  
}
  .set_plot("p2_upset", p)
  p})

  output$p2_infer_txt <- renderText({
    r <- p2_res(); if (is.null(r)) return(NULL)
    infer_gene(r$gene, r$cas, r$d, lang())
  })

  output$p2_download <- downloadHandler(
    filename = function() paste0("ubi_single_gene_scores_", Sys.Date(), ".csv"),
    content = function(con) {
      r <- p2_res(); if (is.null(r)) return(NULL)
      write.csv(.scores_summary_gene(r$d, r$method, include_ubi = TRUE), con, row.names = FALSE)
    })
  output$p2_detail_csv <- downloadHandler(
    filename = function() paste0("ubi_single_gene_detail_", Sys.Date(), ".csv"),
    content = function(con) {
      r <- p2_res(); if (is.null(r)) return(NULL)
      write.csv(.detail_table(r$d, r$method, r$dims), con, row.names = FALSE)
    })
  output$p2_detail_dl <- downloadHandler(
    filename = function() paste0("ubi_single_gene_detail_", Sys.Date(), ".xlsx"),
    content = function(con) {
      r <- p2_res(); if (is.null(r)) return(NULL)
      .write_detail_xlsx(con, .detail_table(r$d, r$method, r$dims),
                         .scores_summary_gene(r$d, r$method, include_ubi = TRUE))
    })

  # ==========================================================================
  # Page 3: custom gene-set scoring
  # ==========================================================================
  p3_res <- reactiveVal(NULL)

  parse_custom_input <- reactive({
    if (identical(input$p3_mode, "file")) {
      f <- input$p3_file
      if (is.null(f) || !file.exists(f$datapath))
        return(data.frame(Gene = character(0), Score = numeric(0)))
      d <- .read_genes_file(f$datapath, tools::file_ext(f$name))
    } else {
      d <- .parse_genes(input$p3_paste %||% "", isTRUE(input$p3_has_score))
    }
    d <- d[!is.na(d$Gene) & nzchar(d$Gene), , drop = FALSE]
    d$Score <- pmin(1, pmax(0, d$Score))
    d
  })

  observeEvent(input$p3_run, {
    pd <- parse_custom_input()
    if (nrow(pd) == 0) { showNotification(tr("p3_no_gene", lang()), type = "warning"); return() }
    genes_raw <- unique(pd$Gene)
    nm <- vapply(genes_raw, .normalize_gene, character(1))
    renamed <- data.frame(orig = genes_raw[!is.na(nm) & genes_raw != nm],
                          matched = nm[!is.na(nm) & genes_raw != nm],
                          stringsAsFactors = FALSE)
    lib_missing <- genes_raw[is.na(nm)]
    pd$Gene <- nm[match(pd$Gene, genes_raw)]
    pd <- pd[!is.na(pd$Gene), , drop = FALSE]
    if (nrow(pd) == 0) {
      showNotification(tr("p3_no_gene", lang()), type = "warning"); return()
    }
    genes <- unique(pd$Gene)
    feat  <- if (any(!is.na(pd$Score))) stats::setNames(pd$Score, pd$Gene) else NULL
    cas   <- input$p3_cancers %||% character(0)
    dims  <- input$p3_dims %||% character(0)
    if (length(cas) == 0 || length(dims) == 0) {
      showNotification(tr("p1_no_sel", lang()), type = "warning", duration = 10); return()
    }
    if ("feature" %in% dims && is.null(feat)) dims <- setdiff(dims, "feature")
    if (length(dims) == 0) { showNotification(tr("p3_no_gene", lang()), type = "warning"); return() }
    method <- input$p3_imm %||% "timer"
    agg    <- input$p3_agg %||% "equal"
    gw <- get_weights("p3", dims); if (!gw$ok) return()
    weights <- if (gw$equal) .eq_weights(dims) else gw$weights
    bw <- get_basic_weights("p3"); if (is.null(bw)) return()
    cw <- get_cell_weights("p3", method)
    mw <- get_meta_weights("p3")
    if (is.null(cw) && length(imm_traits_for(method)) > 0) return()

    withProgress(message = tr("progress", lang()), value = 0, {
      per <- lapply(seq_along(cas), function(i) {
        incProgress(0.9 / length(cas))
        ca <- cas[i]
        rows <- lapply(genes, function(g) {
          if (!g %in% web_genes) return(NULL)
          pc <- get_gene_data(g)$per_cancer[[ca]]
          if (is.null(pc)) return(NULL)
          sc <- pc$scores; sc$Gene <- g; sc
        })
        rows <- Filter(Negate(is.null), rows)
        if (length(rows) == 0) return(NULL)
        sc <- do.call(rbind, rows)
        imm_method <- .imm_method_avail(sc, method)
        sc <- .recompute_dims(sc, dims, bw, imm_method, cw, mw)
        if ("feature" %in% dims) sc$Feature_Score <- feat[match(sc$Gene, names(feat))]
        sc$Cancer <- ca
        sc <- .pct_map_dims(sc)
        sc$Combined_Score <- .combined_scores(sc, dims, weights)
        sc
      })
      per <- Filter(Negate(is.null), per)
      if (length(per) == 0) { showNotification(tr("p1_no_data", lang()), type = "warning"); return() }
      long <- .bind_rows_fill(per)
      if (!"Missing" %in% colnames(long)) long$Missing <- NA_character_
      long <- .add_ranks(long)
      sm <- .gene_summary(long, long, agg, cas, dims, genes)
      if ("feature" %in% dims)
        sm <- merge(sm, unique(long[, c("Gene", "Feature_Score")]), by = "Gene", all.x = TRUE)
      p3_res(list(summ = sm, detail = long, dims = dims, method = method, agg = agg,
                  weights = weights, genes = genes, cas = cas,
                  missing = setdiff(genes, unique(long$Gene)),
                  lib_missing = lib_missing, renamed = renamed,
                  has_feature = "feature" %in% dims))
    })
  })

  output$p3_topnote <- renderText({
    r <- p3_res(); if (is.null(r)) return("")
    n <- length(r$genes); topn <- min(input$p3_topn %||% 20, 20)
    if (n <= topn) return("")
    sprintf(tr("p1_topnote", lang()), n, topn)
  })
  output$p3_notice <- renderUI({
    r <- p3_res(); if (is.null(r)) return(NULL)
    .notice_html(lang(), r$cas)
  })

  output$p3_missing_txt <- renderUI({
    r <- p3_res(); if (is.null(r)) return(NULL)
    lng <- lang()
    parts <- character()
    if (nrow(r$renamed) > 0)
      parts <- c(parts, paste0(tr("p3_renamed", lng), ": ",
                               paste(sprintf("%s -> %s", r$renamed$orig, r$renamed$matched),
                                     collapse = ", "), "."))
    if (length(r$lib_missing) > 0)
      parts <- c(parts, paste0(tr("p3_not_in_lib", lng), ": ",
                               paste(r$lib_missing, collapse = ", "), "."))
    if (length(r$missing) > 0)
      parts <- c(parts, paste0(tr("p3_missing", lng), ": ",
                               paste(r$missing, collapse = ", "), "."))
    if (length(parts) == 0) return(NULL)
    HTML(paste(parts, collapse = "<br>"))
  })

  output$p3_table <- renderDT({
    r <- p3_res(); if (is.null(r)) return(NULL)
    lng <- lang(); d <- r$summ; sel <- r$dims
    tab <- data.frame(Gene = d$Gene, stringsAsFactors = FALSE)
    if ("basic" %in% sel)     tab$Basic <- .fmt3(d$Basic_Score)
    if ("immune" %in% sel)    tab$Immune <- .fmt3(d$Immune_Score)
    if ("metabolic" %in% sel) tab$Metabolic <- .fmt3(d$Metabolic_Score)
    if ("feature" %in% sel)   tab$Feature <- .fmt3(d$Feature_Score)
    tab$Covered <- paste0(d$n_cancers, "/", length(r$cas))
    tab$Combined <- .fmt3(d$Combined_Score)
    tab$Rank <- d$Rank
    cn <- c(tr("c_gene", lng),
            ifelse("basic" %in% sel, tr("c_basic", lng), ""),
            ifelse("immune" %in% sel, tr("c_immune", lng), ""),
            ifelse("metabolic" %in% sel, tr("c_meta", lng), ""),
            ifelse("feature" %in% sel, tr("c_feature", lng), ""),
            tr("c_cov", lng), tr("c_comb", lng), tr("c_rank", lng))
    colnames(tab) <- cn[nzchar(cn)]
    datatable(tab, rownames = FALSE, extensions = "Buttons",
              options = list(pageLength = 15, scrollX = TRUE, dom = "Bfrtip",
                             buttons = c("copy", "csv", "excel")))
  })

  for (gid in c("dl_genes3", "dl_genes4"))
    output[[gid]] <- downloadHandler(
      filename = function() paste0("ubi_pantriage_genes_", Sys.Date(), ".csv"),
      content = function(con) {
        gi <- gene_index
        gi$n_cancers <- lengths(strsplit(gi$Cancers, ","))
        write.csv(gi[, c("Gene", "n_cancers", "Cancers")], con, row.names = FALSE)
      })

  output$p3_4d <- renderPlot({ p <- {

    r <- p3_res(); if (is.null(r)) return(NULL)
    plot_3d_bubble(r$summ)
  
}
  .set_plot("p3_4d", p)
  p})

  output$p3_facet <- renderPlot({ p <- {

    r <- p3_res(); if (is.null(r)) return(NULL)
    topn <- min(input$p3_topn %||% 20, nrow(r$summ))
    top <- r$summ$Gene[seq_len(topn)]
    fl <- .facet_long(r$detail[r$detail$Gene %in% top, , drop = FALSE], r$dims)
    plot_facet_heatmap(fl)
  
}
  .set_plot("p3_facet", p)
  p})

  output$p3_ranking <- renderPlot({ p <- {

    r <- p3_res(); if (is.null(r)) return(NULL)
    topn <- min(input$p3_topn %||% 20, nrow(r$summ))
    plot_ranking_grid(r$summ, top_n = topn)
  
}
  .set_plot("p3_ranking", p)
  p})

  output$p3_download <- downloadHandler(
    filename = function() paste0("custom_gene_scores_", Sys.Date(), ".csv"),
    content = function(con) {
      r <- p3_res(); if (is.null(r)) return(NULL)
      write.csv(.fmt_tbl(r$summ), con, row.names = FALSE)
    })
  output$p3_detail_csv <- downloadHandler(
    filename = function() paste0("custom_gene_detail_", Sys.Date(), ".csv"),
    content = function(con) {
      r <- p3_res(); if (is.null(r)) return(NULL)
      write.csv(.detail_table(r$detail, r$method, r$dims), con, row.names = FALSE)
    })
  output$p3_detail_dl <- downloadHandler(
    filename = function() paste0("custom_gene_detail_", Sys.Date(), ".xlsx"),
    content = function(con) {
      r <- p3_res(); if (is.null(r)) return(NULL)
      .write_detail_xlsx(con, .detail_table(r$detail, r$method, r$dims), .fmt_tbl(r$summ))
    })

  # ==========================================================================
  # Page 4: custom single-gene query
  # ==========================================================================
  p4_res <- reactiveVal(NULL)

  observeEvent(input$p4_run, {
    gene <- trimws(input$p4_gene %||% "")
    cas  <- input$p4_cancers %||% character(0)
    if (!nzchar(gene) || length(cas) == 0) {
      showNotification(tr("p4_no_input", lang()), type = "warning", duration = 10); return()
    }
    ng <- .normalize_gene(gene)
    if (is.na(ng)) {
      showNotification(tr("p4_gene_unknown", lang()), type = "warning", duration = 12); return()
    }
    if (ng != gene) {
      showNotification(sprintf(tr("p4_gene_renamed", lang()), gene, ng),
                       type = "message", duration = 8)
      gene <- ng
    }
    feat <- input$p4_feature
    if (length(feat) == 0 || is.na(feat)) feat <- NA_real_
    method <- input$p4_imm %||% "timer"
    dims   <- input$p4_dims %||% c("basic", "immune", "metabolic")
    if ("feature" %in% dims && is.na(feat)) dims <- setdiff(dims, "feature")
    if (length(dims) == 0) { showNotification(tr("p1_no_sel", lang()), type = "warning", duration = 10); return() }
    unit <- input$p4_expr_unit %||% "tpm"
    show_normal <- identical(input$p4_expr_normal, "both")
    withProgress(message = tr("p4_progress", lang()), value = 0.05, {
      rows <- lapply(seq_along(cas), function(i) {
        incProgress(0.1 + 0.4 * i / length(cas))
        ca <- cas[i]
        pc <- get_gene_data(gene)$per_cancer[[ca]]
        if (is.null(pc)) return(NULL)
        sc <- pc$scores
        sc$Gene <- gene
        imm_method <- .imm_method_avail(sc, method)
        sc <- .recompute_dims(sc, dims, imm_method = imm_method)
        if ("feature" %in% dims) sc$Feature_Score <- feat
        sc$Combined_Score <- .combined_scores(sc, dims, .eq_weights(dims))
        sc$Cancer <- ca
        sc
      })
      rows <- Filter(Negate(is.null), rows)
      if (length(rows) == 0) { showNotification(tr("p4_not_found", lang()), type = "warning"); return() }
      d <- .bind_rows_fill(rows)
      if (!"Missing" %in% colnames(d)) d$Missing <- NA_character_
      plots <- .collect_plots(cas, gene, method, unit, show_normal,
                              prog = function(f) incProgress(0.5 + f * 0.5))
      incProgress(1)
      p4_res(list(d = d, gene = gene, cas = cas, method = method, dims = dims,
                  unit = unit, show_normal = show_normal, feature = feat,
                  caches = stats::setNames(lapply(cas, function(ca) .gene_cache(gene, ca)), cas),
                  imm = plots$imm, meta = plots$meta, box = plots$box,
                  roc = plots$roc, km = plots$km))
    })
  })

  output$p4_notice <- renderUI({
    r <- p4_res(); if (is.null(r)) return(NULL)
    .notice_html(lang(), r$cas)
  })

  output$p4_scores <- renderDT({
    r <- p4_res(); if (is.null(r)) return(NULL)
    lng <- lang()
    s <- .scores_summary_gene(r$d, r$method, include_ubi = FALSE,
                              feature = "feature" %in% r$dims)
    tab <- data.frame(Cancer = s$Cancer,
                      Basic = ifelse(!is.na(s$Missing) & nzchar(s$Missing),
                                     paste0(.fmt3(s$Basic_Score), "*"), .fmt3(s$Basic_Score)),
                      Immune = .fmt3(s$Immune_Score),
                      TopImm = s$Top_immune_cell,
                      Metabolic = .fmt3(s$Metabolic_Score),
                      TopMeta = s$Top_metabolic_trait,
                      stringsAsFactors = FALSE)
    if ("Feature_Score" %in% colnames(s)) tab$Feature <- .fmt3(s$Feature_Score)
    tab$Combined <- .fmt3(s$Combined_Score)
    cn <- c("Cancer", tr("c_basic", lng), tr("c_immune", lng), tr("c_top_imm", lng),
            tr("c_meta", lng), tr("c_top_meta", lng),
            ifelse("Feature_Score" %in% colnames(s), tr("c_feature", lng), ""),
            tr("c_comb", lng))
    colnames(tab) <- cn[nzchar(cn)]
    datatable(tab, rownames = FALSE, extensions = "Buttons",
              options = list(pageLength = 15, scrollX = TRUE, dom = "Bfrtip",
                             buttons = c("copy", "csv", "excel")))
  })

  output$p4_radar <- renderPlot({
    r <- p4_res(); if (is.null(r)) return(NULL)
    if (!"radar" %in% (input$p4_plots %||% character(0))) return(NULL)
    map <- c(basic = "Basic_Score", immune = "Immune_Score",
             metabolic = "Metabolic_Score", feature = "Feature_Score")
    labs <- c(basic = "Basic", immune = "Immune", metabolic = "Metabolic", feature = "Feature")
    vals <- sapply(r$dims, function(x) mean(r$d[[map[x]]], na.rm = TRUE))
    names(vals) <- labs[names(vals)]
    .set_plot("p4_radar", function() plot_radar_fmsb(vals, gene = r$gene))
    plot_radar_fmsb(vals, gene = r$gene)
  })

  output$p4_imm_heat <- renderPlot({ p <- {

    r <- p4_res(); if (is.null(r)) return(NULL)
    if (!"imm" %in% (input$p4_plots %||% character(0))) return(NULL)
    if (is.null(r$imm)) return(.app_na(tr("p2_no_plot", lang())))
    plot_cor_heatmap(r$imm, title = paste0(r$gene, " - immune infiltration correlation (", r$method, ")"))
  
}
  .set_plot("p4_imm_heat", p)
  p})

  output$p4_meta_heat <- renderPlot({ p <- {

    r <- p4_res(); if (is.null(r)) return(NULL)
    if (!"meta" %in% (input$p4_plots %||% character(0))) return(NULL)
    if (is.null(r$meta)) return(.app_na(tr("p2_no_plot", lang())))
    plot_cor_heatmap(r$meta, title = paste0(r$gene, " - metabolic pathway correlation"))
  
}
  .set_plot("p4_meta_heat", p)
  p})

  output$p4_expr <- renderPlot({ p <- {

    r <- p4_res(); if (is.null(r)) return(NULL)
    if (!"expr" %in% (input$p4_plots %||% character(0))) return(NULL)
    plot_expr_pancancer(r$caches, r$gene, r$unit, r$show_normal)
  
}
  .set_plot("p4_expr", p)
  p}, height = 640)

  output$p4_box <- renderImage({
    r <- p4_res()
    if (is.null(r) || !"box" %in% (input$p4_plots %||% character(0)))
      return(.grid_img(NULL, "box", tr("p2_show_box", lang()), NULL, "p4_box"))
    .grid_img(r, "box", paste0(r$gene, " - tumor vs normal expression"),
              function(cch, g) plot_expr_box(cch, g), "p4_box")
  }, deleteFile = FALSE)

  output$p4_roc <- renderImage({
    r <- p4_res()
    if (is.null(r) || !"roc" %in% (input$p4_plots %||% character(0)))
      return(.grid_img(NULL, "roc", tr("p2_show_roc", lang()), NULL, "p4_roc"))
    .grid_img(r, "roc", paste0(r$gene, " - diagnostic ROC (tumor vs normal)"),
              function(cch, g) plot_roc_diag(cch, g), "p4_roc")
  }, deleteFile = FALSE)

  output$p4_km <- renderImage({
    r <- p4_res()
    if (is.null(r) || !"km" %in% (input$p4_plots %||% character(0)))
      return(.grid_img(NULL, "km", tr("p2_show_km", lang()), NULL, "p4_km"))
    .grid_img(r, "km", paste0(r$gene, " - KM survival curves (median split)"),
              function(cch, g) { p <- plot_km(cch, g); p$plot }, "p4_km")
  }, deleteFile = FALSE)

  output$p4_upset <- renderPlot({ p <- {

    r <- p4_res(); if (is.null(r)) return(NULL)
    if (!"upset" %in% (input$p4_plots %||% character(0))) return(NULL)
    thr <- if (identical(input$p4_thr_mode, "fixed")) (input$p4_thr %||% 0.5) else 0.5
    plot_upset_3dim(r$d, threshold = thr)
  
}
  .set_plot("p4_upset", p)
  p})

  output$p4_infer_txt <- renderText({
    r <- p4_res(); if (is.null(r)) return(NULL)
    infer_gene_custom(r$gene, r$cas, r$d, lang(), r$feature)
  })

  output$p4_download <- downloadHandler(
    filename = function() paste0("custom_single_gene_scores_", Sys.Date(), ".csv"),
    content = function(con) {
      r <- p4_res(); if (is.null(r)) return(NULL)
      write.csv(.scores_summary_gene(r$d, r$method, include_ubi = FALSE,
                                     feature = "feature" %in% r$dims),
                con, row.names = FALSE)
    })
  output$p4_detail_csv <- downloadHandler(
    filename = function() paste0("custom_single_gene_detail_", Sys.Date(), ".csv"),
    content = function(con) {
      r <- p4_res(); if (is.null(r)) return(NULL)
      write.csv(.detail_table(r$d, r$method, r$dims), con, row.names = FALSE)
    })
  output$p4_detail_dl <- downloadHandler(
    filename = function() paste0("custom_single_gene_detail_", Sys.Date(), ".xlsx"),
    content = function(con) {
      r <- p4_res(); if (is.null(r)) return(NULL)
      .write_detail_xlsx(con, .detail_table(r$d, r$method, r$dims),
                         .scores_summary_gene(r$d, r$method, include_ubi = FALSE,
                                              feature = "feature" %in% r$dims))
    })
}

# ----------------------------------------------------------------------------
# shared plot grids + bilingual inference (used by the server)
# ----------------------------------------------------------------------------
grid_plot <- function(plots, title, ncol = 3) {
  if (length(plots) == 0) return(.app_na(title))
  sub <- if (grepl("KM", title))
    "Red = High expression, Blue = Low expression" else NULL
  wrap_plots(plots, ncol = ncol, guides = "collect") +
    plot_annotation(title = title, subtitle = sub,
                    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5)))
}
grid_height <- function(n, ncol = 3) max(420, ceiling(n / ncol) * 320)

# ----------------------------------------------------------------------------
# box / roc / km web grids: the per-cancer panels are PNG files; we compose
# them into ONE image with pure grid (no ggplot), so the web grid renders
# identically on every R / ggplot2 version. (Older ggplot2 on the server
# silently fails to draw rasterGrobs inside annotation_custom, leaving an
# empty plot with only the title - the bug this replaces.)
# ----------------------------------------------------------------------------
compose_png_grid <- function(paths, out_file, ncol = 3L) {
  imgs <- lapply(paths, function(f) tryCatch(png::readPNG(f), error = function(e) NULL))
  keep <- !vapply(imgs, is.null, logical(1))
  if (!all(keep))
    message("[plotcache] unreadable png: ", sum(!keep), "/", length(paths))
  imgs <- imgs[keep]
  if (length(imgs) == 0) return(FALSE)
  hh <- max(vapply(imgs, function(z) dim(z)[1L], integer(1)))
  ww <- max(vapply(imgs, function(z) dim(z)[2L], integer(1)))
  n <- length(imgs); rows <- ceiling(n / ncol)
  canvas <- array(1, dim = c(rows * hh, ncol * ww, 3))
  for (i in seq_len(n)) {
    im <- imgs[[i]]; h <- dim(im)[1L]; w <- dim(im)[2L]
    r <- ceiling(i / ncol); c <- ((i - 1L) %% ncol) + 1L
    y0 <- (r - 1L) * hh + 1L; x0 <- (c - 1L) * ww + 1L
    canvas[y0:(y0 + h - 1L), x0:(x0 + w - 1L), ] <- im
  }
  tryCatch({ png::writePNG(canvas, out_file); TRUE },
           error = function(e) { message("[plotcache] writePNG failed: ", conditionMessage(e)); FALSE })
}

# per-cancer PNG paths for a plot kind; if the query produced none (e.g. an
# older cache without pre-rendered panels), compute them live and write through.
.grid_png_paths <- function(r, kind, fun) {
  paths <- r[[kind]]
  if (length(paths) == 0 && !is.null(r) && length(r$caches) && !is.null(fun)) {
    for (ca in names(r$caches)) {
      f <- .plot_cached_png(r$gene, r$caches[[ca]], ca, kind, fun)
      if (!is.null(f)) paths[[ca]] <- f
    }
    if (length(paths))
      message("[plotcache] live-computed ", length(paths), " ", kind, " panels for ", r$gene)
  }
  paths
}

# placeholder image written to a file (shown when no data / no query yet)
.grid_na_write <- function(out, msg) {
  grDevices::png(out, width = 760, height = 420)
  grid::grid.newpage()
  grid::grid.text(msg, gp = grid::gpar(fontsize = 15, col = "grey40"))
  grDevices::dev.off()
}

infer_gene <- function(gene, cas, d, lng) {
  en <- lng == "en"
  ann <- pcs$annotation[pcs$annotation$Gene == gene, , drop = FALSE]
  lines <- if (en) paste0("[", gene, " pan-cancer inference]") else paste0("[", gene, " 泛癌推论]")
  if (nrow(ann) > 0 && !is.na(ann$Ubi_Type_Full[1]))
    lines <- c(lines, if (en) paste0("Ubiquitin class: ", ann$Ubi_Type_Full[1], ".")
               else paste0("泛素类型：", ann$Ubi_Type_Full[1], "。"))
  pick_top <- function(col, n = 3) {
    dd <- d[!is.na(d[[col]]), , drop = FALSE]
    if (nrow(dd) == 0) return(if (en) "none" else "无")
    dd <- dd[order(dd[[col]], decreasing = TRUE), , drop = FALSE]
    paste(head(dd$Cancer, n), collapse = ", ")
  }
  lines <- c(lines, if (en) paste0("Top cancers by basic score: ", pick_top("Basic_Score"), ".")
             else paste0("基础分最高的癌种：", pick_top("Basic_Score"), "。"))
  lines <- c(lines, if (en) paste0("Top cancers by immune score: ", pick_top("Immune_Score"), ".")
             else paste0("免疫分最高的癌种：", pick_top("Immune_Score"), "。"))
  lines <- c(lines, if (en) paste0("Top cancers by metabolic score: ", pick_top("Metabolic_Score"), ".")
             else paste0("代谢分最高的癌种：", pick_top("Metabolic_Score"), "。"))
  # significant Cox survival directions
  st <- .surv_stats(gene, cas)
  dirs <- character()
  for (i in seq_len(nrow(st))) {
    hr <- st$HR[i]; hp <- st$HR_p[i]; dirn <- st$direction[i]
    if (!is.na(hp) && hp < 0.05 && !is.na(hr) && !is.na(dirn) && hr > 0) {
      role <- if (dirn == "prot") (if (en) "protective" else "抑癌") else (if (en) "risk" else "促癌")
      dirs <- c(dirs, paste0(st$Cancer[i], ": ", role))
    }
  }
  lines <- c(lines, if (length(dirs) > 0)
    (if (en) paste0("Significant survival associations (Cox p < 0.05): ", paste(dirs, collapse = "; "), ".")
     else paste0("显著的生存关联（Cox p < 0.05）：", paste(dirs, collapse = "；"), "。"))
    else (if (en) "No significant Cox survival association in the selected cancers."
          else "所选癌种中未发现显著的 Cox 生存关联。"))
  # conservative recommendation (hypothesis only)
  rec <- d$Cancer[which.max(d$Combined_Score)]
  tips <- character()
  hi <- function(col) any(!is.na(d[[col]]) & d[[col]] > 0.7)
  if (hi("Basic_Score"))
    tips <- c(tips, if (en)
      "high basic score (expression / differential / prognosis evidence is strong); experiment design still needs literature and clinical judgment"
      else "基础分较高（表达、差异、预后证据充分），实验设计仍需结合文献与临床判断")
  if (hi("Immune_Score"))
    tips <- c(tips, if (en) "high immune score favors immune-microenvironment research"
             else "免疫分较高，适合免疫微环境研究")
  if (hi("Metabolic_Score"))
    tips <- c(tips, if (en) "high metabolic score favors metabolic research"
             else "代谢分较高，适合代谢相关研究")
  if (length(tips) == 0) tips <- if (en) "baseline expression should be validated first"
    else "建议先验证基础表达"
  lines <- c(lines, if (en)
    paste0("Recommendation (hypothesis only): validate in ", rec, " first (highest combined score; ",
           paste(tips, collapse = "; "), ").")
    else paste0("建议（仅为推论）：优先在 ", rec, " 中验证（综合分最高；", paste(tips, collapse = "；"), "）。"))
  miss <- intersect(cas, cancers_incomplete)
  if (length(miss) > 0)
    lines <- c(lines, if (en)
      paste0("Note: ", paste(miss, collapse = ", "),
             " lack normal tissue / DEG data; their basic scores are re-normalized (*).")
      else paste0("注：", paste(miss, collapse = "、"),
                  " 缺少正常组织/差异数据，其基础分已重归一化（*）。"))
  paste(lines, collapse = "\n")
}

infer_gene_custom <- function(gene, cas, d, lng, feature = NA_real_) {
  en <- lng == "en"
  lines <- if (en) paste0("[", gene, " pan-cancer inference (custom gene)]")
    else paste0("[", gene, " 泛癌推论（自定义基因）]")
  if (!is.na(feature))
    lines <- c(lines, if (en) paste0("User feature score: ", sprintf("%.3f", feature), ".")
               else paste0("用户特征分：", sprintf("%.3f", feature), "。"))
  pick_top <- function(col, n = 3) {
    dd <- d[!is.na(d[[col]]), , drop = FALSE]
    if (nrow(dd) == 0) return(if (en) "none" else "无")
    dd <- dd[order(dd[[col]], decreasing = TRUE), , drop = FALSE]
    paste(head(dd$Cancer, n), collapse = ", ")
  }
  lines <- c(lines, if (en) paste0("Top cancers by basic score: ", pick_top("Basic_Score"), ".")
             else paste0("基础分最高的癌种：", pick_top("Basic_Score"), "。"))
  lines <- c(lines, if (en) paste0("Top cancers by immune score: ", pick_top("Immune_Score"), ".")
             else paste0("免疫分最高的癌种：", pick_top("Immune_Score"), "。"))
  lines <- c(lines, if (en) paste0("Top cancers by metabolic score: ", pick_top("Metabolic_Score"), ".")
             else paste0("代谢分最高的癌种：", pick_top("Metabolic_Score"), "。"))
  st <- .surv_stats(gene, cas)
  dirs <- character()
  for (i in seq_len(nrow(st))) {
    hr <- st$HR[i]; hp <- st$HR_p[i]; dirn <- st$direction[i]
    if (!is.na(hp) && hp < 0.05 && !is.na(hr) && !is.na(dirn) && hr > 0) {
      role <- if (dirn == "prot") (if (en) "protective" else "抑癌") else (if (en) "risk" else "促癌")
      dirs <- c(dirs, paste0(st$Cancer[i], ": ", role))
    }
  }
  lines <- c(lines, if (length(dirs) > 0)
    (if (en) paste0("Significant survival associations (Cox p < 0.05): ", paste(dirs, collapse = "; "), ".")
     else paste0("显著的生存关联（Cox p < 0.05）：", paste(dirs, collapse = "；"), "。"))
    else (if (en) "No significant Cox survival association in the selected cancers."
          else "所选癌种中未发现显著的 Cox 生存关联。"))
  rec <- d$Cancer[which.max(d$Combined_Score)]
  tips <- character()
  hi <- function(col) any(!is.na(d[[col]]) & d[[col]] > 0.7)
  if (hi("Basic_Score"))
    tips <- c(tips, if (en)
      "high basic score (expression / differential / prognosis evidence is strong); experiment design still needs literature and clinical judgment"
      else "基础分较高（表达、差异、预后证据充分），实验设计仍需结合文献与临床判断")
  if (hi("Immune_Score"))
    tips <- c(tips, if (en) "high immune score favors immune-microenvironment research"
             else "免疫分较高，适合免疫微环境研究")
  if (hi("Metabolic_Score"))
    tips <- c(tips, if (en) "high metabolic score favors metabolic research"
             else "代谢分较高，适合代谢相关研究")
  if (length(tips) == 0) tips <- if (en) "baseline expression should be validated first"
    else "建议先验证基础表达"
  lines <- c(lines, if (en)
    paste0("Recommendation (hypothesis only): validate in ", rec, " first (highest combined score; ",
           paste(tips, collapse = "; "), ").")
    else paste0("建议（仅为推论）：优先在 ", rec, " 中验证（综合分最高；", paste(tips, collapse = "；"), "）。"))
  miss <- intersect(cas, cancers_incomplete)
  if (length(miss) > 0)
    lines <- c(lines, if (en)
      paste0("Note: ", paste(miss, collapse = ", "),
             " lack normal tissue / DEG data; their basic scores are re-normalized (*).")
      else paste0("注：", paste(miss, collapse = "、"),
                  " 缺少正常组织/差异数据，其基础分已重归一化（*）。"))
  paste(lines, collapse = "\n")
}

shinyApp(ui, server)










