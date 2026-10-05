## 0.5.10 (2026-10-06)
- README: append a version query (`?v=20261006`) to the jsDelivr figure URLs so
  GitHub's camo image proxy and client browsers re-fetch the pictures instead of
  reusing a cached (possibly failed) copy. Documentation only.

## 0.5.9 (2026-10-05)
- Fix the Shiny app refusing to start on a correctly built server image with
  "No pre-computed cancer caches found ... extdata". The deployment image
  deliberately deletes the legacy per-cancer `<CANCER>_cache.rds` files (the web
  app reads the per-gene `lookup/` library instead, which keeps the image ~3.8
  GB smaller), but `run_shiny_app()` still required those caches. The launch
  check now accepts the `lookup/` library, or the legacy caches.

## 0.5.8 (2026-10-05)
- README: serve the example figures through the jsDelivr CDN
  (https://cdn.jsdelivr.net/gh/WangYY-666/UbiPanTriage-R@main/man/figures/...)
  instead of repo-relative paths. GitHub renders relative README images from
  raw.githubusercontent.com, which is frequently unreachable in mainland China,
  so the pictures failed to load (broken-image icons); jsDelivr is reachable
  there and GitHub proxies it through camo. No package code changed.

## 0.5.7 (2026-10-05)
- Plots: the shared theme now draws the x / y axis lines and tick marks (the
  light grid background is kept), so every figure clearly shows its axes like a
  classic journal figure. Tile heatmaps (gene x cancer score / correlation) and
  the polar radar chart are excluded, since they have no cartesian axes.
- The change applies to the R package figures and to the Shiny web app; the
  pre-rendered box / ROC / KM PNG cache (inst/extdata/plotcache) and the README
  example gallery were rebuilt with the new theme.

## 0.5.6 (2026-10-05)
- Web app: every result table now offers a "Download details (CSV)" button in
  addition to the Excel export. The detail file lists, for each gene x cancer,
  the basic sub-items (expression / differential / survival / ROC, plus HR and
  AUC), the correlation with every immune cell of the selected method, the 7
  GSVA metabolic pathways, the top correlated immune cell / metabolic pathway,
  the per-dimension ranks and (for gene sets) the high-score-cancer lists. The
  plain "Download summary (CSV)" keeps the concise table shown on screen.
- Web app single-gene tables (ubiquitin single gene, page 2; custom single gene,
  page 4) gained two columns: the most correlated immune cell type and the most
  correlated metabolic pathway, derived from the selected immune method.
- Pages 3 and 4 (custom gene set / custom single gene) follow the same detail
  download and extra-column changes.

## 0.5.5 (2026-09-24)
- Fix "The output directory '/tmp/RtmpXXXX/bslib-XXXX' does not exist" on
  long-running web servers: bslib compiles the theme CSS into a sub-directory
  of tempdir() and serves it from there, and the OS eventually deletes
  stale /tmp entries (systemd-tmpfiles / tmpwatch / cloud security agents),
  after which every page request fails. run_shiny_app() now refreshes the
  modification time of the session temp files every 6 hours (`later`), so
  age-based /tmp cleanup no longer picks them up. Applies to both the Shiny
  app shipped in inst/shiny and run_shiny_app().
- The bundled deployment guide/Dockerfile now set TMPDIR=/app/tmp (a path
  outside /tmp that nothing cleans) and install a daily container restart,
  so a self-hosted server stays healthy without manual maintenance.

## 0.5.4 (2026-09-01)
- download_extdata() tries mirror sources automatically when the primary
  raw.githubusercontent.com source fails (e.g. DNS blocking in mainland
  China): user-supplied `mirror` bases first, then built-in gh-proxy
  mirrors. Add `download_extdata(mirror = "https://...")` to force one.
- Failure message now lists every tried URL and the last error, with clear
  guidance (VPN/proxy, mirror, or local_dir offline install).

## 0.5.3 (2026-08-31)
- Re-split the pan-cancer data into 68 parts of 64 MiB (67.1 MB) each
  (`inst/extdata_manifest.csv` updated). The previous 46 parts were 100.66 MB
  each, above GitHub's 100 MB per-file limit, so the `data` branch push was
  rejected by GitHub (HTTP 500) and `download_extdata()` returned 404 for all
  users. The new parts reconstruct the byte-identical archive (verified MD5).

## 0.5.2 (2026-08-31)
- Publish the complete 46-part pan-cancer data on the GitHub `data` branch so
  `download_extdata()` works for all users (previously the branch was not
  pushed and downloads failed with HTTP 404)
- `download_extdata()` gains a `release` argument: download the data parts
  from a GitHub Release (e.g. `download_extdata(release = "data-v1")`) as an
  alternative to the `data` branch; ships the parts manifest
  `inst/extdata_release_manifest.csv`
- More robust parts cache: fall back to the session temp directory when the R
  user cache dir is not writable, so downloads no longer fail on locked-down
  cache locations
- README / FAQ updated with network guidance for mainland China

## 0.5.1 (2026-08-30)
- Remove plot_alluvial() (gene -> cancer -> direction alluvial diagram) and its
  README section; ggalluvial dropped from Suggests
- score_genes() / score_genes_multi() gain a feature_scores argument: a named
  numeric vector (0-1) that replaces the annotation-tier ubiquitin score with
  user-supplied feature scores for custom gene sets
- download_extdata() gains local_dir: install the pre-computed caches from an
  existing local extdata folder (offline install) and prints a clear hint when
  the data-branch download returns HTTP 404
- README: add Example 3 (custom gene set with/without feature scores) and
  Example 4 (custom single gene)
## 0.5.0 (2026-08-30)
- README fully rewritten in English (SCI style): web portal link
  (http://129.211.3.138/) for quick checks; the R package positioned as the
  full-autonomy version with method selection, weight tuning and downloadable
  pre-computed data
- Example gallery rebuilt with the current package code and real TCGA data:
  ubiquitination-related gene set (9 genes x 5 cancers) + single-gene MUL1
  query, matching the Shiny outputs; reproducible script added at
  data-raw/make_readme_figures.R

## 0.4.0 (2026-08-30)
- GitHub 发布：源码已上传至 https://github.com/WangYY-666/UbiPanTriage-R，支持 `remotes::install_github("WangYY-666/UbiPanTriage-R")` 一键安装
- 数据分发：33 癌种预计算缓存（约 4.2 GB）不再随源码仓库推送，改由 `download_extdata()` 从 GitHub `data` 分支分片下载、校验（MD5）、拼接并安装到本地；`check_extdata()` 可随时检查数据是否就绪
- 缺失数据友好提示：`available_cancers()` / `load_cancer_data()` / `run_shiny_app()` 在数据未安装时给出明确指引
- README 全面重写：安装、快速开始、函数速查表、评分原理、示例代码与可视化图
- 新增函数：`download_extdata()`、`check_extdata()`
# UbiPanTriage 更新日志 (Changelog)

## 0.3.0 (2026-08-26)
- 泛癌数据：新增 33 个 TCGA 癌种缓存（全基因 count/TPM/FPKM 整数编码矩阵 + 预计算子分）与泛癌预计算表 `pancancer_ubi_scores.rds`（33 × 1553 泛素基因、182 列全子分）
- 界面重构为 5 页：① 泛素基因集评分排序 ② 泛素单基因查询 ③ 自定义基因集 ④ 自定义单基因 ⑤ 使用说明；默认英文，右上角一键切换简体中文
- 新增维度：泛素功能层级分（E1/E2/E3/DUB/UBD/ULD/Other）作为第 1 维，类型全称展示（如 Ubiquitin ligase (E3)）
- 新增可视化：四维合一点图、分面热图、基因→癌种→推荐方向桑基图、UpSet（基础+免疫+代谢高分癌种）、表达比较图（TPM/count/FPKM、肿瘤±正常）
- 权重体系：维度权重 0.1–1.0、一位小数、总和为 1；基础子项 / 免疫细胞 / 代谢途径相对权重可调；跨癌种汇总支持等权平均（默认）与按样本量加权；免疫方法默认 TIMER、五种方法内部等权
- 缺失数据处理：无正常组织/差异癌种（ACC/DLBC/LAML/LGG/MESO/OV/TGCT/UCS/UVM）基础分重归一化并以 * 标注；LAML 免疫分回退 CIBERSORT；选中缺失癌种时弹窗/顶部提示
- 所有评分表统一三位小数，评分表/详情（xlsx）均可下载
- 启动器：`run_shiny_app()` 浏览器打开失败不再崩溃（打印 URL 继续运行）；新增 `LaunchUbiPanTriage.bat` 双击启动
- 测试：新增 5 页 Shiny 回归测试（testServer），修复 `stats::tapply`（本机 R 中 tapply 位于 base）、UpSet 单集合/缺列、相关性热图显著性星号、fmsb 雷达图数据框等问题

## 0.1.0 (2026-08-24)
- 初始版本：基于 LUAD 示例数据的单癌种三维评分系统
- 评分函数：calc_basic_score / calc_immune_score / calc_metabolic_score / score_genes / score_genes_multi
- 可视化：plot_radar / plot_immune_cor / plot_metabolic_cor / plot_expr_box / plot_roc_diag / plot_km / plot_ranking / plot_score_heatmap
- 推论：gene_inference / gene_inference_pancancer
- 交互界面：run_shiny_app（基因集排序 / 单基因查询 / 自定义基因集 / 方法说明）
- 数据：内置 LUAD 预计算缓存 inst/extdata/LUAD_cache.rds
