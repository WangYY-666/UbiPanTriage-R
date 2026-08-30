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
