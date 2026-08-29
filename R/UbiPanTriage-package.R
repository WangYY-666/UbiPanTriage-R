#' UbiPanTriage: Pan-Cancer Ubiquitin Immunometabolic Scoring
#'
#' Computes a four-dimension (Ubiquitin / Basic / Immune / Metabolic)
#' scoring system for user-defined gene sets in TCGA pan-cancer
#' transcriptomes, ranks genes across cancers, generates publication-style
#' figures (radar, correlation heatmaps, boxplots, ROC, KM survival,
#' alluvial, UpSet), and provides an interactive Shiny interface plus
#' automatic bilingual (EN/ZH) research-hypothesis inference.
#' 泛癌泛素免疫代谢四维评分系统：评分、排序、可视化与交互式界面。
#'
#' @keywords internal
#' @importFrom bslib page_navbar nav_panel navset_bar page_sidebar layout_sidebar card sidebar
#' @importFrom DT DTOutput datatable renderDT
#' @importFrom dplyr filter mutate select arrange group_by summarise bind_rows left_join right_join inner_join full_join rename distinct pull case_when n_distinct if_else desc across
#' @importFrom grid grid.newpage grid.draw gpar textGrob unit viewport pushViewport popViewport
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook read.xlsx write.xlsx
#' @importFrom png readPNG
#' @importFrom readxl read_excel excel_sheets
#' @importFrom scales percent label_percent
#' @importFrom stringr str_detect str_replace str_remove str_extract str_split str_to_upper str_to_lower str_trim
#' @importFrom tibble tibble as_tibble
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom timeROC timeROC
#' @importFrom stats aggregate wilcox.test
#' @importFrom utils head
"_PACKAGE"

# Declare tidy-eval pronouns and other unquoted symbols used in NSE code
# (声明 tidy-eval 代词等非标准求值符号，供 R CMD check 识别)
utils::globalVariables(c(".data", "stratum"))
