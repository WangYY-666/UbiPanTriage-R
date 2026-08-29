# ============================================================================
# inference.R - automatic gene inference (bilingual, conservative)
# ============================================================================

#' Single-cancer gene inference (single-cancer inference, bilingual)
#'
#' A conservative, evidence-based summary for one gene in one cancer:
#' dimension strengths, significant survival association, diagnostic value
#' and a first validation suggestion.
#' @param cache cancer cache
#' @param gene gene symbol
#' @param result output of \code{\link{score_genes}}
#' @param lang language: "en" (default) or "zh"
#' @return character string
#' @export
gene_inference <- function(cache, gene, result, lang = c("en", "zh")) {
  lang <- match.arg(lang)
  sc <- result$scores[result$scores$Gene == gene, , drop = FALSE]
  if (nrow(sc) == 0) return("")
  en <- lang == "en"
  L <- list()
  L$title <- if (en) paste0("[", gene, " inference in ", cache$cancer, "]") else
    paste0("[", gene, " 在 ", cache$cancer, " 中的推论]")
  ty <- sc$Ubi_Type_Full[1]
  L$type <- if (!is.na(ty)) {
    if (en) paste0("Ubiquitin class: ", ty, ".") else paste0("泛素类型：", ty, "。")
  } else ""
  if ("Basic_Score" %in% colnames(sc)) {
    L$basic <- if (en) paste0("Basic score = ", round(sc$Basic_Score, 3),
                              " (expression / differential / survival / diagnostic evidence).")
    else paste0("基础分 = ", round(sc$Basic_Score, 3), "（表达 / 差异 / 生存 / 诊断证据）。")
  }
  if ("Immune_Score" %in% colnames(sc)) {
    L$immune <- if (en) paste0("Immune score = ", round(sc$Immune_Score, 3),
                               "; top immune cell: ", sc$Top_immune_cell, " (r = ",
                               round(sc$Cor_value1, 3), ").")
    else paste0("免疫分 = ", round(sc$Immune_Score, 3), "；主要相关免疫细胞：",
                sc$Top_immune_cell, "（r = ", round(sc$Cor_value1, 3), "）。")
  }
  if ("Metabolic_Score" %in% colnames(sc)) {
    L$metab <- if (en) paste0("Metabolic score = ", round(sc$Metabolic_Score, 3),
                              "; top pathway: ", sc$Top_metabolic_trait, " (r = ",
                              round(sc$Cor_value2, 3), ").")
    else paste0("代谢分 = ", round(sc$Metabolic_Score, 3), "；主要相关代谢途径：",
                sc$Top_metabolic_trait, "（r = ", round(sc$Cor_value2, 3), "）。")
  }
  if ("HR" %in% colnames(sc) && !is.na(sc$HR) && !is.na(sc$HR_p) && sc$HR_p < 0.05) {
    role <- if (sc$direction == "prot") (if (en) "protective" else "抑癌") else (if (en) "risk" else "促癌")
    L$surv <- if (en) paste0("Significant Cox survival association (HR = ", round(sc$HR, 3),
                             ", p = ", format.pval(sc$HR_p, digits = 3), ", ", role, ").")
    else paste0("存在显著 Cox 生存关联（HR = ", round(sc$HR, 3), "，p = ",
                format.pval(sc$HR_p, digits = 3), "，", role, "）。")
  }
  if ("AUC_diag" %in% colnames(sc) && !is.na(sc$AUC_diag)) {
    L$diag <- if (en) paste0("Diagnostic AUC (tumor vs normal) = ", round(sc$AUC_diag, 3), ".")
    else paste0("诊断 AUC（癌 vs 正常）= ", round(sc$AUC_diag, 3), "。")
  }
  if (!is.null(sc$Missing) && !is.na(sc$Missing)) {
    L$missing <- if (en) paste0("Missing evidence in this cancer: ", sc$Missing, " (basic score re-normalized).")
    else paste0("该癌种缺失证据：", sc$Missing, "（基础分已重归一化）。")
  }
  L$rec <- if (en) paste0("Suggestion: validate the gene in ", cache$cancer,
                          " first; confirm the baseline expression before functional assays.")
  else paste0("建议：优先在 ", cache$cancer, " 开展验证实验，先确认基础表达再设计功能实验。")

  out <- c(L$title, L$type, L$basic, L$immune, L$metab, L$surv, L$diag, L$missing, L$rec)
  paste(out[nzchar(out)], collapse = "\n")
}

#' Pan-cancer gene inference (pan-cancer inference, bilingual)
#'
#' Summarizes one gene across multiple cancers: dimension strengths,
#' expression/survival directions per cancer, and the recommended cancer
#' for the first validation experiment.
#' @param caches named list of cancer caches
#' @param gene gene symbol
#' @param result_multi output of \code{\link{score_genes_multi}}
#' @param lang language: "en" (default) or "zh"
#' @return character string
#' @export
gene_inference_pancancer <- function(caches, gene, result_multi, lang = c("en", "zh")) {
  lang <- match.arg(lang)
  long <- result_multi$scores_long
  d <- long[long$Gene == gene, , drop = FALSE]
  if (nrow(d) == 0) {
    return(if (lang == "en") paste0("Gene ", gene, " was not detected in any selected cancer.")
           else paste0("基因 ", gene, " 在所选的癌种中均未检出。"))
  }
  has_u <- "Ubi_Score" %in% colnames(d)
  has_b <- "Basic_Score" %in% colnames(d)
  has_i <- "Immune_Score" %in% colnames(d)
  has_m <- "Metabolic_Score" %in% colnames(d)

  pick_top <- function(col, n = 3) {
    dd <- d[!is.na(d[[col]]), , drop = FALSE]
    if (nrow(dd) == 0) return(if (lang == "en") "none" else "无")
    dd <- dd[order(dd[[col]], decreasing = TRUE), , drop = FALSE]
    paste(head(dd$Cancer, n), collapse = ", ")
  }
  en <- lang == "en"
  out <- if (en) paste0("[", gene, " pan-cancer inference]") else paste0("[", gene, " 泛癌推论]")
  if (has_u) {
    ty <- d$Ubi_Type_Full[!is.na(d$Ubi_Type_Full)][1]
    if (length(ty) == 0 || is.na(ty)) ty <- if (en) "not annotated" else "未注释"
    out <- c(out, if (en) paste0("Ubiquitin class: ", ty, ".") else paste0("泛素类型：", ty, "。"))
  }

  # survival direction per cancer (significant Cox only)
  dir_txt <- character()
  for (ca in d$Cancer) {
    res <- result_multi$per_cancer[[ca]]
    sc <- res$scores[res$scores$Gene == gene, , drop = FALSE]
    if (nrow(sc) == 0) next
    hr <- if (!is.null(sc$HR)) sc$HR[1] else NA_real_
    hp <- if (!is.null(sc$HR_p)) sc$HR_p[1] else NA_real_
    dirn <- if (!is.null(sc$direction)) sc$direction[1] else NA_character_
    if (!is.na(hp) && hp < 0.05 && !is.na(hr) && !is.na(dirn)) {
      role <- if (dirn == "prot") (if (en) "protective" else "抑癌") else (if (en) "risk" else "促癌")
      dir_txt <- c(dir_txt, paste0(ca, ": ", role))
    }
  }
  out <- c(out, if (length(dir_txt) > 0) {
    if (en) paste0("Significant survival associations (Cox p < 0.05): ", paste(dir_txt, collapse = "; "), ".")
    else paste0("显著的生存关联（Cox p < 0.05）：", paste(dir_txt, collapse = "；"), "。")
  } else if (en) {
    "No significant Cox survival association was found in the selected cancers."
  } else {
    "所选癌种中未发现显著的 Cox 生存关联。"
  })

  if (has_b) out <- c(out, if (en) paste0("Top cancers by basic score: ", pick_top("Basic_Score"), ".")
                      else paste0("基础分最高的癌种：", pick_top("Basic_Score"), "。"))
  if (has_i) out <- c(out, if (en) paste0("Top cancers by immune score: ", pick_top("Immune_Score"), ".")
                      else paste0("免疫分最高的癌种：", pick_top("Immune_Score"), "。"))
  if (has_m) out <- c(out, if (en) paste0("Top cancers by metabolic score: ", pick_top("Metabolic_Score"), ".")
                      else paste0("代谢分最高的癌种：", pick_top("Metabolic_Score"), "。"))

  # recommended cancer
  rec <- d$Cancer[which.max(d$Combined_Score)]
  tips <- character()
  if (has_b && max(d$Basic_Pct, na.rm = TRUE) > 0.7)
    tips <- c(tips, if (en)
      "the high basic score suggests strong expression / differential / prognosis evidence (experiment design still needs literature and clinical judgment)"
      else "基础分较高，提示表达、差异与预后证据充分（实验设计仍需结合文献与临床判断）")
  if (has_i && max(d$Immune_Pct, na.rm = TRUE) > 0.7)
    tips <- c(tips, if (en) "the high immune score favors immune-microenvironment research"
             else "免疫分较高，适合开展免疫微环境相关研究")
  if (has_m && max(d$Metabolic_Pct, na.rm = TRUE) > 0.7)
    tips <- c(tips, if (en) "the high metabolic score favors metabolic research" else "代谢分较高，适合开展代谢相关研究")
  if (length(tips) == 0)
    tips <- if (en) "all dimensions need baseline validation first" else "各维度均需先验证基础表达"
  out <- c(out, if (en) paste0("Recommendation: ", rec, " (highest combined score; ", paste(tips, collapse = "; "),
                               "). Validate the gene in ", rec, " first.")
           else paste0("建议：综合分最高的癌种为 ", rec, "（", paste(tips, collapse = "；"),
                       "）。建议优先在 ", rec, " 中验证。"))
  paste(out, collapse = "\n")
}