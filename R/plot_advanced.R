# ============================================================================
# plot_advanced.R - SCI-style plots for the pan-cancer interface
# (expression comparison, 4D dot plot, facet heatmap, UpSet, radar)
# ============================================================================

#' Format numbers with three decimals (统一保留三位小数)
#' @param x numeric vector (数值向量)
#' @return character vector with exactly 3 decimals (保留三位小数的字符)
#' @noRd
.fmt3 <- function(x) sprintf("%.3f", as.numeric(x))

#' Decode integer-encoded expression (log2(x+1)*100 -> x)
#' @noRd
.decode_expr <- function(v) {
  out <- 2^(as.numeric(v) / 100) - 1
  names(out) <- names(v)
  out
}

#' Raw per-sample gene expression in a chosen unit (原始表达矩阵取值)
#' @param cache cancer cache (癌种缓存)
#' @param gene gene symbol (基因符号)
#' @param unit "tpm" (default) / "count" / "fpkm"
#' @return named numeric vector in the raw scale (原始量纲的命名向量)
#' @export
gene_expr_raw <- function(cache, gene, unit = c("tpm", "count", "fpkm")) {
  unit <- match.arg(unit)
  m <- switch(unit,
              tpm   = cache$expr_tpm_full,
              count = cache$expr_count,
              fpkm  = cache$expr_fpkm)
  if (is.null(m) || !gene %in% rownames(m))
    return(stats::setNames(rep(NA_real_, ncol(m)), colnames(m)))
  .decode_expr(m[gene, , drop = TRUE])
}

#' log2(x+1) gene expression for display (log2(x+1) 显示尺度)
#' @param cache cancer cache (癌种缓存)
#' @param gene gene symbol (基因符号)
#' @param unit expression unit: "tpm" (default) / "count" / "fpkm"
#' @return named numeric vector of log2(expression + 1) (log2 表达命名向量)
#' @export
gene_expr_log2 <- function(cache, gene, unit = c("tpm", "count", "fpkm")) {
  log2(gene_expr_raw(cache, gene, unit) + 1)
}

#' Expression data.frame for the tumor / normal comparison plot
#' @noRd
expr_compare_df <- function(cache, gene, unit = "tpm", show_normal = TRUE) {
  v <- gene_expr_log2(cache, gene, unit)
  d <- data.frame(sample = names(v), expr = as.numeric(v), stringsAsFactors = FALSE)
  d <- merge(d, cache$sample_info[, c("sample", "type"), drop = FALSE], by = "sample")
  d <- d[d$type %in% c("Tumor", "Normal"), , drop = FALSE]
  if (!show_normal) d <- d[d$type == "Tumor", , drop = FALSE]
  d[!is.na(d$expr), , drop = FALSE]
}

#' Tumor vs normal expression boxplot (ggpubr style, Wilcoxon significance)
#'
#' SCI-style boxplot with jittered points and a Wilcoxon rank-sum test;
#' significance is shown with the standard star notation (ns / * / ** / ***).
#' @param cache cancer cache (癌种缓存)
#' @param gene gene symbol (基因符号)
#' @param unit expression unit: "tpm" (default) / "count" / "fpkm"
#' @param show_normal whether to include normal tissue samples (是否包含正常组织)
#' @return a ggplot object
#' @export
plot_expr_compare <- function(cache, gene, unit = c("tpm", "count", "fpkm"),
                              show_normal = TRUE) {
  unit <- match.arg(unit)
  d <- expr_compare_df(cache, gene, unit, show_normal)
  unit_lab <- switch(unit, tpm = "log2(TPM+1)", count = "log2(Count+1)", fpkm = "log2(FPKM+1)")
  title <- paste0(gene, " - expression (", unit, ") in ", cache$cancer)
  if (nrow(d) == 0)
    return(.plot_na(paste0(gene, " not detected in ", cache$cancer, " (", unit, ")")))
  d$type <- factor(d$type, levels = c("Normal", "Tumor"))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$type, y = .data$expr, fill = .data$type)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.55) +
    ggplot2::geom_jitter(width = 0.18, size = 0.4, alpha = 0.35, color = "grey30") +
    ggplot2::scale_fill_manual(values = c(Tumor = "#C0392B", Normal = "#5DADE2"), guide = "none") +
    .theme_ubi(base_size = 12) +
    ggplot2::labs(title = title, x = NULL, y = unit_lab) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12),
                   axis.title.y = ggplot2::element_text(size = 11))
  if (show_normal && length(unique(d$type)) == 2) {
    pv <- tryCatch(stats::wilcox.test(expr ~ type, data = d)$p.value, error = function(e) NA_real_)
    star <- if (is.na(pv)) "" else if (pv < 0.001) "***" else if (pv < 0.01) "**" else if (pv < 0.05) "*" else "ns"
    if (nzchar(star)) {
      ymax <- max(d$expr, na.rm = TRUE)
      p <- p + ggplot2::annotate("segment", x = 1, xend = 2, y = ymax * 1.12, yend = ymax * 1.12,
                                 linewidth = 0.5, color = "grey30") +
        ggplot2::annotate("text", x = 1.5, y = ymax * 1.22, label = star, size = 5.5, color = "grey20") +
        ggplot2::coord_cartesian(ylim = c(min(d$expr, na.rm = TRUE), ymax * 1.35))
    } else {
      p <- p + ggplot2::labs(caption = if (is.na(pv)) "p = NA" else paste0("Wilcoxon p = ", format.pval(pv, digits = 3)))
    }
  }
  p
}

#' Four-dimension combined dot plot (四维合一点图)
#'
#' One row per gene, one column per dimension; point color and size encode
#' the score, genes are ordered by the combined score.
#' @param d long data.frame with columns Gene, Dim, Score, Combined
#' @param gene_col column name of the gene (基因列名)
#' @param dim_col column name of the dimension (维度列名)
#' @param score_col column name of the score (评分列名)
#' @param comb_col column name of the combined score used for ordering (综合分列名)
#' @return a ggplot object
#' @export
plot_4d_dot <- function(d, gene_col = "Gene", dim_col = "Dim",
                        score_col = "Score", comb_col = "Combined") {
  d <- d[!is.na(d[[score_col]]), , drop = FALSE]
  if (nrow(d) == 0) return(.plot_na("No scores available"))
  ord <- if (comb_col %in% colnames(d) && any(!is.na(d[[comb_col]]))) {
    v <- stats::aggregate(d[[comb_col]], by = list(g = d[[gene_col]]), mean, na.rm = TRUE)
    as.character(v$g[order(-v$x)])
  } else {
    sort(unique(d[[gene_col]]))
  }
  d[[gene_col]] <- factor(d[[gene_col]], levels = rev(ord))
  d[[dim_col]] <- factor(d[[dim_col]], levels = unique(d[[dim_col]]))
  ggplot2::ggplot(d, ggplot2::aes(x = .data[[dim_col]], y = .data[[gene_col]],
                                  color = .data[[score_col]], size = .data[[score_col]])) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::scale_color_gradient2(low = "#2E86C1", mid = "#F7F7F7", high = "#C0392B",
                                   midpoint = 0.5, limits = c(0, 1), name = "Score") +
    ggplot2::scale_size_continuous(range = c(2, 7), limits = c(0, 1), guide = "none") +
    .theme_ubi(base_size = 12) +
    ggplot2::labs(title = "Dimension scores across selected cancers",
                  x = NULL, y = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
                   legend.position = "right")
}

#' Gene x cancer score heatmap, one facet per dimension (分面热图)
#' @param d long data.frame with columns Gene, Cancer, Dim, Score
#' @param score_col column name of the score (评分列名)
#' @param dim_col column name of the dimension (维度列名)
#' @param gene_col column name of the gene (基因列名)
#' @param cancer_col column name of the cancer (癌种列名)
#' @return a ggplot object
#' @export
plot_facet_heatmap <- function(d, score_col = "Score", dim_col = "Dim",
                               gene_col = "Gene", cancer_col = "Cancer") {
  d <- d[!is.na(d[[score_col]]), , drop = FALSE]
  if (nrow(d) == 0) return(.plot_na("No scores available"))
  d[[gene_col]] <- factor(d[[gene_col]], levels = sort(unique(d[[gene_col]])))
  d[[cancer_col]] <- factor(d[[cancer_col]], levels = sort(unique(d[[cancer_col]])))
  ggplot2::ggplot(d, ggplot2::aes(x = .data[[cancer_col]], y = .data[[gene_col]],
                                  fill = .data[[score_col]])) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::scale_fill_gradient2(low = "#2E86C1", mid = "white", high = "#C0392B",
                                  midpoint = 0.5, limits = c(0, 1), name = "Score") +
    ggplot2::facet_wrap(stats::as.formula(paste0("~ ", dim_col)), nrow = 1) +
    .theme_ubi(base_size = 11, axes = FALSE) +
    ggplot2::labs(title = "Gene x cancer scores", x = NULL, y = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
                   strip.background = ggplot2::element_rect(fill = "#EEF3F8", color = NA),
                   strip.text = ggplot2::element_text(face = "bold", size = 11))
}


#' UpSet plot for one gene: cancers with high basic / immune / metabolic (单基因 UpSet)
#'
#' Three sets: cancers where the gene's basic / immune / metabolic score is
#' above the threshold. Default threshold is the cancer-specific universe
#' median; a fixed 0-1 cutoff can be supplied instead.
#' @param d one row per cancer: columns Cancer, Basic_Score, Immune_Score,
#'   Metabolic_Score
#' @param medians named numeric vector (per cancer, per dimension) used when
#'   threshold is NULL; names like "LUAD.Basic"
#' @param threshold optional fixed cutoff; when NULL the median rule is used
#' @return an UpSetR plot object (printed by renderPlot)
#' @export
plot_upset_3dim <- function(d, medians = NULL, threshold = NULL) {
  if (!requireNamespace("UpSetR", quietly = TRUE))
    return(.plot_na("package 'UpSetR' is required"))
  dims <- c("Basic_Score", "Immune_Score", "Metabolic_Score")
  sets <- matrix(FALSE, nrow = nrow(d), ncol = 3,
                 dimnames = list(NULL, c("Basic", "Immune", "Metabolic")))
  for (j in seq_along(dims)) {
    nm <- sub("_Score$", "", dims[j])
    if (!is.null(threshold)) {
      sets[, nm] <- !is.na(d[[dims[j]]]) & d[[dims[j]]] > threshold
    } else if (!is.null(medians)) {
      key <- paste0(d$Cancer, ".", nm)
      sets[, nm] <- !is.na(d[[dims[j]]]) & d[[dims[j]]] > medians[key]
    }
  }
  combos <- apply(sets, 1, function(r) paste(names(which(r)), collapse = "&"))
  combos[combos == ""] <- "none"
  tt <- table(combos)
  expr <- stats::setNames(as.numeric(tt), names(tt))
  all6 <- c("Basic", "Immune", "Metabolic", "Basic&Immune", "Basic&Metabolic",
            "Immune&Metabolic", "Basic&Immune&Metabolic")
  expr <- expr[intersect(all6, names(expr))]
  if (length(expr) == 0) return(.plot_na("No cancer passes the high-score rule"))
  fe <- UpSetR::fromExpression(expr)
  sets_used <- intersect(c("Basic", "Immune", "Metabolic"), colnames(fe))
  if (length(sets_used) < 2)
    return(.plot_na("Fewer than two high dimensions among the selected cancers"))
  UpSetR::upset(fe,
                sets = sets_used, order.by = "freq",
                main.bar.color = "#2C5FA8", sets.bar.color = "#8B9EB7",
                mainbar.y.label = "Cancers", sets.x.label = "Set size",
                point.size = 3.2, line.size = 0.9)
}

#' Filled radar chart (fmsb style) of the scoring dimensions (fmsb 雷达图)
#' @param values named numeric vector, e.g. c(Basic = 0.8, Immune = 0.5,
#'   Metabolic = 0.6, Ubiquitin = 0.9)
#' @param gene gene name used in the title (基因名)
#' @param max_val outer scale maximum (外圈最大值)
#' @return invisible NULL; draws to the active device (绘制到当前设备)
#' @export
plot_radar_fmsb <- function(values, gene = "", max_val = 1) {
  if (!requireNamespace("fmsb", quietly = TRUE)) {
    warning("package 'fmsb' is required for the radar chart; using the ggplot fallback")
    return(plot_radar(values, gene, max_val))
  }
  v <- as.numeric(values); names(v) <- names(values)
  v[is.na(v)] <- 0
  if (length(v) < 3) v <- c(v, stats::setNames(rep(0, 3 - length(v)),
                                               paste0("Dim", 1:(3 - length(v)))))
  dat <- as.data.frame(rbind(rep(max_val, length(v)), rep(0, length(v)), v))
  colnames(dat) <- names(v)
  graphics::par(mar = c(1.2, 1.2, 2.2, 1.2))
  fmsb::radarchart(dat, axistype = 1, seg = 4, pcol = "#C0392B", pfcol = "#C0392B33",
                   plwd = 2.2, cglcol = "grey75", cglty = 1, axislabcol = "grey50",
                   caxislabels = seq(0, max_val, length.out = 5), vlcex = 1.05)
  if (nzchar(gene))
    graphics::title(main = paste0(gene, " - dimension scores"), font.main = 2)
  graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))
  invisible(NULL)
}

#' Score bubble plot (评分气泡图: x = immune, y = metabolic, size = basic)
#'
#' SCI-style bubble chart: x = immune score, y = metabolic score, bubble size =
#' basic score; point color optionally encodes the ubiquitin class
#' (E1 / E2 / E3 / DUB / UBD / ULD / other, shown with full names).
#' 评分气泡图：横轴为免疫分、纵轴为代谢分、气泡大小为基础分；
#' 气泡颜色可表示泛素酶类别（全称显示）。
#' @param d data.frame with Basic_Score, Immune_Score, Metabolic_Score
#'   (含三维评分的汇总表)
#' @param type_col optional column name of the ubiquitin class full name
#'   (可选：泛素类别全称列名)
#' @return a ggplot object
#' @export
plot_3d_bubble <- function(d, type_col = NULL) {
  need <- c("Basic_Score", "Immune_Score", "Metabolic_Score")
  have <- need[vapply(need, function(n) n %in% colnames(d) && any(!is.na(d[[n]])), logical(1))]
  if (length(have) == 0) return(.plot_na("No scores available"))
  dd <- d[stats::complete.cases(d[, have, drop = FALSE]), , drop = FALSE]
  if (nrow(dd) == 0) return(.plot_na("No scores available"))
  n_drop <- nrow(d) - nrow(dd)
  missing_dim <- setdiff(need, have)
  cap <- character()
  if (length(missing_dim) > 0)
    cap <- c(cap, paste0("Missing dimension: ", paste(missing_dim, collapse = ", "),
                         " (select all three dimensions for the 3D bubble)"))
  if (n_drop > 0)
    cap <- c(cap, paste0(n_drop, " gene(s) with incomplete scores not shown"))
  cap <- paste(cap, collapse = "; ")
  cap <- if (nzchar(cap)) cap else NULL
  x <- if ("Immune_Score" %in% have) "Immune_Score" else have[1]
  y <- if ("Metabolic_Score" %in% have) "Metabolic_Score" else
         if ("Immune_Score" %in% have) "Immune_Score" else have[1]
  size_col <- if ("Basic_Score" %in% have) "Basic_Score" else NULL
  if (length(have) == 1) {
    v <- have[1]
    if (!is.null(type_col) && type_col %in% colnames(dd)) {
      dd[[type_col]][is.na(dd[[type_col]]) | !nzchar(dd[[type_col]])] <- "Other"
      p <- ggplot2::ggplot(dd, ggplot2::aes(x = stats::reorder(.data[[v]], .data[[v]]),
                                            y = .data[[v]], color = .data[[type_col]])) +
        ggplot2::geom_point(size = 3, alpha = 0.85) + ggplot2::coord_flip() +
        ggplot2::scale_color_brewer(palette = "Set2", name = "Ubiquitin class")
    } else {
      p <- ggplot2::ggplot(dd, ggplot2::aes(x = stats::reorder(.data[[v]], .data[[v]]),
                                            y = .data[[v]])) +
        ggplot2::geom_point(color = "#3B7DD8", size = 3, alpha = 0.85) +
        ggplot2::coord_flip()
    }
    return(p + ggplot2::labs(x = "Gene", y = v, title = paste0(v, " scores"),
                             caption = cap) + .theme_ubi(base_size = 12))
  }
  p <- ggplot2::ggplot(dd, ggplot2::aes(x = .data[[x]], y = .data[[y]]))
  if (!is.null(type_col) && type_col %in% colnames(dd)) {
    dd[[type_col]][is.na(dd[[type_col]]) | !nzchar(dd[[type_col]])] <- "Other"
    p <- p + ggplot2::aes(color = .data[[type_col]]) +
      ggplot2::geom_point(alpha = 0.75) +
      ggplot2::scale_color_brewer(palette = "Set2", name = "Ubiquitin class")
  } else {
    p <- p + ggplot2::geom_point(color = "#3B7DD8", alpha = 0.75)
  }
  if (!is.null(size_col))
    p <- p + ggplot2::aes(size = .data[[size_col]]) +
      ggplot2::scale_size_continuous(range = c(2, 11), name = "Basic score")
  p + ggplot2::labs(x = x, y = y, title = paste0(paste(have, collapse = " / "), " scores"),
                    caption = cap) +
    .theme_ubi(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

#' Ranking bar chart with score labels (带得分标注的排名条形图)
#'
#' Genes are ordered by the score column; the score is printed at the end of
#' each bar so the exact value is readable (3 decimals).
#' 按评分排序的条形图，每条 bar 尾部标注三位小数得分。
#' @param d data.frame containing the gene and score columns (含基因列与评分列)
#' @param value score column used for ranking (用于排序的评分列)
#' @param gene_col gene column (基因列)
#' @param top_n number of genes to show (显示基因数)
#' @param label label text for the score axis and title (评分轴/标题标签)
#' @return a ggplot object
#' @export
plot_ranking_bar <- function(d, value = "Combined_Score", gene_col = "Gene",
                             top_n = 20, label = NULL) {
  if (!value %in% colnames(d) || !gene_col %in% colnames(d))
    return(.plot_na("Score column not available"))
  if (is.null(label)) label <- value
  dd <- d[!is.na(d[[value]]), , drop = FALSE]
  dd <- dd[order(dd[[value]], decreasing = TRUE), , drop = FALSE]
  dd <- head(dd, top_n)
  if (nrow(dd) == 0) return(.plot_na("No scores available"))
  dd[[gene_col]] <- factor(dd[[gene_col]], levels = rev(dd[[gene_col]]))
  ggplot2::ggplot(dd, ggplot2::aes(x = .data[[gene_col]], y = .data[[value]])) +
    ggplot2::geom_col(fill = "#3B7DD8", alpha = 0.85, width = 0.68) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", .data[[value]])),
                       hjust = -0.15, size = 3, color = "grey25") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.13))) +
    ggplot2::coord_flip() +
    .theme_ubi(base_size = 12) +
    ggplot2::labs(title = paste0("Top ", nrow(dd), " genes by ", label),
                  x = "Gene", y = label) +
    ggplot2::theme(legend.position = "none",
                   panel.grid.major.y = ggplot2::element_blank())
}

#' Ranking bar charts for all four score dimensions (四维排名条形图组合)
#'
#' One horizontal bar chart per dimension (basic / metabolic / immune /
#' combined), each ranked by its own score, arranged as a 2 x 2 patchwork.
#' 基础 / 代谢 / 免疫 / 综合四个维度各绘制一张横向条形图，按各自评分排序，
#' 以 2 x 2 拼图展示。
#' @param d data.frame containing the gene and the four score columns
#'   (含基因列与四维评分列的汇总表)
#' @param top_n number of genes to show per chart (每个子图显示的基因数)
#' @return a ggplot / patchwork object
#' @export
plot_ranking_grid <- function(d, top_n = 20) {
  cols <- c(Basic = "Basic_Score", Metabolic = "Metabolic_Score",
            Immune = "Immune_Score", Combined = "Combined_Score")
  cols <- cols[cols %in% colnames(d)]
  if (length(cols) == 0) return(.plot_na("No ranking scores available"))
  ps <- lapply(names(cols), function(nm) {
    plot_ranking_bar(d, value = cols[[nm]], gene_col = "Gene",
                     top_n = top_n, label = nm)
  })
  if (length(ps) == 1) return(ps[[1]])
  patchwork::wrap_plots(ps, ncol = 2)
}

#' Pan-cancer expression violin plot (跨癌种表达小提琴图)
#'
#' One panel for all selected cancers; x = cancer, y = log2(unit + 1).
#' When normal samples are included, tumor and normal are drawn side by
#' side per cancer with Wilcoxon significance stars (* p < 0.05,
#' ** p < 0.01, *** p < 0.001).
#' 所有所选癌种画在一张图上：横轴为癌种、纵轴为 log2(单位+1) 表达值；
#' 若包含正常组织，每个癌种内肿瘤与正常并排显示并标注 Wilcoxon 显著性星号。
#' @param caches named list of cancer caches (癌种缓存命名列表)
#' @param gene gene symbol (基因符号)
#' @param unit "tpm" (default) / "count" / "fpkm"
#' @param show_normal include normal samples (是否包含正常组织)
#' @return a ggplot object
#' @export
plot_expr_pancancer <- function(caches, gene, unit = c("tpm", "count", "fpkm"),
                                show_normal = TRUE) {
  unit <- match.arg(unit)
  rows <- lapply(names(caches), function(ca) {
    cache <- caches[[ca]]
    v <- tryCatch(gene_expr_log2(cache, gene, unit), error = function(e) NULL)
    if (is.null(v) || all(is.na(v))) return(NULL)
    d <- data.frame(sample = names(v), expr = as.numeric(v), stringsAsFactors = FALSE)
    d <- merge(d, cache$sample_info[, c("sample", "type"), drop = FALSE], by = "sample")
    d <- d[d$type %in% c("Tumor", "Normal"), , drop = FALSE]
    if (!show_normal) d <- d[d$type == "Tumor", , drop = FALSE]
    if (nrow(d) == 0) return(NULL)
    d$Cancer <- ca
    d
  })
  d <- do.call(rbind, rows)
  if (is.null(d) || nrow(d) == 0)
    return(.plot_na(paste0(gene, " not detected in the selected cancers")))
  unit_lab <- switch(unit, tpm = "log2(TPM+1)", count = "log2(Count+1)",
                     fpkm = "log2(FPKM+1)")
  d$Cancer <- factor(d$Cancer, levels = sort(unique(d$Cancer)))
  ymax <- max(d$expr, na.rm = TRUE)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$Cancer, y = .data$expr))
  if (show_normal && any(d$type == "Normal")) {
    d$type <- factor(d$type, levels = c("Normal", "Tumor"))
    p <- p + ggplot2::aes(fill = .data$type) +
      ggplot2::geom_violin(position = ggplot2::position_dodge(0.8), trim = TRUE,
                           width = 0.9, alpha = 0.65, linewidth = 0.3) +
      ggplot2::geom_boxplot(position = ggplot2::position_dodge(0.8), width = 0.13,
                            outlier.shape = NA, alpha = 0.9, linewidth = 0.3) +
      ggplot2::scale_fill_manual(values = c(Tumor = "#C0392B", Normal = "#5DADE2"),
                                 name = NULL)
    stars <- vapply(levels(d$Cancer), function(ca) {
      sub <- d[d$Cancer == ca, , drop = FALSE]
      if (length(unique(sub$type)) < 2 || sum(sub$type == "Tumor") < 3 ||
          sum(sub$type == "Normal") < 3) return("")
      pv <- tryCatch(stats::wilcox.test(expr ~ type, data = sub)$p.value,
                     error = function(e) NA_real_)
      if (is.na(pv)) "" else if (pv < 0.001) "***" else if (pv < 0.01) "**"
      else if (pv < 0.05) "*" else ""
    }, character(1))
    sd <- data.frame(Cancer = names(stars)[nzchar(stars)], star = stars[nzchar(stars)],
                     stringsAsFactors = FALSE)
    if (nrow(sd) > 0) {
      sd$x <- match(sd$Cancer, levels(d$Cancer))
      p <- p + ggplot2::annotate("text", x = sd$x, y = ymax * 1.08, label = sd$star,
                                 size = 3.8, color = "grey20") +
        ggplot2::coord_cartesian(ylim = c(min(d$expr, na.rm = TRUE), ymax * 1.22))
    }
  } else {
    p <- p + ggplot2::geom_violin(fill = "#C0392B", alpha = 0.6, trim = TRUE,
                                  width = 0.85, linewidth = 0.3) +
      ggplot2::geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white",
                            color = "#8B1A1A", alpha = 0.9, linewidth = 0.3)
  }
  p + .theme_ubi(base_size = 11) +
    ggplot2::labs(title = paste0(gene, " - expression across cancers (", unit, ")"),
                  x = NULL, y = unit_lab) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
                   legend.position = "top")
}
