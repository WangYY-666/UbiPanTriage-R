# ============================================================================
# plots.R - visualization functions (SCI-style, English labels)
# ============================================================================

#' Minimal plot theme used by the package
#' @noRd
.theme_ubi <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                   axis.text = ggplot2::element_text(color = "grey20"),
                   legend.position = "right")
}

#' Placeholder plot for missing data
#' @noRd
.plot_na <- function(msg) {
  ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_text(label = msg, size = 4.5, color = "grey40") +
    ggplot2::theme_void() +
    ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1)
}

#' Radar chart of the scoring dimensions
#'
#' @param values named numeric vector, e.g. c(Ubiquitin = 0.8, Basic = 0.5,
#'   Immune = 0.6, Metabolic = 0.5)
#' @param gene gene name used in the title
#' @param max_val outer scale maximum
#' @return a ggplot object
#' @examples
#' plot_radar(c(Basic = 0.8, Immune = 0.5, Metabolic = 0.6), gene = "MDM2")
#' @export
plot_radar <- function(values, gene = "", max_val = 1) {
  v <- as.numeric(values)
  names(v) <- names(values)
  v[is.na(v)] <- 0
  if (length(v) < 3) v <- c(v, stats::setNames(rep(0, 3 - length(v)), LETTERS[1:(3 - length(v))]))
  df <- data.frame(
    dim = factor(c(names(v), names(v)[1]), levels = names(v)),
    val = c(v, v[1])
  )
  ggplot2::ggplot(df, ggplot2::aes(x = .data$dim, y = .data$val, group = 1)) +
    ggplot2::geom_polygon(fill = "#3B7DD8", alpha = 0.35, color = "#2C5FA8", linewidth = 1) +
    ggplot2::geom_point(size = 2.5, color = "#2C5FA8") +
    ggplot2::coord_polar() +
    ggplot2::ylim(0, max_val) +
    .theme_ubi() +
    ggplot2::labs(title = if (nzchar(gene)) paste0(gene, " - dimension scores") else "Dimension scores",
                  x = NULL, y = NULL) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

#' Correlation bar plot of a gene against traits
#' @noRd
.cor_bar <- function(cors, title = NULL, fill_colors = c("#C0392B", "#2E86C1")) {
  d <- data.frame(trait = names(cors), cor = as.numeric(cors), stringsAsFactors = FALSE)
  d <- d[!is.na(d$cor), , drop = FALSE]
  if (nrow(d) == 0) return(.plot_na("No correlation data available"))
  d <- d[order(abs(d$cor), decreasing = TRUE), , drop = FALSE]
  d$trait <- factor(d$trait, levels = rev(d$trait))
  d$sign <- ifelse(d$cor >= 0, "Positive", "Negative")
  ggplot2::ggplot(d, ggplot2::aes(x = .data$trait, y = .data$cor, fill = .data$sign)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(values = stats::setNames(fill_colors, c("Positive", "Negative"))) +
    ggplot2::coord_flip() +
    .theme_ubi() +
    ggplot2::labs(title = title, x = NULL, y = "Spearman rho") +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3)
}

#' Immune cell infiltration correlation plot
#'
#' @param cache cancer cache
#' @param gene gene symbol
#' @param method immune method: timer / cibersort / mcp / ssgsea / xcell
#' @return a ggplot object
#' @export
plot_immune_cor <- function(cache, gene, method = "timer") {
  cm <- immune_cor(cache, gene, method)
  if (is.null(cm) || ncol(cm) == 0)
    stop("no immune infiltration data available for method: ", method)
  .cor_bar(cm[1, ], title = paste0(gene, " - immune infiltration correlation (", method, ")"))
}

#' Metabolic pathway correlation plot
#'
#' @param cache cancer cache
#' @param gene gene symbol
#' @return a ggplot object
#' @export
plot_metabolic_cor <- function(cache, gene) {
  cm <- metabolic_cor(cache, gene)
  if (is.null(cm) || ncol(cm) == 0)
    stop("no metabolic pathway data available")
  .cor_bar(cm[1, ], title = paste0(gene, " - metabolic pathway correlation"))
}

#' Tumor vs normal expression boxplot
#'
#' @param cache cancer cache
#' @param gene gene symbol
#' @return a ggplot object
#' @export
plot_expr_box <- function(cache, gene) {
  d <- expr_for_plot(cache, gene)
  d <- d[!is.na(d$expression), , drop = FALSE]
  if (length(unique(d$type)) < 2)
    stop("no normal tissue available for this cancer")
  pv <- tryCatch(wilcox.test(expression ~ type, data = d)$p.value, error = function(e) NA_real_)
  star <- if (is.na(pv)) "" else if (pv < 0.001) "***" else if (pv < 0.01) "**"
  else if (pv < 0.05) "*" else "ns"
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$type, y = .data$expression, fill = .data$type)) +
    ggplot2::geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
    ggplot2::scale_fill_manual(values = c(Tumor = "#C0392B", Normal = "#2E86C1")) +
    .theme_ubi() +
    ggplot2::labs(title = paste0(gene, " - tumor vs normal expression (log2 TPM+1)"),
                  x = NULL, y = "log2(TPM+1)") +
    ggplot2::theme(legend.position = "none")
  if (nzchar(star)) {
    ymax <- max(d$expression, na.rm = TRUE)
    p <- p + ggplot2::annotate("segment", x = 1, xend = 2, y = ymax * 1.1, yend = ymax * 1.1,
                               linewidth = 0.5, color = "grey30") +
      ggplot2::annotate("text", x = 1.5, y = ymax * 1.22, label = star, size = 5,
                        color = "grey20") +
      ggplot2::coord_cartesian(ylim = c(min(d$expression, na.rm = TRUE), ymax * 1.35))
  }
  p
}

#' Diagnostic ROC curve (tumor vs normal)
#'
#' @param cache cancer cache
#' @param gene gene symbol
#' @return a ggplot object
#' @export
plot_roc_diag <- function(cache, gene) {
  d <- expr_for_plot(cache, gene)
  d <- d[!is.na(d$expression) & d$type %in% c("Tumor", "Normal"), , drop = FALSE]
  if (nrow(d) < 10 || length(unique(d$type)) < 2)
    stop("not enough tumor/normal samples for a diagnostic ROC")
  roc_obj <- pROC::roc(d$type ~ d$expression, levels = c("Normal", "Tumor"),
                       direction = ">", quiet = TRUE)
  auc <- as.numeric(pROC::auc(roc_obj))
  ci <- tryCatch(as.numeric(pROC::ci.auc(roc_obj)), error = function(e) NULL)
  auc_lab <- if (!is.null(ci) && length(ci) == 3)
    paste0("AUC = ", sprintf("%.3f", auc), " (95% CI: ", sprintf("%.3f", ci[1]),
           " - ", sprintf("%.3f", ci[3]), ")")
  else paste0("AUC = ", sprintf("%.3f", auc))
  pROC::ggroc(roc_obj) +
    ggplot2::geom_abline(slope = 1, intercept = 1, linetype = 2, color = "grey50") +
    .theme_ubi() +
    ggplot2::labs(title = paste0(gene, " - diagnostic ROC (tumor vs normal)"),
                  x = "1 - Specificity", y = "Sensitivity") +
    ggplot2::annotate("text", x = 0.7, y = 0.2, size = 4.5,
                      label = auc_lab)
}

#' Kaplan-Meier survival curves by high/low expression
#'
#' Samples are split by the median expression (log2 TPM+1) in tumors.
#'
#' @param cache cancer cache
#' @param gene gene symbol
#' @param palette colors for high/low groups
#' @return a survminer ggsurvplot object
#' @export
plot_km <- function(cache, gene, palette = c("#C0392B", "#2E86C1")) {
  v <- gene_expression(cache, gene)
  d <- data.frame(sample = substr(names(v), 1, 12), expr = as.numeric(v),
                  stringsAsFactors = FALSE)
  d <- d[!duplicated(d$sample), , drop = FALSE]
  d <- merge(d, cache$clinical, by = "sample")
  d <- d[!is.na(d$expr) & !is.na(d$OS.time) & d$OS.time > 0, , drop = FALSE]
  if (nrow(d) < 20 || sum(d$OS, na.rm = TRUE) < 5)
    stop("not enough survival events for a KM analysis")
  med <- stats::median(d$expr, na.rm = TRUE)
  d$group <- ifelse(d$expr >= med, "High", "Low")
  fit <- survival::survfit(survival::Surv(OS.time, OS) ~ group, data = d)
  survminer::ggsurvplot(fit, data = d, pval = TRUE, risk.table = TRUE,
                        palette = palette, legend = "none",
                        legend.labs = c("High", "Low"), ggtheme = .theme_ubi(),
                        title = paste0(gene, " - KM survival curves"),
                        xlab = "Time (days)", ylab = "Survival probability")
}

#' Ranking bar chart of a score column
#'
#' @param scores data.frame from score_genes
#' @param value column to rank by, e.g. "Combined_Score"
#' @param top_n number of genes to show
#' @return a ggplot object
#' @export
plot_ranking <- function(scores, value = "Combined_Score", top_n = 20) {
  d <- scores[!is.na(scores[[value]]), , drop = FALSE]
  d <- d[order(d[[value]], decreasing = TRUE), , drop = FALSE]
  d <- head(d, top_n)
  if (nrow(d) == 0) return(.plot_na("No scores available"))
  d$Gene <- factor(d$Gene, levels = rev(d$Gene))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$Gene, y = .data[[value]])) +
    ggplot2::geom_col(fill = "#3B7DD8", alpha = 0.85, width = 0.7) +
    ggplot2::coord_flip() +
    .theme_ubi() +
    ggplot2::labs(title = paste0("Top ", nrow(d), " genes by ", value),
                  x = NULL, y = value) +
    ggplot2::theme(legend.position = "none")
}

#' Gene x cancer score heatmap
#'
#' @param scores_long long table from score_genes_multi
#' @param value column to plot
#' @return a ggplot object
#' @export
plot_score_heatmap <- function(scores_long, value = "Combined_Score") {
  d <- scores_long[!is.na(scores_long[[value]]), c("Gene", "Cancer", value), drop = FALSE]
  if (nrow(d) == 0) return(.plot_na("No scores available"))
  d$Gene <- factor(d$Gene, levels = sort(unique(d$Gene)))
  d$Cancer <- factor(d$Cancer, levels = sort(unique(d$Cancer)))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$Cancer, y = .data$Gene, fill = .data[[value]])) +
    ggplot2::geom_tile(color = "white", linewidth = 0.6) +
    ggplot2::scale_fill_gradient2(low = "#2E86C1", mid = "white", high = "#C0392B",
                                  midpoint = 0.5, limits = c(0, 1)) +
    .theme_ubi() +
    ggplot2::labs(title = paste0("Gene x cancer ", value, " heatmap"),
                  x = NULL, y = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Gene x cancer correlation heatmap (immune or metabolic)
#'
#' @param d long data.frame with columns Cancer, Trait, Cor
#' @param title plot title (图标题)
#' @param low low-end gradient color (低值渐变颜色)
#' @param high high-end gradient color (高值渐变颜色)
#' @param p_col column name of the p value used for significance stars (P 值列名)
#' @return a ggplot object
#' @export
plot_cor_heatmap <- function(d, title = "Correlation heatmap",
                             low = "#2E86C1", high = "#C0392B", p_col = "P") {
  d <- d[!is.na(d$Cor), , drop = FALSE]
  if (nrow(d) == 0) return(.plot_na("No correlation data available"))
  d$Cancer <- factor(d$Cancer, levels = sort(unique(d$Cancer)))
  d$Trait <- factor(d$Trait, levels = sort(unique(d$Trait)))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$Trait, y = .data$Cancer, fill = .data$Cor)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::scale_fill_gradient2(low = low, mid = "white", high = high,
                                  midpoint = 0, limits = c(-1, 1)) +
    .theme_ubi() +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  if (p_col %in% colnames(d) && any(!is.na(d[[p_col]]))) {
    d$star <- ifelse(is.na(d[[p_col]]), "",
              ifelse(d[[p_col]] < 0.001, "***",
              ifelse(d[[p_col]] < 0.01, "**",
              ifelse(d[[p_col]] < 0.05, "*", ""))))
    p <- p + ggplot2::geom_text(data = d, ggplot2::aes(label = .data$star), size = 3.2, color = "grey25")
  }
  p
}
