# ============================================================================
# Shiny app regression tests (5 pages)
# Loads the bundled app.R exactly as shiny::runApp does (UTF-8, locale-safe),
# then exercises all four scoring pages with testServer().
# NOTE: `invisible(eval(...))` is required -- app.R ends with `shinyApp(ui, server)`
# and an auto-printed app object would launch a real server via print.shiny.appobj.
# ============================================================================
skip_if_not_installed("shiny")
suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
})

# --- load the app exactly as shiny::runApp does (UTF-8, locale-independent) ----
app_file <- system.file("shiny", "app.R", package = "UbiPanTriage")
skip_if(!file.exists(app_file), "shiny app not bundled in installed package")
env <- new.env(parent = globalenv())
lines <- readLines(app_file, encoding = "UTF-8", warn = FALSE)
invisible(eval(parse(text = lines, encoding = "UTF-8"), envir = env))
server <- env$server

# byte-level substring search (robust under a C locale)
has_bytes <- function(x, target) {
  s <- charToRaw(paste(as.character(x), collapse = ""))
  t <- charToRaw(target)
  if (length(t) == 0 || length(s) < length(t)) return(FALSE)
  any(vapply(seq_len(length(s) - length(t) + 1),
             function(i) identical(s[i:(i + length(t) - 1)], t), logical(1)))
}

# sub-item weight inputs (TIMER cells + 7 metabolic pathways), value = 1 each
sub_weights <- function(prefix) {
  stats::setNames(rep(1, 17),
                  c(paste0(prefix, "_cell_", c("B_cell", "T_cell_CD4", "T_cell_CD8",
                                               "Neutrophil", "Macrophage", "DC")),
                    paste0(prefix, "_meta_", c("Amino_acid", "Carbohydrate", "Energy",
                                               "Lipid", "Nucleotide", "TCA_cycle",
                                               "Vitamin_cofactor")),
                    paste0(prefix, "_bw_expression"), paste0(prefix, "_bw_diff"),
                    paste0(prefix, "_bw_survive"), paste0(prefix, "_bw_roc")))
}

test_that("main UI renders in English and Chinese", {
  testServer(server, {
    expect_false(is.null(output$app_ui))
    expect_true(has_bytes(output$app_ui, "Ubiquitin gene-set ranking"))
    session$setInputs(lang = "zh")
    session$flushReact()
    expect_false(is.null(output$app_ui))
    expect_true(has_bytes(output$app_ui, "泛素基因集评分排序"))
    session$setInputs(lang = "en")
    session$flushReact()
    expect_false(is.null(output$app_ui))
  })
})

test_that("page 1 ubiquitin gene-set ranking computes and renders", {
  testServer(server, {
    do.call(session$setInputs, c(list(
      p1_genes = c("MDM2", "USP7", "TRIM44"),
      p1_cancers = c("LUAD", "BRCA", "LAML"),
      p1_dims = c("ubiquitin", "basic", "immune", "metabolic"),
      p1_imm = "timer", p1_agg = "equal", p1_wmode = "equal"), as.list(sub_weights("p1"))))
    session$setInputs(p1_run = 1)
    session$flushReact()
    r <- p1_res()
    expect_false(is.null(r))
    expect_equal(nrow(r$summ), 3L)
    expect_true(all(c("Gene", "Ubi_Type_Full", "Combined_Score", "Rank", "Pct") %in%
                    colnames(r$summ)))
    # LAML (no normal tissue) must be flagged as missing for at least one gene
    expect_true(any(r$summ$n_missing > 0))
    expect_false(is.null(output$p1_notice))
    expect_false(is.null(output$p1_table))
    expect_false(is.null(output$p1_4d))
    expect_false(is.null(output$p1_facet))
    expect_false(is.null(output$p1_ranking))
  })
})

test_that("page 2 ubiquitin single-gene query renders all outputs", {
  testServer(server, {
    session$setInputs(p2_gene = "MDM2", p2_cancers = c("LUAD", "BRCA"),
                      p2_imm = "timer",
                      p2_dims = c("ubiquitin", "basic", "immune", "metabolic"),
                      p2_plots = c("radar", "imm", "meta", "expr", "box", "roc",
                                   "km", "upset"),
                      p2_expr_unit = "tpm", p2_expr_normal = "both",
                      p2_thr_mode = "median", p2_thr = 0.5)
    session$setInputs(p2_run = 1)
    session$flushReact()
    r <- p2_res()
    expect_false(is.null(r))
    expect_equal(nrow(r$d), 2L)
    expect_true("Combined_Rank" %in% colnames(r$d))
    expect_false(is.null(output$p2_scores))
    expect_false(is.null(output$p2_radar))
    expect_false(is.null(output$p2_imm_heat))
    expect_false(is.null(output$p2_meta_heat))
    expect_false(is.null(output$p2_expr))
    expect_false(is.null(output$p2_box))
    expect_false(is.null(output$p2_roc))
    expect_false(is.null(output$p2_km))
    expect_false(is.null(output$p2_upset))
    expect_true(is.character(output$p2_infer_txt))
    expect_true(nzchar(output$p2_infer_txt))
  })
})

test_that("page 3 custom gene set computes and renders", {
  testServer(server, {
    do.call(session$setInputs, c(list(
      p3_mode = "paste", p3_paste = "CDK1\nAURKA\nTOP2A",
      p3_has_score = FALSE, p3_cancers = c("LUAD", "BRCA"),
      p3_dims = c("basic", "immune", "metabolic"),
      p3_imm = "timer", p3_agg = "equal", p3_wmode = "equal"), as.list(sub_weights("p3"))))
    session$setInputs(p3_run = 1)
    session$flushReact()
    r <- p3_res()
    expect_false(is.null(r))
    expect_equal(nrow(r$summ), 3L)
    expect_false(is.null(output$p3_table))
    expect_false(is.null(output$p3_4d))
    expect_false(is.null(output$p3_facet))
    expect_false(is.null(output$p3_ranking))
  })
})

test_that("page 4 custom single-gene query renders all outputs", {
  testServer(server, {
    session$setInputs(p4_gene = "CDK1", p4_cancers = c("LUAD", "BRCA"),
                      p4_imm = "timer", p4_dims = c("basic", "immune", "metabolic"),
                      p4_plots = c("radar", "imm", "meta", "expr", "box", "roc",
                                   "km", "upset"),
                      p4_expr_unit = "tpm", p4_expr_normal = "both",
                      p4_thr_mode = "median", p4_thr = 0.5)
    session$setInputs(p4_run = 1)
    session$flushReact()
    r <- p4_res()
    expect_false(is.null(r))
    expect_equal(nrow(r$d), 2L)
    expect_false(is.null(output$p4_scores))
    expect_false(is.null(output$p4_radar))
    expect_false(is.null(output$p4_imm_heat))
    expect_false(is.null(output$p4_meta_heat))
    expect_false(is.null(output$p4_expr))
    expect_false(is.null(output$p4_box))
    expect_false(is.null(output$p4_roc))
    expect_false(is.null(output$p4_km))
    expect_false(is.null(output$p4_upset))
    expect_true(is.character(output$p4_infer_txt))
    expect_true(nzchar(output$p4_infer_txt))
  })
})
