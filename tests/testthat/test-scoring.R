# scoring / plotting smoke tests (scoring and plotting smoke tests)
library(UbiPanTriage)

cache_file <- system.file("extdata", "LUAD_cache.rds", package = "UbiPanTriage")
skip_if_not(file.exists(cache_file), "LUAD cache not bundled")

luad <- load_cancer_data("LUAD")
genes <- c("MDM2", "TRIM44", "UBE2C", "FBXW7", "EGFR")

test_that("score_genes returns four-dimension columns and ranges", {
  res <- score_genes(luad, genes)
  sc <- res$scores
  expect_true(all(c("Gene", "Ubi_Score", "Basic_Score", "Immune_Score",
                    "Metabolic_Score", "Combined_Score", "Combined_Rank") %in% colnames(sc)))
  expect_true(all(sc$Ubi_Score >= 0 & sc$Ubi_Score <= 1, na.rm = TRUE))
  expect_true(all(sc$Basic_Score >= 0 & sc$Basic_Score <= 1))
  expect_true(all(sc$Immune_Score >= 0 & sc$Immune_Score <= 1))
  expect_true(all(sc$Metabolic_Score >= 0 & sc$Metabolic_Score <= 1))
  expect_true(all(sc$Combined_Score >= 0 & sc$Combined_Score <= 1))
  expect_equal(nrow(sc), length(genes))
  expect_true(all(c("Ubi_Type", "Ubi_Type_Full") %in% colnames(sc)))
  expect_true(any(sc$Ubi_Type == "E3"))
})

test_that("ubiquitin tier scores follow the confirmed hierarchy", {
  ann <- annotate_genes(luad, c("MDM2", "TRIM44", "UBE2C", "USP7", "RPS27A"))
  tiers <- ann$ubi_score[match(c("MDM2", "TRIM44", "UBE2C", "USP7", "RPS27A"), ann$gene)]
  expect_true(all(tiers[!is.na(tiers)] >= 0.6 & tiers[!is.na(tiers)] <= 1))
})

test_that("weight validation enforces 0.1-1.0 / one decimal / sum = 1", {
  expect_true(check_dim_weights(c(ubiquitin = .3, basic = .3, immune = .2, metabolic = .2))$ok)
  expect_false(check_dim_weights(c(ubiquitin = .5, basic = .5, immune = .5, metabolic = .5))$ok)
  expect_false(check_dim_weights(c(ubiquitin = .35, basic = .25, immune = .2, metabolic = .2))$ok)
  expect_false(check_dim_weights(c(ubiquitin = .8, basic = .1, immune = .05, metabolic = .05))$ok)
})

test_that("weight changes alter scores", {
  bp <- default_basic_params()
  bp$weights <- c(expression = 0, diff = 0, survive = 0, roc = 1)
  r1 <- score_genes(luad, "MDM2", basic_params = bp)$scores$Basic_Score
  r2 <- score_genes(luad, "MDM2")$scores$Basic_Score
  expect_false(isTRUE(all.equal(r1, r2)))
})

test_that("annotation is attached", {
  ann <- annotate_genes(luad, c("MDM2", "EGFR"))
  expect_true(all(c("family", "primary", "feature_score", "ubi_type", "ubi_type_full") %in% colnames(ann)))
  expect_true(any(grepl("E3", ann$primary[ann$gene == "MDM2"])))
})

test_that("cancers without normal tissue get Missing markers", {
  ov_file <- system.file("extdata", "OV_cache.rds", package = "UbiPanTriage")
  skip_if_not(file.exists(ov_file), "OV cache not bundled")
  ov <- load_cancer_data("OV")
  r <- score_genes(ov, "MDM2")
  expect_true(any(grepl("diff", r$scores$Missing)))
  expect_true(any(grepl("diag", r$scores$Missing)))
})

test_that("plots return ggplot objects", {
  skip_if_not(requireNamespace("ggplot2", quietly = TRUE))
  expect_s3_class(plot_radar(c(Basic = .8, Immune = .5, Metabolic = .6), "MDM2"), "ggplot")
  expect_s3_class(plot_immune_cor(luad, "MDM2", "timer"), "ggplot")
  expect_s3_class(plot_metabolic_cor(luad, "MDM2"), "ggplot")
  expect_s3_class(plot_expr_box(luad, "MDM2"), "ggplot")
  expect_s3_class(plot_roc_diag(luad, "MDM2"), "ggplot")
})

test_that("inference produces text in both languages", {
  res <- score_genes(luad, "MDM2")
  en <- gene_inference(luad, "MDM2", res, "en")
  zh <- gene_inference(luad, "MDM2", res, "zh")
  expect_type(en, "character")
  expect_type(zh, "character")
  expect_true(nchar(en) > 50)
  expect_true(nchar(zh) > 50)
})

test_that("multi-cancer scoring works", {
  rm <- score_genes_multi(list(LUAD = luad), c("MDM2", "TRIM44"))
  expect_true(all(c("Gene", "Cancer", "Ubi_Score", "Combined_Score") %in% colnames(rm$scores_long)))
  expect_true("Pan_Rank" %in% colnames(rm$summary))
})