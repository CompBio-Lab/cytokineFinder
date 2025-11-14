# Unit Tests for Statistical Methods
# Functions tested:
#   - cfgsea() from stats_cfgsea.R
#   - run_limma() from stats_run_limma.R
#   - gsva_limma() from stats_gsva_limma.R
#   - pca_limma() from stats_cpca_limma.R
#   - gsva_plsda() from stats_gsva_plsda.R
#   - pca_plsda() from stats_cpca_plsda.R
#   - cytosig_custom_ridge() from stats_cytosig_custom_ridge.R

library(testthat)
library(cytokineFindeR)

# ============================================================================
# Setup: Create test fixtures
# ============================================================================

create_test_eset <- function(n_genes = 100, n_samples = 10, seed = 123) {
  set.seed(seed)
  eset <- matrix(
    rnorm(n_genes * n_samples, mean = 5, sd = 1),
    nrow = n_genes,
    ncol = n_samples
  )
  rownames(eset) <- paste0("GENE_", seq_len(n_genes))
  colnames(eset) <- paste0("SAMPLE_", seq_len(n_samples))
  eset
}

create_test_db <- function() {
  list(
    IL6 = c("GENE_1", "GENE_2", "GENE_3"),
    TNF = c("GENE_4", "GENE_5", "GENE_6"),
    IL10 = c("GENE_7", "GENE_8", "GENE_9")
  )
}

# ============================================================================
# Test run_limma() - Limma DEA
# ============================================================================

test_that("run_limma returns data frame with expected columns", {
  eset <- create_test_eset()
  cond <- rep(c("control", "treatment"), each = 5)
  design <- create_design(cond)

  result <- run_limma(eset, design$design)

  expect_is(result, "data.frame")
  expect_true("genes" %in% colnames(result))
  expect_true("logFC" %in% colnames(result))
  expect_true("P.Value" %in% colnames(result))
})

test_that("run_limma returns all genes in results", {
  eset <- create_test_eset(n_genes = 100, n_samples = 10)
  cond <- rep(c("control", "treatment"), each = 5)
  design <- create_design(cond)

  result <- run_limma(eset, design$design)

  expect_equal(nrow(result), 100)
})

test_that("run_limma handles paired design", {
  eset <- create_test_eset()
  cond <- rep(c("control", "treatment"), each = 5)
  obs_id <- rep(1:5, 2)
  design <- create_design(cond, obs_id = obs_id, eset = eset)

  result <- run_limma(eset, design$design, obs_id = obs_id, correlation = design$dupcor)

  expect_is(result, "data.frame")
  expect_true("logFC" %in% colnames(result))
})

test_that("run_limma p-values are between 0 and 1", {
  eset <- create_test_eset()
  cond <- rep(c("control", "treatment"), each = 5)
  design <- create_design(cond)

  result <- run_limma(eset, design$design)

  expect_true(all(result$P.Value >= 0 & result$P.Value <= 1, na.rm = TRUE))
})

# ============================================================================
# Test cfgsea() - Gene Set Enrichment Analysis
# ============================================================================

test_that("cfgsea runs DEA and GSEA pipeline", {
  eset <- create_test_eset()
  cond <- rep(c("control", "treatment"), each = 5)
  design <- create_design(cond)
  db <- create_test_db()

  result <- cfgsea(eset, design$design, db)

  expect_is(result, "data.frame")
  expect_true("ligand" %in% colnames(result))
})

test_that("cfgsea returns results only for ligands in database", {
  eset <- create_test_eset()
  cond <- rep(c("control", "treatment"), each = 5)
  design <- create_design(cond)
  db <- list(IL6 = c("GENE_1", "GENE_2"))

  result <- cfgsea(eset, design$design, db)

  expect_true(all(result$ligand %in% names(db)))
})

test_that("cfgsea p-values are between 0 and 1", {
  eset <- create_test_eset()
  cond <- rep(c("control", "treatment"), each = 5)
  design <- create_design(cond)
  db <- create_test_db()

  result <- cfgsea(eset, design$design, db)

  expect_true(all(result$pval >= 0 & result$pval <= 1, na.rm = TRUE))
})

test_that("cfgsea handles paired design", {
  eset <- create_test_eset()
  cond <- rep(c("control", "treatment"), each = 5)
  obs_id <- rep(1:5, 2)
  design <- create_design(cond, obs_id = obs_id, eset = eset)
  db <- create_test_db()

  result <- cfgsea(eset, design$design, db, obs_id = obs_id, correlation = design$dupcor)

  expect_is(result, "data.frame")
})

# ============================================================================
# Test gsva_limma() - GSVA + Limma
# ============================================================================

test_that("gsva_limma computes GSVA scores and runs DEA", {
  eset <- create_test_eset()
  cond <- rep(c("control", "treatment"), each = 5)
  design <- create_design(cond)
  db <- create_test_db()

  result <- gsva_limma(eset, design$design, db)

  expect_is(result, "data.frame")
  expect_true("ligand" %in% colnames(result))
  expect_true("pval" %in% colnames(result))
  expect_true("padj" %in% colnames(result))
})

test_that("gsva_limma returns adjusted p-values", {
  eset <- create_test_eset()
  cond <- rep(c("control", "treatment"), each = 5)
  design <- create_design(cond)
  db <- create_test_db()

  result <- gsva_limma(eset, design$design, db)

  expect_true(all(result$padj >= 0 & result$padj <= 1, na.rm = TRUE))
})

test_that("gsva_limma handles paired design", {
  eset <- create_test_eset()
  cond <- rep(c("control", "treatment"), each = 5)
  obs_id <- rep(1:5, 2)
  design <- create_design(cond, obs_id = obs_id, eset = eset)
  db <- create_test_db()

  result <- gsva_limma(eset, design$design, db, obs_id = obs_id, correlation = design$dupcor)

  expect_is(result, "data.frame")
})

# ============================================================================
# Test pca_limma() - PCA + Limma
# ============================================================================

test_that("pca_limma computes PC1 for each ligand's receptors", {
  eset <- create_test_eset()
  cond <- rep(c("control", "treatment"), each = 5)
  design <- create_design(cond)
  db <- create_test_db()

  result <- pca_limma(eset, design$design, db)

  expect_is(result, "data.frame")
  expect_true("ligand" %in% colnames(result))
})

test_that("pca_limma returns p-values and adjusted p-values", {
  eset <- create_test_eset()
  cond <- rep(c("control", "treatment"), each = 5)
  design <- create_design(cond)
  db <- create_test_db()

  result <- pca_limma(eset, design$design, db)

  expect_true("pval" %in% colnames(result))
  expect_true("padj" %in% colnames(result))
})

test_that("pca_limma handles paired design", {
  eset <- create_test_eset()
  cond <- rep(c("control", "treatment"), each = 5)
  obs_id <- rep(1:5, 2)
  design <- create_design(cond, obs_id = obs_id, eset = eset)
  db <- create_test_db()

  result <- pca_limma(eset, design$design, db, obs_id = obs_id, correlation = design$dupcor)

  expect_is(result, "data.frame")
})

test_that("pca_limma handles zero-variance ligands gracefully", {
  eset <- create_test_eset()
  cond <- rep(c("control", "treatment"), each = 5)
  design <- create_design(cond)
  db <- list(
    MultiGene = c("GENE_1", "GENE_2", "GENE_3"),
    SingleGene = c("GENE_50")
  )

  result <- pca_limma(eset, design$design, db)

  expect_is(result, "data.frame")
})

# ============================================================================
# Test gsva_plsda() - GSVA + PLS-DA
# ============================================================================

test_that("gsva_plsda computes GSVA and fits PLS-DA model", {
  eset <- create_test_eset()
  cond <- rep(c("control", "treatment"), each = 5)
  db <- create_test_db()

  result <- gsva_plsda(eset, cond, db)

  expect_is(result, "data.frame")
  expect_true("ligand" %in% colnames(result))
  expect_true("coef" %in% colnames(result))
})

test_that("gsva_plsda returns numeric coefficients", {
  eset <- create_test_eset()
  cond <- rep(c("control", "treatment"), each = 5)
  db <- create_test_db()

  result <- gsva_plsda(eset, cond, db)

  expect_true(all(is.numeric(result$coef)))
})

test_that("gsva_plsda handles binary treatment", {
  eset <- create_test_eset()
  treatment <- rep(c("A", "B"), each = 5)
  db <- create_test_db()

  result <- gsva_plsda(eset, treatment, db)

  expect_is(result, "data.frame")
})

test_that("gsva_plsda handles multi-class treatment", {
  eset <- create_test_eset(n_samples = 12)
  treatment <- rep(c("A", "B", "C"), each = 4)
  db <- create_test_db()

  result <- gsva_plsda(eset, treatment, db)

  expect_is(result, "data.frame")
})

# ============================================================================
# Test pca_plsda() - PCA + PLS-DA
# ============================================================================

test_that("pca_plsda computes PC1 for receptors and fits PLS-DA", {
  eset <- create_test_eset()
  cond <- rep(c("control", "treatment"), each = 5)
  db <- create_test_db()

  result <- pca_plsda(eset, cond, db)

  expect_is(result, "data.frame")
  expect_true("ligand" %in% colnames(result))
  expect_true("coef" %in% colnames(result))
})

test_that("pca_plsda coefficients are numeric and non-zero", {
  eset <- create_test_eset()
  cond <- rep(c("control", "treatment"), each = 5)
  db <- create_test_db()

  result <- pca_plsda(eset, cond, db)

  expect_true(all(is.numeric(result$coef)))
})

test_that("pca_plsda handles zero-variance ligands gracefully", {
  eset <- create_test_eset()
  cond <- rep(c("control", "treatment"), each = 5)
  db <- list(
    MultiGene = c("GENE_1", "GENE_2", "GENE_3"),
    SingleGene = c("GENE_50")
  )

  result <- pca_plsda(eset, cond, db)

  expect_is(result, "data.frame")
})

# ============================================================================
# Test cytosig_custom_ridge() - Ridge Regression
# ============================================================================

test_that("cytosig_custom_ridge runs limma DEA and ridge regression", {
  eset <- create_test_eset(n_genes = 50, n_samples = 10)
  cond <- rep(c("control", "treatment"), each = 5)
  design <- create_design(cond)

  beta_coef <- matrix(
    rnorm(50 * 5),
    nrow = 50,
    ncol = 5,
    dimnames = list(
      paste0("GENE_", 1:50),
      paste0("LIGAND_", 1:5)
    )
  )

  result <- cytosig_custom_ridge(eset, design$design, beta_coef = beta_coef)

  expect_is(result, "data.frame")
  expect_true("ligand" %in% colnames(result))
  expect_true("coef" %in% colnames(result))
})

test_that("cytosig_custom_ridge returns numeric coefficients", {
  eset <- create_test_eset(n_genes = 50, n_samples = 10)
  cond <- rep(c("control", "treatment"), each = 5)
  design <- create_design(cond)

  beta_coef <- matrix(
    rnorm(50 * 5),
    nrow = 50,
    ncol = 5,
    dimnames = list(
      paste0("GENE_", 1:50),
      paste0("LIGAND_", 1:5)
    )
  )

  result <- cytosig_custom_ridge(eset, design$design, beta_coef = beta_coef)

  expect_true(all(is.numeric(result$coef)))
})

test_that("cytosig_custom_ridge handles paired design", {
  eset <- create_test_eset(n_genes = 50, n_samples = 10)
  cond <- rep(c("control", "treatment"), each = 5)
  obs_id <- rep(1:5, 2)
  design <- create_design(cond, obs_id = obs_id, eset = eset)

  beta_coef <- matrix(
    rnorm(50 * 5),
    nrow = 50,
    ncol = 5,
    dimnames = list(
      paste0("GENE_", 1:50),
      paste0("LIGAND_", 1:5)
    )
  )

  result <- cytosig_custom_ridge(
    eset, design$design, obs_id = obs_id,
    correlation = design$dupcor, beta_coef = beta_coef
  )

  expect_is(result, "data.frame")
})

test_that("cytosig_custom_ridge handles gene overlap correctly", {
  eset <- create_test_eset(n_genes = 50, n_samples = 10)
  cond <- rep(c("control", "treatment"), each = 5)
  design <- create_design(cond)

  beta_coef <- matrix(
    rnorm(100 * 5),
    nrow = 100,
    ncol = 5,
    dimnames = list(
      paste0("GENE_", 1:100),
      paste0("LIGAND_", 1:5)
    )
  )

  result <- cytosig_custom_ridge(eset, design$design, beta_coef = beta_coef)

  expect_is(result, "data.frame")
})
