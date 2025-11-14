# Unit Tests for run_lri_methods() - Parallel Execution Engine
# Location: R/run_lri_methods.R

library(testthat)
library(cytokineFindeR)

# ============================================================================
# Setup: Create test fixtures and helper functions
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
    database1 = list(
      IL6 = c("GENE_1", "GENE_2", "GENE_3"),
      TNF = c("GENE_4", "GENE_5", "GENE_6"),
      IL10 = c("GENE_7", "GENE_8", "GENE_9")
    ),
    database2 = list(
      IL6 = c("GENE_2", "GENE_3", "GENE_10"),
      TNF = c("GENE_5", "GENE_6", "GENE_11")
    )
  )
}

# ============================================================================
# Test 1.2.1: Method Execution
# ============================================================================

test_that("run_lri_methods executes single method on single database", {
  eset <- create_test_eset(n_genes = 100, n_samples = 10)
  cond <- rep(c("control", "treatment"), each = 5)
  design <- create_design(cond)
  databases <- list(
    db1 = list(IL6 = c("GENE_1", "GENE_2", "GENE_3"))
  )
  methods <- c("cfgsea")

  result <- run_lri_methods(eset, design$design, databases, methods)

  expect_is(result, "BenchmarkResults")
  expect_true("cfgsea" %in% names(result))
  expect_true("db1" %in% names(result$cfgsea))
})

test_that("run_lri_methods executes multiple methods", {
  eset <- create_test_eset(n_genes = 100, n_samples = 10)
  cond <- rep(c("control", "treatment"), each = 5)
  design <- create_design(cond)
  databases <- list(
    db1 = list(IL6 = c("GENE_1", "GENE_2", "GENE_3"))
  )
  methods <- c("cfgsea", "run_limma")

  result <- run_lri_methods(eset, design$design, databases, methods)

  expect_is(result, "BenchmarkResults")
  expect_true("cfgsea" %in% names(result))
  expect_true("run_limma" %in% names(result))
})

test_that("run_lri_methods processes multiple databases", {
  eset <- create_test_eset(n_genes = 100, n_samples = 10)
  cond <- rep(c("control", "treatment"), each = 5)
  design <- create_design(cond)
  databases <- create_test_db()
  methods <- c("cfgsea")

  result <- run_lri_methods(eset, design$design, databases, methods)

  expect_true("database1" %in% names(result$cfgsea))
  expect_true("database2" %in% names(result$cfgsea))
})

# ============================================================================
# Test 1.2.2: PLSDA Method Handling
# ============================================================================

test_that("run_lri_methods applies PLSDA-specific parameters", {
  eset <- create_test_eset(n_genes = 100, n_samples = 10)
  cond <- rep(c("control", "treatment"), each = 5)
  design <- create_design(cond)
  databases <- list(
    db1 = list(IL6 = c("GENE_1", "GENE_2", "GENE_3"))
  )
  methods <- c("gsva_plsda")

  result <- run_lri_methods(
    eset, design$design, databases, methods,
    treatment = cond
  )

  expect_is(result, "BenchmarkResults")
  expect_true("gsva_plsda" %in% names(result))
})

# ============================================================================
# Test 1.2.3: Paired Design Support
# ============================================================================

test_that("run_lri_methods handles unpaired design", {
  eset <- create_test_eset(n_genes = 100, n_samples = 10)
  cond <- rep(c("control", "treatment"), each = 5)
  design <- create_design(cond)
  databases <- list(
    db1 = list(IL6 = c("GENE_1", "GENE_2", "GENE_3"))
  )
  methods <- c("cfgsea")

  result <- run_lri_methods(eset, design$design, databases, methods)

  expect_is(result, "BenchmarkResults")
  expect_is(result$cfgsea$db1, "data.frame")
})

test_that("run_lri_methods handles paired design with dupcor", {
  eset <- create_test_eset(n_genes = 100, n_samples = 10)
  cond <- rep(c("control", "treatment"), each = 5)
  obs_id <- rep(1:5, 2)
  design <- create_design(cond, obs_id = obs_id, eset = eset)

  databases <- list(
    db1 = list(IL6 = c("GENE_1", "GENE_2", "GENE_3"))
  )
  methods <- c("cfgsea")

  result <- run_lri_methods(
    eset, design$design, databases, methods,
    obs_id = obs_id,
    correlation = design$dupcor
  )

  expect_is(result, "BenchmarkResults")
})

# ============================================================================
# Test 1.2.4: Output Structure
# ============================================================================

test_that("run_lri_methods returns BenchmarkResults class", {
  eset <- create_test_eset(n_genes = 100, n_samples = 10)
  cond <- rep(c("control", "treatment"), each = 5)
  design <- create_design(cond)
  databases <- list(
    db1 = list(IL6 = c("GENE_1", "GENE_2", "GENE_3"))
  )
  methods <- c("cfgsea")

  result <- run_lri_methods(eset, design$design, databases, methods)

  expect_s3_class(result, "BenchmarkResults")
})

test_that("run_lri_methods results are data frames", {
  eset <- create_test_eset(n_genes = 100, n_samples = 10)
  cond <- rep(c("control", "treatment"), each = 5)
  design <- create_design(cond)
  databases <- list(
    db1 = list(IL6 = c("GENE_1", "GENE_2", "GENE_3"))
  )
  methods <- c("cfgsea")

  result <- run_lri_methods(eset, design$design, databases, methods)

  expect_is(result$cfgsea$db1, "data.frame")
})
