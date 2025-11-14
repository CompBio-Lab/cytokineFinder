# Unit Tests for run_cytokinefinder() - Main Workflow Function
# Location: R/run_cytokinefinder_withList.R

library(testthat)
library(cytokineFinder)

# ============================================================================
# Setup: Create test fixtures and helper functions
# ============================================================================

# Helper: Create minimal valid test eset
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

# Helper: Create minimal LRI database
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

# Helper: Create study_data list
create_study_data <- function(eset, condition, obs_id = NULL) {
  study_data <- list(
    qc_eset = eset,
    cond = condition
  )
  if (!is.null(obs_id)) {
    study_data$obs_id <- obs_id
  }
  study_data
}

# ============================================================================
# Test 1.1.1: Input Validation Tests
# ============================================================================

test_that("run_cytokinefinder validates required qc_eset parameter", {
  # Setup: Missing qc_eset
  study_data <- list(cond = rep(c("control", "treatment"), each = 5))
  databases <- create_test_db()
  methods <- c("cfgsea")

  # Execute & Verify: Should throw error
  expect_error(
    run_cytokinefinder(study_data, databases, methods),
    "Missing required elements"
  )
  expect_error(
    run_cytokinefinder(study_data, databases, methods),
    "qc_eset"
  )
})

test_that("run_cytokinefinder validates required cond parameter", {
  # Setup: Missing cond
  eset <- create_test_eset(n_genes = 100, n_samples = 10)
  study_data <- list(qc_eset = eset)
  databases <- create_test_db()
  methods <- c("cfgsea")

  # Execute & Verify: Should throw error
  expect_error(
    run_cytokinefinder(study_data, databases, methods),
    "Missing required elements"
  )
  expect_error(
    run_cytokinefinder(study_data, databases, methods),
    "cond"
  )
})

test_that("run_cytokinefinder validates cond length matches sample count", {
  # Setup: cond length != ncol(eset)
  eset <- create_test_eset(n_genes = 100, n_samples = 10)
  cond <- rep(c("control", "treatment"), each = 3)  # Length 6, not 10
  study_data <- create_study_data(eset, cond)
  databases <- create_test_db()
  methods <- c("cfgsea")

  # Execute & Verify: Should throw error
  expect_error(
    run_cytokinefinder(study_data, databases, methods),
    "non-conformable"  # model.matrix error
  )
})

# ============================================================================
# Test 1.1.2: Output Structure Validation Tests
# ============================================================================

test_that("run_cytokinefinder returns correctly structured output", {
  # Setup: Create minimal valid inputs
  eset <- create_test_eset(n_genes = 100, n_samples = 10)
  cond <- rep(c("control", "treatment"), each = 5)
  study_data <- create_study_data(eset, cond)
  databases <- create_test_db()
  methods <- c("cfgsea")

  # Execute
  result <- run_cytokinefinder(study_data, databases, methods)

  # Verify: Output contains expected elements
  expect_true(is.list(result))
  expect_true("qc_eset" %in% names(result))
  expect_true("cond" %in% names(result))
  expect_true("design" %in% names(result))
  expect_true("benchmarks" %in% names(result))
})

test_that("run_cytokinefinder preserves original qc_eset", {
  # Setup
  eset <- create_test_eset(n_genes = 100, n_samples = 10)
  cond <- rep(c("control", "treatment"), each = 5)
  study_data <- create_study_data(eset, cond)
  databases <- create_test_db()
  methods <- c("cfgsea")

  # Execute
  result <- run_cytokinefinder(study_data, databases, methods)

  # Verify: Original eset is preserved
  expect_identical(result$qc_eset, eset)
  expect_equal(ncol(result$qc_eset), ncol(eset))
})

test_that("run_cytokinefinder output design is a list with design matrix", {
  # Setup
  eset <- create_test_eset(n_genes = 100, n_samples = 10)
  cond <- rep(c("control", "treatment"), each = 5)
  study_data <- create_study_data(eset, cond)
  databases <- create_test_db()
  methods <- c("cfgsea")

  # Execute
  result <- run_cytokinefinder(study_data, databases, methods)

  # Verify: design is a list with design matrix
  expect_is(result$design, "list")
  expect_true("design" %in% names(result$design))
  expect_is(result$design$design, "matrix")
})

test_that("run_cytokinefinder benchmarks is BenchmarkResults class", {
  # Setup
  eset <- create_test_eset(n_genes = 100, n_samples = 10)
  cond <- rep(c("control", "treatment"), each = 5)
  study_data <- create_study_data(eset, cond)
  databases <- create_test_db()
  methods <- c("cfgsea")

  # Execute
  result <- run_cytokinefinder(study_data, databases, methods)

  # Verify: benchmarks is BenchmarkResults object
  expect_is(result$benchmarks, "BenchmarkResults")
  expect_is(result$benchmarks, "list")
})

# ============================================================================
# Test 1.1.3: Paired vs Unpaired Design Handling
# ============================================================================

test_that("run_cytokinefinder handles unpaired design correctly", {
  # Setup: No obs_id (unpaired)
  eset <- create_test_eset(n_genes = 100, n_samples = 10)
  cond <- rep(c("control", "treatment"), each = 5)
  study_data <- create_study_data(eset, cond, obs_id = NULL)
  databases <- create_test_db()
  methods <- c("cfgsea")

  # Execute
  result <- run_cytokinefinder(study_data, databases, methods)

  # Verify: design matrix exists
  expect_is(result$design$design, "matrix")
  expect_equal(nrow(result$design$design), 10)  # Same as samples

  # Verify: dupcor is NULL for unpaired
  expect_null(result$design$dupcor)
})

test_that("run_cytokinefinder handles paired design correctly", {
  # Setup: With obs_id (paired)
  eset <- create_test_eset(n_genes = 100, n_samples = 10)
  cond <- rep(c("control", "treatment"), each = 5)
  obs_id <- rep(1:5, 2)  # 5 pairs
  study_data <- create_study_data(eset, cond, obs_id = obs_id)
  databases <- create_test_db()
  methods <- c("cfgsea")

  # Execute
  result <- run_cytokinefinder(study_data, databases, methods)

  # Verify: design matrix exists
  expect_is(result$design$design, "matrix")

  # Verify: dupcor is calculated for paired
  expect_false(is.null(result$design$dupcor))
  expect_is(result$design$dupcor, "list")
  expect_true("consensus.correlation" %in% names(result$design$dupcor))
})

# ============================================================================
# Test 1.1.4: Database Filtering Integration
# ============================================================================

test_that("run_cytokinefinder filters databases to available genes", {
  # Setup: eset has limited genes
  eset <- create_test_eset(n_genes = 50, n_samples = 10)  # Genes 1-50
  cond <- rep(c("control", "treatment"), each = 5)
  study_data <- create_study_data(eset, cond)

  # Database with genes beyond what's in eset
  databases <- list(
    database1 = list(
      IL6 = c("GENE_1", "GENE_2", "GENE_3", "GENE_1000"),  # GENE_1000 not in eset
      TNF = c("GENE_4", "GENE_5", "GENE_6")
    )
  )
  methods <- c("cfgsea")

  # Execute
  result <- run_cytokinefinder(study_data, databases, methods)

  # Verify: Results don't error (filtering happened)
  expect_is(result$benchmarks, "BenchmarkResults")
})

# ============================================================================
# Integration Test: Full Workflow
# ============================================================================

test_that("run_cytokinefinder completes full workflow without error", {
  # Setup: Full realistic scenario
  eset <- create_test_eset(n_genes = 150, n_samples = 12)
  cond <- rep(c("control", "treatment"), each = 6)
  study_data <- create_study_data(eset, cond)
  databases <- create_test_db()
  methods <- c("cfgsea")

  # Execute
  expect_no_error(
    result <- run_cytokinefinder(study_data, databases, methods)
  )

  # Verify: All expected components present
  expect_true(all(c("qc_eset", "cond", "design", "benchmarks") %in% names(result)))
})

test_that("run_cytokinefinder preserves study_data original elements", {
  # Setup
  eset <- create_test_eset(n_genes = 100, n_samples = 10)
  cond <- rep(c("control", "treatment"), each = 5)
  study_data <- create_study_data(eset, cond)
  study_data$custom_element <- "my_value"  # Add custom element

  databases <- create_test_db()
  methods <- c("cfgsea")

  # Execute
  result <- run_cytokinefinder(study_data, databases, methods)

  # Verify: Custom element preserved
  expect_true("custom_element" %in% names(result))
  expect_equal(result$custom_element, "my_value")
})
