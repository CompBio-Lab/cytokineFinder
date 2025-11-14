# Unit Tests for Results Processing Functions
# Functions tested:
#   - benchlist_to_tbl() from benchmark_list_to_tbl.R
#   - create_ensemble_results() from benchmark_with_ensemble.R
#   - extract_ligands() from utils_extract_ligands.R
#   - process_method_db() from utils_extract_ligands.R

library(testthat)
library(cytokineFinder)
library(tibble)
library(dplyr)

# ============================================================================
# Setup: Create test fixtures
# ============================================================================

create_test_benchmark_results <- function() {
  list(
    cfgsea = list(
      db1 = data.frame(
        ligand = c("IL6", "TNF", "IL10"),
        pval = c(0.01, 0.05, 0.1),
        padj = c(0.02, 0.1, 0.2)
      ),
      db2 = data.frame(
        ligand = c("IL6", "TNF"),
        pval = c(0.02, 0.03),
        padj = c(0.04, 0.06)
      )
    ),
    run_limma = list(
      db1 = data.frame(
        ligand = c("IL6", "TNF", "IL10"),
        coef = c(1.5, -0.8, 0.3),
        P.Value = c(0.001, 0.05, 0.3)
      )
    )
  )
}

create_test_nested_results <- function() {
  list(
    study1 = list(
      qc_eset = matrix(rnorm(100), nrow = 10),
      cond = rep(c("A", "B"), each = 5),
      benchmarks = create_test_benchmark_results()
    )
  )
}

create_test_cytosig_results <- function() {
  list(
    study1 = list(
      cfgsea = list(
        db1 = data.frame(
          ligand = c("IL6", "TNF", "IL10"),
          pval = c(0.02, 0.06, 0.15)
        )
      )
    )
  )
}

# ============================================================================
# Test benchlist_to_tbl() - Convert Results to Tibble
# ============================================================================

test_that("benchlist_to_tbl flattens nested benchmark results", {
  results <- create_test_nested_results()

  result_tbl <- benchlist_to_tbl(results, study_type = "test_study", has_benchmarks_layer = TRUE)

  expect_is(result_tbl, "data.frame")
  expect_true("study_type" %in% colnames(result_tbl))
  expect_true("method" %in% colnames(result_tbl))
  expect_true("database" %in% colnames(result_tbl))
  expect_true("ligand_tables" %in% colnames(result_tbl))
})

test_that("benchlist_to_tbl handles structure without benchmarks layer", {
  results <- create_test_cytosig_results()

  result_tbl <- benchlist_to_tbl(results, study_type = "test_study", has_benchmarks_layer = FALSE)

  expect_is(result_tbl, "data.frame")
  expect_true("study_type" %in% colnames(result_tbl))
})

test_that("benchlist_to_tbl preserves method and database names", {
  results <- create_test_nested_results()

  result_tbl <- benchlist_to_tbl(results, study_type = "test_study", has_benchmarks_layer = TRUE)

  expect_true("cfgsea" %in% result_tbl$method)
  expect_true("run_limma" %in% result_tbl$method)
  expect_true("db1" %in% result_tbl$database)
})

test_that("benchlist_to_tbl sets correct class for LRI methods", {
  results <- create_test_nested_results()

  result_tbl <- benchlist_to_tbl(results, study_type = "test_study", has_benchmarks_layer = TRUE)

  expect_true(all(result_tbl$class == "LRI"))
})

test_that("benchlist_to_tbl sets correct class for CytoSig methods", {
  results <- create_test_cytosig_results()

  result_tbl <- benchlist_to_tbl(results, study_type = "test_study", has_benchmarks_layer = FALSE)

  expect_true(all(result_tbl$class == "CytoSig_Web"))
})

test_that("benchlist_to_tbl ligand_tables contains data frames", {
  results <- create_test_nested_results()

  result_tbl <- benchlist_to_tbl(results, study_type = "test_study", has_benchmarks_layer = TRUE)

  for (i in seq_len(nrow(result_tbl))) {
    expect_is(result_tbl$ligand_tables[[i]], "data.frame")
  }
})

# ============================================================================
# Test extract_ligands() - Extract and Rank Ligands
# ============================================================================

test_that("extract_ligands extracts specified ligands", {
  benchmark_results <- create_test_benchmark_results()

  result <- extract_ligands(
    benchmark_results,
    ligands = c("IL6", "TNF"),
    metrics = c("padj", "coef")
  )

  expect_is(result, "data.frame")
  expect_true(all(result$ligand %in% c("IL6", "TNF")))
})

test_that("extract_ligands returns data frame with required columns", {
  benchmark_results <- create_test_benchmark_results()

  result <- extract_ligands(
    benchmark_results,
    ligands = c("IL6"),
    metrics = c("padj")
  )

  expect_true("ligand" %in% colnames(result))
  expect_true("value" %in% colnames(result))
  expect_true("method" %in% colnames(result))
  expect_true("database" %in% colnames(result))
  expect_true("rank" %in% colnames(result))
  expect_true("metric_type" %in% colnames(result))
})

test_that("extract_ligands handles coefficient metric correctly", {
  benchmark_results <- create_test_benchmark_results()

  result <- extract_ligands(
    benchmark_results,
    ligands = c("IL6"),
    metrics = c("coef")
  )

  expect_true(all(!is.na(result$value)))
})

test_that("extract_ligands handles p-value metric correctly", {
  benchmark_results <- create_test_benchmark_results()

  result <- extract_ligands(
    benchmark_results,
    ligands = c("IL6", "TNF"),
    metrics = c("padj")
  )

  expect_true(all(result$value >= 0 & result$value <= 1, na.rm = TRUE))
})

test_that("extract_ligands processes all method/database combinations", {
  benchmark_results <- create_test_benchmark_results()

  result <- extract_ligands(
    benchmark_results,
    ligands = c("IL6", "TNF"),
    metrics = c("coef", "padj")
  )

  expect_is(result, "data.frame")
  expect_gt(nrow(result), 0)
})

# ============================================================================
# Test process_method_db() - Process Method/Database Combination
# ============================================================================

test_that("process_method_db returns data frame with rank column", {
  test_df <- data.frame(
    ligand = c("IL6", "TNF", "IL10"),
    padj = c(0.01, 0.05, 0.1)
  )

  result <- process_method_db(test_df, "test_method", "test_db", metrics = "padj")

  expect_is(result, "data.frame")
  expect_true("rank" %in% colnames(result))
})

test_that("process_method_db calculates correct ranks for p-values", {
  test_df <- data.frame(
    ligand = c("IL6", "TNF", "IL10"),
    padj = c(0.01, 0.05, 0.1)
  )

  result <- process_method_db(test_df, "test_method", "test_db", metrics = "padj")

  expect_true(all(result$rank >= 0 & result$rank <= 100))
})

test_that("process_method_db ranks coefficients in descending order", {
  test_df <- data.frame(
    ligand = c("IL6", "TNF", "IL10"),
    coef = c(1.5, 0.5, 2.0)
  )

  result <- process_method_db(test_df, "test_method", "test_db", metrics = "coef")

  ranked <- result[order(result$rank, decreasing = TRUE), ]
  expect_equal(ranked$ligand[1], "IL10")
})

test_that("process_method_db adds method and database names", {
  test_df <- data.frame(
    ligand = c("IL6", "TNF"),
    padj = c(0.01, 0.05)
  )

  result <- process_method_db(test_df, "cfgsea", "db1", metrics = "padj")

  expect_true(all(result$method == "cfgsea"))
  expect_true(all(result$database == "db1"))
})

# ============================================================================
# Test create_ensemble_results() - Rank-Based Ensemble
# ============================================================================

test_that("create_ensemble_results filters for LRI methods", {
  results <- create_test_nested_results()
  result_tbl <- benchlist_to_tbl(results, study_type = "test_study", has_benchmarks_layer = TRUE)

  ensemble <- create_ensemble_results(result_tbl)

  expect_true(all(ensemble$class == "LRI"))
})

test_that("create_ensemble_results returns tibble with ensemble columns", {
  results <- create_test_nested_results()
  result_tbl <- benchlist_to_tbl(results, study_type = "test_study", has_benchmarks_layer = TRUE)

  ensemble <- create_ensemble_results(result_tbl)

  expect_true("ensemble_rank" %in% colnames(ensemble))
  expect_true("lri_rank" %in% colnames(ensemble))
})

test_that("create_ensemble_results computes numeric ranks", {
  results <- create_test_nested_results()
  result_tbl <- benchlist_to_tbl(results, study_type = "test_study", has_benchmarks_layer = TRUE)

  ensemble <- create_ensemble_results(result_tbl)

  expect_true(all(is.numeric(ensemble$lri_rank) | is.na(ensemble$lri_rank)))
})

test_that("create_ensemble_results handles missing CytoSig data", {
  results <- create_test_nested_results()
  result_tbl <- benchlist_to_tbl(results, study_type = "test_study", has_benchmarks_layer = TRUE)

  ensemble <- create_ensemble_results(result_tbl)

  expect_is(ensemble, "data.frame")
})

test_that("create_ensemble_results computes overlap_count correctly", {
  results <- create_test_nested_results()
  result_tbl <- benchlist_to_tbl(results, study_type = "test_study", has_benchmarks_layer = TRUE)

  ensemble <- create_ensemble_results(result_tbl)

  expect_true("overlap_count" %in% colnames(ensemble))
  expect_true(all(is.numeric(ensemble$overlap_count)))
})
