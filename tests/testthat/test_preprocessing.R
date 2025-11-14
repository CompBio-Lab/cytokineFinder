# Unit Tests for Preprocessing Functions
# Functions tested:
#   - create_design() from utils_create_design.R
#   - preprocess_eset() from utils_preprocess_eset.R
#   - filter_db_against_eset() from utils_preprocess_eset.R
#   - remove_zero_variance_ligands() from utils_preprocess_eset.R
#   - clean_eset() from utils_clean_data.R
#   - create_db_space() from utils_create_db_space.R

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
# Test create_design() - Design Matrix Creation
# ============================================================================

test_that("create_design generates design matrix from treatment vector", {
  treatment <- rep(c("control", "treatment"), each = 5)
  result <- create_design(treatment)

  expect_is(result, "list")
  expect_true("design" %in% names(result))
  expect_is(result$design, "matrix")
})

test_that("create_design handles unpaired design (no obs_id)", {
  treatment <- rep(c("control", "treatment"), each = 5)
  result <- create_design(treatment, obs_id = NULL)

  expect_is(result$design, "matrix")
  expect_equal(nrow(result$design), 10)
  expect_null(result$dupcor)
})

test_that("create_design handles paired design with obs_id", {
  treatment <- rep(c("control", "treatment"), each = 5)
  obs_id <- rep(1:5, 2)
  eset <- create_test_eset()

  result <- create_design(treatment, obs_id = obs_id, eset = eset)

  expect_is(result$design, "matrix")
  expect_false(is.null(result$dupcor))
  expect_is(result$dupcor, "list")
})

test_that("create_design computes dupcor correctly for paired samples", {
  treatment <- rep(c("control", "treatment"), each = 5)
  obs_id <- rep(1:5, 2)
  eset <- create_test_eset()

  result <- create_design(treatment, obs_id = obs_id, eset = eset)

  expect_true("consensus.correlation" %in% names(result$dupcor))
  expect_is(result$dupcor$consensus.correlation, "numeric")
})

test_that("create_design handles multiple treatment levels", {
  treatment <- c("A", "A", "B", "B", "C", "C")
  result <- create_design(treatment)

  expect_is(result$design, "matrix")
  expect_equal(ncol(result$design), 3)
})

# ============================================================================
# Test filter_db_against_eset() - Database Filtering
# ============================================================================

test_that("filter_db_against_eset retains only genes in eset", {
  eset <- create_test_eset(n_genes = 50, n_samples = 10)
  databases <- list(
    db1 = list(
      IL6 = c("GENE_1", "GENE_2", "GENE_1000"),
      TNF = c("GENE_4", "GENE_5", "GENE_6")
    )
  )

  result <- filter_db_against_eset(eset, databases)

  expect_false("GENE_1000" %in% result$db1$IL6)
  expect_true("GENE_1" %in% result$db1$IL6)
})

test_that("filter_db_against_eset processes all database ligands", {
  eset <- create_test_eset()
  databases <- create_test_db()

  result <- filter_db_against_eset(eset, databases)

  expect_true("IL6" %in% names(result$database1))
  expect_true("TNF" %in% names(result$database1))
  expect_true("IL10" %in% names(result$database1))
})

test_that("filter_db_against_eset removes ligands with no matching receptors", {
  eset <- create_test_eset(n_genes = 50, n_samples = 10)
  databases <- list(
    db1 = list(
      IL6 = c("GENE_1000", "GENE_1001"),
      TNF = c("GENE_4", "GENE_5")
    )
  )

  result <- filter_db_against_eset(eset, databases)

  expect_false("IL6" %in% names(result$db1))
  expect_true("TNF" %in% names(result$db1))
})

# ============================================================================
# Test remove_zero_variance_ligands() - Zero Variance Filtering
# ============================================================================

test_that("remove_zero_variance_ligands processes database", {
  eset <- create_test_eset()
  databases <- create_test_db()

  result <- remove_zero_variance_ligands(eset, databases)

  expect_is(result, "list")
  expect_true("database1" %in% names(result))
})

test_that("remove_zero_variance_ligands maintains ligand structure", {
  eset <- create_test_eset()
  databases <- create_test_db()

  result <- remove_zero_variance_ligands(eset, databases)

  expect_is(result$database1, "list")
})

test_that("remove_zero_variance_ligands handles single-gene ligands", {
  eset <- create_test_eset()
  databases <- list(
    db1 = list(
      SingleGene = c("GENE_1"),
      MultiGene = c("GENE_2", "GENE_3")
    )
  )

  result <- remove_zero_variance_ligands(eset, databases)

  expect_is(result$db1, "list")
})

# ============================================================================
# Test preprocess_eset() - Full Preprocessing Pipeline
# ============================================================================

test_that("preprocess_eset returns filtered eset and databases", {
  eset <- create_test_eset()
  databases <- create_test_db()

  result <- preprocess_eset(eset, databases)

  expect_is(result, "list")
  expect_true("eset_f" %in% names(result))
  expect_true("dbs_f" %in% names(result))
})

test_that("preprocess_eset returns matrix for filtered eset", {
  eset <- create_test_eset()
  databases <- create_test_db()

  result <- preprocess_eset(eset, databases)

  expect_is(result$eset_f, "matrix")
})

test_that("preprocess_eset maintains sample integrity", {
  eset <- create_test_eset(n_genes = 100, n_samples = 10)
  databases <- create_test_db()

  result <- preprocess_eset(eset, databases)

  expect_equal(ncol(result$eset_f), ncol(eset))
})

test_that("preprocess_eset filters databases to genes in eset", {
  eset <- create_test_eset(n_genes = 50, n_samples = 10)
  databases <- list(
    db1 = list(
      IL6 = c("GENE_1", "GENE_2", "GENE_1000")
    )
  )

  result <- preprocess_eset(eset, databases)

  all_genes <- unique(unlist(result$dbs_f))
  expect_true(all(all_genes %in% rownames(eset)))
})

# ============================================================================
# Test clean_eset() - Gene Name Standardization
# ============================================================================

test_that("clean_eset maps probe IDs to gene symbols", {
  eset <- create_test_eset(n_genes = 10)
  gene_list <- data.frame(
    probeids = rownames(eset)[1:5],
    gensym = c("TP53", "TP53", "EGFR", "EGFR", "MYC")
  )

  result <- clean_eset(eset, gene_list)

  expect_true("TP53" %in% rownames(result))
  expect_true("EGFR" %in% rownames(result))
  expect_true("MYC" %in% rownames(result))
})

test_that("clean_eset collapses multiple probes to same gene", {
  eset <- create_test_eset(n_genes = 10)
  gene_list <- data.frame(
    probeids = rownames(eset)[1:2],
    gensym = c("TP53", "TP53")
  )

  result <- clean_eset(eset, gene_list)

  tp53_count <- sum(rownames(result) == "TP53")
  expect_equal(tp53_count, 1)
})

test_that("clean_eset returns matrix", {
  eset <- create_test_eset(n_genes = 10)
  gene_list <- data.frame(
    probeids = rownames(eset)[1:5],
    gensym = c("TP53", "EGFR", "MYC", "BRCA1", "PTEN")
  )

  result <- clean_eset(eset, gene_list)

  expect_is(result, "matrix")
})

# ============================================================================
# Test create_db_space() - Database Storage
# ============================================================================

test_that("create_db_space does not save when saveToFile=FALSE", {
  temp_db <- list(IL6 = c("GENE_1", "GENE_2"))
  temp_path <- tempfile(fileext = ".rda")
  temp_dir <- dirname(temp_path)

  result <- create_db_space(
    dbs_test = temp_db,
    filePath = temp_path,
    saveToFile = FALSE
  )

  expect_false(file.exists(temp_path))
})

test_that("create_db_space saves file when saveToFile=TRUE", {
  temp_db <- list(IL6 = c("GENE_1", "GENE_2"))
  temp_path <- tempfile(fileext = ".rda")

  create_db_space(
    dbs_test = temp_db,
    filePath = temp_path,
    saveToFile = TRUE
  )

  expect_true(file.exists(temp_path))

  if (file.exists(temp_path)) {
    file.remove(temp_path)
  }
})

test_that("create_db_space errors with non-existent directory", {
  temp_db <- list(IL6 = c("GENE_1", "GENE_2"))
  bad_path <- "/nonexistent/directory/db.rda"

  expect_error(
    create_db_space(
      dbs_test = temp_db,
      filePath = bad_path,
      saveToFile = TRUE
    ),
    "does not exist"
  )
})
