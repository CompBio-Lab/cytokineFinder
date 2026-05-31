#' Rank ligands using PCA on receptor gene sets followed by limma DE
#'
#' For each ligand in \code{db}, extracts its receptor genes from \code{eset},
#' computes PC1 across samples, then runs limma to test whether PC1 differs
#' between conditions. Ligands with fewer than 2 receptor genes with non-zero
#' variance are silently dropped.
#'
#' @param eset Numeric matrix (genes x samples).
#' @param design Design matrix from \code{create_design()}.
#' @param db Named list: ligand -> character vector of receptor gene names.
#' @param obs_id Optional donor ID vector for paired analysis.
#' @param correlation Optional within-block correlation.
#'
#' @return data.frame with columns: ligand, logFC, AveExpr, t, pval, padj, B.
#'   Returns \code{NULL} if no ligands survive filtering.
#' @export
#'
#' @importFrom stats prcomp
pca_limma <- function(eset, design, db, obs_id = NULL, correlation = NULL) {

  if (!is.data.frame(design) && !is.matrix(design))
    stop("design must be a data.frame or matrix.")

  # Compute PC1 for each ligand's receptor gene set
  pc_list <- lapply(db, function(receptors) {
    tryCatch({
      genexp <- t(eset[intersect(rownames(eset), receptors), , drop = FALSE])
      # Bug fix: pre-filter zero-variance receptor genes before PCA.
      # prcomp with scale.=TRUE errors on zero-variance columns, which previously
      # caused the entire ligand to be dropped even when most receptors were
      # informative. Strip zero-variance genes first; the length>1 check below
      # then correctly drops only ligands with genuinely insufficient signal.
      nonzero_cols <- apply(genexp, 2, var) > 0
      genexp <- genexp[, nonzero_cols, drop = FALSE]
      prcomp(genexp, center = TRUE, scale. = TRUE, rank. = 1)$x[, "PC1"]
    }, error = function(e) NA)
  })

  # Keep only ligands where PC1 has >1 value (i.e., >=2 receptor genes with non-zero variance)
  keep    <- sapply(pc_list, function(x) length(x) > 1)
  pc_list <- pc_list[keep]
  if (length(pc_list) == 0) return(NULL)

  # BUG FIX 1: save ligand names before rbind.
  # do.call(rbind, single-element-list) drops the name and produces rowname "1".
  lig_names <- names(pc_list)
  pc        <- do.call(rbind, pc_list)
  rownames(pc) <- lig_names   # restore correct ligand names

  # Delegate to run_limma() — consistent paired/unpaired handling (refactor R2)
  top <- run_limma(pc, design, obs_id = obs_id, correlation = correlation)

  # run_limma returns column "genes"; rename to "ligand"
  top <- dplyr::rename(top, ligand = genes, pval = P.Value, padj = adj.P.Val)
  return(top)
}
