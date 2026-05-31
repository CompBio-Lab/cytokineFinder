#' Rank ligands using GSVA scores followed by limma DE
#'
#' Computes GSVA enrichment scores for each ligand's receptor gene set, then
#' runs limma to test whether scores differ between conditions.
#'
#' @param eset Numeric matrix (genes x samples).
#' @param design Design matrix from \code{create_design()}.
#' @param db Named list: ligand -> character vector of receptor gene names.
#' @param obs_id Optional donor ID vector for paired analysis.
#' @param correlation Optional within-block correlation.
#'
#' @return data.frame with columns: ligand, logFC, AveExpr, t, pval, padj, B.
#' @export
#'
#' @importFrom GSVA gsva gsvaParam
gsva_limma <- function(eset, design, db, obs_id = NULL, correlation = NULL) {

  # Compute GSVA enrichment scores (ligands x samples)
  params      <- GSVA::gsvaParam(eset, db)
  gsva_scores <- GSVA::gsva(params, verbose = FALSE)

  # Delegate to run_limma() — consistent paired/unpaired handling (refactor R2)
  top <- run_limma(gsva_scores, design, obs_id = obs_id, correlation = correlation)

  top <- dplyr::rename(top, ligand = genes, pval = P.Value, padj = adj.P.Val)
  return(top)
}
