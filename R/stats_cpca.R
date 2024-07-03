#' Calculate top ligands using a PCA approach and use for weights to 
#'
#' @param eset Expression Set object containing gene expression data.
#' @param y 
#' @param obs_id 
#' @param db 
#'
#' @return List of differentially expressed ligands ordered by p-values
#' @export
#'
#' @examples
#' 
#' @importFrom stats prcomp
#' @importFrom limma eBayes
#' @importFrom limma lmFit
#' @importFrom limma topTable

cpca <- function(eset, y, obs_id, db){
  pc <- t(sapply(db, function(ligand){
    genexp <- t(eset[intersect(rownames(eset), ligand), , drop=FALSE])
    prcomp(genexp, center = TRUE, scale. = TRUE, rank. = 1)$x[, "PC1"]
  }))
  design <- create_design(y, obs_id)
  fit <- eBayes(lmFit(pc, design))
  top <- topTable(fit, coef = 2, n = nrow(fit))
  pval <- top$P.Value
  names(pval) <- rownames(top)
  return(pval[order(pval)])
}