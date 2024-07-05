#' Title
#'
#' @param eset Expression Set object containing gene expression data.
#' @param y Response variable
#' @param db ligand-receptor database
#'
#' @return 
#' @export
#'
#' @examples

cpcr <- function(eset, treatment, db){
  pcs <- sapply(db, function(ligand){
    genexp <- t(eset[intersect(rownames(eset), ligand), , drop=FALSE])
    stats::prcomp(genexp, center = TRUE, scale. = TRUE, rank. = 1)$x[, "PC1"]
  })
  fit <- mixOmics::plsda(pcs, treatment)
  coef <- abs(mixOmics::selectVar(fit, comp=1)$value$value.var)
  names(coef) <- rownames(mixOmics::selectVar(fit, comp=1)$value)
  return(coef[order(coef, decreasing = TRUE)])
}