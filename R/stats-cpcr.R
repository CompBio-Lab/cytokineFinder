#' Title
#'
#' @param eset 
#' @param y 
#' @param obs_id 
#' @param db 
#'
#' @return
#' @export
#'
#' @examples

cpcr <- function(eset, y, obs_id, db){
  pcs <- sapply(db, function(ligand){
    genexp <- t(eset[intersect(rownames(eset), ligand), , drop=FALSE])
    prcomp(genexp, center = TRUE, scale. = TRUE, rank. = 1)$x[, "PC1"]
  })
  fit <- mixOmics::plsda(pcs, y)
  coef <- abs(mixOmics::selectVar(fit, comp=1)$value$value.var)
  names(coef) <- rownames(mixOmics::selectVar(fit, comp=1)$value)
  coef[order(coef, decreasing = TRUE)]
}