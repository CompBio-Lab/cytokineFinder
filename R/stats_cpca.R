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

cpca <- function(eset, y, obs_id, db){
  pc <- t(sapply(db, function(ligand){
    genexp <- t(eset[intersect(rownames(eset), ligand), , drop=FALSE])
    stats::prcomp(genexp, center = TRUE, scale. = TRUE, rank. = 1)$x[, "PC1"]
  }))
  design <- create_design(y, obs_id)
  fit <- limma::eBayes(limma::lmFit(pc, design))
  top <- limma::topTable(fit, coef = 2, n = nrow(fit))
  pval <- top$P.Value
  names(pval) <- rownames(top)
  pval[order(pval)]
}