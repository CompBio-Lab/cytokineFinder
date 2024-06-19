#' Gene Set Variation Analysis for cytokines (cGSVA)
#'
#' @param eset 
#' @param design 
#' @param db 
#'
#' @return
#' @export
#'
#' @examples

cgsva <- function(eset, design, db){
  gsvapar <- GSVA::gsvaParam(eset, db, maxDiff = TRUE)
  gsva_eset <- GSVA::gsva(gsvapar)
  fit <- limma::eBayes(limma::lmFit(gsva_eset, design))
  top <- limma::topTable(fit, coef = 2, n = nrow(fit))
  pval <- top$P.Value
  names(pval) <- rownames(top)
  return(pval[order(pval)])
}