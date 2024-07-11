#' Gene Set Variation Analysis for cytokines (cGSVA)
#'
#' @param eset Expression Set object containing gene expression data.
#' @param design Design matrix generated from create_design()
#' @param db Ligand-receptor database
#'
#' @return
#' @export
#'
#' @importFrom limma lmFit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' @importFrom GSVA gsvaParam
#' @importFrom GSVA gsva
#' 
#' @examples

cgsva <- function(eset, design, db) {
  gsvapar <- gsvaParam(eset, db, maxDiff = TRUE)
  gsva_eset <- gsva(gsvapar)
  
  # Run DEA
  fit <- eBayes(lmFit(gsva_eset, design))
  top <- topTable(fit, coef = 2, n = nrow(fit))
  pval <- top$P.Value
  names(pval) <- rownames(top)
  return(pval[order(pval)])
}