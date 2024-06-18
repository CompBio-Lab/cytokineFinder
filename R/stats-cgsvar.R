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

cgsvar <- function(eset, y, obs_id, db){
  gsvapar <- GSVA::gsvaParam(eset, db, maxDiff = TRUE)
  gsva_eset <- GSVA::gsva(gsvapar)
  fit <- mixOmics::plsda(gsva_eset, y)
  coef <- abs(mixOmics::selectVar(fit, comp=1)$value$value.var)
  names(coef) <- rownames(mixOmics::selectVar(fit, comp=1)$value)
  coef[order(coef, decreasing = TRUE)]
}