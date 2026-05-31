#' Run limma in unpaired or paired mode
#'
#' @param eset Numeric matrix (genes x samples).
#' @param design Design matrix from \code{create_design()}.
#' @param obs_id Optional donor ID vector for paired analysis.
#' @param correlation Optional within-block correlation from
#'   \code{limma::duplicateCorrelation()}.
#'
#' @return data.frame with columns: genes, logFC, AveExpr, t, P.Value, adj.P.Val, B
#' @export
#'
#' @importFrom limma lmFit eBayes topTable
#' @importFrom tibble rownames_to_column
run_limma <- function(eset, design, obs_id = NULL, correlation = NULL) {
  if (!is.null(obs_id)) {
    fit <- limma::lmFit(eset, design, block = obs_id, correlation = correlation)
  } else {
    fit <- limma::lmFit(eset, design)
  }
  efit <- limma::eBayes(fit)
  top  <- limma::topTable(efit, coef = 2, number = nrow(efit)) %>%
    tibble::rownames_to_column(var = "genes")
  return(top)
}
