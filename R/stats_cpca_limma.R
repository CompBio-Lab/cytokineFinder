#' Calculate top ligands using a PCA approach from receptor genes given a database 
#'
#' @param eset Expression Set object containing gene expression data.
#' @param design Design matrix generated from create_design()
#' @param db Ligand-receptor database
#' @param obs_id Optional: provide a vector of sample IDs making sure the order matches with the eset
#' @param correlation Optional: input the correlation consensus between the samples to evaluate if it is paired data
#' 
#' @return a data frame of differentially expressed ligands ordered by p-values
#' @export
#'
#' @examples
#' 
#' @importFrom stats prcomp
#' @importFrom limma eBayes
#' @importFrom limma lmFit
#' @importFrom limma topTable
#' @importFrom tibble enframe
#' @importFrom tibble rownames_to_column

pca_limma <- function(eset, 
                      design, 
                      db, 
                      obs_id = NULL, 
                      correlation = NULL) {
  # Check if the design matrix is a data frame or matrix
  if (!is.data.frame(design) && !is.matrix(design)) {
    stop("The design argument must be a data frame or matrix.")
  }
  
  # Run PCA to get the first PC
  pc <- lapply(db, function(ligand){
    tryCatch({
      genexp <- t(eset[intersect(rownames(eset), ligand), , drop=FALSE])
      prcomp(genexp, center = TRUE, scale. = TRUE, rank. = 1)$x[, "PC1"]  
    }, error = function(e) NA)
    })
  # remove 0 variability ligands
  pc <- pc[sapply(pc, length) > 1] %>% 
    do.call(rbind, .)
  
  # run DEA using the design matrix integrated from previous create_design()
  # Check if paired experiment
  if (!is.null(obs_id)) {
    fit <- lmFit(pc, design, block = obs_id, correlation = correlation)
  } else {
    fit <- lmFit(pc, design)
  }
  efit <- eBayes(fit)
  top <- topTable(efit, coef = 2, number = nrow(efit)) %>%
    rownames_to_column("ligand") %>%
    dplyr::rename(pval = P.Value, padj = adj.P.Val)
  
  return(top)
}