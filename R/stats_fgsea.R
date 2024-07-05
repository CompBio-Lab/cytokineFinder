#' Fast Gene Set Enrichment Analysis for cytokines (cFGSEA)
#'
#' @description
#' Calculate regulatory activity of a list of receptors using FGSEA given some
#' cytokine of interest
#'
#' @details
#' GSEA (Aravind et al., 2005) starts by transforming the input expression
#' matrix to ranks for each sample. Then, an enrichment score `fgsea` is
#' calculated by walking down the list of features, increasing
#' a running-sum statistic when a feature in the target feature set is
#' encountered and decreasing it when it is not. The final score is the maximum
#' deviation from zero encountered in the random walk. Finally, a normalized
#' score `norm_fgsea`, can be obtained by computing the z-score of the estimate
#' compared to a null distribution obtained from N random permutations. The used
#' implementation is taken from the package `fgsea` (Korotkevich et al., 2021).
#'
#' Aravind S. et al. (2005) Gene set enrichment analysis: A knowledge-based
#' approach for interpreting genome-wide expression profiles. PNAS. 102, 43.
#'
#' Korotkevich G. et al. (2021) Fast gene set enrichment analysis. bioRxiv.
#' DOI: https://doi.org/10.1101/060012.
#'
#' @param eset Expression Set object containing gene expression data.
#' @param treatment Optional Response variable
#' @param db ligand-receptor database
#' 
#' @return GSVA result for set of ligands that are significantly enriched
#' ordered by p-value
#' 
#' @importFrom limma eBayes
#' @importFrom limma lmFit
#' @importFrom limma topTable
#' @importFrom fgsea fgsea
#' 

cfgsea <- function(eset, treatment, db){
  # Initialize design matrix
  design <- create_design(treatment, obs_id)
  
  # run DEA
  fit <- eBayes(lmFit(eset, design))
  top <- topTable(fit, coef = 2, n = nrow(fit))
  
  # run fgsea
  stats = top$t
  names(stats) <- rownames(top)
  run_fgsea <- fgsea(pathways = db,
                     stats    = stats,
                     minSize  = min(sapply(db, length)),
                     maxSize  = max(sapply(db, length)))
  result <- run_fgsea$pval
  names(result) <- run_fgsea$pathway
  return(result[order(result)])
}