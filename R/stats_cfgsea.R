#' Rank ligands using fGSEA on limma t-statistics
#'
#' Runs limma on the full eset to obtain per-gene t-statistics, then uses
#' \code{fgsea} to score each ligand's receptor gene set.
#'
#' @param eset Numeric matrix (genes x samples).
#' @param design Design matrix from \code{create_design()}.
#' @param db Named list: ligand -> character vector of receptor gene names.
#' @param obs_id Optional donor ID vector for paired analysis.
#' @param correlation Optional within-block correlation.
#'
#' @return data.frame with fgsea results; one row per ligand. Rows where fgsea
#'   returns \code{pval = NA} (gene set too small to permute) are retained with
#'   \code{pval = NA}. The ranking step sets \code{pct_rank = NA} for these
#'   rows rather than silently ranking them last.
#' @export
#'
#' @importFrom fgsea fgsea
cfgsea <- function(eset, design, db, obs_id = NULL, correlation = NULL) {

  top   <- run_limma(eset, design, obs_id, correlation)
  stats <- setNames(top$t, top$genes)

  length_receptors <- sapply(db, length)

  fgsea_results <- fgsea::fgsea(
    pathways = db,
    stats    = stats,
    minSize  = max(1L, min(length_receptors)),
    maxSize  = max(length_receptors)
  )
  fgsea_results <- dplyr::rename(fgsea_results, ligand = pathway)
  return(fgsea_results)
}
