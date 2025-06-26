#' Run cytokine workflow on multiple studies
#'
#' @param multi_study_data Named list where each element is a study containing qc_eset, cond, obs_id
#' @param databases List of LRI databases for benchmarking  
#' @param methods Vector of method names to benchmark
#'
#' @return multi_study_data with added results for each study
#' @export
#'
#' @importFrom rlang %||%
#'
#' @examples
#' # Single study example
#' study_data <- list(
#'   qc_eset = matrix(rnorm(1000), nrow = 100, ncol = 10),
#'   cond = rep(c("control", "treatment"), each = 5)
#' )
#' 
#' # With paired samples
#' study_data_paired <- list(
#'   qc_eset = matrix(rnorm(1000), nrow = 100, ncol = 10),
#'   cond = rep(c("control", "treatment"), each = 5),
#'   obs_id = rep(1:5, 2)
#' )
#' 
#' \dontrun{
#' # Load databases and methods
#' data(dbs_all)  # or however you load your databases
#' methods <- c("gsea", "ssgsea", "cytosig_custom_ridge")
#' 
#' # Run workflow
#' results <- run_cytokine_workflow(study_data, dbs_all, methods)
#' 
#' # Access results
#' benchmark_results <- results$benchmarks
#' design_matrix <- results$design
#' }

run_cytokinefinder_workflow_batch <- function(multi_study_data, databases, methods) {
  for (study_name in names(multi_study_data)) {
    message(paste("Processing study:", study_name))
    multi_study_data[[study_name]] <- run_cytokinefinder_workflow(
      multi_study_data[[study_name]], 
      databases, 
      methods
    )
  }
  return(multi_study_data)
}