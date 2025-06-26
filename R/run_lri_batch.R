#' Wrapper function to run complete cytokinefinder benchmarking workflow with preprocessing
#'
#' @param study_data Named list containing study data with required elements:
#'   \itemize{
#'     \item qc_eset: Quality-controlled expression set matrix
#'     \item cond: Treatment/condition vector  
#'     \item obs_id: Sample IDs for paired experiments (optional)
#'   }
#' @param databases List of LRI databases for benchmarking
#' @param methods Vector of method names to benchmark
#'
#' @return Original study_data with added 'benchmarks', 'design', and 'preprocessing' elements
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
#' results <- run_cytokinefinder_workflow(study_data, dbs_all, methods)
#' 
#' # Access results
#' benchmark_results <- results$benchmarks
#' design_matrix <- results$design
#' }
run_lri_methods_batch <- function(study_data, databases, methods) {
  # Validate required elements
  required_elements <- c("qc_eset", "cond")
  missing <- setdiff(required_elements, names(study_data))
  if(length(missing) > 0) {
    stop("Missing required elements: ", paste(missing, collapse = ", "))
  }
  
  # Extract data
  eset <- study_data$qc_eset
  obs_id <- study_data$obs_id %||% NULL
  
  # preprocess
  preprocess <- preprocess_eset(eset, databases)
  eset_f <- preprocess$eset_f
  dbs_f <- preprocess$dbs_f
  
  # Create design matrix for original data (model-based methods need this)
  if (!is.null(obs_id)) {
    design_results <- create_design(study_data$cond, obs_id = obs_id, eset = eset)
    lri_method_results <- run_lri_methods(eset = eset_f,
                                          design = design_results$design,
                                          dbs = dbs_f,
                                          methods = methods,
                                          treatment = study_data$cond,
                                          obs_id = obs_id,
                                          correlation = design_results$dupcor$consensus)
  } else {
    design_results <- create_design(study_data$cond, eset = eset)
    benchmark_results <- cytokinefinder(eset = eset_f,
                                        design = design_results$design,
                                        dbs = dbs_f, 
                                        methods = methods, 
                                        treatment = study_data$cond)
  }
  
  # Add results to study data
  study_data$design <- design_results
  study_data$benchmarks <- benchmark_results
  
  return(study_data)
}