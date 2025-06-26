#' The function to run all methods for benchmarking
#'
#' @param eset An expression set (numeric matrix) of genes x samples
#' @param design The design matrix for the data set used to generate the model
#' @param dbs The databases in a list of list format (for package default, use dbs_all)
#' @param methods A vector of methods contained in this package
#' @param treatment A vector containing the treatment (specific to the demo 
#' data set, this is to analyze differentially expressed genes between week 0 
#' and week 6) with the drug golimumab on Ulcerative colitis patients)
#' @param dupCor A logical TRUE/FALSE to determine whether or not duplicate correlate exists. Default is set to FALSE.
#' @param obs_id A vector of sample IDs if biological replicates are tied to a sample of origin
#' @param correlation the average estimated inter-duplicate correlation for paired sample relationships
#'
#' @return A large BenchmarkResults object containing a nested list of methods
#' and the results 
#' @export
#'
#' @importFrom future plan multicore
#' @importFrom future.apply future_lapply future_sapply
#' @examples
cytokinefinder <- function(eset, design, dbs, methods, 
                           treatment = NULL, obs_id = NULL, correlation = NULL) { 
  future::plan(future::multicore)
  results <- list()
  
  # Define which methods don't require databases
  modelbased_methods <- c("cytosig_custom_ridge")  # Add more as needed
  
  for (method_name in methods) {
    method <- get(method_name)
    
    if (method_name %in% modelbased_methods) {
      # Handle database-independent methods
      results[[method_name]] <- run_model_method(
        method, method_name, eset, design, treatment, obs_id, correlation
      )
    } else {
      # Handle database-dependent methods (existing logic)
      results[[method_name]] <- run_lri_method(
        method, method_name, eset, design, dbs, treatment, obs_id, correlation
      )
    }
  }
  return(structure(results, class = "BenchmarkResults"))
}

#' Run database-independent methods
#' @keywords internal
run_model_method <- function(method, method_name, eset, design,
                             treatment, obs_id, correlation) {
  message(paste("Processing database-independent method:", method_name))
  
  if (grepl("plsda", method_name)) {
    result <- method(eset, treatment)
  } else {
    if (!is.null(obs_id)) {
      result <- method(eset, design, obs_id = obs_id, correlation = correlation)
    } else {
      result <- method(eset, design)
    }
  }
  
  message(paste("Finished processing method:", method_name))
  
  # Return in consistent format - single result labeled as "standalone"
  return(list(standalone = result))
}

#' Run database-dependent methods
#' @keywords internal
run_lri_method <- function(method, method_name, eset, design, dbs, 
                                    treatment, obs_id, correlation) {
  method_results <- future_lapply(names(dbs), function(database) {
    message(paste("Processing method:", method_name, "with database:", database))
    
    if (grepl("plsda", method_name)) {
      result <- method(eset, treatment, dbs[[database]])
    } else {
      if (!is.null(obs_id)) {
        result <- method(eset, design, dbs[[database]], 
                         obs_id = obs_id, correlation = correlation)
      } else {
        result <- method(eset, design, dbs[[database]])
      }
    }
    
    message(paste("Finished processing method:", method_name, "with database:", database))
    return(list(database = database, result = result))
  }, future.seed = TRUE)
  
  # Organize results into named list
  message("Combining into one BenchmarkResults Object")
  method_results_named <- setNames(
    lapply(method_results, `[[`, "result"),
    sapply(method_results, `[[`, "database")
  )
  
  return(method_results_named)
}

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
run_cytokinefinder_workflow <- function(study_data, databases, methods) {
  # Validate required elements
  required_elements <- c("qc_eset", "cond")
  missing <- setdiff(required_elements, names(study_data))
  if(length(missing) > 0) {
    stop("Missing required elements: ", paste(missing, collapse = ", "))
  }
  
  # Extract and preprocess data
  eset <- study_data$qc_eset
  preprocess_results <- preprocess_eset(eset = eset, dbs = databases)
  
  # Check if obs_id exists to determine if paired
  obs_id <- study_data$obs_id %||% NULL
  
  if (!is.null(obs_id)) {
    # Paired experiment workflow
    design_results <- create_design(study_data$cond, 
                                    obs_id = obs_id, 
                                    eset = preprocess_results$eset_f)
    
    # Run benchmark using paired args
    benchmark_results <- cytokinefinder(eset = preprocess_results$eset_f, 
                                        design = design_results$design, 
                                        dbs = preprocess_results$dbs_f, 
                                        methods = methods,
                                        treatment = study_data$cond, 
                                        obs_id = obs_id, 
                                        correlation = design_results$dupcor$consensus)
  } else {
    # Unpaired experiment workflow
    design_results <- create_design(study_data$cond, 
                                    eset = preprocess_results$eset_f)
    
    # Run benchmark function without correlation or obs_id
    benchmark_results <- cytokinefinder(eset = preprocess_results$eset_f,
                                        design = design_results$design,
                                        dbs = preprocess_results$dbs_f, 
                                        methods = methods, 
                                        treatment = study_data$cond)
  }
  
  # Add results to study data
  study_data$design <- design_results
  study_data$benchmarks <- benchmark_results
  study_data$preprocessing <- preprocess_results
  
  return(study_data)
}