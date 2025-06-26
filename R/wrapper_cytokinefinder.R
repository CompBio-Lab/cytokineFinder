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
#' 
run_lri_methods <- function(eset, design, dbs, methods,
                           treatment = NULL, obs_id = NULL, correlation = NULL) {
  # Set up the future plan
  future::plan(future::multicore)  # Set up multicore parallelism
  
  # Initialize an empty list to store results
  results <- list()
  
  # Iterate over each method in the methods list
  for (method_name in methods) {
    method <- get(method_name)
    
    # Parallel execution for the current method
    method_results <- future_lapply(names(dbs), function(database) {
      message(paste("Processing method:",
                    method_name,
                    "with database:",
                    database)
      )  # Debug statement
      
      # Check if the method is PLSDA-based:
      if (grepl("plsda", method_name)) {
        # Run the plsda method
        result <- method(eset, treatment, dbs[[database]])
      } else {
        # Run non-PLSDA method
        ## Check if it is a paired experiment by checking obs_id
        if (!is.null(obs_id)) {
          result <- method(eset, design,
                           dbs[[database]],
                           obs_id = obs_id,
                           correlation = correlation)
        } else{
          result <- method(eset, design, dbs[[database]])
        }
      }
      
      message(paste("Finished processing method:",
                    method_name,
                    "with database:",
                    database)
      )  # Debug statement
      return(list(database = database, result = result))
    }, future.seed = TRUE)  # Set future.seed to ensure reproducibility
    
    # Organize results into a named list
    message(paste("Combining into one BenchmarkResults Object"))
    method_results_named <- setNames(
      lapply(method_results, `[[`, "result"), # extract method results
      sapply(method_results, `[[`, "database") # extract database name
    )
    # Store the results in the main results list
    results[[method_name]] <- method_results_named
  }
  
  # Create and return results object
  return(structure(results, class = "BenchmarkResults"))
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
#' methods <- c("gsea", "gva_limma", "cytosig_custom_ridge")
#' 
#' # Run workflow
#' results <- run_cytokinefinder_workflow(study_data, dbs_all, methods)
#' 
#' # Access results
#' benchmark_results <- results$benchmarks
#' design_matrix <- results$design
#' }
#' Wrapper function to run complete cytokinefinder benchmarking workflow
#'
#' @param study_data Named list containing study data with required elements
#' @param databases List of LRI databases for benchmarking  
#' @param methods Vector of method names to benchmark
#'
#' @return Original study_data with added results
#' @export
run_cytokinefinder_workflow <- function(study_data, databases, methods) {
  # Validate required elements
  required_elements <- c("qc_eset", "cond")
  missing <- setdiff(required_elements, names(study_data))
  if(length(missing) > 0) {
    stop("Missing required elements: ", paste(missing, collapse = ", "))
  }
  
  # Separate methods by type
  modelbased_methods <- c("cytosig_custom_ridge")
  lri_methods <- setdiff(methods, modelbased_methods)
  
  # Initialize results
  all_results <- list()
  
  # Process LRI methods (need preprocessing)
  if (length(lri_methods) > 0) {
    preprocess_results <- preprocess_eset(eset = study_data$qc_eset, dbs = databases)
    
    obs_id <- study_data$obs_id %||% NULL
    
    if (!is.null(obs_id)) {
      design_results <- create_design(study_data$cond, 
                                      obs_id = obs_id, 
                                      eset = preprocess_results$eset_f)
      
      lri_results <- cytokinefinder(eset = preprocess_results$eset_f, 
                                    design = design_results$design, 
                                    dbs = preprocess_results$dbs_f, 
                                    methods = lri_methods,
                                    treatment = study_data$cond, 
                                    obs_id = obs_id, 
                                    correlation = design_results$dupcor$consensus)
    } else {
      design_results <- create_design(study_data$cond, 
                                      eset = preprocess_results$eset_f)
      
      lri_results <- cytokinefinder(eset = preprocess_results$eset_f,
                                    design = design_results$design,
                                    dbs = preprocess_results$dbs_f, 
                                    methods = lri_methods, 
                                    treatment = study_data$cond)
    }
    
    all_results <- c(all_results, lri_results)
    study_data$preprocessing <- preprocess_results
    study_data$design <- design_results
  }
  
  # Process model-based methods (use original qc_eset)
  if (length(modelbased_methods) > 0) {
    obs_id <- study_data$obs_id %||% NULL
    
    if (!is.null(obs_id)) {
      design_results_original <- create_design(study_data$cond, 
                                               obs_id = obs_id, 
                                               eset = study_data$qc_eset)
      
      model_results <- cytokinefinder(eset = study_data$qc_eset, 
                                      design = design_results_original$design, 
                                      dbs = NULL,  # Not needed for model methods
                                      methods = modelbased_methods,
                                      treatment = study_data$cond, 
                                      obs_id = obs_id, 
                                      correlation = design_results_original$dupcor$consensus)
    } else {
      design_results_original <- create_design(study_data$cond, 
                                               eset = study_data$qc_eset)
      
      model_results <- cytokinefinder(eset = study_data$qc_eset,
                                      design = design_results_original$design,
                                      dbs = NULL,  # Not needed
                                      methods = modelbased_methods, 
                                      treatment = study_data$cond)
    }
    
    all_results <- c(all_results, model_results)
    
    # Store original design if no LRI methods were run
    if (length(lri_methods) == 0) {
      study_data$design <- design_results_original
    }
  }
  
  study_data$benchmarks <- structure(all_results, class = "BenchmarkResults")
  return(study_data)
}