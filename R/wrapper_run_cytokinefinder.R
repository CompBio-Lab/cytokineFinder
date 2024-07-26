#' The wrapper function to run the benchmarking
#'
#' @param eset 
#' @param design 
#' @param dbs 
#' @param methods
#'
#' @return A large BenchmarkResults object containing a nested list of methods
#' and the results 
#' @export
#'
#' @importFrom future plan multicore
#' @importFrom future.apply future_lapply future_sapply
#' @examples

cytokinefinder <- function(eset, design, dbs, methods) { 
  # Set up the future plan
  future::plan(future::multicore)  # Set up multicore parallelism
  
  # Initialize an empty list to store results
  results <- list()
  
  # Iterate over each method in the methods list
  for (method_name in methods) {
    method <- get(method_name)
    
    # Parallel execution for the current method
    method_results <- future_lapply(names(dbs), function(database) {
      print(paste("Processing method:", 
                  method_name, 
                  "with database:", 
                  database)
            )  # Debug statement
      result <- method(eset, design, dbs[[database]])
      print(paste("Finished processing method:", 
                  method_name, 
                  "with database:", 
                  database)
            )  # Debug statement
      return(list(database = database, result = result))
    }, future.seed = TRUE)  # Set future.seed to ensure reproducibility
    
    # Organize results into a named list
    print(paste("Combining into one BenchmarkResults Object"))
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