# Example usage:
# Assuming 'res' is the output from extract_ligand
# plot_ligand_summary(res)

#' Helper function to extract specific ligands and merge results
#'
#' This function extracts specified ligands from a `BenchmarkResults` object 
#' and merges the results into a single data frame.
#'
#' @param benchmark_results a BenchmarkResults object that contains a nested list of results
#' @param ligands a vector of ligands to extract 
#'
#' @return a data frame containing the extracted results for the specified ligands.
#' @export
#' @name extract_ligands
#' @import dplyr
#' @import purrr
#' 
#' @examples
#' # Example usage requires running cytokinefinder() first, 
#' # for details check ?cytokinefinder():
#' \dontrun{
#' results_df <- extract_ligands(
#'     benchmark_results = results, 
#'     ligands = c("LigandA", "LigandB"))
#' }

# Function to extract specific ligands and merge results
extract_ligands <- function(benchmark_results, ligands) {
  # Convert ligands to a vector if it's not already
  ligands <- as.vector(ligands)
  
  # Use map_dfr to iterate over methods and databases
  results_df <- map_dfr(
    names(benchmark_results), # Iterate over method names
    function(method_name) {
      map_dfr(
        names(benchmark_results[[method_name]]), # Iterate over database names
        function(db_name) {
          process_method_db(
            benchmark_results[[method_name]][[db_name]],
            method_name,
            db_name
          )
        }
      )
    }
  )
  
  return(results_df)
}

#' Define a function to process each method and database
#'
#' @param df 
#' @param method_name 
#' @param db_name 
#'
#' @return
#' @export
#'
#' @examples

process_method_db <- function(df, method_name, db_name) {
  # Extract and filter data for the specified ligands
  reindex_rank_order <- df %>%
    # add column and populate value with method and database
    mutate(database = db_name, method = method_name) %>%
    # reorder by lowest p-value 
    arrange(pval) %>%
    mutate(rank = 100*(1 - (row_number()/n()))) %>%
    filter(ligand %in% ligands)
  
  return(reindex_rank_order)
}