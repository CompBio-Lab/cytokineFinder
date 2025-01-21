#' extract specific ligands and merge results
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
  # Look at only the ligands of interest
  summary_df <- summarize_df(results_df) %>%
    filter(ligand %in% ligands)
  
  return(summary_df)
}

#' Define a function to process each method and database
#'
#' @param df 
#' @param method_name 
#' @param db_name 
#'
#' @return
#' @export
#' @import dplyr
#' @examples

process_method_db <- function(df, method_name, db_name) {
  # Determine if the dataframe contains 'padj' or 'coef'
  if ("padj" %in% colnames(df)) {
    metric_col <- "padj"
  } else if ("coef" %in% colnames(df)) {
    metric_col <- "coef"
  } else {
    stop("Data frame must contain either 'padj' or 'coef' column.")
  }
  # Extract and filter data for the specified ligands
  reindex_rank_order <- df %>%
    # add column and populate value with method and database
    mutate(database = db_name, method = method_name) %>% 
    # Dynamically sort by metric_col
    arrange(if (metric_col == "padj") {
      padj
    } else {
      desc(coef)
    }) %>%
    mutate(rank = 100*(1 - (row_number()/n()))) 
  
  return(reindex_rank_order)
}

#' Primary helper function to create the summary DataFrame
#'
#' @param results_df 
#'
#' @return
#' @export
#' @import dplyr
#' @import tidyr
#' @examples

summarize_df <- function(results_df) {
  # Reshape both p-values and coefficients using the helper function
  pval_summary <- reshape_metric(results_df, "padj", "padj")
  #coef_summary <- reshape_metric(results_df, "coef", "coef")
  
  # Combine both summaries into one final DataFrame
  final_summary_df <- pval_summary %>%
    arrange(method, type, rank)  # Sort if needed
  
  return(final_summary_df)
}


#' Secondary helper function to reshape a specific metric into the desired format
#'
#' @param df 
#' @param metric 
#' @param type 
#'
#' @return
#' @export
#' 
#' @import dplyr
#' @import tidyr
#'
#' @examples

reshape_metric <- function(df, metric, type) {
  reshaped_df <- df %>%
    select(ligand, !!sym(metric), method, database, rank) %>%  # Use !!sym to dynamically select the metric column
    rename(value = !!sym(metric)) %>%
    mutate(type = type) %>%
    filter(!is.na(value))  # Filter out NA values
  
  return(reshaped_df)
}
