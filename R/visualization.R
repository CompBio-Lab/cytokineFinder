#' Plot Ligand of Interest summary result across methods and databases
#'
#' @param data 
#'
#' @return
#' @export
#' @name plot_ligand_summary
#' @import dplyr
#' @import ggplot2
#' @examples

plot_ligand_summary <- function(data) {
  # Check if the required columns are present in the dataframe
  required_columns <- c("method", "ligand", "pval", "database")
  missing_columns <- setdiff(required_columns, colnames(data))
  
  if (length(missing_columns) > 0) {
    stop("The following required columns are missing: ", paste(missing_columns, collapse = ", "))
  }
  
  # Summarize the data
  summary_data <- data %>%
    group_by(method, ligand) %>%
    mutate(sig = -log10(pval))
  
  # Create the plot
  ggplot(summary_data, aes(x = reorder(method, sig), y = sig, fill = database)) +
    geom_bar(stat = "identity", position = "dodge") +
    #geom_errorbar(aes(ymin = mean_pval - sd_pval, ymax = mean_pval + sd_pval), width = .2, position = position_dodge(.9)) +
    ylab("Significance") +
    xlab("Method + Database Annotation") +
    facet_wrap(~ligand, ncol = 1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

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
#' Example usage:
#' # results_df <- extract_ligands(benchmark_results, c("LigandA", "LigandB"))


# Function to extract specific ligands and merge results
extract_ligands <- function(benchmark_results, ligands) {
  # Convert ligands to a vector if it's not already
  ligands <- as.vector(ligands)
  
  # Define a function to process each method and database
  process_method_db <- function(df, method_name, db_name) {
    # Extract and filter data for the specified ligands
    filtered_results <- df %>%
      filter(ligand %in% ligands) %>%
      mutate(database = db_name, method = method_name)
    
    return(filtered_results)
  }
  
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

