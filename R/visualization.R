#' Plot Ligand of Interest summary result across methods and databases
#'
#' @param data 
#'
#' @return
#' @export
#' @import dplyr
#' @import ggplot2
#' @examples

plot_ligand_summary <- function(data) {
  # Check if the required columns are present in the dataframe
  required_columns <- c("method", "pathway", "pval", "database")
  missing_columns <- setdiff(required_columns, colnames(data))
  
  if (length(missing_columns) > 0) {
    stop("The following required columns are missing: ", paste(missing_columns, collapse = ", "))
  }
  
  # Summarize the data
  summary_data <- data %>%
    group_by(method, pathway)
  
  # Create the plot
  ggplot(summary_data, aes(x = reorder(method, pval), y = pval, fill = database)) +
    geom_bar(stat = "identity", position = "dodge") +
    #geom_errorbar(aes(ymin = mean_pval - sd_pval, ymax = mean_pval + sd_pval), width = .2, position = position_dodge(.9)) +
    ylab("p-value") +
    xlab("Method + Database Annotation") +
    facet_wrap(~pathway, ncol = 1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

# Example usage:
# Assuming 'res' is the output from extract_ligand
# plot_ligand_summary(res)

#' Helper function to extract specific ligands and merge results
#'
#' @param benchmark_results 
#' @param ligands 
#'
#' @return
#' @export
#' 
#' @import dplyr
#' @import map
#' 
#' @examples

library(purrr)
library(dplyr)

# Function to extract specific ligands and merge results
extract_ligands <- function(benchmark_results, ligands) {
  # Convert ligands to a vector if it's not already
  ligands <- as.vector(ligands)
  
  # Define a function to process each method and database
  process_method_db <- function(df, method_name, db_name) {
    # Extract and filter data for the specified ligands
    filtered_results <- df %>%
      filter(pathway %in% ligands) %>%
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

# Example usage
# results_df <- extract_ligands(benchmark_results, c("LigandA", "LigandB"))
