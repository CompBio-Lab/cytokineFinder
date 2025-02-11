#' Plot Ligand of Interest summary result across methods and databases
#'
#' @param data 
#'
#' @return A bar plot figure showing the different methods+database combinations
#' @export
#' @name plot_ligand_summary
#' @import dplyr
#' @import ggplot2
#' @examples
#' #Example usage:
#' #Assuming 'res' is the output from extract_ligand
#' \dontrun{
#' plot_ligand_summary(res)
#' }

plot_ligand_summary <- function(data) {
  # Check if the required columns are present in the dataframe
  required_columns <- c("method", "ligand", "database", "rank")
  missing_columns <- setdiff(required_columns, colnames(data))
  
  if (length(missing_columns) > 0) {
    stop("The following required columns are missing: ", paste(missing_columns, collapse = ", "))
  }
  
  # Summarize the data
  summary_data <- data %>%
    group_by(method, ligand) #%>%
    #mutate(sig = -log10(pval))
  
  # Create the plot
  plot <- ggplot(summary_data, aes(x = reorder(method, rank), y = rank, fill = database)) +
    geom_bar(stat = "identity", position = "dodge") +
    #geom_errorbar(aes(ymin = mean_pval - sd_pval, ymax = mean_pval + sd_pval), width = .2, position = position_dodge(.9)) +
    ylab("Percentile Rank") +
    xlab("Method + Database Annotation") +
    facet_wrap(~ligand, ncol = 1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank())
  
  return(plot)
}