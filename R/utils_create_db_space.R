
#' Title
#'
#' @param ... 
#' @param filePath 
#' @param saveToFile 
#'
#' @return
#' @export
#'
#' @examples

create_db_space <- function(..., 
                            filePath = "data/ligand_receptor_db.rda", 
                            saveToFile = TRUE) {  
  # Check if the data directory exists might need to come back to this
  data_dir <- dirname(filePath)
  if (!dir.exists(data_dir)) {
    stop(paste("The directory", 
               data_dir, 
               "does not exist. Please create it first."))
  }
  
  # Save the database list to the file if the save_to_file parameter is TRUE
  if (saveToFile) {
    # Save the list to the specified file path
    save(..., file = filePath)
    message("Database list saved to ", filePath)
  } else {
    message("Database list not saved.")
  }
}