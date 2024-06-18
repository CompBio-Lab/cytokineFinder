
#' Title
#'
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples

create_db_space <- function(...) {
  # Check if the required package is installed
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package 'here' is required but not installed.")
  }
  # Create the file path
  file_path <- here::here("data", "ligand_receptor_db.RData")
  
  save(..., file = file_path)
}