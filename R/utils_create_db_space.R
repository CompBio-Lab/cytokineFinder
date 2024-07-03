
#' Title
#'
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples

create_db_space <- function(..., file_path) {
  # Check if the required package is installed
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package 'here' is required but not installed.")
  }
  # Save to file path
  save(..., file = file_path)
}