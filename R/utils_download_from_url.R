#' Download a file from a URL if it doesn't exist locally.
#' This function checks if the specified destination file exists locally.
#' If not, it downloads the file from the provided URL to the destination file.
#' Optionally, it can extract the contents if the file is a zip archive.
#'
#' @param url The URL of the file to download.
#' @param destFile Destination file path where the downloaded file should be saved.
#' @param mode Optional. Mode for downloading the file. Default is NULL (automatic detection).
#'             Use "wb" for binary files like zip archives.
#' @param extractZip Logical, whether to extract the file if it's a zip archive (default is FALSE).
#' @param extractDir Optional. Directory where the contents of the zip file should be extracted.
#'                   If NULL, extracts to the current working directory.
#'
#'
#' @return NULL
#' @export
#'
#' @examples
download_from_url <- function(url, 
                              destFile, 
                              mode = NULL, 
                              extractZip = FALSE, 
                              extractDir = NULL) {
  # Check if file exists
  if (file.exists(destFile)) { 
    message(paste0("File", destFile, "already exists. Skipping download."))
  } else {
    if (is.null(mode)) {
      # Determine mode if not specified and evaluate if it is a zip file
      mode <- ifelse(grepl("\\.zip$", url, ignore.case = TRUE), "wb", "auto")
    }
    # Download File
    utils::download.file(url = url, destfile = destFile, mode = mode)
    
    # Extract if requested
    if (extractZip && grepl("\\.zip$", ignore.case = TRUE)) {
      if (is.null(extractDir)) {
        utils::unzip(destFile, exdir = dirname(destFile))
      } else {
        utils::unzip(destFile, exdir = extractDir)
      }
    }
  }
  
}

