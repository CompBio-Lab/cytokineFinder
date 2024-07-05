#' Helper fun to create the design matrix
#'
#' @param treatment Treatment response variable
#' @param obs_id Observation ID: some samples may have unique IDs but come from
#' the same tissue of origin, if that exists, provide a vector of this to make
#' sure the Expression matrix accounts for this to avoid incorrect DEA input.
#' @return
#' @export
#'
#' @examples

create_design <- function(treatment, obs_id=NULL){
  if(is.null(obs_id)) {
    #for unpaired datasets
    design <- model.matrix(~treatment)
  } else{
    design <- model.matrix(~treatment+obs_id)
  }
  return(design)
}
