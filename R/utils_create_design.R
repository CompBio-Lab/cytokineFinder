#' Helper fun to create the design matrix
#'
#' @param treatment Treatment response variable
#' @param obs_id Observation ID or sample if looking there are biological replicates
#'
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
