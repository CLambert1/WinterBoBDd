#' Adding observation on effort.
#'
#' Adding the number of individuals and the group size on prepared effort data
#' (e.g. on legdata and segdata).
#'
#' @param countdata_leg data.frame containing the detection number and the total number of
#' individuals at leg scale
#' @param legdata \code{data.frame} with infos at leg scale. Built with \code{prepare_effort_dsm}
#' @param countdata_seg data.frame containing the detection number and the total number of
#' individuals at leg scale
#' @param segdata \code{data.frame} with infos at segment scale. Built with \code{prepare_effort_dsm}
#'
#' @return This function return a list containing :
#'         \itemize{
#'           \item{legdata_obs :}{\code{data.frame} corresponding to "legdata" with observation information}
#'           \item{segdata_obs :}{\code{data.frame} corresponding to "segdata" with observation information}
#'         }
#' @export
#' @examples
#' library(pelaDSM)
#'
#' legdata <- SPEE_data_dsm$effort_output$legdata
#' segdata <- SPEE_data_dsm$effort_output$segdata
#'
#' # loop over MOLMOL and JELLY to see n_detection per session for each taxa
#'
#' taxa <- c("JELLY","MOLMOL")
#'
#' for(s in seq_len(length(taxa))) {
#'
#'   # add JELLY "n_ind" and "n_detection" on "segdata" and "legdata"
#'   countdata_leg <- SPEE_data_dsm$list_prepare_obs_by_sp[[s]]$countdata_leg
#'   countdata_seg <- SPEE_data_dsm$list_prepare_obs_by_sp[[s]]$countdata_seg
#'
#'   effort_w_obs <- add_obs_dsm(
#'     countdata_leg = countdata_leg,
#'     legdata = legdata,
#'     countdata_seg,
#'     segdata = segdata
#'   )
#'
#'   detection_by_session <- aggregate(
#'     n_detection ~ session,
#'     FUN = sum,
#'     data = effort_w_obs$legdata_obs
#'   )
#'
#'   barplot(
#'     n_detection ~ session,
#'     data = detection_by_session,
#'     main = paste0("Number of ",taxa[s]," detection per session")
#'   )
#' }
add_obs_dsm <- function(countdata_leg, legdata, countdata_seg, segdata){

  ### merge Legdata avec obs pour n_detection et n_ind ###
  legdata_obs <- merge(legdata, countdata_leg, by = c("Transect.Label","leg_id"),
                       all.x = TRUE)
  legdata_obs$n_detection[is.na(legdata_obs$n_detection)] <- 0
  legdata_obs$n_ind[is.na(legdata_obs$n_ind)] <- 0

  ### merge Segdata avec obs pour n_detected et n_ind ###
  segdata_obs <- merge(segdata, countdata_seg, by = c("Transect.Label","leg_id","Sample.Label"),
                       all.x = TRUE)
  segdata_obs$n_detection[is.na(segdata_obs$n_detection)] <- 0
  segdata_obs$n_ind[is.na(segdata_obs$n_ind)] <- 0

  return(list(segdata_obs = segdata_obs, legdata_obs = legdata_obs))
}

