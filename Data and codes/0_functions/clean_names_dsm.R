#' Standardization of effort data column names.
#'
#' Attribute standard names to column names of effort data.frame
#' for the next steps of the analysis.
#'
#' @param effort_base sf object containing column names to standardize
#' @param multi_survey logical, whether several surveys are included in effort_base
#' @param unit_km Are "seg_length" and "leg_length" in *km* ?
#' \itemize{
#'   \item If \code{TRUE}, there are considered in *km* and set in *km* with \code{units::set_units()}
#'   \item If \code{FALSE}, there are considered in *m* and set to *m* with \code{units::set_units()}
#' }
#' @return This function return the same sf object as in input with the standardized column name of class *cleaned_effort_dsm*
#' @note When "region_label" is not found at the end of the function, it creates a region_label
#'       column automatically corresponding to the merge of "sub_region" and "strate". There are merged by "_".
#'       Example, for a row, if \code{strate = N1} and \code{sub_region = ATL} we will get
#'       \code{region_label = ATL_N1}.
#' @import assertthat sf
#' @importFrom units set_units
#' @importFrom usethis ui_done
#' @importFrom lubridate int_end interval ymd_hms
#' @export
#' @examples
#' library(pelaDSM)
#'
#' raw_effort <- SPEE_data_dsm$effort_raw
#'
#' cleaned_effort <- clean_effort_dsm(
#'   effort_base = raw_effort,
#'   unit_km = TRUE # "seg_length" and "leg_length" are in km
#' )
clean_effort_dsm <- function(effort_base, multi_survey = FALSE, unit_km = TRUE){

  assert_that(any(is(effort_base) == "sf"))

  col_name_neces <- c("sea_state","subjective","sub_region","strate", "survey",
                      "region_label","transect","leg_id","leg_length","seg_id",
                      "seg_length","left","right", "center_time")

  if(all(col_name_neces %in% colnames(effort_base))){

    return(effort_base)

  } else {

    # Put everything tolowercase and without "_" and check if there are matches with colnames ----
    lower_colnames <- tolower(colnames(effort_base))
    lower_no_under_colnames <- gsub("_", "", lower_colnames)
    lower_neces <- tolower(col_name_neces)

    for (n in 1:length(lower_neces)) {
      if (lower_neces[n] %in% lower_no_under_colnames) {
        colnames(effort_base)[lower_no_under_colnames %in% lower_neces[n]] <- col_name_neces[n]
      }
    }

    # Are there all the needed columns ? ----
    missing_needed_col <- which(!(col_name_neces %in% colnames(effort_base)))

    if (!("sea_state" %in% colnames(effort_base))) {
      colnames(effort_base)[lower_no_under_colnames %in% c("seastate","beaufort")] <- "sea_state"
    }
    if (!("sub_region" %in% colnames(effort_base))) {
      colnames(effort_base)[lower_no_under_colnames %in% c("subregion")] <- "sub_region"
    }
    if (!("region_label" %in% colnames(effort_base))) {
      colnames(effort_base)[lower_no_under_colnames %in% c("regionlabel","stratesec")] <- "region_label"
    }
    if (!("leg_id" %in% colnames(effort_base))) {
      colnames(effort_base)[lower_no_under_colnames %in% c("legid")] <- "leg_id"
    }
    if (!("leg_length" %in% colnames(effort_base))) {
      colnames(effort_base)[lower_no_under_colnames %in% c("leglengthkm","leglengkm")] <- "leg_length"
    }
    if(!any(is(effort_base$leg_length)=="units")) {
      if(unit_km == TRUE) {
        effort_base$leg_length <- units::set_units(effort_base$leg_length, km)
        ui_done("\"leg_length\" set to [km]")
      } else {
        effort_base$leg_length <- units::set_units(effort_base$leg_length, m)
        ui_done("\"leg_length\" set to [m]")
      }
    } else {
      effort_base$leg_length <- set_units(effort_base$leg_length, km)
      ui_done("\"leg_length\" set to [km]")
    }
    if (!("seg_id" %in% colnames(effort_base))) {
      colnames(effort_base)[lower_no_under_colnames %in% c("segid")] <- "seg_id"
    }
    if (!("seg_length" %in% colnames(effort_base))) {
      colnames(effort_base)[lower_no_under_colnames %in% c("seglengthkm","seglengkm")] <- "seg_length"
    }
    if(!any(is(effort_base$seg_length)=="units")) {
      if(unit_km == TRUE) {
        effort_base$seg_length <- units::set_units(effort_base$seg_length, km)
        ui_done("\"seg_length\" set to [km]")
      } else {
        effort_base$seg_length <- units::set_units(effort_base$seg_length, m)
        ui_done("\"seg_length\" set to [m]")
      }
    } else {
      effort_base$seg_length <- units::set_units(effort_base$seg_length, km)
      ui_done("\"seg_length\" set to [km]")
    }

    # Replace facultative colnames ----
    if (!("session" %in% colnames(effort_base))) {
      colnames(effort_base)[lower_no_under_colnames %in% c("session")] <- "session"
    }
    if (!("center_time" %in% colnames(effort_base))) {
      colnames(effort_base)[lower_no_under_colnames %in% c("centertime")] <- "center_time"
    }

    # Build region_label if it doesn't exist ----
    if (!("region_label" %in% colnames(effort_base))) {
      effort_base$region_label <- paste(effort_base$sub_region, effort_base$strate, sep = "_")
    }

    # Build center time : divise l'interval de temps du leg en 2, et center time = la limite de droite
    if (!("centertime" %in% lower_no_under_colnames) &
        ("dtstart" %in% lower_no_under_colnames) &
        ("dtend" %in% lower_no_under_colnames)) {

      df_effort_base <- effort_base %>%
        as.data.frame()

      effort_base$center_time <- lubridate::int_end(
        lubridate::interval(
          start = ymd_hms(
            df_effort_base[,lower_no_under_colnames %in% c("dtstart")]
          ),
          end = ymd_hms(
            df_effort_base[,lower_no_under_colnames %in% c("dtend")]
          )
        ) /2
      )
    }

    # If multi_survey == T, ajoute le nom de campagne dans les id region, leg et seg, session
    if(multi_survey == TRUE){
      effort_base$region_label <- paste(effort_base$survey, effort_base$region_label, sep = "_")
      effort_base$leg_id <- paste(effort_base$survey, effort_base$leg_id, sep = "_")
      effort_base$seg_id <- paste(effort_base$survey, effort_base$seg_id, sep = "_")
      if ("session" %in% colnames(effort_base)) {
        effort_base$session <- paste(effort_base$survey, effort_base$session, sep = "_")
      }
    }

    # Error message if needed columns are missing ----
    if (!all(col_name_neces %in% colnames(effort_base))) {
      var_still_alone <- col_name_neces[!(col_name_neces %in% colnames(effort_base))]
      stop(paste(c("For effort data, can't find equivalent name for : ", var_still_alone), collapse="\n"))
    }

    # # Delete "lon" and "lat" columns (we build them with prepare_functions)
    # if("lon" %in% colnames(effort_base)) {
    #   effort_base <- effort_base %>% dplyr::select(-lon) # ne modifie pas dans gitlab, car normalement les packages annexes sont chargés quand charge pelaDSM (ce qui ne marche pas ici pour le moment
    # }
    # if("lat" %in% colnames(effort_base)) {
    #   effort_base <- effort_base %>% dplyr::select(-lat)
    # }

    class(effort_base) <- c("cleaned_effort_dsm","sf","data.frame")
    return(effort_base)
  }

}


#' Standardization of observation data column names
#'
#' Attribute standard names to column names of observation sf object
#' for the next steps of the analysis
#'
#' @param obs_base sf object containing column names to standardize
#' @param multi_survey logical, whether several surveys are included in "obs_base"
#' @param unit_km Is \code{distance} column of "obs_base" in km ?
#' \itemize{
#'   \item If \code{TRUE}, it is considered in *km* and set in *km* with \code{units::set_units()}
#'   \item If \code{FALSE}, it is considered in *m* and set to *m* with \code{units::set_units()}
#' }
#' @return This function return the same data.frame as in input with the standardized column name of class *cleaned_obs_dsm*
#' @import sf assertthat
#' @importFrom units set_units
#' @importFrom usethis ui_done
#' @export
#' @examples
#' library(pelaDSM)
#'
#' raw_obs <- SPEE_data_dsm$obs_with_sp
#'
#' cleaned_obs <- clean_obs_dsm(
#'   obs_base = raw_obs,
#'   unit_km = FALSE # "distance" column of raw_obs is in meter
#' )
clean_obs_dsm <- function(obs_base, multi_survey = FALSE, unit_km = TRUE) {

  assert_that(any(is(obs_base) == "sf"))

  col_name_neces <- c("strate","sub_region","transect","taxon","group",
                      "region_label","family","species","pod_size","dec_angle",
                      "leg_id","seg_id","distance", "survey")

  if(all(col_name_neces %in% colnames(obs_base))){

    return(obs_base)

  } else {

    # Put everything tolowercase and without "_" and check if there are matches with colnames ----
    lower_colnames <- tolower(colnames(obs_base))
    lower_no_under_colnames <- gsub("_", "", lower_colnames)
    lower_neces <- tolower(col_name_neces)

    for (n in 1:length(lower_neces)) {
      if (lower_neces[n] %in% lower_no_under_colnames) {
        colnames(obs_base)[lower_no_under_colnames %in% lower_neces[n]] <- col_name_neces[n]
      }
    }

    # Are there all the needed columns ? ----
    missing_needed_col <- which(!(col_name_neces %in% colnames(obs_base)))

    if (!("region_label" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("regionlabel","stratesec")] <- "region_label"
    }
    if (!("sub_region" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("subregion","secteur")] <- "sub_region"
    }
    if (!("taxon" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("taxon","taxonfr")] <- "taxon"
    }
    if (!("group" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("group","groupe","groupefr")] <- "group"
    }
    if (!("family" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("family","famille","famillefr")] <- "family"
    }
    if (!("species" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("species","codespecies","codesp","speciescode")] <- "species"
    }
    if (!("pod_size" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("podsize","size")] <- "pod_size"
    }
    if (!("dec_angle" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("decangle")] <- "dec_angle"
    }
    if (!("leg_id" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("legid")] <- "leg_id"
    }
    if (!("distance" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("perpdist")] <- "distance"
    }
    if(!any(is(obs_base$distance)=="units")) {
      if(unit_km == TRUE) {
        obs_base$distance <- units::set_units(obs_base$distance, km)
        ui_done("\"distance\" set to [km]")
      } else {
        obs_base$distance <- units::set_units(obs_base$distance, m)
        ui_done("\"distance\" set to [m]")
      }
    } else {
      obs_base$distance <- units::set_units(obs_base$distance, km)
      ui_done("\"distance\" set to [km]")
    }
    if (!("seg_id" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("segid")] <- "seg_id"
    }

    # Replace facultative colnames ----
    if (!("session" %in% colnames(obs_base))) {
      colnames(obs_base)[lower_no_under_colnames %in% c("session")] <- "session"
    }


    # Build region_label if it doesn't exist ----
    if (!("region_label" %in% colnames(obs_base))) {
      obs_base$region_label <- paste(obs_base$sub_region, obs_base$strate, sep = "_")
    }

    # If multi_survey == T, ajoute le nom de campagne dans les id region, leg, session
    if(multi_survey == TRUE){
      obs_base$region_label <- paste(obs_base$survey, obs_base$region_label, sep = "_")
      obs_base$leg_id <- paste(obs_base$survey, obs_base$leg_id, sep = "_")
      obs_base$seg_id <- paste(obs_base$survey, obs_base$seg_id, sep = "_")
      if ("session" %in% colnames(obs_base)) {
        obs_base$session <- paste(obs_base$survey, obs_base$session, sep = "_")
      }
    }
    # Transform observer into observerId ----
    if("observer" %in% colnames(obs_base)) {
      obs_base$observer_id <- obs_base$observer
      obs_base$observer <- NULL
    }
    # Error message if needed columns are missing ----
    if (!all(col_name_neces %in% colnames(obs_base))) {
      var_still_alone <- col_name_neces[!(col_name_neces %in% colnames(obs_base))]
      stop(paste(c("For observation data, can't find equivalent name for : ", var_still_alone), collapse="\n "))
    }

    # # Delete "lon" and "lat" columns (we build them with prepare_functions)
    # if("lon" %in% colnames(obs_base)) {
    #   obs_base <- obs_base %>% dplyr::select(-lon)
    # }
    # if("lat" %in% colnames(obs_base)) {
    #   obs_base <- obs_base %>% dplyr::select(-lat)
    # }

    class(obs_base) <- c("cleaned_obs_dsm","sf","data.frame")
    return(obs_base)
  }
}
