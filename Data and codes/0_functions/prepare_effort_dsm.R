#' prepare effort for dsm cleaned data
#'
#' @param effort_base sf object with standardized columns effort data with [clean_effort_dsm].
#' @param predictor \code{vector} of segment-scale predictor names to keep in output of the function.
#' @param block_area \code{data.frame} with 2 colnames :
#'                     \itemize{
#'                       \item{\code{Block} :}{names that matches region_label}
#'                       \item{\code{Area} :}{area of each block in *km^2* or *m^2* can have units set}
#'                     }
#' @param unit_Area_km  Is the unit of \code{Area} in "block_area" set in km ?
#' \itemize{
#'   \item If \code{TRUE}, it is considered in *km^2* and set to *km^2* with \code{units::set_units()}
#'   \item If \code{FALSE}, it is considered in *m^2* and set to *m^2* with \code{units::set_units()}
#' }
#' @param projection coordinate reference system: integer with the EPSG code, or character with proj4string
#' @param optimal Argument which allows to keep data sampled in optimal conditions.
#'        Default setting is only data sampled in good conditions are kept (\code{optimal = T}), indexes
#'        \code{"c("GG", "GM", "MG", "EG", "GE", "EE", "ME", "EM", "MM")"}
#' @param col2keep \code{vector} of leg-scale column names to keep from "effort_base" in output of the function.
#' @return This function return a list containing:
#'         \enumerate{
#'           \item legdata : \code{data.frame} with infos at leg scale.
#'           \item segdata : \code{data.frame} with infos at segment scale.
#'         }
#' @import dplyr units
#' @importFrom lubridate ymd
#' @importFrom pelaCDS prepare_effort
#' @method prepare_effort cleaned_effort_dsm
#' @export
#' @examples
#' library(pelaDSM)
#'
#' Block_area_sq_km <- data.frame(
#'   Block = "ATL_N1",
#'   Area = 14947.03 # no unit but it is in [km^2]
#' )
#'
#' prepared_effort <- prepare_effort_dsm(
#'   effort_base = SPEE_data_dsm$standard_effort,
#'   predictor = c("SST_month","depth"),
#'   block_area = Block_area_sq_km,
#'   unit_Area_km = TRUE, # Block_area_sq_km$Area is in [km^2]
#'   projection = 4326,
#'   optimal = TRUE,      # we want to keep observation sampled in good condition
#'   col2keep = NULL
#' )
prepare_effort.cleaned_effort_dsm <- function(effort_base,
                                              predictor = NULL,
                                              block_area,
                                              unit_Area_km = TRUE,
                                              projection,
                                              optimal = TRUE,
                                              col2keep = NULL,
                                              quiet = FALSE) {


  effort <- effort_base %>%
    st_transform(crs = projection) %>%
    # st_centroid() %>%
    as_Spatial() %>%
    suppressWarnings() %>%
    as.data.frame() #%>%
    # rename(lon = coords.x1, lat = coords.x2)

  # verifier si les colonnes sont bien dans le DF effort
  col_name_neces <- c("sea_state","subjective","lon","lat",
                      "survey","region_label","transect","leg_id","leg_length",
                      "seg_id","seg_length",
                      "left","right","center_time")


  if(!all(col_name_neces %in% colnames(effort))){
    var_alone <- setdiff(col_name_neces, colnames(effort))

    stop(paste("Variables: ",var_alone, "are not in the effort table.", sep="\n",collapse = ", "),
         paste("use clean_effort_name function","\n",sep=""))
  }


  ### Block et Surface
  # bon nom pour block et area
  assert_that(is.data.frame(block_area))
  assert_that(block_area %has_name% c("Block","Area"))

  # correspondance entre region_label et block area
  if(!any(block_area$Block %in% effort$region_label)){
    stop(paste("The values in \"Block\" variable from \"block_area\" are different
               from \"region_label\" values of \"effort_base\":",
               "- \"region_label\" created with clean_effort_dsm have incorrect format",
               "- \"Block\" values are not in \"region_label\" of \"effort_base\"",sep="\n"))
  }

  # units of block_area$Area
  if(!any(is(block_area$Area) == "units")) {

    if(unit_Area_km == TRUE) {

      block_area$Area <- units::set_units(block_area$Area, "km^2")

    } else {

      block_area$Area <- units::set_units(block_area$Area, "m^2")

    }

  }

  if(units::deparse_unit(effort_base$seg_length) == "km") {

    block_area$Area <- units::set_units(block_area$Area, "km^2")

    if(quiet == FALSE) {
      ui_done("\"Area\" set to [km^2] as \"seg_length\" and \"leg_length\" are set to [km]")
    }

  } else {

    block_area$Area <- units::set_units(block_area$Area, "m^2")

    if(quiet == FALSE) {
      ui_done("\"Area\" set to [m^2] as \"seg_length\" and \"leg_length\" are set to [m]")
    }
  }


  # predictor
  allvar = predictor


  # selection effort et obs en bonnes conditions
  if(optimal==T) {
    effort <- subset(effort, sea_state <= 3 & subjective %in% c("GG", "GM", "MG", "EG", "GE", "EE", "ME", "EM", "MM"))
  }


  # creation des tableaux necessaires a l'analyse  ---------------

  # est-ce la peine de conserver legdata ?? oui au cas où

  #-- legdata --#
  #-------------#
  if("session" %in% colnames(effort)) {
    legdata <- effort %>%
      group_by(survey, region_label, transect, leg_id, left, right, session, leg_length, sea_state, subjective) %>%
      summarise(.groups = "drop") %>%
      rename(Effort = leg_length)
  } else {
    legdata <- effort %>%
      group_by(survey, region_label, transect, leg_id, left, right, leg_length, sea_state, subjective) %>%
      summarise(.groups = "drop") %>%
      rename(Effort = leg_length)
  }

  legdata <- as.data.frame(legdata)

  names(legdata)[which(names(legdata) %in% c("region_label", "transect", "leg_id"))] <- c("Region.Label", "Transect.Label",
                                                                                          "leg_id")

  # merge col2keep
  if(!is.null(col2keep)){

    col2keep <- col2keep[!(col2keep %in% colnames(legdata))]

    legdata_col_joined <- legdata %>%
      left_join(effort %>%
                  select(c(all_of(col2keep),"leg_id")) %>%
                  unique(), # remove duplicated lines
                by = "leg_id")

    if(nrow(legdata) < nrow(legdata_col_joined)){
      stop("col2keep are not unique at leg scale")
    }

    legdata <- legdata_col_joined

  }

  # Assigner area à legdata en fonction du nom du block commun avec block_area
  legdata <- legdata %>%
    left_join(block_area, by = c("Region.Label"="Block"))




  #-- segdata --#
  #-------------#
  if("session" %in% colnames(effort)) {
    segdata <- data.frame(effort[, c("center_time", "survey", "transect", "leg_id", "seg_id", "left", "right", # add left and right to keep observer info in segdata
                                     "seg_length", "lon", "lat", "region_label",
                                     "sea_state", "subjective","session", allvar, col2keep)
    ])
  } else {
    segdata <- data.frame(effort[, c("center_time", "survey", "transect", "leg_id", "seg_id","left", "right",
                                     "seg_length", "lon", "lat", "region_label",
                                     "sea_state", "subjective", allvar, col2keep)
    ])
  }

  segdata$center_time <- ymd_hms(segdata$center_time)
  names(segdata)[1:11] <- c("date", "survey", "Transect.Label", "leg_id", "Sample.Label", "left", "right",
                            "Effort", "lon", "lat", "Region.Label")

  # Assigner area à segdata en fonction du nom du block commun avec block_area
  # convert Area in [km] or [m] depending on Effort unit
  segdata <- segdata %>%
    left_join(block_area, by = c("Region.Label"="Block"))




  # renvoyer les outputs
  return(list(legdata = legdata, segdata = segdata))

}
