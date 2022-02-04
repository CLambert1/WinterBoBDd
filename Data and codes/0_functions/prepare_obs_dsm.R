#' @export
prepare_obs <- function(x, ...) {
  UseMethod("prepare_obs")
}

#' Preparation of observation data
#'
#' Transform raw observation data into multiple sub data.frame for next analysis of
#' other functions of the package
#'
#' @param sp name of taxon, group, family or species
#' @param obs_base observation database \code{sf} object
#' @param segdata segdata data.frame, built with \code{\link{prepare_effort_dsm}}
#' @param projection coordinate reference system: integer with the EPSG code, or character with proj4string
#' @param predictors covariates to keep in distdata. It should be a \code{vector} of columns
#' @param truncation distance of truncation of the sampling. If \code{NULL},
#' a distance of truncation will be set, defined by the 95th quantile of \code{distance} column.
#' Default = \code{NULL}.
#' @param unit_truncation_km Default = \code{FALSE}. Is the unit of truncation set in km ?
#' @param quiet logical, wether to print information on species, strip-transcet, units of leg and seg
#' @return This function return a list containing :
#'         \enumerate{
#'           \item distdata : \code{data.frame} with distance information, which aim to be an input for
#'           \code{\link[pelaDSM]{cds_detection4dsm}}.
#'           \item obsdata : Data.frame containing information at observation scale.
#'           \item countdata_leg : A data.frame that merge leg scale effort informations
#'           with number of sightings and number of individuals (n_detection and n_ind).
#'           \item countdata_seg : A data.frame that merge segment scale effort informations
#'           with number of sightings and number of individuals (n_detection and n_ind).
#'         }
#' @import dplyr sf cli
#' @importFrom usethis ui_done
#' @importFrom units deparse_unit set_units
#' @export
#' @examples
#' library(pelaDSM)
#'
#' cleaned_obs <- SPEE_data_dsm$standard_obs # output of clean_obs_dsm()
#'
#' segdata <- SPEE_data_dsm$effort_output$segdata # output of prepare_effort_dsm()
#'
#' # strip-transect in case of Jellyfish, all value of distance are already set to 0.2
#' observation_output_JELLY <- prepare_obs_dsm(
#'   sp = "JELLY",
#'   obs_base = cleaned_obs,
#'   truncation = NULL,         # strip-transect set distance to 0.2
#'   segdata = segdata,
#'   projection = 4326,        # WGS84 projection of lon and lat
#'   unit_truncation_km = TRUE # truncation is in [km]
#' )
#'
#' str(observation_output_JELLY, max.level = 2)
#' hist(observation_output_JELLY$distdata$distance,
#'      col = "grey",
#'      breaks = seq(0,0.5,0.05),
#'      xlab = paste0(
#'        "distance in [",
#'        units(observation_output_JELLY$distdata$distance)[1],
#'        "]"
#'      ),
#'      main = "Distance histogram for JELLY")
#'
#' # observation distance are sampled for Sunfish
#' observation_output_MOLMOL <- prepare_obs_dsm(
#'   sp = "MOLMOL",
#'   obs_base = cleaned_obs,
#'   truncation = 0.3, # truncation set to 0.3 km
#'   segdata = segdata,
#'   projection = 4326, # WGS84 projection of lon and lat
#'   unit_truncation_km = TRUE
#' )
#'
#' str(observation_output_MOLMOL, max.level = 2)
#' hist(observation_output_MOLMOL$distdata$distance,
#'      col = "grey",
#'      breaks = seq(0,0.5,0.05),
#'      xlab = paste0(
#'        "distance in [",
#'        units(observation_output_MOLMOL$distdata$distance)[1],
#'        "]"
#'      ),
#'      main = "Distance histogram for MOLMOL")
prepare_obs.cleaned_obs_dsm <- function(obs_base,
                                        sp,
                                        segdata,
                                        projection,
                                        predictors = NULL,
                                        truncation = NULL,
                                        unit_truncation_km = TRUE,
                                        quiet = FALSE) {

  # Declare species
  if(quiet == FALSE) {
    cli::cli_h1("{.sp {sp}}")
  }


  # Prepare sp_data filtering -----------------------------------------------

  # Create "lon" and "lat" from geometry and projection arg
  raw_obs <- obs_base %>%
    st_transform(crs = projection) %>%
    as_Spatial() %>%
    suppressWarnings() %>%
    as.data.frame() #%>%
    # rename(lon = coords.x1, lat = coords.x2)


  # check if all necessary columns are present
  col_name_neces <- c("survey","strate","sub_region","region_label","transect","leg_id","seg_id",
                      "side", "taxon","group", "family","species","pod_size","dec_angle","distance")

  if(!all(col_name_neces %in% colnames(raw_obs))){
    var_alone <- setdiff(col_name_neces, colnames(raw_obs))

    stop(paste("Variables: ",var_alone, "are not in obs. table", sep="\n",collapse = ", "),
         paste("use clean_obs_name function","\n",sep=""))
  }


  # remove center observation
  raw_obs <- subset(raw_obs, side != "CENTER")


  # ne prendre que l'espece choisie ###TEST###
  # cas pour plusieurs espèces en même temps

  sp_data <- raw_obs[0,]

  if (any(sp %in% unique(raw_obs$group))) {
    match_group <- sp[which(sp %in% unique(raw_obs$group))]
    sp_data_group <- subset(raw_obs, group %in% match_group)
    sp_data <- rbind(sp_data, sp_data_group)
  }
  if (any(sp %in% unique(raw_obs$taxon))) {
    match_taxon <- sp[which(sp %in% unique(raw_obs$taxon))]
    sp_data_taxon <- subset(raw_obs, taxon %in% match_taxon)
    sp_data <- rbind(sp_data, sp_data_taxon)
  }
  if (any(sp %in% unique(raw_obs$family))) {
    match_family <- sp[which(sp %in% unique(raw_obs$family))]
    sp_data_family <- subset(raw_obs, family %in% match_family)
    sp_data <- rbind(sp_data, sp_data_family)
  }
  if (any(sp %in% unique(raw_obs$species))) {
    match_species <- sp[which(sp %in% unique(raw_obs$species))]
    sp_data_species <- subset(raw_obs, species %in% match_species)
    sp_data <- rbind(sp_data, sp_data_species)
  }


  # ne garder que les obs dans la bande pour les oiseaux
  if(all(sp_data$taxon %in% c("Oiseau marin","Seabird")) && quiet == FALSE) {
    if(quiet == FALSE) {
      cli_alert_info("Keeping only observations within the 200 m around the transect")
    }
    sp_data <- subset(sp_data, dec_angle %in% c(1, 3))
  }

  # convert distance in segdata effort unit
  if(any(units(segdata$Effort) %in% "km")) {

    sp_data$distance <- sp_data$distance %>%
      units::set_units(km)

    if(quiet == FALSE) {
      usethis::ui_done("\"distance\" set to [km] as \"Effort\" is set to [km]")
    }

  } else {
    sp_data$distance <- sp_data$distance %>%
      units::set_units(m)

    if(quiet == FALSE) {
      usethis::ui_done("\"distance\" set to [m] as \"Effort\" is set to [m]")
    }
  }

  # truncation
  if(is.null(truncation)) {

    # set truncation to 0.95 quantile
    pas <- 0.05
    dist <- sp_data$distance %>%
      units::set_units(km)
    wa <- as.numeric(quantile(dist,  prob = 0.95))
    wa <- pas*floor(wa/pas) + ifelse(wa %% pas == 0, 0, pas)
    wa <- set_units(wa, km)

  } else {

    if(all(is(truncation) != "units")) {

      if(unit_truncation_km == TRUE) {
        wa <- set_units(truncation, km)
      } else if(unit_truncation_km == FALSE) {
        wa <- set_units(truncation, m)
      }

    } else {
      wa <- truncation
    }

  }

  # convert wa in distance units
  if(units::deparse_unit(segdata$Effort) == "km") {
    wa <- set_units(wa, km)
    if(quiet == FALSE) {
      usethis::ui_done("\"truncation\" set to [km] as \"Effort\" is set to [km]")
    }
  } else {
    wa <- set_units(wa, m)
    if(quiet == FALSE) {
      usethis::ui_done("\"truncation\" set to [m] as \"Effort\" is set to [m]")
    }
  }

  sp_data <- subset(sp_data, distance <= wa)





  # Build countdata_leg and countdata_seg -----------------------------------

  if("session" %in% colnames(sp_data)) {
    countdata_seg <- as.data.frame(sp_data[, c("transect", "leg_id", "seg_id","pod_size","session")] %>%
                                     group_by(transect, leg_id, seg_id) %>%
                                     summarize(n = n(),
                                               count = sum(pod_size),
                                               .groups = "drop"))
  } else {
    countdata_seg <- as.data.frame(sp_data[, c("transect", "leg_id", "seg_id","pod_size")] %>%
                                     group_by(transect, leg_id, seg_id) %>%
                                     summarize(n = n(),
                                               count = sum(pod_size),
                                               .groups = "drop"))
  }

  colnames(countdata_seg) <- c("Transect.Label", "leg_id", "Sample.Label","n_detection", "n_ind")

  countdata_leg <- as.data.frame(countdata_seg %>%
                                   group_by(Transect.Label, leg_id) %>%
                                   summarize(n_detection = sum(n_detection), n_ind = sum(n_ind)))




  # Build distdata with truncation ------------------------------------------

  if("session" %in% colnames(sp_data)) {
    distdata <- sp_data[, c("transect", "region_label", "leg_id", "seg_id",
                            "pod_size", "distance", "observer_id","session",
                            "lon", "lat")]
    names(distdata) <- c("Transect.Label", "Region.Label", "leg_id", "Sample.Label",
                         "size", "distance","observer_id","session",
                         "lon", "lat")
  } else {
    distdata <- sp_data[, c("transect", "region_label", "leg_id", "seg_id",
                            "pod_size", "distance", "observer_id", "lon", "lat")]
    names(distdata) <- c("Transect.Label", "Region.Label", "leg_id", "Sample.Label",
                          "size", "distance","observer_id","lon","lat")
  }


  distdata$object <- as.numeric(row.names(tibble(distdata)))
  distdata$detected <- 1
  distdata <- subset(distdata, Sample.Label %in% segdata$Sample.Label) # passe au segdata pour ajuster detection sur les segments pour le dsm


  if("session" %in% colnames(sp_data)) {
    distdata <- left_join(dplyr::select(segdata, date:session, Area, all_of(predictors)), # enl?ve les pr?dictors de segdata
                          dplyr::select(distdata, -Transect.Label, -Region.Label, -session, -leg_id, -lon, -lat),
                          by="Sample.Label")
  } else {
    distdata <- left_join(dplyr::select(segdata, date:subjective, Area, all_of(predictors)),
                          dplyr::select(distdata, -Transect.Label, -Region.Label, -leg_id, -lon, -lat),
                          by="Sample.Label")
  }


  distdata$detected[is.na(distdata$detected)] <- 0 # à verifier




  # Build obsdata ---------------------------------------------------------

  if("session" %in% colnames(sp_data)) {
    obsdata <- subset(distdata, detected == 1)[, c("distance", "size", "Transect.Label", "Region.Label",
                                                   "Sample.Label","observer_id","session")]
  } else {
    obsdata <- subset(distdata, detected == 1)[, c("distance", "size", "Transect.Label", "Region.Label",
                                                   "Sample.Label","observer_id")]
  }
  distObject <- distdata$object
  obsdata$object <- distObject[!is.na(distObject)]


  class(distdata) <- c("distdata_dsm","data.frame")

  return(list(distdata = distdata,
              obsdata = obsdata,
              countdata_leg = countdata_leg,
              countdata_seg = countdata_seg,
              trunc = wa))
}


