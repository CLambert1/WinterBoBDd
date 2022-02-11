##--------------------------------------------------------------------------------------------------------
## PROJECT : Infra-seasonal variations of common dolphin distribution in the Bay of Biscay
## SCRIPT : Prepare data analysis (format tables, extract covariates)
##
## Author : Charlotte Lambert
## Last update : Feb 2022
## R version 4.0.4 (2021-02-15) -- "Lost Library Book"
##--------------------------------------------------------------------------------------------------------


library(rgdal)
library(rgeos)
options(stringsAsFactors = FALSE)
library(dplyr)
library(tidyverse)
library(raster)
library(maps)
library(maptools)
library(ggplot2)
library(ggthemes)
library(ggspatial)
library(ggpubr)
library(RColorBrewer)
library(reshape2)
library(lubridate)
library(Distance)
library(coda)
library(dsm)
library(ggdist)
library(pelaCDS)
library(lubridate)
library(pelaStan)
library(assertthat)
library(sf)
library(birk)


#### ----------------------- ####
#### ---- SOURCES FILES ---- ####
#### ----------------------- ####
outdir <- "~/3_script_results"
sourcedir <- "~/1_data"
funcdir <- "~/0_functions"
setwd(outdir)

L93 <- 2154
WGS <- 4326


### 1 - Load land, isobath -----------------------------------------------------
## A - Open, clip and project land
land <- st_read(dsn = paste(sourcedir, "Map_files", sep = "/"), layer = "Land")
bbx <- as(raster::extent(-4.8,0.5,44.3,47.1), "SpatialPolygons")
raster::crs(bbx) <- raster::crs(land)
bbx_L93 <- bbx %>% st_as_sf() %>% st_transform(crs = st_crs(L93)) 
land_L93 <- land %>% st_transform(crs = st_crs(L93))


## B - Open, clip and project isob
isob <- st_read(dsn = paste(sourcedir, "Map_files", sep = "/"), layer = "Isobath")
isob_L93 <- isob %>% st_transform(crs = st_crs(L93))


### 2 - Open CAPECET files -----------------------------------------------------
## seg
lin_cap <- st_read(dsn = paste(sourcedir, "Surveys", sep = "/"), layer = "CAPECET_Seg5km_wgs84")
lin_cap$day <- date(lin_cap$date) 
lin_cap$week <- week(lin_cap$day)
lin_cap$year <- year(lin_cap$date)
lin_cap$month <- month(lin_cap$datehhmmss)
lin_cap$YearWeek <- str_c(lin_cap$year,lin_cap$week,sep = "_")

## obs
obs_cap <- st_read(dsn = paste(sourcedir, "Surveys", sep = "/"), layer = "CAPECET_Obs_wgs84")
obs_cap$lon <- parse_number(obs_cap$lon)
obs_cap$day <- date(obs_cap$date) 
obs_cap$week <- week(obs_cap$day)
obs_cap$year <- year(obs_cap$date)
obs_cap$month <- month(obs_cap$datehhmmss)
obs_cap$YearWeek <- str_c(obs_cap$year,obs_cap$week,sep = "_")


### 3 - Open SPEE files --------------------------------------------------------
## seg
lin_spee <- st_read(dsn = paste(sourcedir, "Surveys", sep = "/"), layer = "SPEE_Seg5km_wgs84")
lin_spee$day <- date(lin_spee$date) 
lin_spee$week <- week(lin_spee$day)
lin_spee$year <- year(lin_spee$date)
lin_spee$month <- month(lin_spee$datehhmmss)
lin_spee$YearWeek <- str_c(lin_spee$year,lin_spee$week,sep = "_")

## obs
obs_spee <- st_read(dsn = paste(sourcedir, "Surveys", sep = "/"), layer = "SPEE_Obs_wgs84")
obs_spee$day <- date(obs_spee$date) 
obs_spee$week <- week(obs_spee$day)
obs_spee$year <- year(obs_spee$date)
obs_spee$month <- month(obs_spee$datehhmmss)
obs_spee$YearWeek <- str_c(obs_spee$year,obs_spee$week,sep = "_")


### 4 - Open SAMM2 files -------------------------------------------------------
## seg
lin_samm <- st_read(dsn = paste(sourcedir, "Surveys", sep = "/"), layer = "SAMM_Seg5km_wgs84") %>%
  rename(Date_Time1 = DtTmBgn, datehhmmss = DtTmCnt, DATE_TIME2 = DtTmEnd, 
         routeType = routTyp, strateType = strtTyp, altitude = altitud, subRegion = subRegn,
         seaState = seastat, turbidity = turbdty, skyGlint = skyGlnt, glareFrom = glarFrm, glareSever = glarSvr, 
         glareUnder = glrUndr, cloudCover = clodCvr, subjective = subjctv, transect = IdTrnsc, 
         legLengKm = LgLngtK, segLengKm = SgLngtK, lon = Longitd, lat = Latitud, legID = LegID, segID = SegID)
xx <- st_within(x = st_transform(lin_samm, crs = st_crs(L93)), y = bbx_L93, sparse = FALSE) 
lin_samm <- lin_samm[xx[,1] == T,]
lin_samm$day <- date(lin_samm$datehhmmss) 
lin_samm$week <- week(lin_samm$day)
lin_samm$year <- year(lin_samm$datehhmmss)
lin_samm$month <- month(lin_samm$datehhmmss)
lin_samm$YearWeek <- str_c(lin_samm$year,lin_samm$week,sep = "_")
lin_samm$center_time <- lin_samm$datehhmmss

## obs
obs_samm <- st_read(dsn = paste(sourcedir, "Surveys", sep = "/"), layer = "SAMM_Obs_wgs84") %>%
  subset(species %in% c("STEDEL", "DELDEL", "SMADEL")) %>% 
  rename(DateTimeBegin = DtTmBgn, DateTimeCentre = DtTmCnt, DateTimeEnd = DtTmEnd, altitude = altitud,
         transect = IdTrnsct_x, subRegion = subRegn, strateType = strtTyp, computer = computr,
         routeType = routTyp, sighting = sightng, decAngle= decAngl, behaviour = behavir,
         gpsDelay = gpsDely, aircraft = aircrft, datehhmmss = dthhmms, perpDist = perpDst,
         IdTransect_y = IdTrnsct_y, observer = observr, legID = LegID, segID = SegID) 
xx <- st_within(x = st_transform(obs_samm, crs = st_crs(L93)), y = bbx_L93, sparse = FALSE)
obs_samm <- obs_samm[xx[,1] == T,]
obs_samm$day <- date(obs_samm$date)
obs_samm$week <- week(obs_samm$day)
obs_samm$year <- year(obs_samm$date)
obs_samm$month <- month(obs_samm$datehhmmss)
obs_samm$YearWeek <- str_c(obs_samm$year,obs_samm$week,sep = "_")
obs_samm$center_time <- obs_samm$DateTimeCentre


### 5 - Open SPEE9 files -------------------------------------------------------
## seg
lin_spee9 <- st_read(dsn = paste(sourcedir, "Surveys", sep = "/"), layer = "SPEE9_Seg5km_wgs84") %>%
  rename(Date_Time1 = DtTmBgn, datehhmmss = DtTmCnt, DATE_TIME2 = DtTmEnd,  
         routeType = routTyp, strateType = strtTyp, altitude = altitud, subRegion = subRegn,
         seaState = seastat, turbidity = turbdty, skyGlint = skyGlnt, glareFrom = glarFrm, glareSever = glarSvr, 
         glareUnder = glrUndr, cloudCover = clodCvr, subjective = subjctv, transect = IdTrnsc, 
         legLengKm = LgLngtK, segLengKm = SgLngtK, lon = Longitd, lat = Latitud, legID = LegID, segID = SegID)
lin_spee9$day <- date(lin_spee9$datehhmmss) 
lin_spee9$week <- week(lin_spee9$day)
lin_spee9$year <- year(lin_spee9$datehhmmss)
lin_spee9$month <- month(lin_spee9$datehhmmss)
lin_spee9$YearWeek <- str_c(lin_spee9$year,lin_spee9$week,sep = "_")
lin_spee9$center_time <- lin_spee9$datehhmmss

## obs
obs_spee9 <- st_read(dsn = paste(sourcedir, "Surveys", sep = "/"), layer = "SPEE9_Obs_wgs84") %>%
  subset(species %in% c("STEDEL", "DELDEL", "SMADEL")) %>% 
  rename(DateTimeBegin = DtTmBgn, DateTimeCentre = DtTmCnt, DateTimeEnd = DtTmEnd, altitude = altitud,
         transect = IdTrnsct_x, subRegion = subRegn, strateType = strtTyp, computer = computr,
         routeType = routTyp, sighting = sightng, decAngle= decAngl, behaviour = behavir,
         gpsDelay = gpsDely, aircraft = aircrft, datehhmmss = dthhmms, datehhmmss_adjusted = dthhmm_,
         perpDist = perpDst, IdTransect_y = IdTrnsct_y, observer = observr, legID = LegID, segID = SegID)
obs_spee9$day <- date(obs_spee9$date)
obs_spee9$week <- week(obs_spee9$day)
obs_spee9$year <- year(obs_spee9$date)
obs_spee9$month <- month(obs_spee9$datehhmmss)
obs_spee9$YearWeek <- str_c(obs_spee9$year,obs_spee9$week,sep = "_")
obs_spee9$center_time <- obs_spee9$DateTimeCentre


### 6 - Merge all surveys ------------------------------------------------------
surv_all <- bind(as_Spatial(lin_cap), as_Spatial(lin_spee), as_Spatial(lin_samm), as_Spatial(lin_spee9))
surv_all$YearWeek <- factor(surv_all$YearWeek,
                            levels = c("2020_2", "2020_4", "2020_6", "2020_7", "2020_8", "2020_11",
                                       "2021_2",  "2021_3", "2021_4", "2021_6", "2021_7", "2021_8", "2021_9", "2021_10", "2021_11", "2021_12"))
surv_all <- surv_all %>% st_as_sf() %>%
  mutate(Month = case_when(
    month == 1 ~ "January", 
    month == 2 ~ "February", 
    month == 3 ~ "March"),
    Period = case_when(year == 2020 ~ "2020", year == 2021 ~ "2021"))
surv_all$Month <- factor(surv_all$Month, levels = c("January", "February", "March"))


obs_all <- bind(as_Spatial(obs_cap), as_Spatial(obs_spee), as_Spatial(obs_samm), as_Spatial(obs_spee9))
obs_all$YearWeek <- factor(obs_all$YearWeek,
                           levels = c("2020_2", "2020_4", "2020_6", "2020_7", "2020_8", "2020_11",
                                      "2021_2",  "2021_3", "2021_4", "2021_6", "2021_7", "2021_8", "2021_9", "2021_10", "2021_11", "2021_12"))
obs_all <- obs_all %>% st_as_sf() %>%
  mutate(Month = case_when(
    month == 1 ~ "January", 
    month == 2 ~ "February", 
    month == 3 ~ "March"),
    Period = case_when(year == 2020 ~ "2020", year == 2021 ~ "2021"))
obs_all$Month <- factor(obs_all$Month, levels = c("January", "February", "March"))


### 7 - Create shapes of surveys -----------------------------------------------
## CAPECET - SPEE - 2019-2020
# CAPECET
shp_capecet <- as(raster::extent(lin_cap), "SpatialPolygons") 
raster::crs(shp_capecet) <- raster::crs(lin_cap)
shp_capecet <- gDifference(shp_capecet, as_Spatial(land), byid = TRUE)
shp_capecet_L93 <- shp_capecet %>% st_as_sf() %>% st_transform(crs = st_crs(L93)) %>% as_Spatial()
shp_capecet$Area <- gArea(shp_capecet_L93)
shp_capecet$Block <- "CAPECET_ATL_N1"


# SPEE
shp_spee <- as(raster::extent(lin_spee), "SpatialPolygons")
raster::crs(shp_spee) <- raster::crs(lin_spee)
shp_spee <- gDifference(shp_spee, as_Spatial(land), byid = TRUE)
shp_spee_L93 <- shp_spee %>% st_as_sf() %>% st_transform(crs = st_crs(L93)) %>% as_Spatial()
shp_spee$Area <- gArea(shp_spee_L93)
shp_spee$Block <- "SPEE_ATL_N1"


## SAMM2 - SPEE9 - 2020-2021
# SAMM2
shp_samm <- as(raster::extent(lin_samm), "SpatialPolygons")
raster::crs(shp_samm) <- raster::crs(lin_samm)
shp_samm <- gDifference(shp_samm, as_Spatial(land), byid = TRUE)
shp_samm_L93 <- shp_samm %>% st_as_sf() %>% st_transform(crs = st_crs(L93)) %>% as_Spatial()
shp_samm$Area <- gArea(shp_samm_L93)
shp_samm$Block <- "SAMM_ATL_N2"


### 8 - Load custom functions --------------------------------------------------
library(pelaCDS)

# these functions are today incorporated into the pelaDSM package
source(paste(funcdir, "clean_names_dsm.R", sep = "/"))
source(paste(funcdir, "prepare_effort_dsm.R", sep = "/"))
source(paste(funcdir, "prepare_obs_dsm.R", sep = "/"))
source(paste(funcdir, "cds_detection4dsm.R", sep = "/"))
source(paste(funcdir, "add_obs_dsm.R", sep = "/"))
library(usethis)
library(units)
library(cli)


#### ---------------------- ####
#### ---- Prepare data ---- ####
#### ---------------------- ####
# Association of region with its surface
Block_Area <- rbind(shp_capecet@data, shp_spee@data, shp_samm@data)


### 1 - Prepare CAPECET data ---------------------------------------------------
# prepare effort
eff_cap <- lin_cap %>%
  mutate(Sample.Label = paste(survey, segID, sep = "_")) %>%
  clean_effort_dsm(., multi_survey = T)
eff_cap2 <- prepare_effort.cleaned_effort_dsm(
  effort_base = eff_cap,
  optimal = T,
  block_area = Block_Area,
  unit_Area_km = FALSE,
  projection = L93)


# prepare obs
obs_cap$taxon <- "Marine Mammal"
obs_cap$group <- "Marine Mammal"
obs_cap$family <- "Cetacea"
sp <- c("DELDEL", "STEDEL")
obs_cap2 <- obs_cap %>%
  mutate(Sample.Label = paste(survey, segID, sep = "_")) %>%
  clean_obs_dsm(., multi_survey = T, unit_km = F)
obs_cap2 <- prepare_obs.cleaned_obs_dsm(
  sp = sp,
  obs_base = obs_cap2,
  truncation = 0.5,
  segdata = eff_cap2$segdata,
  projection = 2154,
  unit_truncation_km = T)


# Add observation on effort
cap <- add_obs_dsm(countdata_seg = obs_cap2$countdata_seg,
                   segdata = eff_cap2$segdata,
                   countdata_leg = obs_cap2$countdata_leg,
                   legdata = eff_cap2$legdata)



### 2 - Prepare SPEE data ------------------------------------------------------
# prepare effort
eff_spee <- lin_spee %>%
  mutate(Sample.Label = paste(survey, segID, sep = "_")) %>%
  clean_effort_dsm(., multi_survey = T)
eff_spee2 <- prepare_effort.cleaned_effort_dsm(
  effort_base = eff_spee,
  optimal = T,
  block_area = Block_Area,
  unit_Area_km = FALSE,
  projection = L93)


# prepare obs
obs_spee$taxon <- "Marine Mammal"
obs_spee$group <- "Marine Mammal"
obs_spee$family <- "Cetacea"
sp <- c("DELDEL", "STEDEL")
obs_spee2 <- obs_spee %>%
  mutate(Sample.Label = paste(survey, segID, sep = "_")) %>%
  clean_obs_dsm(., multi_survey = T, unit_km = F)
obs_spee2 <- prepare_obs.cleaned_obs_dsm(
  sp = sp,
  obs_base = obs_spee2,
  truncation = 0.5,
  segdata = eff_spee2$segdata,
  projection = 2154,
  unit_truncation_km = T)


# Add observation on effort
spee <- add_obs_dsm(countdata_seg = obs_spee2$countdata_seg,
                    segdata = eff_spee2$segdata,
                    countdata_leg = obs_spee2$countdata_leg,
                    legdata = eff_spee2$legdata)



### 3 - Prepare SAMM data ------------------------------------------------------
# prepare effort
eff_samm <- lin_samm %>%
  mutate(Sample.Label = paste(survey, segID, sep = "_")) %>%
  clean_effort_dsm(., multi_survey = T)
eff_samm2 <- prepare_effort.cleaned_effort_dsm(
  effort_base = eff_samm,
  optimal = T,
  block_area = Block_Area,
  unit_Area_km = FALSE,
  projection = L93)


# prepare obs
obs_samm$taxon <- "Marine Mammal"
obs_samm$group <- "Marine Mammal"
obs_samm$family <- "Cetacea"
sp <- c("DELDEL", "STEDEL")
obs_samm2 <- obs_samm %>%
  mutate(Sample.Label = paste(survey, segID, sep = "_")) %>%
  clean_obs_dsm(., multi_survey = T, unit_km = F)
obs_samm2 <- prepare_obs.cleaned_obs_dsm(
  sp = sp,
  obs_base = obs_samm2,
  truncation = 0.5,
  segdata = eff_samm2$segdata,
  projection = 2154,
  unit_truncation_km = T)


# Add observation on effort
samm <- add_obs_dsm(countdata_seg = obs_samm2$countdata_seg,
                    segdata = eff_samm2$segdata,
                    countdata_leg = obs_samm2$countdata_leg,
                    legdata = eff_samm2$legdata)



### 4 - Prepare SPEE9 data -----------------------------------------------------
# prepare effort
eff_spee9 <- lin_spee9 %>%
  mutate(Sample.Label = paste(survey, segID, sep = "_")) %>%
  clean_effort_dsm(., multi_survey = T)
eff_spee92 <- prepare_effort.cleaned_effort_dsm(
  effort_base = eff_spee9,
  optimal = T,
  block_area = Block_Area,
  unit_Area_km = FALSE,
  projection = L93)


# prepare obs
obs_spee9$taxon <- "Marine Mammal"
obs_spee9$group <- "Marine Mammal"
obs_spee9$family <- "Cetacea"
sp <- c("DELDEL", "STEDEL")
obs_spee92 <- obs_spee9 %>%
  mutate(Sample.Label = paste(survey, segID, sep = "_")) %>%
  clean_obs_dsm(., multi_survey = T, unit_km = F)
obs_spee92 <- prepare_obs.cleaned_obs_dsm(
  sp = sp,
  obs_base = obs_spee92,
  truncation = 0.5,
  segdata = eff_spee92$segdata,
  projection = 2154,
  unit_truncation_km = T)


# Add observation on effort
spee9 <- add_obs_dsm(countdata_seg = obs_spee92$countdata_seg,
                     segdata = eff_spee92$segdata,
                     countdata_leg = obs_spee92$countdata_leg,
                     legdata = eff_spee92$legdata)



### 5 - Merge effort tables ----------------------------------------------------
## merge
effort <- rbind(cap$segdata_obs, spee$segdata_obs, samm$segdata_obs, spee9$segdata_obs)

## add ESW
effort <- effort %>% mutate(
  ESW = case_when(
    survey == "CAPECET" ~0.172,
    survey == "SPEE" & year(date) == 2020 ~0.156,
    survey == "SAMM" ~0.198,
    survey == "SPEE" & year(date) == 2021 ~0.181
  )
)



#### -------------------------- ####
#### ---- COVAR EXTRACTION ---- ####
#### -------------------------- ####

### 1 - List of covariates -----------------------------------------------------
## list
vars = c("Chl", "gradChl", "NPP", "gradNPP", "DissIC", "SPCO2", "Phyto", "ZEu", "Temp", "gradTemp", "CurrentSpeed",
         "EKE", "Salinity", "gradSal", "SSH", "gradSSH", "MLD", "gradMLD",
         "gradChl_front_distance", "gradNPP_front_distance", "gradMLD_front_distance", "gradTemp_front_distance")

## open bathy
bathy <- raster(paste(sourcedir, "/Variables/Bathy.tif", sep = ""))


### 2 - Join static variables --------------------------------------------------
## a - Join Bathy
data <- effort
data$Bathy <- raster::extract(x = bathy, y = data[, c("lon", "lat")], method = "simple")


## b - Join yearly static variables
extractStaticVar <- function(effort, dateField, variable, coordLonLat){
  # open rasters
  persist_1920 <- stack(paste(sourcedir, "/Variables/", variable, "_frontal_persistence-2019-12-01_2020-04-13.envi", sep =""))
  avg_dist_1920 <- stack(paste(sourcedir, "/Variables/", variable, "_front_distance-averaged-2019-12-01_2020-04-13.envi", sep =""))
  avg_dist_1920 <- avg_dist_1920 * 202 # there is a typo in the raster processing, this is for correction
  
  persist_2021 <- stack(paste(sourcedir, "/Variables/", variable, "_frontal_persistence-2020-12-01_2021-04-13.envi", sep =""))
  avg_dist_2021 <- stack(paste(sourcedir, "/Variables/", variable, "_front_distance-averaged-2020-12-01_2021-04-13.envi", sep =""))
  avg_dist_2021 <- avg_dist_2021 * 201 # there is a typo in the raster processing, this is for correction
  
  # time interval for each period
  int_1920 <- interval(ymd("2019-12-01"), ymd("2020-04-13"))
  int_2021 <- interval(ymd("2020-12-01"), ymd("2021-04-13"))

  # extract persistence
  effort[which( date(ymd_hms(effort[, dateField]))  %within% int_1920 ), 
         paste(variable, "_front_persistence", sep = "")] <- raster::extract(x = persist_1920, 
                                                                                     y = effort[which( date(ymd_hms(effort[, dateField]))  %within% int_1920 ),
                                                                                                coordLonLat],
                                                                                     method = "simple")
  effort[which( date(ymd_hms(effort[, dateField]))  %within% int_2021 ), 
         paste(variable, "_front_persistence", sep = "")] <- raster::extract(x = persist_2021, 
                                                                                     y = effort[which( date(ymd_hms(effort[, dateField]))  %within% int_2021 ),
                                                                                                coordLonLat],
                                                                                     method = "simple")
  # extract average distance
  effort[which( date(ymd_hms(effort[, dateField]))  %within% int_1920 ), 
         paste(variable, "_avg_dist_to_front", sep = "")] <- raster::extract(x = avg_dist_1920, 
                                                                                     y = effort[which( date(ymd_hms(effort[, dateField]))  %within% int_1920 ),
                                                                                                coordLonLat],
                                                                                     method = "simple")
  effort[which( date(ymd_hms(effort[, dateField]))  %within% int_2021 ), 
         paste(variable, "_avg_dist_to_front", sep = "")] <- raster::extract(x = avg_dist_2021, 
                                                                                     y = effort[which( date(ymd_hms(effort[, dateField]))  %within% int_2021 ),
                                                                                                coordLonLat],
                                                                                     method = "simple")
  return(effort)
}

for(i in c(2,4,10,18)){
data <- extractStaticVar(effort = data, dateField = "date", coordLonLat = c(9,10),
                      variable = vars[i])
}



### 3 - Join dynamic variables -------------------------------------------------
## create function to extract dynamic variable
extractDynVar <- function(effort, dateField, variable, coordLonLat, t_res) {
  ## extract dates from effort
  dates <- date(ymd_hms(effort[, dateField]))
  dates <- dates[order(dates)]
  
  ## merge raster stacks
  rasters_1920 <- stack(paste(datadir, "/", variable, "-daily-2019-12-01_2020-04-13.envi", sep =""))
  rasters_2021 <- stack(paste(datadir, "/", variable, "-daily-2020-12-01_2021-04-13.envi", sep =""))
  rasters <- stack(rasters_1920, rasters_2021)
  
  ## extract based on dates
  if(t_res == 0){
    
    for (i in 1:length(dates)){
      date <- dates[i]
      closestDate <- ifelse(date <= ymd(str_sub(names(rasters), start = -10))[which.closest(ymd(str_sub(names(rasters), start = -10)), date)], 
                            which.closest(ymd(str_sub(names(rasters), start = -10)), date), 
                            which.closest(ymd(str_sub(names(rasters), start = -10)), date) + 1)
      effort[which(date(ymd_hms(effort[, dateField])) == date), variable] <- raster::extract(rasters[[closestDate]], 
                                                                                             effort[which(date(ymd_hms(effort[, dateField])) == date), coordLonLat], 
                                                                                             method = "simple", df = T)[,2] ## or method = "bilinear"
    }
    
  } else { 
    
    for (i in 1:length(dates)){
      date <- dates[i]
      closestDate <- ifelse(date <= ymd(str_sub(names(rasters), start = -10))[which.closest(ymd(str_sub(names(rasters), start = -10)), date)], 
                            which.closest(ymd(str_sub(names(rasters), start = -10)), date), 
                            which.closest(ymd(str_sub(names(rasters), start = -10)), date) + 1)
      effort[which(date(ymd_hms(effort[, dateField])) == date), paste(variable, "_", t_res, "d", sep = "")] <- raster::extract(rasters[[closestDate -t_res]], 
                                                                                                                               effort[which(date(ymd_hms(effort[, dateField])) == date), coordLonLat], 
                                                                                                                               method = "simple", df = T)[,2] ## or method = "bilinear"
    }
  }
  return(effort)
}


## extract dynamic variables
resolutions <- c(1,2,4,7,10,30)

for (j in 1:length(vars)){
  writeLines(vars[j])
  # variable jour j
  data <- extractDynVar(effort = data, dateField = "date", coordLonLat = c(9,10),
                        variable = vars[j], t_res = 0)
  # les autres résolutions
  for(i in 1:length(resolutions)){
    data <- extractDynVar(effort = data, dateField = "date", coordLonLat = c(9,10),
                          variable = vars[j], t_res = resolutions[i])
  }
}
beepr::beep("fanfare")



#### ---------------------------------- ####
#### ----- Saturate variables --------- ####
#### ---------------------------------- ####

## function
saturate <- function(x, threshold = NULL, lower = TRUE, alpha = 0.01, both = FALSE) {
  if(both) {
    if(is.null(threshold)) {
      threshold <- quantile(x, c(alpha/2, 1-alpha/2), na.rm = TRUE) 
    }
    if(length(threshold) != 2){
      stop("must provide a lower and upper thresholds")
    }
    x <- ifelse(x < min(threshold), min(threshold), x)
    x <- ifelse(x > max(threshold), max(threshold), x)
  }
  else{
    if(is.null(threshold)) {
      threshold <- quantile(x, alpha, na.rm = TRUE) 
    }
    if(lower) { x <- ifelse(x < threshold, threshold, x) }
    else { x <- ifelse(x > threshold, threshold, x) }
  }
  return(x)
}

## generate vector of variable names
var_name <- c("Bathy", 
              "gradChl_front_persistence", "gradChl_avg_dist_to_front",
              "gradNPP_front_persistence", "gradNPP_avg_dist_to_front",
              "gradMLD_front_persistence", "gradMLD_avg_dist_to_front",
              "gradTemp_front_persistence", "gradTemp_avg_dist_to_front")
for(j in 1:length(vars)){
  var_name <- c(var_name, vars[j])
  
  for(t in 1:length(resolutions)){
    var_name <- c(var_name, paste(vars[j], "_", resolutions[t], "d", sep = ""))
  }
}

### 1 - Check distributions of variables and saturate --------------------------
seg <- data
seg <- seg[-which(is.na(seg[,var_name]),arr.ind=T),]

# Chl
boxplot(seg[, str_which(names(seg), pattern = "^Chl(?!_front_distance)") ])
xx <- str_which(names(seg), pattern = "^Chl(?!_front_distance)")
for(j in 1:length(xx)) { seg[, xx[j]] <- saturate(seg[, xx[j]], both = T, alpha = 0.05) }

# NPP
boxplot(seg[, str_which(names(seg), pattern = "^NPP(?!_front_distance)") ])
xx <- str_which(names(seg), pattern = "^NPP(?!_front_distance)")
for(j in 1:length(xx)) { seg[, xx[j]] <- saturate(seg[, xx[j]], both = T, alpha = 0.05) }

# Temp 
boxplot(seg[, str_which(names(seg), pattern = "^Temp(?!_front_distance)") ])
xx <- str_which(names(seg), pattern = "^Temp(?!_front_distance)")
seg[, xx[7]] <- saturate(seg[, xx[7]], both = T, alpha = 0.035) 

# SSH 
boxplot(seg[, str_which(names(seg), pattern = "^SSH(?!_front_distance)") ])
xx <- str_which(names(seg), pattern = "^SSH(?!_front_distance)")
for(j in 1:length(xx)) { seg[, xx[j]] <- saturate(seg[, xx[j]], both = T, alpha = 0.02) }

# Salinity 
boxplot(seg[, str_which(names(seg), pattern = "^Salinity(?!_front_distance)") ])
xx <- str_which(names(seg), pattern = "^Salinity(?!_front_distance)")
for(j in 1:length(xx)) { seg[, xx[j]] <- saturate(seg[, xx[j]], lower = T, alpha = 0.03) }

# EKE 
boxplot(seg[, str_which(names(seg), pattern = "^EKE(?!_front_distance)") ])
xx <- str_which(names(seg), pattern = "^EKE(?!_front_distance)")
for(j in c(1:7)) { seg[, xx[j]] <- saturate(seg[, xx[j]], both = T, alpha = 0.07) }

# CurrentSpeed 
boxplot(seg[, str_which(names(seg), pattern = "^CurrentSpeed(?!_front_distance)") ])
xx <- str_which(names(seg), pattern = "^CurrentSpeed(?!_front_distance)")
for(j in 1:length(xx)) { seg[, xx[j]] <- saturate(seg[, xx[j]], both = T, alpha = 0.02) }
seg[,xx[6]] <- saturate(seg[, xx[6]], both = T, alpha = 0.04) 

# Phyto 
boxplot(seg[, str_which(names(seg), pattern = "^Phyto(?!_front_distance)") ])
xx <- str_which(names(seg), pattern = "^Phyto(?!_front_distance)")
for(j in c(1:5,7)) { seg[, xx[j]] <- saturate(seg[, xx[j]], both = T, alpha = 0.04) }
seg[, xx[6]] <- saturate(seg[, xx[6]], both = T, alpha = 0.06)

# ZEu no need
boxplot(seg[, str_which(names(seg), pattern = "^ZEu(?!_front_distance)") ])

# MLD no need
boxplot(seg[, str_which(names(seg), pattern = "^MLD(?!_front_distance)") ])
xx <- str_which(names(seg), pattern = "^MLD(?!_front_distance)")
for(j in c(3,4,7)) { seg[, xx[j]] <- saturate(seg[, xx[j]], both = T, alpha = 0.04) }

# gradChl 
boxplot(seg[, str_which(names(seg), pattern = "^gradChl(?!_front_distance)") ])
xx <- str_which(names(seg), pattern = "^gradChl(?!_front_distance)")
for(j in 3:length(xx)) { seg[, xx[j]] <- saturate(seg[, xx[j]], both = T, alpha = 0.03) }

# gradNPP 
boxplot(seg[, str_which(names(seg), pattern = "^gradNPP(?!_front_distance)") ])
xx <- str_which(names(seg), pattern = "^gradNPP(?!_front_distance)")
for(j in 3:length(xx)) { seg[, xx[j]] <- saturate(seg[, xx[j]], both = T, alpha = 0.02) }
seg[, xx[2]] <- saturate(seg[, xx[2]], both = T, alpha = 0.025) 

# gradMLD 
boxplot(seg[, str_which(names(seg), pattern = "^gradMLD(?!_front_distance)") ])
xx <- str_which(names(seg), pattern = "^gradMLD(?!_front_distance)")
for(j in c(3,5:9)) { seg[, xx[j]] <- saturate(seg[, xx[j]], both = T, alpha = 0.05) }
seg[, xx[4]] <- saturate(seg[, xx[4]], both = T, alpha = 0.085) 
seg[, xx[2]] <- saturate(seg[, xx[2]], both = T, alpha = 0.04) 

# gradTemp 
boxplot(seg[, str_which(names(seg), pattern = "^gradTemp(?!_front_distance)") ])
xx <- str_which(names(seg), pattern = "^gradTemp(?!_front_distance)")
for(j in 3:length(xx)) { seg[, xx[j]] <- saturate(seg[, xx[j]], both = T, alpha = 0.035) }
seg[, xx[2]] <- saturate(seg[, xx[2]], both = T, alpha = 0.03) 

# gradSSH 
boxplot(seg[, str_which(names(seg), pattern = "^gradSSH(?!_front_distance)") ])
xx <- str_which(names(seg), pattern = "^gradSSH(?!_front_distance)")
for(j in 1:length(xx)) { seg[, xx[j]] <- saturate(seg[, xx[j]], both = T, alpha = 0.03) }
seg[, xx[5]] <- saturate(seg[, xx[5]], both = T, alpha = 0.04) 

# gradSal 
boxplot(seg[, str_which(names(seg), pattern = "^gradSal(?!_front_distance)") ])
xx <- str_which(names(seg), pattern = "^gradSal(?!_front_distance)")
for(j in 1:length(xx)) { seg[, xx[j]] <- saturate(seg[, xx[j]], both = T, alpha = 0.04) }

# DissIC 
boxplot(seg[, str_which(names(seg), pattern = "^DissIC(?!_front_distance)") ])
xx <- str_which(names(seg), pattern = "^DissIC(?!_front_distance)")
for(j in 1:length(xx)) { seg[, xx[j]] <- saturate(seg[, xx[j]], lower = T, alpha = 0.02) }

# SPCO2 
boxplot(seg[, str_which(names(seg), pattern = "^SPCO2(?!_front_distance)") ])
xx <- str_which(names(seg), pattern = "^SPCO2(?!_front_distance)")
for(j in 1:length(xx)) { seg[, xx[j]] <- saturate(seg[, xx[j]], both = T, alpha = 0.03) }

# bathy
boxplot(seg[, "Bathy" ])
seg[, "Bathy"] <- saturate(seg[, "Bathy"], lower = T, alpha = 0.03) 

# gradTemp dist
boxplot(seg[, c("gradChl_front_persistence", "gradNPP_front_persistence", 
                "gradTemp_front_persistence", "gradMLD_front_persistence")])
xx <- str_which(names(seg), pattern = "^gradTemp_front_distance")
for(j in c(1,3:7)) { seg[, xx[j]] <- saturate(seg[, xx[j]], both = T, alpha = 0.04) }
seg[, xx[2]] <- saturate(seg[, xx[2]], both = T, alpha = 0.045)

boxplot(seg[, c("gradChl_avg_dist_to_front", "gradNPP_avg_dist_to_front", 
                "gradTemp_avg_dist_to_front", "gradMLD_avg_dist_to_front")])
seg[, "gradMLD_avg_dist_to_front"] <- saturate(seg[, "gradMLD_avg_dist_to_front"], both = T, alpha = 0.045)



### 3 - Export clean file ------------------------------------------------------
write.table(seg, paste(outdir, "Effort_with_covar.csv", sep = "/"), dec = ".", sep = ";")
save.image(paste(outdir, "Data_processing.RData", sep = "/")) # the environment will be used later






