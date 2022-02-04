##--------------------------------------------------------------------------------------------------------
## PROJECT : Infra-seasonal variations of common dolphin distribution in the Bay of Biscay
## SCRIPT : Compute observed encounter rates and Conventional Distance Sampling analysis
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
library(mvtnorm)


#### ----------------------- ####
#### ---- SOURCES FILES ---- ####
#### ----------------------- ####
outdir <- "F:/CAPECET/paper/data/3_script_results"
sourcedir <- "F:/CAPECET/paper/data/1_data"
setwd(outdir)

L93 <- 2154
WGS <- 4326


### 1 - Load land, isobath -----------------------------------------------------
## A - Open, clip and project land
land <- st_read(dsn = paste(sourcedir, "Map_files", sep = "/"), layer = "DCSMM_EI_terre_carto_wgs84")
bbx <- as(raster::extent(-4.8,0.5,44.3,47.1), "SpatialPolygons")
raster::crs(bbx) <- raster::crs(land)
bbx_L93 <- bbx %>% st_as_sf() %>% st_transform(crs = st_crs(L93)) 
land_L93 <- land %>% st_transform(crs = st_crs(L93))


## B - Open, clip and project isob
isob <- st_read(dsn = paste(sourcedir, "Map_files", sep = "/"), layer = "isobath_reduits")
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



#### ---------------------- ####
#### ---- Prepare data ---- ####
#### ---------------------- ####
# Association of region with its surface
Block_Area <- rbind(shp_capecet@data, shp_spee@data, shp_samm@data)

# NB: this section uses dedicated functions from the pelaCDS package

### 1 - Prepare CAPECET data ---------------------------------------------------
# prepare effort
eff_cap <- lin_cap %>%
  mutate(center_time = datehhmmss,
         session = paste(session, day, sep = "_")) %>% # inclue le jour dans la session pour avoir une estimation par jour par CDS
  clean_effort_cds(., multi_survey = T)
eff_cap2 <- prepare_effort_cds(effort_base = eff_cap,
                                optimal = T,
                                block_area = Block_Area[1,],
                                unit_km = FALSE)
eff_cap2$Area <- Block_Area[1,1]


# prepare obs
obs_cap$taxon <- "Marine Mammal"
obs_cap$group <- "Marine Mammal"
obs_cap$family <- "Cetacea"
sp <- c("DELDEL", "STEDEL")
obs_cap2 <- obs_cap %>%
  mutate(center_time = datehhmmss,
         session = paste(session, day, sep = "_")) %>%
  clean_obs_cds(., multi_survey = T)
obs_cap2 <- prepare_obs_cds(sp = sp,
                             obs_base = obs_cap2,
                             truncation = 0.5,
                             legdata = eff_cap2,
                             projection = 2154,
                             unit_km = FALSE)
obs_cap2$distdata$Area <- Block_Area[1,1]

# Add observation on effort
cap <- add_obs_cds(countdata_leg = obs_cap2$countdata_leg,
                    legdata = eff_cap2)



### 2 - Prepare SPEE data ------------------------------------------------------
# prepare effort
eff_spee <- lin_spee %>%
  mutate(center_time = datehhmmss,
         session = paste(session, day, sep = "_")) %>%
  clean_effort_cds(., multi_survey = T)
eff_speeb <- prepare_effort_cds(effort_base = eff_spee,
                                 optimal = T,
                                 block_area = Block_Area[2,],
                                 unit_km = FALSE)
eff_speeb$Area <- Block_Area[2,1]


# prepare obs
obs_spee$taxon <- "Marine Mammal"
obs_spee$group <- "Marine Mammal"
obs_spee$family <- "Cetacea"
sp <- c("DELDEL", "STEDEL")
obs_speeb <- obs_spee %>%
  mutate(center_time = datehhmmss,
         session = paste(session, day, sep = "_")) %>%
  clean_obs_cds(., multi_survey = T)
obs_speeb <- prepare_obs_cds(sp = sp,
                              obs_base = obs_speeb,
                              truncation = 0.5,
                              legdata = eff_speeb,
                              projection = 2154,
                              unit_km = FALSE)

# Add observation on effort
spee <- add_obs_cds(countdata_leg = obs_speeb$countdata_leg,
                     legdata = eff_speeb)


### 3 - Prepare SAMM data ------------------------------------------------------
# prepare effort
eff_samm <- lin_samm %>%
  mutate(center_time = datehhmmss,
         session = paste(session, day, sep = "_")) %>%
  clean_effort_cds(., multi_survey = T)
eff_samm2 <- prepare_effort_cds(effort_base = eff_samm,
                                optimal = T,
                                block_area = Block_Area[3,],
                                unit_km = FALSE)
eff_samm2$Area <- Block_Area[3,1]


# prepare obs
obs_samm$taxon <- "Marine Mammal"
obs_samm$group <- "Marine Mammal"
obs_samm$family <- "Cetacea"
sp <- c("DELDEL", "STECOE", "STEDEL")
obs_samm2 <- obs_samm %>%
  mutate(center_time = datehhmmss,
         session = paste(session, day, sep = "_")) %>%
  clean_obs_cds(., multi_survey = T)
obs_samm2 <- prepare_obs_cds(sp = sp,
                             obs_base = obs_samm2,
                             truncation = 0.5,
                             legdata = eff_samm2,
                             projection = 2154,
                             unit_km = FALSE)


# Add observation on effort
samm <- add_obs_cds(countdata_leg = obs_samm2$countdata_leg,
                    legdata = eff_samm2)



### 4 - Prepare SPEE9 data -----------------------------------------------------
# prepare effort
eff_spee9 <- lin_spee9 %>%
  mutate(center_time = datehhmmss,
         session = paste(session, day, sep = "_")) %>%
  clean_effort_cds(., multi_survey = T)
eff_spee9b <- prepare_effort_cds(effort_base = eff_spee9,
                                optimal = T,
                                block_area = Block_Area[2,],
                                unit_km = FALSE)
eff_spee9b$Area <- Block_Area[2,1]

# prepare obs
obs_spee9$taxon <- "Marine Mammal"
obs_spee9$group <- "Marine Mammal"
obs_spee9$family <- "Cetacea"
sp <- c("DELDEL", "STEDEL")
obs_spee9b <- obs_spee9 %>%
  mutate(center_time = datehhmmss,
         session = paste(session, day, sep = "_")) %>%
  clean_obs_cds(., multi_survey = T)
obs_spee9b <- prepare_obs_cds(sp = sp,
                              obs_base = obs_spee9b,
                              truncation = 0.5,
                              legdata = eff_spee9b,
                              projection = 2154,
                              unit_km = FALSE)

# Add observation on effort
spee9 <- add_obs_cds(countdata_leg = obs_spee9b$countdata_leg,
                    legdata = eff_spee9b)


### 5 - Merge effort tables ----------------------------------------------------
effort <- rbind(cap, spee, samm, spee9)
distdata <- rbind(obs_cap2$distdata, obs_speeb$distdata, obs_samm2$distdata, obs_spee9b$distdata)
# obsdata <- rbind(obs_cap2$obsdata, obs_speeb$obsdata, obs_samm2$obsdata, obs_spee9b$obsdata)





#### ---------------------- ####
#### ---- SUMMARY & ER ---- #### 
#### ---------------------- ####

## function
my_summary <- function (df, var, what = c("YearWeek", "session", "NONE")) {
  df$var <- df[, var]
  if(what == "YearWeek") {
    tt <- df %>% 
      group_by(YearWeek) %>% 
      summarise(
        n = length(var),
        sum = sum(var, na.rm = TRUE),
        moyenne = mean(var, na.rm = TRUE)
      )
  }
  if(what == "session") {
    tt <- df %>% 
      group_by(session) %>% 
      summarise(
        n = length(var),
        sum = sum(var, na.rm = TRUE),
        moyenne = mean(var, na.rm = TRUE)
      )
  }
  if(what == "NONE") {
    tt <- df  %>% 
      summarise(
        n = length(var),
        sum = sum(var, na.rm = TRUE),
        moyenne = mean(var, na.rm = TRUE)
      )
  }
  return(as.data.frame(tt))
}


## summary effort total (in all conditions)
my_summary(df = as_Spatial(lin_cap)@data, var = "segLengKm", what = "session")
my_summary(df = as_Spatial(lin_spee)@data, var = "segLengKm", what = "session")
my_summary(df = as_Spatial(lin_samm)@data, var = "segLengKm", what = "session")
my_summary(df = as_Spatial(lin_spee9)@data, var = "segLengKm", what = "session")


## Summary obs
summary_sightings <- my_summary(df = effort, var = "n_detection", what = "session")
summary_individuals <- my_summary(df = effort, var = "n_ind", what = "session")
effort$ER <- effort$n_ind/effort$Effort
summary_ER <- my_summary(df = effort, var = "ER", what = "session")

## summary effort
summary_effort <- my_summary(df = effort, var = "Effort", what = "session")
summary_effort <- summary_effort[,1:3]
colnames(summary_effort) <- c("session", "N seg", "Total effort (km)")

## paste summary obs on effort and compute ER
summary_effort$`N individuals` <- summary_individuals$sum
summary_effort$`N detections` <- summary_sightings$sum
summary_effort$ER <- summary_ER$moyenne

## format table
summary_effort$Survey <- sapply(str_split(summary_effort$session, pattern = "_"), "[[", 1)
summary_effort$Session <- sapply(str_split(summary_effort$session, pattern = "_"), "[[", 2)
summary_effort$Day <- date(sapply(str_split(summary_effort$session, pattern = "_"), "[[", 3))
summary_effort$Year <- year(summary_effort$Day)
summary_effort$Period <- ifelse(summary_effort$Year == 2021, "2021", "2020")


## summarize again
summary_effort %>% 
  group_by(Survey, Session) %>% 
  summarise(
    n_leg = length(`Total effort (km)`),
    effort = sum(`Total effort (km)`, na.rm = TRUE) ,
    n_ind = sum(`N individuals`, na.rm = TRUE),
    n_det = sum(`N detections`, na.rm = TRUE),
    ER = mean(ER, na.rm = T))

# values used in the paper
summary_effort %>% 
  group_by(Survey, Session) %>% 
  summarise(
    n_leg = length(`Total effort (km)`),
    effort = sum(`Total effort (km)`, na.rm = TRUE) ,
    n_ind = sum(`N individuals`, na.rm = TRUE),
    n_det = sum(`N detections`, na.rm = TRUE))

effort %>% 
  group_by(survey) %>% 
  summarise(
    ER = mean(ER, na.rm = T))

effort %>% subset(survey == "SPEE") %>%
  mutate(year = year(as_date(str_sub(session, 8)))) %>%
  group_by(year) %>% 
  summarise(
    ER = mean(ER, na.rm = T))



#### ----------------------------------- ####
#### ---- Adjust detection function ---- ####
#### ----------------------------------- ####

### 1 - CAPECET ----------------------------------------------------------------
obs_cap2$distdata <- obs_cap2$distdata %>% mutate(Region.Label = session) # use session as region label to get an estimate per day
hn_cap <- cds_detection4cds(distdata = obs_cap2$distdata,
                            bin = seq(0.0, 0.4, 0.05),
                            key = "halfnorm",
                            upper = 0.4,
                            is_seabird = FALSE)
hn_cap$graph
summary(hn_cap$distFit)


### 2 - SPEE ----------------------------------------------------------------
obs_speeb$distdata <- obs_speeb$distdata %>% mutate(Region.Label = session)
hn_spee <- cds_detection4cds(distdata = obs_speeb$distdata,
                            bin = seq(0.0, 0.4, 0.05),
                            key = "halfnorm",
                            upper = 0.4,
                            is_seabird = FALSE) 
hn_spee$graph
summary(hn_spee$distFit)


### 3 - SAMM ----------------------------------------------------------------
obs_samm2$distdata <- obs_samm2$distdata %>% mutate(Region.Label = session) 
hn_samm <- cds_detection4cds(distdata = obs_samm2$distdata,
                             bin = seq(0.0, 0.4, 0.05),
                             key = "halfnorm",
                             upper = 0.4,
                             is_seabird = FALSE)
hn_samm$graph
summary(hn_samm$distFit)


### 4 - SPEE9 ----------------------------------------------------------------
obs_spee9b$distdata <- obs_spee9b$distdata %>% mutate(Region.Label = session)
hn_spee9 <- cds_detection4cds(distdata = obs_spee9b$distdata,
                             bin = seq(0.0, 0.4, 0.05),
                             key = "halfnorm",
                             upper = 0.4,
                             is_seabird = FALSE)
hn_spee9$graph
summary(hn_spee9$distFit)


### 5 - Format dataset ---------------------------------------------------------
## A - Extract values from CDS summaries
# estimate encounter rates with mean group size per day
CDS_table <- rbind(hn_cap$distFit$dht$individuals$summary,
                   hn_spee$distFit$dht$individuals$summary,
                   hn_samm$distFit$dht$individuals$summary,
                   hn_spee9$distFit$dht$individuals$summary)
colnames(CDS_table)[c(1,6)] <- c("session", "ER.hat")

# estimated density per day
D_table <- rbind(hn_cap$distFit$dht$individuals$D,
                 hn_spee$distFit$dht$individuals$D,
                 hn_samm$distFit$dht$individuals$D,
                 hn_spee9$distFit$dht$individuals$D)
colnames(D_table)[c(1:6)] <- c("session", "D.hat", "D.se", "D.cv", "D.lcl", "D.ucl")

# estimated density per leg (density = Dhat)
leg_table <- rbind(hn_cap$distFit$dht$individuals$bysample,
                   hn_spee$distFit$dht$individuals$bysample,
                   hn_samm$distFit$dht$individuals$bysample,
                   hn_spee9$distFit$dht$individuals$bysample)
colnames(leg_table)[1:2] <- c("sesion", "Sample.Label")


## B - Paste to summary_effort (estimates of ER and D)
summary_effort <- left_join(summary_effort,
                            dplyr::select(CDS_table, -Area, -CoveredArea),
                            by = "session")
summary_effort <- left_join(summary_effort,
                            dplyr::select(D_table, -df),
                            by = "session")
write.table(summary_effort, paste(outdir, "Encounter_rate_per_day_2Years.txt", sep = "/"))


## C - Paste to legs (estimates of D)
effort <- left_join(effort,
                    dplyr::select(leg_table, -Area),
                    by = "Sample.Label")


