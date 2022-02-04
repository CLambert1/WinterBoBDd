##--------------------------------------------------------------------------------------------------------
## PROJECT : Infra-seasonal variations of common dolphin distribution in the Bay of Biscay
## SCRIPT : Prediction reliability assessment - proper scoring of predictions
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
library(dsm)
library(mgcv)
library(qpcR)
library(gamm4)
library(itsadug)
library(viridis)
library(pbapply)
library(prettymapr)
library(sf)
library(scoringRules)
library(ggdist)


#### ----------------------- ####
#### ---- SOURCES FILES ---- ####
#### ----------------------- ####

### 1 - Load geographic files --------------------------------------------------
outdir <- "F:/CAPECET/paper/data/3_script_results"
load(paste(outdir, "/Data_processing.RData", sep = ""))
varhandle:: rm.all.but(keep = c("shp_capecet", "shp_capecet_L93", "shp_samm", "shp_samm_L93", "shp_spee", "shp_spee_L93",
                                "surv", "surv_all", "land", "land_L93", "isob", "isob_L93", "bbx", "bbx_L93", "bathy",
                                "lin_cap", "lin_samm", "lin_spee", "lin_spee9", "obs_cap", "obs_samm", "obs_spee", "obs_spee9", "obs_all",
                                "L93", "WGS"))
land_clip <- st_intersection(land_L93, bbx_L93) %>% st_transform(crs = st_crs(WGS)) %>% as_Spatial()
isob_clip <- st_intersection(isob_L93, bbx_L93) %>% st_transform(crs = st_crs(WGS)) %>% as_Spatial()

land <- as_Spatial(land)
isob <- as_Spatial(isob)
surv <- as_Spatial(surv)
surv_all <- as_Spatial(surv_all)


### 2 - Load seg file and fit model --------------------------------------------
outdir <- "F:/CAPECET/paper/data/3_script_results"
sourcedir <- "F:/CAPECET/paper/data/1_data"
funcdir <- "F:/CAPECET/paper/data/0_functions"
setwd(outdir)

seg <- read.table(paste(outdir, "Effort_with_covar.csv", sep = "/"), dec = ".", sep = ";", h = T)
seg <- seg %>%
  mutate(day = as_date(date), 
         month = month(day),
         year = year(date),
         YearWeek = str_c(year(date),week(date),sep = "_"),
         Month = case_when(
           month == 1 ~ "January",
           month == 2 ~ "February",
           month == 3 ~ "March"),
         Period = case_when(year == 2020 ~ "2020", year == 2021 ~ "2021")) %>%
  subset(month %in% c(1:3))

seg$Month <- factor(seg$Month, levels = c("January", "February", "March"))
seg$YearWeek <- factor(seg$YearWeek,
                       levels = c("2020_2", "2020_4", "2020_6", "2020_7", "2020_8", "2020_11",
                                  "2021_2",  "2021_3", "2021_4", "2021_6", "2021_7", "2021_8", "2021_9", "2021_10", "2021_11", "2021_12"))


### 3 - Open model -------------------------------------------------------------
mod1 <- gam(n_ind ~ 1 +s(gradMLD_avg_dist_to_front, k = 4, bs = 'tp') +
              s(MLD_10d, k = 4, bs = 'tp') +
              s(Chl_7d, k = 4, bs = 'tp') +
              s(gradChl_front_distance_2d, k = 4, bs = 'tp')
            , data = seg, offset = log(seg$Effort*2*seg$ESW), family = tw())
covar <- c("gradMLD_avg_dist_to_front", "MLD_10d", "Chl_7d", "gradChl_front_distance_2d")


### 4 - Dates ------------------------------------------------------------------
dates <- c(ymd("2020-01-01"):ymd("2020-04-13"), ymd("2021-01-01"):ymd("2021-04-13"))
dates <- as_date(dates)


### 5 - Open predictions -------------------------------------------------------
# prediction stack - daily model
pred_stack <- stack(paste(outdir, "/DailyModel_Predictions.envi", sep =""))

# prediction - seasonal model
pred_win_20 <- raster(paste(outdir, "/SeasonalModel_2020_prediction.tif", sep = ""))
pred_win_21 <- raster(paste(outdir, "/SeasonalModel_2021_prediction.tif", sep = ""))

# interval and prediction dates
int_20 <- interval(ymd("2019-12-01"), ymd("2020-04-13"))
int_21 <- interval(ymd("2020-12-01"), ymd("2021-04-13"))

date20 <- which(ymd(str_sub(names(pred_stack), start = -10)) %within% int_20)
date21 <- which(ymd(str_sub(names(pred_stack), start = -10)) %within% int_21)


### 6 - Generate mask to clip on shelf -----------------------------------------
masking_layer <- raster(paste(sourcedir, "/Variables/Bathy.tif", sep = ""))
masking_layer <- crop(masking_layer, extent(-4.8,-1.2,44.3,47.1))
masking_layer <- resample(masking_layer, raster(paste(sourcedir, "/Variables/gradTemp_frontal_persistence-2019-12-01_2020-04-13.envi", sep ="")))

masking_layer[masking_layer>0] <- NA
masking_layer[masking_layer<(-2500)] <- NA

masking_layer <- mask(masking_layer, land, inverse = T)

crs(pred_win_20) <- crs(masking_layer)
crs(pred_win_21) <- crs(masking_layer)


### 7 - load functions ---------------------------------------------------------
source(paste(funcdir, "fct_simulations.R", sep = "/"))




#### -------------------------- ####
#### ---- Prepare grid and ---- ####
#### ------- densities -------- ####
#### -------------------------- ####


### 1 - Project reference raster extent for coordinates ------------------------
projExt <- projectExtent(masking_layer, crs = crs(shp_capecet_L93))
# visualize to check
pp <- projectRaster(from = masking_layer, to = projExt) 
plot(pp)


### 2 - Compute prediction grid ------------------------------------------------
predgrid <- raster(nrows = 102, ncols = 192, 
                   xmn = 78695.73, xmx = 510417.1, 
                   ymn = 6358856, ymx = 6695834,
                   crs = crs(shp_capecet_L93), resolution = c(2248.549, 3303.698),
                   vals = 1:19584)
plot(predgrid)




#### --------------------------------- ####
#### ---- Simulations predictions ---- ####
#### ------------ modele ------------- ####
#### --------------------------------- ####
# caution: surveys are treated separatly to avoid conflict on SegID. 
# the seg ID out of the simulation functions is structural, not based upon the original segID

## 1 - Compute CAPECET simulations ---------------------------------------------
dates_ech_CAP <- unique(date(lin_cap$date))
lin_cap <- lin_cap %>% rename(seg_id = segID,
                              seg_length = segLengKm,
                              sea_state = seaState)

simulation_CAP <- pbapply::pblapply(seq_along(dates_ech_CAP), function(i){
  # select daily prediction layer
  layer <- pred_stack[[which(as_date(str_sub(names(pred_stack), 5)) == dates_ech_CAP[i])]]
  crs(layer) <- crs(masking_layer)
  
  # select daily sampling design (projected in Lamb93)
  transects <- lin_cap %>% subset(date(date) == dates_ech_CAP[i]) %>% st_transform(crs = L93) 
  
  # project raster layer (in Lamb93)
  layer_df <- as.data.frame(projectRaster(from = layer, to = predgrid), xy = T) # project raster and convert to table
  layer_df <- layer_df %>% drop_na()
  layer_L93 <- rasterize(x = layer_df[,1:2], field = layer_df[,3], y = predgrid) # re assign to regular grid
  
  # simulate
  sim <- simul(predgrid = predgrid,
               truth = layer_L93,
               is.raster = TRUE, 
               N_obs = 10000,
               mean_group_size = 6.6, # extracted from CDS
               projection_espg = L93,
               design_transects = transects,
               strip = FALSE,
               esw = 0.172, # extracted from CDS
               soap = TRUE,
               seed_id = NULL,
               distance.max = 700,
               n_sim = 100)
  return(sim)
})


## 2 - Compute SPEE simulations ------------------------------------------------
dates_ech_SPEE <- unique(lin_spee$day)
lin_spee <- lin_spee %>% rename(seg_id = segID,
                              seg_length = segLengKm,
                              sea_state = seaState)

simulation_SPEE <- pbapply::pblapply(seq_along(dates_ech_SPEE), function(i){
  # select daily prediction layer
  layer <- pred_stack[[which(as_date(str_sub(names(pred_stack), 5)) == dates_ech_SPEE[i])]]
  crs(layer) <- crs(masking_layer)
  
  # select daily sampling design (projected in Lamb93)
  transects <- lin_spee %>% subset(date(date) == dates_ech_SPEE[i]) %>% st_transform(crs = L93) 
  
  # project raster layer (in Lamb93)
  layer_df <- as.data.frame(projectRaster(from = layer, to = predgrid), xy = T) # project raster and convert to table
  layer_df <- layer_df %>% drop_na()
  layer_L93 <- rasterize(x = layer_df[,1:2], field = layer_df[,3], y = predgrid) # re assign to regular grid
  
  # simulate
  sim <- simul(predgrid = predgrid,
               truth = layer_L93,
               is.raster = TRUE, 
               N_obs = 10000,
               mean_group_size = 5.40,
               projection_espg = L93,
               design_transects = transects,
               strip = FALSE,
               esw = 0.156,
               soap = TRUE,
               seed_id = NULL,
               distance.max = 700,
               n_sim = 100)
  return(sim)
})


## 3 - Compute SAMM simulations ------------------------------------------------
dates_ech_SAMM <- unique(lin_samm$day)
lin_samm <- lin_samm %>% rename(seg_id = segID,
                                seg_length = segLengKm,
                                sea_state = seaState)

simulation_SAMM <- pbapply::pblapply(seq_along(dates_ech_SAMM), function(i){
  # select daily prediction layer
  layer <- pred_stack[[which(as_date(str_sub(names(pred_stack), 5)) == dates_ech_SAMM[i])]]
  crs(layer) <- crs(masking_layer)
  
  # select daily sampling design (projected in Lamb93)
  transects <- lin_samm %>% subset(date(day) == dates_ech_SAMM[i]) %>% st_transform(crs = L93) 
  
  # project raster layer (in Lamb93)
  layer_df <- as.data.frame(projectRaster(from = layer, to = predgrid), xy = T) # project raster and convert to table
  layer_df <- layer_df %>% drop_na()
  layer_L93 <- rasterize(x = layer_df[,1:2], field = layer_df[,3], y = predgrid) # re assign to regular grid
  
  # simulate
  sim <- simul(predgrid = predgrid,
               truth = layer_L93,
               is.raster = TRUE, 
               N_obs = 10000,
               mean_group_size = 8.83,
               projection_espg = L93,
               design_transects = transects,
               strip = FALSE,
               esw = 0.198,
               soap = TRUE,
               seed_id = NULL,
               distance.max = 700,
               n_sim = 100)
  return(sim)
})



## 4 - Compute SPEE9 simulations -----------------------------------------------
dates_ech_SPEE9 <- unique(lin_spee9$day)
lin_spee9 <- lin_spee9 %>% rename(seg_id = segID,
                                seg_length = segLengKm,
                                sea_state = seaState)

simulation_SPEE9 <- pbapply::pblapply(seq_along(dates_ech_SPEE9), function(i){
  # select daily prediction layer
  layer <- pred_stack[[which(as_date(str_sub(names(pred_stack), 5)) == dates_ech_SPEE9[i])]]
  crs(layer) <- crs(masking_layer)
  
  # select daily sampling design (projected in Lamb93)
  transects <- lin_spee9 %>% subset(date(day) == dates_ech_SPEE9[i]) %>% st_transform(crs = L93) 
  
  # project raster layer (in Lamb93)
  layer_df <- as.data.frame(projectRaster(from = layer, to = predgrid), xy = T) # project raster and convert to table
  layer_df <- layer_df %>% drop_na()
  layer_L93 <- rasterize(x = layer_df[,1:2], field = layer_df[,3], y = predgrid) # re assign to regular grid
  
  # simulate
  sim <- simul(predgrid = predgrid,
               truth = layer_L93,
               is.raster = TRUE, 
               N_obs = 10000,
               mean_group_size = 7.81,
               projection_espg = L93,
               design_transects = transects,
               strip = FALSE,
               esw = 0.181,
               soap = TRUE,
               seed_id = NULL,
               distance.max = 700,
               n_sim = 100)
  return(sim)
})




#### --------------------------------- ####
#### ---- Simulations predictions ---- ####
#### ----- Uniform distribution ------ ####
#### --------------------------------- ####

## 1 - Generate uniform layer --------------------------------------------------
unif <- raster(nrows = 102, ncols = 192,
               xmn = 78695.73, xmx = 510417.1,
               ymn = 6358856, ymx = 6695834,
               crs = crs(shp_capecet_L93), resolution = c(2248.549, 3303.698),
               vals = 1)


## 2 - CAPECET -----------------------------------------------------------------
simulation_CAP_unif <- simul(predgrid = predgrid,
                             truth = unif,
                             is.raster = TRUE, 
                             N_obs = 10000,
                             mean_group_size = 6.6,
                             projection_espg = L93,
                             design_transects = lin_cap %>% st_transform(crs = L93) ,
                             strip = FALSE,
                             esw = 0.172,
                             soap = TRUE,
                             seed_id = NULL,
                             distance.max = 700,
                             n_sim = 100)


## 2 - SPEE --------------------------------------------------------------------
simulation_SPEE_unif <- simul(predgrid = predgrid,
                              truth = unif,
                              is.raster = TRUE, 
                              N_obs = 10000,
                              mean_group_size = 5.40,
                              projection_espg = L93,
                              design_transects = lin_spee %>% st_transform(crs = L93) ,
                              strip = FALSE,
                              esw = 0.156,
                              soap = TRUE,
                              seed_id = NULL,
                              distance.max = 700,
                              n_sim = 100)


## 3 - SAMM --------------------------------------------------------------------
simulation_SAMM_unif <- simul(predgrid = predgrid,
                              truth = unif,
                              is.raster = TRUE, 
                              N_obs = 10000,
                              mean_group_size = 8.83,
                              projection_espg = L93,
                              design_transects = lin_samm %>% st_transform(crs = L93) ,
                              strip = FALSE,
                              esw = 0.198,
                              soap = TRUE,
                              seed_id = NULL,
                              distance.max = 700,
                              n_sim = 100)


## 4 - SPEE9 -------------------------------------------------------------------
simulation_SPEE9_unif <- simul(predgrid = predgrid,
                               truth = unif,
                               is.raster = TRUE, 
                               N_obs = 10000,
                               mean_group_size = 7.81,
                               projection_espg = L93,
                               design_transects = lin_spee9 %>% st_transform(crs = L93) ,
                               strip = FALSE,
                               esw = 0.181,
                               soap = TRUE,
                               seed_id = NULL,
                               distance.max = 700,
                               n_sim = 100)



#### --------------------------------- ####
#### ---- Simulations predictions ---- ####
#### -------- Seasonal model --------- ####
#### --------------------------------- ####

## 0 - Project raster layers ---------------------------------------------------
# 2020
pred_win_20_df <- as.data.frame(projectRaster(from = pred_win_20, to = predgrid), xy = T) # project raster and convert to table
pred_win_20_df <- pred_win_20_df %>% drop_na()
pred_win_20_L93 <- rasterize(x = pred_win_20_df[,1:2], field = pred_win_20_df[,3], y = predgrid) # re assign to regular grid

# 2021
pred_win_21_df <- as.data.frame(projectRaster(from = pred_win_21, to = predgrid), xy = T) 
pred_win_21_df <- pred_win_21_df %>% drop_na()
pred_win_21_L93 <- rasterize(x = pred_win_21_df[,1:2], field = pred_win_21_df[,3], y = predgrid) 



## 1 - Compute CAPECET simulations ---------------------------------------------
simulation_CAP_win <- simul(predgrid = predgrid,
                             truth = pred_win_20_L93,
                             is.raster = TRUE, 
                             N_obs = 10000,
                             mean_group_size = 6.6,
                             projection_espg = L93,
                             design_transects = lin_cap %>% st_transform(crs = L93) ,
                             strip = FALSE,
                             esw = 0.172,
                             soap = TRUE,
                             seed_id = NULL,
                             distance.max = 700,
                             n_sim = 100)



## 2 - Compute SPEE simulations ------------------------------------------------
simulation_SPEE_win <- simul(predgrid = predgrid,
                              truth = pred_win_20_L93,
                              is.raster = TRUE, 
                              N_obs = 10000,
                              mean_group_size = 5.40,
                              projection_espg = L93,
                              design_transects = lin_spee %>% st_transform(crs = L93) ,
                              strip = FALSE,
                              esw = 0.156,
                              soap = TRUE,
                              seed_id = NULL,
                              distance.max = 700,
                              n_sim = 100)


## 3 - Compute SAMM simulations ------------------------------------------------
simulation_SAMM_win <- simul(predgrid = predgrid,
                             truth = pred_win_21_L93,
                             is.raster = TRUE, 
                             N_obs = 10000,
                             mean_group_size = 8.83,
                             projection_espg = L93,
                             design_transects = lin_samm %>% st_transform(crs = L93) ,
                             strip = FALSE,
                             esw = 0.198,
                             soap = TRUE,
                             seed_id = NULL,
                             distance.max = 700,
                             n_sim = 100)


## 4 - Compute SPEE9 simulations -----------------------------------------------
simulation_SPEE9_win <- simul(predgrid = predgrid,
                               truth = pred_win_21_L93,
                               is.raster = TRUE, 
                               N_obs = 10000,
                               mean_group_size = 7.81,
                               projection_espg = L93,
                               design_transects = lin_spee9 %>% st_transform(crs = L93) ,
                               strip = FALSE,
                               esw = 0.181,
                               soap = TRUE,
                               seed_id = NULL,
                               distance.max = 700,
                               n_sim = 100)





#### ------------------------------- ####
#### ------ Retrieve dataset ------- ####
#### ----- Join observed value ----- ####
#### ------------------------------- ####

### 1 - Daily prediction -------------------------------------------------------
sim_detect_CAP <- simulation_CAP %>% 
  purrr::map(`[[`, 1) %>% 
  data.table::rbindlist() 
sim_detect_SPEE <- simulation_SPEE %>%
  purrr::map(`[[`, 1) %>%
  data.table::rbindlist()
sim_detect_SAMM <- simulation_SAMM %>% 
  purrr::map(`[[`, 1) %>% 
  data.table::rbindlist()
sim_detect_SPEE9 <- simulation_SPEE9 %>%
  purrr::map(`[[`, 1) %>%
  data.table::rbindlist()
Daily_simulations <- rbind(sim_detect_CAP, sim_detect_SPEE, sim_detect_SAMM, sim_detect_SPEE9)


### 2 - Uniform distribution ---------------------------------------------------
sim_detect_CAP_unif <- simulation_CAP_unif$Detections
sim_detect_SPEE_unif <- simulation_SPEE_unif$Detections
sim_detect_SAMM_unif <- simulation_SAMM_unif$Detections
sim_detect_SPEE9_unif <- simulation_SPEE9_unif$Detections
Unif_simulations <- rbind(sim_detect_CAP_unif, sim_detect_SPEE_unif, sim_detect_SAMM_unif, sim_detect_SPEE9_unif)


### 3 - Seasonal prediction ----------------------------------------------------
sim_detect_CAP_win <- simulation_CAP_win$Detections
sim_detect_SPEE_win <- simulation_SPEE_win$Detections
sim_detect_SAMM_win <- simulation_SAMM_win$Detections
sim_detect_SPEE9_win <- simulation_SPEE9_win$Detections
Seasonal_simulations <- rbind(sim_detect_CAP_win, sim_detect_SPEE_win, sim_detect_SAMM_win, sim_detect_SPEE9_win)


### 4 - Join observed values ---------------------------------------------------
Daily_simulations <- left_join(Daily_simulations, 
                              dplyr::select(seg, c("Sample.Label", "n_ind")), 
                              by = "Sample.Label")
names(Daily_simulations)[111] <- c("Obs_n_idv")

Unif_simulations <- left_join(Unif_simulations, 
                              dplyr::select(seg, c("Sample.Label", "n_ind")), 
                              by = "Sample.Label")
names(Unif_simulations)[111] <- c("Obs_n_idv")

Seasonal_simulations <- left_join(Seasonal_simulations, 
                                dplyr::select(seg, c("Sample.Label", "n_ind")), 
                                by = "Sample.Label")
names(Seasonal_simulations)[111] <- c("Obs_n_idv")



#### --------------------------- ####
#### ----- Proper Scoring ------ ####
#### --------- Rules ----------- ####
#### --------------------------- ####


### 1 - Prepare tables ---------------------------------------------------------
Daily_simulations <- Daily_simulations %>% drop_na(Obs_n_idv)

Unif_simulations <- Unif_simulations %>% drop_na(Obs_n_idv)

Seasonal_simulations <- Seasonal_simulations %>% drop_na(Obs_n_idv)


### 2 - Compute scores ---------------------------------------------------------
Daily_simulations$Scores <- crps_sample(y = Daily_simulations$Obs_n_idv,
                                        dat = as.matrix(Daily_simulations[, 10:109]))

Unif_simulations$Scores <- crps_sample(y = Unif_simulations$Obs_n_idv,
                                      dat = as.matrix(Unif_simulations[, 10:109]))

Seasonal_simulations$Scores <- crps_sample(y = Seasonal_simulations$Obs_n_idv,
                                        dat = as.matrix(Seasonal_simulations[, 10:109]))


### 3 - Pivot longer for plotting ----------------------------------------------
Simulations <- Daily_simulations[, c(1:09, 110:111)]
Simulations$Daily_predictions <- Daily_simulations[, "Scores"]
Simulations$Uniform_distribution <- Unif_simulations[, "Scores"]
Simulations$Seasonal_prediction <- Seasonal_simulations[, "Scores"]

Simulations_long <- Simulations %>% pivot_longer(cols = c("Daily_predictions", "Uniform_distribution", "Seasonal_prediction")) %>%
  mutate(Model = case_when(
    name == "Daily_predictions" ~"Daily predictions",
    name == "Uniform_distribution" ~"Uniform distribution",
    name == "Seasonal_prediction" ~"Seasonal prediction"
  ),
  Year = year(day))
write.table(Simulations_long, paste(outdir, "ProperScores.txt", sep = "/"), row.names = F)


### 4 - Plots ------------------------------------------------------------------
ggplot(Simulations_long, aes(y = value, x = Model)) + stat_halfeye() + ylim(c(0,5)) +
  theme_bw() + coord_flip() + xlab("") + ylab("Scores")  + facet_grid(Year~.) 
ggsave(paste(outdir, "ProperScores.png", sep = "/"), dev = "png", dpi = 350, width = 5, height = 6)


