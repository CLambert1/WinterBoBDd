##--------------------------------------------------------------------------------------------------------
## PROJECT : Infra-seasonal variations of common dolphin distribution in the Bay of Biscay
## SCRIPT : Extract seasonal vars, format data and fit the seasonal model
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
library(viridis)
library(mgcv)
library(qpcR)
library(gamm4)
library(itsadug)
library(ggridges)
library(viridis)
library(prettymapr)

lapply(c("dplyr", "tidyr","lubridate", "wdpar", "stringr", "curl", "units", "rgdal", "fields", "birk", "viridis", "rasterVis", 
         "downloader", "directlabels", "grec", "ggplot2", "rgl", "ncdf4", "mapdata", "maptools", "lattice", "reshape2", "raster", "reticulate"), 
       library, character.only = TRUE)


#### ---------------------------- ####
#### ----- OPEN DATA SOURCE ----- ####
#### ---------------------------- ####

### 1 - Open variable mean rasters ---------------------------------------------
load("F:/CAPECET/Mean_covariates.RData")
var_list <- c("Temp", "Salinity", "CurrentSpeed", "Chl", "NPP", "Phyto",
              "EKE", "gradSSH", "gradTemp", "gradMLD", "gradSal",
              "SSH", "gradChl", "gradNPP", "ZEu", "MLD", "DissIC","SPCO2")

for (i in 1:length(var_list)){
  map1 <- get(paste("map_", var_list[i], "_20", sep = ""))
  map2 <- get(paste("map_", var_list[i], "_21", sep = ""))
  map_mn <- stack(map1[[1]], map2[[1]])
  map_sd <- stack(map1[[2]], map2[[2]])
  names(map_mn) <- c("2020", "2021")
  names(map_sd) <- c("2020", "2021")
  assign(paste(var_list[i], "Mean", sep = "_"), map_mn)
  assign(paste(var_list[i], "SD", sep = "_"), map_sd)
}
rm(list = c(ls(pattern = "map"), ls(pattern = "ts"), ls(pattern = "persist"), 
            "shp", "shp0", "ras", "rasters", 
            "surv", "ll", "lin_cap", "lin_spee", "isob", "bbx"))


outdir <- "F:/CAPECET/paper/data/3_script_results"
sourcedir <- "F:/CAPECET/paper/data/1_data"
funcdir <- "F:/CAPECET/paper/data/0_functions"
setwd(outdir)


### 2 - Open effort file -------------------------------------------------------
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

data <- seg[, c(1:27, 182:186)] # drop dynamic variables


### 3 - open land --------------------------------------------------------------------
land <- readOGR(dsn = paste(sourcedir, "Map_files", sep = "/"), layer = "DCSMM_EI_terre_carto_wgs84")
bbx <- as(extent(-4.8,0.5,44.3,47.1), "SpatialPolygons")
crs(bbx) <- crs(land)
land <- gIntersection(land, bbx, byid = TRUE)



#### -------------------------- ####
#### ----- EXTRACT COVARS ----- ####
#### ------- & PRED GRID ------ ####
#### -------------------------- ####

### 1 - Extraction -------------------------------------------------------------
extractVar <- function(effort, variable, coordLonLat, what = "Mean") {
  if(what == "Mean"){
    ras <- get(paste(variable, "Mean", sep = "_"))
    ras_20 <- ras[["X2020"]]
    ras_21 <- ras[["X2021"]]
    
    effort[which(effort$Period == "2020"), 
           paste(variable, "mean", sep = "_")] <- raster::extract(ras_20, 
                                                                  effort[which(effort$Period == "2020"), 
                                                                         coordLonLat], 
                                                                  method = "simple", df = T)[,2] 
    effort[which(effort$Period == "2020-2021"), 
           paste(variable, "mean", sep = "_")] <- raster::extract(ras_21, 
                                                                  effort[which(effort$Period == "2021"), 
                                                                         coordLonLat], 
                                                                  method = "simple", df = T)[,2] 
    return(effort)
  }
  if(what == "SD"){
    ras <- get(paste(variable, "SD", sep = "_"))
    ras_20 <- ras[["X2020"]]
    ras_21 <- ras[["X2021"]]
    
    effort[which(effort$Period == "2020"), 
           paste(variable, "sd", sep = "_")] <- raster::extract(ras_20, 
                                                                  effort[which(effort$Period == "2020"), 
                                                                         coordLonLat], 
                                                                  method = "simple", df = T)[,2] 
    effort[which(effort$Period == "2021"), 
           paste(variable, "sd", sep = "_")] <- raster::extract(ras_21, 
                                                                  effort[which(effort$Period == "2021"), 
                                                                         coordLonLat], 
                                                                  method = "simple", df = T)[,2] 
    return(effort)
  }
}

for (j in 1:length(var_list)){
data <- extractVar(effort = data, coordLonLat = c(9,10), variable = var_list[j], what = "Mean")
data <- extractVar(effort = data, coordLonLat = c(9,10), variable = var_list[j], what = "SD")
}



### 2 - Rescale and saturate ---------------------------------------------------
# functions
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

# prep data
var_to_test <- names(data)[c(19:27, 33:68)]
seg <- data

# gradMLD_avg_dist_to_front
boxplot(seg[, "gradMLD_avg_dist_to_front" ])

# gradTemp_avg_dist_to_front
boxplot(seg[, "gradTemp_avg_dist_to_front" ])
seg[, "gradTemp_avg_dist_to_front"] <- saturate(seg[, "gradTemp_avg_dist_to_front"], both = T, alpha = 0.04) 

# gradChl_avg_dist_to_front
boxplot(seg[, "gradChl_avg_dist_to_front" ])

# gradNPP_avg_dist_to_front
boxplot(seg[, "gradNPP_avg_dist_to_front" ])
seg[, "gradNPP_avg_dist_to_front"] <- saturate(seg[, "gradNPP_avg_dist_to_front"], both = T, alpha = 0.03) 

# gradMLD_front_persistence
boxplot(seg[, "gradMLD_front_persistence" ])
seg[, "gradMLD_front_persistence"] <- saturate(seg[, "gradMLD_front_persistence"], both = T, alpha = 0.02) 

# gradTemp_front_persistence no need
boxplot(seg[, "gradTemp_front_persistence" ])

# gradChl_front_persistence  no need
boxplot(seg[, "gradChl_front_persistence" ])

# gradNPP_front_persistence  no need
boxplot(seg[, "gradNPP_front_persistence" ])
seg[, "gradNPP_front_persistence"] <- saturate(seg[, "gradNPP_front_persistence"], both = T, alpha = 0.04) 

# Bathy
boxplot(seg[, "Bathy" ]) 

# Temp
boxplot(seg[, c("Temp_mean", "Temp_sd") ])
seg[, "Temp_mean"] <- saturate(seg[, "Temp_mean"], lower = T, alpha = 0.03) 
seg[, "Temp_sd"] <- saturate(seg[, "Temp_sd"], both = T, alpha = 0.04) 

# Salinity
boxplot(seg[, c("Salinity_mean", "Salinity_sd") ])
seg[, "Salinity_mean"] <- saturate(seg[, "Salinity_mean"], lower = T, alpha = 0.05) 

# CurrentSpeed
boxplot(seg[, c("CurrentSpeed_mean", "CurrentSpeed_sd") ])
seg[, "CurrentSpeed_mean"] <- saturate(seg[, "CurrentSpeed_mean"], both = T, alpha = 0.07) 
seg[, "CurrentSpeed_sd"] <- saturate(seg[, "CurrentSpeed_sd"], both = T, alpha = 0.05) 

# Chl
boxplot(seg[, c("Chl_mean", "Chl_sd") ])
seg[, "Chl_mean"] <- saturate(seg[, "Chl_mean"], both = T, alpha = 0.05) 
seg[, "Chl_sd"] <- saturate(seg[, "Chl_sd"], both = T, alpha = 0.05) 

# NPP
boxplot(seg[, c("NPP_mean", "NPP_sd") ])
seg[, "NPP_mean"] <- saturate(seg[, "NPP_mean"], both = T, alpha = 0.05) 
seg[, "NPP_sd"] <- saturate(seg[, "NPP_sd"], both = T, alpha = 0.05) 

# Phyto
boxplot(seg[, c("Phyto_mean", "Phyto_sd") ])
seg[, "Phyto_mean"] <- saturate(seg[, "Phyto_mean"], both = T, alpha = 0.05) 
seg[, "Phyto_sd"] <- saturate(seg[, "Phyto_sd"], both = T, alpha = 0.02) 

# EKE
boxplot(seg[, c("EKE_mean", "EKE_sd") ])
seg[, "EKE_mean"] <- saturate(seg[, "EKE_mean"], both = T, alpha = 0.08) 
seg[, "EKE_sd"] <- saturate(seg[, "EKE_sd"], both = T, alpha = 0.075) 

# gradSSH
boxplot(seg[, c("gradSSH_mean", "gradSSH_sd") ])
seg[, "gradSSH_mean"] <- saturate(seg[, "gradSSH_mean"], both = T, alpha = 0.05) 
seg[, "gradSSH_sd"] <- saturate(seg[, "gradSSH_sd"], both = T, alpha = 0.05) 

# gradTemp
boxplot(seg[, c("gradTemp_mean", "gradTemp_sd") ])
seg[, "gradTemp_mean"] <- saturate(seg[, "gradTemp_mean"], both = T, alpha = 0.05) 
seg[, "gradTemp_sd"] <- saturate(seg[, "gradTemp_sd"], both = T, alpha = 0.05) 

# gradMLD
boxplot(seg[, c("gradMLD_mean", "gradMLD_sd") ])
seg[, "gradMLD_mean"] <- saturate(seg[, "gradMLD_mean"], both = T, alpha = 0.02) 

# gradSal
boxplot(seg[, c("gradSal_mean", "gradSal_sd") ])
seg[, "gradSal_mean"] <- saturate(seg[, "gradSal_mean"], both = T, alpha = 0.06) 
seg[, "gradSal_sd"] <- saturate(seg[, "gradSal_sd"], both = T, alpha = 0.03) 

# SSH
boxplot(seg[, c("SSH_mean", "SSH_sd") ])
seg[, "SSH_mean"] <- saturate(seg[, "SSH_mean"], both = T, alpha = 0.02) 
seg[, "SSH_sd"] <- saturate(seg[, "SSH_sd"], both = T, alpha = 0.02) 

# gradChl
boxplot(seg[, c("gradChl_mean", "gradChl_sd") ])
seg[, "gradChl_mean"] <- saturate(seg[, "gradChl_mean"], both = T, alpha = 0.04) 
seg[, "gradChl_sd"] <- saturate(seg[, "gradChl_sd"], both = T, alpha = 0.04) 

# gradNPP
boxplot(seg[, c("gradNPP_mean", "gradNPP_sd") ])
seg[, "gradNPP_mean"] <- saturate(seg[, "gradNPP_mean"], both = T, alpha = 0.04) 
seg[, "gradNPP_sd"] <- saturate(seg[, "gradNPP_sd"], both = T, alpha = 0.05) 

# ZEu no need
boxplot(seg[, c("ZEu_mean", "ZEu_sd") ])

# MLD no need
boxplot(seg[, c("MLD_mean", "MLD_sd") ])

# DissIC
boxplot(seg[, c("DissIC_mean", "DissIC_sd") ])
seg[, "DissIC_mean"] <- saturate(seg[, "DissIC_mean"], lower = T, alpha = 0.05) 

# SPCO2
boxplot(seg[, c("SPCO2_mean", "SPCO2_sd") ])
seg[, "SPCO2_mean"] <- saturate(seg[, "SPCO2_mean"], both = T, alpha = 0.05) 
seg[, "SPCO2_sd"] <- saturate(seg[, "SPCO2_sd"], both = T, alpha = 0.05) 

# plot(gam(n_ind~s(gradNPP_mean, k = 4, bs = 'tp'), data = seg), rug =T)




### 3 - Prediction grid --------------------------------------------------------
raster_ref <- raster(paste(sourcedir, "/Variables/gradTemp_frontal_persistence-2019-12-01_2020-04-13.envi", sep =""))
bathy_rspl <- raster::resample(x = bathy, y = raster_ref)

# table 2020
pred_tab_20 <- cbind(as.data.frame(bathy_rspl, xy = T),
                       
                       as.data.frame(raster(paste(datadir, "/gradChl_frontal_persistence-2019-12-01_2020-04-13.envi", sep =""))),
                       as.data.frame(raster(paste(datadir, "/gradChl_front_distance-averaged-2019-12-01_2020-04-13.envi", sep ="")))*202,
                       
                       as.data.frame(raster(paste(datadir, "/gradNPP_frontal_persistence-2019-12-01_2020-04-13.envi", sep =""))),
                       as.data.frame(raster(paste(datadir, "/gradNPP_front_distance-averaged-2019-12-01_2020-04-13.envi", sep ="")))*202,
                       
                       as.data.frame(raster(paste(datadir, "/gradMLD_frontal_persistence-2019-12-01_2020-04-13.envi", sep =""))),
                       as.data.frame(raster(paste(datadir, "/gradMLD_front_distance-averaged-2019-12-01_2020-04-13.envi", sep ="")))*202,
                       
                       as.data.frame(raster(paste(datadir, "/gradTemp_frontal_persistence-2019-12-01_2020-04-13.envi", sep =""))),
                       as.data.frame(raster(paste(datadir, "/gradTemp_front_distance-averaged-2019-12-01_2020-04-13.envi", sep ="")))*202
)
names(pred_tab_20) <- c("lon", "lat", "Bathy", 
                          "gradChl_front_persistence", "gradChl_avg_dist_to_front",
                          "gradNPP_front_persistence", "gradNPP_avg_dist_to_front",
                          "gradMLD_front_persistence", "gradMLD_avg_dist_to_front",
                          "gradTemp_front_persistence", "gradTemp_avg_dist_to_front")
for(i in 1:length(var_list)){
  # mean
  ras_mean <- get(paste(var_list[i], "Mean", sep = "_"))
  pred_tab_20[, paste(var_list[i], "mean", sep = "_")] <- as.data.frame(ras_mean[["X2020"]], xy = F)
  # sd
  ras_sd <- get(paste(var_list[i], "SD", sep = "_"))
  pred_tab_20[, paste(var_list[i], "sd", sep = "_")] <- as.data.frame(ras_sd[["X2020"]], xy = F)
}


# table 2021
pred_tab_21 <- cbind(as.data.frame(bathy_rspl, xy = T),
                       
                       as.data.frame(raster(paste(datadir, "/gradChl_frontal_persistence-2020-12-01_2021-04-13.envi", sep =""))),
                       as.data.frame(raster(paste(datadir, "/gradChl_front_distance-averaged-2020-12-01_2021-04-13.envi", sep ="")))*201,
                       
                       as.data.frame(raster(paste(datadir, "/gradNPP_frontal_persistence-2020-12-01_2021-04-13.envi", sep =""))),
                       as.data.frame(raster(paste(datadir, "/gradNPP_front_distance-averaged-2020-12-01_2021-04-13.envi", sep ="")))*201,
                       
                       as.data.frame(raster(paste(datadir, "/gradMLD_frontal_persistence-2020-12-01_2021-04-13.envi", sep =""))),
                       as.data.frame(raster(paste(datadir, "/gradMLD_front_distance-averaged-2020-12-01_2021-04-13.envi", sep ="")))*201,
                       
                       as.data.frame(raster(paste(datadir, "/gradTemp_frontal_persistence-2020-12-01_2021-04-13.envi", sep =""))),
                       as.data.frame(raster(paste(datadir, "/gradTemp_front_distance-averaged-2020-12-01_2021-04-13.envi", sep ="")))*201
)
names(pred_tab_21) <- c("lon", "lat", "Bathy", 
                          "gradChl_front_persistence", "gradChl_avg_dist_to_front",
                          "gradNPP_front_persistence", "gradNPP_avg_dist_to_front",
                          "gradMLD_front_persistence", "gradMLD_avg_dist_to_front",
                          "gradTemp_front_persistence", "gradTemp_avg_dist_to_front")
for(i in 1:length(var_list)){
  # mean
  ras_mean <- get(paste(var_list[i], "Mean", sep = "_"))
  pred_tab_21[, paste(var_list[i], "mean", sep = "_")] <- as.data.frame(ras_mean[["X2021"]], xy = F)
  # sd
  ras_sd <- get(paste(var_list[i], "SD", sep = "_"))
  pred_tab_21[, paste(var_list[i], "sd", sep = "_")] <- as.data.frame(ras_sd[["X2021"]], xy = F)
}




#### ------------------------- ####
#### ------- MODELLING ------- ####
#### ------------------------- ####
# functions
source(paste(funcdir, "fctGAM_20210316.R", sep = "/"))
correlations <- function(predictors, effort) {
  corr <- Hmisc::rcorr(as.matrix(effort[!is.na(effort[,predictors]), predictors]), type = "spearman")
  # symnum(corr, abbr.colnames=FALSE)
  corrplot::corrplot(corr$r, type="upper", order="hclust", tl.col="black", tl.srt=45, method = "number", diag = F,
                     addCoefasPercent = TRUE, tl.cex = 0.7, is.corr = TRUE, number.cex = .5)
}


### 1 - Inspect correlations ----------------------------------------------------
correlations(effort = seg, predictors = var_to_test)



### 2 - Fit all models ---------------------------------------------------------
## Fitting
fit_best <- fit_all_gam(envdata = seg, esw = seg$ESW, outcome = "n_ind",
                        predictors = var_to_test,
                        rescale = F, weight_by_g_naught = 1, family = "tweedie", 
                        max_cor = 0.6, nb_min_pred = 1, nb_max_pred = 4, complexity = 4)
fit_best_ord <- fit_best[order(fit_best$AIC),]
aics <- akaike.weights(fit_best_ord[,3])
fit_best_ord$Delta.AIC <- aics$deltaAIC
fit_best_ord$rel.likelihood <- aics$rel.LL 
fit_best_ord$Akaike.weight <- aics$weights

## Plots
plot(fit_best_ord$Akaike.weight, xlim = c(0, 50))
fit_best_ord$Number <- 1:nrow(fit_best_ord)
ggplot(fit_best_ord) + geom_point(aes(y = Akaike.weight, x = Number)) + xlim(0, 50) + theme_bw()
ggplot(fit_best_ord) + geom_point(aes(y = Delta.AIC, x = Number)) + xlim(0, 50) + theme_bw()



### 3 - Average AIC weight of variables ----------------------------------------
library(data.table)
table <- data.frame(var = NA, count = NA, percent = NA, akaike.weight = NA)
for (i in 1:length(var_to_test)) {
  select <- fit_best_ord[fit_best_ord$model %like% var_to_test[i], ]
  table[i,1] <- var_to_test[i]
  table[i,2] <- nrow(select) 
  table[i,3] <- nrow(select)/nrow(fit_best_ord)*100
  table[i,4] <- sum(select$Akaike.weight)*100
}
table <- table[order(-table$akaike.weight),]

plot(table$akaike.weight)



### 4 - Inspect best models ----------------------------------------------------
mod1 <- gam(n_ind~ 1 +s(Temp_sd  , k = 4, bs = 'tp') + 
              s(gradTemp_sd   , k = 4, bs = 'tp') + 
              s(gradMLD_avg_dist_to_front, k = 4, bs = 'tp') + 
              s(ZEu_sd, k = 4, bs = 'tp')
            , data = seg, offset = log(seg$Effort*2*seg$ESW), family = tw())
plot(mod1, pages = 1, scale = 0, res = T)
summary(mod1)
AIC(mod1) 



### 5 - Predict from best model ------------------------------------------------
# generate mask to clip on shelf
masking_layer <- raster(paste(sourcedir, "/Variables/Bathy.tif", sep = ""))
masking_layer <- crop(masking_layer, extent(-4.8,-1.2,44.3,47.1))
masking_layer <- resample(masking_layer, raster(paste(sourcedir, "/Variables/gradTemp_frontal_persistence-2019-12-01_2020-04-13.envi", sep ="")))

masking_layer[masking_layer>0] <- NA
masking_layer[masking_layer<(-2500)] <- NA

masking_layer <- mask(masking_layer, land, inverse = T)



# fit prediction to all models to explore
fit_pred <- function(mod, newdata, se.fit = T){
  if(se.fit == T){
    pred <- predict.gam(object = mod, newdata = newdata, type = "response", se.fit = T)
    
    newdata$fit <- pred$fit
    pred_fit <- rasterFromXYZ(newdata[, c("lon", "lat", "fit")])
    pred_fit <- mask(pred_fit, masking_layer)
    
    newdata$se.fit <- pred$se.fit
    pred_sefit <- rasterFromXYZ(newdata[, c("lon", "lat", "se.fit")])
    pred_sefit <- mask(pred_sefit, masking_layer)
    return(list(fit = pred_fit, se = pred_sefit))
  }
  if(se.fit == F){
    pred <- predict.gam(object = mod, newdata = newdata, type = "response", se.fit = F)
    
    newdata$fit <- pred
    pred_fit <- rasterFromXYZ(newdata[, c("lon", "lat", "fit")])
    pred_fit <- mask(pred_fit, masking_layer)
    return(pred_fit)
  }
}

par(mfrow = c(1,2), mar = c(1.5,1.5,1.5,1.5))
plot(fit_pred(mod = mod1, newdata = pred_tab_20, se.fit = F), 
     col = viridis(20), main = "mod1",
     addfun = c(plot(land, add = T, col = "grey")))
plot(fit_pred(mod = mod1, newdata = pred_tab_21, se.fit = F), 
     col = viridis(20), main = "mod1",
     addfun = c(plot(land, add = T, col = "grey")))


# extrat prediction et incertitude du best model
pred_20 <- fit_pred(mod = mod1, newdata = pred_tab_20, se.fit = T)
pred_20_fit <- pred_20$fit
pred_20_sefit <- pred_20$se

pred_21 <- fit_pred(mod = mod1, newdata = pred_tab_21, se.fit = T)
pred_21_fit <- pred_21$fit
pred_21_sefit <- pred_21$se


### 6 - Extraction des elements a utiliser pour la suite -----------------------
## Prediction
png(paste(outdir, "SeasonalModel_prediction.png", sep = "/"), height = 15, width = 15, units = "cm", res = 300)
par(mfrow = c(2,2), mar = c(0.8,1,2,0.3),oma = c(0,0,0.5,0))

plot(pred_20_fit, col = viridis(20),  axes = F, xlim = c(-4.8, -1),
     addfun = c(plot(land, add = T, col = "grey"), box(),
                scalebar(type = "bar", below = "kilom.", cex = 0.7, d = 30)))
addnortharrow(scale = 0.4, padin = c(0.05, 0.05))
mtext("Mean prediction", side = 3, adj = 0, cex = 0.9)
mtext("A - 2020", side = 3, adj = 0, cex = 0.9, line = 1.2, font = 2, col = "#da4308")
plot(pred_20_sefit, col = viridis(20),  axes = F, xlim = c(-4.8, -1),
     addfun = c(plot(land, add = T, col = "grey"), box(),
                scalebar(type = "bar", below = "kilom.", cex = 0.7, d = 30)))
addnortharrow(scale = 0.4, padin = c(0.05, 0.05))
mtext("Uncertainty", side = 3, adj = 0, cex = 0.9)

plot(pred_21_fit, col = viridis(20),  axes = F, xlim = c(-4.8, -1),
     addfun = c(plot(land, add = T, col = "grey"), box(),
                scalebar(type = "bar", below = "kilom.", cex = 0.7, d = 30)))
addnortharrow(scale = 0.4, padin = c(0.05, 0.05))
mtext("Mean prediction", side = 3, adj = 0, cex = 0.9)
mtext("B - 2021", side = 3, adj = 0, cex = 0.9, line = 1.2, font = 2, col = "#da4308")
plot(pred_21_sefit, col = viridis(20),  axes = F, xlim = c(-4.8, -1),
     addfun = c(plot(land, add = T, col = "grey"), box(),
                scalebar(type = "bar", below = "kilom.", cex = 0.7, d = 30)))
addnortharrow(scale = 0.4, padin = c(0.05, 0.05))
mtext("Uncertainty", side = 3, adj = 0, cex = 0.9)
dev.off()

## pred 2
png(paste(outdir, "SeasonalModel_prediction2.png", sep = "/"), height = 7, width = 15, units = "cm", res = 300)
par(mfrow = c(1,2), mar = c(0.8,1,1.2,0.3))
plot(pred_20_fit, col = viridis(20),  axes = F, xlim = c(-4.8, -1),
     addfun = c(plot(land, add = T, col = "grey"), box(),
                scalebar(type = "bar", below = "kilom.", cex = 0.7, d = 30)))
addnortharrow(scale = 0.4, padin = c(0.05, 0.05))
mtext("A - 2020", side = 3, adj = 0, cex = 0.9, font = 2)

plot(pred_21_fit, col = viridis(20),  axes = F, xlim = c(-4.8, -1),
     addfun = c(plot(land, add = T, col = "grey"), box(),
                scalebar(type = "bar", below = "kilom.", cex = 0.7, d = 30)))
addnortharrow(scale = 0.4, padin = c(0.05, 0.05))
mtext("B - 2021", side = 3, adj = 0, cex = 0.9, font = 2)
dev.off()


## exporte pred
writeRaster(x = pred_20_fit, f = paste(outdir, "/SeasonalModel_2020_prediction.tif", sep = ""), format = "GTiff", overwrite = T)
writeRaster(x = pred_20_sefit, f = paste(outdir, "/SeasonalModel_2020_SE.tif", sep = ""), format = "GTiff", overwrite = T)

writeRaster(x = pred_21_fit, f = paste(outdir, "/SeasonalModel_2021_prediction.tif", sep = ""), format = "GTiff", overwrite = T)
writeRaster(x = pred_21_sefit, f = paste(outdir, "/SeasonalModel_2021_SE.tif", sep = ""), format = "GTiff", overwrite = T)

## Model
png(paste(outdir, "SeasonalModel_response_curves.png", sep = "/"), height = 5, width = 20, res = 400, units = "cm")
par(mfrow = c(1,4), mar = c(2,2,2,1), oma = c(0,1.5,0.8,0.5), mgp = c(2.5,0.5,0), tck = -0.04)
plot_smooth(mod1, view = "gradMLD_avg_dist_to_front", sim.ci = T, set.seed(1e4), ylab = "", main = "", hide.label = T, h0 = 1, transform = exp) ; box()
mtext(expression(ind.km^2), side = 2, line = 2, cex = 0.8)
mtext("Average distance\n to MLD front", side = 3, line = 0.3, cex = 0.8)
plot_smooth(mod1, view = "Temp_sd", sim.ci = T, set.seed(1e4), ylab = "", main = "", hide.label = T, h0 = 1, transform = exp) ; box() #, transform = exp
mtext("Temp (sd)", side = 3, line = 0.3, cex = 0.8)
plot_smooth(mod1, view = "gradTemp_sd", sim.ci = T, set.seed(1e4), ylab = "", main = "", hide.label = T, h0 = 1, transform = exp) ; box()
mtext("Temp front (sd)", side = 3, line = 0.3, cex = 0.8)
plot_smooth(mod1, view = "ZEu_sd", sim.ci = T, set.seed(1e4), ylab = "", main = "", hide.label = T, h0 = 1, transform = exp) ; box()
mtext("ZEu (sd)", side = 3, line = 0.3, cex = 0.8)
dev.off()




#### ----------------------------- ####
#### ------- EVALUATE PRED ------- ####
#### ----------------------------- ####
source("F:/CAPECET/functions/make_cfact_2.R")
source("F:/CAPECET/functions/rescale2.R")
library(assertthat)
library(WhatIf)

### 1 - Quantify the proportion of extrapolation in prediction grids -----
covar <- c("gradMLD_avg_dist_to_front", "gradTemp_sd", "Temp_sd", "ZEu_sd")

# 2020
newdata20 <- pred_tab_20
newdata20 <- newdata20[,c("lon", "lat", covar, "Bathy")]
newdata20[which(newdata20$Bathy > 0),] <- NA # masking layer (sur bathy et sur long pour enlever estuaire)
newdata20[which(newdata20$Bathy < (-2500)),] <- NA
newdata20[which(newdata20$x > -1.2),] <- NA
newdata20 <- newdata20 %>% drop_na(all_of(covar)) # remove NA bc make_cfact_2 removes it
make_cfact_2(calibration_data = seg, 
             test_data = newdata20, 
             var_name = covar)
# 0.402


# 2021
newdata21 <- pred_tab_21
newdata21 <- newdata21[,c("lon", "lat", covar, "Bathy")]
newdata21[which(newdata21$Bathy > 0),] <- NA # masking layer (sur bathy et sur long pour enlever estuaire)
newdata21[which(newdata21$Bathy < (-2500)),] <- NA
newdata21[which(newdata21$x > -1.2),] <- NA
newdata21 <- newdata21 %>% drop_na(all_of(covar)) # remove NA bc make_cfact_2 removes it
make_cfact_2(calibration_data = seg, 
             test_data = newdata21, 
             var_name = covar)
# 0.251


### 3 - Map daily extrapolation ------------------------------------------------
## identify extrapolation points in daily prediction grids (map)
pts20 <- make_cfact_2(calibration_data = seg, 
                        test_data = newdata20, 
                        var_name = covar,
                        percent = FALSE)
extrap_map_20 <- rasterFromXYZ( cbind(newdata20[,c("lon", "lat")], pts20) )

pts21 <- make_cfact_2(calibration_data = seg, 
                        test_data = newdata21, 
                        var_name = covar,
                        percent = FALSE)
extrap_map_21 <- rasterFromXYZ( cbind(newdata21[,c("lon", "lat")], pts21) )

extrap_map <- stack(extrap_map_20, extrap_map_21)
names(extrap_map) <- c("2020", "2021")


### 4 - Export rasters ---------------------------------------------------------
writeRaster(x = extrap_map_20, f = paste(outdir, "/SeasonalModel_extrap_2020.tif", sep = ""), format = "GTiff", overwrite = T)
writeRaster(x = extrap_map_21, f = paste(outdir, "/SeasonalModel_extrap_2021.tif", sep = ""), format = "GTiff", overwrite = T)


### 5 - Plots totaux -----------------------------------------------------------

png(paste(outdir, "Seasonal_model_maps.png", sep = "/"), height = 12, width = 18.5, units = "cm", res = 400)
par(mfrow = c(2,3), mar = c(0.8,1.5,3,1),oma = c(0,0,0,1))
# 2020
plot(pred_20_fit, col = viridis(20),  axes = F, xlim = c(-4.8, -1), ylim = c(44.4,47),
     addfun = c(plot(land, add = T, col = "grey"),box(),
                scalebar(type = "bar", below = "kilometers", cex = 0.8, d = 50)))
addnortharrow(scale = 0.4, padin = c(0.05, 0.05))
mtext("A - 2020", side = 3, line = 1.5, font = 2, col = "#da4308", adj = 0, at = -5.1)
mtext("Prediction", side = 3, adj = 0, cex = 0.95)
plot(pred_20_sefit, col = viridis(20, option = "C"),  axes = F, xlim = c(-4.8, -1), ylim = c(44.4,47),
     addfun = c(plot(land, add = T, col = "grey"),box(),
                scalebar(type = "bar", below = "kilometers", cex = 0.8, d = 50)))
addnortharrow(scale = 0.4, padin = c(0.05, 0.05))
mtext("Uncertainty", side = 3, adj = 0, cex = 0.95)
plot(extrap_map_20, col = c("grey60", "grey20"),  axes = F, xlim = c(-4.8, -1), ylim = c(44.4,47), legend = F,
     addfun = c(plot(land, add = T, col = "grey"),box(),
                scalebar(type = "bar", below = "kilometers", cex = 0.8, d = 50)))
legend("topright", cex = 1, legend = c("0", "1"), fill = c("grey60", "grey20"))
mtext("Extrapolation", side = 3, adj = 0, cex = 0.95)
# 2021
plot(pred_21_fit, col = viridis(20),   axes = F, xlim = c(-4.8, -1), ylim = c(44.4,47),
     addfun = c(plot(land, add = T, col = "grey"),box(),
                scalebar(type = "bar", below = "kilometers", cex = 0.8, d = 50)))
addnortharrow(scale = 0.4, padin = c(0.05, 0.05))
mtext("B - 2021", side = 3, line = 1.5, font = 2, col = "#da4308", adj = 0, at = -5.1)
mtext("Prediction", side = 3, adj = 0, cex = 0.95)
plot(pred_21_sefit, col = viridis(20, option = "C"),   axes = F, xlim = c(-4.8, -1), ylim = c(44.4,47),
     addfun = c(plot(land, add = T, col = "grey"), box(),
                scalebar(type = "bar", below = "kilometers", cex = 0.8, d = 50)))
addnortharrow(scale = 0.4, padin = c(0.05, 0.05))
mtext("Uncertainty", side = 3, adj = 0, cex = 0.95)
plot(extrap_map_21, col = c("grey60", "grey20"), axes = F, xlim = c(-4.8, -1), ylim = c(44.4,47), legend = F,
     addfun = c(plot(land, add = T, col = "grey"),box(),
                scalebar(type = "bar", below = "kilometers", cex = 0.8, d = 50)))
legend("topright", cex = 1, legend = c("0", "1"), fill = c("grey60", "grey20"))
mtext("Extrapolation", side = 3, adj = 0, cex = 0.95)
dev.off()


