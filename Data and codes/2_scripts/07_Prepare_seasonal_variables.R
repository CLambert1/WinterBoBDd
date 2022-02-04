##--------------------------------------------------------------------------------------------------------
## PROJECT : Variation hivernale de la distribution du dauphin commun dans le GoG
## SCRIPT : Compute prediction tables for the seasonal model
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

#### ----------------------- ####
#### ---- SOURCES FILES ---- ####
#### ----------------------- ####
outdir <- "F:/CAPECET/paper/data/3_script_results"
sourcedir <- "F:/CAPECET/paper/data/1_data"
setwd(outdir)

### 1 - Open map files
## open bathy
bathy <- raster(paste(sourcedir, "/Variables/Bathy.tif", sep = ""))
bathy <- crop(bathy, extent(-4.8,0.5,44.3,47.1))


### 2 - Dates 
dates_1920 <- ymd("2019-12-01"):ymd("2020-04-13")
dates_1920 <- as_date(dates_1920)
int_1920 <- interval(ymd("2019-12-01"), ymd("2020-04-13"))

dates_2021 <- ymd("2020-12-01"):ymd("2021-04-13")
dates_2021 <- as_date(dates_2021)
int_2021 <- interval(ymd("2020-12-01"), ymd("2021-04-13"))



#### ------------------------- ####
#### --- Average variables --- ####
#### ------------------------- ####


var_list <- c("Temp", "Salinity", "CurrentSpeed", "Chl", "NPP", "Phyto",
"EKE", "gradSSH", "gradTemp", "gradMLD", "gradSal",
"SSH", "gradChl", "gradNPP", "ZEu", "MLD", "DissIC","SPCO2",
"gradMLD_front_distance", "gradTemp_front_distance", "gradChl_front_distance", "gradNPP_front_distance",
"gradMLD_front_persistence", "gradTemp_front_persistence", "gradChl_front_persistence", "gradNPP_front_persistence")


## functions
mean_ras <- function(variable, what = "Rasters", year){
  if(year == 2020){
    rasters <- stack(paste(datadir, "/", variable, "-daily-2019-12-01_2020-04-13.envi", sep =""))
    if(variable == "SPCO2") rasters <- rasters[[which(as_date(str_sub(names(rasters), 6)) %within% int_1920)]]
    if(variable != "SPCO2") rasters <- rasters[[which(as_date(str_sub(names(rasters), 5)) %within% int_1920)]]
    
    if(what == "Rasters"){
      avg <- overlay(rasters, fun = "mean", na.rm = T) # average in each cell across layers
      sd <- overlay(rasters, fun = "sd", na.rm = T)
      return(list(mean = avg, sd = sd))
    }
    if(what == "TimeSeries"){
      Zavg <- cellStats(rasters, stat = "mean", na.rm = T)   # summary stat of each layer (each date)
      Zsd <- cellStats(rasters, stat = "sd", na.rm = T) 
      ts <- as.data.frame(cbind(dates_1920, Zavg, Zsd))
      names(ts) <- c("Date", "Zavg", "Zsd")
      ts$Variable <- variable
      return(ts = ts)
    }
  }
  if(year == 2021){
    rasters <- stack(paste(datadir, "/", variable, "-daily-2020-12-01_2021-04-13.envi", sep =""))
    if(variable == "SPCO2") rasters <- rasters[[which(as_date(str_sub(names(rasters), 6)) %within% int_2021)]]
    if(variable != "SPCO2") rasters <- rasters[[which(as_date(str_sub(names(rasters), 5)) %within% int_2021)]]
    
    if(what == "Rasters"){
      avg <- overlay(rasters, fun = "mean", na.rm = T) # average in each cell across layers
      sd <- overlay(rasters, fun = "sd", na.rm = T)
      return(list(mean = avg, sd = sd))
    }
    if(what == "TimeSeries"){
      Zavg <- cellStats(rasters, stat = "mean", na.rm = T)   # summary stat of each layer (each date)
      Zsd <- cellStats(rasters, stat = "sd", na.rm = T) 
      ts <- as.data.frame(cbind(dates_2021, Zavg, Zsd))
      names(ts) <- c("Date", "Zavg", "Zsd")
      ts$Variable <- variable
      return(ts = ts)
    }
  }
}


## summarize rasters
# all vars
for(i in 1:22){
  smry_20 <- mean_ras(variable = var_list[i], what = "Rasters", year = 2020)
  assign(paste("map", var_list[i], "20", sep = "_"), smry_20)

  smry_21 <- mean_ras(variable = var_list[i], what = "Rasters", year = 2021)
  assign(paste("map", var_list[i], "21", sep = "_"), smry_21)
}

# front persistence
Chl_persist_20 <- stack(paste(datadir, "/", "gradChl_frontal_persistence-2019-12-01_2020-04-13.envi", sep =""))
NPP_persist_20 <- stack(paste(datadir, "/", "gradNPP_frontal_persistence-2019-12-01_2020-04-13.envi", sep =""))
MLD_persist_20 <- stack(paste(datadir, "/", "gradMLD_frontal_persistence-2019-12-01_2020-04-13.envi", sep =""))
Temp_persist_20 <- stack(paste(datadir, "/", "gradTemp_frontal_persistence-2019-12-01_2020-04-13.envi", sep =""))

Chl_persist_21 <- stack(paste(datadir, "/", "gradChl_frontal_persistence-2020-12-01_2021-04-13.envi", sep =""))
NPP_persist_21 <- stack(paste(datadir, "/", "gradNPP_frontal_persistence-2020-12-01_2021-04-13.envi", sep =""))
MLD_persist_21 <- stack(paste(datadir, "/", "gradMLD_frontal_persistence-2020-12-01_2021-04-13.envi", sep =""))
Temp_persist_21 <- stack(paste(datadir, "/", "gradTemp_frontal_persistence-2020-12-01_2021-04-13.envi", sep =""))


## summarize time series
for(i in 1:22){
  smry_20 <- mean_ras(variable = var_list[i], what = "TimeSeries", year = 2020)
  assign(paste("ts", var_list[i], "20", sep = "_"), smry_20)

  smry_21 <- mean_ras(variable = var_list[i], what = "TimeSeries", year = 2021)
  assign(paste("ts", var_list[i], "21", sep = "_"), smry_21)
}


save.image(paste(outdir, "/Mean_covariates.RData", sep = ""))




