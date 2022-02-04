##--------------------------------------------------------------------------------------------------------
## PROJECT : Infra-seasonal variations of common dolphin distribution in the Bay of Biscay
## SCRIPT : Predict from the daily model
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


#### ----------------------- ####
#### ---- SOURCES FILES ---- ####
#### ----------------------- ####

### 1 - Load geographic files --------------------------------------------------------
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


### 2 - Load seg file and fit model --------------------------------------------------
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


mod1 <- gam(n_ind ~ 1 +s(gradMLD_avg_dist_to_front, k = 4, bs = 'tp') +
              s(MLD_10d, k = 4, bs = 'tp') +
              s(Chl_7d, k = 4, bs = 'tp') +
              s(gradChl_front_distance_2d, k = 4, bs = 'tp')
            , data = seg, offset = log(seg$Effort*2*seg$ESW), family = tw())
covar <- c("gradMLD_avg_dist_to_front", "MLD_10d", "Chl_7d", "gradChl_front_distance_2d")


### 3 - Date vector ------------------------------------------------------------
dates <- as_date(c(ymd("2020-01-01"):ymd("2020-04-13"), ymd("2021-01-01"):ymd("2021-04-13")))
dates <- as_date(dates)


### 4 - Generate mask to clip on shelf -----------------------------------------
masking_layer <- raster(paste(sourcedir, "/Variables/Bathy.tif", sep = ""))
masking_layer <- crop(masking_layer, extent(-4.8,-1.2,44.3,47.1))
masking_layer <- resample(masking_layer, raster(paste(sourcedir, "/Variables/gradTemp_frontal_persistence-2019-12-01_2020-04-13.envi", sep ="")))

masking_layer[masking_layer>0] <- NA
masking_layer[masking_layer<(-2500)] <- NA

masking_layer <- mask(masking_layer, land, inverse = T)


### 5 - Open CDS table ---------------------------------------------------------
cds <- read.table(paste(outdir, "Encounter_rate_per_day_2Years.txt", sep = "/"), h = T) %>%
  subset(Year != 2019) %>%
  mutate(date2 = case_when( # trick for plotting on the same levels
      Period == "2020" ~Day,
      Period == "2021" ~as_date(Day - years(1))
    )
  )    
dates_ech <- cds[,c("Day", "date2", "Period")]



#### --------------------------------- ####
#### ---- PREDICTION FULL DATASET ---- ####
#### --------------------------------- ####

### 1 - Create functions -------------------------------------------------------
prediction_stack <- function(model, dates, datadir){
  list_pred <- pbapply::pblapply(1:length(dates), function(i){
    newdata <- read.table(paste(sourcedir, "/Variables/PredictionGrids/PredictionGrid_", dates[i], ".txt", sep = ""), h = T, sep = "\t", dec = ".")
    pred <- predict.gam(object = model, newdata = newdata, type = "response", se.fit = T)
    
    newdata$fit <- pred$fit
    pred_fit <- rasterFromXYZ(newdata[, c("x", "y", "fit")])
    names(pred_fit) <- paste("fit", dates[i], sep = "_")
    pred_fit <- mask(x = pred_fit, mask = masking_layer)
    
    newdata$se.fit <- pred$se.fit
    pred_sefit <- rasterFromXYZ(newdata[, c("x", "y", "se.fit")])
    names(pred_sefit) <- paste("se.fit", dates[i], sep = "_")
    pred_sefit <- mask(x = pred_sefit, mask = masking_layer)
    
    list(pred_fit, pred_sefit)
  })
  pred_stack <- stack(sapply(list_pred, "[[", 1))
  se_stack <- stack(sapply(list_pred, "[[", 2))
  
  return(list(prediction_fit = pred_stack,
              prediction_se = se_stack)) 
}

writeENVI <- function(x, f, ...) {
  writeRaster(x, f, overwrite=TRUE, ...)
  cat("band names = {", paste(names(x),collapse=","), "}", "\n", file=extension(f, "hdr"), append=TRUE)
}


### 2 - Predict mean and se ----------------------------------------------------
## daily pred 
raw_stack <- prediction_stack(model = mod1, dates = dates, datadir = datadir)

pred_stack <- raw_stack[["prediction_fit"]]
pred_stack[(pred_stack <0)] <- 0
writeENVI(x = pred_stack, f = paste(outdir, "/DailyModel_Predictions.envi", sep =""))

se_stack <- raw_stack[["prediction_se"]]
writeENVI(x = se_stack, f = paste(outdir, "/DailyModel_PredictionsSE.envi", sep =""))

save.image(paste(outdir, "/DailyModel_Prediction.RData", sep = ""))


#### --------------------------- ####
#### ---- PLOTS PREDICTIONS ---- ####
#### --------------------------- ####

## 1 - reopen raster stacks ----------------------------------------------------
pred_stack <- stack(paste(outdir, "/DailyModel_Predictions.envi", sep =""))
se_stack <- stack(paste(outdir, "/DailyModel_PredictionsSE.envi", sep =""))



### 2 - Time series plots  -----------------------------------------------------
## a - generate pred in abundance
surf <- raster::area(pred_stack[[1]])
pred_stack_ab <- pred_stack * surf
names(pred_stack_ab) <- names(pred_stack)
surf_tot <- cellStats(surf, stat = "sum", na.rm = T)

## b - generate table with one value per day
pred_smry <- data.frame(date = dates) 
pred_smry <- pred_smry %>% 
  mutate(pred_smry,
         Mean = cellStats(pred_stack, stat = "mean", na.rm = T),
         SD = cellStats(pred_stack, stat = "sd", na.rm = T), 
         Min = cellStats(pred_stack, stat = "min", na.rm = T),
         Max = cellStats(pred_stack, stat = "max", na.rm = T),
         Sum = cellStats(pred_stack_ab, stat = "sum", na.rm = T)
  )%>% 
  pivot_longer(cols = c(2:6), names_to = "variable") %>%
  mutate(What = "Prediction")
se_smry <- data.frame(date = dates) 
se_smry <- se_smry %>% 
  mutate(se_smry,
         Mean = cellStats(se_stack, stat = "mean", na.rm = T),
         SD = cellStats(se_stack, stat = "sd", na.rm = T), 
         Min = cellStats(se_stack, stat = "min", na.rm = T),
         Max = cellStats(se_stack, stat = "max", na.rm = T)
  ) %>% 
  pivot_longer(cols = c(2:5), names_to = "variable") %>%
  mutate(What = "Uncertainty")
smry <- full_join(pred_smry, se_smry)
write.table(smry, paste(outdir, "DailyModel_Pred_TimeSeries.csv", sep = "/"), row.names = F)


Sys.setlocale("LC_TIME", "English")


## c - join tables
ts <-  smry %>% mutate(
  Year = year(date),
  Period = case_when(
    date %within% int_1920 ~"2020",
    date %within% int_2021 ~"2021"
  ),
  date2 = case_when(
    Period == "2020" ~date,
    Period == "2021" ~as_date(date - years(1))
  )
)

## d - pivot wider
ts2 <-  ts %>% 
  pivot_wider(values_from = value, id_cols = c("date", "Period", "What"), names_from = variable) %>% 
  mutate(
    date2 = case_when(
      Period == "2020" ~date,
      Period == "2021" ~as_date(date - years(1))
    )
  )
ts2b <-  ts %>% 
  pivot_wider(values_from = value, id_cols = c("date", "Period", "variable"), names_from = What) %>% 
  mutate(
    date2 = case_when(
      Period == "2020" ~date,
      Period == "2021" ~as_date(date - years(1))
    )
  )




#### ----------------------- ####
#### ---- Extrapolation ---- ####
#### ----------------------- ####

### 0 - Functions --------------------------------------------------------------
source("F:/CAPECET/functions/make_cfact_2.R")
source("F:/CAPECET/functions/rescale2.R")
library(assertthat)
library(WhatIf)
'%!in%' <- function(x,y)!('%in%'(x,y))
'%not_in%' <- purrr::negate(`%in%`)


### 1 - Quantify daily extrapolation -------------------------------------------
extrap <- data.frame(Date = dates, Extrapolation = NA, InformativeData = NA)

## quantify the proportion of extrapolation in daily prediction grids
extrap$Extrapolation <- unlist(
  pblapply(1:length(dates), function(i){
    newdata <- read.table(paste(sourcedir, "/Variables/PredictionGrids/PredictionGrid_", dates[i], ".txt", sep = ""), h = T, sep = "\t", dec = ".")
    newdata <- newdata[,c("Date", "x", "y", covar, "Bathy")]
    newdata[which(newdata$Bathy > 0),] <- NA # masking layer (on bathy & longitude)
    newdata[which(newdata$Bathy < (-2500)),] <- NA
    newdata[which(newdata$x > -1.2),] <- NA
    newdata <- newdata %>% drop_na(all_of(covar)) # remove NA bc make_cfact_2 removes it
    if(i == 110){
      val <- NA
    } else{
      val <- suppressMessages(
             make_cfact_2(calibration_data = seg,
                          test_data = newdata,
                          var_name = covar)   )

    }
    return(val)
  })
)

## plot
extrap <-  extrap %>% mutate(
  Year = year(Date),
  Period = case_when(
    Date %within% int_1920 ~"2020",
    Date %within% int_2021 ~"2021"
  ),
  date2 = case_when(
    Period == "2020" ~Date,
    Period == "2021" ~as_date(Date - years(1))
  )
)
extrap %>% group_by(Period) %>% summarise(mean = mean(Extrapolation, na.rm = T),
                                            min = min(Extrapolation, na.rm = T),
                                            max = max(Extrapolation, na.rm = T),
                                            sd = sd(Extrapolation, na.rm = T))
# Period       mean   min   max    sd
# 1 2019 - 2020 0.267 0.091 0.666 0.123
# 2 2020 - 2021 0.346 0.137 0.81  0.167
write.table(extrap, paste(outdir, "DailyModel_ExtrapolationTable.csv", sep = "/"), row.names = F)


ggplot(extrap) +
  geom_vline(data = dates_ech2, aes(xintercept = as_date(date2)), col = "grey95") +
  geom_line(data = ts2[which(ts2$What == "Prediction"),], aes(x = date2, y = Mean), size = 0.8) +
  geom_ribbon(data = ts2b[which(ts2b$variable == "Mean"),], fill = "grey65",
              aes(x = date2, y = Prediction, ymin = Prediction - Uncertainty, ymax = Prediction + Uncertainty), alpha = 0.5) +
  geom_line(aes(x = date2, y = Extrapolation), color = "#da4308", size = 0.8) +
  geom_pointrange(data = cds, aes(x = as_date(date2), y = D.hat, ymin = D.lcl, ymax = D.ucl), size = 0.2, color = "#14875e") +
  scale_x_date(date_labels = "%b") + scale_color_colorblind() + coord_cartesian(ylim = c(0,2.5)) +
  theme_bw() + layout + facet_grid(Period~., scales = "free") + 
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8), 
        title = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  labs(title = "Estimated density (green), predicted density (black) and extrapolation (red)")
ggsave(paste(outdir, "DailyModel_PredictedDensity_Extrap2b.png", sep = "/"), width = 6, height = 6, dpi = 450)




### 2 - Map daily extrapolation ------------------------------------------------
## identify extrapolation points in daily prediction grids (map)
extrap_stack <- stack(sapply(
  pblapply(1:length(dates), function(i){
    newdata <- read.table(paste(sourcedir, "/Variables/PredictionGrids/PredictionGrid_", dates[i], ".txt", sep = ""), h = T, sep = "\t", dec = ".")
    newdata <- newdata[,c("Date", "x", "y", covar, "Bathy")]
    newdata[which(newdata$Bathy > 0),] <- NA # masking layer
    newdata[which(newdata$Bathy < (-2500)),] <- NA
    newdata[which(newdata$x > -1.2),] <- NA
    newdata <- newdata %>% drop_na(all_of(covar)) # remove NA bc make_cfact_2 removes it
    if(i == 110){
      ras <- rasterFromXYZ( cbind(newdata[,c("x", "y")], NA) )
      names(ras) <- dates[i]
    }else{
      pts <- suppressMessages(
        make_cfact_2(calibration_data = seg, 
                     test_data = newdata, 
                     var_name = covar,
                     percent = FALSE)   )
      ras <- rasterFromXYZ( cbind(newdata[,c("x", "y")], pts) )
      names(ras) <- dates[i]
    }
    list(ras)
  }), "[[", 1 ))
writeENVI(x = extrap_stack, f = paste(outdir, "/DailyModel_Extrapolation.envi", sep =""))




#### ------------------------- ####
#### ---- GIF PREDICTIONS ---- ####
#### ------------------------- ####
Sys.setlocale("LC_TIME", "English")

loop.animate <- function(prediction, uncertainty, ts, extrapolation) {
  ts_pred_1920 <- ts %>% subset(date %within% int_1920 & variable == "Mean" & What == "Prediction")
  ts_pred_1920 <- ts_pred_1920[which(ts_pred_1920$date != "2020-02-29"),]
  se_pred_1920 <- ts %>% subset(date %within% int_1920 & variable == "Mean" & What == "Uncertainty")
  se_pred_1920 <- se_pred_1920[which(se_pred_1920$date != "2020-02-29"),]
  
  ts_pred_2021 <- ts %>% subset(date %within% int_2021 & variable == "Mean" & What == "Prediction")
  se_pred_2021 <- ts %>% subset(date %within% int_2021 & variable == "Mean" & What == "Uncertainty")
  
  extrap_1920 <- extrap %>% subset(Period == "2019 - 2020")
  extrap_2021 <- extrap %>% subset(Period == "2020 - 2021")
  
  pred_1920 <- prediction[[ which(ymd(str_sub(names(prediction), start = -10)) %within% int_1920) ]]
  pred_1920 <- dropLayer(pred_1920, which(ymd(str_sub(names(pred_1920), start = -10))== "2020-02-29"))
  unc_1920 <- uncertainty[[ which(ymd(str_sub(names(uncertainty), start = -10)) %within% int_1920) ]]
  unc_1920 <- dropLayer(unc_1920, which(ymd(str_sub(names(unc_1920), start = -10))== "2020-02-29"))
  
  pred_2021 <- prediction[[ which(ymd(str_sub(names(prediction), start = -10)) %within% int_2021) ]]
  unc_2021 <- uncertainty[[ which(ymd(str_sub(names(uncertainty), start = -10)) %within% int_2021) ]]
  
  p_rmin <- min(cellStats(prediction, "min"))
  p_rmax <- max(cellStats(prediction, "max"))
  
  u_rmin <- min(cellStats(uncertainty, "min"))
  u_rmax <- max(cellStats(uncertainty, "max"))

  pbapply::pblapply(1:nlayers(pred_1920), function(i) { 
    par(mfrow = c(2,3), mar = c(2,2.5,4.2,2), oma = c(2,0,1,1), mgp = c(3,0.8,0))
    dat <- as_date(str_sub(names(pred_1920[[i]]), start = -10))
    # 19-20
    plot(pred_1920[[i]], col = viridis(20, option = "H"), axes = T,
         zlim = c(p_rmin, p_rmax), xlim = c(-4.8, -1),
         addfun = c(plot(land, add = T, col = "grey"), 
                    lines(isob))) ; box()
    rect(xleft = -4.7, xright = -3.5, ybottom = 44.4, ytop = 44.6, col = "white", border = NA)
    text(x = -4.1, y = 44.5, font = 2, cex = 1.5,
         as_date(str_sub(names(pred_1920[[i]]), start = -10)))
    mtext("Prediction", side = 3, adj = 0, line = 0.5, cex = 1.5)
    mtext("A - 2020", side = 3, adj = 0, line = 2.2, cex = 1.6, col = "#da4308")
    
    plot(ts_pred_1920$date, ts_pred_1920$value, col = NA,
         type = "o", xlab = "", ylab = "Mean density", 
         xlim = c(min(ts_pred_1920$date), max(ts_pred_1920$date)), 
         ylim = c(0, max( ts_pred_1920$value)))
    points(ts_pred_1920$date[which(ts_pred_1920$date <= dat)],
           ts_pred_1920$value[which(ts_pred_1920$date <= dat)], type = "l")
    points(se_pred_1920$date[which(se_pred_1920$date <= dat)],
           se_pred_1920$value[which(se_pred_1920$date <= dat)], type = "l", lty = 2, col = "black")
    points(extrap_1920$Date[which(extrap_1920$Date <= dat)],
           extrap_1920$Extrapolation[which(extrap_1920$Date <= dat)], type = "l", col = "red")
    legend("topleft", lty = c(1, 2, 1), cex = 1.3,
           legend = c("Prediction", "Uncertainty", "Extrapolation"), col = c("black", "black", "red"))
    
    plot(unc_1920[[i]], col = viridis(20, option = "H"), axes = T,
         zlim = c(u_rmin, u_rmax), xlim = c(-4.8, -1),
         addfun = c(plot(land, add = T, col = "grey"), 
                    lines(isob))) ; box()
    rect(xleft = -4.7, xright = -3.5, ybottom = 44.4, ytop = 44.6, col = "white", border = NA)
    text(x = -4.1, y = 44.5, font = 2, cex = 1.5,
         as_date(str_sub(names(pred_1920[[i]]), start = -10)))
    mtext("Uncertainty", side = 3, adj = 0, line = 0.5, cex = 1.5)
    
    # 20-21
    plot(pred_2021[[i]], col = viridis(20, option = "H"), axes = T,
         zlim = c(p_rmin, p_rmax), xlim = c(-4.8, -1),
         addfun = c(plot(land, add = T, col = "grey"), 
                    lines(isob))) ; box()
    rect(xleft = -4.7, xright = -3.5, ybottom = 44.4, ytop = 44.6, col = "white", border = NA)
    text(x = -4.1, y = 44.5, font = 2, cex = 1.5,
         as_date(str_sub(names(pred_2021[[i]]), start = -10)))
    mtext("Prediction", side = 3, adj = 0, line = 0.5, cex = 1.5)
    mtext("B - 2021", side = 3, adj = 0, line = 2.2, cex = 1.6, col = "#da4308")
    
    plot(ts_pred_2021$date, ts_pred_2021$value, col = NA,
         type = "o", xlab = "", ylab = "Mean density", 
         xlim = c(min(ts_pred_2021$date), max(ts_pred_2021$date)), 
         ylim = c(0, max(ts_pred_2021$value)))
    points(ts_pred_2021$date[which(ts_pred_2021$date <= dat+years(1))],
           ts_pred_2021$value[which(ts_pred_2021$date <= dat+years(1))], type = "l")
    points(se_pred_2021$date[which(se_pred_2021$date <= dat+years(1))],
           se_pred_2021$value[which(se_pred_2021$date <= dat+years(1))], type = "l", lty = 2, col = "black")
    points(extrap_2021$Date[which(extrap_2021$Date <= dat+years(1))],
           extrap_2021$Extrapolation[which(extrap_2021$Date <= dat+years(1))], type = "l", col = "red")
    legend("topleft", lty = c(1, 2, 1), cex = 1.3,
           legend = c("Prediction", "Uncertainty", "Extrapolation"), col = c("black", "black", "red"))
    
    plot(unc_2021[[i]], col = viridis(20, option = "H"), axes = T,
         zlim = c(u_rmin, u_rmax), xlim = c(-4.8, -1),
         addfun = c(plot(land, add = T, col = "grey"), 
                    lines(isob))) ; box()
    rect(xleft = -4.7, xright = -3.5, ybottom = 44.4, ytop = 44.6, col = "white", border = NA)
    text(x = -4.1, y = 44.5, font = 2, cex = 1.5,
         as_date(str_sub(names(pred_2021[[i]]), start = -10)))
    mtext("Uncertainty", side = 3, adj = 0, line = 0.5, cex = 1.5)
  })
}


animation::saveGIF(loop.animate(prediction = pred_stack, uncertainty = se_stack, ts = ts, extrapolation = extrap), 
                   interval = .5, ani.width = 800, ani.height = 600,
                   movie.name = paste(outdir, "/Animation_Prediction_Uncertainty_Winter.gif", sep = "") )

