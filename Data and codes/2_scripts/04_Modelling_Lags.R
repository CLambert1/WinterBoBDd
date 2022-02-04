##--------------------------------------------------------------------------------------------------------
## PROJECT : Infra-seasonal variations of common dolphin distribution in the Bay of Biscay
## SCRIPT : Modelling the lags between the dolphin distributions and the env. conditions
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
library(qpcR)
library(itsadug)


#### ----------------------- ####
#### ---- SOURCES FILES ---- ####
#### ----------------------- ####
outdir <- "F:/CAPECET/paper/data/3_script_results"
sourcedir <- "F:/CAPECET/paper/data/1_data"
funcdir <- "F:/CAPECET/paper/data/0_functions"
setwd(outdir)

### 1 - Open campaign files ----------------------------------------------------
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


### 2 - Vectors of variables ---------------------------------------------------
vars <- c("Chl", "gradChl", "NPP", "gradNPP", "DissIC", "SPCO2", "Phyto", "ZEu", "Temp", "gradTemp", "CurrentSpeed",
          "EKE", "Salinity", "gradSal", "SSH", "gradSSH", "MLD", "gradMLD",
          "gradChl_front_distance", "gradNPP_front_distance", "gradMLD_front_distance", "gradTemp_front_distance")
resolutions <- c(1, 2, 4, 7, 10, 30)
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



#### ------------------------------- ####
#### ---- Variable correlations ---- ####
#### ------------------------------- ####

### 0 - Function ---------------------------------------------------------------
correlations <- function(predictors, effort) {
  corr <- Hmisc::rcorr(as.matrix(effort[!is.na(effort[,predictors]), predictors]), type = "spearman")
  # symnum(corr, abbr.colnames=FALSE)
  corrplot::corrplot(corr$r, type="upper", order="hclust", tl.col="black", tl.srt=45, method = "number", diag = F,
                     addCoefasPercent = TRUE, tl.cex = 0.7, is.corr = TRUE, number.cex = .5)
}

### 1 - Between variable categories, by temporal resolution --------------------
## daily 
correlations(effort = seg, 
             predictors = var_name[c(str_which(var_name, "_d"),
                                     which(var_name %in% c("Bathy",
                                                           "gradChl_front_persistence", "gradChl_avg_dist_to_front",
                                                           "gradNPP_front_persistence", "gradNPP_avg_dist_to_front",
                                                           "gradMLD_front_persistence", "gradMLD_avg_dist_to_front",
                                                           "gradTemp_front_persistence", "gradTemp_avg_dist_to_front")))])

## others
for(i in 1:length(resolutions)){
  correlations(effort = seg, 
               predictors = var_name[c(str_which(var_name, paste("_", resolutions[i], "d", sep = "")),
                                       which(var_name %in% c("Bathy",
                                                             "gradChl_front_persistence", "gradChl_avg_dist_to_front",
                                                             "gradNPP_front_persistence", "gradNPP_avg_dist_to_front",
                                                             "gradMLD_front_persistence", "gradMLD_avg_dist_to_front",
                                                             "gradTemp_front_persistence", "gradTemp_avg_dist_to_front")))])
}


### 2 - Between temporal resolution, by variable categories --------------------
for(i in 1:length(vars)){
  correlations(effort = seg, 
               predictors = var_name[c(str_which(var_name, vars[i]))])
}




#### ------------------------------------------- ####
#### ---- Testing temporal lag of variables ---- ####
#### ------------------------------------------- ####


### 0 - Functions --------------------------------------------------------------
source(paste("F:/CAPECET/functions/fctGAM_20210316.R"))


### 1 - testing resolutions for each variables ---------------------------------
dyn_var <- var_name[10:163]

## apply fit_all_gam to each dynamic variable, to see which resolution performs best
fit_single <- NULL
for(i in 1:length(vars)){
  # select predictors to test
  pred_to_test <- dyn_var[c(str_which(dyn_var, 
                                      pattern = paste("^", vars[i],"(?!_front_distance)", sep ="")
  ) )] 
  # fit models
  dd <- fit_all_gam(envdata = seg, esw = seg$ESW, outcome = "n_ind",
                    predictors = pred_to_test,
                    rescale = F, weight_by_g_naught = 1, family = "tweedie", 
                    max_cor = 0.7, nb_max_pred = 1, complexity = 4)
  dd.ord <- dd[order(dd$AIC),]
  aics <- akaike.weights(dd.ord[,3])
  dd.ord$Delta.AIC <- aics$deltaAIC
  dd.ord$rel.likelihood <- aics$rel.LL
  dd.ord$Akaike.weight <- aics$weights
  fit_single[[vars[i]]] <- cbind(variable = vars[i], dd.ord)
}

## same with to static var
for(i in 1:9){
  dd <- fit_all_gam(envdata = seg, esw = seg$ESW, outcome = "n_ind",
                    predictors = var_name[i], 
                    rescale = F, weight_by_g_naught = 1, family = "tweedie", 
                    max_cor = 0.7, nb_max_pred = 1, complexity = 4)
  dd.ord <- dd[order(dd$AIC),]
  aics <- akaike.weights(dd.ord[,3])
  
  dd.ord$Delta.AIC <- aics$deltaAIC  
  dd.ord$rel.likelihood <- aics$rel.LL
  dd.ord$Akaike.weight <- aics$weights
  fit_single[[length(vars) + i]] <- cbind(variable = var_name[i], dd.ord)
}

## transforme table
fit_single2 <- do.call("rbind", fit_single)

# retrieve mod number
fit_single2$mod_nbr <- str_sub(rownames(fit_single2), start = -1)
fit_single2$mod_nbr[177:194] <- rep(c(1,2),18)

# retrieve resolution (= lag)
fit_single2$resolution <- NA
for(i in 1:nrow(fit_single2)){
  strg <- unlist( str_extract_all(fit_single2$model[ i ], "[0-9]+") )
  x <- NULL
  if(fit_single2$variable[i] != "SPCO2"){
    if(length(strg) == 3){ x = strg[2] } # model with lagged variables (type "Chl_1d")
    if(length(strg) == 2){ x = 0 }       # model with daily variables (type "Chl")
    if(length(strg) == 1){ x = NA }      # model with intercept only
  } else {
    if(length(strg) == 4){ x = strg[3] }
    if(length(strg) == 3){ x = 0 }      
    if(length(strg) == 1){ x = NA }     
  }
  fit_single2$resolution[i] <- x
}
fit_single2$resolution <- factor(fit_single2$resolution, levels = c("NA", "0", "1", "2", "4", "7", "10", "30"))

# type of variable
fit_single2$type <- ifelse(fit_single2$variable %in% vars, "Dynamic", 
                           ifelse(fit_single2$variable %in% c("Bathy"), "Static", "Seasonal"))



### 2 - plot results -----------------------------------------------------------
# Akaike weight
g1 <- ggplot(fit_single2, aes(x = resolution, y = variable)) + theme_bw() +
  scale_x_discrete(name = "Lag", labels = c("Intercept", "None", "1d", "2d", "4d", "7d", "10d", "30d"),
                   limits = c(NA, "0", "1", "2", "4", "7", "10", "30"), position = "top") +
  scale_y_discrete(name = "Environmental variable",
                   limits = rev(
                     c("Temp", "Salinity", "CurrentSpeed", "Chl", "NPP", "Phyto",
                       "EKE", "gradSSH", "gradTemp", "gradMLD", "gradSal",
                       "SSH", "gradChl", "gradNPP", "ZEu", "MLD", "DissIC","SPCO2",
                       "gradMLD_front_distance", "gradTemp_front_distance", "gradChl_front_distance", "gradNPP_front_distance",
                       "Bathy",
                       "gradMLD_avg_dist_to_front", "gradTemp_avg_dist_to_front", "gradChl_avg_dist_to_front", "gradNPP_avg_dist_to_front",
                       "gradMLD_front_persistence", "gradTemp_front_persistence", "gradChl_front_persistence", "gradNPP_front_persistence"
                     )),
                   labels = rev(
                     c("Temp", "Salinity", "CurrentSpeed", "Chl", "NPP", "Phyto",
                       "EKE", "gradSSH", "gradTemp", "gradMLD", "gradSal",
                       "SSH", "gradChl", "gradNPP", "ZEu", "MLD", "DissIC","SPCO2",
                       "Distance MLD front", "Distance Temp front", "Distance Chl front", "Distance NPP front",
                       "Bathy",
                       "Average dist MLD front", "Average dist Temp front", "Average dist Chl front", "Average dist NPP front",
                       "Persistence MLD front", "Persistence Temp front", "Persistence Chl front", "Persistence NPP front"
                     )))
g1 + geom_point(aes(size = Akaike.weight, color = round(Delta.AIC,1))) + 
  scale_color_gradientn(colors = rev(brewer.pal(5, "YlOrRd")), values = c(0, exp(seq(-5, 0, length.out = 8)))) +
  labs(color = "Delta.AIC") 
ggsave(paste(outdir, "/Single_variable_models_AICWeight.png", sep = ""), dpi = 350, height = 8, width = 7)



### 3 - Plot gam curves for each lag  ------------------------------------------
var <- c("Temp", "Salinity", "CurrentSpeed", "Chl", "NPP", "Phyto", 
         "EKE", "gradSSH", "gradTemp", "gradMLD", "gradSal", 
         "SSH", "gradChl", "gradNPP", "ZEu", "MLD", "DissIC", "SPCO2",
         "gradMLD_front_distance", "gradTemp_front_distance", "gradChl_front_distance", "gradNPP_front_distance", "Bathy",
         "gradMLD_avg_dist_to_front", "gradTemp_avg_dist_to_front", "gradChl_avg_dist_to_front", "gradNPP_avg_dist_to_front",
         "gradMLD_front_persistence", "gradTemp_front_persistence", "gradChl_front_persistence", "gradNPP_front_persistence"
)
# var instead of vars just to reorder levels for plotting purpose

png(paste(outdir, "Single_var_modelsA.png", sep = "/"), res = 400, height = 25, width = 20, units = "cm")
par(mfrow = c(11, 7), mar = c(2,2,2,.5), oma = c(0,0,3,0), cex = 0.5, mgp = c(2, 0.5, 0))
# plot curves
suppressMessages(
  lapply(1:11, function(i){
    var_plot <- dyn_var[ c(str_which(dyn_var, pattern = paste("^", var[i],"(?!_front_distance)", sep =""))) ]
    
    lapply(1:length(var_plot), function(i){
      form <- paste("n_ind ~ s(", var_plot[i], ", k = 4, bs = 'tp')", sep = "")
      mod <- gam(as.formula(form), data = seg, offset = log(seg$Effort*2*seg$ESW), family = tw() )
      
      plot_smooth(mod, view = mod$smooth[[1]]$term, sim.ci = T, set.seed(1e4), ylab = "log(individuals)", main = "",
                  print.summary = F, hide.label = T) ; box() 
      if(i == 1) {  mtext(var_plot[1], side = 3, font = 3, line = .5, cex = 0.7, adj = 0, col = "#da4308") }
    })
  }) )
# add outer labels - lags
mtext("No lag", side = 3, outer = TRUE, at = 0.08, cex = 0.8)
mtext("1d lag", side = 3, outer = TRUE, at = 0.218, cex = 0.8)
mtext("2d lag", side = 3, outer = TRUE, at = 0.362, cex = 0.8)
mtext("4d lag", side = 3, outer = TRUE, at = 0.505, cex = 0.8)
mtext("7d lag", side = 3, outer = TRUE, at = 0.648, cex = 0.8)
mtext("10d lag", side = 3, outer = TRUE, at = 0.79, cex = 0.8)
mtext("30d lag", side = 3, outer = TRUE, at = 0.938, cex = 0.8)
dev.off()


png(paste(outdir, "Single_var_modelsB.png", sep = "/"), res = 400, height = 25, width = 20, units = "cm")
par(mfrow = c(11, 7), mar = c(2,2,2,.5), oma = c(0,0,3,0), cex = 0.5, mgp = c(2, 0.5, 0))
# plot curves
suppressMessages(
  lapply(12:22, function(i){
    var_plot <- dyn_var[ c(str_which(dyn_var, pattern = paste("^", var[i],"(?!_front_distance)", sep =""))) ]
    
    lapply(1:length(var_plot), function(i){
      form <- paste("n_ind ~ s(", var_plot[i], ", k = 4, bs = 'tp')", sep = "")
      mod <- gam(as.formula(form), data = seg, offset = log(seg$Effort*2*seg$ESW), family = tw() )
      
      plot_smooth(mod, view = mod$smooth[[1]]$term, sim.ci = T, set.seed(1e4), ylab = "log(individuals)", main = "",
                  print.summary = F, hide.label = T) ; box() 
      if(i == 1) {  mtext(var_plot[1], side = 3, font = 3, line = .5, cex = 0.7, adj = 0, col = "#da4308") }
    })
  }) )
# add outer labels - lags
mtext("No lag", side = 3, outer = TRUE, at = 0.08, cex = 0.8)
mtext("1d lag", side = 3, outer = TRUE, at = 0.218, cex = 0.8)
mtext("2d lag", side = 3, outer = TRUE, at = 0.362, cex = 0.8)
mtext("4d lag", side = 3, outer = TRUE, at = 0.505, cex = 0.8)
mtext("7d lag", side = 3, outer = TRUE, at = 0.648, cex = 0.8)
mtext("10d lag", side = 3, outer = TRUE, at = 0.79, cex = 0.8)
mtext("30d lag", side = 3, outer = TRUE, at = 0.938, cex = 0.8)
dev.off()


png(paste(outdir, "Single_var_modelsC.png", sep = "/"), res = 400, height = 12, width = 15, units = "cm")
par(mfrow = c(3, 3), mar = c(2,2,2,.5), oma = c(0,0,3,0), cex = 0.5, mgp = c(2, 0.5, 0))
# plot curves
suppressMessages(
  lapply(23:31, function(i){
    var_plot <- var[ i ]
    
    lapply(1:length(var_plot), function(i){
      form <- paste("n_ind ~ s(", var_plot[i], ", k = 4, bs = 'tp')", sep = "")
      mod <- gam(as.formula(form), data = seg, offset = log(seg$Effort*2*seg$ESW), family = tw() )
      
      plot_smooth(mod, view = mod$smooth[[1]]$term, sim.ci = T, set.seed(1e4), ylab = "log(individuals)", main = "",
                  print.summary = F, hide.label = T) ; box() 
      if(i == 1) {  mtext(var_plot[1], side = 3, font = 3, line = .5, cex = 0.7, adj = 0, col = "#da4308") }
    })
  }) )
# add outer labels - lags
mtext("No lag", side = 3, outer = TRUE, at = 0.01, adj = 0, cex = 0.8)
dev.off()


