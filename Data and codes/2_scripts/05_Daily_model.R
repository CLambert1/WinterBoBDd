##--------------------------------------------------------------------------------------------------------
## PROJECT : Infra-seasonal variations of common dolphin distribution in the Bay of Biscay
## SCRIPT : Adjust the daily model 
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
outdir <- "~/3_script_results"
sourcedir <- "~/1_data"
funcdir <- "~/0_functions"
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



### 3 - open land --------------------------------------------------------------------
land <- st_read(dsn = paste(sourcedir, "Map_files", sep = "/"), layer = "Land")
bbx <- as(raster::extent(-4.8,0.5,44.3,47.1), "SpatialPolygons")
crs(bbx) <- crs(land)
land <- gIntersection(as_Spatial(land), bbx, byid = TRUE)



#### ------------------------- ####
#### ---- Fit daily model ---- ####
#### ------------------------- ####

### 0 - Functions --------------------------------------------------------------
source("F:/CAPECET/functions/fctGAM_20210316.R")
correlations <- function(predictors, effort) {
  corr <- Hmisc::rcorr(as.matrix(effort[!is.na(effort[,predictors]), predictors]), type = "spearman")
  corrplot::corrplot(corr$r, type="upper", order="hclust", tl.col="black", tl.srt=45, method = "number", diag = F,
                     addCoefasPercent = TRUE, tl.cex = 0.7, is.corr = TRUE, number.cex = .5)
}


### 1 - Model selection with variable resolutions pre-selected from script 04 ----
## select, for each variable, resolutions with DeltaAIC < 2 & max AICw <=0.5
var_to_test <- var_name[var_name %in% c("gradChl_front_persistence", "gradChl_avg_dist_to_front",
                                        "gradNPP_front_persistence", "gradNPP_avg_dist_to_front",
                                        "gradMLD_front_persistence", "gradMLD_avg_dist_to_front",
                                        "gradTemp_front_persistence",
                                        "Chl_7d",
                                        "NPP_10d",
                                        "DissIC_10d",
                                        "SPCO2_10d",
                                        "Phyto_10d", 
                                        "ZEu_7d",
                                        "Temp_10d", 
                                        "EKE_10d","EKE_4d",
                                        "Salinity_10d",
                                        "gradSal_10d",
                                        "gradSSH_10d",
                                        "MLD_10d",
                                        "gradMLD_10d", "gradMLD_30d", "gradMLD_2d", 
                                        "gradChl_front_distance_2d",
                                        "gradNPP_front_distance_4d",
                                        "gradMLD_front_distance_4d", "gradMLD_front_distance_30d",
                                        "gradTemp_front_distance_30d"
)]


## correlations
correlations(effort = seg, predictors = var_to_test)


## fit all 
fit_best <- fit_all_gam(envdata = seg, esw = seg$ESW, outcome = "n_ind",
                        predictors = var_to_test, 
                        rescale = F, weight_by_g_naught = 1, family = "tweedie", 
                        max_cor = 0.6, nb_min_pred = 1, nb_max_pred = 4, complexity = 4)

# compute AIC weights and likelihoods
fit_best_ord <- fit_best[order(fit_best$AIC),]
aics <- akaike.weights(fit_best_ord[,3])

fit_best_ord$Delta.AIC <- aics$deltaAIC  # < 2: probable model; 4 - 7: less probable; > 10: not probable

fit_best_ord$rel.likelihood <- aics$rel.LL # relative likelihood : exp((AICmin-AICi)/2)

fit_best_ord$Akaike.weight <- aics$weights # AICw = the probability a model is the best among all tested models


## save
write.table(fit_best_ord, paste(outdir, "/Summary_model_selection.csv", sep = ""), sep = ";", dec = ".")

save.image(paste(outdir, "/2winters_JanMar_Modelling_weeks.RData", sep = ""))



### 2 - Average AIC weight of variables ----
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
write.table(table, paste(outdir, "/Summary_model_selection_average_AICweight_variables.csv", sep = ""), sep = ";", dec = ".")



### 3 - Fit and check best models ----
## Fit
# first model
mod1 <- gam(n_ind ~ 1 +s(gradMLD_avg_dist_to_front, k = 4, bs = 'tp') + 
              s(gradTemp_front_persistence, k = 4, bs = 'tp') + 
              s(Chl_7d, k = 4, bs = 'tp') + 
              s(gradChl_front_distance_2d, k = 4, bs = 'tp') 
            , data = seg, offset = log(seg$Effort*2*seg$ESW), family = tw())
plot(mod1, pages = 1, scale = 0, res = T)
summary(mod1)
mod1$aic

# try first model with lon, lat
mod1_b <- gam(n_ind ~ 1 +s(gradMLD_avg_dist_to_front, k = 4, bs = 'tp') + 
                s(gradTemp_front_persistence, k = 4, bs = 'tp') + 
                s(Chl_7d, k = 4, bs = 'tp') + 
                s(gradChl_front_distance_2d, k = 4, bs = 'tp') 
              + te(lon, lat, bs = "tp", k = 4)
              , data = seg, offset = log(seg$Effort*2*seg$ESW), family = tw())
plot(mod1_b, pages = 1, scale = 0, res = T)
summary(mod1_b)
mod1_b$aic


# second model
mod2 <- gam(n_ind ~ 1 +s(gradMLD_avg_dist_to_front, k = 4, bs = 'tp') + 
              s(Chl_7d, k = 4, bs = 'tp') + 
              s(MLD_10d, k = 4, bs = 'tp') + 
              s(gradChl_front_distance_2d, k = 4, bs = 'tp') 
            , data = seg, offset = log(seg$Effort*2*seg$ESW), family = tw())
plot(mod2, pages = 1, scale = 0, res = T)
summary(mod2)
mod2$aic

mod2_b <- gam(n_ind ~ 1 +s(gradMLD_avg_dist_to_front, k = 4, bs = 'tp') + 
                s(Chl_7d, k = 4, bs = 'tp') + 
                s(MLD_10d, k = 4, bs = 'tp') + 
                s(gradChl_front_distance_2d, k = 4, bs = 'tp')   
              + te(lon, lat, bs = "tp", k = 4)
              , data = seg, offset = log(seg$Effort*2*seg$ESW), family = tw())
plot(mod2_b, pages = 1, scale = 0, res = T)
summary(mod2_b)
mod2_b$aic



### 3 - Plot best model  ----
## plot summed effects rather than partial effects (based on predictions)
png(paste(outdir, "/Best_model_smoothplot.png", sep = ""), height = 5, width = 20, res = 400, units = "cm")
par(mfrow = c(1,4), mar = c(2,2,2,1), oma = c(0,1.5,0.8,0.5), mgp = c(2.5,0.5,0), tck = -0.04)
plot_smooth(mod2, view = "Chl_7d", sim.ci = T, set.seed(1e4), ylab = "", main = "", hide.label = T, h0 = 1, transform = exp) ; box()
mtext(expression(ind.km^2), side = 2, line = 2, cex = 0.8)
mtext("Chl. concentration\n(7d lag)", side = 3, line = 0.3, cex = 0.8)
plot_smooth(mod2, view = "MLD_10d", sim.ci = T, set.seed(1e4), ylab = "", main = "", hide.label = T, h0 = 1, transform = exp, ylim = c(0,4)) ; box()
mtext("MLD (10d lag)", side = 3, line = 0.7, cex = 0.8)
plot_smooth(mod2, view = "gradMLD_avg_dist_to_front", sim.ci = T, set.seed(1e4), ylab = "", main = "", hide.label = T, h0 = 1, transform = exp) ; box()
mtext("Averaged distance\nto MLD front", side = 3, line = 0.3, cex = 0.8)
plot_smooth(mod2, view = "gradChl_front_distance_2d", sim.ci = T, set.seed(1e4), ylab = "", main = "", hide.label = T, h0 = 1, transform = exp) ; box()
mtext("Dist. to Chl fronts (2d lag)", side = 3, line = 0.7, cex = 0.8)
dev.off()


