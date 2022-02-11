##--------------------------------------------------------------------------------------------------------
## PROJECT : Infra-seasonal variations of common dolphin distribution in the Bay of Biscay
## SCRIPT : Prepare predictions grids for the daily model
##
## Author : Charlotte Lambert
## Last update : Feb 2022
## R version 4.0.4 (2021-02-15) -- "Lost Library Book"
##--------------------------------------------------------------------------------------------------------


library(rgdal)
library(rgeos)
options(stringsAsFactors = FALSE)
library(raster)

outdir <- "~/3_script_results"
sourcedir <- "~/1_data"

## 1 - Useful vectors ----------------------------------------------------------
static_vars <- c("Bathy", 
              "gradChl_front_persistence", "gradChl_avg_dist_to_front",
              "gradNPP_front_persistence", "gradNPP_avg_dist_to_front",
              "gradMLD_front_persistence", "gradMLD_avg_dist_to_front",
              "gradTemp_front_persistence", "gradTemp_avg_dist_to_front")

dyn_vars = c("Chl", "gradChl", "NPP", "gradNPP", "DissIC", "SPCO2", "Phyto", "ZEu", "Temp", "gradTemp", "CurrentSpeed",
         "EKE", "Salinity", "gradSal", "SSH", "gradSSH", "MLD", "gradMLD",
         "gradChl_front_distance", "gradNPP_front_distance", "gradMLD_front_distance", "gradTemp_front_distance")
lags <- c(1,2,4,7,10,30)

dates <- as_date(c(ymd("2020-01-01"):ymd("2020-04-13"), ymd("2021-01-01"):ymd("2021-04-13")))


## 2 - Resample static layers 
raster_ref <- raster(paste(sourcedir, "/Variables/gradTemp_frontal_persistence-2019-12-01_2020-04-13.envi", sep =""))

bathy <- raster(paste(sourcedir, "/Variables/Bathy.tif", sep = ""))
bathy_rspl <- raster::resample(x = bathy, y = raster_ref)


## 3 - Compose prediction tables with all covariates ---------------------------
pbapply::pboptions(type = "timer", style = 4)
pbapply::pblapply(1:length(dates), function(i){
  # static variables
  tab <- cbind(as.data.frame(bathy_rspl, xy = T),
               
               as.data.frame(raster(paste(sourcedir, "/Variables/gradChl_frontal_persistence-", year(dates[i]-1), "-12-01_", year(dates[i]), "-04-13.envi", sep =""))),
               as.data.frame(raster(paste(sourcedir, "/Variables/gradChl_front_distance-averaged-", year(dates[i]-1), "-12-01_", year(dates[i]), "-04-13.envi", sep =""))),
               
               as.data.frame(raster(paste(sourcedir, "/Variables/gradNPP_frontal_persistence-", year(dates[i]-1), "-12-01_", year(dates[i]), "-04-13.envi", sep =""))),
               as.data.frame(raster(paste(sourcedir, "/Variables/gradNPP_front_distance-averaged-", year(dates[i]-1), "-12-01_", year(dates[i]), "-04-13.envi", sep =""))),
               
               as.data.frame(raster(paste(sourcedir, "/Variables/gradMLD_frontal_persistence-", year(dates[i]-1), "-12-01_", year(dates[i]), "-04-13.envi", sep =""))),
               as.data.frame(raster(paste(sourcedir, "/Variables/gradMLD_front_distance-averaged-", year(dates[i]-1), "-12-01_", year(dates[i]), "-04-13.envi", sep =""))),
               
               as.data.frame(raster(paste(sourcedir, "/Variables/gradTemp_frontal_persistence-", year(dates[i]-1), "-12-01_", year(dates[i]), "-04-13.envi", sep =""))),
               as.data.frame(raster(paste(sourcedir, "/Variables/gradTemp_front_distance-averaged-", year(dates[i]-1), "-12-01_", year(dates[i]), "-04-13.envi", sep ="")))
  )
  names(tab) <- c("x", "y", static_vars) # be careful to check the order of names
  if(year(dates[i]) == 2020){ # correct for a typo when creating the front rasters
    tab[, "gradChl_avg_dist_to_front"] <- tab[, "gradChl_avg_dist_to_front"] * 202
    tab[, "gradNPP_avg_dist_to_front"] <- tab[, "gradNPP_avg_dist_to_front"] * 202
    tab[, "gradMLD_avg_dist_to_front"] <- tab[, "gradMLD_avg_dist_to_front"] * 202
    tab[, "gradTemp_avg_dist_to_front"] <- tab[, "gradTemp_avg_dist_to_front"] * 202
  }
  if(year(dates[i]) == 2021){ # correct for a typo when creating the front rasters
    tab[, "gradChl_avg_dist_to_front"] <- tab[, "gradChl_avg_dist_to_front"] * 201
    tab[, "gradNPP_avg_dist_to_front"] <- tab[, "gradNPP_avg_dist_to_front"] * 201
    tab[, "gradMLD_avg_dist_to_front"] <- tab[, "gradMLD_avg_dist_to_front"] * 201
    tab[, "gradTemp_avg_dist_to_front"] <- tab[, "gradTemp_avg_dist_to_front"] * 201
  }
 
  # dynamic variables
  for(j in 1:length(dyn_vars)){
    # open source file
    ras <- stack(paste(sourcedir, "/Variables/", dyn_vars[j], "-daily-", year(dates[i]-1), "-12-01_", year(dates[i]), "-04-13.envi", sep =""))
    # isolate target layer
    ras_d <- ras[[ which(ymd(str_sub(names(ras), start = -10)) == dates[i]) ]]
    names(ras_d) <- dyn_vars[j]
    # add dayly value (0d lag)
    tab2 <- as.data.frame(ras_d, xy = T)
    # add other lags
    for(k in 1:length(lags)){
      ras_tres <- ras[[ which(ymd(str_sub(names(ras), start = -10)) == dates[i]) - lags[k] ]]
      names(ras_tres) <- paste(dyn_vars[j], "_", lags[k], "d", sep = "")
      suppressMessages( tab2 <- full_join(tab2, as.data.frame(ras_tres, xy = T)) )
    }
    suppressMessages( tab <- full_join(tab, tab2) )
  }
  
  # export as table
  tab <- cbind(Date = Date, tab)
  write.table(tab, paste(sourcedir, "/PredictionGrids/PredictionGrid_", dates[i], ".txt", sep = ""), sep = "\t", dec = ".", row.names = F)

  return(tab)
})
beepr::beep("fanfare")


