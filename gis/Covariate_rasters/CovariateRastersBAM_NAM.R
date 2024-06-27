# ---
# title: National Models 4.1 - covariate raster GIS
# author: Elly Knight
# created: October 11, 2023
# ---

#NOTES################################

#This script wrangles covariate rasters prior to buffer extraction for the national models:
#1. merges 175 annual landcover raster stacks for Alaska and saves them as separate wall to wall tifs for each year of data.

#This script uses CONUS Albers equal-area conic (EPSG: 5072)

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(terra) #basic raster handling

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_NationalModels5"

#ABOVE#####

#1. Get list of tifs----
files <- data.frame(path=list.files("G:/Shared drives/BAM_NationalModels/NationalModels4.1/CovariateRasters/ABoVE/rawdata", pattern="*tif", full.names = TRUE)) %>% 
  dplyr::filter(str_sub(path, -22, -13)!="Simplified")

#2. Read them in as a list----
rlist <- list()
for(i in 1:nrow(files)){
  
  rlist[[i]] <- rast(files$path[i])
  
  print(paste0("Read in raster ", i, " of ", nrow(files)))
}

#3. Make a spatraster collection----
rsrc <- sprc(rlist)

#4. Merge----
m <- merge(rsrc)

#5. Separate bands and save----
years <- seq(1984, 2014, 1)
for(i in 1:length(years)){
  
  m.i <- m[[1]]
  
  names(m.i) <- "landcover"
  
  writeRaster(m.i, file.path(root, "CovariateRasters", "ABoVE", "merged", paste0("ABoVE_", years[i], ".tif")), overwrite=TRUE)
  
  print(paste0("Saved raster ", i, " of ", length(years)))
  
}

