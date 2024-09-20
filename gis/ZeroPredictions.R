#-----------------------------------
#  Create zero predictions for raster mosaicking
#  Author: Elly Knight
#  Created: September 2024
#-----------------------------------

#NOTES###########

#This code produces a layer of zero bird density for each BCR modelling area that is used to supplement the actual model predictions in `09.MosaicPredictions.R` for subunits that aren't modelled for a particular species

#PREAMBLE##########

#1. Load packages ----
library(tidyverse)
library(terra)
library(sf)

#2. NAD83(NSRS2007)/Conus Albers projection (epsg:5072) ----
crs <- "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"

#3. Set file paths ----
root <- "G:/Shared drives/BAM_NationalModels5"

#4. Read in buffered and merged region shapefile----
bcr <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel_Buffered.shp"))

#ZERO PREDICTIONS#####

#1. Get the list of stacks for one year----
stacks <- data.frame(file=list.files(file.path(root, "gis", "stacks"), pattern="*.tif"),
                     path=list.files(file.path(root, "gis", "stacks"), pattern="*.tif", full.names = TRUE)) |> 
  separate(file, into=c("bcr", "year", "tif"), remove=FALSE) |>
  mutate(year = as.numeric(year)) |> 
  dplyr::filter(str_sub(file, -3, -1)!="xml",
                year==2020)

#2. Set up loop----
for(i in 1:nrow(stacks)){
  
  #3. Read in the stack and take just the first one-----
  stack.i <- rast(stacks$path[i])[[1]]
  
  #4. Convert numbers to 0s-----
  zero.i <- classify(stack.i, cbind(NA, NA, NA))
  zero.i[!is.na(zero.i)] <- 0
  
  #5. Save----
  terra::writeRaster(zero.i, file.path(root, "gis", "zeros", paste0(bcr, ".tif")), overwrite=TRUE)
  
  cat("Finished", i, "\n")
  
}
