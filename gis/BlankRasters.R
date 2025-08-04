# ---
# title: National Models 5.0 - create blank rasters
# author: Elly Knight
# created: July 21, 2025
# ---

#NOTES################################

# This script creates blank template rasters for each model subunit and the mosaic

#PREAMBLE############################

#1. Load packages----
library(tidyverse)
library(terra)
library(sf)

#2. Set root path----
root <- "G:/Shared drives/BAM_NationalModels5"

#3. Subunit polygons----
bcr.country <- read_sf(file.path(root, "gis", "BAM_BCR_NationalModel_Unbuffered.shp")) |> 
  mutate(bcr= paste0(country, subUnit))

#BLANK RASTERS#######################

#1. Set up loop----
for(i in 1:nrow(bcr.country)){
  
  #2. Load the raster stack----
  rast.i <- rast(file.path(root, "gis", "stacks", paste0(bcr.country$bcr[i], "_2000.tif")))
  
  #3. Get the BCR boundary----
  sf.i <- bcr.country |> 
    dplyr::filter(bcr==bcr.country$bcr[i]) |> 
    st_transform(crs=crs(rast.i)) |> 
    vect()
  
  #4. Crop to BCR boundary----
  blanks.i <- rast.i$method |> 
    crop(sf.i, mask=TRUE)
  
  #5. Set non NAs to zero----
  blanks.i <- ifel(!is.na(blanks.i), 0, NA)
  
  #5. Save----
  terra::writeRaster(blanks.i, file.path(root, "gis", "blanks", paste0(bcr.country$bcr[i], ".tif")), overwrite=TRUE)
  
  cat("Finished raster", i, "of", nrow(bcr.country), "\n")
  
}

#6. Load a mosaic raster----
rast.mosaic <- rast(file.path(root, "output", "08_mosaics", "OSFL_2020.tif"))

#7. Add it to the list----
blanks.m <- ifel(!is.na(rast.mosaic[[1]]), 0, NA)

#9. Save it ----
terra::writeRaster(blanks.m, file.path(root, "gis", "blanks", "mosaic.tif"), overwrite=TRUE)
