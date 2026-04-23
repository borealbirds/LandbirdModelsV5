# ---
# Title: Packaging BCR model reliability layers
# Author: Elly Knight
# Date: April 2026
# ---

#NOTES################################

#PREAMBLE####

#1. Load packages----
library(sf)
library(tidyverse)
library(terra)

#2. Root file path----
root <- "G:/Shared drives/BAM_NationalModels5"

#3. Subunit polygons----
bcr.all <- read_sf(file.path(root, "gis", "Subregions_unbuffered.shp")) |> 
  st_transform("EPSG:3978")

#4. Make a todo list ----
loop <- expand.grid(bcr = bcr.all$bcr, year=unique(extrap$year))

#PACKAGE############

#1. Set up loop ----
for(i in 1:nrow(loop)){
  
  bcr.i <- loop$bcr[i]
  year.i <- loop$year[i]
  
  #2. Get the bcr boundary ----
  bcr <- bcr.all |> 
    dplyr::filter(bcr==bcr.i)
  
  #3. Get the 3 layers ----
  extrap.i <- rast(file.path(root, "gis", "extrapolation", paste0(bcr.i, "_", year.i, ".tif")))
  train.i <- rast(file.path(root, "gis", "samplingdensity", "train", paste0(bcr.i, "_", year.i, ".tif")))
  test.i <- rast(file.path(root, "gis", "samplingdensity", "test", paste0(bcr.i, "_", year.i, ".tif")))
  
  #4. Calculate means ----
  extrap.mean <- mean(extrap.i)
  train.mean <- mean(train.i)
  test.mean <- mean(test.i)
  
  #3. Stack & tidy ----
  out <- c(extrap.mean, train.mean, test.mean) |> 
    project("EPSG:3978") |> 
    crop(bcr, mask=TRUE)
  
  #4. Fix names ----
  names(out) <- c("extrapolation", "samplingdensity_train", "samplingdensity_test")
  
  #5. Make folder as needed ----
  if(!(file.exists(file.path(root, "output", "10_packaged","Sampling", bcr.i)))){
    dir.create(file.path(root, "output", "10_packaged", "Sampling", bcr.i))
  }
  
  #6. Save
  writeRaster(out, file.path(root, "output", "10_packaged", "Sampling", bcr.i, paste0(bcr.i, "_", year.i, ".tif")), overwrite=T)
  
  cat(i, "  ")
  
}