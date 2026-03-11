# ---
# Title: Calculating extrapolation of models
# Author: Anna Drake, adapted by Elly Knight
# Date: June 2025
# ---

#NOTES################################

#PREAMBLE####

#1. Load packages----
print("* Loading packages on master *")
library(sf)
library(tidyverse)
library(terra)

#2. Root file path----
root <- "G:/Shared drives/BAM_NationalModels5"

#3. Get overlap raster----
MosaicOverlap <- rast(file.path(root, "gis", "ModelOverlap.tif"))

#4. List of rasters----
extrap <- data.frame(file = list.files(file.path(root, "gis", "extrapolation"), pattern="*.tif", recursive = TRUE)) |> 
  separate(file, into=c("bcr", "year", "filetype"), remove=FALSE) |> 
  mutate(year = as.numeric(year),
         path = file.path(root, "gis", "extrapolation", file),
         layer = "extrapolation") |> 
  dplyr::select(-file)

train <- data.frame(file = list.files(file.path(root, "gis", "samplingdensity", "train"), pattern="*.tif", recursive = TRUE)) |> 
  separate(file, into=c("bcr", "year", "filetype"), remove=FALSE) |> 
  mutate(year = as.numeric(year),
         path = file.path(root, "gis", "samplingdensity", "train", file),
         layer = "train") |> 
  dplyr::select(-file)

test <- data.frame(file = list.files(file.path(root, "gis", "samplingdensity", "test"), pattern="*.tif", recursive = TRUE)) |> 
  separate(file, into=c("bcr", "year", "filetype"), remove=FALSE) |> 
  mutate(year = as.numeric(year),
         path = file.path(root, "gis", "samplingdensity", "test", file),
         layer = "test") |> 
  dplyr::select(-file)

files <- rbind(extrap, train, test)

#6. Subunit polygons----
bcr.country <- read_sf(file.path(root, "gis", "Subregions_unbuffered.shp"))

#7. Mosaic polygons----
print("* Mosaicing bcrs *")

akbox <- st_as_sfc(st_bbox(c(xmin= -3693930,
                             ymin = 2187083,
                             xmax = -1617275,
                             ymax = 4562389),
                           crs = st_crs(bcr.country)))

bcr.can <- bcr.country |> 
  dplyr::filter(country=="can") |> 
  st_union() |> 
  st_transform("EPSG:3978") |> 
  vect() |> 
  aggregate() |> 
  fillHoles() |> 
  st_as_sf()

bcr.ak <- bcr.country |> 
  dplyr::filter(bcr %in% c("usa41423", "usa2", "usa40", "usa43", "usa5")) |> 
  st_crop(akbox) |> 
  st_union() |> 
  st_transform("EPSG:3978") |> 
  vect() |> 
  aggregate() |> 
  fillHoles() |> 
  st_as_sf()

bcr.48 <- bcr.country |> 
  dplyr::filter(bcr %in% c("usa5", "usa9", "usa10", "usa11", "usa13", "usa14", "usa23", "usa28")) |> 
  st_difference(akbox) |> 
  st_union() |> 
  st_transform("EPSG:3978") |> 
  vect() |> 
  aggregate() |> 
  fillHoles() |> 
  st_as_sf()

#8. Get the water layer----
water <- read_sf(file.path(root, "gis", "WaterMask.shp")) |> 
  st_transform("EPSG:3978") |> 
  vect()

#MOSAIC#######

#1. Set up to loop through years ----
loop <- expand.grid(year = seq(1985, 2020, 5),
                     bcr = c("Canada", "Alaska", "Lower48"))

loop <- expand.grid(year = c(2020),
                     bcr = c("Canada", "Alaska", "Lower48"))

for(i in 1:length(loop)){
  
  year.i <- loop$year[i]
  country.i <- loop$bcr[i]
  
  #2. Filter to the year ----
  files.year <- dplyr::filter(files, year==year.i)
  
  #filter to the right list for country
  if(country.i=="Canada"){
    files.bcr <- files.year |> 
      dplyr::filter(str_sub(bcr, 1, 3)=="can")
    bcr <- bcr.can
  }
  
  if(country.i=="Alaska"){
    files.bcr <- files.year |> 
      dplyr::filter(bcr %in% c("usa41423", "usa2", "usa40", "usa43", "usa5"))
    bcr <- bcr.ak
  }
  
  if(country.i=="Lower48"){
    files.bcr <- files.year |> 
      dplyr::filter(bcr %in% c("usa5", "usa9", "usa10", "usa11", "usa13", "usa14", "usa23", "usa28"))
    bcr <- bcr.48
  }
  
  #3. Get the files for each layer ----
  files.extrap <- files.bcr |> 
    dplyr::filter(layer == "extrapolation") 
  files.train <- files.bcr |> 
    dplyr::filter(layer=="train")
  files.test <- files.bcr |> 
    dplyr::filter(layer=="test")
  
  #4. Import, mosaic, and sum rasters-------
  mosaic.extrap <- lapply(files.extrap$path, rast) |> 
    sprc() |> 
    mosaic(fun="sum")
  mosaic.train <- lapply(files.train$path, rast) |> 
    sprc() |> 
    mosaic(fun="sum")
  mosaic.test <- lapply(files.test$path, rast) |> 
    sprc() |> 
    mosaic(fun="sum")
  
  #5. Fix extent of overlap raster----
  overlap.i <- resample(MosaicOverlap, mosaic.extrap[[1]])

  #6. Divide extrapolation layers by number of overlaps, average -----
  mean.extrap <- mean(mosaic.extrap/overlap.i, na.rm=TRUE)
  mean.train <- mean(mosaic.train/overlap.i, na.rm=TRUE)
  mean.test <- mean(mosaic.test/overlap.i, na.rm=TRUE)
  
  #7. Stack and crop ----
  out <- c(mean.extrap, mean.train, mean.test) |> 
    project("EPSG:3978") |> 
    crop(bcr, mask=TRUE) |> 
    mask(water, inverse=TRUE)
  
  #8. Fix names ----
  names(out) <- c("extrapolation", "samplingdensity_train", "samplingdensity_test")
  
  #8. Make folder as needed ----
  if(!(file.exists(file.path(root, "output", "10_packaged", "SamplingReliability", country.i)))){
    dir.create(file.path(root, "output", "10_packaged", "SamplingReliability", country.i))
  }
  
  #8. Save
  writeRaster(out, file.path(root, "output", "10_packaged", "SamplingReliability", country.i, paste0(country.i, "_", year.i, ".tif")), overwrite=T)
  
  cat(i, "  ")
  
}