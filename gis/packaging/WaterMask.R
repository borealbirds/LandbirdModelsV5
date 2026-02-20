# ---
# title: National Models 5.0 - water mask
# author: Elly Knight
# created: Feb 14, 2026
# ---

library(tidyverse)
library(terra)
library(sf)

root <- "G:/Shared drives/BAM_NationalModels5"

#1. Read in the 1 km file for North America ----
raw <- rast(file.path(root, "gis", "GlobalSurfaceWater_Occurrence_1km.tif"))
occur <- raw |> 
  project("EPSG:3978")

#2. Remove occurrence < 1
water <- ifel(occur < 80, NA, 1)

#2. Get the BCRs ----
bcr.country <- read_sf(file.path(root, "gis", "Subregions_unbuffered.shp")) |> 
  st_transform("EPSG:3978") 

#3. Get the mosaic extents ----
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

#3. Get a mosaic of each extent ----
m.can <- rast(file.path(root, "output", "08_mosaics", "Canada", "ALFL", "ALFL_2020.tif"))
m.ak <- rast(file.path(root, "output", "08_mosaics", "Alaska", "ALFL", "ALFL_2020.tif"))
m.48 <- rast(file.path(root, "output", "08_mosaics", "Lower48", "ALFL", "ALFL_2020.tif"))

#4. Resample ----
can <- water |> 
  crop(bcr.can, mask=TRUE) |> 
  resample(m.can)

ak <- water |> 
  crop(bcr.ak, mask=TRUE) |> 
  resample(m.ak)

l48 <- water |> 
  crop(bcr.48, mask=TRUE) |> 
  resample(m.48)
  
#5. Save ----
writeRaster(can, file.path(root, "gis", "WaterMask_Canada.tif"))
writeRaster(ak, file.path(root, "gis", "WaterMask_Alaska.tif"))
writeRaster(l48, file.path(root, "gis", "WaterMask_Lower48.tif"))
