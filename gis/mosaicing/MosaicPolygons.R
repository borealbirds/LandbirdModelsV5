# ---
# Title: Mosaic polygons
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

#3. Get region file ----
bcr.all <- read_sf(file.path(root, "gis", "Subregions_unbuffered.shp"))

#MAKE POLYGONS############

#1. Make an alaska boundary box ----
akbox <- st_as_sfc(st_bbox(c(xmin= -3693930,
                             ymin = 2187083,
                             xmax = -1617275,
                             ymax = 4562389),
                           crs = st_crs(bcr.all)))

#2. Make the mosaic regions ----
bcr.can <- bcr.all |> 
  dplyr::filter(country=="can") |> 
  st_union()  |> 
  st_transform("EPSG:3978") |> 
  vect() |> 
  aggregate() |> 
  fillHoles() |> 
  st_as_sf() |> 
  mutate(bcr="Canada",
         subUnit=NA,
         country="can")

bcr.ak <- bcr.all |> 
  dplyr::filter(bcr %in% c("usa41423", "usa2", "usa40", "usa43", "usa5")) |> 
  st_crop(akbox) |> 
  st_union()  |> 
  st_transform("EPSG:3978") |> 
  vect() |> 
  aggregate() |> 
  fillHoles() |> 
  st_as_sf() |> 
  mutate(bcr="Alaska",
         subUnit=NA,
         country="usa")

bcr.48 <- bcr.all |> 
  dplyr::filter(bcr %in% c("usa5", "usa9", "usa10", "usa11", "usa12", "usa13", "usa14", "usa23", "usa28")) |> 
  st_difference(akbox) |> 
  st_union() |> 
  st_transform("EPSG:3978") |> 
  vect() |> 
  aggregate() |> 
  fillHoles() |> 
  st_as_sf() |> 
  mutate(bcr="Lower48",
         subUnit=NA,
         country="usa")

#3. Put back together ----
bcr.out <- bcr.all |> 
  st_transform("EPSG:3978") |> 
  rbind(bcr.can, bcr.ak, bcr.48)

#4. Save ----
write_sf(bcr.out, file.path(root, "gis", "Subregions_Mosaics_EPSG3978.shp"))
