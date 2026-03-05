# ---
# title: National Models 5.0 - water mask
# author: Elly Knight
# created: Feb 14, 2026
# ---

library(tidyverse)
library(terra)
library(sf)

root <- "G:/Shared drives/BAM_NationalModels5"

#WATER POLYGON LAYER ####

#1. Read in the raw shapefile ---- 
raw <- read_sf(file.path(root, "Regions", "Lakes_and_Rivers", "lhy_000c16a_e.shp"))

#2. Estimate area ----
water <- raw |> 
  mutate(area = as.numeric(st_area(raw))) |> 
  dplyr::filter(area > 1000000) |> 
   st_transform("EPSG:3978")

write_sf(water, file.path(root, "gis", "WaterMask.shp"))