# ---
# title: National Models 5.0 - water mask
# author: Elly Knight
# created: Feb 14, 2026
# ---

library(tidyverse)
library(terra)
library(sf)

root <- "G:/Shared drives/BAM_NationalModels5"

#CANADA WATER POLYGON LAYER ####

#1. Read in the raw shapefile ---- 
raw_ca <- read_sf(file.path(root, "Regions", "WaterMask", "Lakes_and_Rivers", "lhy_000c16a_e.shp"))

#2. Estimate area ----
water_ca <- raw_ca |> 
  mutate(area = as.numeric(st_area(raw))) |> 
  dplyr::filter(area > 1000000) |> 
  st_transform("EPSG:3978")

write_sf(water_ca, file.path(root, "gis", "WaterMask_Canada.shp"))

#US WATER POLYGON LAYER #########

#1. Read the raw AK shapefile ----
path_ak <- file.path(root, "Regions", "WaterMask", "3dhp_all_Alaska_20260112_GDB", "3dhp_all_Alaska_20260112_GDB.gdb")
st_layers(path_ak)
raw_ak <- read_sf(path_ak, "hydro_3dhp_all_waterbody")

#2. Filter by area ----
use_ak <- raw_ak |> 
  dplyr::filter(areasqkm >= 1, featuretypelabel %in% c("Lake", "River")) |> 
  st_transform("EPSG:3978")
rm(raw_ak)

#3. Read the raw conus shapefile ----
path_us <- file.path(root, "Regions", "WaterMask", "3dhp_all_CONUS_20260112_GDB", "3dhp_all_CONUS_20260112_GDB.gdb")
st_layers(path_us)
raw_us <- read_sf(path_us, "hydro_3dhp_all_waterbody")

#4. Filter by area ----
use_us <- raw_us |> 
  dplyr::filter(areasqkm >= 1, featuretypelabel %in% c("Lake", "River")) |> 
  st_transform("EPSG:3978")
rm(raw_us)

#5. Put together ----
water_us <- rbind(use_ak, use_us) |> 
  rename(id = id3dhp, name = gnisidlabel, type = featuretypelabel) |> 
  dplyr::select(id, name, type, areasqkm) |> 
  st_zm(drop=TRUE, what = "ZM")

write_sf(water_us, file.path(root, "gis", "WaterMask_US.shp"))
