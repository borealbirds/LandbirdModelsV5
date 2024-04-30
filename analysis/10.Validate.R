# ---
# title: National Models 5.0 - validate predictions
# author: Elly Knight
# created: April 29, 2024
# ---

#NOTES################################

# This script uses the withheld data from each bootstrap to validate the spatial predictions of the models at two levels:
# 1. Study area: uses the tiff output from '10.MosaicPredictions.R' 
# 2. BCR: uses the tiff output from '08.Predict.R'

#PREAMBLE############################

#1. Load packages----
library(tidyverse)
library(terra)
library(sf)

#2. Set root path----
root <- "G:/Shared drives/BAM_NationalModels/NationalModels5.0"

#3. Load data file----

#SUBUNIT POLYGONS#####################

#1. Read in country shapefiles----
can <- read_sf(file.path(root, "Regions", "CAN_adm", "CAN_adm0.shp")) |> 
  st_transform(crs=5072) 
usa <- read_sf(file.path(root, "Regions", "USA_adm", "USA_adm0.shp")) |> 
  st_transform(crs=5072) 

#2. Read in BCR shapefile----
#Remove subunit 1
bcr <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel.shp")) |> 
  dplyr::filter(subUnit!=1) |> 
  st_transform(crs=5072)

ggplot(bcr) +
  geom_sf(aes(fill=factor(subUnit)))

#3. Identify BCRs for each country----
bcr.ca <- bcr |> 
  st_intersection(can) |> 
  mutate(country="can") |> 
  dplyr::select(subUnit, country)

bcr.usa <- bcr |> 
  st_intersection(usa) |> 
  mutate(country="usa") |> 
  dplyr::select(subUnit, country)

#4. Make merged subunits----
bcr.can4142 <- bcr.ca |> 
  dplyr::filter(subUnit %in% c(41, 42)) |> 
  st_union() |> 
  nngeo::st_remove_holes() |> 
  st_sf() |> 
  mutate(country="can",
         subUnit = 4142)
st_geometry(bcr.can4142) <- "geometry"

bcr.usa41423 <- bcr.ca |> 
  dplyr::filter(subUnit %in% c(41, 42, 3)) |> 
  st_union() |> 
  nngeo::st_remove_holes() |> 
  st_sf() |> 
  mutate(country="usa",
         subUnit = 41423)
st_geometry(bcr.usa41423) <- "geometry"

bcr.usa414232 <- bcr.usa |> 
  dplyr::filter(subUnit %in% c(41, 42, 3, 2)) |> 
  st_union() |> 
  nngeo::st_remove_holes() |> 
  st_sf() |> 
  mutate(country="usa",
         subUnit = 414232)
st_geometry(bcr.usa414232) <- "geometry"

bcr.can8182 <- bcr.ca |> 
  dplyr::filter(subUnit %in% c(81, 82)) |> 
  st_union() |> 
  nngeo::st_remove_holes() |> 
  st_sf() |> 
  mutate(country="can",
         subUnit = 8182)
st_geometry(bcr.can8182) <- "geometry"

#polygons to remove
bcr.remove <- data.frame(country=c("can", "usa", "usa", "can", "usa"),
                         subUnit=c(41, 42, 3, 42, 41))

#5. Put together----
bcr.country <- rbind(bcr.ca, bcr.usa, bcr.can4142, bcr.usa41423, bcr.usa414232, bcr.can8182) |> 
  anti_join(bcr.remove)

#BCR###############

#1. Get list of predictions----
preds <- data.frame(predpath = list.files(file.path(root, "output", "predictions"), pattern="*.tiff", full.names=TRUE),
                    predfile = list.files(file.path(root, "output", "predictions"), pattern="*.tiff")) |> 
  separate(predfile, into=c("spp", "bcr", "boot", "year"), sep="_", remove=FALSE) |>  
  mutate(year = as.numeric(str_sub(year, -100, -6)),
         boot = as.numeric(boot))

#2. Get list of models----
mods <- data.frame(modpath = list.files(file.path(root, "output", "bootstraps"), pattern="*.R", full.names=TRUE),
                   modfile = list.files(file.path(root, "output", "bootstraps"), pattern="*.R")) |> 
  separate(modfile, into=c("spp", "bcr", "boot"), sep="_", remove=FALSE) |>  
  mutate(boot = as.numeric(str_sub(boot, -100, -3)))

#3. Put together----
todo <- inner_join(preds, mods)

#Check
nrow(preds)==nrow(todo)

#4. Set up loop----
for(i in 1:nrow(todo)){
  
  #5. Read in files----
  pred.i <- rast(preds$predpath[i])
  mod.i <- load(mods$modpath[i])
  
  
  
}
