# ---
# title: National Models 5.0 - calculate means
# author: Elly Knight
# created: April 29, 2024
# ---

#NOTES################################

# This script calculates means and standard deviation across bootstraps of model predictions for packaging

# Means and sds are calculated for two extends:
#1. Subunit: this the output of `07.Predict.R` run on compute canada, with the predictions trimmed to the BCR boundary.
#2. Mosaic: this is the output of `09.MosaicPredictions.R`

# Each averaged raster is also masked by water bodies

#PREAMBLE############################

#1. Load packages----
library(tidyverse)
library(terra)
library(sf)

#2. Set root path----
root <- "G:/Shared drives/BAM_NationalModels5"

#3. Get the water layer----
water <- read_sf(file.path(root, "Regions", "Lakes_and_Rivers", "hydrography_p_lakes_v2.shp")) |> 
  dplyr::filter(TYPE %in% c(16, 18)) |> 
  st_transform(crs=5072) |> 
  vect()

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

bcr.usa41423 <- bcr.usa |> 
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
  anti_join(bcr.remove) |> 
  mutate(bcr = paste0(country, subUnit))

#SUBUNIT AVERAGING#########

#1. Get list of predictions----
predicted <- data.frame(path.pred = list.files(file.path(root, "output", "predictions"), pattern="*.tiff", full.names=TRUE),
                        file.pred = list.files(file.path(root, "output", "predictions"), pattern="*.tiff")) |> 
  separate(file.pred, into=c("spp", "bcr", "boot", "year"), sep="_", remove=FALSE) |> 
  mutate(year = as.numeric(str_sub(year, -100, -5)))

#2. Remove ones without a full set of bootstraps----
boots <- 10

todo <- predicted |> 
  group_by(spp, bcr, year) |> 
  summarize(n=n()) |> 
  ungroup() |> 
  dplyr::filter(n==boots)

#3. Check which have been run----
done <- data.frame(file.mean = list.files(file.path(root, "output", "averages", "subunit"), pattern="*.tiff"))  |> 
  separate(file.mean, into=c("spp", "bcr", "year"), sep="_", remove=FALSE) |> 
  mutate(year = as.numeric(str_sub(year, -100, -5)))

#4. Make the todo list----
loop <- anti_join(todo, done)

#5. Set up the loop----
for(i in 1:nrow(loop)){
  
  spp.i <- loop$spp[i]
  bcr.i <- loop$bcr[i]
  year.i <- loop$year[i]
  
  #6. Get the list of prediction bootstraps----
  pred.i <- predicted |> 
    dplyr::filter(spp==spp.i,
                  bcr==bcr.i,
                  year==year.i)
  
  #7. Read them in----
  stack.i <- rast(pred.i$path.pred)
  
  #8. Get the BCR boundary----
  sf.i <- bcr.country |> 
    dplyr::filter(bcr==bcr.i) |> 
    vect()
  
  #9. Calculate mean----
  #Also crop and mask
  mean.i <- mean(stack.i, na.rm=TRUE) |> 
    crop(sf.i, mask=TRUE) |> 
    mask(water, inverse=TRUE)
  
  #10. Calculate sd----
  sd.i <- stdev(stack.i, na.rm=TRUE) |> 
    crop(sf.i, mask=TRUE) |> 
    mask(water, inverse=TRUE)
  
  #11. Stack----
  out.i <- c(mean.i, sd.i)
  
  #12. Save----
  writeRaster(out.i, filename = file.path(root, "output", "averages", "subunit", paste0(spp.i, "_", bcr.i, "_", year.i, ".tiff")), overwrite=TRUE)
  
  cat("Finished", i, "of", nrow(loop), "prediction sets\n")
  
}

#MOSAIC AVERAGING#########

#1. Get list of predictions----
predicted <- data.frame(path.pred = list.files(file.path(root, "output", "mosaics"), pattern="*.tiff", full.names=TRUE),
                        file.pred = list.files(file.path(root, "output", "mosaics"), pattern="*.tiff")) |> 
  separate(file.pred, into=c("spp", "boot", "year"), sep="_", remove=FALSE) |> 
  mutate(year = as.numeric(str_sub(year, -100, -5)))

#2. Remove ones without a full set of bootstraps----
boots <- 10

todo <- predicted |> 
  group_by(spp, bcr, year) |> 
  summarize(n=n()) |> 
  ungroup() |> 
  dplyr::filter(n==boots)

#3. Check which have been run----
done <- data.frame(file.mean = list.files(file.path(root, "output", "averages", "mosaic"), pattern="*.tiff"))  |> 
  separate(file.mean, into=c("spp", "year"), sep="_", remove=FALSE) |> 
  mutate(year = as.numeric(str_sub(year, -100, -5)))

#4. Make the todo list----
loop <- anti_join(todo, done)

#5. Set up the loop----
for(i in 1:nrow(loop)){
  
  spp.i <- loop$spp[i]
  year.i <- loop$year[i]
  
  #6. Get the list of prediction bootstraps----
  pred.i <- predicted |> 
    dplyr::filter(spp==spp.i,
                  year==year.i)
  
  #7. Read them in----
  stack.i <- rast(pred.i$path.pred)
  
  #8. Calculate mean----
  #Also crop and mask
  mean.i <- mean(stack.i, na.rm=TRUE) |> 
    crop(sf.i, mask=TRUE) |> 
    mask(water, inverse=TRUE)
  
  #9. Calculate sd----
  #Also crop and mask
  sd.i <- stdev(stack.i, na.rm=TRUE) |> 
    crop(sf.i, mask=TRUE) |> 
    mask(water, inverse=TRUE)
  
  #10. Stack----
  out.i <- c(mean.i, sd.i)
  
  #11. Save----
  writeRaster(out.i, filename = file.path(root, "output", "averages", "mosaic", paste0(spp.i, "_", year.i, ".tiff")), overwrite=TRUE)
  
  cat("Finished", i, "of", nrow(loop), "prediction sets")
  
}
