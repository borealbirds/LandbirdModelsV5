#-----------------------------------
#  Produce weighting rasters for
#  Prediction Mosaicking function
#  Author: Anna Drake, edited by Elly Knight
#  Created: April 2024
#-----------------------------------

#NOTES###########

# This code produces 3 outputs: 
# 1 & 2. A scored raster of BCR overlap for US and Canada + a shapefile of overlap zones 
# 3. Distance-weighting rasters for every BCR (decreasing weight toward edge)
# These outputs are generic and are used within species-prediction mosaic code 
# Distance weighting uses approach by Dominic Roye (dominicroye.github.io)
# All outputs go into "MosaicWeighting" folder
# Note: use latest version of Terra. terra::mosaic function has bug between 7-29 & 7-39 (https://github.com/rspatial/terra/issues/1262)

#PREAMBLE##########

#1. Load packages ----
library(sf)
library(tidyverse)
library(terra)

#2. NAD83(NSRS2007)/Conus Albers projection (epsg:5072) ----
crs <- "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"

#3. Set file paths ----
root <- "G:/Shared drives/BAM_NationalModels5"

#4. Read in buffered and merged region shapefile----
bcr <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel_Buffered.shp"))

#OVERLAP#######

#1. Get extent -----
#can use any prediction raster
Baselayer<-terra::rast(file.path(root, "PredictionRasters/ClimateNormal/AHM_1km.tif"))
Baselayer[!is.na(Baselayer)] <- 1
e <- ext(Baselayer) 

#2. Set up BCR loop----
overlap <- c()
for (i in (1:nrow(bcr))) {
  
  #3. Get the BCR shapefile----
  bcr.i <- bcr[i,]
  
  #4. Make it a raster----
  overlap[[i]] <- Baselayer |>
    project(crs)  |>
    crop(e) |>
    mask(bcr.i) 
  
}

#5. Mosaic----
# == amount of prediction layer overlap
overlap_sum <- sprc(overlap) |> 
  mosaic(fun="sum")  

#6. Combine the two----
writeRaster(overlap_sum, file.path(root, "gis","ModelOverlap.tif"), overwrite=TRUE)

#DISTANCE TO EDGE######

#1. Set up loop for BCR buffering----
for(i in 1:nrow(bcr)){
  
  #2. Get the BCR shapefile----
  bcr.i <- bcr[i,]
  
  #3. Make fishnet----
  #(10km2 to speed processing - could make this finer)
  grid <- st_make_grid(bcr.i, cellsize = 10000, what = "centers") |> 
    st_intersection(bcr.i)
  
  #6. Polygon -> line----
  bcr.line <- st_cast(bcr.i, "MULTILINESTRING")
  
  #7. Get distance from the edge for each point----
  dist <- st_distance(bcr.line, grid)
  
  #8. Dataframe as proportion of max distance----
  df <- data.frame(dist = as.vector(dist)/as.vector(max(dist)),
                   st_coordinates(grid))
  
  #9. Rasterize + match resolution----
  r <- terra::rast(resolution = 10000, ext = ext(bcr.line), crs = crs)
  
  dist_sf <- st_as_sf(df, coords = c("X", "Y")) |> 
    st_set_crs(crs)
  
  dist_raster <- rasterize(dist_sf, r, "dist", fun = mean) |> 
    terra::disagg(fact=10, method="bilinear")
  
  #10. Save----
  writeRaster(dist_raster, file.path(root, "gis", "edgeweights", paste0(bcr.i$bcr, ".tif")), overwrite=TRUE)
  
  cat("Finished", i, "of", nrow(bcr.country), "BCRs \n")
  
}

#11. Check it----
check <- list.files(file.path(root, "gis", "edgeweights"), full.names = TRUE)[-c(17, 30)] |> 
  lapply(rast) |> 
  sprc() |> 
  mosaic(fun="sum")
plot(check)

###################### END OF CODE #####################################
                
