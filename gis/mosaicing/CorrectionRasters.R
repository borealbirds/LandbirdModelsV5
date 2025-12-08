# ---
# Title: Using extapolation rasters to produce correction masks
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

#4. List of extrapolation rasters----
extrap <- data.frame(file = list.files(file.path(root, "gis", "extrapolation"), pattern="*.tif", recursive = TRUE)) |> 
  separate(file, into=c("bcr", "year", "filetype"), remove=FALSE) |> 
  mutate(year = as.numeric(year),
         path = file.path(root, "gis", "extrapolation", file)) |> 
  dplyr::select(-file)

#5. BCR perimeter ----
bcr <- read_sf(file.path(root, "gis", "BAM_BCR_NationalModel_Unbuffered.shp")) |> 
  mutate(bcr = paste0(country, subUnit)) |> 
  st_transform(crs="ESRI:102001") |> 
  vect()

#1. Set up to loop through years ----
years <- seq(1985, 2020, 5)
for(i in 1:length(years)){
  
  #2. Filter to the year ----
  extrap.i <- dplyr::filter(extrap, year==years[i])
  
  #6. Import, mosaic, and sum extrapolation rasters-------
  f.mosaic <- lapply(loop.i$path.extrap, rast) |> 
    sprc() |> 
    mosaic(fun="sum") |> 
    crop(ext(MosaicOverlap))
  
  #7. Fix extent of overlap raster----
  overlap.i <- crop(MosaicOverlap, ext(f.mosaic))
  
  #8. Divide extrapolation layers by available layers to determine if alternate exists-----
  # If == 1 there is no potential raster, if == 0.25-0.75 then a potential raster exists
  Overlay <- f.mosaic/overlap.i
  
  #9. Produce region-wide extrapolation raster----
  #Areas we retain extrapolation because we have no alternative == 1
  Extrapolation <- Overlay
  Extrapolation[Extrapolation < 1] <- 0
  Extrapolation[Extrapolation > 0] <- 1
  # writeRaster(Extrapolation, file.path(root, "output", "mosaics", "extrapolation", paste0(spp.i, "_", boot.i, "_", year.i, ".tif")), overwrite=T)
  
  #10. Produce correction raster----
  Correction <- Overlay
  Correction[Correction==1] <- 0 #ignore areas where nothing *or* everything is missing 
  Correction[Correction>0] <- 1 #flag areas where non-extrapolated predictions exist in any layer
  
  #8. Save
  writeRaster(Extrap.mn, file.path(root, "output", "10_packaged", "extrapolation", paste0("Extrapolation_", years[i], ".tif")), overwrite=T)
  
  cat(i, "  ")
  
}


