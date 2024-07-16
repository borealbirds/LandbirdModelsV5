# ---
# Title: Mosaicking prediction rasters by bootstrap with extrapolation mask & border blend
# Author: Anna Drake, adapted by Elly Knight
# Date: March 2024
# ---

#NOTES################################

# Extrapolation analysis is run for each species-BCR-bootstrap. Because of time-varying
# prediction rasters (e.g. biomass variables), year needs to be specified as well.

#------------------------------------------------
# This code does the following:
# 1) Identify which of these areas have alternate prediction layers that are 
# *not* extrapolated (originating from the buffer predictions of neighbouring BCRs)
# 2) Where a better option exists, clip out the extrapolated areas at the BCR level 
# 3) Where no better option exists retain these areas - for non-overlapping areas this will be the extrapolation of 1 model, for overlap areas where all models have no non-extrapolated data, it will be the weighted average output
# 4) Mosaic and output species-and-year-specific region-wide "Extrapolated area" raster that shows areas where the only option was to retain extrapolated predictions
# 5) Mosaic and output the species-specific region-wide predictions using distance weighting to blend BCR predictions and using the clipping from (3) to retain only non-extrapolated predictions in overlap zones that have multiple outputs

# --------------------
# This code requires:
# 1) Downloading the "dsmextra" package from https://github.com/densitymodelling/dsmextra
# 2) USA/CAN BCR overlap rasters and individual BCR weighting rasters produced 
# in "gis/WeightingRasters.R"
# 3) Buffered BCR subunit shapefile produced in "gis/BufferedSubunits.R"
# -------------------

#PREAMBLE####

#1. Load packages -------- 
library(sf)
library(tidyverse)
library(terra)

#2. NAD83(NSRS2007)/Conus Albers projection (epsg:5072) ----
crs <- "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"

#3. Set file paths ----
root <- "G:/Shared drives/BAM_NationalModels5"

#4. Buffered region shapefile----
bcr <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel_Buffered.shp")) |> 
  mutate(bcr = paste0(country, subUnit))

#5. Get overlap raster----
MosaicOverlap <- rast(file.path(root, "gis", "ModelOverlap.tif"))

#6. Data package----
load(file.path(root, "data", "04_NM5.0_data_stratify.R"))

#INVENTORY####

#Take out the BCR aggregation option that's not being used

#1. Get list of predictions----
predicted <- data.frame(path.pred = list.files(file.path(root, "output", "predictions"), pattern="*.tiff", full.names=TRUE),
                     file.pred = list.files(file.path(root, "output", "predictions"), pattern="*.tiff")) |> 
  separate(file.pred, into=c("spp", "bcr", "boot", "year"), sep="_", remove=FALSE) |> 
  mutate(year = as.numeric(str_sub(year, -100, -5))) |> 
  dplyr::filter(!bcr %in% c("can8182", "usa41423", "usa2"))

#2. Get list of extrapolations----
extrapolated <- data.frame(path.extrap = list.files(file.path(root, "output", "extrapolation"), pattern="*.tiff", full.names=TRUE),
                        file.extrap = list.files(file.path(root, "output", "extrapolation"), pattern="*.tiff")) |>
  separate(file.extrap, into=c("spp", "bcr", "boot", "year"), sep="_", remove=FALSE) |>
  mutate(year = as.numeric(str_sub(year, -100, -5))) |> 
  dplyr::filter(!bcr %in% c("can8182", "usa41423", "usa2"))

#3. Get list of mosaics completed----
mosaicked <- data.frame(path.mosaic = list.files(file.path(root, "output", "mosaics"), pattern="*.tiff", full.names=TRUE),
                        file.mosaic = list.files(file.path(root, "output", "mosaics"), pattern="*.tiff")) |> 
  separate(file.mosaic, into=c("spp", "bcr", "boot"), sep="_", remove=FALSE) |> 
  mutate(boot = str_sub(boot, -100, -3))

#4. Make the to-do list----
all <- inner_join(predicted, extrapolated) |> 
  anti_join(mosaicked)

todo <- all |> 
  dplyr::select(spp, boot, year) |> 
  unique()
  
#5. Check against bcr list per species----
spp <- birdlist |> 
  pivot_longer(AGOS:YTVI, names_to="spp", values_to="use") |> 
  dplyr::filter(use==TRUE)  

loop <- data.frame()
for(i in 1:nrow(todo)){
  
  all.i <- all |> 
    dplyr::filter(spp==todo$spp[i],
                  boot==todo$boot[i],
                  year==todo$year[i])
  
  spp.i <- spp |> 
    dplyr::filter(spp==todo$spp[i])
  
  if(nrow(all.i)==nrow(spp.i)){ loop <- rbind(loop, todo[i])}
  
}

#MOSAIC####

#1. Set up the loop----
for(i in 1:nrow(loop)){
  
  start <- Sys.time()
  
  #2. Loop settings----
  spp.i <- loop$spp[i]
  boot.i <- loop$boot[i]
  year.i <- loop$year[i]
  
  #3. Loop file lists----
  predicted.i <- dplyr::filter(predicted, spp==spp.i, boot==boot.i, year==year.i)
  extrapolated.i <- dplyr::filter(extrapolated, spp==spp.i, boot==boot.i, year==year.i)
  
  #Skip loop if rows don't match
  if(nrow(predicted.i)!=nrow(extrapolated.i)){ next } else { loop.i <- full_join(predicted.i, extrapolated.i)}
  
  #2. Import, mosaic, and sum extrapolation rasters-------
  f.mosaic <- lapply(extrapolated.i$path.extrap, rast) |> 
    sprc() |> 
    mosaic(fun="sum") |> 
    crop(ext(MosaicOverlap))
  
  #3. Fix extent of overlap raster----
  overlap.i <- crop(MosaicOverlap, ext(f.mosaic))
  
  #4. Divide extrapolation layers by available prediction layers to determine if alternate exists-----
  # If ==1 there is no suitable raster, if ==0.25-0.75 then a suitable raster exists
  Overlay <- f.mosaic/overlap.i
  
  #5. Produce region-wide extrapolation raster----
  #Areas we retain extrapolation because we have no alternative == 1
  Extrapolation <- Overlay
  Extrapolation[Extrapolation < 1] <- 0
  Extrapolation[Extrapolation > 0] <- 1
  writeRaster(Extrapolation, file.path(root, "output", "mosaics", "extrapolation", paste0(spp.i, "_", boot.i, "_", year.i, ".tiff")), overwrite=T)
  
  #6. Produce correction raster----
  Correction <- Overlay
  Correction[Correction==1] <- 0 #ignore areas where nothing *or* everything is missing 
  Correction[Correction>0] <- 1 #flag areas where non-extrapolated predictions exist in any layer
  
  #7. Set up bcr loop to correct rasters----
  
  MosaicStack<-list() # hold the weighting rasters
  SppStack<-list() #hold the species prediction rasters
  
  for(j in c(1:nrow(predicted.i))){
    
    #8. Read in masking raster----
    mask.j <- rast(loop.i$path.extrap[j])|>
      crop(Correction) #correct edges (Newfoundland is off)
    
    #9. Crop correction raster to bcr----
    cor <- Correction |> 
      crop(mask.j, mask=TRUE)
    
    #10. Multiply by masking raster----
    #if missing and substitute exists (1*1=1), if not missing or no substitute (0*1/1*0=0)
    cor2 <- cor*mask.j
    #invert values to give regions with substitute values no weight
    cori <- classify(cor2, cbind(1,0), others=1)
    
    #11. Get the edge weighting raster----
    w <- rast(file.path(root, "gis", "edgeweights", paste0(loop.i$bcr[j], ".tif"))) |>
      resample(cori)
    
    #12. Modify edge weighting raster to account for removed sections----
    # 0-out areas of extrapolated values where alternate predictions exist
    MosaicStack[[j]] <- w*cori
    
    #13. Read in the prediction----
    p <- rast(loop.i$path.pred[j]) |> 
      terra::project(crs(w)) |> 
      crop(w)
    
    #14. Weight the prediction----
    SppStack[[j]] <- p*w  #apply weighting to the predictions
    
    cat("Finished BCR", j, "of", nrow(loop.i), "\n")
  }
  
  # !sapply... deals with missing layers, can remove ultimately
  CANwideW <- MosaicStack[!sapply(MosaicStack,is.null)]|>
    sprc()|>
    mosaic(fun="sum") #sum weighting (divisor)
  
  CANwideSp <- SppStack[!sapply(SppStack,is.null)]|>
    sprc()|>
    mosaic(fun="sum") # sum weighted predictions
  
  # correct the weighting -------------
  FinalOut <- CANwideSp/CANwideW |> 
    mask(st_transform(bcr, crs(CANwideSp)))
  
  writeRaster(FinalOut, file.path(root, "output", "mosaics", "predictions", paste0(spp.i, "_", boot.i, "_", year.i, ".tiff")), overwrite=T)
  
  duration <- Sys.time() - start
  
  cat("Finished loop", i, "of", nrow(loop), "in", duration, attr(duration, "units"))
  
}
