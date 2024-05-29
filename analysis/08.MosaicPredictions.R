#-----------------------------------------------------
# Title: Mosaicking prediction rasters by bootstrap with extrapolation mask & border blend
# Author: Anna Drake, Elly Knight
# Date: March 2024           
#-----------------------------------------------------------------------------

# Extrapolation analysis is run for each species-BCR-bootstrap. Because of time-varying
# prediction rasters (e.g. biomass variables), year needs to be specified as well.

#------------------------------------------------
# This code does the following:
# 1) Identify areas of a BCR where the BCR-, species- & year-specific suite of predictor variables
# will cause model extrapolation (Type I and Type II) 
# 2) Identify which of these areas have alternate prediction layers that are 
# *not* extrapolated (originating from the buffer predictions of neighbouring BCRs)
# 3) Where a better option exists, clip out the extrapolated areas at the BCR level 
# 4) Where no better option exists retain these areas - for non-overlapping areas this will be
# the extrapolation of 1 model, for overlap areas where all models have no non-extrapolated
# data, it will be the weighted average output
# 5) Mosaic and output species-and-year-specific region-wide "Extrapolated area" raster that 
# shows areas where the only option was to retain extrapolated predictions
# 6) Mosaic and output the species-specific region-wide predictions using distance weighting
# to blend BCR predictions and using the clipping from (3) to retain only non-extrapolated
# predictions in overlap zones that have multiple outputs

# --------------------
# This code requires:
# 1) Downloading the "dsmextra" package from https://github.com/densitymodelling/dsmextra
# 2) USA/CAN BCR overlap rasters and individual BCR weighting rasters produced 
# in "gis/WeightingRasters.R"
# 3) Buffered BCR subunit shapefile produced in "gis/BufferedSubunits.R"
# -------------------

#TO DO: COLLAPSE INTO ONE LOOP####
#TO DO: CONVERT TO COMPUTE CANADA SCRIPT####
#TO DO: MOVE OUTPUT TO OUTPUT FOLDER####
#TO DO: ADD MEAN AND SD ACROSS BOOTSTRAPS#####
#TO DO: ADD EXTRAPOLATION STUDY AREA RASTERS#########

#PREAMBLE####

#1. Load packages -------- 
library(sf)
library(tidyverse)
library(terra)
library(leaflet)
library(dsmextra)

#2. NAD83(NSRS2007)/Conus Albers projection (epsg:5072) ----
crs <- "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"

#3. Set file paths ----
root<-"G:/Shared drives/BAM_NationalModels/NationalModels5.0"

#4. Buffered region shapefile----
bcr <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel_Buffered.shp")) |> 
  mutate(bcr = paste0(country, subUnit))

#5. Load data package with covariates and BCR lists----
load(file.path(root, "Data", "04_NM5.0_data_stratify.R"))

#EXTRAPOLATION ANALYSIS####

#1. Get list of covariates to include----
cov_clean <- cov |> 
  dplyr::select(where(is.numeric))
cov_clean <- cov_clean[names(cov_clean)!="hli3cl_1km"] #remove hli - it is categorical
cov_clean <- cov_clean[names(cov_clean)!="TRI_1km"] #remove tri - not using

#2. Get list of species and years to process----
SPP<-"OVEN"
YR <-1985

#3. Select regions that were modelled----
BCR_lookup <- birdlist[,colnames(birdlist) %in% c("bcr", SPP)]
BCR_lookup <- BCR_lookup[BCR_lookup[,2]==TRUE,]

#TO DO: MOVE STEP 5 & 6 OUTSIDE J LOOP#########
#TO DO: MAKE LIST FOR AT THE BOOTSTRAP LEVEL########

#4. Set up BCR loop----
for(j in 1:nrow(BCR_lookup)){  # select BCR  
  
  #5. Get list of bootstraps and check for completion----
  file <- list.files(file.path(root, "Output","bootstraps"), pattern=paste(SPP,"_", BCR_lookup$bcr[j], ".+", sep=""))
  
  ck <- list.files(file.path(root,"MosaicWeighting", "MaskingRasters"),
                pattern=paste("ExtrapolatedArea_", BCR_lookup$bcr[j], ".*", SPP, "_", YR, ".tif", sep=""))
  
  if(length(ck) == length(file)) {
    cat("extrapolation complete, skipping...\n")
    next
  }
  
  #6. Retrieve list of continuous covariates for that species---
  #Note: covariates are the same across boots, retrieve from one...
  load(file.path(root, "Output", "bootstraps", file[1]))
  
  Variables <- b.i$var.names[b.i$var.names %in% names(cov_clean)]
  
  #clean up
  rm(visit.i, out.i, b.i)

  #7. Load buffered BCR----
  bcr.i <- bcr |> 
    dplyr::filter(bcr==BCR_lookup$bcr[j])
  
  #8. Get the dataframe of raster values for that BCR*YR combination----
  Dataframes <- as.data.frame(rast(file.path(root, "stacks", paste0(BCR_lookup$bcr[j], "_", YR, ".tif"))), xy=TRUE) 
  
  #9. Set up bootstrap loop----
  for (k in 1:length(file)){  #run through each bootstrap
    
    load(file.path(root, "Output","bootstraps",file[k]))

    #10. Get training data for that bootstrap----
    sample <- cov_clean |> 
      dplyr::filter(id %in% visit.i$id) |> 
      dplyr::select(all_of(Variables))

    #11. Create target dataframes----
    target <- Dataframes |> 
      dplyr::select(c("x","y",all_of(Variables))) 
    
    #12. Compute extrapolation----
    Extrapol <- compute_extrapolation(samples = sample,
                                      covariate.names = names(sample),
                                      prediction.grid = target,
                                      coordinate.system = crs,
                                      verbose=F)
    
    #13. Produce binary raster of extrapolation area----
    raster <- as(Extrapol$rasters$mic$all, "SpatRaster")
    raster[raster>0] <- 1  # extrapolated locations in each boot
    raster[raster!=1|is.na(raster)] <- 0 #eliminate NA areas
    raster <- mask(raster,bcr.i) # crop to BCR + buffer
    
    #14. Write raster----
    writeRaster(raster, file.path(root, "MosaicWeighting", "MaskingRasters",
                paste0("ExtrapolatedArea_", BCR_lookup$bcr[j], "_boot", k, "_", SPP, "_", YR, ".tif")), 
                overwrite=T)
    
    # Print progress
    cat("Finished", k, "of", length(file), "bootstraps \n")
    
  }
  
  cat("Processed", j, "of", nrow(BCR_lookup), "BCRs \n")

}

#MOSAIC##########

#using edge weighting & removing extrapolation areas where alternate predictions exist

#1. Get overlap raster----
MosaicOverlap <- rast(file.path(root,"MosaicWeighting","ModelOverlap.tif"))

cat("Mosaicking", YR, "predictions for", SPP,"\n")

#2. Set up bootstrap loop-----
for (k in c(1:10)){
  
  #3. Import, mosaic, and sum extrapolation rasters-------
  f.mosaic <- lapply(list.files(path=file.path(root, "MosaicWeighting", "MaskingRasters"), pattern=paste0("_boot", k, "_", SPP, "_", YR, "*"), full.names = TRUE), terra::rast) |> 
    sprc() |> 
    mosaic(fun="sum") |> 
    crop(ext(MosaicOverlap))
  
  #4. Divide extrapolation layers by available prediction layers to determine if alternate exists-----
  
  #Fix extent of overlap raster
  overlap.i <- crop(MosaicOverlap, ext(f.mosaic))
  
  # If ==1 there is no suitable raster, if ==0.25-0.75 then a suitable raster exists
  Overlay <- f.mosaic/overlap.i
  plot(Overlay)
  
  #5. Produce region-wide extrapolation raster----

  #Areas we retain extrapolation because we have no alternative == 1
  Extrapolation <- Overlay
  Extrapolation[Extrapolation < 1] <- 0
  plot(Extrapolation)
  writeRaster(Extrapolation, file.path(root, "MosaicWeighting", paste0("ExtrapolatedArea_", "boot", k, "_", SPP, "_", YR, ".tif")), overwrite=T)
  
  #6. Produce correction raster----
  Correction <- Overlay
  Correction[Correction==1] <- 0 #ignore areas where nothing *or* everything is missing 
  Correction[Correction>0] <- 1 #flag areas where non-extrapolated predictions exist in any layer
  plot(Correction)
  
  #7. Set up bcr loop to correct rasters----
  MosaicStack<-list() # hold the weighting rasters
  SppStack<-list() #hold the species prediction rasters
  for(j in c(1:nrow(BCR_lookup))){  # select BCR  
    
    c <- list.files(file.path(root, "MosaicWeighting","MaskingRasters"),
                  paste0("ExtrapolatedArea_", BCR_lookup$bcr[j],"_boot", k, "_", SPP, "_", YR, ".tif")) 
    
    if(length(c)==0){ 
      cat("Missing", BCR_lookup$Names[j], "extrapolation output for", YR, SPP, "boot",k,"\n")
      next}
    
    #8. Read in masking raster----
    rast.j <- terra::rast(file.path(root, "MosaicWeighting", "MaskingRasters", c))|>
      crop(Correction) #correct edges (Newfoundland is off)
    
    #9. Crop correction raster to bcr----
    cor <- Correction |> 
      crop(rast.j) |> 
      mask(rast.j)
    
    #10. Multiply by masking raster----
    #if missing and substitute exists (1*1=1), if not missing or no substitute (0*1/1*0=0)
    cor2 <- cor*rast.j
    #invert values to give regions with substitute values no weight
    cori <- classify(cor2, cbind(1,0), others=1)
    
    #11. Modify weighting raster to account for removed sections---- 
    w <- terra::rast(file.path(root, "MosaicWeighting",  "BCR_Weighting", paste0("EdgeWeighting_", BCR_lookup$bcr[j], ".tif"))) |> 
      resample(cori)
    
    # 0-out areas of extrapolated values where alternate predictions exist
    MosaicStack[[j]] <- w*cori
    
    # Apply BCR weighting raster to Prediction raster 
    plist <- list.files(path=file.path(root, "Output", "predictions"), pattern=paste(SPP, "_", BCR_lookup$bcr[j], "_", k, "_", YR, ".tiff", sep=""), full.names = TRUE)
    
    if (length(plist)==0) {cat("Missing bootstrap", k, "predictions for", SPP, BCR_lookup$Names[j], "in", YR,"\n")
      next} 
    
    p <- terra::rast(plist) |> 
      terra::project(crs(w)) |> 
      crop(w)
    SppStack[[j]] <- p*w  #apply weighting to the predictions
    
    cat("Finished BCR", j, "of", nrow(BCR_lookup), "\n")
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

  writeRaster(FinalOut, file.path(root, "MosaicWeighting", paste0("PredictionMosaic_boot", k, "_", "_", SPP, "_", YR, ".tif")), overwrite=T)
  
}
