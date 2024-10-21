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
# 4) Zeroed BCR subunit predictions tifs produced in "gis/ZeroPredictions.R"
# -------------------

#TO DO: FILL IN REST OF BCRS WITH ZERO PREDICTIONS##########

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

#1. Get list of predictions----
predicted <- data.frame(path.pred = list.files(file.path(root, "output", "predictions"), pattern="*.tiff", full.names=TRUE, recursive=TRUE),
                        file.pred = list.files(file.path(root, "output", "predictions"), pattern="*.tiff", recursive = TRUE)) |> 
  separate(file.pred, into=c("folder", "spp", "bcr", "boot", "year", "file"), remove=FALSE) |>  
  mutate(year = as.numeric(year),
         boot = as.numeric(boot)) |> 
  dplyr::select(-folder, -file)

#2. Get list of extrapolations----
extrapolated <- data.frame(path.extrap = list.files(file.path(root, "output", "extrapolation"), pattern="*.tiff", full.names=TRUE, recursive=TRUE),
                           file.extrap = list.files(file.path(root, "output", "extrapolation"), pattern="*.tiff", recursive = TRUE)) |> 
    separate(file.extrap, into=c("folder", "spp", "bcr", "boot", "year", "file"), remove=FALSE) |>  
    mutate(year = as.numeric(year),
           boot = as.numeric(boot)) |> 
    dplyr::select(-folder, -file) |> 
  dplyr::filter(!bcr %in% c("can8182", "usa41423", "usa2"))

#3. Get list of mosaics completed----
mosaicked <- data.frame(path.mosaic = list.files(file.path(root, "output", "mosaics", "predictions"), pattern="*.tiff", full.names=TRUE),
                        file.mosaic = list.files(file.path(root, "output", "mosaics", "predictions"), pattern="*.tiff")) |> 
  separate(file.mosaic, into=c("spp", "boot", "year"), sep="_", remove=FALSE) |> 
  mutate(year = as.numeric(str_sub(year, -100, -5)),
         boot = as.numeric(boot))

#4. Figure out which spp*boot*year*bcr are predicted and mosaicked----
all <- inner_join(predicted, extrapolated)

#5. Make the to-do list----
sppuse <- read.csv(file.path(root, "data", "priority_spp_with_model_performance.csv"))$species_code

todo <- birdlist |> 
  pivot_longer(AGOS:YTVI, names_to="spp", values_to="use") |> 
  dplyr::filter(use==TRUE, spp %in% sppuse) |> 
  expand_grid(year = seq(1985, 2020, 5),
              boot = seq(1, 10, 1))

#6. Remove spp*boot*year combinations that don't have all BCRs----
loop <- full_join(all, todo) |> 
    mutate(na = ifelse(is.na(file.pred), 1, 0)) |>
    group_by(spp, boot, year) |>
    summarize(nas = sum(na)) |>
    ungroup() |>
    dplyr::filter(nas==0)

#5. Check against bcr list per species----

# #remove spp*bcr combinations never hit a satisfactory learning rate in model tuning
# norun <- read.csv(file.path(root, "output", "SpeciesBCRCombos_NotTuned.csv"))
# 
# #remove spp that had more than 4 bcr combiantions that never hit a satisfactory learning rate in model tuning
# norun_spp <- norun |> 
#   group_by(spp) |> 
#   summarize(n = n()) |> 
#   dplyr::filter(n > 4)
# 
# spp <- birdlist |> 
#   pivot_longer(AGOS:YTVI, names_to="spp", values_to="use") |> 
#   dplyr::filter(use==TRUE) |> 
#   anti_join(norun) |> 
#   anti_join(norun_spp)

# #This takes a couple minutes to run
# loop <- data.frame()
# for(i in 1:nrow(todo)){
#   
#   all.i <- all |> 
#     dplyr::filter(spp==todo$spp[i],
#                   boot==todo$boot[i],
#                   year==todo$year[i])
#   
#   spp.i <- spp |> 
#     dplyr::filter(spp==todo$spp[i])
#   
#   if(nrow(all.i)==nrow(spp.i)){ loop <- rbind(loop, todo[i,])}
#   
# }

#MOSAIC####

#1. Get the full list of zeroed bcr raster----
zeros <- data.frame(path = list.files(file.path(root, "gis", "zeros"), full.names = TRUE),
                    file = list.files(file.path(root, "gis", "zeros"))) |>
  separate(file, into=c("bcr", "tif")) |> 
  dplyr::select(path, bcr)

#2. Make an object to catch corrupt file errors----
corrupt <- data.frame()

#3. Set up the loop----
for(i in 1:nrow(loop)){
  
  start <- Sys.time()
  
  #4. Loop settings----
  spp.i <- loop$spp[i]
  boot.i <- loop$boot[i]
  year.i <- loop$year[i]
  
  #5. Loop file lists----
  predicted.i <- dplyr::filter(predicted, spp==spp.i, boot==boot.i, year==year.i) |> 
    inner_join(spp)
  extrapolated.i <- dplyr::filter(extrapolated, spp==spp.i, boot==boot.i, year==year.i) |> 
    inner_join(spp)
  
  #6. Skip loop if rows don't match-----
  if(nrow(predicted.i)!=nrow(extrapolated.i)){ next } 
  
  #7. Determine required blank predictions for unmodelled BCRs----
  zeros.i <- zeros |> 
    anti_join(predicted.i)
  
  #8. Add zero files to available files----
  loop.i <- full_join(predicted.i, extrapolated.i) |> 
    dplyr::select(-file.pred, -file.extrap, -use) |> 
    rbind(zeros.i |> 
            rename(path.pred = path) |> 
            mutate(path.extrap = path.pred,
                   spp = spp.i,
                   boot = boot.i,
                   year = year.i))
    
  #9. Import, mosaic, and sum extrapolation rasters-------
  f.mosaic <- lapply(loop.i$path.extrap, rast) |> 
    sprc() |> 
    mosaic(fun="sum") |> 
    crop(ext(MosaicOverlap))
  
  #10. Fix extent of overlap raster----
  overlap.i <- crop(MosaicOverlap, ext(f.mosaic))
  
  #11. Divide extrapolation layers by available layers to determine if alternate exists-----
  # If == 1 there is no potential raster, if == 0.25-0.75 then a potential raster exists
  Overlay <- f.mosaic/overlap.i
  
  #12. Produce region-wide extrapolation raster----
  #Areas we retain extrapolation because we have no alternative == 1
  Extrapolation <- Overlay
  Extrapolation[Extrapolation < 1] <- 0
  Extrapolation[Extrapolation > 0] <- 1
  writeRaster(Extrapolation, file.path(root, "output", "mosaics", "extrapolation", paste0(spp.i, "_", boot.i, "_", year.i, ".tiff")), overwrite=T)
  
  #13. Produce correction raster----
  Correction <- Overlay
  Correction[Correction==1] <- 0 #ignore areas where nothing *or* everything is missing 
  Correction[Correction>0] <- 1 #flag areas where non-extrapolated predictions exist in any layer
  
  #14. Set up bcr loop to correct rasters----
  MosaicStack<-list() # hold the weighting rasters
  SppStack<-list() #hold the species prediction rasters
  
  for(j in c(1:nrow(loop.i))){
    
    #15. Read in masking raster----
    mask.j <- try(rast(loop.i$path.extrap[j])|>
      crop(Correction)) #correct edges (Newfoundland is off)
    
    if(inherits(mask.j, "try-error")){break}
    
    #16. Crop correction raster to bcr----
    cor <- Correction |> 
      crop(mask.j, mask=TRUE)
    
    #17. Multiply by masking raster----
    #if missing and substitute exists (1*1=1), if not missing or no substitute (0*1/1*0=0)
    cor2 <- cor*mask.j
    #invert values to give regions with substitute values no weight
    cori <- classify(cor2, cbind(1,0), others=1)
    
    #18. Get the edge weighting raster----
    w <- rast(file.path(root, "gis", "edgeweights", paste0(loop.i$bcr[j], ".tif"))) |>
      resample(cori)
    
    #19. Modify edge weighting raster to account for removed sections----
    # 0-out areas of extrapolated values where alternate predictions exist
    MosaicStack[[j]] <- w*cori
    
    #20. Read in the prediction----
    p <- try(rast(loop.i$path.pred[j]) |> 
      terra::project(crs(w)) |> 
      crop(w))
    
    if(inherits(p, "try-error")){break}
    
    #21. Weight the prediction----
    SppStack[[j]] <- p*w  #apply weighting to the predictions
    
    cat("Finished BCR", j, "of", nrow(loop.i), "\n")
  }
  
  #22. Skip to next if there were corrupt files----
  if(inherits(p, "try-error") | inherits(mask.j, "try-error")){
    
    print("INPUT RASTER ERROR")
    
    corrupt <- rbind(corrupt, loop[i,])
    
    next
    
  }
  
  #23. Sum the two stacks----
  CANwideW <- MosaicStack |> 
    sprc()|>
    mosaic(fun="sum") #sum weighting (divisor)
  
  CANwideSp <- SppStack |> 
    sprc()|>
    mosaic(fun="sum") # sum weighted predictions
  
  #24. Correct the weighting ----
  FinalOut <- CANwideSp/(CANwideW) |> 
    mask(st_transform(bcr, crs(CANwideSp)))
  
  #25. Save-----
  writeRaster(FinalOut, file.path(root, "output", "mosaics", "predictions", paste0(spp.i, "_", boot.i, "_", year.i, ".tiff")), overwrite=T)
  
  duration <- Sys.time() - start
  
  cat("Finished loop", i, "of", nrow(loop), "in", duration, attr(duration, "units"))
  
}

#19. Delete the corrupt files----

