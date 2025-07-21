# ---
# Title: Mosaicking prediction rasters by bootstrap with extrapolation mask & border blend
# Author: Anna Drake, adapted by Elly Knight
# Date: March 2024
# ---

#NOTES################################

# Extrapolation analysis is run for each species-BCR-bootstrap. Because of time-varying prediction rasters (e.g. biomass variables), year needs to be specified as well.

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
# 3) Zeroed BCR subunit predictions tifs produced in "gis/ZeroPredictions.R"
# -------------------

# This script collects a list of prediction model objects that will not load and deletes them at the end of the script. The prediction script will need to be rerun for those files.

#TO DO FOR V6: Fix crs as EPSG:5072. Is not perfectly consistent with CRS used in previous scripts.

#PREAMBLE####

#1. Load packages -------- 
print("* Loading packages on master *")
library(sf)
library(tidyverse)
library(terra)
library(parallel)

#2. Determine if testing and on local or cluster----
test <- TRUE
cc <- FALSE

#3. Set nodes for local vs cluster----
if(cc){ cores <- 32}
if(!cc | test){ cores <- 1}

#4. Create and register clusters----
print("* Creating clusters *")
print(table(cores))
cl <- makePSOCKcluster(cores, type="PSOCK")
length(clusterCall(cl, function() Sys.info()[c("nodename", "machine")]))

#5. Set root path----
print("* Setting root file path *")
if(cc){root <- "/scratch/ecknight/NationalModels"}
if(!cc){root <- "G:/Shared drives/BAM_NationalModels5"}

#6. NAD83(NSRS2007)/Conus Albers projection (epsg:5072) ----
crs <- "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"

#8. Data package----
load(file.path(root, "data", "04_NM5.0_data_stratify.Rdata"))
rm(offsets, cov, bcrlist, visit, bird)

#9. BCR perimeter ----
bcr <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel_Unbuffered.shp")) |> 
  mutate(bcr = paste0(country, subUnit))

#10. Load packages on clusters----
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(sf))
tmpcl <- clusterEvalQ(cl, library(tidyverse))
tmpcl <- clusterEvalQ(cl, library(terra))

#WRITE FUNCTION##########

#1. Set up----
brt_mosaic <- function(i){
  
  #2. Loop settings----
  spp.i <- loop$spp[i]
  year.i <- loop$year[i]
  
  #7. Get overlap raster----
  MosaicOverlap <- rast(file.path(root, "gis", "ModelOverlap.tif"))
  
  #3. Loop file lists----
  predicted.i <- dplyr::filter(predicted, spp==spp.i, year==year.i)
  
  #4. Determine required blank predictions for unmodelled BCRs----
  zeros.i <- zeros |> 
    anti_join(predicted.i)
  
  #5. Add zero files to available files----
  loop.i <- predicted.i |> 
    dplyr::select(-file.pred) |> 
    mutate(path.extrap = file.path(root, "MosaicWeighting", "Extrapolation", paste0(bcr, "_", year.i, ".tif"))) |> 
    rbind(zeros.i |> 
            rename(path.pred = path.zero) |> 
            mutate(path.extrap = path.pred,
                   spp = spp.i,
                   year = year.i))
  
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
  
  #11. Set up bcr loop to correct rasters----
  MosaicStack<-list() # hold the weighting rasters
  SppStack<-list() #hold the species prediction rasters
  
  for(j in c(1:nrow(loop.i))){
    
    #12. Read in masking raster----
    mask.j <- try(rast(loop.i$path.extrap[j])|>
                    crop(Correction)) #correct edges (Newfoundland is off)
    
    if(inherits(mask.j, "try-error")){break}
    
    #13. Crop correction raster to bcr----
    cor <- Correction |> 
      crop(mask.j, mask=TRUE)
    
    #14. Multiply by masking raster----
    #if missing and substitute exists (1*1=1), if not missing or no substitute (0*1/1*0=0)
    cor2 <- cor*mask.j
    #invert values to give regions with substitute values no weight
    cori <- classify(cor2, cbind(1,0), others=1)
    
    #15. Get the edge weighting raster----
    w <- rast(file.path(root, "gis", "edgeweights", paste0(loop.i$bcr[j], ".tif"))) |>
      resample(cori)
    
    #16. Modify edge weighting raster to account for removed sections----
    # 0-out areas of extrapolated values where alternate predictions exist
    MosaicStack[[j]] <- w*cori
    
    #17. Read in the prediction----
    p <- try(rast(loop.i$path.pred[j]) |> 
               terra::project(crs(w)) |> 
               crop(w))
    
    if(inherits(p, "try-error")){break}
    
    #18. Weight the prediction----
    SppStack[[j]] <- p*w  #apply weighting to the predictions
    
    cat("Finished BCR", j, "of", nrow(loop.i), "\n")
  }
  
  #19. Skip to next if there were corrupt files----
  if(inherits(p, "try-error") | inherits(mask.j, "try-error")){
    
    return(NULL)
    
  }
  
  #20. Clean up ----
  rm(f.mosaic, overlap.i, Overlay, Extrapolation, Correction, w, p)
  gc()
  
  #21. Sum the two stacks----
  CANwideW <- MosaicStack |> 
    sprc()|>
    mosaic(fun="sum")  #sum weighting (divisor)
  
  CANwideSp <- SppStack |> 
    sprc()|>
    mosaic(fun="sum") # sum weighted predictions
  
  #22. Correct the weighting ----
  Weighted <- CANwideSp/(CANwideW)
  
  #23. Fill in NAs----
  #some from hard border transition
  FinalOut <- cover(Weighted, focal(Weighted, w=9, fun=mean, na.policy="only", na.rm=TRUE)) |> 
    mask(st_transform(bcr, crs(CANwideSp)))
  #remaining NAs are waterbodies that will be masked out anyway
  
  #24. Save----- 
  writeRaster(FinalOut, file.path(root, "output", "08_mosaics",paste0(spp.i, "_", year.i, ".tif")), overwrite=T)
  
}

#INVENTORY####

#1. Get list of predictions----
predicted <- data.frame(path.pred = list.files(file.path(root, "output", "07_predictions"), pattern="*.tif", full.names=TRUE, recursive=TRUE),
                        file.pred = list.files(file.path(root, "output", "07_predictions"), pattern="*.tif", recursive = TRUE)) |> 
  separate(file.pred, into=c("folder", "spp", "bcr", "year", "file"), remove=FALSE) |> 
  mutate(year = as.numeric(year)) |> 
  dplyr::select(-folder, -file)

#2. Get list of mosaics completed----
mosaicked <- data.frame(path.mosaic = list.files(file.path(root, "output", "08_mosaics"), pattern="*.tif", full.names=TRUE),
                        file.mosaic = list.files(file.path(root, "output", "08_mosaics"), pattern="*.tif")) |> 
  separate(file.mosaic, into=c("spp", "year"), sep="_", remove=FALSE) |> 
  mutate(year = as.numeric(str_sub(year, -100, -5)))

#3. Make the to-do list----
sppuse <- read.csv(file.path(root, "data", "priority_spp_with_model_performance.csv"))$species_code

todo <- birdlist |> 
  pivot_longer(-bcr, names_to="spp", values_to="use") |> 
  dplyr::filter(use==TRUE, spp %in% sppuse,
                str_sub(bcr, 1, 3)=="can") |> 
  expand_grid(year = seq(1985, 2020, 5))

#4. Remove spp*boot*year combinations that don't have all BCRs----
loop <- full_join(predicted, todo) |> 
  mutate(na = ifelse(is.na(path.pred), 1, 0)) |> 
  group_by(spp, year) |>
  summarize(nas = sum(na)) |> 
  ungroup() |> 
  dplyr::filter(nas==0) |> 
  anti_join(mosaicked)

#5. Get the full list of zeroed bcr raster----
zeros <- data.frame(path.zero = list.files(file.path(root, "gis", "zeros"), full.names = TRUE),
                    file.zero = list.files(file.path(root, "gis", "zeros"))) |>
  separate(file.zero, into=c("bcr", "tif")) |> 
  dplyr::select(path.zero, bcr)

#MOSAIC####

#1. Export objects to clusters----
tmpcl <- clusterExport(cl, c("loop", "zeros", "predicted", "bcr", "crs", "brt_mosaic", "root"))

#2. Run BRT function in parallel----
print("* Fitting models *")
mosaics <- parLapply(cl,
                  X=1:nrow(loop),
                  fun=brt_mosaic)
