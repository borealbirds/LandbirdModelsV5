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

#TO DO FOR V6: Fix crs. Is not perfectly consistent with CRS used in previous scripts.

#PREAMBLE####

#1. Load packages -------- 
print("* Loading packages on master *")
library(sf)
library(tidyverse)
library(terra)
library(parallel)

#2. Determine if testing and on local or cluster----
cc <- FALSE

#3. Set nodes for local vs cluster----
if(cc){ cores <- 16}
if(!cc){ cores <- 1} #can't run on more than 1 local core without swamping RAM

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
bcr <- read_sf(file.path(root, "gis", "Subregions_unbuffered.shp"))

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
  
  #3. Loop file lists----
  predicted.i <- dplyr::filter(predicted, spp==spp.i, year==year.i)
  
  #4. Determine required blank predictions for unmodelled BCRs----
  zeros.i <- zeros |> 
    anti_join(predicted.i)
  
  #5. Add zero files to available files----
  loop.i <- predicted.i |> 
    dplyr::select(-file.pred) |> 
    mutate(path.extrap = file.path(root, "gis", "extrapolation", paste0(bcr, "_", year.i, ".tif"))) |> 
    rbind(zeros.i |> 
            rename(path.pred = path.zero) |> 
            mutate(path.extrap = path.pred,
                   spp = spp.i,
                   year = year.i))
  
  #6. Set up bcr loop to correct rasters----
  SppStack <- list() #hold the species prediction rasters
  MosaicStack <- list() #hold the weighting rasters
  for(j in c(1:nrow(loop.i))){
    
    #7. Get the edge weighting raster----
    w <- rast(file.path(root, "gis", "edgeweights", paste0(loop.i$bcr[j], ".tif")))

    #Add to the list for division
    MosaicStack[[j]] <- w
    
    #8. Read in the prediction----
    p <- try(rast(loop.i$path.pred[j]))
    
    if(inherits(p, "try-error")){break}
    
    if(ext(p)!=ext(w)){
      p <- p |> 
        resample(w) |> 
        crop(w)
    }
    
    #9. Weight the prediction----
    SppStack[[j]] <- p*w  #apply weighting to the predictions
    
    cat(j, "  ")

  }
  
  #10. Skip to next if there were corrupt files----
  if(inherits(p, "try-error")){
    
    return(NULL)
    
  }

  #11. Sum the two stacks----
  CANwideW <- MosaicStack |> 
    sprc()|>
    mosaic(fun="sum")  #sum weighting (divisor)
  
  CANwideSp <- SppStack |> 
    sprc()|>
    mosaic(fun="sum") # sum weighted predictions
  
  #12. Correct the weighting ----
  FinalOut <- CANwideSp/(CANwideW)
  
  #14. Save----- 
  writeRaster(FinalOut, file.path(root, "output", "08_mosaics_can",paste0(spp.i, "_", year.i, ".tif")), overwrite=T)
  
}

#INVENTORY####

#1. Get list of predictions----
predicted <- data.frame(file.pred = list.files(file.path(root, "output", "07_predictions"), pattern="*.tif", recursive = TRUE)) |> 
  separate(file.pred, into=c("folder", "spp", "bcr", "year", "file"), remove=FALSE) |> 
  mutate(year = as.numeric(year),
         path.pred = file.path(root, "output", "07_predictions", file.pred)) |> 
  dplyr::select(-folder, -file) |> 
  dplyr::filter(str_sub(bcr, 1, 3)=="can")

#2. Get list of mosaics completed----
mosaiced <- data.frame(file.mosaic = list.files(file.path(root, "output", "08_mosaics_can"), pattern="*.tif")) |> 
  separate(file.mosaic, into=c("spp", "year"), sep="_", remove=FALSE) |> 
  mutate(year = as.numeric(str_sub(year, -100, -5)),
         path.mosaic = file.path(root, "output", "08_mosaics_can", file.mosaic))

#3. Make the to-do list----
sppuse <- read.csv(file.path(root, "data", "priority_spp_with_model_performance.csv"))$species_code

todo <- birdlist |> 
  pivot_longer(-bcr, names_to="spp", values_to="use") |> 
  dplyr::filter(use==TRUE, spp %in% sppuse,
                str_sub(bcr, 1, 3)=="can") |> 
  expand_grid(year = seq(1985, 2020, 5))

#4. Remove spp*boot*year combinations that don't have all BCRs----
loop <- inner_join(predicted, todo) |> 
  mutate(na = ifelse(is.na(path.pred), 1, 0)) |> 
  group_by(spp, year) |>
  summarize(nas = sum(na)) |> 
  ungroup() |> 
  dplyr::filter(nas==0) |> 
  anti_join(mosaiced) |> 
  arrange(-year) |> 
  dplyr::filter(year > 1985)

#5. Get the full list of zeroed bcr raster----
zeros <- data.frame(path.zero = list.files(file.path(root, "gis", "zeros"), full.names = TRUE),
                    file.zero = list.files(file.path(root, "gis", "zeros"))) |>
  separate(file.zero, into=c("bcr", "tif")) |> 
  dplyr::select(path.zero, bcr)

#MOSAIC####

#1. Export objects to clusters----
tmpcl <- clusterExport(cl, c("loop", "zeros", "predicted", "bcr", "crs", "brt_mosaic", "root"))

#2. Run BRT function in parallel----
print("* Mosaicing predictions *")
mosaics <- parLapply(cl,
                  X=1:nrow(loop),
                  fun=brt_mosaic)
