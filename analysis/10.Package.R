# ---
# title: National Models 5.0 - package predictions
# author: Elly Knight
# created: April 29, 2024
# ---

#NOTES################################

# This script calculates means and coefficients of variation (CV) across bootstraps of model predictions, limits the predictions by range limits calculated from the dataset, and packages with the sampling distance layer.

# Means and cvs are calculated for two extents:
#1. Subunit: this the output of `07.Predict.R` run on compute canada, with the predictions trimmed to the BCR boundary.
#2. Mosaic: this is the output of `08.MosaicPredictions.R`

# Range limitation is done using the output of `gis/Probable_Range_Extent.R`

# Products are stacked, attributed, and saved in a hierarchical structure that is accessed by the BAMExploreR package. Do not change this file structure or the attribution approach for future versions.

# This script also standardizes all projects to ESRI:102001 (NAD83 Canada Alberta Equal Conic)

# Each averaged raster is also masked by water bodies

# Products from 1985 are excluded from packaging for 3 reasons: 1) There are very few sampling points in the dataset prior to 1993, the 1985 SCANFI products used for prediction are less reliable, 3) there are very few covariate layers that are available for those years of prediction

#PREAMBLE############################

#1. Load packages----
print("* Loading packages on master *")
library(tidyverse)
library(terra)
library(sf)
library(parallel)

#2. Determine if testing and on local or cluster----
cc <- FALSE

#3. Set nodes for local vs cluster----
if(cc){ cores <- 32}
if(!cc){ cores <- 6}

#4. Create and register clusters----
print("* Creating clusters *")
print(table(cores))
cl <- makePSOCKcluster(cores, type="PSOCK")
length(clusterCall(cl, function() Sys.info()[c("nodename", "machine")]))

#5. Set root path----
print("* Setting root file path *")
if(cc){root <- "/scratch/ecknight/NationalModels"}
if(!cc){root <- "G:/Shared drives/BAM_NationalModels5"}

#6. Get the water layer----
print("* Getting water layer *")
water <- read_sf(file.path(root, "gis", "Lakes_and_Rivers.shp"))

#7. Subunit polygons----
print("* Getting bcrs *")
bcr.country <- read_sf(file.path(root, "gis", "Subregions_unbuffered.shp")) |> 
  dplyr::filter(country=="can")

#8. Mosaic polygon----
print("* Mosaicing bcrs *")
bcr.all <- st_union(bcr.country)

#9. Load packages on clusters----
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(sf))
tmpcl <- clusterEvalQ(cl, library(tidyverse))
tmpcl <- clusterEvalQ(cl, library(terra))

#FUNCTION###########

#1. Set up the loop----
brt_package <- function(i){
  
  spp.i <- loop$spp[i]
  year.i <- loop$year[i]
  bcr.i <- loop$bcr[i]
  
  #2. Get the BCR boundary----
  if(bcr.i!="mosaic"){
    sf.i <- bcr.country |> 
      dplyr::filter(bcr==bcr.i) |> 
      st_transform("ESRI:102001") |> 
      vect()
  } else {
    sf.i <- bcr.all |> 
      st_transform("ESRI:102001") |> 
      vect()
  }
  
  #3. Get list of predictions, extrapolations, sampling layers----
  if(bcr.i!="mosaic"){
    
    files.i <- sampled |> 
      dplyr::filter(spp==spp.i,
                    year==year.i,
                    bcr==bcr.i) |> 
      mutate(predpath = file.path(root, "output", "07_predictions", spp,
                                  paste0(spp, "_", bcr.i, "_", year, ".tif")),
             samplepath = file.path(root, "output", "09_sampling", spp, 
                                    paste0(spp, "_", bcr.i, "_", year, ".tif")))
    
  } else {
    
    files.i <- sampled |> 
      dplyr::filter(spp==spp.i,
                    year==year.i,
                    bcr==bcr.i) |> 
      rename(samplepath = path) |> 
      mutate(predpath = file.path(root, "output", "08_mosaics_can", 
                                    paste0(spp, "_", year, ".tif")))
    
  }
  
  #4. Read in the predictions----
  stack.i <- try(rast(files.i$predpath) |> 
                   project("ESRI:102001"))
  
  if(inherits(stack.i, "try-error")){return(NULL)}
  
  #5. Calculate mean----
  mean.i <- mean(stack.i, na.rm=TRUE) |> 
    crop(sf.i, mask=TRUE) |> 
    mask(water, inverse=TRUE)
  
  #Truncat to 99.8% quantile----
  q99mn <- global(mean.i, fun=function(x) quantile(x, 0.998, na.rm=TRUE))
  mean.i[[1]][values(mean.i[[1]]) > q99mn[,1]] <- q99mn[,1]
  
  #6. Calculate sd----
  #FIX DIVISION BY ZERO
  sd.i <- stdev(stack.i, na.rm=TRUE) |> 
    crop(sf.i, mask=TRUE) |> 
    mask(water, inverse=TRUE)
  
  #7. Calculate cv----
  cv.i <- sd.i/(mean.i+0.0000001)
  
  #Truncate to 99.8% quantile----
  q99cv <- global(cv.i, fun=function(x) quantile(x, 0.998, na.rm=TRUE))
  cv.i[[1]][values(cv.i[[1]]) > q99cv] <- q99cv[,1]
  
  # #8. Read in the sampling distance layers----
  # sample.i <- try(rast(files.i$samplepath) |> 
  #                   project("ESRI:102001"))
  # 
  # if(inherits(sample.i, "try-error")){return(NULL)}
  # 
  # #9. Calculate mean sampling distance----
  # samplemn.i <- mean(sample.i, na.rm=TRUE) |> 
  #   resample(mean.i) |>  
  #   crop(sf.i, mask=TRUE) |> 
  #   mask(water, inverse=TRUE)
  
  #10. Stack----
  stack.i <- c(mean.i, cv.i)
  
  #11. Mask outside range----
  limit <- c("CAWA", "BAWW", "BBWA", "BHVI", "BRCR", "CONW", "EVGR", "GWWA", "HAFL", "LCSP", "LEFL", "PUFI", "VATH", "VESP", "WIWR")
  
  if(spp.i %in% limit){
    limit.i <- try(read_sf(file.path(root, "gis", "ranges", paste0(spp.i, "_rangelimit.shp"))) |>
                     dplyr::filter(Limit=="0.1% limit") |>
                     st_transform(crs="ESRI:102001"))
    
    if(inherits(limit.i, "try-error")){return(NULL)}
    
    range.i <- mask(stack.i, limit.i)
    range.i[is.na(range.i)] <- 0
  } else {range.i <- stack.i}

  #12. Mask by water and NA layer ----
  na <- rast(file.path(root, "gis", "DataLimitationsMask.tif"))
  na2 <- resample(na, range.i)
  out.i <- range.i |> 
    mask(water, inverse=TRUE) |> 
    mask(na2)
  
  names(out.i) <- c("mean", "cv")
  
  #13. Add some attributes----
  attr(out.i, "species") <- spp.i
  attr(out.i, "subunit") <- bcr.i
  attr(out.i, "year") <- year.i
  
  #14. Make folders as needed-----
  if(!(file.exists(file.path(root, "output", "10_packaged_can", spp.i)))){
    dir.create(file.path(root, "output", "10_packaged_can", spp.i))
  }
  
  if(!(file.exists(file.path(root, "output", "10_packaged_can", spp.i, bcr.i)))){
    dir.create(file.path(root, "output", "10_packaged_can", spp.i, bcr.i))
  }
  
  #15. Save----
  writeRaster(out.i, filename = file.path(root, "output", "10_packaged_can", spp.i, bcr.i, paste0(spp.i, "_", bcr.i, "_", year.i, ".tif")), overwrite=TRUE)
  
}

#INVENTORY#########

#1. Get the list of sampling layers----
#remove 1985 - we're not providing those predictions
# sampled <- data.frame(file = list.files(file.path(root, "output", "09_sampling"), pattern="*.tif", recursive = TRUE)) |> 
#   separate(file, into=c("spf", "spp", "bcr", "year", "filetype"), remove=FALSE) |>  
#   mutate(year = as.numeric(year),
#          path = file.path(root, "output", "09_sampling", file)) |> 
#   dplyr::filter(year!=1985) |> 
#   dplyr::select(-filetype, -file, -spf)

predicted <-data.frame(file = list.files(file.path(root, "output", "07_predictions"), pattern="*.tif", recursive = TRUE)) |> 
  separate(file, into=c("folder", "spp", "bcr", "year", "file"), remove=FALSE) |> 
  mutate(year = as.numeric(year),
         path = file.path(root, "output", "07_predictions", file)) |> 
  dplyr::select(-folder, -file) |> 
  dplyr::filter(str_sub(bcr, 1, 3)=="can")

mosaiced <- data.frame(file = list.files(file.path(root, "output", "08_mosaics_can"), pattern="*.tif", recursive = TRUE)) |> 
  separate(file, into=c("spp", "year", "filetype"), remove=FALSE) |>  
  mutate(year = as.numeric(year),
         path = file.path(root, "output", "09_sampling", file)) |> 
  dplyr::filter(year!=1985) |> 
  dplyr::select(-filetype, -file) |> 
  mutate(bcr="mosaic")

sampled <- rbind(predicted, mosaiced)

#2. Check which have been run----
done <- data.frame(file = list.files(file.path(root, "output", "10_packaged_can"), pattern="*.tif", recursive=TRUE))  |> 
  separate(file, into=c("sppfolder", "bcrfolder", "spp", "bcr", "year", "filetype"), remove=FALSE) |> 
  mutate(year = as.numeric(year)) |> 
  dplyr::select(-filetype) |> 
  dplyr::filter(sppfolder!="extrapolation")

#3. Make the todo list----
#remove species that we are omitting for now
loop <- sampled |> 
  anti_join(done) |> 
  dplyr::filter(!spp %in% c("BEKI", "SPGR", "SPSA", "CMWA"),
                year > 1985) |> 
  dplyr::filter(bcr=="mosaic") |> 
  arrange(-year, spp, bcr)

#PACKAGE########

#1. Export objects to clusters----
tmpcl <- clusterExport(cl, c("loop", "sampled", "bcr.all", "bcr.country", "brt_package", "root", "water"))

#2. Run BRT function in parallel----
print("* Packaging *")
packaged <- parLapply(cl,
                      X=1:nrow(loop),
                      fun=brt_package)
