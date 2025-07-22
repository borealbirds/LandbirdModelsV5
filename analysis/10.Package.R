# ---
# title: National Models 5.0 - calculate means
# author: Elly Knight
# created: April 29, 2024
# ---

#NOTES################################

# This script calculates means and coefficients of variation (CV) across bootstraps of model predictions, limits the predictions by range limits calculated from the dataset, and packages with other existing layers (sampling & extrapolation)

# Means and cvs are calculated for two extents:
#1. Subunit: this the output of `07.Predict.R` run on compute canada, with the predictions trimmed to the BCR boundary.
#2. Mosaic: this is the output of `09.MosaicPredictions.R`

# Range limitation is done using the output of `gis/Probable_Range_Extent.R`

# Products are stacked, attributed, and saved in a hierarchical structure that is accessed by the BAMExploreR package. Do not change this file structure or the attribution approach for future versions.

# This script also standardizes all projects to ESRI:102001 (NAD83 Canada Alberta Equal Conic)

# Each averaged raster is also masked by water bodies

# Products from 1985 are excluded from packaging for 3 reasons: 1) There are very few sampling points in the dataset prior to 1993, the 1985 SCANFI products used for prediction are less reliable, 3) there are very few covariate layers that are available for those years of prediction

#TO DO: PARALLELIZE

#PREAMBLE############################

#1. Load packages----
library(tidyverse)
library(terra)
library(sf)

#2. Set root path----
root <- "G:/Shared drives/BAM_NationalModels5"

#3. Get the water layer----
water <- read_sf(file.path(root, "gis", "Lakes_and_Rivers", "hydrography_p_lakes_v2.shp")) |> 
  dplyr::filter(TYPE %in% c(16, 18)) |> 
  st_transform(crs="ESRI:102001") |> 
  vect()

#4. Subunit polygons----
bcr.country <- read_sf(file.path(root, "gis", "BAM_BCR_NationalModel_Unbuffered.shp")) |> 
  mutate(bcr= paste0(country, subUnit)) |> 
  st_transform(crs="ESRI:102001")

#5. Mosaic polygon----
bcr.all <- st_union(bcr.country)

#6. Data package----
load(file.path(root, "data", "04_NM5.0_data_stratify.Rdata"))

#INVENTORY#########

#1. Get the list of sampling layers----
#remove 1985 - we're not providing those predictions
sampled <- data.frame(file = list.files(file.path(root, "output", "09_sampling"), pattern="*.tif", recursive = TRUE)) |> 
  separate(file, into=c("spf", "spp", "bcr", "year", "filetype"), remove=FALSE) |>  
  mutate(year = as.numeric(year),
         path = file.path(root, "output", "09_sampling", file)) |> 
  dplyr::filter(year!=1985) |> 
  dplyr::select(-filetype, -path, -file, -spf)

#2. Check which have been run----
done <- data.frame(file = list.files(file.path(root, "output", "10_packaged"), pattern="*.tif", recursive=TRUE))  |> 
  separate(file, into=c("sppfolder", "bcrfolder", "spp", "bcr", "year", "filetype"), remove=FALSE) |> 
  mutate(year = as.numeric(year)) |> 
  dplyr::select(-filetype)

#3. Make the todo list----
loop <- sampled |> 
  anti_join(done) |> 
  dplyr::filter(bcr %in% c("mosaic"),
                year==2020)

#PACKAGE###########

#1. Set up the loop----
for(i in 1:nrow(loop)){
  
  start <- Sys.time()
  
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
    
    files.i <- mosaics |> 
      dplyr::filter(spp==spp.i,
                    year==year.i) |> 
      mutate(predpath = file.path(root, "output", "07_predictions", spp,
                                  paste0(spp, "_", bcr.i, "_", year, ".tif")),
             samplepath = file.path(root, "output", "09_sampling", spp, 
                                    paste0(spp, "_", bcr.i, "_", year, ".tif")))
    
  } else {
    
    files.i <- mosaics |> 
      dplyr::filter(spp==spp.i,
                    year==year.i) |> 
      rename(predpath = path) |> 
      mutate(samplepath = file.path(root, "output", "09_sampling", spp, 
                                    paste0(spp, "_", bcr.i, "_", year, ".tif")))
    
  }
  
  #4. Read in the predictions----
  stack.i <- try(rast(files.i$predpath) |> 
                   project("ESRI:102001"))
  
  if(inherits(stack.i, "try-error")) {next}
  
  #5. Truncate to 99.8% quantile----
  #99.9% still gives Inf for some species
  q99 <- global(stack.i, fun=function(x) quantile(x, 0.998, na.rm=TRUE))
  
  truncate.i <- stack.i
  for(k in 1:nlyr(truncate.i)){
    truncate.i[[k]][values(truncate.i[[k]]) > q99[k,1]] <- q99[k,1]
  }
  
  #6. Calculate mean----
  mean.i <- mean(truncate.i, na.rm=TRUE) |> 
    crop(sf.i, mask=TRUE) |> 
    mask(water, inverse=TRUE)
  
  #7. Calculate sd----
  sd.i <- stdev(truncate.i, na.rm=TRUE) |> 
    crop(sf.i, mask=TRUE) |> 
    mask(water, inverse=TRUE)
  
  #8. Calculate cv----
  cv.i <- sd.i/mean.i
  
  #10. Read in the sampling distance layers----
  sample.i <- try(rast(files.i$samplepath) |> 
                    project("ESRI:102001"))
  
  if(inherits(sample.i, "try-error")) {next}
  
  #11. Calculate mean sampling distance----
  samplemn.i <- mean(sample.i, na.rm=TRUE) |> 
    resample(mean.i) |>  
    crop(sf.i, mask=TRUE) |> 
    mask(water, inverse=TRUE)
  
  #12. Read in range limit shp----
  range.i <- try(read_sf(file.path(root, "gis", "ranges", paste0(spp.i, "_rangelimit.shp"))) |> 
    dplyr::filter(Limit=="0.1% limit") |> 
    st_transform(crs="ESRI:102001"))
  
  if(inherits(range.i, "try-error")) {next}
  
  #13. Zero out mean prediction outside range----
  mask.i <- mask(mean.i, range.i)
  mask.i[is.na(mask.i)] <- 0
  range.i <- crop(mask.i, sf.i, mask=TRUE) |> 
    mask(water, inverse=TRUE)
  
  #14. Stack----
  out.i <- c(range.i, mean.i, cv.i, samplemn.i)
  names(out.i) <- c( "range-limited mean", "mean", "cv", "detections")
  
  #16. Add some attributes----
  attr(out.i, "species") <- spp.i
  attr(out.i, "subunit") <- bcr.i
  attr(out.i, "year") <- year.i
  
  #17. Make folders as needed-----
  if(!(file.exists(file.path(root, "output", "10_packaged", spp.i)))){
    dir.create(file.path(root, "output", "10_packaged", spp.i))
  }
  
  if(!(file.exists(file.path(root, "output", "10_packaged", spp.i, bcr.i)))){
    dir.create(file.path(root, "output", "10_packaged", spp.i, bcr.i))
  }
  
  #18. Save----
  writeRaster(out.i, filename = file.path(root, "output", "10_packaged", "extrapolation", spp.i, bcr.i, paste0(spp.i, "_", bcr.i, "_", year.i, ".tif")), overwrite=TRUE)
  
  end <- Sys.time()
  
  cat("FINISHED", i, "of", nrow(loop), "PACKAGES\n")
  difftime(end, start)
  
}
