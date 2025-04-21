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
water <- read_sf(file.path(root, "Regions", "Lakes_and_Rivers", "hydrography_p_lakes_v2.shp")) |> 
  dplyr::filter(TYPE %in% c(16, 18)) |> 
  st_transform(crs="ESRI:102001") |> 
  vect()

#4. Subunit polygons----
bcr.country <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel_Unbuffered.shp")) |> 
  mutate(bcr= paste0(country, subUnit)) |> 
  st_transform(crs="ESRI:102001")

#5. Mosaic polygon----
bcr.all <- st_union(bcr.country)

#6. Data package----
load(file.path(root, "data", "04_NM5.0_data_stratify.Rdata"))

#INVENTORY#########

#1. Get list of sets that have been mosaicked----
#We use the mosaics as input because we know they have all subunits and the mosaics and extrapolation complete
#remove 1985 - we're not providing those predictions
mosaics <- data.frame(path = list.files(file.path(root, "output", "mosaics", "predictions"), pattern="*.tiff", full.names=TRUE, recursive = TRUE),
                        file = list.files(file.path(root, "output", "mosaics", "predictions"), pattern="*.tiff", recursive = TRUE)) |> 
  separate(file, into=c("spp", "boot", "year", "filetype"), remove=FALSE) |>  
  mutate(year = as.numeric(year),
         boot = as.numeric(boot)) |> 
  dplyr::filter(year!=1985) |> 
  dplyr::select(-filetype)

#2. Get the list of sampling layers----
sampled <- data.frame(path = list.files(file.path(root, "output", "sampling"), pattern="*.tiff", full.names=TRUE, recursive = TRUE),
                      file = list.files(file.path(root, "output", "sampling"), pattern="*.tiff", recursive = TRUE)) |> 
  separate(file, into=c("spf", "spp", "bcr", "boot", "year", "filetype"), remove=FALSE) |>  
  mutate(year = as.numeric(year),
         boot = as.numeric(boot)) |> 
  dplyr::filter(year!=1985) |> 
  dplyr::select(-filetype, -path, -file, -spf)

#2. Make the todo list----
#Add all 10 bootstraps
todo <- inner_join(mosaics, sampled) |> 
  group_by(spp, year, bcr) |> 
  summarize(boots = n()) |> 
  ungroup() |> 
  dplyr::filter(boots==10)

missing <- inner_join(mosaics, sampled) |> 
  group_by(spp, year, bcr) |> 
  summarize(boots = n()) |> 
  ungroup() |> 
  dplyr::filter(boots < 10) |> 
  dplyr::filter(bcr %in% c("can71", "can81"))

#3. Check which have been run----
done <- data.frame(file = list.files(file.path(root, "output", "packaged"), pattern="*.tiff", recursive=TRUE))  |> 
  separate(file, into=c("sppfolder", "bcrfolder", "spp", "bcr", "year", "filetype"), remove=FALSE) |> 
  mutate(year = as.numeric(year)) |> 
  dplyr::select(-filetype) |> 
  dplyr::filter(bcrfolder!="old")

#4. Make the todo list----
loop <- anti_join(todo, done) |> 
  dplyr::filter(bcr %in% c("mosaic"),
                year==2020,
                spp %in% c("TEWA", "OSFL"))

#PACKAGE###########

#1. Set up the loop----
for(i in 1:nrow(loop)){
  
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
      mutate(predpath = file.path(root, "output", "predictions", spp,
                                  paste0(spp, "_", bcr.i, "_", boot, "_", year, ".tiff")),
             extrappath = file.path(root, "output", "extrapolation", spp,
                                    paste0(spp, "_", bcr.i, "_", boot, "_", year, ".tiff")),
             samplepath = file.path(root, "output", "sampling", spp, 
                                    paste0(spp, "_", bcr.i, "_", boot, "_", year, ".tiff")))
    
  } else {
    
    files.i <- mosaics |> 
      dplyr::filter(spp==spp.i,
                    year==year.i) |> 
      rename(predpath = path) |> 
      mutate(extrappath = file.path(root, "output", "mosaics", "extrapolation", file),
             samplepath = file.path(root, "output", "sampling", spp, 
                                    paste0(spp, "_", bcr.i, "_", boot, "_", year, ".tiff")))
    
  }
  
  #4. Read in the predictions----
  stack.i <- try(rast(files.i$predpath) |> 
                   project("ESRI:102001"))
  
  if(inherits(stack.i, "try-error")) {next}
  
  #5. Truncate to 99.9% quantile----
  q99 <- global(stack.i, fun=function(x) quantile(x, 0.999, na.rm=TRUE))
  
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
  
  #9. Read in the extrapolations----
  extrap.i <- try(rast(files.i$extrappath) |> 
                    project("ESRI:102001"))
  
  if(inherits(extrap.i, "try-error")) {next}
  
  #10. Calculate mean extrapolation----
  extrapmn.i <- mean(extrap.i, na.rm=TRUE) |> 
    resample(mean.i) |>  
    crop(sf.i, mask=TRUE) |> 
    mask(water, inverse=TRUE)
  
  #11. Read in the sampling distance layers----
  sample.i <- try(rast(files.i$samplepath) |> 
                    project("ESRI:102001"))
  
  if(inherits(sample.i, "try-error")) {next}
  
  #12. Calculate mean sampling distance----
  samplemn.i <- mean(sample.i, na.rm=TRUE) |> 
    resample(mean.i) |>  
    crop(sf.i, mask=TRUE) |> 
    mask(water, inverse=TRUE)
  
  #13. Read in range limit shp----
  range.i <- try(read_sf(file.path(root, "gis", "ranges", paste0(spp.i, "_rangelimit.shp"))) |> 
    dplyr::filter(Limit=="0.1% limit") |> 
    st_transform(crs="ESRI:102001"))
  
  if(inherits(range.i, "try-error")) {next}
  
  #14. Zero out mean prediction outside range----
  mask.i <- mask(mean.i, range.i)
  mask.i[is.na(mask.i)] <- 0
  range.i <- crop(mask.i, sf.i, mask=TRUE) |> 
    mask(water, inverse=TRUE)
  
  #15. Stack----
  out.i <- c(mean.i, cv.i, extrapmn.i, samplemn.i, range.i)
  names(out.i) <- c("mean", "cv", "extrapolation", "detections", "range-limited mean")
  
  #16. Add some attributes----
  attr(out.i, "species") <- spp.i
  attr(out.i, "subunit") <- bcr.i
  attr(out.i, "year") <- year.i
  
  #17. Make folders as needed-----
  if(!(file.exists(file.path(root, "output", "packaged", spp.i)))){
    dir.create(file.path(root, "output", "packaged", spp.i))
  }
  
  if(!(file.exists(file.path(root, "output", "packaged", spp.i, bcr.i)))){
    dir.create(file.path(root, "output", "packaged", spp.i, bcr.i))
  }
  
  #18. Save----
  writeRaster(out.i, filename = file.path(root, "output", "packaged", spp.i, bcr.i, paste0(spp.i, "_", bcr.i, "_", year.i, ".tiff")), overwrite=TRUE)
  
  cat("FINISHED", i, "of", nrow(loop), "PACKAGES\n")
  
}
