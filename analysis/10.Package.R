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
if(!cc){ cores <- 4}

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
  dplyr::filter(country=="can") |> 
  st_union() |> 
  st_area()

#8. Mosaic polygons----
print("* Mosaicing bcrs *")

akbox <- st_as_sfc(st_bbox(c(xmin= -3693930,
                             ymin = 2187083,
                             xmax = -1617275,
                             ymax = 4562389),
                           crs = st_crs(bcr.country)))

bcr.can <- bcr.country |> 
  dplyr::filter(country=="can") |> 
  st_union() |> 
  st_transform("ESRI:102001")

bcr.ak <- bcr.country |> 
  dplyr::filter(bcr %in% c("usa41423", "usa2", "usa40", "usa43", "usa5")) |> 
  st_crop(akbox) |> 
  st_union() |> 
  st_transform("ESRI:102001") |> 
  vect()

bcr.48 <- bcr.country |> 
  dplyr::filter(bcr %in% c("usa5", "usa9", "usa10", "usa11", "usa13", "usa14", "usa23", "usa28")) |> 
  st_difference(akbox) |> 
  st_union() |> 
  st_transform("ESRI:102001") |> 
  vect()

#9. Load packages on clusters----
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(sf))
tmpcl <- clusterEvalQ(cl, library(tidyverse))
tmpcl <- clusterEvalQ(cl, library(terra))
# 
# #10. List species for range limitation ----
# limit <- c("CAWA", "BAWW", "BBWA", "BHVI", "BRCR", "CONW", "EVGR", "GWWA", "HAFL", "LCSP", "LEFL", "PUFI", "VATH", "VESP", "WIWR")

#10. Data limit mask ----
limit <- rast(file.path(root, "gis", "DataExtentMask.tif"))
#FUNCTION###########

#1. Set up the loop----
brt_package <- function(i){
  
  spp.i <- loop$spp[i]
  year.i <- loop$year[i]
  bcr.i <- loop$bcr[i]
  
  #2. Get the BCR boundary----
  if(!bcr.i %in% c("Canada", "Alaska", "Lower48")){
    
    sf.i <- bcr.country |> 
      dplyr::filter(bcr==bcr.i) |> 
      vect()
    
  } else {
    
    if(bcr.i=="Canada"){sf.i <- vect(bcr.can)}
    if(bcr.i=="Alaska"){sf.i <- vect(bcr.ak)}
    if(bcr.i=="Lower48"){sf.i <- vect(bcr.48)}

  }
  
  #3. Get list of predictions, extrapolations, sampling layers----
  if(!bcr.i %in% c("Canada", "Alaska", "Lower48")){
    
    files.i <- sampled |> 
      dplyr::filter(spp==spp.i,
                    year==year.i,
                    bcr==bcr.i) |> 
      rename(samplepath = path) |> 
      mutate(predpath = file.path(root, "output", "07_predictions", spp,
                                  paste0(spp, "_", bcr.i, "_", year, ".tif")))
    
  } else {
    
    files.i <- sampled |> 
      dplyr::filter(spp==spp.i,
                    year==year.i,
                    bcr==bcr.i) |> 
      rename(samplepath = path) |> 
      mutate(predpath = file.path(root, "output", "08_mosaics", bcr.i, spp,
                                    paste0(spp, "_", year, ".tif")))
    
  }
  
  #4. Read in the predictions----
  stack.i <- try(rast(files.i$predpath) |> 
                   project("ESRI:102001"))
  
  if(inherits(stack.i, "try-error")){return(NULL)}
  
  #5. Calculate mean----
  mean.i <- mean(stack.i, na.rm=TRUE) |> 
    crop(sf.i, mask=TRUE)
  
  #Truncate to 99.8% quantile----
  q99mn <- global(mean.i, fun=function(x) quantile(x, 0.998, na.rm=TRUE))
  mean.i[[1]][values(mean.i[[1]]) > q99mn[,1]] <- q99mn[,1]
  
  #6. Calculate sd----
  sd.i <- stdev(stack.i, na.rm=TRUE) |> 
    crop(sf.i, mask=TRUE)
  
  #7. Calculate cv----
  cv.i <- sd.i/(mean.i+0.0000001)
  
  #Truncate to 99.8% quantile----
  q99cv <- global(cv.i, fun=function(x) quantile(x, 0.998, na.rm=TRUE))
  cv.i[[1]][values(cv.i[[1]]) > q99cv] <- q99cv[,1]
  
  #8. Read in the sampling distance layers----
  sample.i <- try(rast(files.i$samplepath) |>
                    project("ESRI:102001"))

  if(inherits(sample.i, "try-error")){return(NULL)}

  #9. Calculate mean sampling distance----
  samplemn.i <- mean(sample.i, na.rm=TRUE) |>
    resample(mean.i) |>
    crop(sf.i, mask=TRUE) |>
    mask(water, inverse=TRUE)
  
  #10. Stack----
  stack.i <- c(mean.i, cv.i, samplemn.i)
  names(stack.i) <- c("mean", "cv", "detectiondistance")
  
  #11. Mask outside range----
  limit.i <- rast(file.path(root, "gis", "ranges", paste0(spp.i, ".tif"))) |>
    project("ESRI:102001") |> 
    resample(stack.i)
  
  mask.i <- stack.i * limit.i
  mask.i[is.na(mask.i)] <- 0

  #12. Mask by water and NA layer ----
  na <- resample(limit, mask.i)
  
  out.i <- mask.i |> 
    crop(sf.i, mask=TRUE) |> 
    mask(na) |> 
    mask(water, inverse=TRUE)
  
  #13. Add some attributes----
  attr(out.i, "species") <- spp.i
  attr(out.i, "subunit") <- bcr.i
  attr(out.i, "year") <- year.i
  
  #14. Make folders as needed-----
  if(!(file.exists(file.path(root, "output", "10_packaged", spp.i)))){
    dir.create(file.path(root, "output", "10_packaged", spp.i))
  }
  
  if(!(file.exists(file.path(root, "output", "10_packaged", spp.i, bcr.i)))){
    dir.create(file.path(root, "output", "10_packaged", spp.i, bcr.i))
  }
  
  #15. Save----
  writeRaster(out.i, filename = file.path(root, "output", "10_packaged", spp.i, bcr.i, paste0(spp.i, "_", bcr.i, "_", year.i, ".tif")), overwrite=TRUE)
  
}

#INVENTORY#########

#1. Get the list of sampling layers----
#remove 1985 - we're not providing those predictions
sampled <- data.frame(file = list.files(file.path(root, "output", "09_sampling"), pattern="*.tif", recursive = TRUE)) |>
  separate(file, into=c("spf", "spp", "bcr", "year", "filetype"), remove=FALSE) |>
  mutate(year = as.numeric(year),
         path = file.path(root, "output", "09_sampling", file)) |>
  dplyr::filter(year!=1985) |>
  dplyr::select(-filetype, -file, -spf)

predicted <-data.frame(file = list.files(file.path(root, "output", "07_predictions"), pattern="*.tif", recursive = TRUE)) |> 
  separate(file, into=c("folder", "spp", "bcr", "year", "file"), remove=FALSE) |> 
  mutate(year = as.numeric(year),
         path = file.path(root, "output", "07_predictions", file)) |> 
  dplyr::select(-folder, -file) |> 
  dplyr::filter(str_sub(bcr, 1, 3)=="can")

mosaiced <- data.frame(file = list.files(file.path(root, "output", "08_mosaics"), pattern="*.tif", recursive = TRUE)) |> 
  separate(file, into=c("country", "sppfolder", "spp", "year", "filetype"), remove=FALSE) |>  
  mutate(year = as.numeric(year),
         path = file.path(root, "output", "09_sampling", file)) |> 
  dplyr::filter(year!=1985,
                country=="Canada") |> 
  mutate(bcr="mosaic") |> 
  dplyr::select(all_of(colnames(predicted)))

#2. Check which have been run----
done <- data.frame(file = list.files(file.path(root, "output", "10_packaged"), pattern="*.tif", recursive=TRUE))  |> 
  separate(file, into=c("sppfolder", "bcrfolder", "spp", "bcr", "year", "filetype"), remove=FALSE) |> 
  mutate(year = as.numeric(year)) |> 
  dplyr::select(-filetype) |> 
  dplyr::filter(sppfolder!="extrapolation")

#3. Make the todo list----
#remove species that we are omitting for now
loop <- sampled |> 
  anti_join(done) |> 
  dplyr::filter(!spp %in% c("BEKI", "SPGR", "SPSA", "CMWA"),
                year >= 2000,
                bcr=="Canada") |> 
  arrange(-year, spp, bcr)

#PACKAGE########

#1. Export objects to clusters----
tmpcl <- clusterExport(cl, c("loop", "sampled", "bcr.all", "bcr.country", "brt_package", "root", "water", "limit"))

#2. Run BRT function in parallel----
print("* Packaging *")
packaged <- parLapply(cl,
                      X=1:nrow(loop),
                      fun=brt_package)
