# ---
# title: National Models 5.0 - package predictions
# author: Elly Knight
# created: April 29, 2024
# ---

#NOTES################################

# This script calculates means and standard error (CV) across bootstraps of model predictions, limits the predictions by range limits calculated from the dataset, and packages with the sampling distance layer.

# Means and ses are calculated for two extents:
#1. Subunit: this the output of `07.Predict.R` run on compute canada, with the predictions trimmed to the BCR boundary.
#2. Mosaic: this is the output of `08.MosaicPredictions.R`

# We communicate model prediction variance via standard error instead of coefficient of variation or confidence intervals for several reasons: 1) we only want to provide one raster band within the stack for starage size and download efficiency, 2) while CV is better for visualizing variance, it does not have interpretable units, and so SE was prioritized as the better option after conversation with data users. R functions will be made available in the BAMexploreR package to derive CV and CIs from the SE band. CIs from the SE band will be approximations that in some cases may overestimate the lower CI due to a nonnormal distribution of bootstraps, but this is uncommon and is preferred to storing CI rasters.

# Range limitation is done using the output of `gis/Probable_Range_Extent.R`

# Products are stacked, attributed, and saved in a hierarchical structure that is accessed by the BAMExploreR package. Do not change this file structure or the attribution approach for future versions.

# This script also standardizes all projects to EPSG:3978 NAD 83 / Canada Atlas Lambert

# Each averaged raster is also masked by water bodies

# Products from 1985 are excluded from packaging for 3 reasons: 1) There are very few sampling points in the dataset prior to 1993, the 1985 SCANFI products used for prediction are less reliable, 3) there are very few covariate layers that are available for those years of prediction

# This script uses a template to resample from 5072 to 3978 due to changes in projection decisions part way through the modelling process. Future versions will not require this step because they should use 3978 from the onset.


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
water_ca <- read_sf(file.path(root, "gis", "WaterMask_Canada.shp"))
water_us <- read_sf(file.path(root, "gis", "WaterMask_US.shp"))

#7. Subunit polygons----
print("* Getting bcrs *")
bcr.all <- read_sf(file.path(root, "gis", "Subregions_unbuffered.shp"))
bcr.country <- bcr.all |> 
  st_transform("EPSG:3978") 

#8. Mosaic polygons----
print("* Mosaicing bcrs *")

akbox <- st_as_sfc(st_bbox(c(xmin= -3693930,
                             ymin = 2187083,
                             xmax = -1617275,
                             ymax = 4562389),
                           crs = st_crs(bcr.all)))

bcr.can <- bcr.all |> 
  dplyr::filter(country=="can") |> 
  st_union()  |> 
  st_transform("EPSG:3978") |> 
  vect() |> 
  aggregate() |> 
  fillHoles() |> 
  st_as_sf()

bcr.ak <- bcr.all |> 
  dplyr::filter(bcr %in% c("usa41423", "usa2", "usa40", "usa43", "usa5")) |> 
  st_crop(akbox) |> 
  st_union()  |> 
  st_transform("EPSG:3978") |> 
  vect() |> 
  aggregate() |> 
  fillHoles() |> 
  st_as_sf()

bcr.48 <- bcr.all |> 
  dplyr::filter(bcr %in% c("usa5", "usa9", "usa10", "usa11", "usa13", "usa14", "usa23", "usa28")) |> 
  st_difference(akbox) |> 
  st_union() |> 
  st_transform("EPSG:3978") |> 
  vect() |> 
  aggregate() |> 
  fillHoles() |> 
  st_as_sf() 

#9. Load packages on clusters----
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(sf))
tmpcl <- clusterEvalQ(cl, library(tidyverse))
tmpcl <- clusterEvalQ(cl, library(terra))

#10. Data limit mask ----
limit <- read_sf(file.path(root, "gis", "DataLimitationsMask.shp")) |> 
  st_transform("EPSG:3978")

#11. Data package ----
load(file.path(root, "data", "04_NM5.0_data_stratify.Rdata"))
rm(cov, covlist, bcrlist, birdlist, bootlist)

#12. Truncation values ----
q <- read.csv(file.path(root, "data", "SpeciesPredictionTruncationValues.csv"))

#FUNCTION###########

#1. Set up the loop----
brt_package <- function(i){
  
  spp.i <- loop$spp[i]
  year.i <- loop$year[i]
  bcr.i <- loop$bcr[i]
  
  #2. Get the BCR boundary & correct water mask ----
  if(!bcr.i %in% c("Canada", "Alaska", "Lower48")){
    
    sf.i <- bcr.country |> 
      dplyr::filter(bcr==bcr.i) |> 
      vect()
    
  } else {
    
    if(bcr.i=="Canada"){sf.i <- vect(bcr.can)}
    if(bcr.i=="Alaska"){sf.i <- vect(bcr.ak)}
    if(bcr.i=="Lower48"){sf.i <- vect(bcr.48)}

  }
  
  if(bcr.i=="Canada" | str_sub(bcr.i, 1, 3)=="can"){
    water <- water_ca
  } else {water <- water_us}
  
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
  rast.i <- try(rast(files.i$predpath))
  
  if(inherits(rast.i, "try-error")){return(NULL)}
  
  #5. Truncate by count quantile ----
  
  #Get truncation
  qsp <- q[q$spp==spp.i,]$q
  
  #truncate
  truncate.i <- clamp(rast.i, upper = qsp, values=TRUE)
  
  #6. Calculate mean----
  mn.i <- app(truncate.i, mean, na.rm=TRUE)
  
  #7. Secondary upper truncation ----
  q99 <- global(mn.i, quantile, probs=0.999, na.rm=TRUE)[1,1]
  
  truncate2.i <- clamp(truncate.i, upper=q99, values=TRUE)
  mn2.i <- clamp(mn.i, upper=q99, values=TRUE)
    
  #8. Truncate values below 99.9% of the population estimate to zero ----
  mndf.i <- as.data.frame(mn2.i, xy=TRUE) |> 
    arrange(desc(mean)) |>
    mutate(mean_prop = mean/sum(mean, na.rm=TRUE),
           mean_cum = cumsum(mean_prop),
           mean = ifelse(mean_cum > 0.999, 0, mean))
  
  q0 <- min(mndf.i[mndf.i$mean > 0,]$mean)
  
  mean.i <- mndf.i |> 
    dplyr::select(x, y, mean) |> 
    rast(type="xyz", crs(mn2.i)) |> 
    project("EPSG:3978", res=1000)
  
  truncate3.i <- ifel(truncate2.i < q0, 0, truncate2.i)
  
  #9. Calculate sd----
  sd.i <- app(truncate3.i, sd, na.rm=TRUE) |> 
    project("EPSG:3978", res=1000) |> 
    resample(mean.i)

  #10. Stack----
  stack.i <- c(mean.i, sd.i)

  #11. Mask outside range----
  range.i <- rast(file.path(root, "gis", "ranges", paste0(spp.i, ".tif"))) |>
    resample(stack.i)
  
  mask.i <- stack.i * range.i
  mask.i[is.na(mask.i)] <- 0
  
  #12. Read in the sampling distance layers----
  sample.i <- try(rast(files.i$samplepath))
  
  if(inherits(sample.i, "try-error")){return(NULL)}
  
  
  #13. Calculate mean sampling distance----
  samplemn.i <- app(sample.i, mean, na.rm=TRUE) |> 
    project("EPSG:3978", res=1000) |> 
    resample(mean.i)
  
  #14. Stack again ----
  stack2.i <- c(mask.i, samplemn.i)
  names(stack2.i) <- c("mean", "standard_deviation", "detection_distance")

  #15. Mask by water and NA layer ----
  out.i <- stack2.i |> 
    crop(vect(limit), mask=TRUE) |> 
    crop(sf.i, mask=TRUE) |> 
    mask(vect(water), inverse=TRUE)
  
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
  
  #1.8 Save----
  writeRaster(out.i, filename = file.path(root, "output", "10_packaged", spp.i, bcr.i, paste0(spp.i, "_", bcr.i, "_", year.i, ".tif")), overwrite=TRUE, datatype = "FLT4S")
  
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

#2. Check which have been run----
done <- data.frame(file = list.files(file.path(root, "output", "10_packaged"), pattern="*.tif", recursive=TRUE))  |> 
  separate(file, into=c("sppfolder", "bcrfolder", "spp", "bcr", "year", "filetype"), remove=FALSE) |> 
  mutate(year = as.numeric(year)) |> 
  dplyr::select(-filetype) |> 
  dplyr::filter(sppfolder!="SamplingReliability")

#3. Make the todo list----
#remove species that we are omitting for now
loop <- sampled |> 
  anti_join(done) |> 
  arrange(-year, spp, bcr)

#PACKAGE########

#1. Export objects to clusters----
tmpcl <- clusterExport(cl, c("loop", "sampled", "bcr.can", "bcr.ak", "bcr.48", "bcr.country", "brt_package", "root", "water_ca", "water_us", "limit", "q"))

#2. Run BRT function in parallel----
print("* Packaging *")
packaged <- parLapply(cl,
                      X=1:nrow(loop),
                      fun=brt_package)

