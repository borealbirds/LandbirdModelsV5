# ---
# title: National Models 5.0 - package predictions
# author: Elly Knight
# created: April 29, 2024
# ---

#NOTES################################

# This script calculates means and standard deviation across bootstraps of model predictions, limits the predictions by range limits calculated from the dataset, and packages with the sampling distance layer.

# Means and ses are calculated for two extents:
#1. Subunit: this the output of `07.Predict.R`
#2. Mosaic: this is the output of `08.MosaicPredictions.R`

# We communicate model prediction variance via standard deviation instead of coefficient of variation or confidence intervals for several reasons: 1) we only want to provide one raster band within the stack for starage size and download efficiency, 2) while CV is better for visualizing variance, it does not have interpretable units, and so SE was prioritized as the better option after conversation with data users. R functions will be made available in the BAMexploreR package to derive CV and CIs from the SE band. CIs from the SE band will be approximations that in some cases may overestimate the lower CI due to a nonnormal distribution of bootstraps, but this is uncommon and is preferred to storing CI rasters.

# Products are stacked, attributed, and saved in a hierarchical structure that is accessed by the BAMExploreR package. Do not change this file structure or the attribution approach for future versions.

# Products from 1985 are excluded from mosaicing (and other scripts) for 3 reasons: 1) There are very few sampling points in the dataset prior to 1993, the 1985 SCANFI products used for prediction are less reliable, 3) there are very few covariate layers that are available for those years of prediction

#PREAMBLE############################

#1. Load packages----
print("* Loading packages on master *")
library(tidyverse)
library(terra)
library(sf)
library(parallel)

#2. Determine if on local or cluster----
cc <- FALSE

#3. Set nodes for local vs cluster----
if(cc){ cores <- 32}
if(!cc){ cores <- 2}

#4. Create and register clusters----
print("* Creating clusters *")
print(table(cores))
cl <- makePSOCKcluster(cores, type="PSOCK")
length(clusterCall(cl, function() Sys.info()[c("nodename", "machine")]))

#5. Set root path----
print("* Setting root file path *")
if(cc){root <- "/scratch/ecknight/NationalModels"}
if(!cc){root <- "G:/Shared drives/BAM_NationalModels5"}

#6. Load packages on clusters----
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(sf))
tmpcl <- clusterEvalQ(cl, library(tidyverse))
tmpcl <- clusterEvalQ(cl, library(terra))

#7. Subunit polygons----
bcr.country <- read_sf(file.path(root, "gis", "Subregions_Mosaics_EPSG3978.shp"))

#FUNCTION###########

#1. Set up the loop----
brt_package <- function(i){
  
  spp.i <- loop$spp[i]
  year.i <- loop$year[i]
  bcr.i <- loop$bcr[i]
  
  #2. Get the BCR boundary & correct water mask ----
  sf.i <- bcr.country |> 
    dplyr::filter(bcr==bcr.i)
  
  #4. Read in the predictions----
  rast.i <- try(rast(file.path(root, "output", "10_truncated", spp.i, bcr.i,
                               paste0(spp.i, "_", bcr.i, "_", year.i, ".tif"))))
  
  if(inherits(rast.i, "try-error")){return(NULL)}

  #6. Calculate mean----
  mn.i <- app(rast.i, mean, na.rm=TRUE)
  
  #8. Calculate sd----
  sd.i <- app(rast.i, sd, na.rm=TRUE)

  #12. Read in the sampling distance layers----
  sample.i <- try(rast(file.path(root, "output", "09_sampling", spp.i,
                                 paste0(spp.i, "_", bcr.i, "_", year.i, ".tif"))))
  
  if(inherits(sample.i, "try-error")){return(NULL)}
  
  #13. Calculate mean sampling distance----
  samplemn.i <- app(sample.i, mean, na.rm=TRUE) |> 
    project("EPSG:3978", res=1000) |> 
    crop(vect(sf.i), mask=TRUE) |> 
    resample(mn.i)
  
  #14. Stack again ----
  out.i <- c(mn.i, sd.i, samplemn.i)
  names(out.i) <- c("mean", "standard_deviation", "detection_distance")
  
  #16. Add some attributes----
  attr(out.i, "species") <- spp.i
  attr(out.i, "subunit") <- bcr.i
  attr(out.i, "year") <- year.i
  
  #17. Make folders as needed-----
  if(!(file.exists(file.path(root, "output", "11_packaged", spp.i)))){
    dir.create(file.path(root, "output", "11_packaged", spp.i))
  }
  
  if(!(file.exists(file.path(root, "output", "11_packaged", spp.i, bcr.i)))){
    dir.create(file.path(root, "output", "11_packaged", spp.i, bcr.i))
  }
  
  #1.8 Save----
  writeRaster(out.i, filename = file.path(root, "output", "11_packaged", spp.i, bcr.i, paste0(spp.i, "_", bcr.i, "_", year.i, ".tif")), overwrite=TRUE, datatype = "FLT4S")
  
}

#INVENTORY#########

#1. Get the list of sampling layers----
sampled <- data.frame(file = list.files(file.path(root, "output", "09_sampling"), pattern="*.tif", recursive = TRUE)) |>
  separate(file, into=c("spf", "spp", "bcr", "year", "filetype"), remove=FALSE) |>
  mutate(year = as.numeric(year),
         path = file.path(root, "output", "09_sampling", file)) |>
  dplyr::filter(year!=1985) |>
  dplyr::select(-filetype, -file, -spf)

#2. Get the list of truncated models ----
truncated <- data.frame(file = list.files(file.path(root, "output", "10_truncated"), pattern="*.tif", recursive=TRUE))  |> 
  separate(file, into=c("sppfolder", "bcrfolder", "spp", "bcr", "year", "filetype"), remove=FALSE) |> 
  mutate(year = as.numeric(year)) |> 
  dplyr::select(-filetype)

#2. Check which have been run----
done <- data.frame(file = list.files(file.path(root, "output", "11_packaged"), pattern="*.tif", recursive=TRUE))  |> 
  separate(file, into=c("sppfolder", "bcrfolder", "spp", "bcr", "year", "filetype"), remove=FALSE) |> 
  mutate(year = as.numeric(year)) |> 
  dplyr::select(-filetype) |> 
  dplyr::filter(sppfolder!="SamplingReliability")

#3. Make the todo list----
#remove species that we are omitting for now
loop <- sampled |> 
  inner_join(truncated) |> 
  anti_join(done) |> 
  arrange(-year, spp, bcr)

#4. Shut down if nothing left to do----
if(nrow(loop)==0){
  print("* Shutting down clusters *")
  stopCluster(cl)
  
  if(cc){ q() }
}

#PACKAGE########

#1. Export objects to clusters----
tmpcl <- clusterExport(cl, c("root", "loop", "bcr.country"))

#2. Run BRT function in parallel----
print("* Packaging *")
packaged <- parLapply(cl,
                       X=1:nrow(loop),
                       fun=brt_package)

#CONCLUDE####

#1. Close clusters----
print("* Shutting down clusters *")
stopCluster(cl)

if(cc){ q() }
