# ---
# title: National Models 5.0 - make predictions
# author: Elly Knight
# created: March 5, 2024
# ---

#NOTES################################

#PREAMBLE############################

#1. Load packages----
print("* Loading packages on master *")
library(tidyverse)
library(gbm)
library(parallel)
library(Matrix)
library(terra)
library(sf)

#2. Determine if testing and on local or cluster----
test <- TRUE
cc <- FALSE

#3. Set nodes for local vs cluster----
if(cc){ nodes <- 32}
if(!cc){ nodes <- 2}

#4. Set root path for data on google drive (for local testing)----
root <- "G:/Shared drives/BAM_NationalModels/NationalModels5.0"

#5. Create nodes list----
print("* Creating nodes list *")

#For running on cluster
nodeslist <- unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))

#For testing on local
if(test){ nodeslist <- nodes }

print(nodeslist)

#6. Create and register clusters----
print("* Creating clusters *")
cl <- makePSOCKcluster(nodeslist, type="PSOCK")

#7. Load packages on clusters----
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(gbm))
tmpcl <- clusterEvalQ(cl, library(tidyverse))
tmpcl <- clusterEvalQ(cl, library(Matrix))
tmpcl <- clusterEvalQ(cl, library(terra))
tmpcl <- clusterEvalQ(cl, library(sf))

#8. Load the region layer----
bcr <- read_sf(file.path("data", "BAM_BCR_NationalModel_Buffered.shp"))

#9. Load the water layer----
water <- read_sf(file.path("data", "hydrography_p_lakes_v2.shp")) |> 
  dplyr::filter(TYPE %in% c(16, 18)) |> 
  st_transform(crs=crs(bcr)) |> 
  vect()

#10. Load shapefiles on clusters----
print("* Loading regions on workers *")

if(cc){ tmpcl <- clusterEvalQ(cl, setwd("/home/ecknight/NationalModels")) }
tmpcl <- clusterExport(cl, c("bcr", "water"))

#WRITE FUNCTION##########

brt_predict <- function(i){
  
  t0 <- proc.time()
  
  #1. Get model settings---
  bcr.i <- loop$bcr[i]
  spp.i <- loop$spp[i]
  boot.i <- loop$boot[i]
  year.i <- loop$year[i]
  
  #2. Load model----
  load(paste0("output/bootstraps/", spp.i, "_", bcr.i, "_", boot.i, ".R"))
  
  #3. Load raster stack----
  stack.i <- rast(paste0("stacks/", bcr.i, "_", year.i, ".tif"))
  names(stack.i) <- c("meth.i", names(stack.i)[2:dim(stack.i)[3]])

  #4. Predict----
  pred.i <- terra::predict(model=b.i, object=stack.i)
  rm(stack.i)

  #5. Get the region polygon----
  shp.i <- dplyr::filter(bcr, country==str_sub(bcr.i, 1, 3), subUnit==as.numeric(str_sub(bcr.i, 4, 10))) |>
    vect()

  #6. Mask----
  mask.i <- mask(pred.i, shp.i) |>
    mask(water, inverse=TRUE)
  
  # #7. Clear for RAM----
  # rm(stack.i, shp.i)
  # 
  # #7. Save----
  # writeRaster(mask.i, file=file.path("output/predictions", paste0(spp.i, "_", bcr.i, "_", boot.i, "_", year.i, ".tiff")))
  
}

#8. Export to clusters----
print("* Loading function on workers *")

tmpcl <- clusterExport(cl, c("brt_predict"))

#RUN MODELS#########

#1. Set desired years----
#years <- seq(1985, 2020, 5)
years <- 2020

#2. Get list of models that are tuned----
booted <- data.frame(path = list.files("output/bootstraps", pattern="*.R", full.names=TRUE),
                    file = list.files("output/bootstraps", pattern="*.R")) |> 
  separate(file, into=c("spp", "bcr", "boot"), sep="_", remove=FALSE) |> 
  mutate(boot = as.numeric(str_sub(boot, -100, -3)))

#3. Create to do list----
#Sort longest to shortest duration to get the big models going first
todo <- booted |> 
  dplyr::select(bcr, spp, boot) |> 
  expand_grid(year=years)

#4. Determine which are already done----
done <- data.frame(path = list.files("output/predictions", pattern="*.csv", full.names=TRUE),
                   file = list.files("output/predictions", pattern="*.csv")) |> 
  separate(file, into=c("spp", "bcr", "boot"), sep="_", remove=FALSE) |> 
  mutate(boot = as.numeric(str_sub(boot, -100, -3)))

#5. Create final to do list----
if(nrow(done) > 0){
  
  loop <- todo |> 
    anti_join(done)
  
} else { loop <- todo }

#For testing - take the shortest duration models
if(test) {loop <- loop[1:2,]}

print("* Loading model loop on workers *")
tmpcl <- clusterExport(cl, c("loop"))

#6. Run BRT function in parallel----
print("* Fitting models *")
mods <- parLapply(cl,
                  X=1:nrow(loop),
                  fun=brt_predict)

#CONCLUDE####

#1. Close clusters----
print("* Shutting down clusters *")
stopCluster(cl)

if(cc){ q() }