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

#2. Determine if testing and on local or cluster----
test <- TRUE
cc <- TRUE

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

  #7. Save----
  #set options to not write the .aux.xml file
  rgdal::setCPLConfigOption("GDAL_PAM_ENABLED", "FALSE")
  writeRaster(pred.i, file=file.path("output/predictions", paste0(spp.i, "_", bcr.i, "_", boot.i, "_", year.i, ".tiff")), overwrite=TRUE)
  
}

#8. Export to clusters----
print("* Loading function on workers *")

tmpcl <- clusterExport(cl, c("brt_predict"))

#RUN MODELS#########

#1. Set desired years----
years <- seq(1985, 2020, 5)

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