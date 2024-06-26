# ---
# title: National Models 5.0 - make predictions
# author: Elly Knight
# created: March 5, 2024
# ---

#NOTES################################

# This script uses the model output from 07.Bootstrap.R to make spatial predictions for each bootstrap of each species by BCR combination.

# The script requires that the "gis/SubunitStacks.R" script is run prior to build stacks of rasters for each combination of BCR and year of prediction (currently any 5-year period between 1985 and 2020).

# The years for prediction can be set within this script.

# The script then builds a to-do list for the bootstrapped models and desired years of prediction and checks against the predictions that have already been run prior to building the to-do list so that it can be run multiple times.

# The output is a tif for each bootstrap of each species by BCR combination on the response scale (i.e., # of male birds per ha).

#PREAMBLE############################

#1. Load packages----
print("* Loading packages on master *")
library(tidyverse)
library(gbm)
library(parallel)
library(Matrix)
library(terra)

#2. Determine if testing and on local or cluster----
test <- FALSE
cc <- TRUE

#3. Set nodes for local vs cluster----
if(cc){ nodes <- 32}
if(!cc | test){ nodes <- 2}

#4. Create nodes list----
print("* Creating nodes list *")

#For running on cluster
nodeslist <- unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))

#For testing on local
if(!cc){ nodeslist <- nodes }

print(nodeslist)

#5. Create and register clusters----
print("* Creating clusters *")
cl <- makePSOCKcluster(nodeslist, type="PSOCK")

#6. Set root path----
print("* Setting root file path *")
if(cc){root <- "/scratch/ecknight"}
if(!cc){root <- "G:/Shared drives/BAM_NationalModels5"}

tmpcl <- clusterExport(cl, c("root"))

#7. Load packages on clusters----
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(gbm))
tmpcl <- clusterEvalQ(cl, library(tidyverse))
tmpcl <- clusterEvalQ(cl, library(Matrix))
tmpcl <- clusterEvalQ(cl, library(terra))

#WRITE FUNCTION##########

brt_predict <- function(i){
  
  #1. Get model settings---
  bcr.i <- loop$bcr[i]
  spp.i <- loop$spp[i]
  boot.i <- loop$boot[i]
  year.i <- loop$year[i]
  
  #2. Load model----
  load(file.path(root, "output", "bootstraps", paste0(spp.i, "_", bcr.i, "_", boot.i, ".R")))
  
  #3. Load raster stack----
  stack.i <- rast(file.path(root, "gis", "stacks", paste0(bcr.i, "_", year.i, ".tif")))
  stack.i$meth.i <- stack.i$method

  #4. Predict----
  pred.i <- try(terra::predict(model=b.i, object=stack.i, type="response"))

  #7. Save----
  if(class(pred.i)[1]=="SpatRaster"){
    writeRaster(pred.i, file=file.path(root, "output", "predictions", paste0(spp.i, "_", bcr.i, "_", boot.i, "_", year.i, ".tiff")), overwrite=TRUE)
  }
  
}

#8. Export to clusters----
print("* Loading function on workers *")

tmpcl <- clusterExport(cl, c("brt_predict"))

#RUN MODELS#########

#1. Set desired years----
years <- seq(1985, 2020, 5)

#2. Get list of models that are bootstrapped----
booted <- data.frame(path = list.files(file.path(root, "output", "bootstraps"), pattern="*.R", full.names=TRUE),
                    file = list.files(file.path(root, "output", "bootstraps"), pattern="*.R")) |> 
  separate(file, into=c("spp", "bcr", "boot"), sep="_", remove=FALSE) |> 
  mutate(boot = as.numeric(str_sub(boot, -100, -3)))

#3. Create to do list----
#Sort longest to shortest duration to get the big models going first
todo <- booted |> 
  dplyr::select(bcr, spp, boot) |> 
  expand_grid(year=years)

#4. Determine which are already done----
done <- data.frame(path = list.files(file.path(root, "output", "predictions"), pattern="*.tiff", full.names=TRUE),
                   file = list.files(file.path(root, "output", "predictions"), pattern="*.tiff")) |> 
  separate(file, into=c("spp", "bcr", "boot", "year"), sep="_", remove=FALSE) |>  
  mutate(year = as.numeric(str_sub(year, -100, -6)),
         boot = as.numeric(boot))

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
print("* Making predictions *")
mods <- parLapply(cl,
                  X=1:nrow(loop),
                  fun=brt_predict)

#CONCLUDE####

#1. Close clusters----
print("* Shutting down clusters *")
stopCluster(cl)

if(cc){ q() }
