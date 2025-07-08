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

#2. Set species subset ----
set <- c(1)
set_sp <- "TEWA"
set_yr <- 2020

#3. Determine if testing and on local or cluster----
test <- FALSE
cc <- TRUE

#4. Set nodes for local vs cluster----
if(cc){ cores <- 8}
if(!cc | test){ cores <- 1}

#5. Create and register clusters----
print("* Creating clusters *")
print(table(cores))
cl <- makePSOCKcluster(cores, type="PSOCK")
length(clusterCall(cl, function() Sys.info()[c("nodename", "machine")]))

#6. Set root path----
print("* Setting root file path *")
if(cc){root <- "/scratch/ecknight/NationalModels"}
if(!cc){root <- "G:/Shared drives/BAM_NationalModels5"}

#7. Get the species list----
sppuse <- read.csv(file.path(root, "data", "priority_spp_with_model_performance.csv")) |> 
  rename(spp = species_code) |> 
  dplyr::select(spp, rerun) |> 
  dplyr::filter(rerun %in% set)

#8. Load packages on clusters----
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(gbm))
tmpcl <- clusterEvalQ(cl, library(tidyverse))
tmpcl <- clusterEvalQ(cl, library(Matrix))
tmpcl <- clusterEvalQ(cl, library(terra))

#WRITE FUNCTION##########

brt_predict <- function(i){
  
  #1. Load bootstraps----
  trymod <- try(load(file=file.path(root, "output", "06_bootstraps", spp.i, paste0(spp.i, "_", bcr.i, ".Rdata"))))
  if(inherits(trymod, "try-error")){ return(NULL) }
  
  #2. Get the model----
  b.i <- b.list[[i]]
  remove(b.list)
  
  #3. Load the raster stack----
  #Apparently faster than exporting, gets around requirement to wrap & unwrap
  stack <- try(rast(file.path(root, "gis", "stacks", paste0(bcr.i, "_", year.i, ".tif"))))
  if(inherits(stack, "try-error")){ return(NULL) }

  #4. Reduce to just the required layers----
  #to save RAM
  stack.i <- try(stack[[b.i$var.names]])
  if(inherits(stack.i, "try-error")){ return(NULL)}
  rm(stack)

  #5. Set prediction file path ----
  #We write to a temp file because it is apparently most efficient and avoids having to wrap, plus it returns a list of file paths, which we can read in and save as a multiband raster
  #We use scratch on the cluster to avoid blowing out the RAM
  tf <- file.path(root, "tempfiles", paste0("pred_", spp.i, "_", bcr.i, "_", year.i, "_", i, ".tif"))
  
  #6. Predict----
  pred.i <- try(terra::predict(model=b.i, object=stack.i, type="response", filename = tf, overwrite=TRUE))
  
  #8. Return output
  return(tf)

}

tmpcl <- clusterExport(cl, c("brt_predict"))

#GET INVENTORY###############

#1. Set desired years----
years <- seq(1985, 2020, 5)

#2. Get list of models that are bootstrapped----
booted <- data.frame(path = list.files(file.path(root, "output", "06_bootstraps"), pattern="*.Rdata", full.names=TRUE, recursive = TRUE),
                    file = list.files(file.path(root, "output", "06_bootstraps"), pattern="*.Rdata", recursive = TRUE)) |> 
  separate(file, into=c("folder", "spp", "bcr", "file"), remove=FALSE)

#3. Create to do list----
#currently set to prioritize
todo <- booted |>
  dplyr::select(bcr, spp) |>
  expand_grid(year=years) |>
  arrange(spp, year, bcr) |> 
  inner_join(sppuse |> 
               dplyr::select(spp, rerun))

#RUN MODELS#########

#1. Determine which are already done----
done <- data.frame(file = list.files(file.path(root, "output", "07_predictions"), pattern="*.tif", recursive = TRUE)) |> 
  separate(file, into=c("folder", "spp", "bcr", "year", "file"), remove=FALSE) |>  
  mutate(year = as.numeric(year)) |> 
  dplyr::select(-folder, -file)

#2. Create final to do list----
if(nrow(done) > 0){
  
  loop <- todo |> 
    anti_join(done) |> 
    dplyr::filter(spp == set_sp)
#    dplyr::filter(year == set_yr)
  
} else { loop <- todo }

#For testing
if(test) {loop <- loop[1:2,]}

#3. Shut down if nothing left to do----
if(nrow(loop)==0){
  print("* Shutting down clusters *")
  stopCluster(cl)
  
  if(cc){ q() }
}

#4. Otherwise set up while loop ----
while(nrow(loop) > 0){
  
  #5. Get model settings ----
  bcr.i <- loop$bcr[1]
  spp.i <- loop$spp[1]
  year.i <- loop$year[1]
  
  #6. Export to cores ----
  print("* Loading function requirements on workers *")
  tmpcl <- clusterExport(cl, c("loop", "brt_predict", "bcr.i", "spp.i", "year.i", "cc", "root", "cores"))
  
  #7. Run the models ----
  print("* Making predictions *")
  if(test){file.list <- parLapply(cl,
                                    X=1:cores,
                                    fun=brt_predict)}
  if(!test){file.list <- parLapply(cl,
                                     X=c(1:32),
                                     fun=brt_predict)}
  
  #8. Tidy up----
  gc()
  
  #9. Error handling ----
  if(is.null(file.list[[1]])){
    
    if(!(file.exists(file.path(root, "output", "07_predictions", spp.i)))){
      dir.create(file.path(root, "output", "07_predictions", spp.i))
    }
    
    write.csv(loop[1,], file=file.path(root, "output", "07_predictions", spp.i, paste0(spp.i, "_", bcr.i, "_", year.i, "_fail.csv")), row.names = FALSE)
    
  } else {
    
    #10. Read in the temp files and name----
    print("* Stacking predictions *")
    pred <- terra::rast(unlist(file.list))
    names(pred) <- paste0("b", seq(1:length(file.list)))
    
    #11. Save model----
    print("* Saving predictions *")
    if(!(file.exists(file.path(root, "output", "07_predictions", spp.i)))){
      dir.create(file.path(root, "output", "07_predictions", spp.i))
    }
    
    terra::writeRaster(pred, file=file.path(root, "output", "07_predictions", spp.i, paste0(spp.i, "_", bcr.i, "_", year.i, ".tif")),
                       overwrite=TRUE)
    
  }
  
  #12. Update the list ----
  done <- data.frame(file = list.files(file.path(root, "output", "07_predictions"), pattern="*.tif", recursive = TRUE)) |> 
    separate(file, into=c("folder", "spp", "bcr", "year", "file"), remove=FALSE) |>  
    mutate(year = as.numeric(year)) |> 
    dplyr::select(-folder, -file)
  
  loop <- loop |> 
    anti_join(done, by=c("bcr", "spp", "year"))
  
}

#CONCLUDE####

#1. Close clusters----
print("* Shutting down clusters *")
stopCluster(cl)

if(cc){ q() }
