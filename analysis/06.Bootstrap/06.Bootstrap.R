# ---
# title: National Models 5.0 - bootstrap models
# author: Elly Knight
# created: March 5, 2024
# ---

#NOTES################################

#This script uses the model output from 06.Tune.R to determine the best learning rate and number of trees for the gbm function. It reads in the summary output from those models and selects the best parameters for each bcr*spp model. It also eliminates bcr*spp models that did not have enough data to achieve at least 1000 trees with a sufficiently high learning rate (lr.stop >= 1e-6 from the previous script).

#This script also uses the model output from 05.Tune.R to remove covariates that explained < 0.1% of the deviation in the model as part of the effort to continue to simplify the models to improve efficiency of the prediction step.

# The script then builds a to-do list for all the tuned combinations of bcr and spp and runs the simpler gbm::gbm() function using the selected learning rate and number of trees.

#The script checks the combinations of bcr & spp that have already been run prior to building the to-do list so that it can be run multiple times.

#The output is an Rdata object with a list of the bootstrapped raw models for each bcr & spp combination and a dataframe of model paramters for all bootstaps.

#Although the `gbm.fit` function might provide faster performance, we use `gbm` here to ensure the model terms match the terms in the prediction raster stack. The prediction and extrapolation scripts are by far the most time-intensive steps in the model building process, and so we prioritize redundancy over speed in the `06.Bootstrap.R` script to ensure predictions are correct. Otherwise, the order of variables in the raster stacks must match those in the model building.

#TO DO FOR V6: Consider the most efficient place to reduce the models by filtering out covs < 0.01 relative influence

#PREAMBLE############################

#1. Load packages----
print("* Loading packages on master *")
library(tidyverse)
library(gbm)
library(parallel)
library(Matrix)

#2. Determine if testing and on local or cluster----
test <- FALSE
cc <- FALSE

#3. Set cores for local vs cluster----
if(cc){ cores <- 32 }
if(!cc | test){ cores <- 2}

#4. Create and register clusters----
print("* Creating clusters *")
print(table(cores))
cl <- makePSOCKcluster(cores, type="PSOCK")
length(clusterCall(cl, function() Sys.info()[c("nodename", "machine")]))

#5. Set root path----
print("* Setting root file path *")
if(cc){root <- "/scratch/ecknight/NationalModels"}
if(!cc){root <- "G:/Shared drives/BAM_NationalModels5"}

tmpcl <- clusterExport(cl, c("root"))

#6. Get the species list----
sppuse <- read.csv(file.path(root, "data", "priority_spp_with_model_performance.csv")) |> 
  dplyr::filter(rerun==2) |> 
  rename(spp = species_code) |> 
  dplyr::select(spp, rerun)

#7. Load packages on clusters----
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(gbm))
tmpcl <- clusterEvalQ(cl, library(tidyverse))
tmpcl <- clusterEvalQ(cl, library(Matrix))
tmpcl <- clusterEvalQ(cl, library(stats))

#8. Load data package----
print("* Loading data on master *")

load(file.path(root, "data", "04_NM5.0_data_stratify.Rdata"))

#9. Load data objects----
print("* Loading data on workers *")

tmpcl <- clusterExport(cl, c("bird", "offsets", "cov", "covlist", "bcrlist", "bootlist", "visit"))

#WRITE FUNCTION##########

brt_boot <- function(i){
  
  t0 <- proc.time()
  
  #6. Get visits to include----
  visit.i <- bcrlist[bcrlist[,bcr.i]>0, c("id", bcr.i)] |> 
    dplyr::filter(id %in% bootlist[[i + 2]])
  
  #7. Get response data (bird data)----
  bird.i <- bird[as.character(visit.i$id), spp.i]
  
  #8. Get covariates and remove the nonsignificant ones----
  cov.i <- cov[cov$id %in% visit.i$id, colnames(cov) %in% covsnew.i$var]
  
  #9. Add year----
  year.i <- visit[visit$id %in% visit.i$id, "year"]
  
  #10. Put together data object----
  dat.i <- cbind(bird.i, year.i, cov.i) |> 
    rename(count = bird.i)
  
  #11. Get offsets----
  off.i <- offsets[offsets$id %in% visit.i$id, spp.i]
  
  #12. Clean up to save space----
  rm(bird.i, year.i, cov.i)
  
  #13. Run model----
  set.seed(i)
  b.i <- try(gbm::gbm(dat.i$count ~ . + offset(off.i),
                      data = dat.i[, -1],
                      n.trees = trees.i,
                      interaction.depth = 3,
                      shrinkage = lr.i,
                      distribution="poisson",
                      keep.data = FALSE,
                      n.cores=1))
  
  if(inherits(b.i, "try-error")){ return(NULL) }
  
  #14. Add performance as attributes----
  attr(b.i, "bcr") <- bcr.i
  attr(b.i, "spp") <- spp.i
  attr(b.i, "time") <- (proc.time()-t0)[3]
  attr(b.i, "n") <- nrow(dat.i)
  attr(b.i, "ncount") <- nrow(dplyr::filter(dat.i, count > 0))
  
  return(b.i)
}

#16. Export to clusters----
print("* Loading function on workers *")
tmpcl <- clusterExport(cl, c("brt_boot"))

#GET INVENTORY###############

#1. Get list of models that are tuned----
tuned <- data.frame(path = list.files(file.path(root, "output", "05_tuning"), pattern="*.csv", full.names=TRUE, recursive = TRUE),
                    file = list.files(file.path(root, "output", "05_tuning"), pattern="*.csv", recursive=TRUE)) |> 
  dplyr::filter(file!="SpeciesBCRCombos_NotTuned.csv") |> 
  separate(file, into=c("step", "spp", "bcr"), sep="_", remove=FALSE) |> 
  dplyr::select(-step) |>
  mutate(bcr = str_remove(bcr, ".csv")) |> 
  inner_join(data.frame(file = list.files(file.path(root, "output", "05_fullmodels"), pattern="*.Rdata", recursive=TRUE)) |> 
               separate(file, into=c("folder", "bcr"), sep="_", remove=TRUE) |> 
               separate(folder, into=c("folder", "spp")) |> 
               dplyr::select(-folder) |> 
               mutate(bcr = str_remove(bcr, ".Rdata")),
             by=c("spp", "bcr")) |> 
  inner_join(sppuse)

#2. Get learning rates----
perf <- map_dfr(read.csv, .x=tuned$path) |> 
  dplyr::filter(trees >= 1000 & trees < 10000)

#3. Get the list of ones that need forcing----
notune <- read.csv(file.path(root, "output", "05_tuning", "SpeciesBCRCombos_NotTuned.csv")) |> 
  dplyr::select(bcr, spp) |> 
  mutate(lr = 1e-05,
         trees = 9999,
         time = mean(perf$time),
         tuned = FALSE)

#4. Create to do list----
#Sort longest to shortest duration to get the big models going first
todo <- perf |> 
  dplyr::select(bcr, spp, lr, trees, time) |> 
  mutate(tuned = TRUE) |> 
  rbind(notune) |> 
  arrange(-time) |> 
  unique() |> 
  inner_join(sppuse)

#RUN MODELS############

#1. Determine which are already done----
done <- data.frame(file = list.files(file.path(root, "output", "06_bootstraps"), pattern="*.Rdata", recursive=TRUE)) |> 
  separate(file, into=c("folder", "spp", "bcr", "filetype"), remove=FALSE) |> 
  inner_join(sppuse )

#2. To do list----
if(nrow(done) > 0){
  
  loop <- todo |> 
    anti_join(done)
  
} else { loop <- todo }

if(test) {loop <- arrange(loop, time)[1:2,]}

#3. Shut down if nothing left to do----
if(nrow(loop)==0){
  print("* Shutting down clusters *")
  stopCluster(cl)
  
  if(cc){ q() }
}

#4. Otherwise set up while loop ----
cat("* Starting modelling for", nrow(loop), "rows *")
while(nrow(loop) > 0){
  
  #5. Get model settings----
  bcr.i <- loop$bcr[1]
  spp.i <- loop$spp[1]
  lr.i <- loop$lr[1]
  trees.i <- loop$trees[1]
  tuned.i <- loop$tuned[1]
  
  #6. Load full model----
  if(tuned.i==TRUE){
    trymod <- try(load(file=file.path(root, "output", "05_fullmodels", spp.i, paste0(spp.i, "_", bcr.i, ".Rdata"))))
    
    if(inherits(trymod, "try-error")){ return(NULL) }
    
    #7. Make the new covlist if already tuned----
    #make sure to retain method and year
    covsnew.i <- m.i[["contributions"]] |> 
      dplyr::filter(rel.inf >= 0.1 |
                      var %in% c("year", "method"))
    
    #remove the model
    rm(m.i)
    
  }
  
  #8. Make the new covlist if not tuned----
  if(tuned.i==FALSE){
    
    covsnew.i <- covlist |> 
      dplyr::filter(bcr==bcr.i) |> 
      pivot_longer(-bcr, names_to="var", values_to="use") |> 
      dplyr::filter(use==TRUE)
    
  }
  
  #9. Export loop info to cores ----
  print("* Loading model loop info on workers *")
  tmpcl <- clusterExport(cl, c("bcr.i", "spp.i", "lr.i", "trees.i", "tuned.i", "covsnew.i"))
  
  #10. Run the models ----
  print("* Fitting models *")
  b.list <- parLapply(cl,
                      X=1:32,
                      fun=brt_boot)
  
  #11. Save model----
  if(!(file.exists(file.path(root, "output", "06_bootstraps", spp.i)))){
    dir.create(file.path(root, "output", "06_bootstraps", spp.i))
  }
  
  save(b.list, file=file.path(root, "output", "06_bootstraps", spp.i, paste0(spp.i, "_", bcr.i, ".Rdata")))
  
  #12. Update the list ----
  done <- data.frame(file = list.files(file.path(root, "output", "06_bootstraps"), pattern="*.Rdata", recursive=TRUE)) |> 
    separate(file, into=c("folder", "spp", "bcr", "filetype"), remove=FALSE)
  
  loop <- loop |> 
    anti_join(done)
  
}

#CONCLUDE####

#1. Close clusters----
print("* Shutting down clusters *")
stopCluster(cl)

if(cc){ q() }
