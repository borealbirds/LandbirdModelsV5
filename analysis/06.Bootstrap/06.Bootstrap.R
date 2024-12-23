# ---
# title: National Models 5.0 - bootstrap models
# author: Elly Knight
# created: March 5, 2024
# ---

#NOTES################################

#This script uses the model output from 06.Tune.R to determine the best learning rate and number of trees for the gbm function. It reads in the summary output from those models and selects the best parameters for each bcr*spp model. It also eliminates bcr*spp models that did not have enough data to achieve at least 1000 trees with a sufficiently high learning rate (lr.stop >= 1e-6 from the previous script).

#This script also uses the model output from 06.Tune.R to remove covariates that explained < 0.1% of the deviation in the model as part of the effort to continue to simplify the models to improve efficiency of the prediction step.

# The script then builds a to-do list for all the combinations of bcr and spp and the requested number of bootstraps and runs the simpler gbm::gbm() function using the selected learning rate and number of trees.

#The script checks the combinations of bcr, spp, and bootstrap that have already been run prior to building the to-do list so that it can be run multiple times.

#The output is an R object with the raw model, a record of the bcr, spp, bootstrap, parameters, and runtime, and a list of the visits that were included for that bootstrap.

#The list of visits in the output is crucial for determining the training and testing data when we do model evaluation.

#Although the `gbm.fit` function might provided faster performance, we use `gbm` here to ensure the model terms match the terms in the prediction raster stack. The prediction and extrapolation scripts are by far the most time-intensive steps in the model building process, and so we prioritize redundancy over speed in the `06.Bootstrap.R` script to ensure predictions are correct. Otherwise, the order of variables in the raster stacks must match those in the model building.

#PREAMBLE############################

#1. Load packages----
print("* Loading packages on master *")
library(tidyverse)
library(gbm)
library(parallel)
library(Matrix)

#2. Determine if testing and on local or cluster----
test <- FALSE
cc <- TRUE

#3. Set nodes for local vs cluster----
if(cc){ nodes <- 48}
if(!cc | test){ nodes <- 2}

#4. Create and register clusters----
print("* Creating clusters *")
cl <- makePSOCKcluster(nodes, type="PSOCK")

#5. Set root path----
print("* Setting root file path *")
if(cc){root <- "/scratch/ecknight/NationalModels"}
if(!cc){root <- "G:/Shared drives/BAM_NationalModels5"}

tmpcl <- clusterExport(cl, c("root"))

#6. Load packages on clusters----
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(gbm))
tmpcl <- clusterEvalQ(cl, library(tidyverse))
tmpcl <- clusterEvalQ(cl, library(Matrix))
tmpcl <- clusterEvalQ(cl, library(stats))

#7. Load data package----
print("* Loading data on master *")

load(file.path(root, "data", "04_NM5.0_data_stratify.R"))

#8. Load data objects----
print("* Loading data on workers *")

tmpcl <- clusterExport(cl, c("bird", "offsets", "cov", "covlist", "bcrlist", "gridlist", "visit"))

#WRITE FUNCTION##########

brt_boot <- function(i){

  t0 <- proc.time()
  
  #1. Get model settings----
  bcr.i <- loop$bcr[i]
  spp.i <- loop$spp[i]
  boot.i <- loop$boot[i]
  lr.i <- loop$lr[i]
  trees.i <- loop$trees[i]
  
  #2. Load full model----
  trymod <- try(load(file=file.path(root, "output", "fullmodels", spp.i, paste0(loop$spp[i], "_", loop$bcr[i], "_", loop$lr[i], ".R"))))
  
  if(inherits(trymod, "try-error")){ return(NULL) }
  
  #3. Make the new covlist----
  #make sure to retain method and year
  covsnew.i <- m.i[["contributions"]] |> 
    dplyr::filter(rel.inf >= 0.1 |
                    var %in% c("year", "method"))
  
  #remove the model
  rm(b.i)
  
  #4. Get visits to include----
  set.seed(boot.i)
  visit.i <- gridlist[bcrlist[,bcr.i],] |> 
    group_by(year, cell) |> 
    mutate(rowid = row_number(),
           use = sample(1:max(rowid), 1)) |> 
    ungroup() |> 
    dplyr::filter(rowid==use)
  
  #5. Get response data (bird data)----
  bird.i <- bird[as.character(visit.i$id), spp.i]
  
  #6. Get covariates and remove the nonsignificant ones----
  cov.i <- cov[cov$id %in% visit.i$id, colnames(cov) %in% covsnew.i$var]
  
  #7. Add year----
  year.i <- visit[visit$id %in% visit.i$id, "year"]
  
  #8. Put together data object----
  dat.i <- cbind(bird.i, year.i, cov.i) |> 
    rename(count = bird.i)
  
  #9. Get offsets----
  off.i <- offsets[offsets$id %in% visit.i$id, spp.i]
  
  #10. Clean up to save space----
  rm(bird.i, year.i, cov.i)
  
  #11. Run model----
  set.seed(boot.i)
  b.i <- try(gbm::gbm(dat.i$count ~ . + offset(off.i),
                  data = dat.i[, -1],
                  n.trees = trees.i,
                  interaction.depth = 3,
                  shrinkage = lr.i,
                  distribution="poisson",
                  keep.data = FALSE,
                  n.cores=1))
  
  #11. Get performance metrics----
  out.i <- loop[i,] |> 
    cbind(data.frame(n = nrow(dat.i),
                     ncount = nrow(dplyr::filter(dat.i, count > 0)),
                     time = (proc.time()-t0)[3]))
  
  #12. Save model----
  if(!(file.exists(file.path(root, "output", "predictions", spp.i)))){
    dir.create(file.path(root, "output", "predictions", spp.i))
  }
  
  save(b.i, out.i, visit.i, file=file.path(root, "output", "bootstraps", spp.i, paste0(spp.i, "_", bcr.i, "_", boot.i, ".R")))
  
}

#13. Export to clusters----
print("* Loading function on workers *")

tmpcl <- clusterExport(cl, c("brt_boot"))

#RUN MODELS###############

#1. Get list of models that are tuned----
tuned <- data.frame(path = list.files(file.path(root, "output", "tuning"), pattern="*.csv", full.names=TRUE, recursive = TRUE),
                    file = list.files(file.path(root, "output", "tuning"), pattern="*.csv", recursive=TRUE)) |> 
  separate(file, into=c("step", "spp", "bcr", "lr"), sep="_", remove=FALSE) |> 
  mutate(lr = as.numeric(str_sub(lr, -100, -5)))  |> 
  dplyr::select(-step) |> 
  inner_join(data.frame(file = list.files(file.path(root, "output", "fullmodels"), pattern="*.R", recursive=TRUE)) |> 
               separate(file, into=c("folder", "bcr", "lr"), sep="_", remove=TRUE) |> 
               separate(folder, into=c("folder", "spp")) |> 
               mutate(lr = as.numeric(str_sub(lr, -100, -3))) |> 
               dplyr::select(-folder))

#2. Get learning rates----
perf <- map_dfr(read.csv, .x=tuned$path) |> 
  dplyr::filter(trees >= 1000 & trees < 10000)

#3. Set number of bootstraps----
boots <- 10

#4. Get the list of ones that need forcing----
notune <- read.csv(file.path(root, "output", "tuning", "SpeciesBCRCombos_NotTuned.csv")) |> 
  dplyr::select(bcr, spp) |> 
  mutate(lr = 1e-05,
         trees = 9999,
         time = mean(perf$time))

#5. Create to do list----
#Sort longest to shortest duration to get the big models going first
todo <- perf |> 
  dplyr::select(bcr, spp, lr, trees, time) |> 
  rbind(notune) |> 
  expand_grid(boot=c(1:boots)) |> 
  arrange(-time)

#6. Determine which are already done----
done <- data.frame(path = list.files(file.path(root, "output", "bootstraps"), pattern="*.R", full.names=TRUE, recursive=TRUE),
                   file = list.files(file.path(root, "output", "bootstraps"), pattern="*.R", recursive=TRUE)) |> 
  separate(file, into=c("folder", "spp", "bcr", "boot"), sep="_", remove=FALSE) |> 
  mutate(boot = as.numeric(str_sub(boot, -100, -3)))

#6. Create final to do list----
if(nrow(done) > 0){
  
  loop <- todo |> 
    anti_join(done)
  
} else { loop <- todo }

#7. Shut down if nothing left to do----
if(nrow(loop)==0){
  print("* Shutting down clusters *")
  stopCluster(cl)
  
  if(cc){ q() }
}

#For testing - take the shortest duration models
if(test) {loop <- arrange(loop, time)[1:2,]}

print("* Loading model loop on workers *")
tmpcl <- clusterExport(cl, c("loop"))

#9. Run BRT function in parallel----
print("* Fitting models *")
mods <- parLapply(cl,
                  X=1:nrow(loop),
                  fun=brt_boot)

#CONCLUDE####

#1. Close clusters----
print("* Shutting down clusters *")
stopCluster(cl)

if(cc){ q() }
