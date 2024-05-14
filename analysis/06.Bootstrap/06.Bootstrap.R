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

#PREAMBLE############################

#1. Load packages----
print("* Loading packages on master *")
library(tidyverse)
library(gbm)
library(parallel)
library(Matrix)

#2. Determine if testing and on local or cluster----
test <- TRUE
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
if(cc){root <- "/home/ecknight/NationalModels"}
if(!cc){root <- "G:/Shared drives/BAM_NationalModels/NationalModels5.0"}

tmpcl <- clusterExport(cl, c("root"))

#7. Load packages on clusters----
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(gbm))
tmpcl <- clusterEvalQ(cl, library(tidyverse))
tmpcl <- clusterEvalQ(cl, library(Matrix))
tmpcl <- clusterEvalQ(cl, library(stats))

#8. Load data package----
print("* Loading data on master *")

load(file.path(root, "data", "04_NM5.0_data_stratify.R"))

#9. Load data objects----
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

  #2. Get visits to include----
  set.seed(i)
  visit.i <- gridlist[bcrlist[,bcr.i],] %>% 
    group_by(year, cell) %>% 
    mutate(rowid = row_number(),
           use = sample(1:max(rowid), 1)) %>% 
    ungroup() %>% 
    dplyr::filter(rowid==use)
  
  #3. Get response data (bird data)----
  bird.i <- bird[as.character(visit.i$id), spp.i]
  
  #4. Get covariates and remove the nonsignificant ones----
  cov.i <- cov[cov$id %in% visit.i$id, colnames(cov) %in% covsnew[[i]]$var]
  
  #5. Get PC vs SPT vs SPM vs eBird----
  meth.i <- visit[visit$id %in% visit.i$id, "method"]
  
  #6. Get year----
  year.i <- visit[visit$id %in% visit.i$id, "year"]
  
  #7. Put together data object----
  dat.i <- cbind(bird.i, year.i, meth.i, cov.i) %>% 
    rename(count = bird.i)
  
  #8. Get offsets----
  off.i <- offsets[offsets$id %in% visit.i$id, spp.i]
  
  #9. Clean up to save space----
  rm(visit.i, bird.i, year.i, meth.i, cov.i)
  
  #10. Run model----
  set.seed(i)
  b.i <- try(gbm::gbm(dat.i$count ~ . + offset(off.i),
                  data = dat.i[, -1],
                  n.trees = trees.i,
                  interaction.depth = 3,
                  shrinkage = lr.i,
                  distribution="poisson",
                  keep.data = FALSE,
                  n.cores=1))
  
  #11. Get performance metrics----
  out.i <- loop[i,] %>% 
    cbind(data.frame(n = nrow(dat.i),
                     ncount = nrow(dplyr::filter(dat.i, count > 0)),
                     time = (proc.time()-t0)[3]))
  
  #12. Save----
  save(b.i, out.i, visit.i, file=file.path(root, "output", "bootstraps", paste0(spp.i, "_", bcr.i, "_", boot.i, ".R")))
  
}

#13. Export to clusters----
print("* Loading function on workers *")

tmpcl <- clusterExport(cl, c("brt_boot"))

#RUN MODELS###############

#1. Get list of models that are tuned----
tuned <- data.frame(path = list.files(file.path(root, "output", "tuning"), pattern="*.csv", full.names=TRUE),
                    file = list.files(file.path(root, "output", "tuning"), pattern="*.csv")) %>% 
  separate(file, into=c("step", "spp", "bcr", "lr"), sep="_", remove=FALSE) %>% 
  mutate(lr = as.numeric(str_sub(lr, -100, -5)))  |> 
  inner_join(data.frame(file = list.files(file.path(root, "output", "fullmodels"), pattern="*.R")) |> 
               separate(file, into=c("spp", "bcr", "lr"), sep="_", remove=TRUE) |> 
               mutate(lr = as.numeric(str_sub(lr, -100, -3))))

#2. Get learning rates----
perf <- map_dfr(read.csv, .x=tuned$path) %>% 
  dplyr::filter(trees >= 1000 & trees < 10000)

#3. Set number of bootstraps----
boots <- 10

#4. Create to do list----
#Sort longest to shortest duration to get the big models going first
todo <- perf %>% 
  dplyr::select(bcr, spp, lr, trees, time) %>% 
  expand_grid(boot=c(1:boots)) %>% 
  arrange(-time)

#5. Determine which are already done----
done <- data.frame(path = list.files(file.path(root, "output", "bootstraps"), pattern="*.R", full.names=TRUE),
                   file = list.files(file.path(root, "output", "bootstraps"), pattern="*.R")) %>% 
  separate(file, into=c("spp", "bcr", "boot"), sep="_", remove=FALSE) %>% 
  mutate(boot = as.numeric(str_sub(boot, -100, -3)))

#6. Create final to do list----
if(nrow(done) > 0){
  
  loop <- todo %>% 
    anti_join(done)
  
} else { loop <- todo }

#For testing - take the shortest duration models
if(test) {loop <- arrange(loop, time)[1:2,]}

print("* Loading model loop on workers *")
tmpcl <- clusterExport(cl, c("loop"))

#7. Update the covariate lists (remove covs that explain < 0.01 % of deviance)----
covsnew <- list()
for(i in 1:nrow(loop)){
  
  load(file=file.path(root, "output", "fullmodels", paste0(loop$spp[i], "_", loop$bcr[i], "_", loop$lr[i], ".R")))
  
  covsnew[[i]] <- m.i[["contributions"]] %>% 
    dplyr::filter(rel.inf >= 0.1)
  
}

print("* Loading new covariate lists *")
tmpcl <- clusterExport(cl, c("covsnew"))

#8. Run BRT function in parallel----
print("* Fitting models *")
mods <- parLapply(cl,
                  X=1:nrow(loop),
                  fun=brt_boot)

#CONCLUDE####

#1. Close clusters----
print("* Shutting down clusters *")
stopCluster(cl)

if(cc){ q() }
