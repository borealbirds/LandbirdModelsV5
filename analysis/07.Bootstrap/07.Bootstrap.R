# ---
# title: National Models 5.0 - tune models
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
if(!cc){ nodes <- 2}

#3. Set species subset if desired----
sppuse <- c("OVEN", "OSFL", "CONW", "PAWA")

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
tmpcl <- clusterEvalQ(cl, library(dismo))
tmpcl <- clusterEvalQ(cl, library(tidyverse))
tmpcl <- clusterEvalQ(cl, library(Matrix))

#8. Load data package----
print("* Loading data on master *")

#For running on cluster
if(cc){ load(file.path("04_NM5.0_data_stratify.R")) }

#For testing on local
if(!cc){ load(file.path(root, "Data", "04_NM5.0_data_stratify.R")) }

#9. Wrangle method----
meth <- cbind(cov %>% dplyr::select(id, tagMethod),
              visit %>% dplyr::select(source)) %>% 
  mutate(method = ifelse(source=="eBird", "eBird", as.character(tagMethod)),
         method = factor(method, levels=c("PC", "eBird", "1SPM", "1SPT")))

#10. Load data objects----
print("* Loading data on workers *")

if(cc){ tmpcl <- clusterEvalQ(cl, setwd("/home/ecknight/NationalModels")) }
tmpcl <- clusterExport(cl, c("bird", "offsets", "cov", "birdlist", "covlist", "bcrlist", "gridlist", "meth", "visit"))

#WRITE FUNCTION##########

brt_boot <- function(i){
  
  t0 <- proc.time()
  
  #1. Get model settings---
  bcr.i <- loop$bcr[i]
  spp.i <- loop$spp[i]
  boot.i <- loop$boot[i]
  lr.i <- loop$lr[i]
  trees.i <- loop$trees[i]
  
  #2. Load full model from tuning----
  load(file=file.path("output/fullmodels", paste0(spp.i, "_", bcr.i, "_", lr.i, ".R")))
  
  #3. Identify covariates that explain > 0.1% of variance----
  covlist.i <- m.i[["contributions"]] %>% 
    dplyr::filter(rel.inf >= 0.1)
  
  #4. Get visits to include----
  set.seed(i)
  visit.i <- gridlist[bcrlist[,bcr.i],] %>% 
    group_by(year, cell) %>% 
    mutate(rowid = row_number(),
           use = sample(1:max(rowid), 1)) %>% 
    ungroup() %>% 
    dplyr::filter(rowid==use)
  
  #5. Get response data (bird data)----
  bird.i <- bird[as.character(visit.i$id), spp.i]
  
  #6. Get covariates and remove the nonsignificant ones----
  cov.i <- cov[cov$id %in% visit.i$id, colnames(cov) %in% covlist.i$var]
  
  #7. Get PC vs SPT vs SPM----
  meth.i <- cov[cov$id %in% visit.i$id, "tagMethod"]
  
  #8. Put together data object----
  dat.i <- cbind(bird.i, meth.i, cov.i) %>% 
    rename(count = bird.i)
  
  #9. Get offsets----
  off.i <- offsets[offsets$id %in% visit.i$id, spp.i]
  
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
  saveRDS(b.i, out.i, visit.i, file=file.path("bootstraps", paste0(spp.i, "_", bcr.i, "_", boot.i, ".R")))
  
  #11. Tidy up----
  rm(bcr.i, spp.i, boot.i, lr.i, id.i, visit.i, bird.i, covlist.i, cov.i, meth.i, dat.i, off.i, m.i, b.i)
  
}

#9. Export to clusters----
print("* Loading function on workers *")

tmpcl <- clusterExport(cl, c("brt_boot"))

#RUN MODELS###############

#1. Get BCR & bird combo list----
#filter by list of desired species
bcr.spp <- birdlist %>% 
  pivot_longer(AGOS:YTVI, names_to="spp", values_to="use") %>% 
  dplyr::filter(use==TRUE) %>% 
  dplyr::filter(spp %in% sppuse) %>% 
  dplyr::select(-use)

#2. Reformat covariate list----
bcr.cov <- covlist %>% 
  pivot_longer(ERAMAP_1km:mTPI_1km, names_to="cov", values_to="use") %>% 
  dplyr::filter(use==TRUE)

print("* Loading covariate list on workers *")
tmpcl <- clusterExport(cl, c("bcr.cov"))

#3. Get list of models that are tuned----
tuned <- data.frame(path = list.files("output/tuning", pattern="*.csv", full.names=TRUE),
                    file = list.files("output/tuning", pattern="*.csv")) %>% 
  separate(file, into=c("step", "spp", "bcr", "lr"), sep="_", remove=FALSE) %>% 
  mutate(lr = as.numeric(str_sub(lr, -100, -5)))

#4. Get learning rates----
perf <- map_dfr(read.csv, .x=tuned$path) %>% 
  dplyr::filter(trees >= 1000 & trees < 10000)

#5. Set number of bootstraps----
boots <- 25

#6. Create to do list----
#Sort longest to shortest duration to get the big models going first
todo <- perf %>% 
  dplyr::select(bcr, spp, lr, trees, time) %>% 
  expand_grid(boot=c(1:boots)) %>% 
  arrange(-time)

#7. Determine which are already done----
done <- data.frame(path = list.files("output/bootstraps", pattern="*.csv", full.names=TRUE),
                   file = list.files("output/bootstraps", pattern="*.csv")) %>% 
  separate(file, into=c("step", "spp", "bcr", "boot"), sep="_", remove=FALSE) %>% 
  mutate(boot = as.numeric(str_sub(boot, -100, -5)))

#8. Create final to do list----
if(nrow(done) > 0){
  
  loop <- todo %>% 
    anti_join(done)
  
} else { loop <- todo }

#For testing - take the shortest duration models
if(test) {loop <- loop[nrow(loop)-2:nrow(loop),]}

print("* Loading model loop on workers *")
tmpcl <- clusterExport(cl, c("loop"))

#10. Run BRT function in parallel----
print("* Fitting models *")
mods <- parLapply(cl,
                  X=1:nrow(loop),
                  fun=brt_boot)

#CONCLUDE####

#1. Close clusters----
print("* Shutting down clusters *")
stopCluster(cl)

if(cc){ q() }