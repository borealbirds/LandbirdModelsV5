# ---
# title: National Models 5.0 - tune models
# author: Elly Knight
# created: November 8, 2023
# ---

#NOTES################################

#This script uses a grid search of tuning parameters to determine ideal BRT parameters for each model subunit (BCR*country)

#Parameters to tune include interaction depth & learning rate

#Best parameters are defined as the highest deviation explained that produces a model with 2000-10000 trees

#Tuning process is only run once using bootstrap #1, assuming that tuning parameters are relatively insensitive to bootstrap.

#This script is written to be run on the compute canada cluster in parallel. Messages are therefore printed in the output for each step in the script to help with debugging. Use the .sh object to run the script on the cluster.

#There is an object at the beginning of the script that can be set to TRUE if running this script on a local machine, which will then overwrite objects to local testing settings. Set this object to FALSE before running on compute canada.

#Files must be transferred between local computer and compute canada's servers using Globus Connect.

#PREAMBLE############################

#1. Load packages----
print("* Loading packages on master *")
library(tidyverse)
library(dismo)
library(parallel)

#2. Determine if testing and on local or cluster----
test <- TRUE
cc <- TRUE

#3. Set root path for data on google drive (for local testing)----
root <- "G:/Shared drives/BAM_NationalModels/NationalModels5.0"

#4. Create nodes list----
print("* Creating nodes list *")

#For running on cluster
nodeslist <- unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))

#For testing on local
if(test){ nodeslist <- 32 }

print(nodeslist)

#5. Create and register clusters----
print("* Creating clusters *")
cl <- makePSOCKcluster(nodeslist, type="PSOCK")

#6. Load packages on clusters----
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(dismo))
tmpcl <- clusterEvalQ(cl, library(tidyverse))

#7. Load data----
print("* Loading data on master *")

#For running on cluster
if(cc){ load(file.path("04_NM5.0_data_stratify.R")) }

#For testing on local
if(!cc){ load(file.path(root, "Data", "04_NM5.0_data_stratify.R")) }

print("* Loading data on workers *")

if(cc){ tmpcl <- clusterEvalQ(cl, setwd("/home/ecknight/NationalModels")) }
tmpcl <- clusterExport(cl, c("visit", "bird", "offsets", "birdlist", "covlist"))

#WRITE FUNCTION##########

brt_tune <- function(i){
  
  t0 <- proc.time()
  
  #1. Get model settings---
  bcr.i <- loop$bcr[i]
  spp.i <- loop$spp[i]
  boot.i <- loop$boot[i]
  lr.i <- loop$lr[i]
  id.i <- loop$id[i]
  
  #2. Get visits to include----
  visit.i <- visit[bootlist[[bcr.i]][,boot.i+1],]
  
  #3. Get response data (bird data)----
  bird.i <- bird[bird$id %in% visit.i$id, spp.i]
  
  #4. Get covariates----
  covlist.i <- bcr.cov %>% 
    dplyr::filter(bcr==bcr.i)
  cov.i <- visit.i[,colnames(visit.i) %in% covlist.i$cov]
  
  #5. Get PC vs SPT vs SPM----
  meth.i <- factor(visit.i$tagMethod)
  
  #6. Put together data object----
  dat.i <- cbind(bird.i, meth.i, cov.i)
  
  #7. Get offsets----
  off.i <- offsets[offsets$id %in% visit.i$id, spp.i]
  
  #8. Run model----
  m.i <- dismo::gbm.step(data=dat.i,
                         gbm.x=c(2:ncol(dat.i)),
                         gbm.y=1,
                         offset=off.i,
                         tree.complexity = id.i,
                         learning.rate = lr.i,
                         family="poisson")
  
  #9. Get performance metrics----
  out[[i]] <- loop[i,] %>% 
    cbind(data.frame(trees = m.i$n.trees,
                     deviance.mean = m.i$cv.statistics$deviance.mean,
                     deviance.se = m.i$cv.statistics$deviance.se,
                     null = m.i$self.statistics$mean.null,
                     resid = m.i$self.statistics$mean.resid,
                     correlation = m.i$self.statistics$correlation,
                     correlation.mean = m.i$cv.statistics$correlation.mean,
                     correlation.se = m.i$cv.statistics$correlation.se,
                     time = (proc.time()-t0)[3]))
  
  #10. Save----
  save(out, file=file.path("results", "ModelTuning.Rdata"))
  
  #11. Tidy up----
  rm(bcr.i, spp.i, boot.i, lr.i, id.i, visit.i, bird.i, covlist.i, cov.i, meth.i, dat.i, off.i, m.i)
  
}

#9. Export to clusters----
print("* Loading function on workers *")

tmpcl <- clusterExport(cl, c("brt_tune"))

#RUN MODELS###############

#1. Set grid search parameters----

#Learning rate
lr <- c(0.01, 0.005, 0.001, 0.0005, 0.0001)

#Interaction depth
id <- c(3)

#2. Set number of bootstraps----
boot <- 1

#3. Get BCR & bird combo list----
bcr.spp <- birdlist %>% 
  pivot_longer(ABDU:YTVI, names_to="spp", values_to="use") %>% 
  dplyr::filter(use==TRUE)

#4. Reformat covariate list----
bcr.cov <- covlist %>% 
  pivot_longer(ERAMAP_1km:mTPI_1km, names_to="cov", values_to="use") %>% 
  dplyr::filter(use==TRUE)

print("* Loading covariate list on workers *")
tmpcl <- clusterExport(cl, c("bcr.cov"))

#5. Make dataframe of models to run----

#Full combinations
loop <- expand.grid(lr=lr, id=id, boot=boot) %>% 
  merge(bcr.spp) %>% 
  dplyr::select(-use) %>% 
  arrange(bcr, lr, id, spp)

#For testing
if(test) {loop <- loop[1:32,]}

print("* Loading model loop on workers *")
tmpcl <- clusterExport(cl, c("loop"))

#6. Make object to store output----
out <- list()

print("* Loading output object on workers *")
tmpcl <- clusterExport(cl, c("out"))

#7. Run BRT function in parallel----
print("* Fitting models *")
mods <- parLapply(cl,
                  1:nrow(loop),
                  fun=brt_tune)

#CONCLUDE####

#1. Close clusters----
print("* Shutting down clusters *")
stopCluster(cl)

if(cc){ q() }