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
tmpcl <- clusterEvalQ(cl, library(terra))

#8. Load the stacks----


#tpmcl <- clusterEvalQ()



#WRITE FUNCTION##########

brt_predict <- function(i){
  
  t0 <- proc.time()
  
  #1. Get model settings---
  bcr.i <- loop$bcr[i]
  spp.i <- loop$spp[i]
  boot.i <- loop$boot[i]
  
  #2. Load model----
  
  #3. Load raster stack----
  
  #4. Predict----
  
  #5. Save----
  
  
  
  
  
  
}

#RUN MODELS#########

#1. Set desired years----
#years <- seq(1985, 2020, 5)
years <- 2020

#2. Get list of models that are tuned----
tuned <- data.frame(path = list.files("output/bootstraps", pattern="*.csv", full.names=TRUE),
                    file = list.files("output/bootstraps", pattern="*.csv")) %>% 
  separate(file, into=c("step", "spp", "bcr", "lr", "boot"), sep="_", remove=FALSE)

#5. Create to do list----
#Sort longest to shortest duration to get the big models going first
todo <- perf %>% 
  dplyr::select(bcr, spp, lr, trees, time) %>% 
  expand_grid(boot=c(1:boots)) %>% 
  arrange(-time)

#6. Determine which are already done----
done <- data.frame(path = list.files("output/bootstraps", pattern="*.csv", full.names=TRUE),
                   file = list.files("output/bootstraps", pattern="*.csv")) %>% 
  separate(file, into=c("step", "spp", "bcr", "boot"), sep="_", remove=FALSE) %>% 
  mutate(boot = as.numeric(str_sub(boot, -100, -5)))

#7. Create final to do list----
if(nrow(done) > 0){
  
  loop <- todo %>% 
    anti_join(done)
  
} else { loop <- todo }

#For testing - take the shortest duration models
if(test) {loop <- arrange(loop, time)[1:2,]}

print("* Loading model loop on workers *")
tmpcl <- clusterExport(cl, c("loop"))

#8. Update the covariate lists (remove covs that explain < 0.01 % of deviance)----
covsnew <- list()
for(i in 1:nrow(loop)){
  
  load(file=file.path("output/fullmodels", paste0(loop$spp[i], "_", loop$bcr[i], "_", loop$lr[i], ".R")))
  
  covsnew[[i]] <- m.i[["contributions"]] %>% 
    dplyr::filter(rel.inf >= 0.1)
  
}

print("* Loading new covariate lists *")
tmpcl <- clusterExport(cl, c("covsnew"))

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