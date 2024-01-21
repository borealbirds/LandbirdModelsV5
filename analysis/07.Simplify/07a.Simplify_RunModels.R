# ---
# title: National Models 5.0 - tune models
# author: Elly Knight
# created: January 2, 2024
# ---

#NOTES################################


#PREAMBLE############################

#1. Load packages----
print("* Loading packages on master *")
library(tidyverse)
library(dismo)
library(parallel)
library(Matrix)

#2. Determine if testing and on local or cluster----
test <- FALSE
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
tmpcl <- clusterEvalQ(cl, library(Matrix))

#7. Load data----
print("* Loading data on master *")

#For running on cluster
if(cc){ load(file.path("06_NM5.0_data_tune.R")) }

#For testing on local
if(!cc){ load(file.path(root, "Data", "06_NM5.0_data_tune.R")) }

print("* Loading data on workers *")

if(cc){ tmpcl <- clusterEvalQ(cl, setwd("/home/ecknight/NationalModels")) }
tmpcl <- clusterExport(cl, c("bird", "offsets", "cov", "birdlist", "covlist", "bcrlist", "gridlist", "lrlist"))

#WRITE FUNCTION##########

brt_simplify <- function(i){
  
  t0 <- proc.time()
  
  #1. Read in the model----
  m.i <- readRDS("SOMETHING HERE")
  
  #2. Simplify model----
  s.i <- dismo::gbm.simplify(m.i)
  
  #3. Get performance metrics----
  out.i <- ????
  #TO DO: DECIDE WHAT TO SAVE####
  
  #4. Save----
  write.csv(out.i, file=file.path("results", paste0("ModelSimplification_", spp.i, "_", bcr.i, ".csv")))
  
  #5. Tidy up----
  rm(m.i, s.i, out.i)
  
}

#9. Export to clusters----
print("* Loading function on workers *")

tmpcl <- clusterExport(cl, c("brt_tune"))

#RUN MODELS###############

#1. Get model list----


#For testing
if(test) {loop <- loop[1:32,]}

print("* Loading model loop on workers *")
tmpcl <- clusterExport(cl, c("loop"))

#7. Run BRT function in parallel----
print("* Fitting models *")
mods <- parLapply(cl,
                  1:nrow(loop),
                  fun=brt_simplify)

#CONCLUDE####

#1. Close clusters----
print("* Shutting down clusters *")
stopCluster(cl)

if(cc){ q() }