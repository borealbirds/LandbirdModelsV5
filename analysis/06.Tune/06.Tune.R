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

#PREAMBLE############################

#1. Load packages----
print("* Loading packages on master *")
library(tidyverse)
library(dismo)
library(parallel)

#2. Determine if testing on local----
test <- TRUE

#3. Set root path for data on google drive (for local testing)----
root <- "G:/Shared drives/BAM_NationalModels/NationalModels5.0"

#4. Create nodes list----
print("* Creating nodes list *")

#For running on cluster
nodeslist <- unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))

#For testing on local
if(test){ nodeslist <- 2 }

print(nodeslist)

#5. Create and register clusters----
print("* Creating clusters *")
cl <- makePSOCKcluster(nodeslist, type="PSOCK")

#6. Load packages on clusters----
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(dismo))

#7. Load data----
print("* Loading data on master *")

#For running on cluster
if(!test){ load(file.path("04_NM5.0_data_stratify.Rdata")) }

#For testing on local
if(test){ load(file.path(root, "Data", "04_NM5.0_data_stratify.Rdata")) }

print("* Loading data on workers *")

if(!test){ tmpcl <- clusterEvalQ(cl, setwd("/home/ecknight/NM")) }
tmpcl <- clusterExport(cl, c("visit", "bird", "offsets", "bootstraps", "birdlist"))

#RUN MODELS###############

#1. Set grid search parameters----

#Learning rate
lr <- c(0.01, 0.005, 0.001, 0.0005, 0.0001)

#Interaction depth
id <- c(2, 3, 4)

#2. Set number of bootstraps----
boot <- 1

#3. Get BCR & bird combo list----
bcr.spp <- birdlist %>% 
  rename(bcr=id) %>% 
  pivot_longer(ABDU:YTVI, names_to="spp", values_to="use") %>% 
  dplyr::filter(use==TRUE)

#4. Make dataframe of models to run----

#For running on cluster
loop <- expand.grid(lr=lr, id=id, boot=boot) %>% 
  merge(bcr.spp) %>% 
  dplyr::select(-use)
tmpcl <- clusterExport(cl, c("loop"))

#For testing on local
if(test) {loop <- loop[1:2,]}

#5. Load BRT function----
load("00.BRTFunction.R")

#6. Run BRT function in parallel----
print("* Fitting models *")
mods <- parLapply(cl,
                  1:nrow(loop),
                  fun=brt.tune)


#SELECT PARAMETERS####



#CONCLUDE####

#1. Save----

#2. Close clusters----
print("* Shutting down clusters *")
stopCluster(cl)

if(!test){ q() }


