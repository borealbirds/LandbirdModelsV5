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

#1. Determine if testing on local----
test <- FALSE

#2. Set root path for data on google drive (for local testing)----
root <- "G:/Shared drives/BAM_NationalModels/NationalModels5.0"

#3. Create nodes list----
print("* Creating nodes list *")

#For running on cluster
nodeslist <- unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))

#For testing on local
if(test){ nodeslist <- 2 }

print(nodeslist)

#4. Create and register clusters----
print("* Creating clusters *")
cl <- makePSOCKcluster(nodeslist, type="PSOCK")

#5. Load packages----
print("* Loading packages on master *")
library(tidyverse)
library(dismo)
library(parallel)

print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(dismo))

#6. Load data----
print("* Loading data on master *")

#For running on cluster
if(!test){ load(file.path("04_NM5.0_data_stratify.Rdata")) }

#For testing on local
if(test){ load(file.path(root, "Data", "04_NM5.0_data_stratify.Rdata")) }

print("* Loading data on workers *")

if(!test){ tmpcl <- clusterEvalQ(cl, setwd("/home/ecknight/NM")) }
tmpcl <- clusterExport(cl, c("DATA OBECTS HERE"))

#RUN MODELS###############

#1. Set grid search parameters----

#Learning rate
lr <- c(0.01, 0.005, 0.001, 0.0005, 0.0001)

#Interaction depth
id <- c(2, 3, 4)

#2. Get BCR & bird lists----
spp <- colnames(bird)[2:ncol(bird)]

bcrs <- rownames(bootstraps)

#For running on cluster
loop <- expand.grid(lr, id)

#For testing on local
if(test) {loop <- loop[,1:2]}

#2. Load BRT function----
load("00.BRTFunction.R")

#3. Run BRT function in parallel----


#SELECT PARAMETERS####



#CONCLUDE####

#1. Save----

#2. Close clusters----
print("* Shutting down clusters *")
stopCluster(cl)

if(!test){ q() }


