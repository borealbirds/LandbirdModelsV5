# ---
# title: National Models 5.0 - tune models
# author: Elly Knight
# created: November 8, 2023
# ---

#NOTES################################

#This script finds the best learning rate for each combination of model subunt (BCR*country) and bird species.

#The species list can be all the species available in the dataset or can be determined with a custom list.

#It starts with a learning rate of 0.001 and then drops up or down one order of magnitude if the model has less than 1000 or 10000 trees (the max). Interaction depth is held constant at 3 (i.e., 3-way interactions are biologicaly plausible). Models that cannot achieve 1000-10000 trees with a learning rate of at least 1e-06 are dropped.

#The other approach to tuning boosted regression trees is using a grid search of interaction depth & learning rate and picking the best parameters as the highest deviation explained. This approach was not used for the national models in the interest of efficiency.

#Tuning process is only run once using bootstrap #1, assuming that tuning parameters are relatively insensitive to bootstrap.

#This script is written to be run on the compute canada cluster in parallel. Messages are therefore printed in the output for each step in the script to help with debugging. Use the .sh object to run the script on the cluster.
#There are two object at the beginning of the script that can be set to TRUE if running tests. One controls whether you are running on a subset of model iterations (i.e., testing), and one controls whether you are running on your local machine or compute canada. Ideally you would:
#1. set test to true and cc to false (test on local)
#2. set cc to true (test on compute canada)
#3. set test to false (run full model set on compute canada)

#The script is set to inventory the models already run and remove them from the to-do list ("loop" object) every time the script is run. Make sure you keep the "results" folder with all the output in it until you are finished running everything, because this is what is being used to inventory the models already run. Note this script can only one on one 48-core node at a time due to the inability to distribute tasks across multiple nodes due to ssh limitations. If you want to run more than that at a time, batch your species list into chunks and run different chunks of the species list concurrently.

#The steps for running on compute canada (for this and other scripts are):
#1. Transfer this script and your data object between local computer and compute canada's servers using Globus Connect.
#2. Create any folders you need for saving output into (as per your script). You'll want to work in the scratch space, given the volume of output we will produce.
#3. Log in to compute canada in the terminal with the ssh function.
#4. Figure out which modules you need to run your R packages (if making any changes) and load them.
#5. Load an instance of R in the test node and install the packages you need. Close it.
#6. Use the nano function to create a shell script that tells the slurm the resources you need, the modules you need, and the script to run. NOTE: You have to do this via the terminal or with a Linux machine. The slurm can't read the encoding if you create it on a windows machine and transfer it over via Globus.
#7. Use the cd function to navigate to the write working directory within copute canada's file servers (i.e., where you put your files with globus connect)
#8. Run your script with sbatch and the name of your sh file!
#9. Use squeue -u and your username to check on the status of your job.
#10. Rerun the script as many times as you need to work through all the required models (see note above).

#see medium.com/the-nature-of-food/how-to-run-your-r-script-with-compute-canada-c325c0ab2973 for the most straightforward tutorial I found for compute canada

#BAM has previously used the Graham cluster; however, we hit major obstacles with innaccessibility of the scratch space during V5 and moved our work to Cedar.

#Although this and the subsequent scripts for the cluster are written to run on more than one node at once, it is typically fastest to request one node for 24 hours at a time. Instead, we use a lookup table (e.g., sppuse below) with subgroupings to run more than group per modelling step at a time, each on a single node.

#PREAMBLE############################

#1. Load packages----
print("* Loading packages on master *")
library(tidyverse)
library(dismo)
library(parallel)
library(Matrix)

#2. Set species subset ----
set <- c(1:5)

#3. Determine if testing and on local or cluster----
test <- FALSE
cc <- FALSE

#4. Set cores for local vs cluster----
if(cc){ cores <- 24 }
if(!cc | test){ cores <- 1}

#5. Create and register clusters----
print("* Creating clusters *")
print(table(cores))
cl <- makePSOCKcluster(cores, type="PSOCK")

#6. Set root path----
print("* Setting root file path *")
if(cc){root <- "/scratch/ecknight/NationalModels"}
if(!cc){root <- "G:/Shared drives/BAM_NationalModels5"}

tmpcl <- clusterExport(cl, c("root"))

#7. Load packages on clusters----
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(dismo))
tmpcl <- clusterEvalQ(cl, library(tidyverse))
tmpcl <- clusterEvalQ(cl, library(Matrix))

#8. Load data package----
print("* Loading data on master *")

load(file.path(root, "data", "04_NM5.0_data_stratify.Rdata"))

#9. Set species subset----
sppuse <- read.csv(file.path(root, "data", "priority_spp_with_model_performance.csv")) |> 
  rename(spp = species_code) |> 
  dplyr::select(spp, rerun) |> 
  dplyr::filter(rerun %in% set)

#10. Load data objects----
print("* Loading data on workers *")

tmpcl <- clusterExport(cl, c("bird", "offsets", "cov", "covlist", "bcrlist", "bootlist", "visit"))

#WRITE FUNCTION##########

brt_tune <- function(i){
  
  t0 <- proc.time()
  
  #1. Get model settings----
  bcr.i <- loop$bcr[i]
  spp.i <- loop$spp[i]
  boot.i <- 1
  lr.i <- loop$lr[i]
  id.i <- 3
  
  #2. Get visits to include----
  #use bootstrap #1
  visit.i <- bcrlist[bcrlist[,bcr.i]>0, c("id", bcr.i)] |> 
    dplyr::filter(id %in% bootlist$b1)
  
  #3. Get response data (bird data)----
  bird.i <- bird[as.character(visit.i$id), spp.i]
  
  #4. Get covariates----
  covlist.i <- covlist |> 
    dplyr::filter(bcr==bcr.i) |> 
    pivot_longer(-bcr, names_to="cov", values_to="use") |> 
    dplyr::filter(use==TRUE)
  cov.i <- cov[cov$id %in% visit.i$id, colnames(cov) %in% covlist.i$cov]
  
  #5. Get PC vs SPT vs SPM vs eBird----
  meth.i <- visit[visit$id %in% visit.i$id, "method"]

  #6. Get year----
  year.i <- visit[visit$id %in% visit.i$id, "year"]
  
  #7. Put together data object----
  dat.i <- cbind(bird.i, year.i, meth.i, cov.i) |> 
    rename(count = bird.i)
  
  #8. Get offsets----
  off.i <- offsets[offsets$id %in% visit.i$id, spp.i]
  
  #9. Clean up to save space----
  rm(visit.i, bird.i, year.i, meth.i, cov.i)
  
  #10. Run model----
  set.seed(1234)
  m.i <- try(dismo::gbm.step(data=dat.i,
                         gbm.x=c(2:ncol(dat.i)),
                         gbm.y=1,
                         offset=off.i,
                         tree.complexity = id.i,
                         learning.rate = lr.i,
                         family="poisson"))
  
  trees.i <- ifelse(class(m.i)%in% c("NULL", "try-error"), 0, m.i$n.trees)
  
  #rerun if lr too high (small sample size e.g., can3)
  #stop at a certain lr and remove that model from the to-do list
  while((trees.i < 1000 & lr.i > lr.min) |
        (trees.i==10000 & lr.i < lr.max)){
    
    if(trees.i < 1000){ lr.i <- lr.i/10}
    if(trees.i==10000 & lr.i < lr.max){ lr.i <- lr.i*10}
    
    set.seed(1234)
    m.i <- try(dismo::gbm.step(data=dat.i,
                               gbm.x=c(2:ncol(dat.i)),
                               gbm.y=1,
                               offset=off.i,
                               tree.complexity = id.i,
                               learning.rate = lr.i,
                               family="poisson"))
    
    # divide lr by 2 if stuck in eternal loop of 10000 and < 1000 
    if(trees.i==10000 & m.i$n.trees < 1000){
      
      lr.i <- lr.i/2
      
      set.seed(1234)
      m.i <- try(dismo::gbm.step(data=dat.i,
                                 gbm.x=c(2:ncol(dat.i)),
                                 gbm.y=1,
                                 offset=off.i,
                                 tree.complexity = id.i,
                                 learning.rate = lr.i,
                                 family="poisson"))
      
    }
    
    trees.i <- ifelse(class(m.i)%in% c("NULL", "try-error"), 0, m.i$n.trees)
    
  }
  
  #11. Get performance metrics----
  if(class(m.i) %in% c("NULL", "try-error")){
    
    out.i <- loop[i,] |> 
      cbind(data.frame(trees = NA,
                       deviance.mean = NA,
                       deviance.se = NA,
                       null = NA,
                       resid = NA,
                       correlation = NA,
                       correlation.mean = NA,
                       correlation.se = NA,
                       n = nrow(dat.i),
                       ncount = nrow(dplyr::filter(dat.i, count > 0)),
                       time = (proc.time()-t0)[3])) |> 
      mutate(lr = lr.i)
    
    if(!(file.exists(file.path(root, "output", "05_tuning", spp.i)))){
      dir.create(file.path(root, "output", "05_tuning", spp.i))
    }
    
    write.csv(out.i, file=file.path(file.path(root, "output", "05_tuning", spp.i, paste0("ModelTuning_", spp.i, "_", bcr.i, ".csv"))), row.names = FALSE)
    
  } else {
    
    out.i <- loop[i,] |> 
      cbind(data.frame(trees = m.i$n.trees,
                       deviance.mean = m.i$cv.statistics$deviance.mean,
                       deviance.se = m.i$cv.statistics$deviance.se,
                       null = m.i$self.statistics$mean.null,
                       resid = m.i$self.statistics$mean.resid,
                       correlation = m.i$self.statistics$correlation,
                       correlation.mean = m.i$cv.statistics$correlation.mean,
                       correlation.se = m.i$cv.statistics$correlation.se,
                       n = nrow(dat.i),
                       ncount = nrow(dplyr::filter(dat.i, count > 0)),
                       time = (proc.time()-t0)[3])) |> 
      mutate(lr = lr.i)
    
    #12. Save model----
    if(!(file.exists(file.path(root, "output", "05_fullmodels", spp.i)))){
      dir.create(file.path(root, "output", "05_fullmodels", spp.i))
    }
    
    if(!(file.exists(file.path(root, "output", "05_tuning", spp.i)))){
      dir.create(file.path(root, "output", "05_tuning", spp.i))
    }
    
    write.csv(out.i, file=file.path(file.path(root, "output", "05_tuning", spp.i, paste0("ModelTuning_", spp.i, "_", bcr.i, ".csv"))), row.names = FALSE)
    save(m.i, file=file.path(root, "output", "05_fullmodels", spp.i, paste0(spp.i, "_", bcr.i, ".Rdata")))
    
  }
  
}

#13. Export to clusters----
print("* Loading function on workers *")

tmpcl <- clusterExport(cl, c("brt_tune"))

#RUN MODELS###############

#1. Get BCR & bird combo list----
#filter by list of desired species
bcr.spp <- birdlist |> 
  pivot_longer(-bcr, names_to="spp", values_to="use") |> 
  dplyr::filter(use==TRUE) |> 
  dplyr::select(-use)

#2. Reformat covariate list----
bcr.cov <- covlist |> 
  pivot_longer(-bcr, names_to="cov", values_to="use") |> 
  dplyr::filter(use==TRUE)

print("* Loading covariate list on workers *")
tmpcl <- clusterExport(cl, c("bcr.cov"))

#3. Get list of models already run----
files <- data.frame(path = list.files(file.path(root, "output", "05_tuning"), pattern="*.csv", full.names=TRUE, recursive = TRUE),
                    file = list.files(file.path(root, "output", "05_tuning"), pattern="*.csv", recursive = TRUE)) |> 
  separate(file, into=c("step", "spp", "bcr"), sep="_", remove=FALSE) |> 
  mutate(bcr = str_sub(bcr, -100, -5)) |> 
  dplyr::filter(!is.na(bcr))

#4. Set learning rate threshold for dropping a spp*bcr combo----
lr.min <- 1e-5
lr.max <- 0.01

print("* Loading stop threshold on workers *")
tmpcl <- clusterExport(cl, c("lr.min", "lr.max"))

#5. Read in performance of those models----

if(nrow(files) > 0){
  
  perf <- map_dfr(read.csv, .x=files$path)
  
  #6. Determine those that won't run----
  norun <- dplyr::filter(perf, is.na(trees))
  
  write.csv(norun, file.path(root, "output", "05_tuning", "SpeciesBCRCombos_NotTuned.csv"), row.names = FALSE)
  
  #7. Make dataframe of models to run----
  #Full combinations, take out done models
  loop <- bcr.spp |> 
    anti_join(perf) |> 
    mutate(lr = 0.001) |> 
    arrange(spp, bcr) |> 
    inner_join(sppuse)
  
}

#8. Make dataframe of models to run if there are no files yet----
if(nrow(files)==0){
  loop <- bcr.spp |> 
    mutate(lr = 0.001) |> 
    arrange(spp, bcr) |> 
    inner_join(sppuse)
}

#For testing
if(test) {loop <- loop[1:cores,]}

cat("* Loading model loop of", nrow(loop), "rows on workers *")
tmpcl <- clusterExport(cl, c("loop"))

#9. Run BRT function in parallel----
print("* Fitting models *")
mods <- parLapply(cl,
                  X=1:nrow(loop),
                  fun=brt_tune)

#CONCLUDE####

#1. Close clusters----
print("* Shutting down clusters *")
stopCluster(cl)

if(cc){ q() }
