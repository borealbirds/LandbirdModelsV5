# ---
# title: National Models 5.0 - tune models
# author: Elly Knight
# created: November 8, 2023
# ---

#NOTES################################

#This script finds the best learning rate for each combination of model subunt (BCR*country) and bird species.

#It starts with a learning rate of 0.001 and then drops up or down one order of magnitude if the model has less than 1000 or 10000 trees (the max). Interaction depth is held constant at 3 (i.e., 3-way interactions are biologicaly plausible). Models that cannot achieve 1000-10000 trees with a learning rate of at least 1e-06 are dropped.

#The other approach to tuning boosted regression trees is using a grid search of interaction depth & learning rate and picking the best parameters as the highest deviation explained. This approach was not used for the national models in the interest of efficiency.

#Tuning process is only run once using bootstrap #1, assuming that tuning parameters are relatively insensitive to bootstrap.

#This script is written to be run on the compute canada cluster in parallel. Messages are therefore printed in the output for each step in the script to help with debugging. Use the .sh object to run the script on the cluster.

#There are two object at the beginning of the script that can be set to TRUE if running tests. One controls whether you are running on a subset of model iterations (i.e., testing), and one controls whether you are running on your local machine or compute canada. Ideally you would:
#1. set test to true and cc to false (test on local)
#2. set cc to true (test on compute canada)
#3. set test to false (run full model set on compute canada)

#The script is set to inventory the models already run and remove them from the to-do list ("loop" object) every time the script is run. This means you can submit it as a [relatively] small job on compute canada instead of requesting enough resources for the entire list of models (which is a lot). Just make sure you keep the "results" folder with all the output in it until you are finished running everything, because this is what is being used to inventory the models already run.

#The steps for running on compute canada (for this and other scripts are):
#1. Transfer this script and your data object between local computer and compute canada's servers using Globus Connect.
#2. Create any folders you need for saving output into (as per your script).
#3. Log in to compute canada in the terminal with the ssh function.
#4. Figure out which modules you need to run your R packages (if making any changes) and load them.
#5. Load an instance of R in the test node and install the packages you need. Close it.
#6. Use the nano function to create a shell script that tells the slurm the resources you need, the modules you need, and the script to run. NOTE: You have to do this via the terminal or with a Linux machine. The slurm can't read the encoding if you create it on a windows machine and transfer it over via Globus.
#7. Use the cd function to navigate to the write working directory within copute canada's file servers (i.e., where you put your files with globus connect)
#8. Run your script with sbatch and the name of your sh file!
#9. Use squeue -u and your username to check on the status of your job.
#10. Rerun the script as many times as you need to work through all the required models (see note above).

#see medium.com/the-nature-of-food/how-to-run-your-r-script-with-compute-canada-c325c0ab2973 for the most straightforward tutorial I found for compute canada

#Things to change for next version: move eBird tagmethod wrangling to script 1. Move remove offset NA visits to stratification script.

#PREAMBLE############################

#1. Load packages----
print("* Loading packages on master *")
library(tidyverse)
library(dismo)
library(parallel)
library(Matrix)

#2. Determine if testing and on local or cluster----
test <- TRUE
cc <- FALSE

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

#10. Take the NA offsets out----
offsets <- offsets[is.infinite(offsets$ALFL)==FALSE,]
meth <- dplyr::filter(meth, id %in% offsets$id)
visit <- dplyr::filter(visit, id %in% offsets$id)

#10. Load data objects----
print("* Loading data on workers *")

if(cc){ tmpcl <- clusterEvalQ(cl, setwd("/home/ecknight/NationalModels")) }
tmpcl <- clusterExport(cl, c("bird", "offsets", "cov", "birdlist", "covlist", "bcrlist", "gridlist", "meth", "visit"))

#WRITE FUNCTION##########

brt_tune <- function(i){
  
  t0 <- proc.time()
  
  #1. Get model settings---
  bcr.i <- loop$bcr[i]
  spp.i <- loop$spp[i]
  boot.i <- 1
  lr.i <- loop$lr[i]
  id.i <- 3
  
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
  
  #4. Get covariates----
  covlist.i <- covlist %>% 
    dplyr::filter(bcr==bcr.i) %>% 
    pivot_longer(ERAMAP_1km:mTPI_1km, names_to="cov", values_to="use") %>% 
    dplyr::filter(use==TRUE)
  cov.i <- cov[cov$id %in% visit.i$id, colnames(cov) %in% covlist.i$cov]
  
  #5. Get PC vs SPT vs SPM vs eBird----
  meth.i <- meth[meth$id %in% visit.i$id, "method"]

  #6. Get year----
  year.i <- visit[visit$id %in% visit.i$id, "year"]
  
  #6. Put together data object----
  dat.i <- cbind(bird.i, year.i, meth.i, cov.i) %>% 
    rename(count = bird.i)
  
  #7. Get offsets----
  off.i <- offsets[offsets$id %in% visit.i$id, spp.i]
  
  #8. Run model----
  set.seed(i)
  m.i <- try(dismo::gbm.step(data=dat.i,
                         gbm.x=c(2:ncol(dat.i)),
                         gbm.y=1,
                         offset=off.i,
                         tree.complexity = id.i,
                         learning.rate = lr.i,
                         family="poisson"))
  
  trees.i <- ifelse(class(m.i)=="NULL", 0, m.i$n.trees)
  
  #rerun if lr too high (small sample size e.g., can3)
  #stop at a certain lr and remove that model from the to-do list
  while((trees.i < 1000 & lr.i >= lr.stop) |
        (trees.i==10000 & lr.i >= lr.stop)){
    
    if(trees.i < 1000){ lr.i <- lr.i/10}
    if(trees.i==10000){ lr.i <- lr.i*10}
    
    set.seed(i)
    m.i <- try(dismo::gbm.step(data=dat.i,
                               gbm.x=c(2:ncol(dat.i)),
                               gbm.y=1,
                               offset=off.i,
                               tree.complexity = id.i,
                               learning.rate = lr.i,
                               family="poisson"))
    
  }
  
  #9. Get performance metrics----
  if(class(m.i)=="NULL"){
    
    out.i <- loop[i,] %>% 
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
                       time = (proc.time()-t0)[3])) %>% 
      mutate(lr = lr.i)
    
    write.csv(out.i, file=file.path("output/tuning", paste0("ModelTuning_", spp.i, "_", bcr.i, "_", lr.i, ".csv")), row.names = FALSE)
    
  } else {
    
    out.i <- loop[i,] %>% 
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
                       time = (proc.time()-t0)[3])) %>% 
      mutate(lr = lr.i)
    
    #10. Save model----
    write.csv(out.i, file=file.path("output/tuning", paste0("ModelTuning_", spp.i, "_", bcr.i, "_", lr.i, ".csv")), row.names = FALSE)
    save(m.i, file=file.path("output/fullmodels", paste0(spp.i, "_", bcr.i, "_", lr.i, ".R")))
    
  }
  
  #16. Tidy up----
  rm(bcr.i, spp.i, boot.i, lr.i, id.i, visit.i, bird.i, covlist.i, cov.i, meth.i, dat.i, off.i, m.i)
  
}

#17. Export to clusters----
print("* Loading function on workers *")

tmpcl <- clusterExport(cl, c("brt_tune"))

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

#3. Get list of models already run----
files <- data.frame(path = list.files("output/tuning", pattern="*.csv", full.names=TRUE),
                    file = list.files("output/tuning", pattern="*.csv")) %>% 
  separate(file, into=c("step", "spp", "bcr", "lr"), sep="_", remove=FALSE) %>% 
  mutate(lr = as.numeric(str_sub(lr, -100, -5)))

#4. Set learning rate threshold for dropping a spp*bcr combo----
lr.stop <- 1e-06

print("* Loading stop threshold on workers *")
tmpcl <- clusterExport(cl, c("lr.stop"))

#5. Read in performance of those models----

if(nrow(files) > 0){
  
  perf <- map_dfr(read.csv, .x=files$path)
  
  #6. Determine which ones are done or won't run----
  done <- perf %>% 
    dplyr::filter((trees >= 1000 & trees < 10000) | lr<=lr.stop) %>% 
    dplyr::select(spp, bcr) %>% 
    unique()
  
  #7. Determine learning rate for those that need rerunning----
  redo <- perf %>% 
    anti_join(done) %>% 
    group_by(spp, bcr) %>% 
    mutate(last = ifelse(trees==10000, max(lr), min(lr))) %>% 
    dplyr::filter(lr == last) %>% 
    mutate(lr.next = case_when(trees == 10000 ~ lr*10,
                               trees < 1000 ~ lr/10)) %>% 
    ungroup()
  
  #8. Make dataframe of models to run----
  #Full combinations, take out done models and redo models, then add redo models back in
  loop <- bcr.spp %>% 
    anti_join(done) %>% 
    anti_join(redo) %>% 
    mutate(lr = 0.001) %>% 
    rbind(redo %>% 
            dplyr::select(spp, bcr, lr.next) %>% 
            rename(lr = lr.next)) %>% 
    arrange(spp, bcr)
  
}

#9. Make dataframe of models to run if there are no files yet----
if(nrow(files)==0){
  loop <- bcr.spp %>% 
    mutate(lr = 0.001) %>% 
    arrange(spp, bcr)
}

#For testing
if(test) {loop <- loop[1:nodes,]}

print("* Loading model loop on workers *")
tmpcl <- clusterExport(cl, c("loop"))

#10. Run BRT function in parallel----
print("* Fitting models *")
mods <- parLapply(cl,
                  X=1:nrow(loop),
                  fun=brt_tune)

#CONCLUDE####

#1. Close clusters----
print("* Shutting down clusters *")
stopCluster(cl)

if(cc){ q() }