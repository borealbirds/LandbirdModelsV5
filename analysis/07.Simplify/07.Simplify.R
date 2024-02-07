# ---
# title: National Models 5.0 - simplify models
# author: Elly Knight
# created: February 5, 2024
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

#3. Set nodes for local vs cluster----
if(cc){ nodes <- 32}
if(!cc){ nodes <- 2}

#4. Set root path for data on google drive (for local testing)----
root <- "G:/Shared drives/BAM_NationalModels/NationalModels5.0"

#5. Create nodes list----
print("* Creating nodes list *")

#For running on cluster
nodeslist <- unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))

#For testing on local
if(test){ nodeslist <- 32 }

print(nodeslist)

#6. Create and register clusters----
print("* Creating clusters *")
cl <- makePSOCKcluster(nodeslist, type="PSOCK")

#7. Load packages on clusters----
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(dismo))
tmpcl <- clusterEvalQ(cl, library(tidyverse))
tmpcl <- clusterEvalQ(cl, library(Matrix))

#WRITE FUNCTION##########

brt_simplify <- function(i){
  
  t0 <- proc.time()
  
  #1. Read in the model----
  load(paste0("output/fullmodels/", use$spp[i], "_", use$bcr[i], "_", use$lr[i], ".R"))
  
  #2. Simplify model----
  s.i <- dismo::gbm.simplify(m.i)
  
  #3. Get performance metrics----
  out.i <- s.i[["final.drops"]] %>% 
    mutate(spp = use$spp[i],
           bcr = use$bcr[i])
  
  #4. Save----
  write.csv(out.i, file=file.path("output/simplification", paste0("ModelSimplification_", spp.i, "_", bcr.i, ".csv")))
  save(s.i, file=file.path("output/simplifiedmodels", paste0(use$spp[i], "_", use$bcr[i], ".R")))
  
  #5. Tidy up----
  rm(m.i, s.i, out.i)
  
}

#9. Export to clusters----
print("* Loading function on workers *")

tmpcl <- clusterExport(cl, c("brt_simplify"))

#RUN MODELS###############

#1. Get list of models already run----
files <- data.frame(path = list.files("output/tuning", pattern="*.csv", full.names=TRUE),
                    file = list.files("output/tuning", pattern="*.csv")) %>% 
  separate(file, into=c("step", "spp", "bcr", "lr"), sep="_", remove=FALSE) %>% 
  mutate(lr = as.numeric(str_sub(lr, -100, -5)))

#2. Read in performance of those models----
#Take out the mutate after running next time
perf <- map_dfr(read.csv, .x=files$path) %>%
  dplyr::select(-lr, -lr.1) %>% 
  cbind(data.frame(lr = files$lr))

#3. Pick the best model for each spp*bcr----
use <- perf %>% 
  dplyr::filter(trees!=10000,
                trees >= 1000) 

#4. Check you have one model for each spp*bcr----
mods <- perf %>% 
  dplyr::select(spp, bcr) %>% 
  unique()

missing <- anti_join(mods, use) %>% 
  left_join(perf)

#For testing
if(test) {use <- use[1:nodes,]}

print("* Loading model loop on workers *")
tmpcl <- clusterExport(cl, c("use"))

#7. Run BRT function in parallel----
print("* Simplifying models *")
mods <- parLapply(cl,
                  1:nrow(use),
                  fun=brt_simplify)

#CONCLUDE####

#1. Close clusters----
print("* Shutting down clusters *")
stopCluster(cl)

if(cc){ q() }
