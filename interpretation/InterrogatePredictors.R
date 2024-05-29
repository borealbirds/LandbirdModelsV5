# ---
# title: National Models 5.0 - covariate interpretation
# author: Mannfred Boehm
# created: May 22, 2024
# ---


#NOTES################################

# This script extracts covariate contributions to model predictions. 

# The outputs will be analysed for identifying covariates of importance across several metrics (e.g. BCR, ecology, etc.)



#PREAMBLE############################

#1. Attach packages----
library(gbm)
library(tidyverse)

#2. Determine if testing and on local or cluster----
test <- FALSE
cc <- TRUE

#3. Set nodes for local vs cluster----
if(cc){ nodes <- 32}
if(!cc | test){ nodes <- 4}

makePSOCKcluster



#WRITE FUNCTION#######################
# (this will eventually be a function for use in a package/shinyapp)



#1. Connect to BAM Drive and create an index of output attributes----
root <- "G:/Shared drives/BAM_NationalModels/NationalModels5.0/"

prediction_files <- list.files(root)

loop <- 
  prediction_files %>% 
  stringr::str_split_fixed(pattern="_", n=3) %>% 
  dplyr::as_tibble(.name_repair = "universal") %>% 
  magrittr::set_colnames(c("spp", "bcr", "boot")) %>% 
  dplyr::arrange(spp, bcr)



#2. Generate covariate lists----

covs <- list()
for(i in 1:nrow(loop)){
  
  # define spp x bcr x bootstrap permutation
  spp.i <- loop$spp[i]
  bcr.i <- loop$bcr[i]
  boot.i <- loop$boot[i]
  # lr.i <- loop$lr[i]
  # id.i <- 3
  
  # loads a `gbm` object named `b.i`
  load(file=file.path(root, "output", "bootstraps", paste0(spp.i, "_", bcr.i, "_", boot.i)))
  
  # `summary(b.i)` is a `data.frame` with columns `var` and `rel.inf`
  covs[[i]] <- summary(b.i)
  
  # print progress
  cat(paste("\riteration", i))
    Sys.sleep(0.5)
}


tmpcl <- clusterExport(cl, c("covs"))





# birdlist <- data.frame(bcr=bcrs)

# # create `loop` object (see: "05.Tune.R")
# bcr.spp <- birdlist %>% 
#   pivot_longer(AGOS:YTVI, names_to="spp", values_to="use") %>% 
#   dplyr::filter(use==TRUE) %>% 
#   dplyr::filter(spp %in% sppuse) %>% 
#   dplyr::select(-use)
# 
# 
# loop <- bcr.spp %>% 
#   anti_join(done) %>% 
#   anti_join(redo) %>% 
#   mutate(lr = 0.001) %>% 
#   rbind(redo %>% 
#           dplyr::select(spp, bcr, lr.next) %>% 
#           rename(lr = lr.next)) %>% 
#   arrange(spp, bcr)