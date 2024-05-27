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
library(tidyverse)



#WRITE FUNCTION#######################
# (this will eventually be a function for use in a package/shinyapp)



#1. Connect to BAM Drive and create an index of output attributes----
root <- "G:/Shared drives/BAM_NationalModels/NationalModels5.0/output/bootstraps"

prediction_files <- list.files(root)

loop <- 
  prediction_files %>% 
  stringr::str_split_fixed(pattern="_", n=3) %>% 
  dplyr::as_tibble() %>% 
  magrittr::set_colnames(c("spp", "bcr", "boot")) %>% 
  dplyr::arrange(spp, bcr)



#2. Generate covariate lists (remove covs that explain < 0.01 % of deviance)----
# Based on code from "06.Bootstrap.R")

covs <- list()
for(i in 1:nrow(loop)){
  
  # define spp x bcr x bootstrap permutation
  bcr.i <- loop$bcr[i]
  spp.i <- loop$spp[i]
  boot.i <- loop$boot[i]
  # lr.i <- loop$lr[i]
  # id.i <- 3
  
  # loads a `gbm` object named `m.i`. 
  # `m.i$contributions` is a `data.frame` with columns `var` and `rel.inf`
  load(file=file.path(root, "output", "bootstraps", paste0(loop$spp[i], "_", loop$bcr[i], "_", loop$lr[i], ".R")))
  
  #make sure to retain method and year
  covs[[i]] <- m.i[["contributions"]] %>% 
    dplyr::filter(rel.inf >= 0.1 |
                    var %in% c("year", "method"))
  
}

print("* Loading new covariate lists *")
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