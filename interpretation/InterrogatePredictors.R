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
print("* Attaching packages on master *")
library(tidyverse)



#WRITE FUNCTION#######################




#2. Generate covariate lists (remove covs that explain < 0.01 % of deviance)----
# (Based on code from "06.Bootstrap.R")

covsnew <- list()
for(i in 1:nrow(loop)){
  
  load(file=file.path(root, "output", "fullmodels", paste0(loop$spp[i], "_", loop$bcr[i], "_", loop$lr[i], ".R")))
  
  #make sure to retain method and year
  covsnew[[i]] <- m.i[["contributions"]] %>% 
    dplyr::filter(rel.inf >= 0.1 |
                    var %in% c("year", "method"))
  
}

print("* Loading new covariate lists *")
tmpcl <- clusterExport(cl, c("covsnew"))