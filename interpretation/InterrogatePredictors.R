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


#2. Create a list of dataframes containing relative influence per covariate----
#   Every list element represents a bootstrap sample 

# Connect to BAM Drive and find bootstrap files 
root <- "G:/Shared drives/BAM_NationalModels/NationalModels5.0/"

# gbm_objs <- list.files(file.path(root, "output", "bootstraps"))

gbm_objs_test <- gbm_objs[1:3] # for testing


# Create a list of dataframes containing relative influence per covariate
covs <- list()
for(i in 1:length(gbm_objs_test)){
  
   # loads a `gbm` object named `b.i`
  load(file=file.path(root, "output", "bootstraps", gbm_objs_test[i]))
  
  # `summary(b.i)` is a `data.frame` with columns `var` and `rel.inf`
  covs[[i]] <- summary(b.i)
  
  # print progress
  cat(paste("\riteration", i))
    Sys.sleep(0.5)
}


#WRITE FUNCTION(S)#######################

#' for roxygen2:

#' @param species is `"all"` or a `character` with the FLBC indicating the species of interest.
#'
#' @param method not sure yet...
#'
#' @param type either 'ip' for an interpolating spline or 'smooth' for a
#' smoothing spline. Uses \code{stats::spline()} or \code{stats::smooth.spline()}, respectively,
#' for curve fitting and estimating the derivatives. Default is \code{type = 'smooth'}.
#' See: ?spline and ?smooth.spline for details.
#'
#' @return a `data.frame` to be passed to a plotting function...tbd
#'
#' @examples ...tbd

# what a user might want:
# rel.inf by species
# rel.inf by bcr


rel_inf_bcr <- function(species=c("all", ...), method){
  
  
  
  
}
                        
                        
# v.4.0 summary data format
#                      BCR                   varclass           x   total       prop
#                      10 Northern Rockies   climate 21.651179070   115 1.882711e-01
#                      10 Northern Rockies      time  2.806286740   115 2.440249e-02
#                      10 Northern Rockies  veglocal 18.907945015   115 1.644169e-01
#                      11 Prairie Potholes      topo  5.749344637   125 4.599476e-02
#                      11 Prairie Potholes   climate 25.935521062   125 2.074842e-01
#                      11 Prairie Potholes  veglocal 20.333884220   125 1.626711e-01


