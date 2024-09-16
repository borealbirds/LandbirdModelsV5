# ---
# title: National Models 4.0 - estimating and ranking 2-way covariate interactions for CONWA and CAWA
# author: Mannfred Boehm
# created: September 13, 2024
# ---

#NOTES################################

# This script extracts covariate contributions to model predictions, and also synthesizes various trait databases.

# The outputs will be analysed for identifying covariates of importance across BCRs



#1. Attach packages----
library(gbm)
library(RcppAlgos)
library(tidyverse)



#2. Create a list of dataframes containing relative influence per covariate----
#   Every list element represents a bootstrap sample 

# connect to BAM Drive and find bootstrap files 
root <- "G:/Shared drives/BAM_NationalModels4/NationalModels4.0"

# gbm objects for CAWA and CONWA stored on the BAM drive
# 16 folders x 32 bootstraps = 512 gbm models
gbm_conw <- list.files(file.path(root, "Feb2020", "out", "boot", "CONW"), 
                            pattern = "^gnmboot-CONW-BCR_([0-9]|[1-9][0-9])-[0-9]", 
                            recursive=TRUE, full.names = TRUE)

gbm_cawa <- list.files(file.path(root, "Feb2020", "out", "boot", "CAWA"), 
                       pattern = "^gnmboot-CAWA-BCR_([0-9]|[1-9][0-9])-[0-9]", 
                       recursive=TRUE, full.names=TRUE)

gbm_objs <- c(gbm_conw, gbm_cawa)




#3 create an index from `gbm_objs` containing the species (FLBC), BCR, and bootstrap replicate----
# append binomial names to FLBC
sample_id <- 
  gbm_objs |> 
  stringr::str_split_fixed(pattern="-", n=4) |> 
  gsub("\\.RData", "", x = _) |>
  tibble::as_tibble() |> 
  dplyr::select(2:4) |> 
  magrittr::set_colnames(c("spp", "bcr", "boot"))



#4. Create a list of two-way covariate interactions by bcr x species x bootstrap replicate----

# this loop creates a list of lists of `data.frame`s. 
# each top level element is a bootstrap replicate, each second-level element represents a 2-way covariate interaction for that bootstrap.
# each covariate interaction has a corresponding `data.frame`, with columns 1 and 2 being the covariate domains and column 3 being the response (density)
# NOTE: In computing interactions involving discrete variables (e.g. `MODISLCC_1km`) the resultant plot will differ compared to interactions between two continuous variables (heatmap).
# `plot.gbm()` produces faceted line plots instead of a heatmap. Each line plot corresponds to one of the 17 levels of `MODISLCC_1km`, showing how `MODISLCC_1km` affects the response variable at each level.


boot_pts_i2 <- list()
for (q in 1:length(gbm_objs)){ 
  
  # load a bootstrap replicate 
  load(file.path(gbm_objs[q]))
  
  # find evaluation points (`data.frame`) for every covariate permutation of degree 2 (indexed by i,j) 
  pts <- list()
  n <- length(out$var.names) # get the number of variables
  interaction_index <- 1 # starts at 1 and increases for every covariate interaction computed. Resets at one when moving to the next bootstrap model.
  
  # end at n-1 to avoid finding the interaction of variable n x variable n
  for (i in 1:(n-1)) {
    
    # start at i+1 to avoid finding the interaction of variable 1 x variable 1
    for (j in (i+1):n) {
      
      # `continuous.resolution` is defaulted at 100: with two covariates this produces a dataframe of length=100*100 (10000 rows)
      # by lowering to 25, we get 625 rows, which should still be enough resolution to find local maximums
      # we then discard the grid (too much data) and keep the mean and std. dev. of the response for covariates i,j
      grid_ij <- plot.gbm(x = out, return.grid = TRUE, i.var = c(i, j), continuous.resolution=25, type="response")  
      
      pts[[interaction_index]] <- matrix(data=c(mean(grid_ij$y), sd(grid_ij$y)), ncol=2, nrow=1)
      colnames(pts[[interaction_index]]) <- c("y_mean", "y_sd")
      
      names(pts)[interaction_index] <- paste(out$var.names[i], out$var.names[j], sep = ".") # label the interaction
      
      interaction_index <- interaction_index + 1
      
    }
  }
  
  boot_pts_i2[[q]] <- pts
  
  # print progress
  cat(paste("\riteration", q))
  Sys.sleep(0.001)
}

saveRDS(boot_pts_i2, file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boot_pts_i2.rds")

