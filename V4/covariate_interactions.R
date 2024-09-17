# ---
# title: National Models 4.0 - estimating and ranking 2-way covariate interactions for CONWA and CAWA
# author: Mannfred Boehm
# created: September 13, 2024
# ---

#NOTES################################

# This script creates a list of dataframes containing relative influence per covariate
# Every list element represents a bootstrap sample 




#1. attach packages----
print("* attaching packages on master *")
library(gbm)
library(parallel)
library(tidyverse)


#2. define local or cluster
test <- FALSE
cc <- TRUE


#3. set nodes for local vs cluster----
if(cc){ nodes <- 16}
if(!cc | test){ nodes <- 4}


#4. create and register clusters----
# `makePSOCKcluster` creates a set of copies of R running in parallel that communicate over "sockets". 
print("* creating clusters *")
cl <- makePSOCKcluster(nodes, type="PSOCK")


#5. set root path----
print("* setting root file path *")

if(!cc){root <- "G:/Shared drives/BAM_NationalModels4/NationalModels4.0"}
if(cc){root <- "/home/mannfred/scratch"}

tmpcl <- clusterExport(cl, c("root"))



#6. attach packages on clusters----
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(gbm))
tmpcl <- clusterEvalQ(cl, library(tidyverse))



#7. index gbm objects for CONW and CAWA----
print("* indexing gbm objects for CONW and CAWA *")

# gbm objects for CAWA and CONWA have been transfered to the symbolic link file "scratch"
# 16 folders x 32 bootstraps = 512 gbm models per species
gbm_conw <- list.files(file.path(root, "CONW"), 
                            pattern = "^gnmboot-CONW-BCR_([0-9]|[1-9][0-9])-[0-9]", 
                            recursive=TRUE, full.names = TRUE)

gbm_cawa <- list.files(file.path(root, "CAWA"), 
                       pattern = "^gnmboot-CAWA-BCR_([0-9]|[1-9][0-9])-[0-9]", 
                       recursive=TRUE, full.names=TRUE)

gbm_objs <- c(gbm_conw, gbm_cawa)


#8. export the necessary variables and functions to the cluster----
print("* exporting gbm_objs and functions to cluster *")
clusterExport(cl, c("gbm_objs", "plot.gbm"))


#9. create function to process each gbm object (parallelized)----
process_gbm <- function(obj_path) {
  
  # load a bootstrap replicate (gbm object)
  load(file.path(obj_path))  
  
  # find evaluation points (`data.frame`) for every covariate permutation of degree 2 (indexed by i,j) 
  pts <- list()
  n <- length(out$var.names)  # retrieve the number of variables
  interaction_index <- 1 # starts at 1 and increases for every covariate interaction computed. Resets at one when moving to the next bootstrap model.
  
  # end at n-1 to avoid finding the interaction of variable n x variable n
  for (i in 1:(n-1)) {
    
    # start at i+1 to avoid finding the interaction of variable 1 x variable 1
    for (j in (i+1):n) {
      
      # `continuous.resolution` is defaulted at 100: with two covariates this produces a dataframe of length=100*100 (10000 rows)
      # by lowering to 25, we get 625 rows, which should still be enough resolution to find local maximums
      # we then discard the grid (too much data) and keep the mean and std. dev. of the response for covariates i,j
      grid_ij <- plot.gbm(x = out, return.grid = TRUE, i.var = c(i, j), continuous.resolution = 25, type = "response")
      
      # calculate mean and sd of response
      pts[[interaction_index]] <- matrix(data = c(mean(grid_ij$y), sd(grid_ij$y)), ncol = 2, nrow = 1)
      colnames(pts[[interaction_index]]) <- c("y_mean", "y_sd")
      
      # label the interaction
      names(pts)[interaction_index] <- paste(out$var.names[i], out$var.names[j], sep = ".")  # Label the interaction
      interaction_index <- interaction_index + 1
    }
  }
  
  return(pts)  # return the list of interaction points
}

#10. run the function in parallel----
print("* running `process_gbm` in parallel *")
boot_pts_i2 <- parLapply(cl, gbm_objs, process_gbm)


#11. stop the cluster----
print("* stopping cluster *")
stopCluster(cl)


#12. 
print("* saving RDS file *")
saveRDS(boot_pts_i2, file=file.path(root, "boot_pts_i2.rds"))

print("* completed :-) *")
