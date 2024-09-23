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
test <- TRUE
cc <- FALSE


#3. set nodes for local vs cluster----
if(cc){ nodes <- 32}
if(!cc | test){ nodes <- 4}


#4. create and register clusters----
# creates 16 copies of R running in parallel via 16 nodes, on one if the cluster's sockets (processors). 
# Belgua has ~965 nodes
print("* creating clusters *")
cl <- makePSOCKcluster(nodes, type="PSOCK")

# print number of nodes and host name for confirmation
cl


#5. set root path----
print("* setting root file path *")

if(!cc){root <- "G:/Shared drives/BAM_NationalModels4/NationalModels4.0/Feb2020/out/boot"}
if(cc){root <- "/home/mannfred/scratch"}

tmpcl <- clusterExport(cl, c("root"))



#6. attach packages on clusters----
# `clusterEvalQ` evaluates a literal expression on each cluster node. 
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(gbm))
tmpcl <- clusterEvalQ(cl, library(tidyverse))



#7. index gbm objects for CONW and CAWA----
print("* indexing gbm objects for CONW and CAWA *")

# gbm objects for CAWA and CONWA have been transfered to the symbolic link file "scratch"
# 16 folders x 32 bootstraps = 512 gbm models per species
# gbm_conw <- list.files(file.path(root, "CONW"), 
#                             pattern = "^gnmboot-CONW-BCR_([0-9]|[1-9][0-9])-[0-9]", 
#                             recursive=TRUE, full.names = TRUE)
# 
# gbm_cawa <- list.files(file.path(root, "CAWA"), 
#                        pattern = "^gnmboot-CAWA-BCR_([0-9]|[1-9][0-9])-[0-9]", 
#                        recursive=TRUE, full.names=TRUE)

gbm_conw_bcr6 <- list.files(file.path(root, "CONW", "BCR_60"), full.names = TRUE)

gbm_objs <- gbm_conw_bcr6[1:3] 


#8. create function to process each gbm object (parallelized)----
process_gbm <- function(obj_path) {
  
  # load a bootstrap replicate (gbm object)
  try(load(file.path(obj_path)))
  if(inherits(load, "try-error")){ return(NULL) }
  
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


#9. export the necessary variables and functions to the cluster----
print("* exporting gbm_objs and functions to cluster *")
clusterExport(cl, c("gbm_objs", "plot.gbm", "process_gbm"))


#10. run the function in parallel----
print("* running `process_gbm` in parallel *")
boot_pts_i2 <- parLapply(cl = cl, X = gbm_objs, fun = process_gbm)


#11. stop the cluster----
print("* stopping cluster :-)*")
stopCluster(cl)


#12. save list of 2-way interactions per bootstrap----
print("* saving RDS file *")
saveRDS(boot_pts_i2, file=file.path(root, "boot_pts_i2_conw60.rds"))

if(cc){ q() }
