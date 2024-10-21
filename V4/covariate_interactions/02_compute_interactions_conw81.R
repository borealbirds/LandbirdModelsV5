# ---
# title: National Models 4.0 - estimating and ranking 2-way covariate interactions for CONW and CAWA
# author: Mannfred Boehm
# created: September 13, 2024
# ---

#NOTES################################

# This script creates a list of dataframes containing relative influence per covariate
# Every list element represents a bootstrap sample 




#1. attach packages----
print("* attaching packages on master *")
library(gbm)
library(Matrix)
library(parallel)
library(tidyverse)


#2. define local or cluster
test <- FALSE
cc <- TRUE

#3. set number of tasks for local vs cluster----
if(cc){ n_tasks <- 32}
if(!cc | test){ n_tasks <- 4}


#4. create and register clusters----
# creates 16 copies of R running in parallel via 16 tasks, on one of the cluster's sockets (processors). 
# Belgua has ~965 nodes
print("* creating clusters *")
cl <- makePSOCKcluster(n_tasks, type="PSOCK")

# print number of tasks and host name for confirmation
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
tmpcl <- clusterEvalQ(cl, library(Matrix))
tmpcl <- clusterEvalQ(cl, library(tidyverse))



#7. index gbm objects for CONW and CAWA----
print("* indexing gbm objects for CONW BCR81 *")

# gbm objects for CAWA and CONW have been transfered to the symbolic link file "scratch"
# dirs <- c("CAWA/BCR_60", "CAWA/BCR_12", "CONW/BCR_60", "CONW/BCR_81")
# gbm_objs <- unlist(lapply(file.path(root, dirs), list.files, full.names = TRUE))
dirs <- "CONW/BCR_81"
gbm_objs <- list.files(file.path(root, dirs), full.names = TRUE)
message("Found ", length(gbm_objs), " files: ", gbm_objs)


# load in V4 count data
try_load <- suppressWarnings(try(load(file.path(root, "BAMdb-GNMsubset-2020-01-08.RData")), silent=TRUE))

# check if `load` was successful and V4 data exists
if (inherits(try_load, "try-error") || !exists("dd")){
  message("could not V4 data from: ", root)
} else {
  message("successfully loaded V4 data from: ", root)
} #close else()



#8. create function to process each gbm object (parallelized)----
process_gbm <- function(obj_path) {
  
  
  # load a bootstrap replicate (gbm object)
  try_load <- suppressWarnings(try(load(file.path(obj_path)), silent = TRUE))
  
  # check if `load` was successful and `out` exists
  if (inherits(try_load, "try-error") || !exists("out")) {
    message("could not load gbm_obj from: ", obj_path)
  } else {
    message("successfully loaded gbm_obj from: ", obj_path)
  } #close else()
  flush.console()
  
  
  # build `DAT` for the target species x bcr
  # define species and region
  spp <- attr(out, "__settings__")$species
  bcr <- attr(out, "__settings__")$region
  
  # create an index indicating where in `dd` we will find the current BCR's data
  # ss = subset
  ss <- dd[, bcr] == 1L
  
  # reconstruct DAT object (see: https://github.com/borealbirds/GNM/blob/master/R/04-boot.R)
  # `yy` contains the counts data, and is subset using `ss` and `spp` 
  DAT <- 
    data.frame(
      count = as.numeric(yy[ss, spp]),  
      offset = off[ss, spp],
      cyid = dd$cyid[ss],
      YEAR = dd$YEAR[ss],
      ARU = dd$ARU[ss],   
      dd2[ss, out$var.names[-c(1:2)]]  #-c(1:2) because we already have `YEAR` and `ARU` from `dd`
    )
  
  # keep just one observation (row) for each unique "cell x year" combination
  DAT <- DAT[sample.int(nrow(DAT)),]
  DAT <- DAT[!duplicated(DAT$cyid),]
  
  # ensure `DAT` has the covariates found in `out`
  required_vars <- out$var.names
  if (!all(required_vars %in% colnames(DAT))) {
    stop("one or more covariates in `out` are missing from `DAT`")
  } else {
    message("check: all covariates in `out` are in `DAT`")
  } 
  flush.console()

  # find evaluation points (`data.frame`) for every covariate permutation of degree 2 (indexed by i,j) 
  pts <- list()
  n <- length(out$var.names); message("check: `out` has ", n, " covariates")
  interaction_index <- 1 # starts at 1 and increases for every covariate interaction computed. Resets at one when moving to the next bootstrap model.
  flush.console()

  # end at n-1 to avoid finding the interaction of variable n x variable n
  for (i in 1:(n-1)) {
  message("computing all interactions with ", out$var.names[i])
  flush.console()

    # start at i+1 to avoid finding the interaction of variable 1 x variable 1
    for (j in (i+1):n) {
      
      # set seed for reproducibility of the data subsampling procedure
      set.seed(interaction_index)
      
      # check if the response data is non-empty
      if (sum(DAT$count) < 1) {
        pts[[interaction_index]] <- NA  # assign NA if there are no occurrence records
      
        } else {
        
        # subsample the data using `slice_sample` to speed up `interact.gbm`
        DAT_sample <- 
          dplyr::slice_sample(DAT, prop = 0.25) |> 
          dplyr::select(all_of(required_vars)) # select only the necessary covariates (i.e those in `out`)
        
        # check for variability in the selected covariates
        var_i_unique <- length(unique(DAT_sample[[out$var.names[i]]]))
        var_j_unique <- length(unique(DAT_sample[[out$var.names[j]]]))
        
        if (var_i_unique < 2 || var_j_unique < 2) {
          message("skipping interaction for ", out$var.names[i], " and ", out$var.names[j], " due to insufficient variability.")
          pts[[interaction_index]] <- NA
          flush.console()
        } else {
          
        
        pts[[interaction_index]] <- gbm::interact.gbm(x = out, data = DAT_sample, i.var = c(i, j))  # test the interaction between variable i and j
        
        } # close nested else()
        
      } # close top else ()
      
      # label the interaction
      names(pts)[interaction_index] <- paste(out$var.names[i], out$var.names[j], sep = ".")  # Label the interaction
      interaction_index <- interaction_index + 1
    } # close nested loop
    
  } # close top loop
  
  return(pts)  # return the list of interaction points
  
} # close function


#9. export the necessary variables and functions to the cluster----
print("* exporting gbm_objs and functions to cluster *")
clusterExport(cl, c("gbm_objs", "interact.gbm", "process_gbm", "dd", "dd2", "yy", "off"))


#10. run the function in parallel----
print("* running `process_gbm` in parallel *")
boot_pts_i2 <- parLapply(cl = cl, X = gbm_objs, fun = process_gbm)


#11. stop the cluster----
print("* stopping cluster :-)*")
stopCluster(cl)


#12. save list of 2-way interactions per bootstrap----
print("* saving RDS file *")
saveRDS(boot_pts_i2, file=file.path(root, "boot_pts_i2_conw81.rds"))

if(cc){ q() }

