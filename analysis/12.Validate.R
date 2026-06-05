# ---
# title: National Models 5.0 - validate models
# author: Elly Knight
# created: Feb 2, 2026
# ---

#NOTES################################

# This script uses the withheld data from each bootstrap to validate the spatial the models.

# The script works by running only on species and bootstrap combinations that have been fully run for all bcrs and for all years, which is determined by whether they've been mosaicked for all years.

# Validation is done by comparing the withheld data to predictions made from the model object. A previous version of the model script derived predictions from the spatial predictions for each year; however, the evaluation was not comparable to previous versions.

# Only withheld data from within the BCR is used; not data from the 100 km buffer

# Validation is done for two spatial extents: BCR & mosaic. The mosaic extent validation is done by summing the dataset across BCRs for that species 

# The metrics calculated are the same as used for V4 of the models. Additional metrics may be desired for future versions.

#PREAMBLE############################

#1. Load packages----
print("* Loading packages on master *")
library(tidyverse)
library(gbm)
library(epiR)
library(data.table)
library(parallel)

#2. Determine if on local or cluster----
cc <- FALSE

#3. Set nodes for local vs cluster----
if(cc){ cores <- 32}
if(!cc){ cores <- 2}

#4. Create and register clusters----
print("* Creating clusters *")
print(table(cores))
cl <- makePSOCKcluster(cores, type="PSOCK")
length(clusterCall(cl, function() Sys.info()[c("nodename", "machine")]))

#5. Set root path----
print("* Setting root file path *")
if(cc){root <- "/scratch/ecknight/NationalModels"}
if(!cc){root <- "G:/Shared drives/BAM_NationalModels5"}

#6. Load packages on clusters----
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(tidyverse))
tmpcl <- clusterEvalQ(cl, library(gbm))
tmpcl <- clusterEvalQ(cl, library(purrr))
tmpcl <- clusterEvalQ(cl, library(epiR))
tmpcl <- clusterEvalQ(cl, library(data.table))

#7. Load data file----
load(file.path(root, "data", "04_NM5.0_data_stratify.Rdata"))

#8. Get the survey ids for unique site-year data for the OCCC evaluation ----
visit$cyid <- paste0(visit$location, "-", visit$year)
id.occc <- visit[which(!duplicated(visit$cyid)), "id"]

#WRITE FUNCTIONS ###############

#1. Functions for evaluation----
#As per NMV4
pseudo_r2 <- function(observed, fitted, null=NULL, p=0) {
  if (is.null(null))
    null <- mean(observed)
  ll0 <- sum(dpois(observed, null, log=TRUE))
  lls <- sum(dpois(observed, observed, log=TRUE))
  llf <- sum(dpois(observed, fitted, log=TRUE))
  n <- length(observed)
  R2 <- 1 - (lls - llf) / (lls - ll0)
  R2adj <- 1 - (1 - R2) * ((n-1) / (n-(p+1)))
  D0 <- -2 * (ll0 - lls)
  DR <- -2 * (llf - lls)
  p_value <- 1 - pchisq(DR, length(observed)-(p+1))
  c(R2=R2, R2adj=R2adj, Deviance=D0 - DR, Dev0=D0, DevR=DR, p_value=p_value)
}

simple_roc <- function(labels, scores){
  Labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(
    TPR=cumsum(Labels)/sum(Labels),
    FPR=cumsum(!Labels)/sum(!Labels),
    Labels=Labels)
}

simple_auc <- function(ROC) {
  ROC$inv_spec <- 1-ROC$FPR
  dx <- diff(ROC$inv_spec)
  sum(dx * ROC$TPR[-1]) / sum(dx)
}

# Function to get train and test data and the predictions from the model object 
get_data_bcr <- function(k, bcr.j, spp.i, b.list){
  
  #1. Get visits ----
  visit.k <- bcrlist[bcrlist[,bcr.j]>0, c("id", bcr.j)]
  colnames(visit.k) <- c("id", "use")
  
  #2. Get year & coords ----
  year.k <- visit[visit$id %in% visit.k$id, c("year", "lon", "lat")]
  
  #3. Get response data (bird data)----
  bird.k <- bird[as.character(visit.k$id), spp.i]
  
  #4. Get the offset ----
  off.k <- offsets[offsets$id %in% visit.k$id, spp.i]
  
  #5. Put together -----
  #exclude the buffer
  #add species presence/absence
  dat.k <- cbind(visit.k, year.k, bird.k, off.k) |> 
    rename(count = bird.k,
           offset = off.k,
           year1 = year) |> 
    mutate(year = round(year1/5)*5,
           year = ifelse(year==1980, 1985, year)) |> 
    dplyr::filter(use==1) |> 
    mutate(p = ifelse(count > 0, 1, 0))
  
  #6. Get the covariates ----
  cov.k <- dat.k |> 
    left_join(cov, by="id") |> 
    dplyr::select(id, offset, count, p, all_of(b.list[[1]]$var.names))
  
  #7. Get full predictions ----
  cov.k$predraw <- suppressWarnings(predict.gbm(b.list[[k]], cov.k, type="link"))
  cov.k$predfull <- exp(cov.k$predraw + cov.k$offset)
  
  #8. Get intercept only predictions ----
  cov.k$predinit <- exp(b.list[[k]]$initF + cov.k$offset)
  
  #9. Split into train and test and occc objects ----
  train.k <- cov.k |> 
    dplyr::filter(id %in% bootlist[[k + 2]])
  test.k <- cov.k |> 
    dplyr::filter(!id %in% bootlist[[k + 2]])
  occc.k <- cov.k |> 
    dplyr::filter(id %in% id.occc$id)
  
  #10. Get the 
  
  #10. Return ----
  out.k <- list(train.k, test.k, occc.k)
  names(out.k) <- c("train", "test", "occc")
  
  return(out.k)
  
}

# Function to run the evaluation metrics for each bootstrap
evaluate_boot <- function(data){
  
  #1. Calculate prevalence ----
  prevalence <- mean(ifelse(data$train$count >0, 1, 0))
  
  #2. Get the AUCs ----
  AUC_init <- simple_auc(simple_roc(ifelse(data$test$count>0, 1, 0), data$test$predinit))
  AUC_final <- simple_auc(simple_roc(ifelse(data$test$count>0, 1, 0), data$test$predfull))
  
  #3. Pseudo-R2 -----
  pseudo_R2 <- pseudo_r2(data$test$count, data$test$predfull, data$test$predinit)[1]
  
  #4. Return ----
  eval.k <- data.frame(prevalence, AUC_init, AUC_final, pseudo_R2)
  return(eval.k)
  
}

#wrapper function
brt_evaluate <- function(i){
  
  #2. Get loop settings----
  spp.i <- loop$spp[i]
  
  #3. Get the bcrs----
  loop.i <- birdlist[,c("bcr", spp.i)] |> 
    pivot_longer(-bcr, names_to="spp", values_to="use") |> 
    dplyr::filter(use==TRUE) |> 
    dplyr::select(bcr, use)
  
  #4. Set up the BCR loop----
  dat <- list()
  eval <- list()
  preds <- list()
  oc <- list()
  for(j in 1:nrow(loop.i)){
    
    #5. Get the bcr ----
    bcr.j <- loop.i$bcr[j]
    
    #6. Get the models ----
    bload <- try(load(file.path(root, "output", "06_bootstraps", spp.i, paste0(spp.i, "_", bcr.j, ".Rdata"))))
    if(class(bload)=="try-error"){
      corrupt <- rbind(corrupt, loop.i[j,])
      break
    }
    
    #7. Get the data -----
    dat[[j]] <- lapply(.x = seq_along(b.list), function(k) {get_data_bcr(.x, bcr.j=bcr.j, spp.i=spp.i, b.list=b.list)})
    
    #8. Evaluate -----
    eval[[j]] <- do.call(rbind, lapply(dat[[j]], function(x){evaluate_boot}))
    
    #9. Calculate OCCC -----
    #Need to get a dataframe of predictions with a column for each bootstrap
    preds[[j]] <- suppressWarnings(do.call(cbind,
                                           lapply(dat[[j]], function(b){as.numeric(b$occc$predfull)})))
    oc[[j]] <- epi.occc(preds[[j]])
    
    cat(j, " ")
    
  }
  if(class(bload)=="try-error") return(NULL)
  
  #10. Add some names----
  names(dat) <- c(loop.i$bcr)
  names(eval) <- c(loop.i$bcr)
  names(preds) <- c(loop.i$bcr)
  names(oc) <- c(loop.i$bcr)
  
  #11. Get the list of mosaics to do ----
  mosaic.i <- dplyr::filter(mosaiced, spp==spp.i)
  
  #12. Set up the mosaic loop ----
  for(j in 1:nrow(mosaic.i)){
    
    #13. Get the bcrs ----
    bcr.j <- mosaic.i$region[j]
    if(bcr.j=="Canada"){bcrs.j <- colnames(bcrlist)[str_sub(colnames(bcrlist), 1, 3)=="can"]}
    if(bcr.j=="Alaska"){bcrs.j <- c("usa41423", "usa2", "usa40", "usa43", "usa5")}
    if(bcr.j=="Lower48"){bcrs.j <- c("usa5", "usa9", "usa10", "usa11", "usa12", "usa13", "usa14", "usa23", "usa28")}
    
    #14. Get the data ----
    dat.bcr <- dat[names(dat) %in% bcrs.j]
    dat[[nrow(loop.i)+j]] <- lapply(seq_len(length(dat.bcr[[1]])), function(i){
      train_i <- do.call(rbind, lapply(dat.bcr, function(r) r[[i]]$train))
      test_i <- do.call(rbind, lapply(dat.bcr, function(r) r[[i]]$test))
      list(train = train_i, test = test_i)
    })
    
    #15. Evaluate ----
    eval[[nrow(loop.i)+j]] <- purrr::map_dfr(dat[[nrow(loop.i)+j]], evaluate_boot)
    
    #16. Calculate OCCC -----
    #Need to get a dataframe of predictions with a column for each bootstrap
    preds[[nrow(loop.i)+j]] <- suppressWarnings(do.call(cbind,
                                                        lapply(dat[[nrow(loop.i)+j]], function(b){as.numeric(b$test$predfull)})))
    oc[[nrow(loop.i)+j]] <- epi.occc(preds[[nrow(loop.i)+j]])
    
    #17. 
    
    cat(j, " ")
    
  }
  
  #18. Add some more names----
  names(dat) <- c(loop.i$bcr, mosaic.i$region)
  names(eval) <- c(loop.i$bcr, mosaic.i$region)
  names(preds) <- c(loop.i$bcr, mosaic.i$region)
  names(oc) <- c(loop.i$bcr, mosaic.i$region)
  
  #19. Save the data----
  save(eval, oc, file = file.path(root, "output", "11_validation",
                                  paste0(spp.i, ".Rdata")))
  
}


#INVENTORY#########

#1. Get list of models that have been packaged ----
#this ensures completeness
packaged <- data.frame(path = list.files(file.path(root, "output", "10_packaged"), pattern="*.tif", full.names=TRUE, recursive = TRUE),
                       file = list.files(file.path(root, "output", "10_packaged"), pattern="*.tif", recursive = TRUE)) |> 
  separate(file, into=c("spp1", "bcr1", "spp", "bcr", "year", "filetype"), remove=FALSE) |>  
  mutate(year = as.numeric(year)) |> 
  dplyr::filter(spp1!="SamplingReliability")

#2. Make the todo list----
todo <- packaged |> 
  dplyr::select(spp) |> 
  unique()

todo <- data.frame(spp = list.files(file.path(root, "output", "07_predictions")))

#3. Check which have been run----
done <- data.frame(file.mean = list.files(file.path(root, "output", "12_validation"), pattern="*.Rdata"))  |> 
  separate(file.mean, into=c("spp", "filetype"), remove=FALSE)

#4. Make the todo list----
loop <- anti_join(todo, done)

#5. Get the mosaics list ----
mosaiced <- data.frame(file = list.files(file.path(root, "output", "08_mosaics"), pattern="*.tif", recursive=TRUE)) |> 
  separate(file, into=c("region", "sppfolder", "spp", "year", "filetype"), remove=FALSE) |>  
  mutate(year = as.numeric(year)) |> 
  dplyr::select(-sppfolder, -filetype, -year, -file) |> 
  unique()

# EVALUATE #########

#1. Export objects to clusters----
tmpcl <- clusterExport(cl, c("root", "loop", "mosaiced", 
                             "bird", "offsets", "bcrlist", "birdlist", "bootlist", "visit", "cov", "id.occc",
                             "pseudo_r2", "simple_roc", "simple_auc", "get_data_bcr", "evaluate_boot", "brt_evaluate"))

#2. Run BRT function in parallel----
print("* Evaluating *")
packaged <- parLapply(cl,
                      X=1:nrow(loop),
                      fun=brt_evaluate)

# CONCLUDE ####

#1. Close clusters----
print("* Shutting down clusters *")
stopCluster(cl)

if(cc){ q() }