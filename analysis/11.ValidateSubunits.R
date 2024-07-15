# ---
# title: National Models 5.0 - validate subunit predictions
# author: Elly Knight
# created: April 29, 2024
# ---

#NOTES################################

# This script uses the withheld data from each bootstrap to validate the spatial the models.

# Validation is done by using the raw model object to predict values to the withheld data instead of spatial predictions because predictions were made for particular years, but the withheld data is across all years.

# Validation is done for two spatial extents:
# 1. BCR: uses the single bootstrapped model to predict. Points that were used to build the model in the 100 km buffer are removed from validation because the packaged subunit predictions exclude the buffer.
# 2. Study area: The predicted values from each subunit model are saved out in this script and used in the subsequent script `12.ValidateMosaics.R` to evaluate the mosaic predictions.

#PREAMBLE############################

#1. Load packages----
library(tidyverse)
library(terra)
library(sf)
library(dismo)
library(gbm)
library(DescTools)
library(carat)

#2. Set root path----
root <- "G:/Shared drives/BAM_NationalModels5"

#3. Load data file----
load(file.path(root, "data", "04_NM5.0_data_stratify.R"))

#4. Convert bird data to dataframe----
bird.df <- as.data.frame(as.matrix(bird))
bird.df$id <- as.numeric(row.names(bird.df))

#5. Function for pseudoR2----
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

#6. Subunit polygons----
bcr.country <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel_Unbuffered.shp"))

#BCR MODELS###############

#1. Get list of models----
todo <- data.frame(modpath = list.files(file.path(root, "output", "bootstraps"), pattern="*.R", full.names=TRUE),
                   modfile = list.files(file.path(root, "output", "bootstraps"), pattern="*.R")) |> 
  separate(modfile, into=c("spp", "bcr", "boot"), sep="_", remove=FALSE) |>  
  mutate(boot = as.numeric(str_sub(boot, -100, -3)))

#2. Get previously run output----
if(file.exists(file.path(root, "output", "validation", "ModelValidation_InterimOutput.RData"))){
  
  load(file.path(root, "output", "validation", "ModelValidation_InterimOutput.RData"))
  
  start <- length(out.list) + 1
  
} else {
  
    out.list <- list()
    start <- 1
    
    }

for(i in start:nrow(todo)){
  
  start.i <- Sys.time()
  
  #3. Get loop settings----
  bcr.i <- todo$bcr[i]
  spp.i <- todo$spp[i]
  boot.i <- todo$boot[i]
  
  #4. Read in files----
  load(todo$modpath[i])
  
  #5. Get training data----
  train.i <- visit.i |> 
    left_join(visit) |> 
    inner_join(bcrlist |> 
                 dplyr::select(id, todo$bcr[i]) |> 
                 rename(bcr = todo$bcr[i]) |>
                 dplyr::filter(bcr==TRUE)) |> 
    inner_join(bird.df |> 
                 dplyr::select(id, todo$spp[i]) |> 
                 rename(count = todo$spp[i]))
  
  #6. Get withheld test data----
  #add species presence/absence
  #add offset
  withheld.i <- gridlist[bcrlist[,bcr.i],] |> 
    anti_join(visit.i) |> 
    left_join(visit) |>
    left_join(bird.df |> 
                dplyr::select(id, todo$spp[i]) |> 
                rename(count = todo$spp[i])) |> 
    mutate(p = ifelse(count > 0, 1, 0)) |> 
    left_join(offsets |> 
                dplyr::select(id, todo$spp[i]) |> 
                rename(offset = todo$spp[i])) |> 
    left_join(cov) |> 
    mutate(meth.i = method)
  
  #7. Make predictions----
  withheld.i$fitted <- predict(b.i, withheld.i)
  withheld.i$prediction <- exp(withheld.i$fitted + withheld.i$offset)
  
  #8. Clip to BCR boundary for validation----
  test.i <- withheld.i  |> 
    st_as_sf(coords=c("lon", "lat"), crs=4326, remove=FALSE) |> 
    st_transform(crs=5072) |> 
    st_intersection(bcr.country |> dplyr::filter(bcr==bcr.i)) |> 
    st_drop_geometry()

  #9. Calculate total and training residual deviance----
  totaldev.i <- calc.deviance(train.i$count, rep(mean(train.i$count), nrow(train.i)), family="poisson", calc.mean = FALSE)/nrow(train.i)
  
  train.resid.i <- mean(b.i$train.error^2)
  
  #10. Determine whether can evaluate----  
  #Skip to next loop if there's no withheld data or positive detections for that time period
  if(nrow(test.i)==0 | sum(test.i$p) == 0){
    
    out.i <- data.frame(spp=spp.i,
                        bcr=bcr.i,
                        boot=boot.i,
                        trees = b.i$n.trees,
                        n.train = b.i$nTrain,
                        n.train.p = nrow(dplyr::filter(train.i, count > 0)),
                        n.train.a = nrow(dplyr::filter(train.i, count == 0)),
                        total.dev = totaldev.i,
                        train.resid = train.resid.i,
                        train.d2 = (totaldev.i - train.resid.i)/totaldev.i,
                        n.test.p = sum(test.i$p),
                        n.test.a = nrow(test.i) - sum(test.i$p))
    
    next
  }
  
  #10. Estimate SSB (spatial sorting bias; Hijmans 2012 Ecology)----
  p.i <- test.i |> 
    dplyr::filter(p==1) |> 
    dplyr::select(lon, lat) |> 
    as.matrix()
  
  #cheat and add another row if there's only one detection
  if(nrow(p.i)==1){p.i <- rbind(p.i, p.i)}
  
  a.i <- test.i |> 
    dplyr::filter(p==0) |> 
    dplyr::select(lon, lat) |> 
    as.matrix()

  train.p.i <- train.i |> 
    dplyr::filter(count > 0) |> 
    dplyr::select(lon, lat) |> 
    as.matrix()
  
  ssb.i <- ssb(p=p.i, a=a.i, reference=train.p.i, lonlat=TRUE)
  
  #11. Dismo evaluate presence & absence----
  p.i <- test.i |> 
    dplyr::filter(p==1)

  a.i <- test.i |> 
    dplyr::filter(p==0)
  
  eval.i <- dismo::evaluate(p.i$prediction, a.i$prediction)
  
  #12. Brier score----
  brier.i <- BrierScore(resp = test.i$count, pred = test.i$prediction)
  
  #13. Calculate other count metrics----
  #Accuracy, discrimination (spearman, pearson, intercept, slope), precision (Norberg et al. 2019 Ecol Mongr, Waldock et al. 2022 Ecography)
  
  accuracy.i <- test.i |> 
    mutate(diff = abs(count - prediction)) |> 
    summarize(accuracy = mean(diff)/mean(count))
  
  precision.i <- sd(test.i$count)/sd(test.i$prediction)
  
  lm.i <- lm(prediction ~ count, data=test.i)
  
  cor.spearman.i = cor(test.i$prediction, test.i$count, method="spearman")
  cor.pearson.i = cor(test.i$prediction, test.i$count, method="pearson")
  
  #14. Calculate test deviance & residuals----
  test.dev.i <- calc.deviance(test.i$count, test.i$prediction, family="poisson")
  
  test.resid.i <- mean(abs(test.i$count - test.i$prediction))
  
  #15. Calculate pseudo-R2----
  r2.i <- pseudo_r2(test.i$count, test.i$prediction)
  
  #16. Put together----
  out.vals <- data.frame(spp=spp.i,
                      bcr=bcr.i,
                      boot=boot.i,
                      trees = b.i$n.trees,
                      n.train = b.i$nTrain,
                      n.train.p = nrow(dplyr::filter(train.i, count > 0)),
                      n.train.a = nrow(dplyr::filter(train.i, count == 0)),
                      total.dev = totaldev.i,
                      test.dev = test.dev.i,
                      train.resid = train.resid.i,
                      test.resid = test.resid.i,
                      train.d2 = (totaldev.i - train.resid.i)/totaldev.i,
                      test.d2 = (totaldev.i - test.resid.i)/totaldev.i, 
                      n.test.p = sum(test.i$p),
                      n.test.a = nrow(test.i) - sum(test.i$p),
                      ssb.p = ssb.i[1],
                      ssb.a = ssb.i[2],
                      ssb = ssb.i[1]/ssb.i[2],
                      auc = eval.i@auc,
                      cor = eval.i@cor,
                      cor.spearman = cor.spearman.i,
                      cor.pearson = cor.pearson.i,
                      brier = brier.i,
                      accuracy = accuracy.i$accuracy,
                      precision = precision.i,
                      discrim.intercept = lm.i$coefficients[1],
                      discrim.slope = lm.i$coefficients[2],
                      pseudor2 <-r2.i[1],
                      duration = as.numeric(difftime(Sys.time(), start.i, units="mins")))
  
  #17. Save interim object----
  out.list[[i]] <- out.vals
  
  save(out.list, file = file.path(root, "output", "validation", "ModelValidation_InterimOutput.RData"))
  
  #18. Save test data for national evaluation----
  test <- test.i |> 
    dplyr::select(id, year, cell, count, offset, fitted, prediction) |> 
    mutate(spp=spp.i,
           bcr=bcr.i,
           boot=boot.i)
  
  write.csv(test, file.path(root, "output", "validation", "data",
                            paste0(spp.i, "_", bcr.i, "_", boot.i, ".csv")), row.names = FALSE)
  
  print(paste0("Finished evaluation ", i, " of ", nrow(todo)))
  

}

#17. Package and save----
out <- data.table::rbindlist(out.list, fill=TRUE)

write.csv(out, file.path(root, "output", "validation", "ModelValidation_BCR.csv"), row.names = FALSE)
