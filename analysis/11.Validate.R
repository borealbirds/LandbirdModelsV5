# ---
# title: National Models 5.0 - validate subunit predictions
# author: Elly Knight
# created: April 29, 2024
# ---

#NOTES################################

# This script uses the withheld data from each bootstrap to validate the spatial the models.

# Validation is done by rounding the withheld data to the nearest five years and extracting the predicted density values from the spatial prediction for each bootstrap.

# Validation is done for two spatial extents: BCR & study area.

# The script works by running only on species and bootstrap combinations that have been fully run for all bcrs (determined by whether they've been mosaicked) and all years

#TO DO: CONSIDER SPLITTING INTO TWO SCRIPTS##########

#TO DO: PARALLELIZE

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

#UNBUFFERED ATTRIBUTION########

#1. Subunit polygons----
bcr.country <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel_Unbuffered.shp"))

#2. Make visit a vector object----
visit.v <- visit |> 
  dplyr::select(id, lat, lon) |> 
  st_as_sf(coords=c("lon", "lat"), crs=4326) |> 
  st_transform(5072) |> 
  vect()

#3. Set up mosaicked for BCR attribution----
bcrdf <- data.frame(id=visit$id)
for(i in 1:nrow(bcr.country)){
  
  #4. Filter bcr shapefile----
  bcr.i <- bcr.country |> 
    dplyr::filter(row_number()==i)
  
  #5. Convert to raster for fast extraction----
  r <- rast(ext(bcr.i), resolution=1000, crs=crs(bcr.i))
  bcr.r <- rasterize(x=bcr.i, y=r, field="subUnit")
  
  #6. Extract raster value----
  bcr.out <- data.frame(subUnit=extract(x=bcr.r, y=visit.v)[,2]) |> 
    mutate(use = ifelse(!is.na(subUnit), TRUE, FALSE))
  bcrdf[,(i+1)] <- bcr.out$use
  
  print(paste0("Finished bcr ", i, " of ", nrow(bcr.country)))
  
}

colnames(bcrdf) <- c("id", paste0(bcr.country$country, bcr.country$subUnit))

#INVENTORY#########

#1. Get list of models----
booted <- data.frame(path.mod = list.files(file.path(root, "output", "bootstraps"), pattern="*.R", full.names=TRUE),
                   file.mod = list.files(file.path(root, "output", "bootstraps"), pattern="*.R")) |> 
  separate(file.mod, into=c("spp", "bcr", "boot"), sep="_", remove=FALSE) |>  
  mutate(boot = as.numeric(str_sub(boot, -100, -3)))

#2. Get list of predictions----
predicted <- data.frame(path.pred = list.files(file.path(root, "output", "predictions"), pattern="*.tiff", full.names=TRUE, recursive=TRUE),
                        file.pred = list.files(file.path(root, "output", "predictions"), pattern="*.tiff", recursive = TRUE)) |> 
  separate(file.pred, into=c("folder", "spp", "bcr", "boot", "year", "file"), remove=FALSE) |>  
  mutate(year = as.numeric(year),
         boot = as.numeric(boot)) |> 
  dplyr::select(-folder, -file)

#3. Get list of mosaics----
mosaicked <- data.frame(path.mos = list.files(file.path(root, "output", "mosaics", "predictions"), pattern="*.tiff", full.names=TRUE, recursive=TRUE),
                        file.mos = list.files(file.path(root, "output", "mosaics", "predictions"), pattern="*.tiff", recursive = TRUE)) |> 
  separate(file.mos, into=c("spp", "boot", "year", "file"), remove=FALSE) |>  
  mutate(year = as.numeric(year),
         boot = as.numeric(boot))

#4. Get previously run output----

#Test data
if(file.exists(file.path(root, "output", "validation", "Validation_TestData.RData"))){
  
  load(file.path(root, "output", "validation", "Validation_TestData.RData"))
  
  tested <- data.frame(name = names(test)) |> 
    separate(name, into=c("spp", "boot")) |> 
    mutate(boot = as.numeric(boot))
  
}  else { tested <- data.frame(spp = NA, boot = NA)}

#Validations
if(file.exists(file.path(root, "output", "validation", "Validation_Results.csv"))){
  
  evaluations <- read.csv(file.path(root, "output", "validation", "Validation_Results.csv"))
  evaluated <- evaluations |> 
    dplyr::select(spp, boot) |> 
    unique()
  
} else { evaluated <- data.frame(spp = NA, boot = NA)}

#TEST DATA COMPILATION###############

#1. Make to do list----
#only do spp*boot with all years of prediction----
nyears <- length(seq(1985, 2020, 5))

loop <- mosaicked |> 
  group_by(spp, boot) |> 
  summarize(years = n()) |> 
  ungroup() |> 
  dplyr::filter(years==nyears) |> 
  mutate(name = paste0(spp, "_", boot)) |> 
  anti_join(tested)

#2. Set up loop----
test.new <- list()
train.new <- list()
mod.new <- list()
for(i in 1:nrow(loop)){
  
  start.i <- Sys.time()
  
  #3. Get loop settings----
  spp.i <- loop$spp[i]
  boot.i <- loop$boot[i]
  
  #4. Get the bcrs----
  loop.i <- loop[i,] |> 
    left_join(booted, multiple="all")
  
  #5. Set up bcr loop----
  test.i <- list()
  train.i <- list()
  mod.i <- list()
  for(j in 1:nrow(loop.i)){
    
    #6. Get the bcr---
    bcr.j <- loop.i$bcr[j]
    
    #7. Read in files----
    load(loop.i$path.mod[j])
    
    #8. Get training data----
    train.j <- visit.i |> 
      left_join(visit, by=c("id", "year")) |> 
      inner_join(bcrlist |> 
                   dplyr::select(id, loop.i$bcr[i]) |> 
                   rename(bcr = loop.i$bcr[i]) |>
                   dplyr::filter(bcr==TRUE),
                 by = "id") |> 
      inner_join(bird.df |> 
                   dplyr::select(id, loop$spp[i]) |> 
                   rename(count = loop$spp[i]),
                 by = "id")
    
    #9. Get withheld test data----
    #add species presence/absence
    #add offset
    withheld.j <- gridlist[bcrlist[,bcr.j],] |> 
      anti_join(visit.i, by=c("id", "year", "cell")) |> 
      left_join(visit, by=c("id", "year")) |>
      left_join(bird.df |> 
                  dplyr::select(id, loop$spp[i]) |> 
                  rename(count = loop$spp[i]),
                by="id") |> 
      mutate(p = ifelse(count > 0, 1, 0)) |> 
      left_join(offsets |> 
                  dplyr::select(id, loop$spp[i]) |> 
                  rename(offset = loop$spp[i]),
                by="id") |> 
      # left_join(cov,
      #           by=c("id", "tagMethod", "method")) |> 
      rename(year1 = year) |> 
      mutate(meth.i = method,
             year = round(year1/5)*5,
             year = ifelse(year==1980, 1985, year))
    
    #10. Get list of predictions----
    #Get rounded years in data
    loop.j <- loop.i[j,] |> 
      left_join(predicted, multiple="all")
    
    #11. Set up loop to get the predictions for the withheld data----
    
    #Make data frame to hold predictions
    prediction.j <- data.frame()

    for(k in 1:nrow(loop.j)){
      
      #12. Read in prediction raster----
      rast.k <- rast(loop.j$path.pred[k])
      
      #13. filter and project data----
      withheld.k <- withheld.j |> 
        dplyr::filter(year==loop.j$year[k])
      
      vect.k <- withheld.k |> 
        st_as_sf(coords=c("lon", "lat"), crs=4326) |> 
        st_transform(crs=crs(rast.k)) |> 
        vect()
      
      #14. Extract predictions----
      withheld.k$prediction <- extract(rast.k, vect.k)$lyr1
      
      prediction.j <- rbind(prediction.j, withheld.k)
      
      cat("Year", loop.j$year[k], "\n")
      
    }
    
    #15. Use within BCR only for validation----
    #i.e., filter out the points in the buffer
    test.i[[j]] <- prediction.j |> 
      inner_join(bcrdf |> 
                   dplyr::select(id, loop.i$bcr[i]) |> 
                   rename(bcr = loop.i$bcr[i]) |>
                   dplyr::filter(bcr==TRUE),
                 by="id")
    train.i[[j]] <- train.j
    
    #16. Get the model object details----
    mod.i[[j]] <- data.frame(spp=spp.i,
                             bcr=bcr.i,
                             boot=boot.i,
                             trees = b.i$n.trees,
                             n.train = b.i$nTrain)
    
    cat("BCR", loop.i$bcr[j], ":", j, "of", nrow(loop.i), "\n")
    
  }
  
  #16. Add to the big list----
  test.new[[i]] <- test.i
  train.new[[i]] <- train.i
  mod.new[[i]] <- mod.i
  names(test.new[[i]]) <- loop.i$bcr 
  names(train.new[[i]]) <- loop.i$bcr
  names(mod.new[[i]]) <- loop.i$bcr
  
  cat("FINISHED SPECIES*BOOTSTRAP", i, "OF", nrow(loop), "\n")
  
}

#17. Save test data predictions----
names(test.new) <- loop$name
names(train.new) <- loop$name
names(mod.new) <- loop$name

if(exists("test")){test <- c(test, test.new)} else {test <- test.new}
if(exists("train")){train <- c(train, train.new)} else {train <- train.new}
if(exists("mod")){mod <- c(mod, mod.new)} else {mod <- mod.new}

save(test, train, mod, file = file.path(root, "output", "validation", "Validation_TestData.RData"))
load(file.path(root, "output", "validation", "Validation_TestData.RData"))

#VALIDATE###########

#1. Make to do list----
#only do spp*boot with all years of prediction----
loop <- data.frame(name = names(test)) |> 
  separate(name, into=c("spp", "boot")) |> 
  mutate(boot = as.numeric(boot)) |> 
  anti_join(evaluated)

#2. Set up loop-----
if(!exists("evaluations")){evaluations <- data.frame()}
for(i in 1:nrow(loop)){
  
  #3. Loop settings----
  spp.i <- loop$spp[i]
  boot.i <- loop$boot[i]
  
  #4. Get the bcrs to run individually----
  test.i <- test[[i]]
  train.i <- train[[i]]
  mod.i <- mod[[i]]
  
  #5. Add the whole study area to the list----
  test.i[[length(test.i)+1]] <- do.call(rbind, test.i)
  names(test.i)[[length(test.i)]] <- "range"
  
  train.i[[length(train.i)+1]] <- do.call(rbind, train.i)
  names(train.i)[[length(train.i)]] <- "range"
  
  mod.i[[length(mod.i)+1]] <- do.call(rbind, mod.i) |> 
    group_by(spp, boot) |> 
    summarize(trees = mean(trees),
              n.train = mean(n.train)) |> 
    ungroup() |> 
    mutate(bcr = "range")
  names(mod.i)[[length(mod.i)]] <- "range"
  
  evals <- list()
  #6. Set up to loop through the list----
  for(j in 1:length(test.i)){
    
    #7. BCR----
    bcr.j <- names(test.i)[j]
    
    #8. Get the train and test data for this iteration----
    test.j <- test.i[[j]]
    train.j <- train.i[[j]]
    
    #9. Calculate total and training residual deviance----
    totaldev.i <- calc.deviance(train.j$count, rep(mean(train.j$count), nrow(train.j)), family="poisson", calc.mean = FALSE)/nrow(train.j)
    
    train.resid.i <- mean(b.i$train.error^2)
    
    #10. Determine whether can evaluate----  
    #Skip to next loop if there's no withheld data or positive detections for that time period
    if(nrow(test.j)==0 | sum(test.j$p) == 0){
      
      #11. Save what we can----
      evals[[j]] <- data.frame(spp=spp.i,
                          bcr=bcr.j,
                          boot=boot.i,
                          trees = mod.i[[j]]$trees,
                          n.train = mod.i[[j]]$n.train,
                          n.train.p = nrow(dplyr::filter(train.j, count > 0)),
                          n.train.a = nrow(dplyr::filter(train.j, count == 0)),
                          total.dev = totaldev.i,
                          train.resid = train.resid.i,
                          train.d2 = (totaldev.i - train.resid.i)/totaldev.i,
                          n.test.p = sum(test.j$p),
                          n.test.a = nrow(test.j) - sum(test.i$p))
      
      cat("Finished BCR", j, "of", length(train.i), "\n")
      
      next
    }
    
    #12. Estimate SSB (spatial sorting bias; Hijmans 2012 Ecology)----
    p.i <- test.j |> 
      dplyr::filter(p==1) |> 
      dplyr::select(lon, lat) |> 
      as.matrix()
    
    #skip if there's only one detection
    if(nrow(p.i)==1){ssb.i <- list(NA, NA, NA)} else {
      
      a.i <- test.j |> 
        dplyr::filter(p==0) |> 
        dplyr::select(lon, lat) |> 
        as.matrix()
      
      train.p.i <- train.j |> 
        dplyr::filter(count > 0) |> 
        dplyr::select(lon, lat) |> 
        as.matrix()
      
      ssb.i <- ssb(p=p.i, a=a.i, reference=train.p.i, lonlat=TRUE)
      ssb.i[3] <- ssb.i[1]/ssb.i[2]
      
    }
    
    #13. Dismo evaluate presence & absence----
    p.i <- test.j |> 
      dplyr::filter(p==1)
    
    a.i <- test.j |> 
      dplyr::filter(p==0)
    
    eval.i <- dismo::evaluate(p.i$prediction, a.i$prediction)
    
    #14. Brier score----
    brier.i <- BrierScore(resp = test.j$count, pred = test.j$prediction)
    
    #15. Calculate other count metrics----
    #Accuracy, discrimination (spearman, pearson, intercept, slope), precision (Norberg et al. 2019 Ecol Mongr, Waldock et al. 2022 Ecography)
    
    accuracy.i <- test.j |> 
      mutate(diff = abs(count - prediction)) |> 
      summarize(accuracy = mean(diff)/mean(count))
    
    precision.i <- sd(test.j$count)/sd(test.j$prediction)
    
    lm.i <- lm(prediction ~ count, data=test.j)
    
    cor.spearman.i = cor(test.j$prediction, test.j$count, method="spearman")
    cor.pearson.i = cor(test.j$prediction, test.j$count, method="pearson")
    
    #16. Calculate test deviance & residuals----
    test.dev.i <- calc.deviance(test.j$count, test.j$prediction, family="poisson")
    
    test.resid.i <- mean(abs(test.j$count - test.j$prediction))
    
    #17. Calculate pseudo-R2----
    r2.i <- pseudo_r2(test.j$count, test.j$prediction)
    
    #18. Put together----
    evals[[j]] <- data.frame(spp=spp.i,
                           bcr=bcr.j,
                           boot=boot.i,
                           trees = mod.i[[j]]$trees,
                           n.train = mod.i[[j]]$n.train,
                           n.train.p = nrow(dplyr::filter(train.j, count > 0)),
                           n.train.a = nrow(dplyr::filter(train.j, count == 0)),
                           total.dev = totaldev.i,
                           test.dev = test.dev.i,
                           train.resid = train.resid.i,
                           test.resid = test.resid.i,
                           train.d2 = (totaldev.i - train.resid.i)/totaldev.i,
                           test.d2 = (totaldev.i - test.resid.i)/totaldev.i, 
                           n.test.p = sum(test.j$p),
                           n.test.a = nrow(test.j) - sum(test.j$p),
                           ssb.p = ssb.i[1],
                           ssb.a = ssb.i[2],
                           ssb = ssb.i[3],
                           auc = eval.i@auc,
                           cor = eval.i@cor,
                           cor.spearman = cor.spearman.i,
                           cor.pearson = cor.pearson.i,
                           brier = brier.i,
                           accuracy = accuracy.i$accuracy,
                           precision = precision.i,
                           discrim.intercept = lm.i$coefficients[1],
                           discrim.slope = lm.i$coefficients[2],
                           pseudor2 <-r2.i[1])
    
    cat("Finished BCR", j, "of", length(train.i), "\n")
    
  }
  
  #19. Bind to dataframe----
  evaluations <- rbind(evaluations, data.table::rbindlist(evals, fill=TRUE))
  
  cat("Finished evaluation", i, "of", nrow(loop), "\n")
  
}

#20. Save----
write.csv(evaluations, file.path(root, "output", "validation", "Validation_Results.csv"))