# ---
# title: National Models 5.0 - validate subunit predictions
# author: Elly Knight
# created: April 29, 2024
# ---

#NOTES################################

# This script uses the withheld data from each bootstrap to validate the spatial the models.

# The script works by running only on species and bootstrap combinations that have been fully run for all bcrs and for all years, which is determined by whether they've been mosaicked for all years.

# Validation is done by rounding the withheld data to the nearest five years and extracting the predicted density values from the spatial prediction for each bootstrap.

# Validation is done for two spatial extents: BCR & mosaic. The mosaic extent is added to the end of the list of BCRs and compiled from the BCR data.

# This script collects a list of bootstrap model objects that will not load and deletes them at the end of the script. The bootstrap script will need to be rerun for those files.

#TO DO FOR V6: Convert to functions and tidy

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
load(file.path(root, "data", "04_NM5.0_data_stratify.Rdata"))

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

#INVENTORY#########

#1. Get list of sets that have been mosaiced----
#We inventory from the mosaics because we also evaluate the mosaic performance
mosaics <- data.frame(path = list.files(file.path(root, "output", "08_mosaics"), pattern="*.tif", full.names=TRUE, recursive = TRUE),
                      file = list.files(file.path(root, "output", "08_mosaics"), pattern="*.tif", recursive = TRUE)) |> 
  separate(file, into=c("spp", "year", "filetype"), remove=FALSE) |>  
  mutate(year = as.numeric(year))

#2. Make the todo list----
#remove spp*boot combinations that don't have all years, we need them to get the predicted values
todo <- mosaics |> 
  dplyr::select(spp) |> 
  unique() |> 
  expand_grid(year = seq(1985, 2020, 5)) |> 
  full_join(mosaics) |> 
  mutate(na = ifelse(is.na(path), 1, 0)) |>
  group_by(spp) |>
  summarize(nas = sum(na)) |> 
  ungroup() |>
  dplyr::filter(nas==0)

#3. Check which have been run----
done <- data.frame(file.mean = list.files(file.path(root, "output", "11_validation"), pattern="*.Rdata"))  |> 
  separate(file.mean, into=c("spp", "filetype"), remove=FALSE)

#4. Make the todo list----
loop <- anti_join(todo, done)

#WRITE FUNCTIONS###########

#TEST DATA COMPILATION###############

#1. Set up loop----
corrupt <- data.frame()
for(i in 1:nrow(loop)){
  
  start.i <- Sys.time()
  
  #2. Get loop settings----
  spp.i <- loop$spp[i]
  
  #3. Get the bcrs----
  loop.i <- birdlist[,c("bcr", spp.i)] |> 
    pivot_longer(-bcr, names_to="spp", values_to="use") |> 
    dplyr::filter(use==TRUE) |> 
    dplyr::select(bcr, use) |> 
    rbind(data.frame(bcr="mosaic", use=TRUE))
  
  #4. Set up bcr loop----
  train <- list()
  test <- list()
  eval <- list()
  for(j in 1:nrow(loop.i)){
    
    #5. Get the bcr---
    bcr.j <- loop.i$bcr[j]
    
    #GET TRAIN & TEST DATA####
    
    #1. Determine if it's the mosaic loop or not
    if(bcr.j!="mosaic"){
      
      #2. Set up bootstrap loop----
      train.j <- list()
      test.j <- list()
      for(k in 1:32){
        
        #3. Get visits ----
        visit.k <- bcrlist[bcrlist[,bcr.j]>0, c("id", bcr.j)]
        colnames(visit.k) <- c("id", "use")
        
        #4. Get year & coords ----
        year.k <- visit[visit$id %in% visit.k$id, c("year", "lon", "lat")]
        
        #5. Get response data (bird data)----
        bird.k <- bird[as.character(visit.k$id), spp.i]
        
        #6. Get the offset ----
        off.k <- offsets[offsets$id %in% visit.k$id, spp.i]
        
        #7. Put together -----
        dat.k <- cbind(visit.k, year.k, bird.k, off.k) |> 
          rename(count = bird.k,
                 offset = off.k,
                 year1 = year) |> 
          mutate(year = round(year1/5)*5,
                 year = ifelse(year==1980, 1985, year))
        
        #8. Get the training data for this bootstrap ----
        train.j[[k]] <- dat.k |> 
          dplyr::filter(id %in% bootlist[[k + 2]])
        
        #9. Get the with withheld test data ----
        #exclude the buffer
        #add species presence/absence
        test.j[[k]] <- dat.k |> 
          dplyr::filter(!id %in% bootlist[[k + 2]],
                        use==1) |> 
          mutate(p = ifelse(count > 0, 1, 0))
        
      }
      
    } else {
      
      #10. If it is the mosaic loop, bind all the other bcrs together----
      train.j <- do.call(rbind, train)
      test.j <- do.call(rbind, test)
      
    }
    
    #GET PREDICTIONS####
    
    #1. Set up loop to get the predictions for the withheld data----
    
    #years of predictions
    years <- seq(1985, 2020, 5)
    
    #Make data frame to hold predictions
    prediction.j <- list()
    
    for(k in 1:length(years)){
      
      year.k <- years[k]
      
      #2. Read in prediction raster----
      if(bcr.j!="mosaic"){
        rast.k <- try(rast(file.path(root, "output", "07_predictions", spp.i,
                                 paste0(spp.i, "_", bcr.j, "_", year.k, ".tif"))))
      } else {
        rast.k <- try(rast(file.path(root, "output", "08_mosaics", "predictions",
                                 paste0(spp.i, "_", year.k, ".tif"))))
      }
      
      if(inherits(rast.k, "try-error")){break}

      #3. filter and project data----
      test.k <- lapply(test.j, function(df) df[df$year==year.k,])
      
      vect.k <- lapply(test.k, function(df){
        vect(df, geom=c("lon", "lat"), crs="EPSG:4326")
      })
      
      #4. Extract predictions----
      for(l in 1:length(test.k)){
        test.k[[l]]$prediction <- extract(rast.k[[l]], vect.k[[l]])[[2]]
      }
      
      #5. remove NAs and infinite predictions----
      prediction.k <- lapply(test.k, function(df){
        df[is.finite(df$prediction) & !is.na(df$prediction),]
      })
      
      #6. Put the years together ----
      for(l in 1:length(prediction.k)){
        if(k==1){
          prediction.j[[l]] <- prediction.k[[l]]
        } else {
          prediction.j[[l]] <- rbind(prediction.j[[l]], prediction.k[[l]])
        }
        
      }
      
      cat("Year", year.k, "\n")
      
    }

    #VALIDATE#############
    
    #1. Calculate total deviance----
    totaldev.i <- sapply(train.j, function(df){
      calc.deviance(df$count, rep(mean(df$count), nrow(df)), family="poisson", calc.mean = FALSE)/nrow(df)
    })
    
    #2. Get total training residual deviance----
    load(file.path(root, "output", "06_bootstraps", spp.i, paste0(spp.i, "_", bcr.j, ".Rdata")))
    
    train.resid.i <- c()
    for(k in 1:length(train.j)){
      train.resid.i <- c(train.resid.i, mean(b.list[[k]]$train.error^2))
    }
    
    #2. Determine whether can evaluate----  
    #Skip to next loop if there's no withheld data or positive detections for that time period
    
    if(any(sapply(test.j, function(df) sum(df$p)==0)) | 
       any(sapply(test.j, function(df) nrow(df)==0))){
      
      #3. Save what we can----
      eval[[j]] <- data.frame(spp=spp.i,
                               bcr=bcr.j,
                               boot=c(1:32),
                               trees = b.list[[1]]$n.trees,
                               n.train = sapply(b.list, function(m) m$nTrain),
                               n.train.p = sapply(train.j, function(df) nrow(dplyr::filter(df, count > 0))),
                               n.train.a = sapply(train.j, function(df) nrow(dplyr::filter(df, count == 0))),
                               total.dev = totaldev.i,
                               train.resid = train.resid.i,
                               train.d2 = (totaldev.i - train.resid.i)/totaldev.i,
                               n.test.p = sapply(test.j, function(df) sum(df$p)),
                               n.test.a = sapply(test.j, function(df) nrow(df) - sum(df$p)))
      
      train[[j]] <- train.j
      test[[j]] <- prediction.j
      
      cat("BCR", loop.i$bcr[j], ":", j, "of", nrow(loop.i), "\n")
      
      next
    }
    
    #5. Dismo evaluate presence & absence----
    p.i <- lapply(prediction.j, function(df) dplyr::filter(df, p==1))
    a.i <- lapply(prediction.j, function(df) dplyr::filter(df, p==0))
    
    eval.i <- mapply(function(a, p) {
      dismo::evaluate(p = p$prediction, a=a$prediction)
    }, a.i, p.i)
    
    #6. Calculate other count metrics----
    #Accuracy, discrimination (spearman, pearson, intercept, slope), precision (Norberg et al. 2019 Ecol Mongr, Waldock et al. 2022 Ecography)
    
    accuracy.i <- sapply(prediction.j, function(df){
      df |> 
        mutate(diff = abs(count - prediction)) |> 
        summarize(accuracy = mean(diff)/mean(count))
    }) |> 
      unlist() |> 
      unname()
    
    precision.i <- sapply(prediction.j, function(df){
      sd(df$count)/sd(df$prediction)
    })
    
    lm.i <- lapply(prediction.j, function(df){
      lm(prediction ~ count, data = df, na.action="na.exclude")
    })
    
    cor.spearman.i <- sapply(prediction.j, function(df){
      cor(df$prediction, df$count, method="spearman")
    })
    
    cor.pearson.i <- sapply(prediction.j, function(df){
      cor(df$prediction, df$count, method="pearson")
    })
    
    #7. Calculate test deviance & residuals----
    test.dev.i <- sapply(prediction.j, function(df){
      calc.deviance(df$count, df$prediction, family="poisson")
    })
    
    test.resid.i <- sapply(prediction.j, function(df){
      mean(abs(df$count - df$prediction))
    })
    
    #8. Calculate pseudo-R2----
    r2.i <- sapply(prediction.j, function(df){
      pseudo_r2(df$count, df$prediction)
    }) |> 
      t() |> 
      data.frame()
    
    #9. Put together----
    eval[[j]] <- data.frame(spp=spp.i,
                            bcr=bcr.j,
                            boot=c(1:32),
                            trees = b.list[[1]]$n.trees,
                            n.train = sapply(b.list, function(m) m$nTrain),
                            n.train.p = sapply(train.j, function(df) nrow(dplyr::filter(df, count > 0))),
                            n.train.a = sapply(train.j, function(df) nrow(dplyr::filter(df, count == 0))),
                            total.dev = totaldev.i,
                            train.resid = train.resid.i,
                            train.d2 = (totaldev.i - train.resid.i)/totaldev.i,
                            n.test.p = sapply(test.j, function(df) sum(df$p)),
                            n.test.a = sapply(test.j, function(df) nrow(df) - sum(df$p)),
                            auc = sapply(eval.i, function(x) x@auc),
                            cor = sapply(eval.i, function(x) x@cor),
                            cor.spearman = cor.spearman.i,
                            cor.pearson = cor.pearson.i,
                            accuracy = accuracy.i,
                            precision = precision.i,
                            discrim.intercept = sapply(lm.i, function(x) x$coefficients[1]),
                            discrim.slope = sapply(lm.i, function(x) x$coefficients[2]),
                            pseudor2 = r2.i$R2adj)
    
    train[[j]] <- train.j
    test[[j]] <- prediction.j
    
    cat("BCR", loop.i$bcr[j], ":", j, "of", nrow(loop.i), "\n")

  }
  
  if(inherits(rast.k, "try-error")){next}
  if(inherits(mod.i, "try-error")){next}
  
  #15. Add some names----
  names(test) <- loop.i$bcr
  names(train) <- loop.i$bcr
  names(eval) <- loop.i$bcr
  
  #16. Save the data----
  save(test, train, eval, file = file.path(root, "output", "validation",
                                                          paste0(spp.i, "_", boot.i, ".Rdata")))
  
  end.i <- Sys.time()
  
  cat("FINISHED MODEL VALIDATION FOR", i, "OF", nrow(loop), "in", difftime(end.i, start.i, units="mins"), "minutes\n")
  
}

#17. Delete corrupt bootstraps----
remove <- corrupt |> 
  mutate(path = file.path(root, "output", "bootstraps", spp,
                          paste0(spp, "_", bcr, "_", boot, ".Rdata")))

file.remove(unique(remove$path))
