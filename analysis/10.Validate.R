# ---
# title: National Models 5.0 - validate subunit predictions
# author: Elly Knight
# created: April 29, 2024
# ---

#NOTES################################

# This script uses the withheld data from each bootstrap to validate the spatial the models.

# The script works by running only on species and bootstrap combinations that have been fully run for all bcrs and for all years, which is determined by whether they've been mosaicked.

# Validation is done by rounding the withheld data to the nearest five years and extracting the predicted density values from the spatial prediction for each bootstrap.

# Validation is done for two spatial extents: BCR & mosaic. The mosaic extent is added to the end of the list of BCRs and compiled from the BCR data.

# This script collects a list of bootstrap model objects that will not load and deletes them at the end of the script. The bootstrap script will need to be rerun for those files.

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

#1. Get list of sets that have been mosaicked----
#We use the mosaics as input because we know they have all subunits, and the mosaics and extrapolation complete at the bootstrap level
mosaics <- data.frame(path = list.files(file.path(root, "output", "mosaics", "predictions"), pattern="*.tiff", full.names=TRUE, recursive = TRUE),
                      file = list.files(file.path(root, "output", "mosaics", "predictions"), pattern="*.tiff", recursive = TRUE)) |> 
  separate(file, into=c("spp", "boot", "year", "filetype"), remove=FALSE) |>  
  mutate(year = as.numeric(year),
         boot = as.numeric(boot))

#2. Make the todo list----
#remove spp*boot combinations that don't have all years
todo <- mosaics |> 
  dplyr::select(spp, boot) |> 
  unique() |> 
  expand_grid(year = seq(1985, 2020, 5)) |> 
  full_join(mosaics) |> 
  mutate(na = ifelse(is.na(path), 1, 0)) |>
  group_by(spp, boot) |>
  summarize(nas = sum(na)) |> 
  ungroup() |>
  dplyr::filter(nas==0)

#3. Check which have been run----
done <- data.frame(file.mean = list.files(file.path(root, "output", "validation"), pattern="*.Rdata"))  |> 
  separate(file.mean, into=c("spp", "boot", "filetype"), remove=FALSE) |> 
  mutate(boot = as.numeric(boot))

#4. Make the todo list----
loop <- anti_join(todo, done)

#TEST DATA COMPILATION###############

#1. Set up loop----
corrupt <- data.frame()
for(i in 1:nrow(loop)){
  
  start.i <- Sys.time()
  
  #2. Get loop settings----
  spp.i <- loop$spp[i]
  boot.i <- loop$boot[i]
  
  #3. Get the bcrs----
  loop.i <- birdlist[,c("bcr", spp.i)] |> 
    pivot_longer(-bcr, names_to="spp", values_to="use") |> 
    dplyr::filter(use==TRUE,
                  str_sub(bcr, 1, 3)=="can") |> 
    dplyr::select(bcr, use) |> 
    rbind(data.frame(bcr="mosaic", use=TRUE))
  
  #4. Set up bcr loop----
  test <- list()
  train <- list()
  eval <- list()
  for(j in 1:nrow(loop.i)){
    
    #5. Get the bcr---
    bcr.j <- loop.i$bcr[j]
    
    #Determine if it's the mosaic loop or not
    if(bcr.j!="mosaic"){
      
      #6. Read in files----
      mod.i <- try(load(file.path(root, "output", "bootstraps", spp.i,
                     paste0(spp.i, "_", bcr.j, "_", boot.i, ".R"))))
      
      if(inherits(mod.i, "try-error")){
        
        corrupt <- rbind(corrupt, loop[i,] |> 
                           mutate(output = "bootstraps",
                                  bcr = bcr.j))
        break
        
      }
      
      #7. Get training data----
      train.j <- visit.i |> 
        left_join(visit, by=c("id", "year")) |> 
        inner_join(bcrlist |> 
                     dplyr::select(id, loop.i$bcr[j]) |> 
                     rename(bcr = bcr.j) |>
                     dplyr::filter(bcr==TRUE),
                   by = "id") |> 
        inner_join(bird.df |> 
                     dplyr::select(id, loop$spp[i]) |> 
                     rename(count = loop$spp[i]),
                   by = "id")
      
      #8. Get withheld test data----
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
        rename(year1 = year) |> 
        mutate(meth.i = method,
               year = round(year1/5)*5,
               year = ifelse(year==1980, 1985, year))
      
    } else {
      
      #9. If it is the mosaic loop, bind all the other bcrs together----
      train.j <- do.call(rbind, train)
      withheld.j <- do.call(rbind, test)
      
    }
    
    #10. Set up loop to get the predictions for the withheld data----
    
    #years of predictions
    years <- seq(1985, 2020, 5)
    
    #Make data frame to hold predictions
    prediction.j <- data.frame()
    
    for(k in 1:length(years)){
      
      year.k <- years[k]
      
      #10. Read in prediction raster----
      if(bcr.j!="mosaic"){
        rast.k <- try(rast(file.path(root, "output", "predictions", spp.i,
                                 paste0(spp.i, "_", bcr.j, "_", boot.i, "_", year.k, ".tiff"))))
      } else {
        rast.k <- try(rast(file.path(root, "output", "mosaics", "predictions",
                                 paste0(spp.i, "_", boot.i, "_", year.k, ".tiff"))))
      }
      
      if(inherits(rast.k, "try-error")){break}

      
      #11. filter and project data----
      withheld.k <- withheld.j |> 
        dplyr::filter(year==year.k)
      
      vect.k <- withheld.k |> 
        st_as_sf(coords=c("lon", "lat"), crs=4326) |> 
        st_transform(crs=crs(rast.k)) |> 
        vect()
      
      #12. Extract predictions----
      withheld.k$prediction <- extract(rast.k, vect.k)$lyr1
      
      #remove NAs and infinite predictions
      prediction.j <- rbind(prediction.j, withheld.k) |> 
        dplyr::filter(!is.na(prediction),
                      !is.infinite(prediction))
      
      cat("Year", year.k, "\n")
      
    }
    
    if(inherits(rast.k, "try-error")){break}
    
    #13. Use within BCR only for validation----
    #i.e., filter out the points in the buffer (non-mosaic only)
    if(bcr.j!="mosaic"){
      test.j <- prediction.j |> 
        inner_join(bcrdf |> 
                     dplyr::select(id, loop.i$bcr[j]) |> 
                     rename(bcr = bcr.j) |>
                     dplyr::filter(bcr==TRUE),
                   by="id")
    } else {
      test.j <- prediction.j
    }

    #VALIDATE#############
    
    #1. Calculate total and training residual deviance----
    totaldev.i <- calc.deviance(train.j$count, rep(mean(train.j$count), nrow(train.j)), family="poisson", calc.mean = FALSE)/nrow(train.j)
    
    train.resid.i <- mean(b.i$train.error^2)
    
    #2. Determine whether can evaluate----  
    #Skip to next loop if there's no withheld data or positive detections for that time period
    if(nrow(test.j)==0 | sum(test.j$p) == 0){
      
      #3. Save what we can----
      eval[[j]] <- data.frame(spp=spp.i,
                               bcr=bcr.j,
                               boot=boot.i,
                               trees = b.i$n.trees,
                               n.train = b.i$nTrain,
                               n.train.p = nrow(dplyr::filter(train.j, count > 0)),
                               n.train.a = nrow(dplyr::filter(train.j, count == 0)),
                               total.dev = totaldev.i,
                               train.resid = train.resid.i,
                               train.d2 = (totaldev.i - train.resid.i)/totaldev.i,
                               n.test.p = sum(test.j$p),
                               n.test.a = nrow(test.j) - sum(test.j$p))
      
      train[[j]] <- train.j
      test[[j]] <- test.j
      
      cat("BCR", loop.i$bcr[j], ":", j, "of", nrow(loop.i), "\n")
      
      next
    }
    
    #4. Estimate SSB (spatial sorting bias; Hijmans 2012 Ecology)----
    #This is commented out for now because it takes a long time to run
    # p.i <- test.j |> 
    #   dplyr::filter(p==1) |> 
    #   dplyr::select(lon, lat) |> 
    #   as.matrix()
    # 
    # #skip if there's only one detection
    # if(nrow(p.i)==1){ssb.i <- list(NA, NA, NA)} else {
    #   
    #   a.i <- test.j |> 
    #     dplyr::filter(p==0) |> 
    #     dplyr::select(lon, lat) |> 
    #     as.matrix()
    #   
    #   train.p.i <- train.j |> 
    #     dplyr::filter(count > 0) |> 
    #     dplyr::select(lon, lat) |> 
    #     as.matrix()
    #   
    #   ssb.i <- ssb(p=p.i, a=a.i, reference=train.p.i, lonlat=TRUE)
    #   ssb.i[3] <- ssb.i[1]/ssb.i[2]
    #   
    # }
    
    #5. Dismo evaluate presence & absence----
    p.i <- test.j |> 
      dplyr::filter(p==1)
    
    a.i <- test.j |> 
      dplyr::filter(p==0)
    
    eval.i <- dismo::evaluate(p.i$prediction, a.i$prediction)
    
    #6. Calculate other count metrics----
    #Accuracy, discrimination (spearman, pearson, intercept, slope), precision (Norberg et al. 2019 Ecol Mongr, Waldock et al. 2022 Ecography)
    
    accuracy.i <- test.j |> 
      mutate(diff = abs(count - prediction)) |> 
      summarize(accuracy = mean(diff)/mean(count))
    
    precision.i <- sd(test.j$count)/sd(test.j$prediction)
    
    lm.i <- lm(prediction ~ count, data=test.j, na.action = "na.exclude")
    
    cor.spearman.i = cor(test.j$prediction, test.j$count, method="spearman")
    cor.pearson.i = cor(test.j$prediction, test.j$count, method="pearson")
    
    #7. Calculate test deviance & residuals----
    test.dev.i <- calc.deviance(test.j$count, test.j$prediction, family="poisson")
    
    test.resid.i <- mean(abs(test.j$count - test.j$prediction))
    
    #8. Calculate pseudo-R2----
    r2.i <- pseudo_r2(test.j$count, test.j$prediction)
    
    #9. Put together----
    eval[[j]] <- data.frame(spp=spp.i,
                             bcr=bcr.j,
                             boot=boot.i,
                             trees = b.i$n.trees,
                             n.train = b.i$nTrain,
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
                             # ssb.p = ssb.i[1],
                             # ssb.a = ssb.i[2],
                             # ssb = ssb.i[3],
                             auc = eval.i@auc,
                             cor = eval.i@cor,
                             cor.spearman = cor.spearman.i,
                             cor.pearson = cor.pearson.i,
                             accuracy = accuracy.i$accuracy,
                             precision = precision.i,
                             discrim.intercept = lm.i$coefficients[1],
                             discrim.slope = lm.i$coefficients[2],
                             pseudor2 = r2.i[1])
    
    train[[j]] <- train.j
    test[[j]] <- test.j
    
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
                          paste0(spp, "_", bcr, "_", boot, ".R")))

file.remove(unique(remove$path))
