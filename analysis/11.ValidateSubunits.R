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

#UNBUFFERED ATTRIBUTION########

#1. Subunit polygons----
bcr.country <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel_Unbuffered.shp"))

#2. Make visit a vector object----
visit.v <- visit |> 
  dplyr::select(id, lat, lon) |> 
  st_as_sf(coords=c("lon", "lat"), crs=4326) |> 
  st_transform(5072) |> 
  vect()

#3. Set up loop for BCR attribution----
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

#VALIDATE SUBUNITS###############

#1. Get list of models----
todo <- data.frame(modpath = list.files(file.path(root, "output", "bootstraps"), pattern="*.R", full.names=TRUE),
                   modfile = list.files(file.path(root, "output", "bootstraps"), pattern="*.R")) |> 
  separate(modfile, into=c("spp", "bcr", "boot"), sep="_", remove=FALSE) |>  
  mutate(boot = as.numeric(str_sub(boot, -100, -3)))

#2. Get list of predictions----
predicted <- data.frame(path.pred = list.files(file.path(root, "output", "predictions"), pattern="*.tiff", full.names=TRUE),
                        file.pred = list.files(file.path(root, "output", "predictions"), pattern="*.tiff")) |> 
  separate(file.pred, into=c("spp", "bcr", "boot", "year"), sep="_", remove=FALSE) |> 
  mutate(year = as.numeric(str_sub(year, -100, -5)),
         boot = as.numeric(boot)) |> 
  dplyr::filter(!bcr %in% c("can8182", "usa41423", "usa2"))

#3. Make todo list----
years <- seq(1985, 2020, 5)

loop <- predicted |> 
  group_by(spp, bcr, boot) |> 
  summarize(n=n()) |> 
  ungroup() |> 
  dplyr::filter(n==length(years)) |> 
  inner_join(todo)

#4. Get previously run output----
if(file.exists(file.path(root, "output", "validation", "ModelValidation_InterimOutput.RData"))){
  
  load(file.path(root, "output", "validation", "ModelValidation_InterimOutput.RData"))
  
  start <- length(out.list) + 1
  
} else {
  
    out.list <- list()
    start <- 1
    
    }

for(i in start:nrow(loop)){
  
  start.i <- Sys.time()
  
  #3. Get loop settings----
  bcr.i <- loop$bcr[i]
  spp.i <- loop$spp[i]
  boot.i <- loop$boot[i]
  
  #4. Read in files----
  load(loop$modpath[i])
  
  #5. Get training data----
  train.i <- visit.i |> 
    left_join(visit, by=c("id", "year")) |> 
    inner_join(bcrlist |> 
                 dplyr::select(id, loop$bcr[i]) |> 
                 rename(bcr = loop$bcr[i]) |>
                 dplyr::filter(bcr==TRUE),
               by = "id") |> 
    inner_join(bird.df |> 
                 dplyr::select(id, loop$spp[i]) |> 
                 rename(count = loop$spp[i]),
               by = "id")
  
  #6. Get withheld test data----
  #add species presence/absence
  #add offset
  withheld.i <- gridlist[bcrlist[,bcr.i],] |> 
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
    left_join(cov,
              by=c("id", "tagMethod", "method")) |> 
    rename(year1 = year) |> 
    mutate(meth.i = method,
           year = round(year1/5)*5,
           year = ifelse(year==1980, 1985, year)) 
  
  #7. Get predictions----
  
  #Get rounded years in data
  years.i <- sort(unique(withheld.i$year))
  
  #Make data frame to hold predictions
  prediction.i <- data.frame()
  
  #Loop through years
  for(j in 1:length(years.i)){
    
    #Read in prediction raster
    rast.j <- try(rast(file.path(root, "output", "predictions",
                             paste0(spp.i, "_", bcr.i, "_", boot.i, "_", years.i[j], ".tiff"))))
    
    #End loop if raster doesn't load
    if(inherits(rast.j, "try-error")){ break }
    
    #filter and project data
    withheld.j <- withheld.i |> 
      dplyr::filter(year==years.i[j])
    
    vect.j <- withheld.j |> 
      st_as_sf(coords=c("lon", "lat"), crs=4326) |> 
      st_transform(crs=crs(rast.j)) |> 
      vect()
    
    #Extract predictions
    withheld.j$prediction <- extract(rast.j, vect.j)$lyr1
    
    prediction.i <- rbind(prediction.i, withheld.j)
    
    cat("Year", years.i[j], "\n")
    
  }
  
  #Go to next in loop
  if(inherits(rast.j, "try-error")){ next }
  
  #8. Use within BCR only for validation----
  test.i <- prediction.i |> 
     inner_join(bcrdf |> 
                 dplyr::select(id, loop$bcr[i]) |> 
                 rename(bcr = loop$bcr[i]) |>
                 dplyr::filter(bcr==TRUE),
                by="id")

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
                        n.test.a = nrow(test.i) - sum(test.i$p),
                        duration = as.numeric(difftime(Sys.time(), start.i, units="mins")))
    
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
    dplyr::select(id, year, cell, count, offset, prediction) |> 
    mutate(spp=spp.i,
           bcr=bcr.i,
           boot=boot.i)
  
  write.csv(test, file.path(root, "output", "validation", "data",
                            paste0(spp.i, "_", bcr.i, "_", boot.i, ".csv")), row.names = FALSE)
  
  print(paste0("Finished evaluation ", i, " of ", nrow(loop), " in ", out.vals$duration, " minutes"))
  

}

#17. Package and save----
out <- data.table::rbindlist(out.list, fill=TRUE)

write.csv(out, file.path(root, "output", "validation", "ModelValidation_BCR.csv"), row.names = FALSE)
