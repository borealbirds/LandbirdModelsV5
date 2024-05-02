# ---
# title: National Models 5.0 - validate predictions
# author: Elly Knight
# created: April 29, 2024
# ---

#NOTES################################

# This script uses the withheld data from each bootstrap to validate the spatial predictions of the models at two levels:
# 1. Study area: uses the tiff output from '10.MosaicPredictions.R' 
# 2. BCR: uses the tiff output from '08.Predict.R'

#PREAMBLE############################

#1. Load packages----
library(tidyverse)
library(terra)
library(sf)
library(dismo)
library(DescTools)
library(carat)

#2. Set root path----
root <- "G:/Shared drives/BAM_NationalModels/NationalModels5.0"

#3. Load data file----
load(file.path(root, "data", "04_NM5.0_data_stratify.R"))

#4. Convert bird data to dataframe----
bird.df <- as.data.frame(as.matrix(bird))
bird.df$id <- as.numeric(row.names(bird.df))

#SUBUNIT POLYGONS#####################

#1. Read in country shapefiles----
can <- read_sf(file.path(root, "Regions", "CAN_adm", "CAN_adm0.shp")) |> 
  st_transform(crs=5072) 
usa <- read_sf(file.path(root, "Regions", "USA_adm", "USA_adm0.shp")) |> 
  st_transform(crs=5072) 

#2. Read in BCR shapefile----
#Remove subunit 1
bcr <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel.shp")) |> 
  dplyr::filter(subUnit!=1) |> 
  st_transform(crs=5072)

ggplot(bcr) +
  geom_sf(aes(fill=factor(subUnit)))

#3. Identify BCRs for each country----
bcr.ca <- bcr |> 
  st_intersection(can) |> 
  mutate(country="can") |> 
  dplyr::select(subUnit, country)

bcr.usa <- bcr |> 
  st_intersection(usa) |> 
  mutate(country="usa") |> 
  dplyr::select(subUnit, country)

#4. Make merged subunits----
bcr.can4142 <- bcr.ca |> 
  dplyr::filter(subUnit %in% c(41, 42)) |> 
  st_union() |> 
  nngeo::st_remove_holes() |> 
  st_sf() |> 
  mutate(country="can",
         subUnit = 4142)
st_geometry(bcr.can4142) <- "geometry"

bcr.usa41423 <- bcr.ca |> 
  dplyr::filter(subUnit %in% c(41, 42, 3)) |> 
  st_union() |> 
  nngeo::st_remove_holes() |> 
  st_sf() |> 
  mutate(country="usa",
         subUnit = 41423)
st_geometry(bcr.usa41423) <- "geometry"

bcr.usa414232 <- bcr.usa |> 
  dplyr::filter(subUnit %in% c(41, 42, 3, 2)) |> 
  st_union() |> 
  nngeo::st_remove_holes() |> 
  st_sf() |> 
  mutate(country="usa",
         subUnit = 414232)
st_geometry(bcr.usa414232) <- "geometry"

bcr.can8182 <- bcr.ca |> 
  dplyr::filter(subUnit %in% c(81, 82)) |> 
  st_union() |> 
  nngeo::st_remove_holes() |> 
  st_sf() |> 
  mutate(country="can",
         subUnit = 8182)
st_geometry(bcr.can8182) <- "geometry"

#polygons to remove
bcr.remove <- data.frame(country=c("can", "usa", "usa", "can", "usa"),
                         subUnit=c(41, 42, 3, 42, 41))

#5. Put together----
bcr.country <- rbind(bcr.ca, bcr.usa, bcr.can4142, bcr.usa41423, bcr.usa414232, bcr.can8182) |> 
  anti_join(bcr.remove) |> 
  mutate(bcr = paste0(country, subUnit))

#BCR###############

#1. Get list of predictions----
preds <- data.frame(predpath = list.files(file.path(root, "output", "predictions"), pattern="*.tiff", full.names=TRUE),
                    predfile = list.files(file.path(root, "output", "predictions"), pattern="*.tiff")) |> 
  separate(predfile, into=c("spp", "bcr", "boot", "year"), sep="_", remove=FALSE) |>  
  mutate(year = as.numeric(str_sub(year, -100, -6)),
         boot = as.numeric(boot))

#2. Get list of models----
mods <- data.frame(modpath = list.files(file.path(root, "output", "bootstraps"), pattern="*.R", full.names=TRUE),
                   modfile = list.files(file.path(root, "output", "bootstraps"), pattern="*.R")) |> 
  separate(modfile, into=c("spp", "bcr", "boot"), sep="_", remove=FALSE) |>  
  mutate(boot = as.numeric(str_sub(boot, -100, -3)))

#3. Put together----
todo <- inner_join(preds, mods)

#Check
nrow(preds)==nrow(todo)

#4. Set up loop----
for(i in 1:nrow(todo)){
  
  #1. Get loop settings----
  bcr.i <- todo$bcr[i]
  spp.i <- todo$spp[i]
  boot.i <- todo$boot[i]
  year.i <- todo$year[i]
  
  #5. Read in files----
  pred.i <- rast(preds$predpath[i])
  load(mods$modpath[i])
  
  #6. Clip prediction by BCR boundary----
  clip.i <- bcr.country |> 
    dplyr::filter(bcr==bcr.i) |> 
    vect() |> 
    terra::crop(x=pred.i, mask=TRUE)
  
  #7. Get withheld data----
  #filter to 5 year interval of prediction
  #add species presence/absence
  #add predictions
  test.i <- gridlist[bcrlist[,bcr.i],] |> 
    anti_join(visit.i) |> 
    mutate(yearpred = round(year/5)*5) |>
    dplyr::filter(yearpred==year.i) |> 
    left_join(visit) |>
    st_as_sf(coords=c("lon", "lat"), crs=4326, remove=FALSE) |> 
    st_transform(crs=5072) |> 
    st_intersection(bcr.country |> dplyr::filter(bcr==bcr.i)) |> 
    terra::extract(x=clip.i, bind=TRUE) |> 
    data.frame() |> 
    rename(prediction = lyr1) |> 
    left_join(bird.df |> 
                dplyr::select(id, todo$spp[i]) |> 
                rename(count = spp.i)) |> 
    mutate(p = ifelse(count > 0, 1, 0)) 
  
  #8. Get training data----
  train.i <- visit.i |> 
    left_join(visit) |> 
    inner_join(bcrlist |> 
                 dplyr::select(id, todo$bcr[i]) |> 
                 rename(bcr = todo$bcr[i]) |>
                 dplyr::filter(bcr==TRUE)) |> 
    inner_join(bird.df |> 
                 dplyr::select(id, todo$spp[i]) |> 
                 rename(count = todo$spp[i]))
  
  #9. Calculate total and training residual deviance----
  total.dev.i <- calc.deviance(train.i$count, rep(mean(train.i$count), nrow(train.i)), family="poisson", calc.mean = FALSE)/nrow(train.i)
  
  train.resid.i <- mean(b.i$train.error^2)
  
  #10. Determine whether can evaluate fully ----  
  #Skip to next loop if there's no withheld data or positive detections for that time period
  if(nrow(test.i)==0 | sum(test.i$p) == 0){
    
    out.i <- data.frame(spp=spp.i,
                        bcr=bcr.i,
                        boot=boot.i,
                        year=year.i,
                        trees = b.i$n.trees,
                        n.train = b.i$nTrain,
                        n.train.p = nrow(train.p.i),
                        n.train.a = nrow(visit.i) - nrow(train.p.i),
                        total.dev = total.dev.i,
                        train.resid = train.resid.i,
                        train.d2 = (totaldev.i - train.resid.i)/totaldev.i,
                        n.test.p = sum(test.i$p),
                        n.test.a = nrow(test.i) - sum(test.i$p))
    
    next
  }
  
  #11. Estimate SSB (spatial sorting bias; Hijmans 2012 Ecology)----
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
  
  #12. Dismo evaluate presence & absence----
  p.i <- test.i |> 
    dplyr::filter(p==1)

  a.i <- test.i |> 
    dplyr::filter(p==0)
  
  eval.i <- dismo::evaluate(p.i$prediction, a.i$prediction)
  
  #13. Brier score----
  brier.i <- BrierScore(resp = test.i$count, pred = test.i$prediction)
  
  #14. Calculate other count metrics----
  #Accuracy, discrimination (spearman, pearson, intercept, slope), precision (Norberg et al. 2019 Ecol Mongr, Waldock et al. 2022 Ecography)
  
  #15. Calculate test deviance & residuals----
  test.dev.i <- calc.deviance(test.i$count, test.i$prediction, family="poisson")
  
  test.resid.i <- mean(abs(test.i$count - test.i$prediction))
  
  #16. Put together----
  out.i <- data.frame(spp=spp.i,
                      bcr=bcr.i,
                      boot=boot.i,
                      year=year.i,
                      trees = b.i$n.trees,
                      n.train = b.i$nTrain,
                      n.train.p = nrow(train.p.i),
                      n.train.a = nrow(visit.i) - nrow(train.p.i),
                      total.dev = total.dev.i,
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
                      brier = brier.i)
  
  
  
  
  
}
