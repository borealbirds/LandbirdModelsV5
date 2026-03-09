# ---
# Title: Calculating sampling density of models 
# Author: Elly Knight
# Date: March 2026 
# ---

#NOTES################################

# Extrapolation analysis is run for each BCR-bootstrap. Because of time-varying
# prediction rasters (e.g. biomass variables), year needs to be specified as well.

#------------------------------------------------
# This code does the following:
# 1) Identify areas of a BCR where the BCR-, species- & year-specific suite of predictor variables will cause model extrapolation (Type I and Type II) 

# --------------------
# This code requires:
# 1) Downloading the "dsmextra" package from https://github.com/densitymodelling/dsmextra if running on a local or by loading the "dsmextra_fn.R" script if running on the cluster to get around the absence of rgdal on compute canada
# 2) Buffered BCR subunit shapefile produced in "gis/BufferedSubunits.R"
# -------------------

#PREAMBLE####

#1. Load packages----
print("* Loading packages on master *")
library(sf)
library(tidyverse)
library(terra)
library(parallel)

#2. Set nodes for local vs cluster----
cores <- 4

#3. Create and register clusters----
print("* Creating clusters *")
print(table(cores))
cl <- makePSOCKcluster(cores, type="PSOCK")

#4. Load packages on clusters----
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(gbm))
tmpcl <- clusterEvalQ(cl, library(tidyverse))
tmpcl <- clusterEvalQ(cl, library(Matrix))
tmpcl <- clusterEvalQ(cl, library(terra))

#5. Set root path----
root <- "G:/Shared drives/BAM_NationalModels5"

#6. Load data package----
load(file.path(root, "data", "04_NM5.0_data_stratify.Rdata"))
rm(bird, offsets, cov, birdlist, covlist)

#7. Set crs----
#NAD83(NSRS2007)/Conus Albers projection (epsg:5072)
crs <- "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"

#8. Define summary radius ----
radius <- 20000

#8. Export objects to clusters ----
tmpcl <- clusterExport(cl, c("root", "crs", "cov", "bootlist", "bcrlist"))

#INVENTORY#########

#1. Set desired years----
years <- seq(1985, 2020, 5)

#2. Create to do list----
todo <- data.frame(bcr = colnames(bcrlist[,-1])) |> 
  expand_grid(year=years) |> 
  arrange(-year, bcr)

#3. Determine which are already done----
done <- data.frame(path = list.files(file.path(root, "gis", "samplingdensity"), pattern="*.tif", full.names=TRUE, recursive=TRUE),
                   file = list.files(file.path(root, "gis", "samplingdensity"), pattern="*.tif", recursive = TRUE)) |> 
  separate(file, into=c("bcr", "year", "file"), remove=FALSE) |>  
  mutate(year = as.numeric(year)) |> 
  dplyr::select(-file)

#4. Create final to do list----
if(nrow(done) > 0){
  
  loop <- todo |> 
    anti_join(done)
  
} else { loop <- todo }

print("* Loading model loop on workers *")
tmpcl <- clusterExport(cl, c("loop"))

#WRITE FUNCTION##########

calc_extrapolation <- function(i){
  
  #1. Get model settings---
  bcr.i <- loop$bcr[i]
  year.i <- loop$year[i]
  
  #2. Get all data for the BCR ----
  visit.i <- bcrlist[bcrlist[,bcr.i]>0, c("id", bcr.i)] |> 
    left_join(visit) |> 
    dplyr::filter(round(year/2)*2==year.i) |> 
    dplyr::select(id, lon, lat, year)
  
  #4. Set up bootstrap list ----
  for(j in 1:32){
    
    #5. Get training & testing data ----
    train.j <- visit.i |> 
      dplyr::filter(id %in% bootlist[[j + 2]]) |> 
      st_as_sf(coords=c("lon", "lat"), crs=4326) |> 
      st_transform(crs)
    
    test.j <- visit.i |> 
      anti_join(train.j) |> 
      st_as_sf(coords=c("lon", "lat"), crs=4326) |> 
      st_transform(crs)
    
    #6. Get the extrap raster for a template ----
    extrap <- rast(file.path(root, "gis", "extrapolation", paste0(bcr.i, "_", year.i, ".tif")))[[1]]
    
    #7. Focal window ----
    w <- focalMat(extrap, radius, type="circle")
    
    #8. Sampling density ----
    train.r <- rasterize(vect(train.j), extrap, fun="count", background=0)
    train.f <- focal(train.r, w=w, fun="sum", na.policy="omit") |> 
      crop(extrap, mask=TRUE)
    
    #9. Testing density ----
    test.r <- rasterize(vect(test.j), extrap, fun="count", background=0)
    test.f <- focal(test.r, w=w, fun="sum", na.policy="omit") |> 
      crop(extrap, mask=TRUE)
    
    #10. Add to output list ----
    if(j==1){
      train.out <- train.f
      test.out <- test.f
    } else {
        train.out <- c(train.out, train.f)
        test.out <- c(test.out, test.f)
    }
  }
  
  #11. Fix names---
  names(train.out) <- paste0("b", seq(1:dim(train.out)[3]))
  names(test.out) <- paste0("b", seq(1:dim(test.out)[3]))
  
  #12. Write raster----
  writeRaster(train.out, file.path(root, "gis", "samplingdensity", "train", paste0(bcr.i, "_", year.i, ".tif")), overwrite=T)
  writeRaster(test.out, file.path(root, "gis", "samplingdensity", "test", paste0(bcr.i, "_", year.i, ".tif")), overwrite=T)
  
}

#11. Export to clusters----
print("* Loading function on workers *")
tmpcl <- clusterExport(cl, c("calc_extrapolation"))

#RUN EXTRAPOLATION###############

#1. Run function in parallel----
print("* Calculating extrapolation *")
mods <- parLapply(cl,
                  X=1:nrow(loop),
                  fun=calc_extrapolation)

#CONCLUDE####

#1. Close clusters----
print("* Shutting down clusters *")
stopCluster(cl)