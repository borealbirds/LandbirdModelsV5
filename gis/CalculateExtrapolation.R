# ---
# Title: Calculating extrapolation of models
# Author: Anna Drake, adapted by Elly Knight
# Date: March 2024 
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
library(dsmextra)

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
tmpcl <- clusterEvalQ(cl, library(dsmextra))

#5. Set root path----
root <- "G:/Shared drives/BAM_NationalModels5"

#6. Load data package----
load(file.path(root, "data", "04_NM5.0_data_stratify.Rdata"))
rm(bird, offsets, visit, birdlist)

#7. Set crs----
#NAD83(NSRS2007)/Conus Albers projection (epsg:5072)
crs <- "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"

#8. Export objects to clusters ----
tmpcl <- clusterExport(cl, c("root", "crs", "cov", "bootlist", "bcrlist", "covlist"))

#INVENTORY#########

#1. Set desired years----
years <- seq(1985, 2020, 5)

#2. Create to do list----
todo <- data.frame(bcr = colnames(bcrlist[,-1])) |> 
  expand_grid(year=years) |> 
  arrange(-year, bcr)

#3. Determine which are already done----
done <- data.frame(path = list.files(file.path(root, "MosaicWeighting", "Extrapolation"), pattern="*.tif", full.names=TRUE, recursive=TRUE),
                   file = list.files(file.path(root, "MosaicWeighting", "Extrapolation"), pattern="*.tif", recursive = TRUE)) |> 
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
  
  #2. Determine the covs for that bcr ----
  covlist.i <- covlist |> 
    dplyr::filter(bcr==bcr.i) |> 
    pivot_longer(-bcr, names_to="cov", values_to="use") |> 
    dplyr::filter(use==TRUE) 
  
  #3. Create target dataframe of prediction raster values----
  target <- rast(file.path(root, "gis", "stacks", paste0(bcr.i, "_", year.i, ".tif"))) |> 
    as.data.frame(xy=TRUE) |> 
    dplyr::select(c("x", "y", all_of(covlist.i$cov))) |> 
    data.frame()
  
  #4. Set up bootstrap list ----
  for(j in 1:32){
    
    #5. Get visits to include----
    sample_all <- bcrlist[bcrlist[,bcr.i]>0, c("id", bcr.i)] |> 
      dplyr::filter(id %in% bootlist[[j + 2]]) |> 
      dplyr::select(id) |> 
      left_join(cov) |> 
      dplyr::select(all_of(covlist.i$cov)) |> 
      dplyr::select(where(is.numeric)) |> 
      data.frame()
    
    #6. Remove variables that are all zero----
    sample <- sample_all[,colSums(sample_all, na.rm = TRUE)!=0]
    
    #7. Compute extrapolation----
    Extrapol <- try(compute_extrapolation(samples = sample,
                                          covariate.names = names(sample),
                                          prediction.grid = target,
                                          coordinate.system = crs,
                                          verbose=F))
    
    #tidy
    rm(sample, target)
    
    if(!inherits(Extrapol, "try-error")){
      
      #8. Produce binary raster of extrapolation area----
      raster <- as(Extrapol$rasters$mic$all, "SpatRaster")
      raster[raster>0] <- 1  # extrapolated locations in each boot
      raster[raster!=1|is.na(raster)] <- 0 #eliminate NA areas
      
      #tidy
      rm(Extrapol)
      
      #9. Add to output list ----
      if(j==1){out <- raster} else {out <- c(out, raster)}
      
    }
    
  }
  
  #9. Fix names---
  names(out) <- paste0("b", seq(1:dim(out)[3]))
  
  #10. Write raster----
  writeRaster(out, file.path(root, "MosaicWeighting", "Extrapolation", paste0(bcr.i, "_", year.i, ".tif")), overwrite=T)
  
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

if(cc){ q() }