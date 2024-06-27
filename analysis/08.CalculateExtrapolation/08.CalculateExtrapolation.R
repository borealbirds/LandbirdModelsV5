#-----------------------------------------------------
# Title: Calculating extrapolation of models
# Author: Anna Drake, adapted by Elly Knight
# Date: March 2024           
#-----------------------------------------------------------------------------

# Extrapolation analysis is run for each species-BCR-bootstrap. Because of time-varying
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

#2. Determine if testing and on local or cluster----
test <- FALSE
cc <- TRUE

#3. Set nodes for local vs cluster----
if(cc){ nodes <- 32}
if(!cc | test){ nodes <- 2}

#4. Create nodes list----
print("* Creating nodes list *")

#For running on cluster
nodeslist <- unlist(strsplit(Sys.getenv("NODESLIST"), split=" "))

#For testing on local
if(!cc){ nodeslist <- nodes }

print(nodeslist)

#5. Create and register clusters----
print("* Creating clusters *")
cl <- makePSOCKcluster(nodeslist, type="PSOCK")

#6. Set root path----
print("* Setting root file path *")
if(cc){root <- "/scratch/ecknight"}
if(!cc){root <- "G:/Shared drives/BAM_NationalModels5"}

tmpcl <- clusterExport(cl, c("root"))

#7. Load packages on clusters----
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(sf))
tmpcl <- clusterEvalQ(cl, library(tidyverse))
tmpcl <- clusterEvalQ(cl, library(terra))

#8. Load data package----
print("* Loading data on master *")

load(file.path(root, "data", "04_NM5.0_data_stratify.R"))
rm(bird, covlist, bcrlist, offsets, visit, gridlist)

#9. Get list of covariates to include----
cov_clean <- cov |> 
  dplyr::select(where(is.numeric))
cov_clean <- cov_clean[names(cov_clean)!="hli3cl_1km"] #remove hli - it is categorical

#10. Set crs----
#NAD83(NSRS2007)/Conus Albers projection (epsg:5072)
crs <- "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"

#11. Load data objects----
print("* Loading data on workers *")

tmpcl <- clusterExport(cl, c("birdlist", "crs", "cov_clean"))

#12. Load the dsmextra functions----
print("* Loading dsmextra *")
if(cc){source("/home/ecknight/NationalModels/dsmextra_fn.R")
  tmpcl <- clusterExport(cl, c("addLegend_decreasing", "check_crs", "compute_extrapolation", "compute_nearby", "ExDet", "make_X", "map_extrapolation", "n_and_p", "proj_rasters", "rescale_cov", "summarise_extrapolation", "whatif", "whatif.opt"))}

if(!cc){library(dsmextra)
  tmpcl <- clusterEvalQ(cl, library(dsmextra))}

#WRITE FUNCTION##########

calc_extrapolation <- function(i){
  
  #1. Get model settings---
  bcr.i <- loop$bcr[i]
  spp.i <- loop$spp[i]
  boot.i <- loop$boot[i]
  year.i <- loop$year[i]
  
  #2. Load model----
  load(file.path(root, "output", "bootstraps", paste0(spp.i, "_", bcr.i, "_", boot.i, ".R")))
  
  #3. Get list of covariates used----
  Variables <- b.i$var.names[b.i$var.names %in% names(cov_clean)]
  
  #4. Get the dataframe of prediction raster values----
  stack.i <- rast(file.path(root, "gis", "stacks", paste0(bcr.i, "_", year.i, ".tif")))
  Dataframes <- as.data.frame(stack.i, xy=TRUE) 
  
  #10. Get training data for that bootstrap----
  sample <- cov_clean |> 
    dplyr::filter(id %in% visit.i$id) |> 
    dplyr::select(all_of(Variables))
  
  #11. Create target dataframes----
  target <- Dataframes |> 
    dplyr::select(c("x", "y", all_of(Variables))) 
  
  #12. Compute extrapolation----
  Extrapol <- compute_extrapolation(samples = sample,
                                    covariate.names = names(sample),
                                    prediction.grid = target,
                                    coordinate.system = crs,
                                    verbose=F)
  
  #13. Produce binary raster of extrapolation area----
  raster <- as(Extrapol$rasters$mic$all, "SpatRaster")
  raster[raster>0] <- 1  # extrapolated locations in each boot
  raster[raster!=1|is.na(raster)] <- 0 #eliminate NA areas
  
  #14. Write raster----
  writeRaster(raster, file.path(root, "output", "extrapolation", paste0(spp.i, "_", bcr.i, "_", boot.i, "_", year.i, ".tiff")), overwrite=T)
  
}

#15. Export to clusters----
print("* Loading function on workers *")

tmpcl <- clusterExport(cl, c("calc_extrapolation"))


#RUN EXTRAPOLATION###############

#1. Set desired years----
years <- seq(1985, 2020, 5)

#2. Get list of models that are bootstrapped----
booted <- data.frame(path = list.files(file.path(root, "output", "bootstraps"), pattern="*.R", full.names=TRUE),
                     file = list.files(file.path(root, "output", "bootstraps"), pattern="*.R")) |> 
  separate(file, into=c("spp", "bcr", "boot"), sep="_", remove=FALSE) |> 
  mutate(boot = as.numeric(str_sub(boot, -100, -3)))

#3. Create to do list----
#Sort longest to shortest duration to get the big models going first
todo <- booted |> 
  dplyr::select(bcr, spp, boot) |> 
  expand_grid(year=years)

#4. Determine which are already done----
done <- data.frame(path = list.files(file.path(root, "output", "extrapolation"), pattern="*.tiff", full.names=TRUE),
                   file = list.files(file.path(root, "output", "extrapolation"), pattern="*.tiff")) |> 
  separate(file, into=c("spp", "bcr", "boot", "year"), sep="_", remove=FALSE) |>  
  mutate(year = as.numeric(str_sub(year, -100, -6)),
         boot = as.numeric(boot))

#5. Create final to do list----
if(nrow(done) > 0){
  
  loop <- todo |> 
    anti_join(done)
  
} else { loop <- todo }

#For testing - take the shortest duration models
if(test) {loop <- loop[1:2,]}

print("* Loading model loop on workers *")
tmpcl <- clusterExport(cl, c("loop"))

#6. Run function in parallel----
print("* Calculating extrapolation *")
mods <- parLapply(cl,
                  X=1:nrow(loop),
                  fun=calc_extrapolation)

#CONCLUDE####

#1. Close clusters----
print("* Shutting down clusters *")
stopCluster(cl)

if(cc){ q() }