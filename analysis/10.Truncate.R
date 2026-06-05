# ---
# title: National Models 5.0 - truncate predictions
# author: Elly Knight
# created: June 4, 2026
# ---

#NOTES################################

# This script tidies up predictions to remove extreme values and overprediction via four steps.

# 1) Truncation of extreme values. This is a two-stage trunation using 1) species-specific densities based on the counts in the raw data using the output of gis/packaging/TruncationValues.R, and 2) a 99.9% quantile.

# 2) Range limitation is done using the output of `gis/packaging/RangeLimitation.R`

# 3) The data limit based on the output of 'gis/packaging/DataExtent.R'

# 4) Each averaged raster is also masked by water bodies, given that we are modelling landbirds.

# This script also transforms from 5072 to 3978 due to changes in projection decisions part way through the modelling process. Future versions will not require this step because they should use 3978 from the onset.

#PREAMBLE############################

#1. Load packages----
print("* Loading packages on master *")
library(tidyverse)
library(terra)
library(sf)
library(parallel)

#2. Determine if on local or cluster----
cc <- TRUE

#3. Set nodes for local vs cluster----
if(cc){ cores <- 32}
if(!cc){ cores <- 6}

#4. Create and register clusters----
print("* Creating clusters *")
print(table(cores))
cl <- makePSOCKcluster(cores, type="PSOCK")
length(clusterCall(cl, function() Sys.info()[c("nodename", "machine")]))

#5. Set root path----
print("* Setting root file path *")
if(cc){root <- "/scratch/ecknight/NationalModels"}
if(!cc){root <- "G:/Shared drives/BAM_NationalModels5"}

#6. Get the water layer----
print("* Getting water layer *")
water_ca <- read_sf(file.path(root, "gis", "WaterMask_Canada.shp"))
water_us <- read_sf(file.path(root, "gis", "WaterMask_US.shp"))

#7. Subunit polygons----
print("* Getting bcrs *")
bcr.country <- read_sf(file.path(root, "gis", "Subregions_Mosaics_EPSG3978.shp"))

#8. Load packages on clusters----
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(sf))
tmpcl <- clusterEvalQ(cl, library(tidyverse))
tmpcl <- clusterEvalQ(cl, library(terra))

#9. Data limit mask ----
limit <- read_sf(file.path(root, "gis", "DataLimitationsMask.shp")) |> 
  st_transform(3978)

#10. Truncation values ----
load(file.path(root, "data", "SpeciesPredictionTruncationValues.Rdata"))

#FUNCTION###########

#1. Set up the loop----
brt_truncate <- function(i){
  
  spp.i <- loop$spp[i]
  year.i <- loop$year[i]
  bcr.i <- loop$bcr[i]
  
  #2. Get the BCR boundary & correct water mask ----
  sf.i <- bcr.country |> 
    dplyr::filter(bcr==bcr.i)
  
  if(bcr.i=="Canada" | str_sub(bcr.i, 1, 3)=="can"){
    water <- water_ca
  } else {water <- water_us}
  
  #3. Get list of predictions, extrapolations, sampling layers----
  if(!bcr.i %in% c("Canada", "Alaska", "Lower48")){
    
    files.i <- sampled |> 
      dplyr::filter(spp==spp.i,
                    year==year.i,
                    bcr==bcr.i) |> 
      rename(samplepath = path) |> 
      mutate(predpath = file.path(root, "output", "07_predictions", spp,
                                  paste0(spp, "_", bcr.i, "_", year, ".tif")))
    
  } else {
    
    files.i <- sampled |> 
      dplyr::filter(spp==spp.i,
                    year==year.i,
                    bcr==bcr.i) |> 
      rename(samplepath = path) |> 
      mutate(predpath = file.path(root, "output", "08_mosaics", bcr.i, spp,
                                  paste0(spp, "_", year, ".tif")))
    
  }
  
  #4. Read in the predictions----
  rast.i <- try(rast(files.i$predpath)) |> 
    terra::project("EPSG:3978", res=1000)
  
  if(inherits(rast.i, "try-error")){return(NULL)}
  
  #5. Truncate by density quantile ----
  
  #Get truncation values
  qsp <- q.out[q.out$spp==spp.i,]$densmax
  
  #truncate
  truncate.i <- clamp(rast.i, upper = qsp, values=TRUE)
  #rm(rast.i)
  
  #6. Secondary upper truncation ----
  
  #take the mean first to get a global value
  mn.i <- app(truncate.i, mean, na.rm=TRUE)
  
  q99 <- global(mn.i, quantile, probs=0.999, na.rm=TRUE)[1,1]
  if(is.na(q99)){return(NULL)}
  
  truncate2.i <- clamp(truncate.i, upper=q99, values=TRUE)
  
  #7. Mask outside range----
  range.i <- rast(file.path(root, "gis", "ranges", paste0(spp.i, ".tif"))) |>
    resample(truncate2.i)
  
  mask.i <- truncate2.i * range.i
  mask.i[is.na(mask.i)] <- 0
  
  #8. Mask by water and NA layer ----
  out.i <- mask.i |> 
    crop(vect(sf.i), mask=TRUE) |> 
    crop(vect(limit), mask=TRUE) |> 
    mask(vect(water), inverse=TRUE)
  
  #9. Add some attributes----
  attr(out.i, "species") <- spp.i
  attr(out.i, "subunit") <- bcr.i
  attr(out.i, "year") <- year.i
  
  #10. Make folders as needed-----
  if(!(file.exists(file.path(root, "output", "10_truncated", spp.i)))){
    dir.create(file.path(root, "output", "10_truncated", spp.i))
  }
  
  if(!(file.exists(file.path(root, "output", "10_truncated", spp.i, bcr.i)))){
    dir.create(file.path(root, "output", "10_truncated", spp.i, bcr.i))
  }
  
  #11. Save----
  writeRaster(out.i, filename = file.path(root, "output", "10_truncated", spp.i, bcr.i, paste0(spp.i, "_", bcr.i, "_", year.i, ".tif")), overwrite=TRUE, datatype = "FLT4S")
  
}

#INVENTORY#########

#1. Get the list of sampling layers----
#remove 1985 - we're not providing those predictions
sampled <- data.frame(file = list.files(file.path(root, "output", "09_sampling"), pattern="*.tif", recursive = TRUE)) |>
  separate(file, into=c("spf", "spp", "bcr", "year", "filetype"), remove=FALSE) |>
  mutate(year = as.numeric(year),
         path = file.path(root, "output", "09_sampling", file)) |>
  dplyr::filter(year!=1985) |>
  dplyr::select(-filetype, -file, -spf)

#2. Check which have been run----
done <- data.frame(file = list.files(file.path(root, "output", "10_truncated"), pattern="*.tif", recursive=TRUE))  |> 
  separate(file, into=c("sppfolder", "bcrfolder", "spp", "bcr", "year", "filetype"), remove=FALSE) |> 
  mutate(year = as.numeric(year)) |> 
  dplyr::select(-filetype)

#3. Make the todo list----
#remove species that we are omitting for now
loop <- sampled |> 
  anti_join(done) |> 
  arrange(-year, spp, bcr)

#4. Shut down if nothing left to do----
if(nrow(loop)==0){
  print("* Shutting down clusters *")
  stopCluster(cl)
  
  if(cc){ q() }
}

#PACKAGE########

#1. Export objects to clusters----
tmpcl <- clusterExport(cl, c("loop", "sampled", "bcr.country", "root", "water_ca", "water_us", "limit", "q.out"))

#2. Run BRT function in parallel----
print("* Truncating *")
truncated <- parLapply(cl,
                      X=1:nrow(loop),
                      fun=brt_truncate)

#CONCLUDE####

#1. Close clusters----
print("* Shutting down clusters *")
stopCluster(cl)

if(cc){ q() }