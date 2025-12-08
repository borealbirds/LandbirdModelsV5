# ---
# title: National Models 5.0 - plot sampling density
# author: Elly Knight, adapted from Anna Drake
# created: December 23, 2024
# ---

#NOTES################################

# This script uses the test data from each bootstrap to average distance to the nearest detection in the models. The mean output from this script is provided in the packaged raster stacks to inform interpretation of the models.

#PREAMBLE############################

#1. Load packages----
print("* Loading packages on master *")
library(tidyverse)
library(terra)
library(sf)
library(parallel)
library(Matrix)

#2. Determine if testing and on local or cluster----
cc <- FALSE

#3. Set nodes for local vs cluster----
if(cc){ cores <- 32}
if(!cc){ cores <- 4}

#4. Create and register clusters----
print("* Creating clusters *")
print(table(cores))
cl <- makePSOCKcluster(cores, type="PSOCK")
length(clusterCall(cl, function() Sys.info()[c("nodename", "machine")]))

#5. Set root path----
print("* Setting root file path *")
if(cc){root <- "/scratch/ecknight/NationalModels"}
if(!cc){root <- "G:/Shared drives/BAM_NationalModels5"}

#6. Subunit polygons----
bcr.country <- read_sf(file.path(root, "gis", "Subregions_unbuffered.shp"))

#7. Names of bcrs ----
bcrs <- bcr.country$bcr

#8. Mosaic polygon----
bcr.all <- st_union(bcr.country)

#9. Data package ----
load(file.path(root, "data", "04_NM5.0_data_stratify.Rdata"))
rm(cov, covlist, offsets)

#10. Load packages on clusters----
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(sf))
tmpcl <- clusterEvalQ(cl, library(tidyverse))
tmpcl <- clusterEvalQ(cl, library(terra))
tmpcl <- clusterEvalQ(cl, library(Matrix))

#FUNCTION####

#1. Set up loop----
brt_sampling <- function(i){
  
  #4. Get the region ----
  bcr.i <- loop$region[i]
  spp.i <- loop$spp[i]
  year.i <- loop$year[i]
  
  #5. Get the template raster ----
  rast.i <- rast(file.path(root, "gis", "zeros", paste0(bcr.i, ".tif")))
  
  #6. Set up bootstrap list ----
  for(k in 1:32){
    
    #7. Get visits to include----
    if(!bcr.i %in% c("Canada", "Alaska", "Lower48")){
      
      visit.k <- bcrlist[bcrlist[,bcr.i]>0, c("id", bcr.i)] |> 
        dplyr::filter(id %in% bootlist[[k + 2]])
      
    } else {
      
      if(bcr.i=="Canada"){bcrs.i <- bcrs[str_sub(bcrs, 1, 3)=="can"]}
      
      if(bcr.i=="Alaska"){bcrs.i <- c("usa41423", "usa2", "usa40", "usa43", "usa5")}
      
      if(bcr.i=="Lower48"){bcrs.i <- c("usa5", "usa9", "usa10", "usa11", "usa13", "usa14", "usa23", "usa28")}
      
      visit.k <- bcrlist |> 
        dplyr::select(all_of(c("id", bcrs.i))) |> 
        pivot_longer(-id, names_to="bcr", values_to="use") |> 
        dplyr::filter(use > 0) |> 
        dplyr::select(id) |> 
        unique() |> 
        arrange(id)
      }
    
    #8. Get year & coords ----
    year.k <- visit[visit$id %in% visit.k$id, c("year", "lon", "lat")]
    
    #9. Get response data (bird data)----
    bird.k <- bird[as.character(visit.k$id), spp.i]
    
    #10. Filter ----
    detections.k <- cbind(bird.k, year.k) |> 
      rename(count = bird.k,
             year1 = year) |> 
      mutate(year = round(year1/5)*5,
             year = ifelse(year==1980, 1985, year)) |> 
      dplyr::filter(count > 0,
                    year==year.i)
    
    #11. If there's no detections, make a blank raster----
    if(nrow(detections.k)==0){
      
      rast.k <- rast.i
      
    } else {
      
      #12. Convert points to a spatial vector----
      pts.k <- detections.k |> 
        st_as_sf(coords=c("lon", "lat"), crs=4326) |> 
        st_transform(crs(rast.i)) |> 
        vect()
      
      #13. Caculate distance to nearest----
      rast.k <- distance(rast.i, pts.k, unit="km") |> 
        crop(rast.i, mask=TRUE)
      
    }
    
    #14. Make a stack ----
    if(k==1){ rast.out <- rast.k} else { rast.out <- c(rast.out, rast.k) }
    
  }
  
  #15. Name the rasters ----
  names(rast.out) <- paste0("b", seq(1:32))
  
  #16. Save raster----
  if(!(file.exists(file.path(root, "output", "09_sampling", spp.i)))){
    dir.create(file.path(root, "output", "09_sampling", spp.i))
  }
  
  writeRaster(rast.out, filename = file.path(root, "output", "09_sampling", spp.i, paste0(spp.i, "_", bcr.i, "_", year.i, ".tif")), overwrite=TRUE)
  
  cat("Sampling", i, "of", nrow(loop), "\n")
  
}

#INVENTORY#########

#1. Get list of sets that have been mosaiced----
mosaiced <- data.frame(file = list.files(file.path(root, "output", "08_mosaics"), pattern="*.tif", recursive=TRUE)) |> 
  separate(file, into=c("region", "sppfolder", "spp", "year", "filetype"), remove=FALSE) |>  
  mutate(year = as.numeric(year)) |> 
  dplyr::select(-sppfolder, -filetype)

#2. Get the individual models ----
predicted <- data.frame(file = list.files(file.path(root, "output", "07_predictions"), pattern="*.tif", recursive=TRUE)) |> 
  separate(file, into=c("sppfolder", "spp",  "region", "year", "filetype"), remove=FALSE) |>  
  mutate(year = as.numeric(year)) |> 
  dplyr::select(-sppfolder, -filetype)

#3. Make the to-do list----
todo <- rbind(mosaiced, predicted)

#4. Check which have been run----
done <- data.frame(file.mean = list.files(file.path(root, "output", "09_sampling"), pattern="*.tif", recursive=TRUE)) |> 
  separate(file.mean, into=c("sppfolder", "spp", "region", "year", "filetype"), remove=FALSE) |> 
  mutate(year = as.numeric(year))

#4. Make the todo list----
loop <- anti_join(todo, done)

#RUN#####

#1. Export objects to clusters----
tmpcl <- clusterExport(cl, c("loop", "bcr.all", "bcr.country", "brt_sampling", "root", "visit", "bird", "bcrlist", "birdlist", "bootlist", "bcrs"))

#2. Run BRT function in parallel----
print("* Calculating sampling *")
sampling <- parLapply(cl,
                     X=1:nrow(loop),
                     fun=brt_sampling)
