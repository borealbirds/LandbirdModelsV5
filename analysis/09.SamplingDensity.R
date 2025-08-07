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
if(!cc){ cores <- 3}

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

#7. Mosaic polygon----
bcr.all <- st_union(bcr.country)

#8. Data package ----
load(file.path(root, "data", "04_NM5.0_data_stratify.Rdata"))
rm(cov, covlist, offsets)

#9. Load packages on clusters----
print("* Loading packages on workers *")
tmpcl <- clusterEvalQ(cl, library(sf))
tmpcl <- clusterEvalQ(cl, library(tidyverse))
tmpcl <- clusterEvalQ(cl, library(terra))
tmpcl <- clusterEvalQ(cl, library(Matrix))

#FUNCTION####

#1. Set up loop----
brt_sampling <- function(i){
  
  #2. Get the loop settings----
  spp.i <- loop$spp[i]
  year.i <- loop$year[i]
  
  #3. Get the list of subunits----
  bcr.i <- birdlist |> 
    pivot_longer(-bcr, names_to="spp", values_to="use") |> 
    dplyr::filter(use==TRUE, spp==spp.i) |> 
    dplyr::select(-use) |> 
    rbind(data.frame(bcr="mosaic", spp=spp.i))
  
  for(j in 1:nrow(bcr.i)){
    
    #4. Get the BCR ----
    bcr.j <- bcr.i$bcr[j]
    
    #5. Get the template raster ----
    rast.j <- rast(file.path(root, "gis", "blanks", paste0(bcr.j, ".tif")))
    
    #5. Get the BCR boundary----
    if(bcr.j!="mosaic"){
      sf.j <- bcr.country |> 
        dplyr::filter(bcr==bcr.j) |> 
        st_transform(crs=crs(rast.j)) |> 
        vect()
    } else {
      sf.j <- bcr.all |> 
        st_transform(crs=crs(rast.j)) |> 
        vect()
    }
    
    #6. Set up bootstrap list ----
    for(k in 1:32){
      
      #7. Get visits to include----
      if(bcr.j!="mosaic"){
        visit.k <- bcrlist[bcrlist[,bcr.j]>0, c("id", bcr.j)] |> 
          dplyr::filter(id %in% bootlist[[k + 2]])
      } else {
        bcrs <- bcr.i[bcr.i$bcr!="mosaic",]$bcr
        
        visit.k <- bcrlist |> 
          dplyr::select(all_of(c("id", bcrs))) |> 
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
        
        rast.k <- rast.j
        
      } else {
        
        #12. Convert points to a spatial vector----
        pts.k <- detections.k |> 
          st_as_sf(coords=c("lon", "lat"), crs=4326) |> 
          st_transform(crs(rast.j)) |> 
          vect()
        
        #13. Caculate distance to nearest----
        rast.k <- distance(rast.j, pts.k, unit="km") |> 
          crop(sf.j, mask=TRUE)
      
      }
      
      #14. Make a stack ----
      if(k==1){ rast.out <- rast.k} else { rast.out <- c(rast.out, rast.k) }
    
    }
    
    #15. Name the rasters ----
    names(rast.out) <- paste0("b", seq(1:32))
    
    #16. Drop the sources----
    rast.save <- rast(rast.out[])
    
    #16. Save raster----
    if(!(file.exists(file.path(root, "output", "09_sampling", spp.i)))){
      dir.create(file.path(root, "output", "09_sampling", spp.i))
    }
    
    writeRaster(rast.out, filename = file.path(root, "output", "09_sampling", spp.i, paste0(spp.i, "_", bcr.j, "_", year.i, ".tif")), overwrite=TRUE)
    
    cat("BCR ", j, "of", nrow(bcr.i), "\n")
    
  }
  
  cat("Species and year", i, "of", nrow(loop), "\n")
  
}

#INVENTORY#########

#1. Get list of sets that have been mosaiced----
mosaiced <- data.frame(path = list.files(file.path(root, "output", "08_mosaics"), pattern="*.tif", full.names=TRUE),
                       file = list.files(file.path(root, "output", "08_mosaics"), pattern="*.tif")) |> 
  separate(file, into=c("spp", "year", "filetype"), remove=FALSE) |>  
  mutate(year = as.numeric(year))

#2. Check which have been run----
done <- data.frame(file.mean = list.files(file.path(root, "output", "09_sampling"), pattern="*.tif", recursive=TRUE)) |> 
  separate(file.mean, into=c("folder", "spp", "bcr", "year", "filetype"), remove=FALSE) |> 
  mutate(year = as.numeric(year)) |> 
  dplyr::filter(bcr=="mosaic")

#4. Make the todo list----
loop <- anti_join(mosaiced, done)

#RUN#####

#1. Export objects to clusters----
tmpcl <- clusterExport(cl, c("loop", "bcr.all", "bcr.country", "brt_sampling", "root", "visit", "bird", "bcrlist", "birdlist", "bootlist"))

#2. Run BRT function in parallel----
print("* Calculating sampling *")
sampling <- parLapply(cl,
                     X=1:nrow(loop),
                     fun=brt_sampling)
