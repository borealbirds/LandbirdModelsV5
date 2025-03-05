# ---
# title: National Models 5.0 - plot sampling density
# author: Elly Knight, adapted from Anna Drake
# created: December 23, 2024
# ---

#NOTES################################

# This script uses the test data from each bootstrap to visualize the detections included in the models.

# The script works by using the test data saved from the model validation script (`10.Validate.R`).

#PREAMBLE############################

#1. Load packages----
library(tidyverse)
library(terra)
library(sf)

#2. Set root path----
root <- "G:/Shared drives/BAM_NationalModels5"

#3. Subunit polygons----
bcr.country <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel_Unbuffered.shp")) |> 
  mutate(bcr= paste0(country, subUnit)) |> 
  st_transform(crs="ESRI:102001")

#4. Mosaic polygon----
bcr.all <- st_union(bcr.country)

#BLANK RASTERS#######################

#1. Set up loop----
rasts <- list()
for(i in 1:nrow(bcr.country)){
  
  #2. Load the raster stack----
  rast.i <- rast(file.path(root, "gis", "stacks", paste0(bcr.country$bcr[i], "_2000.tif")))
  
  #3. Get the BCR boundary----
  sf.i <- bcr.country |> 
    dplyr::filter(bcr==bcr.country$bcr[i]) |> 
    st_transform(crs=crs(rast.i)) |> 
    vect()
  
  #4. Crop to BCR boundary----
  rasts[[i]] <- rast.i$method |> 
    crop(sf.i, mask=TRUE)
  
  #5. Set non NAs to zero----
  rasts[[i]] <- ifel(!is.na(rasts[[i]]), 0, NA)
  
  cat("Finished raster", i, "of", nrow(bcr.country), "\n")
    
}

#6. Names----
names(rasts) <- bcr.country$bcr

#7. Load a mosaic raster----
rast.mosaic <- rast(file.path(root, "output", "mosaics", "predictions", "OVEN_1_2020.tiff"))

#8. Add it to the list----
rasts[["mosaic"]] <- ifel(!is.na(rast.mosaic), 0, NA)

#INVENTORY#########

#1. Get list of sets that have been validated----
#We use the validated datasets because we can pull the training data out of the output easily
validated <- data.frame(path = list.files(file.path(root, "output", "validation"), pattern="*.Rdata", full.names=TRUE),
                        file = list.files(file.path(root, "output", "validation"), pattern="*.Rdata")) |> 
  separate(file, into=c("spp", "boot", "filetype"), remove=FALSE) |>  
  mutate(boot = as.numeric(boot))

#2. Make the todo list----
todo <- validated |> 
  dplyr::select(spp, boot) |> 
  unique()

#3. Check which have been run----
done <- data.frame(file.mean = list.files(file.path(root, "output", "sampling"), pattern="*.tiff", recursive=TRUE)) |> 
  separate(file.mean, into=c("folder", "spp", "bcr", "boot", "year", "filetype"), remove=FALSE) |> 
  mutate(boot = as.numeric(boot)) |> 
  dplyr::filter(bcr=="mosaic") |> 
  dplyr::select(spp, boot) |> 
  unique()

#4. Make the todo list----
loop <- anti_join(todo, done)

#PLOT SAMPLING DENSITY###########

#1. Set up loop----
for(i in 1:nrow(loop)){
  
  #2. Get the loop settings----
  spp.i <- loop$spp[i]
  boot.i <- loop$boot[i]
  
  #3. Load the validation object----
  load(file.path(root, "output", "validation",
                 paste0(spp.i, "_", boot.i, ".Rdata")))
  
  #4. Set up the bootstrap loop----
  for(j in 1:length(train)){
    
    #5. Get the bcr----
    bcr.j <- names(train)[j]
    
    #6. Get the blank raster----
    rast.j <- rasts[[bcr.j]]
    
    #7. Get the BCR boundary----
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
    
    #8. Get the training data----
    train.j <- train[[j]]
    
    #9. Round to nearest year and take detections only----
    detections.j <- train.j |> 
      rename(year1 = year) |> 
      mutate(year = round(year1/5)*5,
             year = ifelse(year==1980, 1985, year)) |> 
      dplyr::filter(count > 0)
    
    #10. Set up year loop----
    years <- seq(1990, 2020, 5)
    for(k in 1:length(years)){
      
      #11. Get the year----
      year.k <- years[k]
      
      #12. Filter the detections to year----
      detections.k <- detections.j |> 
        dplyr::filter(year==year.k)
      
      #13. If there's no detections, make a blank raster----
      if(nrow(detections.k)==0){
        
        rast.out <- rasts[[bcr.j]]
        
      } else {
        
        #14. Convert points to a spatial vector----
        pts.k <- detections.k |> 
          st_as_sf(coords=c("lon", "lat"), crs=4326) |> 
          st_transform(crs(rast.j)) |> 
          vect()
        
        #15. Caculate distance to nearest----
        rast.out <- distance(rast.j, pts.k, unit="km") |> 
          crop(sf.j, mask=TRUE)
        
      }
      
      #16. Save raster----
      if(!(file.exists(file.path(root, "output", "sampling", spp.i)))){
        dir.create(file.path(root, "output", "sampling", spp.i))
      }
      
      writeRaster(rast.out, filename = file.path(root, "output", "sampling", spp.i, paste0(spp.i, "_", bcr.j, "_", boot.i, "_", year.k, ".tiff")), overwrite=TRUE)
      
      cat("Year", year.k, "of bcr", j, "of", length(train), "\n")
      
    }
    
  }
  
  cat("FINISHED BOOTSTRAP", i, "OF", nrow(loop), "\n")
  
}
