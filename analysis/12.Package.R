# ---
# title: National Models 5.0 - calculate means
# author: Elly Knight
# created: April 29, 2024
# ---

#NOTES################################

# This script calculates means and coefficients of variation (CV) across bootstraps of model predictions for packaging

# Means and cvs are calculated for two extents:
#1. Subunit: this the output of `07.Predict.R` run on compute canada, with the predictions trimmed to the BCR boundary.
#2. Range: this is the output of `09.MosaicPredictions.R`

#This script also standardizes all projects to ESRI:102001 (NAD83 Canada Alberta Equal Conic)

# Each averaged raster is also masked by water bodies

#TO DO: PARALLELIZE

#TO DO: CONSIDER ADDING RANGE TRUNCATION TO MOSAIC PRODUCTS

#TO DO: ADD ZEROS TO THE REST OF THE CV LAYER?

#TO DO: ADD ATTRIBUTES FOR SP ETC

#NOTES FROM PACKAGING SCRIPT:
#LIST OF PRODUCTS PER SPECIES:
#2 extents: boreal, regional
#Products in both: mean, cv, extrapolation, sampling

#LIST OF SAMPLING PRODUCT:
#All the extrapolation
#Sampling density
#Distance to sampling

#PREAMBLE############################

#1. Load packages----
library(tidyverse)
library(terra)
library(sf)

#2. Set root path----
root <- "G:/Shared drives/BAM_NationalModels5"

#3. Get the water layer----
water <- read_sf(file.path(root, "Regions", "Lakes_and_Rivers", "hydrography_p_lakes_v2.shp")) |> 
  dplyr::filter(TYPE %in% c(16, 18)) |> 
  st_transform(crs="ESRI:102001") |> 
  vect()

#4. Subunit polygons----
bcr.country <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel_Unbuffered.shp")) |> 
  mutate(bcr= paste0(country, subUnit)) |> 
  st_transform(crs="ESRI:102001")

#5. Mosaic polygon----
bcr.all <- st_union(bcr.country)

#6. Data package----
load(file.path(root, "data", "04_NM5.0_data_stratify.R"))

#INVENTORY#########

#1. Get list of sets that have been mosaicked----
#We use the mosaics as input because we know they have all subunits, ten bootstraps, and the mosaics and extrapolation complete
mosaics <- data.frame(path = list.files(file.path(root, "output", "mosaics", "predictions"), pattern="*.tiff", full.names=TRUE, recursive = TRUE),
                        file = list.files(file.path(root, "output", "mosaics", "predictions"), pattern="*.tiff", recursive = TRUE)) |> 
  separate(file, into=c("spp", "boot", "year", "filetype"), remove=FALSE) |>  
  mutate(year = as.numeric(year),
         boot = as.numeric(boot))

#2. Make the todo list----
todo <- mosaics |> 
  dplyr::select(spp, year) |> 
  unique()

#3. Check which have been run----
done <- data.frame(file = list.files(file.path(root, "output", "packaged"), pattern="*.tiff", recursive=TRUE))  |> 
  separate(file, into=c("sppfolder", "bcrfolder", "spp", "bcr", "year", "filetype"), remove=FALSE) |> 
  mutate(year = as.numeric(year))

#4. Make the todo list----
loop <- anti_join(todo, done |> 
                    dplyr::filter(bcr=="mosaic") |> 
                    dplyr::select(spp, year)) 

#AVERAGE MOSAICS###########

#1. Set up the loop----
for(i in 1:nrow(loop)){
  
  spp.i <- loop$spp[i]
  year.i <- loop$year[i]
  
  #2. Make the list of BCRs-----
  loop.i <- birdlist[,c("bcr", spp.i)] |> 
    pivot_longer(-bcr, names_to="spp", values_to="use") |> 
    dplyr::filter(use==TRUE) |> 
    dplyr::select(bcr, use) |> 
    rbind(data.frame(bcr="mosaic", use=TRUE))
  
  #3. Set up BCR loop----
  for(j in 1:nrow(loop.i)){
    
    #4. Get the bcr----
    bcr.j <- loop.i$bcr[j]
    
    #5. Get the BCR boundary----
    if(bcr.j!="mosaic"){
      sf.j <- bcr.country |> 
        dplyr::filter(bcr==bcr.j) |> 
        st_transform("ESRI:102001") |> 
        vect()
    } else {
      sf.j <- bcr.all |> 
        st_transform("ESRI:102001") |> 
        vect()
    }
    
    #6. Get list of predictions, extrapolations, sampling layers----
    if(bcr.j!="mosaic"){
      
      files.j <- mosaics |> 
        dplyr::filter(spp==spp.i,
                      year==year.i) |> 
        mutate(predpath = file.path(root, "output", "predictions", spp,
                                    paste0(spp, "_", bcr.j, "_", boot, "_", year, ".tiff")),
               extrappath = file.path(root, "output", "extrapolation", spp,
                                      paste0(spp, "_", bcr.j, "_", boot, "_", year, ".tiff")),
               samplepath = file.path(root, "output", "sampling", spp, 
                                      paste0(spp, "_", bcr.j, "_", boot, "_", year, ".tiff")))
      
    } else {
      
      files.j <- mosaics |> 
        dplyr::filter(spp==spp.i,
                      year==year.i) |> 
        rename(predpath = path) |> 
        mutate(extrappath = file.path(root, "output", "mosaics", "extrapolation", file),
               samplepath = file.path(root, "output", "sampling", spp, 
                                      paste0(spp, "_", bcr.j, "_", boot, "_", year, ".tiff")))
      
    }
    
    #7. Read in the predictions----
    stack.j <- try(rast(files.j$predpath) |> 
                     project("ESRI:102001"))
    
    #8. Truncate to 99% quantile----
    q99 <- global(stack.j, fun=function(x) quantile(x, 0.99, na.rm=TRUE))
    
    truncate.j <- stack.j
    for(k in 1:nlyr(truncate.j)){
      truncate.j[[k]][values(truncate.j[[k]]) > q99[k,1]] <- q99[k,1]
    }
    
    #9. Calculate mean----
    mean.j <- mean(truncate.j, na.rm=TRUE) |> 
      crop(sf.j, mask=TRUE) |> 
      mask(water, inverse=TRUE)
    
    #10. Calculate sd----
    sd.j <- stdev(truncate.j, na.rm=TRUE) |> 
      crop(sf.j, mask=TRUE) |> 
      mask(water, inverse=TRUE)
    
    #11. Calculate cv----
    cv.j <- sd.j/mean.j
    
    #12. Read in the extrapolations----
    extrap.j <- try(rast(files.j$extrappath) |> 
                      project("ESRI:102001"))
    
    #13. Calculate mean extrapolation----
    extrapmn.j <- mean(extrap.j, na.rm=TRUE) |> 
      resample(mean.j) |>  
      crop(sf.j, mask=TRUE) |> 
      mask(water, inverse=TRUE)
    
    #14. Read in the sampling distance layers----
    #TO DO: REMOVE SUBSAMPLE LINE############
    sample.j <- try(rast(files.j$samplepath[c(3,4)]) |> 
                      project("ESRI:102001"))
    
    #15. Calculate mean sampling distance----
    samplemn.j <- mean(sample.j, na.rm=TRUE) |> 
      resample(mean.j) |>  
      crop(sf.j, mask=TRUE) |> 
      mask(water, inverse=TRUE)
    
    #16. Stack----
    out.i <- c(mean.j, cv.j, extrapmn.j, samplemn.j)
    names(out.i) <- c("mean", "cv", "extrapolation", "detections")
    
    #17. Add some attributes----
    attr(out.i, "species") <- spp.i
    attr(out.i, "subunit") <- bcr.j
    attr(out.i, "year") <- year.i
    
    #18. Make folders as needed-----
    if(!(file.exists(file.path(root, "output", "packaged", spp.i)))){
      dir.create(file.path(root, "output", "packaged", spp.i))
    }
    
    if(!(file.exists(file.path(root, "output", "packaged", spp.i, bcr.j)))){
      dir.create(file.path(root, "output", "packaged", spp.i, bcr.j))
    }
    
    #19. Save----
    writeRaster(out.i, filename = file.path(root, "output", "packaged", spp.i, bcr.j, paste0(spp.i, "_", bcr.j, "_", year.i, ".tiff")), overwrite=TRUE)
    
    cat("Finished subunit", j, "of", nrow(loop.i), "for", spp.i, "in", year.i, "\n")
  
  }
  
  cat("FINISHED", i, "of", nrow(loop), "PREDICTION SETS\n")
  
}
