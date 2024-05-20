# ---
# title: National Models 5.0 - stack and crop rasters for subunit prediction
# author: Elly Knight
# created: March 5, 2024
# ---

#NOTES################################

#This script creates separate raster stacks for each of the modelling subunits for prediction.

#Next, we identify the covariates that we use for that region, stack them, and then crop them by the polygon for that region. We save out each stack separately.

#Finally, we add the layers to each stack for the nonspatially explicit variables in the model (year, method)

#The stacks can then be moved over to compute canada for model prediction.

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #basic data wrangling
library(terra) #basic raster handling
library(sf) #basic shapefile handling
library(nngeo) #fill holes in polygons
library(tidyterra) #raster plotting

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_NationalModels/NationalModels5.0"

#3. Load data packages with covariate lookup table----
load(file.path(root, "Data", "04_NM5.0_data_stratify.R"))

#4. Set options to not write the .aux.xml file----
rgdal::setCPLConfigOption("GDAL_PAM_ENABLED", "FALSE")

#5. Get buffered region shapefile----
bcr.out <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel_Buffered.shp"))

#GET PREDICTION RASTER LOCATIONS######

#1. Get list of files in the folder----
#Remove Archived files and a couple weirdos, 100m ones
files <- data.frame(path = list.files(file.path(root, "PredictionRasters"), full.names = TRUE, recursive = TRUE, pattern="*.tif"),
                    file = list.files(file.path(root, "PredictionRasters"), full.names = FALSE, recursive = TRUE, pattern="*.tif")) |> 
  dplyr::filter(!str_sub(file, 1, 7)=="Archive") |> 
  rowwise() |> 
  mutate(int = str_locate_all(file, "/"),
         int = max(int)) |> 
  ungroup() |> 
  mutate(file = str_sub(file, int+1, 100)) |> 
  separate(file, into=c("var", "scale", "year", "filetype"), remove=FALSE) |> 
  mutate(year = as.numeric(ifelse(is.na(filetype), NA, year))) |> 
  dplyr::filter(!var %in% c("Standardized", "VLCE100"),
                scale!="100m") |> 
  mutate(cov = paste0(var, "_", scale)) |> 
  unique()

#2. Add the two that are duplicates with temporal mismatch---
add <- files |> 
  dplyr::filter(cov %in% c("ERAPPTsm_1km", "ERATavesm_1km")) |> 
  mutate(cov = gsub(pattern="sm", replacement="smt", x=cov))

files.use <- rbind(files, add) |> 
  arrange(cov, year) 

#3. Check against the covariate list----
missing <- data.frame(list = names(covlist))  |>  
  dplyr::filter(!list %in% c(files.use$cov, "bcr"))
nrow(missing)

#4. Get the t-1 ones----
meth <- readxl::read_excel(file.path(root, "NationalModels_V5_VariableList.xlsx"), sheet = "ExtractionLookup")

lag <- meth |> 
  dplyr::filter(YearMatch==-1)

#5. Match to each desired year of prediction----

#Get the covariate list
covs <- unique(files.use$cov)

#Set the years for prediction
years <- seq(1985, 2020, 5)

year.out <- data.frame()
for(i in 1:length(covs)){
  
  #Get the available years for that covariate
  cov.i <- files.use |> 
    dplyr::filter(cov==covs[i])
  
  #Identify if is a cov that doesn't have annual layers---
  if(is.na(cov.i$year)[1]){
    
    year.out <- rbind(year.out,
                      data.frame(cov = covs[i],
                                 year = NA,
                                 predyear = years))
  }
  
  #Adjust the time lag as needed
  if(covs[i] %in% lag$Label){ cov.i$year <- cov.i$year + 1}
  
  #Match to nearest year of available data
  dt = data.table::data.table(year=unique(cov.i$year), val = unique(cov.i$year))
  data.table::setattr(dt, "sorted", "year")
  data.table::setkey(dt, year)
  year.i <- dt[J(years), roll = "nearest"]$val
  
  #Readjust the time lag as need
  if(covs[i] %in% lag$Label){year.i <- year.i-1}
  
  #Put it together----
  year.out <- rbind(year.out,
                    data.frame(cov = covs[i],
                               year = year.i,
                               predyear = years))
  
}

#6. Put together with file paths----
files.year <- year.out |> 
  left_join(files.use,
            multiple="all") |> 
  unique() |> 
  dplyr::select(cov, predyear, year, path)

summary(files.year)
  
#MAKE STACKS####

#1. Set up loop for subunits * year----
units <- bcr.out |> 
  st_drop_geometry() |> 
  mutate(bcr = paste0(country, subUnit)) |> 
  unique() |> 
  expand_grid(year = seq(1985, 2020, 5)) 

for(i in 1:nrow(units)){
  
  #2. Get subunit----
  bcr.i <- paste0(units$bcr[i])
  year.i <- paste0(units$year[i])
  
  #3. Get covariate list----
  covlist.i <- covlist |> 
    dplyr::filter(bcr==bcr.i) |> 
    pivot_longer(ERAMAP_1km:mTPI_1km, names_to="cov", values_to="use") |> 
    dplyr::filter(use==TRUE)
  
  #4. Get file paths----
  files.i <- files.year |> 
    dplyr::filter(cov %in% covlist.i$cov,
                  predyear == year.i)
  
  #5. Get the bcr shp----
  shp.i <- bcr.out |> 
    dplyr::filter(country==units$country[i],
                  subUnit==units$subUnit[i]) |> 
    vect()
  
  #6. Set up loop for the rasters----
  for(j in 1:nrow(files.i)){
    
    #7. Read it in----
    rast.i <- rast(files.i$path[j])
    
    print("read")
    
    #8. Reproject as needed----
    if(crs(rast.i)!=crs(shp.i)){rast.i <- project(rast.i, crs(shp.i))}
    
    print("reproject")
    
    #9. Crop----
    crop.i <- crop(rast.i, shp.i) |> 
      mask(shp.i)
    names(crop.i) <- files.i$cov[j]
    
    print("crop")
    
    #10. Resample for extent as needed----
    if(j > 1){
      if(ext(crop.i)!=ext(stack.i)){crop.i <- resample(crop.i, stack.i)}
    } 
    
    print("resample")
    
    #11. Stack----
    if(j==1){stack.i <- crop.i} else {stack.i <- c(stack.i, crop.i)}
    
    print("stack")
    
    #12. Sanity check plot the first one----
    if(j==1) {
      
      plot.i <- ggplot() +
        geom_spatraster(data=crop.i) +
        geom_sf(data=shp.i, fill=NA, colour="black", linewidth = 2)
      
      ggsave(plot.i, filename=file.path(root, "stacks", "CheckPlots", paste0(bcr.i, "_", year.i, ".jpeg")), width=6, height=4)
      
    }
    
    #13. Report----
    print(paste0("Finished raster ", j, " of ", nrow(files.i), " for bcr*year ", i, " of ", nrow(units)))
    
  }
  
  #14. Make a new method raster----
  meth.i <- rast(ext(stack.i), resolution=res(stack.i), crs=crs(stack.i))
  values(meth.i) <- "PC"
  
  #15. Make a new year raster----
  year.r <- rast(ext(stack.i), resolution=res(stack.i), crs=crs(stack.i))
  values(year.r) <- units$year[i]
  
  #16. Put together and mask----
  meth.year.i <- c(meth.i, year.r) |> 
    mask(shp.i)
  
  #17. Rename----
  names(meth.year.i) <- c("method", "year")
  
  #18. Restack----
  stack.out <- c(meth.year.i, stack.i)

  #14. Save----
  terra::writeRaster(stack.out, file.path(root, "stacks", paste0(bcr.i, "_", year.i, ".tif")), overwrite=TRUE)
  
  rm(rast.i, stack.i, stack.out, crop.i, shp.i)
  
  #14. Remove temp files to save RAM----
  tmp.dir <- tempdir()
  tmp.files <- list.files(tmp.dir, pattern="*.tif", full.names = TRUE)
  file.remove(tmp.files)
  
}

#15. Check they're all there----
files.stack <- data.frame(file=list.files(file.path(root, "stacks"), pattern="*.tif")) |> 
  separate(file, into=c("bcr", "year", "tif"), remove=FALSE) |>
  mutate(year = as.numeric(year)) |> 
  dplyr::filter(str_sub(file, -3, -1)!="xml")

todo.stack <- anti_join(units, files.stack)
nrow(todo.stack)






#MOVE THIS TO LATER#######
#MASKING####

#1. Get the water layer----
water <- read_sf(file.path(root, "Regions", "Lakes_and_Rivers", "hydrography_p_lakes_v2.shp")) |> 
  dplyr::filter(TYPE %in% c(16, 18)) |> 
  st_transform(crs=crs(bcr.out)) |> 
  vect()

#2. Set up loop----
for(i in 1:nrow(files.stack)){
  
  #3. Get the stack----
  stack.i <- rast(file.path(root, "PredictionRasters", "SubunitStacks", "Unmasked", files.stack$file[i]))
  
  #4. Mask it----
  mask.i <- mask(stack.i, water, inverse=TRUE)
  
  #5. Save it----
  terra::writeRaster(mask.i, file.path(root, "PredictionRasters", "SubunitStacks", "Masked", paste0(files.stack$bcr[i], "_", files.stack$year[i], ".tif")), overwrite=TRUE)
  
  #6. Report----
  print(paste0("Finished raster ", i, " of ", nrow(files.stack)))
  
  
}
