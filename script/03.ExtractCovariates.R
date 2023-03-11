# ---
# title: National Models 4.1 - extract covariates
# author: Elly Knight
# created: December 22, 2022
# ---

#NOTES################################

#This script uses CONUS Albers equal-area conic (EPSG: 5072)

#This script extracts the list of model covariates agreed upon by the BAM National Model team. The full list of those covariates can be found at https://docs.google.com/spreadsheets/d/1XATuq8BOYC2KkbJObaturaf4NFD2ufvn/edit?usp=sharing&ouid=104837701987164094932&rtpof=true&sd=true

#This script extracts covariates for every combination of location & year. Efficiency could be improved by 15-20% by taking only the unique locations for all rows in the method table with TemporalResolution=="static"

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(rgee) #to extract data from google earth engine
library(terra) #basic raster handling
library(sf) #basic shapefile handling
library(exactextractr) #fast & efficient raster extraction
library(rgee) #access & extract layers stored on Google Earth Engine
library(readxl) #read in lookup google sheet

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_NationalModels/NationalModels4.1"

#3. Get extraction methods lookup table----
meth <- read_excel(file.path(root, "NationalModels_V4.1_VariableList.xlsx"), sheet = "ExtractionLookup")

#A. DATA PREP####

#1. Load data----
load(file.path(root, "BirdData", "02_NM4.1_data_offsets.R"))

#2. Create sf object of just location*year for annual layers----
#transform to 5072
loc.yr <- visit %>% 
  dplyr::select(project, location, lat, lon, year) %>% 
  unique() %>% 
  st_as_sf(coords=c("lon", "lat"), crs=4326, remove=FALSE) %>% 
  st_transform(crs=5072)

#EXTRA. Test sample dataset----
set.seed(1234)
loc.n <- loc.yr %>% 
  sample_n(100)

#3. Buffer location objects----
loc.buff <- st_buffer(loc.n, 200)

#B. EXTRACT COVARIATES FROM GOOGLE DRIVE####

#1. Get list of layers to run----
meth.gd <- dplyr::filter(meth, Source=="Google Drive", Running==1)

#2. Plain dataframe for joining to output----
loc.df <- data.frame(loc.n) %>% 
  dplyr::select(-geometry)

#3. Set up loop----
for(i in 41:nrow(meth.gd)){
  
  #4. Determine if stacking----
  #increase i to skip the stacked rows of data
  if(meth.gd$StackCategory[i]==1){
    meth.gd.i <- dplyr::filter(meth.gd, Category==meth.gd$Category[i])
    i <- i + nrow(meth.gd.i) - 1
  } 
  if(meth.gd$StackCategory[i]==0){
    meth.gd.i <- meth.gd[i,]}
  
  #5. Determine if temporally static----
  if(meth.gd$TemporalResolution[i]=="static"){
    
    #6. Read in raster----
    rast.i <- rast(meth.gd.i$Link)
    
    #7. Extract - determine point or buffer extraction---
    if(meth.gd$Extraction[i]=="point"){
      loc.cov <- loc.n %>% 
        st_transform(crs(rast.i)) %>% 
        terra::extract(x=rast.i, ID=FALSE)
    }
    
    if(meth.gd$Extraction[i]=="radius"){
      loc.cov <- loc.buff %>% 
        st_transform(crs(rast.i)) %>% 
        exact_extract(x=rast.i, "mean")
    }
    
  }
  
  #8. Determine if temporally matching----
  if(meth.gd$TemporalResolution[i]=="match"){
    
    #9. Get list of individual files----
    files.i <- data.frame(Link=list.files(meth.gd.i$Link, full.names = TRUE),
                          file=list.files(meth.gd.i$Link)) %>% 
      mutate(year = as.numeric(str_sub(file, -8, -5))) %>% 
      arrange(year)
    
    #10. Match year of file to year of data----
    #http://adomingues.github.io/2015/09/24/finding-closest-element-to-a-number-in-a-list/
    dt = data.table::data.table(year=files.i$year, val = files.i$year)
    data.table::setattr(dt, "sorted", "year")
    data.table::setkey(dt, year)
    loc.n$year.rd <- dt[J(loc.n$year), roll = "nearest"]$val
    loc.buff$year.rd <- dt[J(loc.buff$year), roll = "nearest"]$val
    
    #11. Set up to loop through years----
    loc.cov <- data.frame()
    yrs <- unique(files.i$year)
    for(j in 1:length(yrs)){
      
      loc.n.j <- dplyr::filter(loc.n, year.rd==yrs[j])
      loc.buff.j <- dplyr::filter(loc.buff, year.rd==yrs[j])
      
      #12. Read in raster----
      rast.i <- rast(files.i$Link[j])
      
      #13. Extract - determine point or buffer extraction---
      if(meth.gd$Extraction[i]=="point"){
        loc.cov <- loc.n.j %>% 
          st_transform(crs(rast.i)) %>% 
          terra::extract(x=rast.i, ID=FALSE) %>% 
          rbind(loc.cov)
      }
      
      if(meth.gd$Extraction[i]=="radius"){
        loc.cov <- loc.buff.j %>% 
          st_transform(crs(rast.i)) %>% 
          exact_extract(x=rast.i, "mean") %>% 
          data.frame() %>% 
          rbind(loc.cov)
      }
      
    }
    
  }

  #14. Add to output----
  #fix column names
  if(!is.na(meth.gd$Priority[i])) meth.gd.i$Name <- paste0(meth.gd.i$Name, "-", meth.gd.i$Priority)
  nms <- c(colnames(loc.df)[1:(ncol(loc.df)-1)], meth.gd.i$Name)
  loc.df <- cbind(loc.df, loc.cov)
  names(loc.df) <- nms
  
  #15. Save----
  save(loc.df, file=file.path(root, "BirdData", "03_NM4.1_data_covariates_GD.R"))
  
  #16. Report status----
  print(paste0("Finished variable ", i, " of ", nrow(meth.gd), " - ", meth.gd$Name[i]))

}

#C. EXTRACT COVARIATES FROM GEE####

#1. Initialize GEE----
ee_Initialize()
ee_check()

#2. Create GCS bucket----
#This only needs to be done once!
# project_id <- ee_get_earthengine_path() %>%
#   list.files(., "\\.json$", full.names = TRUE) %>%
#   jsonlite::read_json() %>%
#   '$'(project_id) # Get the Project ID
# 
# googleCloudStorageR::gcs_create_bucket("national_models", projectId = project_id)

#3. Send polygons to GEE---
poly <- sf_as_ee(loc.buff)
point <- sf_as_ee(loc.n)

#4. Get list of layers to run----
meth.gd <- dplyr::filter(meth, Source=="Google Earth Engine", Running==1)

#4. Get layers----


#D. GET BCR####

#E. STRATIFY####

#F. 

#G. SAVE#####

save(visit, bird,  file=file.path(root, "01_NM4.1_data_clean.R"))