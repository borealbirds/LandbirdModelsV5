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
load(file.path(root, "Data", "02_NM4.1_data_offsets.R"))

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
loc.gd <- data.frame(loc.n) %>% 
  dplyr::select(-geometry)

#3. Set up loop----
#TO DO: FIX COUNT OF LOOPS CALCULATION####
#for(i in 1:nrow(meth.gd)){
for(i in 1:12){
    
  
  #4. Determine if stacking----
  #remove additional stacked layers from layer list
  if(meth.gd$StackCategory[i]==1){
    meth.gd.i <- dplyr::filter(meth.gd, Category==meth.gd$Category[i])
    meth.gd <- anti_join(meth.gd, meth.gd.i[2:nrow(meth.gd.i),])
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
    loc.n.i$year.rd <- dt[J(loc.n$year), roll = "nearest"]$val
    loc.buff.i$year.rd <- dt[J(loc.buff$year), roll = "nearest"]$val
    
    #11. Set up to loop through years----
    loc.cov <- data.frame()
    yrs <- unique(files.i$year)
    for(j in 1:length(yrs)){
      
      loc.n.j <- dplyr::filter(loc.n.i, year.rd==yrs[j])
      loc.buff.j <- dplyr::filter(loc.buff.i, year.rd==yrs[j])
      
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
  nms <- c(colnames(loc.gd)[1:ncol(loc.gd)], meth.gd.i$Name)
  loc.gd <- cbind(loc.gd, loc.cov)
  names(loc.gd) <- nms
  
  #15. Save----
  write.csv(loc.gd, file=file.path(root, "Data", "Covariates", "03_NM4.1_data_covariates_GD.csv"))
  
  #16. Report status----
  print(paste0("Finished variable ", i, " of ", nrow(meth.gd), " - ", meth.gd$Name[i]))

}

#C. EXTRACT COVARIATES FROM GEE####

#1. Initialize GEE----
ee_Initialize(gcs=TRUE)

#2. Create GCS bucket----
#This only needs to be done once!
# project_id <- ee_get_earthengine_path() %>%
#   list.files(., "\\.json$", full.names = TRUE) %>%
#   jsonlite::read_json() %>%
#   '$'(project_id) # Get the Project ID
# 
# googleCloudStorageR::gcs_create_bucket("national_models", projectId = project_id)

#3. Get list of static layers to run----
meth.gee <- dplyr::filter(meth, Source=="Google Earth Engine", Running==1, TemporalResolution=="static")

#4. Set up loop for static layer stacking----
for(i in 1:nrow(meth.gee)){
  
  #5. Get the image----
  if(meth.gee$GEEtype[i]=="image"){
    img.i <- ee$Image(meth.gee$Link[i])$select(meth.gee$Name[i])
  }
  if(meth.gee$GEEtype[i]=="imagecollection"){
    img.i <- ee$ImageCollection(meth.gee$Link[i])$select(meth.gee$Name[i])$toBands()
  }
  
  #6. Stack----
  if(i==1){
    img.stack <- img.i
  }
  
  if(i > 1){
    img.stack <- img.stack$addBands(img.i)
  }
  
}

#7. Set up to do batches of 1000----
loc.buff$loop <- ceiling(row_number(loc.buff)/1000


#8. Send polygons to GEE---
poly <- sf_as_ee(loc.buff.i)

#8. Extract----
img.red <- img.stack$reduceRegions(reducer=ee$Reducer$mean(),
                                   collection=poly)

task_vector <- ee_table_to_gcs(collection=img.red,
                               bucket="national_models",
                               fileFormat = "CSV")
task_vector$start()
ee_monitoring(task_vector, max_attempts=1000)

#9. Save to GD----
ee_gcs_to_local(task = task_vector, dsn=file.path(root, "Data", "Covariates", "03_NM4.1_data_covariates_GEE_static.csv"))

check <- read.csv(file.path(root, "Data", "Covariates", "03_NM4.1_data_covariates_GEE_static.csv"))






#D. GET BCR####

#E. STRATIFY####

#F. 

#G. SAVE#####

save(visit, bird,  file=file.path(root, "01_NM4.1_data_clean.R"))
