# ---
# title: National Models 4.1 - extract covariates
# author: Elly Knight
# created: December 22, 2022
# ---

#NOTES################################

#This script uses CONUS Albers equal-area conic (EPSG: 5072)

#This script extracts the list of model covariates agreed upon by the BAM National Model team. The full list of those covariates can be found at https://docs.google.com/spreadsheets/d/1XATuq8BOYC2KkbJObaturaf4NFD2ufvn/edit?usp=sharing&ouid=104837701987164094932&rtpof=true&sd=true

#This script extracts covariates for every combination of location & year. Efficiency could be improved by ~15% by taking only the unique locations for all rows in the method table with TemporalResolution=="static"

#Records that return NA values from the layers stored in google drive are locations that are outside the area of the raster (i.e., coastlines)

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(terra) #basic raster handling
library(sf) #basic shapefile handling
library(exactextractr) #fast & efficient raster extraction
library(rgee)
ee_Initialize(gcs=TRUE)

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_NationalModels/NationalModels4.1"

#3. Get extraction methods lookup table----
meth <- readxl::read_excel(file.path(root, "NationalModels_V4.1_VariableList.xlsx"), sheet = "ExtractionLookup")

#A. DATA PREP####

#1. Load data----
load(file.path(root, "Data", "02_NM4.1_data_offsets.R"))

#2. Create sf object of just location*year for annual layers----
loc.yr <- visit %>% 
  dplyr::select(project, location, lat, lon, year) %>% 
  unique() %>% 
  st_as_sf(coords=c("lon", "lat"), crs=4326, remove=FALSE) %>% 
  st_transform(crs=5072)

#EXTRA. Test sample dataset----
# set.seed(1234)
# loc.n <- loc.yr %>%
#   sample_n(1000)
loc.n <- loc.yr

#3. Buffer location objects----
#200 m radius
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
        exact_extract(x=rast.i, "mean", force_df=TRUE)
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
    loc.n.i <- loc.n
    loc.buff.i <- loc.buff
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
          exact_extract(x=rast.i, "mean", force_df=TRUE) %>% 
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

#17. Check output----
loc.check <- read.csv(file.path(root, "Data", "Covariates", "03_NM4.1_data_covariates_GD.csv"))
summary(loc.check)

#check missing climate normals - looks like all the coastal data
loc.na <- dplyr::filter(loc.check, is.na(CMD)) %>% 
  dplyr::select(project, year, lat, lon, CMD)
plot(loc.na$lon, loc.na$lat)

#check missing hli3cl - also looks coastal
loc.na <- dplyr::filter(loc.check, is.na(hli3cl))
plot(loc.na$lon, loc.na$lat)

#check missing peat - also looks coastal
loc.na <- dplyr::filter(loc.check, is.na(peat))
plot(loc.na$lon, loc.na$lat)

#check missing global HF - also looks coastal
loc.na <- dplyr::filter(loc.check, is.na(HF_5km.2))
plot(loc.na$lon, loc.na$lat)

#check missing greenup - looks like all of canada?
loc.na <- dplyr::filter(loc.check, is.na(greenup)) %>% 
  dplyr::select(project, year, lat, lon, greenup, dormancy)
plot(loc.na$lon, loc.na$lat)

write.csv(loc.na, file.path(root, "Data", "Covariates", "Greenup_NA.csv"), row.names=FALSE)

#TO DO: COME BACK TO MISSING GREENUP & DORMANCY VALUES####

#C. EXTRACT COVARIATES FROM SCANFI####

#1. Get list of covariates to run----
meth.scanfi <- dplyr::filter(meth, Source=="SCANFI", Running==1)

#2. Get & wrangle list of SCANFI layers----
files.scanfi <- data.frame(Link=list.files("G:/.shortcut-targets-by-id/11nj6IZyUe3EqrOEVEmDfrnkDCEQyinSi/SCANFI_share", full.names = TRUE,recursive=TRUE, pattern="*.tif"),
                      file=list.files("G:/.shortcut-targets-by-id/11nj6IZyUe3EqrOEVEmDfrnkDCEQyinSi/SCANFI_share", recursive=TRUE, pattern="*.tif")) %>% 
  dplyr::filter(str_sub(file, -3, -1)=="tif") %>% 
  mutate(file = str_sub(file, 13, 100)) %>% 
  separate(file, into=c("scanfi", "covtype", "variable", "S", "year", "v0", "filetype")) %>% 
  mutate(year = as.numeric(ifelse(filetype=="v0", v0, year)),
         variable = case_when(variable=="VegTypeClass" ~ "landcover",
                              variable=="prcC" ~ "conifer",
                              variable=="prcD" ~ "deciduous",
                              !is.na(variable) ~ tolower(variable))) %>% 
  dplyr::filter(scanfi=="SCANFI",
                variable %in% meth.scanfi$Name) %>% 
  dplyr::select(Link, variable, year) %>% 
  arrange(year, variable)
#Check everything is there
table(files.scanfi$variable, files.scanfi$year)

#3. Match year of data to year of SCANFI----
years.scanfi <- unique(files.scanfi$year)
dt = data.table::data.table(year=years.scanfi, val=years.scanfi)
data.table::setattr(dt, "sorted", "year")
data.table::setkey(dt, year)
loc.buff.i <- loc.buff
loc.buff.i$year.rd <- dt[J(loc.buff$year), roll = "nearest"]$val

#4. Set up to loop through years of SCANFI----
loc.scanfi <- data.frame()
for(i in 1:length(years.scanfi)){
  
  loc.buff.yr <- dplyr::filter(loc.buff.i, year.rd==years.scanfi[i])

  #6. Read in raster----
  files.i <- dplyr::filter(files.scanfi, year==years.scanfi[i])
  rast.i <- rast(files.i$Link)
  names(rast.i) <- files.i$variable
  
  #7. Extract----
  loc.scanfi <- data.frame(loc.buff.yr) %>% 
    dplyr::select(-geometry, -year.rd) %>% 
    cbind(loc.buff.yr %>% 
            st_transform(crs(rast.i)) %>% 
            exact_extract(x=rast.i, "mean", force_df=TRUE)) %>% 
    rbind(loc.scanfi)
  
  print(paste0("Finished year ", i, " of ", length(years.scanfi), ": ", years.scanfi[i]))
  
}

#8. Fix column names----
colnames(loc.scanfi) <- colnames(loc.df, unique(files.scanfi)) 

#9. Save----
write.csv(loc.scanfi, file=file.path(root, "Data", "Covariates", "03_NM4.1_data_covariates_SCANFI.csv"))

#D. EXTRACT COVARIATES FROM GEE####

#1. Create GCS bucket----
#This only needs to be done once!
# project_id <- ee_get_earthengine_path() %>%
#   list.files(., "\\.json$", full.names = TRUE) %>%
#   jsonlite::read_json() %>%
#   '$'(project_id) # Get the Project ID
# 
# googleCloudStorageR::gcs_create_bucket("national_models", projectId = project_id)

#2. Get list of static layers to run----
meth.gee <- dplyr::filter(meth, Source=="Google Earth Engine", Running==1, TemporalResolution=="static")

#3. Set up loop for static layer stacking----
for(i in 1:nrow(meth.gee)){

  #4. Get the image----
  if(meth.gee$GEEtype[i]=="image"){
    if(!is.na(meth.gee$GEEBand[i]))
    img.i <- ee$Image(meth.gee$Link[i])$select(meth.gee$GEEBand[i]) 
    else img.i <- ee$Image(meth.gee$Link[i]) 
  }
  if(meth.gee$GEEtype[i]=="imagecollection"){
    img.i <- ee$ImageCollection(meth.gee$Link[i])$select(meth.gee$GEEBand[i])$toBands()
  }
  
  #5. Stack----
  if(i==1){
    img.stack <- img.i
  }
  
  if(i > 1){
    img.stack <- img.stack$addBands(img.i)
  }
  
}

#6. Set up to do batches of 1000----
loc.buff$loop <- ceiling(row_number(loc.buff)/1000)

task.list <- list()
for(i in 1:max(loc.buff$loop)){
  
  loc.buff.i <- dplyr::filter(loc.buff, loop==i)
  
  #7. Send polygons to GEE---
  poly <- sf_as_ee(loc.buff.i)
  
  #8. Extract----
  img.red <- img.stack$reduceRegions(reducer=ee$Reducer$mean(),
                                     collection=poly)
  
  task.list[[i]] <- ee_table_to_gcs(collection=img.red,
                                 bucket="national_models",
                                 fileFormat = "CSV")
  task.list[[i]]$start()
  
  #10. Check is running
  if(i==1){
    ee_monitoring(task.list[[i]], max_attempts=1000)
  }

  print(paste0("Finished batch ", i, " of ", max(loc.buff$loop)))
  
}

#11. Save the task list----
save(task.list, file.path(root, "Data", "Covariates", "GEETasklist.Rdata"))

#11. Download from cloud----
for(i in 1:length(task.list)){
  ee_gcs_to_local(task = task.list[[i]], dsn=file.path(root, "Data", "Covariates", "GEE", paste0("03_NM4.1_data_covariates_GEE_static_", i, ".csv")))
}

#12. Collapse to one object----
files.gee <- list.files(file.path(root, "Data", "Covariates", "GEE"), full.names = TRUE)
dat.gee <- purrr::map(files.gee, read.csv) %>% 
  data.table::rbindlist() %>% 
  dplyr::select(c(colnames(loc.yr)[colnames(loc.yr)!="geometry"], meth.gee$GEEName))

#13. Fix column names----
meth.gee$Name <- ifelse(is.na(meth.gee$Priority), meth.gee$Name, paste0(meth.gee$Name, ".", meth.gee$Priority))
colnames(dat.gee) <- c(colnames(loc.yr)[colnames(loc.yr)!="geometry"], meth.gee$Name)

#14. Zerofill----
zerocols <- meth.gee %>% 
  dplyr::filter(Zerofill==1)
dat.gee <- dat.gee %>% 
  dplyr::select(zerocols$Name) %>% 
  mutate()




#D. GET BCR####

#E. STRATIFY####

#F. 

#G. SAVE#####

save(visit, bird,  file=file.path(root, "03_NM4.1_data_covariates.R"))
