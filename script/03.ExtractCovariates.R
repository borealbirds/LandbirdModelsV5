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

#CHECKLIST
#GOOGLE DRIVE - DONE EXCEPT GREENUP VARS, TRI
#SCANFI - NOT DONE - ISSUES WITH RASTER READING - TRY RUNNING IN CLUSTERS OF 1000 - NOPE
#GOOGLE EARTH STATIC - DONE
#GOOGLE EARTH TEMPORAL - NOT DONE
#WRANGLING - NOT DONE


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
#meth.gd <- dplyr::filter(meth, Source=="Google Drive", Running==1)
meth.gd <- dplyr::filter(meth, Name=="TRI")

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
                              variable=="LodepolePine" ~ "lodgepolepine",
                              !is.na(variable) ~ tolower(variable))) %>% 
  dplyr::filter(scanfi%in%c("SCANFI", "CaNFIR"),
                variable %in% meth.scanfi$Name) %>% 
  dplyr::select(Link, variable, year) %>% 
  arrange(year, variable)
#Check everything is there
table(files.scanfi$variable, files.scanfi$year)

#3. Remove points outside of SCANFI coverage----
#do this for SCANFI but not for others because scanfi is so much more time intensive (many high resolution layers)
rast.scanfi <- rast(files.scanfi$Link[1])

loc.scanfi.sa <- loc.n %>% 
  st_transform(crs(rast.scanfi)) %>% 
  terra::extract(x=rast.scanfi, ID=FALSE)
colnames(loc.scanfi.sa) <- "scanfi"

#4. Rebuffer----
loc.scanfi.buff <- cbind(loc.n, loc.scanfi.sa) %>% 
  dplyr::filter(!is.na(scanfi)) %>% 
  st_transform(5072) %>% 
  st_buffer(200) %>% 
  st_transform(crs(rast.scanfi))

#5. Match year of data to year of SCANFI----
#Use only buffer object because all extractions are radius method
years.scanfi <- unique(files.scanfi$year)
dt = data.table::data.table(year=years.scanfi, val=years.scanfi)
data.table::setattr(dt, "sorted", "year")
data.table::setkey(dt, year)
loc.scanfi.buff$year.rd <- dt[J(loc.scanfi.buff$year), roll = "nearest"]$val

#6. Set up to loop through years of SCANFI----
loc.scanfi <- data.frame()
loc.error <- data.frame()
for(i in 1:length(years.scanfi)){
  
  loc.buff.yr <- dplyr::filter(loc.scanfi.buff, year.rd==years.scanfi[i]) %>% 
    arrange(lat, lon) %>% 
    mutate(loop = ceiling(row_number()/1000))

  #7. Read in raster----
  files.i <- dplyr::filter(files.scanfi, year==years.scanfi[i])
  rast.i <- rast(files.i$Link)
  names(rast.i) <- files.i$variable
  
  #8. Set up loops----
  for(j in 1:max(loc.buff.yr$loop)){
    
    loc.i <- dplyr::filter(loc.buff.yr, loop==j)
    
    #9. Extract----
    loc.scanfi.i <- try(exact_extract(x=rast.i, y=loc.i, "mean", force_df=TRUE))
    
    #10. Save output if extraction works----
    if(class(loc.scanfi.i)=="data.frame"){
      loc.scanfi <- rbind(loc.scanfi,
                          cbind(data.frame(loc.i) %>% 
                                  dplyr::select(-geometry, -year.rd, -loop, -scanfi),
                                loc.scanfi.i))
    }
    
    #11. Log if extraction does note work----
    if(class(loc.scanfi.i)!="data.frame"){
      loc.error <- rbind(loc.error, loc.i)
    }
    
    #12. Save----
    write.csv(loc.scanfi, file=file.path(root, "Data", "Covariates", "03_NM4.1_data_covariates_SCANFI.csv"))
    
    #13. Report progress----
    print(paste0("Finished loop ", j, " of ", max(loc.buff.yr$loop), " for year ", i, " of ", length(years.scanfi), ": ", years.scanfi[i]))
    
    flush.console()
    
  }
  
}

write.csv(loc.error, file=file.path(root, "Data", "Covariates", "SCANFIerrors.csv"), row.names = FALSE)

# #10. Fix column names----
# colnames(loc.scanfi) <- colnames(unique(files.scanfi)) 
# 
# #11. Add the data outside of scanfi coverage back in
# 
# #12. Save again----
# write.csv(loc.scanfi, file=file.path(root, "Data", "Covariates", "03_NM4.1_data_covariates_SCANFI.csv"))

#D. EXTRACT COVARIATES FROM GEE - TEMPORALLY STATIC####

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
for(i in 580:max(loc.buff$loop)){
  
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
  
  #9. Monitor----
  ee_monitoring(task.list[[i]], max_attempts=1000)
  
  #10.Download to local----
  ee_gcs_to_local(task = task.list[[i]], dsn=file.path(root, "Data", "Covariates", "GEE", paste0("03_NM4.1_data_covariates_GEE_static_", i, ".csv")))

  print(paste0("Finished batch ", i, " of ", max(loc.buff$loop)))
  
}

#11. Collapse to one object----
files.gee <- list.files(file.path(root, "Data", "Covariates", "GEE"), full.names = TRUE)
loc.gee <- purrr::map(files.gee, read.csv) %>% 
  data.table::rbindlist() %>% 
  dplyr::select(c(colnames(loc.yr)[colnames(loc.yr)!="geometry"], meth.gee$GEEName))

#12. Fix column names----
meth.gee$Name2 <- ifelse(is.na(meth.gee$Priority), meth.gee$Name, paste0(meth.gee$Name, ".", meth.gee$Priority))
colnames(loc.gee) <- c(colnames(loc.yr)[colnames(loc.yr)!="geometry"], meth.gee$Name2)

#13. Zerofill----
zerocols <- meth.gee %>% 
  dplyr::filter(Zerofill==1)
zero.gee <- loc.gee %>% 
  dplyr::select(zerocols$Name2) %>% 
  replace_na(list(occurrence=0, recurrence=0, seasonality=0))
loc.gee2 <- loc.gee %>% 
  dplyr::select(-zerocols$Name2) %>% 
  cbind(zero.gee)

#14. Save----
write.csv(loc.gee2, file=file.path(root, "Data", "Covariates", "03_NM4.1_data_covariates_GEE-static.csv"))

#E. EXTRACT COVARIATES FROM GEE - TEMPORALLY MATCHED####

#1. Get list of static layers to run----
meth.gee <- dplyr::filter(meth, Source=="Google Earth Engine", Running==1, TemporalResolution=="match")

#2. Landcover categories----
lc <- data.frame(landcover = c(1:17),
                 lcclass = c("evergreen_needleleaf",
                             "evergreen_broadleaf",
                             "deciduous_needleleaf",
                             "deciduous_broadleaf",
                             "mixed",
                             "shrub_closed",
                             "shrub_open",
                             "savanna_woody",
                             "savanna_open",
                             "grassland",
                             "wetland",
                             "crop_dense",
                             "urban",
                             "crop_natural",
                             "snow", 
                             "barren",
                             "water"))


#3. Plain dataframe for joining to output----
loc.gee <- data.frame(loc.n) %>% 
  dplyr::select(-geometry)

#4. Set up to loop through the layers----
#not worth stacking because almost all layers have different temporal filtering settings
for(i in 1:nrow(meth.gee)){
  
  #5. Identify years of imagery----
  years.gee <- seq(meth.gee$GEEYearMin[i], meth.gee$GEEYearMax[i])
  
  #6. Match year of data to year of data----
  #Use only point object because all extractions are point method
  dt = data.table::data.table(year=years.gee, val=years.gee)
  data.table::setattr(dt, "sorted", "year")
  data.table::setkey(dt, year)
  loc.n.i <- loc.n
  loc.n.i$yearrd <- dt[J(loc.n$year), roll = "nearest"]$val
  
  #7. Set up to loop through years----
  loc.j <- list()
  for(j in 1:length(years.gee)){
    
    loc.n.yr <- dplyr::filter(loc.n.i, yearrd==years.gee[j]) %>% 
      mutate(loop = ceiling(row_number()/1000))
    
    if(nrow(loc.n.yr) > 0){
      
      #8. Set up to loop through sets----
      loc.k <- list()
      for(k in 1:max(loc.n.yr$loop)){
        
        loc.n.loop <- dplyr::filter(loc.n.yr, loop==k)
        
        #9. Send polygons to GEE---
        point <- sf_as_ee(loc.n.loop)
        
        #10. Set start & end date for image filtering---
        start.k <- paste0(years.gee[j]+meth.gee$YearMatch[i], "-", meth.gee$GEEMonthMin[i], "-01")
        
        if(meth.gee$GEEMonthMax[i] > meth.gee$GEEMonthMin[i]){
          end.k <- paste0(years.gee[j]+meth.gee$YearMatch[i], "-", meth.gee$GEEMonthMax[i], "-28")
        }
        if(meth.gee$GEEMonthMax[i] < meth.gee$GEEMonthMin[i]){
          end.k <- paste0(years.gee[j], "-", meth.gee$GEEMonthMax[i], "-28")
        }

        #11. Get the image----
        img.i <- ee$ImageCollection(meth.gee$Link[i])$filter(ee$Filter$date(start.k, end.k))$select(meth.gee$GEEBand[i])$mean()
        
        #12. Extract----
        loc.k[[k]] <- ee_extract(x=img.i, y=point)
        
        print(paste0("Finished batch ", k, " of ", max(loc.n.yr$loop)))
        
      }
      
      #13. Collapse loops for the year----
      loc.j[[j]] <- data.table::rbindlist(loc.k)
      colnames(loc.j[[j]]) <- c(colnames(loc.j[[j]])[1:ncol(loc.j[[j]])-1],meth.gee$Name[i])
      
    }
    
    print(paste0("Finished year ", j, " of ", length(years.gee)))
    
  }
  
  #14. Collapse data across years----
  loc.i <- data.table::rbindlist(loc.j)
  loc.gee <- cbind(loc.gee, loc.i %>% 
                     dplyr::select(meth.gee$Name[i]))
  
  #15. Join to landcover classes----
  if(meth.gee$Name[i]=="landcover"){
    loc.gee <- left_join(loc.gee, lc)
  }
  
  write.csv(loc.gee, file=file.path(root, "Data", "Covariates", "03_NM4.1_data_covariates_GEE-match.csv"))
  
  print(paste0("FINISHED LAYER ", i, " of ", nrow(meth.gee)))
  
}

#D. GET BCR####

#E. STRATIFY####

#F. 

#G. SAVE#####

save(visit, bird,  file=file.path(root, "03_NM4.1_data_covariates.R"))
