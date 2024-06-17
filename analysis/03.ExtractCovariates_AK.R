# ---
# title: National Models 5.0 - extract covariates
# author: Elly Knight
# created: January 17, 2024
# ---

#NOTES################################

#This script is a replicate of the "03.ExtractCovariates.R" script, but was added after that script was run to add an additional set of data for Alaska to the national models v5.

#Future versions of the National Models should consider more modularity from the initial data acquisition step to make it easier to add covariates part way through the process.

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(terra) #basic raster handling
library(sf) #basic shapefile handling
library(exactextractr) #fast & efficient raster extraction
library(rgee)
ee_Initialize(gcs=TRUE)

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_NationalModels5"

#3. Get extraction methods lookup table----
meth <- readxl::read_excel(file.path(root, "NationalModels_V5_VariableList.xlsx"), sheet = "ExtractionLookup")

#A. DATA PREP####

#1. Load data----
load(file.path(root, "Data", "02_NM5.0_data_offsets.R"))
rm(bird)
rm(offsets)

#2. Subset to just the new Alaska data----
visit.ak <- visit %>% 
  dplyr::filter(project=="ALMS2002-22")

#2. Create sf object of just location*year for annual layers----
loc.yr <- visit.ak %>% 
  dplyr::select(project, location, lat, lon, year) %>% 
  unique() %>% 
  mutate(id=paste(project, location, lat, lon, year, sep="_")) %>% 
  st_as_sf(coords=c("lon", "lat"), crs=4326, remove=FALSE) %>% 
  st_transform(crs=5072)

#EXTRA. Test sample dataset----
# set.seed(1234)
# loc.n <- loc.yr %>%
#   sample_n(100000)
loc.n <- loc.yr

#3. Buffer location objects----
#200 m radius for local extent
#2 km radius to landscape extent
loc.buff <- st_buffer(loc.n, 200)
loc.buff2 <- st_buffer(loc.n, 2000)

#B. EXTRACT COVARIATES FROM GOOGLE DRIVE####

#1. Get list of layers to run----
meth.gd <- dplyr::filter(meth, Source=="Google Drive", AK==1)

#2. Plain dataframe for joining to output----
loc.gd <- data.frame(loc.n) %>%
  dplyr::select(-geometry)

#3. Set up loop----
loop <- sort(unique(meth.gd$StackCategory))
for(i in 1:length(loop)){
  
  #4. Filter to stack category----
  meth.gd.i <- dplyr::filter(meth.gd, StackCategory==loop[i])
  
  #5. Determine if temporally static----
  if(meth.gd.i$TemporalResolution[1]=="static"){
    
    #6. Read in raster----
    rast.i <- rast(meth.gd.i$Link)
    names(rast.i) <- meth.gd.i$Label
    
    #7. Extract - determine point or buffer extraction and buffer extent---
    if(meth.gd.i$Extraction[1]=="point"){
      
      #Just extract for everything because point extraction is fast
      loc.cov <- loc.n %>% 
        st_transform(crs(rast.i)) %>% 
        terra::extract(x=rast.i, ID=FALSE)
      
      loc.cov$id <- loc.n$id
    }
    
    if(meth.gd.i$Extraction[1]=="radius"){
      
      #Extract point value, filter out NAs, reproject to UTM, buffer, reproject to raster CRS, extract value, rejoin to full df
      
      loc.cov.i <- loc.n %>% 
        st_transform(crs(rast.i)) %>% 
        terra::extract(x=rast.i, bind=TRUE) %>% 
        st_as_sf() %>% 
        rename(cov = names(rast.i)[1]) %>% 
        dplyr::filter(!is.na(cov)) %>% 
        dplyr::select(-cov) %>% 
        st_transform(crs=5072) %>% 
        st_buffer(meth.gd.i$RadiusExtent[1]) %>% 
        st_transform(crs(rast.i))
      
      loc.cov <- cbind(loc.cov.i,
                       exact_extract(x=rast.i,
                                     y=loc.cov.i,
                                     meth.gd.i$RadiusFunction[1],
                                     force_df=TRUE)) %>% 
        st_drop_geometry() %>% 
        right_join(loc.n %>% 
                     st_drop_geometry())
      
    }
    
  }
  
  #8. Determine if temporally matching----
  if(meth.gd.i$TemporalResolution[1]=="match"){
    
    #9. Get list of individual files----
    
    #For Landfire file structure
    if(loop[i] %in% c(17:24)){
      files.i <- data.frame(Link=list.files(meth.gd.i$Link, full.names = TRUE, recursive=TRUE, pattern="*.tif"),
                            file=list.files(meth.gd.i$Link, recursive=TRUE, pattern="*.tif")) %>% 
        separate(file, into=c("year", "region", "name", "tif"), remove=FALSE) %>% 
        mutate(name = case_when(name %in% c("105cbd", "130cbd", "140cbd", "CBD") ~ "biomass",
                                name %in% c("105cc", "130cc", "140cc", "CC") ~ "closure",
                                name %in% c("105ch", "130ch", "140ch", "CH") ~ "height"),
               year = as.numeric(year)) %>% 
        arrange(year) %>% 
        unique()
      
      #Just height for CV extraction
      if(loop[i] %in% c(21:24)){
        files.i <- files.i %>% 
          dplyr::filter(name=="height")
      }
      
    }
    
    #For NLCD file structure
    if(loop[i]==25){
      files.i <- data.frame(Link=list.files(meth.gd.i$Link, full.names = TRUE, pattern="*.img"),
                            file=list.files(meth.gd.i$Link, pattern="*.img")) %>% 
        mutate(year = as.numeric(str_sub(file, -32, -29))) %>% 
        arrange(year) %>% 
        unique() %>% 
        dplyr::filter(!is.na(year))
    }
    
    #Everything else
    if(!loop[i] %in% c(17:25)){
      files.i <- data.frame(Link=list.files(meth.gd.i$Link, full.names = TRUE, pattern="*tif", recursive=TRUE),
                            file=list.files(meth.gd.i$Link, pattern="*tif", recursive=TRUE)) %>% 
        mutate(year = as.numeric(str_sub(file, -8, -5))) %>% 
        arrange(year) %>% 
        unique() %>% 
        dplyr::filter(!is.na(year))
    }
    
    #10. Match year of file to year of data----
    #http://adomingues.github.io/2015/09/24/finding-closest-element-to-a-number-in-a-list/
    dt = data.table::data.table(year=unique(files.i$year), val = unique(files.i$year))
    data.table::setattr(dt, "sorted", "year")
    data.table::setkey(dt, year)
    loc.n.i <- loc.n
    loc.buff.i <- loc.buff
    loc.buff2.i <- loc.buff2
    loc.n.i$year.rd <- dt[J(loc.n$year), roll = "nearest"]$val
    loc.buff.i$year.rd <- dt[J(loc.buff$year), roll = "nearest"]$val
    loc.buff2.i$year.rd <- dt[J(loc.buff2$year), roll = "nearest"]$val
    
    #11. Set up to loop through years----
    loc.cov <- data.frame()
    yrs <- unique(files.i$year)
    for(j in 1:length(yrs)){
      
      loc.n.j <- dplyr::filter(loc.n.i, year.rd==yrs[j])
      loc.buff.j <- dplyr::filter(loc.buff.i, year.rd==yrs[j])
      loc.buff2.j <- dplyr::filter(loc.buff2.i, year.rd==yrs[j])
      
      #12. Read in raster----
      files.j <- dplyr::filter(files.i, year==yrs[j])
      rast.i <- rast(files.j$Link)
      names(rast.i) <- meth.gd.i$Label
      
      #13. Extract - determine point or buffer extraction and buffer extent---
      if(meth.gd.i$Extraction[1]=="point"){
        loc.cov <- loc.n.j %>% 
          st_transform(crs(rast.i)) %>% 
          terra::extract(x=rast.i, ID=FALSE) %>% 
          data.table::setnames(meth.gd.i$Label) %>% 
          cbind(loc.n.j %>% 
                  st_drop_geometry() %>% 
                  dplyr::select(id)) %>% 
          rbind(loc.cov)
      }
      
      if(meth.gd.i$Extraction[1]=="radius" & meth.gd.i$RadiusExtent[1]==200){
        loc.cov <- loc.buff.j %>% 
          st_transform(crs(rast.i)) %>% 
          exact_extract(x=rast.i, meth.gd.i$RadiusFunction, force_df=TRUE) %>% 
          data.table::setnames(meth.gd.i$Label) %>% 
          cbind(loc.n.j %>% 
                  st_drop_geometry() %>% 
                  dplyr::select(id)) %>% 
          rbind(loc.cov)
      }
      
      if(meth.gd.i$Extraction[1]=="radius" & meth.gd.i$RadiusExtent[1]==2000){
        loc.cov <- loc.buff2.j %>% 
          st_transform(crs(rast.i)) %>% 
          exact_extract(x=rast.i, meth.gd.i$RadiusFunction, force_df=TRUE) %>% 
          data.table::setnames(meth.gd.i$Label) %>% 
          cbind(loc.n.j %>% 
                  st_drop_geometry() %>% 
                  dplyr::select(id)) %>% 
          rbind(loc.cov) 
      }
      
      print(paste0("Finished year ", j, " of ", length(yrs)))
      
    }
    
  }
  
  #14. Fix column names----
  if(loop[i] %in% c(19, 20, 23, 24)){
    colnames(loc.cov) <- c(paste0(str_sub(meth.gd.i$Label, -100, -2), "_conus"), "id")
    nms <- c(colnames(loc.gd), paste(str_sub(meth.gd.i$Label, -100, -2), "_conus"))
  } else {nms <- c(colnames(loc.gd), meth.gd.i$Label) }
  
  #15. Add output to main file----
  loc.gd <- left_join(loc.gd, loc.cov)
  colnames(loc.gd) <- nms
  
  #16. Save----
  write.csv(loc.gd, file=file.path(root, "Data", "Covariates", "03_NM5.0_data_covariates_GD_AK.csv"), row.names=FALSE)
  
  #17. Report status----
  print(paste0("Finished stack category ", i, " of ", length(loop)))
  
}

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
meth.gee <- dplyr::filter(meth, Source=="Google Earth Engine", AK==1, TemporalResolution=="static")

#3. Make/get dataframe----
loc.gee.static <- data.frame()

#3. Make method loop----
loop <- meth.gee %>% 
  dplyr::select(RadiusFunction, RadiusExtent) %>% 
  unique()

for(h in 1:nrow(loop)){
  
  meth.gee.h <- meth.gee %>% 
    dplyr::filter(RadiusFunction==loop$RadiusFunction[h],
                  RadiusExtent==loop$RadiusExtent[h])
  
  #3. Set up loop for static layer stacking----
  for(i in 1:nrow(meth.gee.h)){
    
    #4. Get the image----
    if(meth.gee.h$GEEtype[i]=="image"){
      if(!is.na(meth.gee.h$GEEBand[i]))
        img.i <- ee$Image(meth.gee.h$Link[i])$select(meth.gee.h$GEEBand[i]) 
      else img.i <- ee$Image(meth.gee.h$Link[i]) 
    }
    if(meth.gee.h$GEEtype[i]=="imagecollection"){
      img.i <- ee$ImageCollection(meth.gee.h$Link[i])$select(meth.gee.h$GEEBand[i])$toBands()
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
  if(loop$RadiusExtent[h]==200){
    loc.gee <- loc.buff %>% 
      mutate(loop = ceiling(row_number()/1000))
  }
  
  if(loop$RadiusExtent[h]==2000){
    loc.gee <- loc.buff2 %>% 
      mutate(loop = ceiling(row_number()/1000))
  }
  
  
  task.list <- list()
  for(i in 1:max(loc.gee$loop)){
    
    set.seed(1234)
    loc.gee.i <- dplyr::filter(loc.gee, loop==i)
    
    #8. Send polygons to GEE---
    poly <- sf_as_ee(loc.gee.i)
    
    #9. Extract----
    if(loop$RadiusFunction[h]=="mean"){
      img.red <- img.stack$reduceRegions(reducer=ee$Reducer$mean(),
                                         collection=poly,
                                         scale=meth.gee$GEEScale[h])
    }
    
    if(loop$RadiusFunction[h]=="cv"){
      img.red <- img.stack$reduceRegions(reducer=ee$Reducer$stdDev(),
                                         collection=poly,
                                         scale=meth.gee$GEEScale[h])
    }
    
    task.list[[i]] <- ee_table_to_gcs(collection=img.red,
                                      bucket="national_models",
                                      fileFormat = "CSV")
    task.list[[i]]$start()
    
    #10. Monitor----
    ee_monitoring(task.list[[i]], max_attempts=1000)
    
    #11.Download to local----
    ee_gcs_to_local(task = task.list[[i]], dsn=file.path(root, "Data", "Covariates", "GEE_AK", "Static", paste0("03_NM5.0_data_covariates_GEE_static_", loop$RadiusFunction[h], "_", loop$RadiusExtent[h], "_", i, ".csv")))
    
    print(paste0("Finished batch ", i, " of ", max(loc.gee$loop)))
    
  }
  
  print(paste0("Finished loop ", h, " of ", nrow(loop)))
  
}

#12. Read in output files----
files.gee <- data.frame(path = list.files(file.path(root, "Data", "Covariates", "GEE_AK", "Static"), full.names = TRUE)) %>%
  separate(path, into=c("j1", "j2", "j3", "j4", "j5", "j6", "j7", "j8", "method", "extent", "n"), sep="_", remove=FALSE)

files.gee.cv.200 <- dplyr::filter(files.gee, method=="cv", extent=="200")
loc.gee.cv.200 <- purrr::map(files.gee.cv.200$path, read.csv) %>% 
  data.table::rbindlist()

files.gee.cv.2000 <- dplyr::filter(files.gee, method=="cv", extent=="2000")
loc.gee.cv.2000 <- purrr::map(files.gee.cv.2000$path, read.csv) %>% 
  data.table::rbindlist()

files.gee.mean.200 <- dplyr::filter(files.gee, method=="mean", extent=="200")
loc.gee.mean.200 <- purrr::map(files.gee.mean.200$path, read.csv) %>% 
  data.table::rbindlist()

files.gee.mean.2000 <- dplyr::filter(files.gee, method=="mean", extent=="2000")
loc.gee.mean.2000 <- purrr::map(files.gee.mean.2000$path, read.csv) %>% 
  data.table::rbindlist()

#13. Fix column names----
lab.cv.200 <- meth.gee %>% 
  dplyr::filter(RadiusFunction=="cv",
                RadiusExtent=="200")

lab.cv.2000 <- meth.gee %>% 
  dplyr::filter(RadiusFunction=="cv",
                RadiusExtent=="2000")

lab.mean.200 <- meth.gee %>% 
  dplyr::filter(RadiusFunction=="mean",
                RadiusExtent=="200")

lab.mean.2000 <- meth.gee %>% 
  dplyr::filter(RadiusFunction=="mean",
                RadiusExtent=="2000")

#14. Fix column names, put together, calculate----
loc.gee <- loc.gee.cv.200 %>% 
  data.table::setnames(c(colnames(loc.gee.cv.200)[1:7],
                         lab.cv.200$Label,
                         colnames(loc.gee.cv.200)[8:9])) %>% 
  dplyr::select(c(id, lab.cv.200$Label)) %>% 
  full_join(loc.gee.cv.2000 %>% 
              data.table::setnames(c(colnames(loc.gee.cv.2000)[1:7],
                                     lab.cv.2000$Label,
                                     colnames(loc.gee.cv.2000)[8:9])) %>% 
              dplyr::select(c(id, lab.cv.2000$Label))) %>% 
  full_join(loc.gee.mean.200 %>% 
              data.table::setnames(c(colnames(loc.gee.mean.200)[1],
                                     lab.mean.200$Label[1:2],
                                     colnames(loc.gee.mean.200)[4:8],
                                     lab.mean.200$Label[3],
                                     colnames(loc.gee.mean.200)[10],
                                     lab.mean.200$Label[4:5],
                                     colnames(loc.gee.mean.200)[13:14])) %>% 
              dplyr::select(c(id, lab.mean.200$Label))) %>% 
  full_join(loc.gee.mean.2000 %>% 
              data.table::setnames(c(colnames(loc.gee.mean.2000)[1],
                                     lab.mean.2000$Label[1:2],
                                     colnames(loc.gee.mean.2000)[4:8],
                                     lab.mean.2000$Label[3],
                                     colnames(loc.gee.mean.2000)[10:12])) %>% 
              dplyr::select(c(id, lab.mean.2000$Label))) %>% 
  mutate(ETHheightcv_1km = ETHheightcv_1km/ETHheight_1km,
         ETHheightcv_5x5 = ETHheightcv_5x5/ETHheight_5x5) %>% 
  mutate(ETHheightcv_1km = ifelse(is.infinite(ETHheightcv_1km), 0, ETHheightcv_1km),
         ETHheightcv_5x5 = ifelse(is.infinite(ETHheightcv_5x5), 0, ETHheightcv_5x5))

#15. Zerofill----
zerocols <- meth.gee %>% 
  dplyr::filter(Zerofill==1)
zero.gee <- loc.gee %>% 
  dplyr::select(zerocols$Label) %>% 
  replace_na(list(occurrence=0, recurrence=0, seasonality=0, occurrence_ls=0))
loc.gee.static <- loc.gee %>% 
  dplyr::select(-zerocols$Label) %>% 
  cbind(zero.gee)

#16. Save----
write.csv(loc.gee.static, file=file.path(root, "Data", "Covariates", "03_NM5.0_data_covariates_GEE-static_AK.csv"), row.names = FALSE)

#E. EXTRACT COVARIATES FROM GEE - TEMPORALLY MATCHED####

#1. Get list of static layers to run----
meth.gee <- dplyr::filter(meth, Source=="Google Earth Engine", TemporalResolution=="match", AK==1)

#2. Plain dataframe for joining to output----
loc.gee <- data.frame(loc.n) %>%
  dplyr::select(-geometry)

#3. Set up to loop through the layers----
#not worth stacking because almost all layers have different temporal filtering settings and scales
for(i in 1:nrow(meth.gee)){
  
  #4. Identify years of imagery----
  years.gee <- seq(meth.gee$GEEYearMin[i], meth.gee$GEEYearMax[i])
  
  #5. Match year of data to year of data----
  #Use only point object because all extractions are point method
  dt = data.table::data.table(year=years.gee, val=years.gee)
  data.table::setattr(dt, "sorted", "year")
  data.table::setkey(dt, year)
  loc.n.i <- loc.n
  loc.buff2.i <- loc.buff2
  loc.n.i$yearrd <- dt[J(loc.n$year), roll = "nearest"]$val
  loc.buff2.i$yearrd <- dt[J(loc.buff2$year), roll = "nearest"]$val
  
  #6 Set up to loop through years----
  loc.j <- list()
  for(j in 1:length(years.gee)){
    
    loc.n.yr <- dplyr::filter(loc.n.i, yearrd==years.gee[j]) %>% 
      mutate(loop = ceiling(row_number()/1000))
    loc.buff2.yr <- dplyr::filter(loc.buff2.i, yearrd==years.gee[j]) %>% 
      mutate(loop = ceiling(row_number()/1000))
    
    if(nrow(loc.n.yr) > 0){
      
      #7. Set up to loop through sets----
      loc.k <- list()
      task.list <- list()
      for(k in 1:max(loc.n.yr$loop)){
        
        loc.n.loop <- dplyr::filter(loc.n.yr, loop==k)
        loc.buff2.loop <- dplyr::filter(loc.buff2.yr, loop==k)
        
        #8. Send polygons to GEE---
        point <- sf_as_ee(loc.n.loop)
        poly <- sf_as_ee(loc.buff2.loop)
        
        #9. Set start & end date for image filtering---
        start.k <- paste0(years.gee[j]+meth.gee$YearMatch[i], "-", meth.gee$GEEMonthMin[i], "-01")
        
        if(meth.gee$GEEMonthMax[i] > meth.gee$GEEMonthMin[i]){
          end.k <- paste0(years.gee[j]+meth.gee$YearMatch[i], "-", meth.gee$GEEMonthMax[i], "-28")
        }
        if(meth.gee$GEEMonthMax[i] < meth.gee$GEEMonthMin[i]){
          end.k <- paste0(years.gee[j], "-", meth.gee$GEEMonthMax[i], "-28")
        }
        
        #10. Get the image----
        img.i <- ee$ImageCollection(meth.gee$Link[i])$filter(ee$Filter$date(start.k, end.k))$select(meth.gee$GEEBand[i])$mean()
        
        #11. Extract----
        if(meth.gee$Extraction[i]=="point"){
          loc.k[[k]] <- ee_extract(x=img.i, y=point, scale=meth$GEEScale[h])
        }
        
        if(meth.gee$Extraction[i]=="radius"){
          
          if(meth.gee$RadiusFunction[i]=="mean"){
            img.red <- img.i$reduceRegions(reducer=ee$Reducer$mean(),
                                           collection=poly,
                                           scale=meth.gee$GEEScale[h])
          }
          
          if(meth.gee$RadiusFunction[i]=="mode"){
            img.red <- img.i$reduceRegions(reducer=ee$Reducer$mode(),
                                           collection=poly, 
                                           scale=meth.gee$GEEScale[h])
          }
          
          task.list[[k]] <- ee_table_to_gcs(collection=img.red,
                                            bucket="national_models",
                                            fileFormat = "CSV")
          task.list[[k]]$start()
          
          ee_monitoring(task.list[[k]], max_attempts=1000)
          
          ee_gcs_to_local(task = task.list[[k]], dsn=file.path(root, "Data", "Covariates", "GEE_AK", "Match", paste0("03_NM5.0_data_covariates_GEE_match_", meth.gee$Label[i], "_", years.gee[j], "_", k, ".csv")))
          
          loc.k[[k]] <- read.csv(file.path(root, "Data", "Covariates", "GEE_AK", "Match", paste0("03_NM5.0_data_covariates_GEE_match_", meth.gee$Label[i], "_", years.gee[j], "_", k, ".csv")))
          
        }
        
        print(paste0("Finished batch ", k, " of ", max(loc.n.yr$loop)))
        
      }
      
      if(meth.gee$Extraction[i]=="point"){
        
        #12. Collapse loops for the year----
        loc.j[[j]] <- data.table::rbindlist(loc.k, fill=TRUE) 
        
        #13. Fix column names----
        colnames(loc.j[[j]]) <- c(colnames(loc.j[[j]])[1:ncol(loc.j[[j]])-1], meth.gee$Label[i])
      }
      
      if(meth.gee$Extraction[i]=="radius"){
        
        #12. Collapse loops for the year----
        loc.j[[j]] <- data.table::rbindlist(loc.k, fill=TRUE) %>% 
          dplyr::select(id, mode)
        
        #13. Fix column names----
        colnames(loc.j[[j]]) <- c("id", meth.gee$Label[i])
      }

      
    }
    
    print(paste0("Finished year ", j, " of ", length(years.gee)))
    
  }
  
  #14. Collapse data across years----
  loc.i <- data.table::rbindlist(loc.j, fill=TRUE)
  loc.gee <- cbind(loc.gee, loc.i %>% 
                     dplyr::select(meth.gee$Label[i]))
  
  write.csv(loc.gee, file=file.path(root, "Data", "Covariates", "03_NM5.0_data_covariates_GEE-match_AK.csv"), row.names=FALSE)
  
  print(paste0("FINISHED LAYER ", i, " of ", nrow(meth.gee)))
  
}

#16. Save again----
write.csv(loc.gee, file=file.path(root, "Data", "Covariates", "03_NM5.0_data_covariates_GEE-match_AK.csv"), row.names=FALSE)

#E. ASSEMBLE####

#1. Load data----
load(file.path(root, "Data", "02_NM5.0_data_offsets.R"))

#2. Load previously extracted covariates and bind with new ones----
loc.gee.match <- read.csv(file.path(root, "Data", "Covariates", "03_NM5.0_data_covariates_GEE-match.csv")) %>% 
  bind_rows(read.csv(file.path(root, "Data", "Covariates", "03_NM5.0_data_covariates_GEE-match_AK.csv")) %>% 
              dplyr::select(-project, -location, -lat, -lon, -year))

loc.gee.static <- read.csv(file.path(root, "Data", "Covariates", "03_NM5.0_data_covariates_GEE-static.csv")) %>% 
  bind_rows(read.csv(file.path(root, "Data", "Covariates", "03_NM5.0_data_covariates_GEE-static_AK.csv")))

loc.gd <- read.csv(file=file.path(root, "Data", "Covariates", "03_NM5.0_data_covariates_GD.csv")) %>% 
  bind_rows(read.csv(file.path(root, "Data", "Covariates", "03_NM5.0_data_covariates_GD_AK.csv")))

loc.scanfi <- read.csv(file.path(root, "Data", "Covariates", "03_NM5.0_data_covariates_SCANFI.csv"))

#3. Add id to visit and join together and wrangle----
meth.use <- meth %>% 
  dplyr::filter(Complete==1)

#Remove lat lon fields due to rounding errors that cause mismatches
#format landcover classes as factors
#zero out heights < 0.1 and NA the height cv values for those
#remove landfire values < 0, not sure what's going on there

visit.covs <- visit %>% 
  rename(row = id) %>% 
  mutate(id=paste(project, location, lat, lon, year, sep="_")) %>% 
  left_join(loc.gd %>% 
              dplyr::select(-lat, -lon)) %>% 
  left_join(loc.scanfi %>% 
              dplyr::select(-lat, -lon)) %>% 
  left_join(loc.gee.static) %>% 
  left_join(loc.gee.match) %>% 
  dplyr::select(-id) %>% 
  rename(id = row) %>% 
  dplyr::select(all_of(colnames(visit)), all_of(meth.use$Label)) %>% 
  mutate(MODISLCC_1km = factor(MODISLCC_1km),
         MODISLCC_5x5 = factor(MODISLCC_5x5),
         SCANFI_1km = factor(SCANFI_1km),
         NLCD_1km = factor(NLCD_1km),
         ABoVE_1km = factor(ABoVE_1km),
         VLCE_1km = factor(VLCE_1km)) %>% 
  mutate(LFheigth_1km = ifelse(LFheigth_1km < 0, NA, LFheigth_1km),
         LFheigth_5x5 = ifelse(LFheigth_5x5 < 0, NA, LFheigth_5x5),
         LFheigthcv_1km = ifelse(LFheigthcv_1km < 0, NA, LFheigthcv_1km),
         LFheigthcv_5x5 = ifelse(LFheigthcv_5x5 < 0, NA, LFheigthcv_5x5),
         LFbiomass_1km = ifelse(LFbiomass_1km < 0, NA, LFbiomass_1km),
         LFbiomass_5x5 = ifelse(LFbiomass_5x5 < 0, NA, LFbiomass_5x5),
         LFcrownclosure_1km = ifelse(LFcrownclosure_1km < 0, NA, LFcrownclosure_1km),
         LFcrownclosure_5x5 = ifelse(LFcrownclosure_5x5 < 0, NA, LFcrownclosure_5x5)) %>% 
  mutate(SCANFIheight_1km = ifelse(SCANFIheight_1km < 0.1, 0, SCANFIheight_1km),
         SCANFIheight_5x5 = ifelse(SCANFIheight_5x5 < 0.1, 0, SCANFIheight_5x5),
         SCANFIheightcv_1km = ifelse(SCANFIheight_1km < 0.1, NA, SCANFIheightcv_1km),
         SCANFIheightcv_5x5 = ifelse(SCANFIheight_1km < 0.1, NA, SCANFIheightcv_5x5),
         ETHheight_1km = ifelse(ETHheight_1km < 0.1, 0, ETHheight_1km),
         ETHheight_5x5 = ifelse(ETHheight_5x5 < 0.1, 0, ETHheight_5x5),
         ETHheightcv_1km = ifelse(ETHheight_1km < 0.1, NA, ETHheightcv_1km),
         ETHheightcv_5x5 = ifelse(ETHheight_1km < 0.1, NA, ETHheightcv_5x5),
         LFheigth_1km = ifelse(LFheigth_1km < 0.1, 0, LFheigth_1km),
         LFheigth_5x5 = ifelse(LFheigth_5x5 < 0.1, 0, LFheigth_5x5),
         LFheigthcv_1km = ifelse(LFheigth_1km < 0.1, NA, LFheigthcv_1km),
         LFheigthcv_5x5 = ifelse(LFheigth_1km < 0.1, NA, LFheigthcv_5x5))

#G. SAVE#####
visit <- visit.covs

save(visit, bird, offsets, file=file.path(root, "Data", "03_NM5.0_data_covariates.R"))
