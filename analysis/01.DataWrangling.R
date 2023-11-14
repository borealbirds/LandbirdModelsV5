# ---
# title: National Models 5.0 - get data
# author: Elly Knight, Melina Houle, Anna Drake
# created: Novemeber 17, 2022
# ---

#NOTES################################

#The "projectInstructions.csv" file is a list of all projects currently in WildTrax should not be used in national (instructions=="DO NOT USE"). This file should be updated for each iteration of in collaboration with Erin Bayne. Future versions of this spreadsheet can hopefully be derived by a combination of organization and a google form poll for consent from other organizations. Note this category also includes all ABMI projects with inaccurate (i.e., buffered) coordinates.

#The "projectInstructions.csv" file also contains information on which ARU projects are processed for a single species or taxa (instructions=="DO NOT USE") and therefore those visits should only be used for models for the appropriate taxa. This file should be updated for each iteration of the national models in collaboration with Erin Bayne. These projects are currently not included in the models.

#There are a handful of projects that are not downloading properly via wildRtrax. An issue is open on this. These projects are listed in the error.log object. These files should be downloaded manually.

#raw eBird data is downloaded from the eBird interface at https://ebird.org/data/download/ebd prior to wrangling with the "auk" package and will require a request for access. Use the custom download tool to download only the datasets for Canada and the US instead of the global dataset. Note you will also need the global sampling file to use the auk package for zero filling.

#raw eBird data omits Great Grey Owl & Northern Hawk Owl as sensitive species (https://support.ebird.org/en/support/solutions/articles/48000803210?b_id=1928&_gl=1*xq054u*_ga*ODczMTUyMjcuMTY2OTE0MDI4Ng..*_ga_QR4NVXZ8BM*MTY2OTE0MDI4NS4xLjEuMTY2OTE0MDM3OC4zNS4wLjA.&_ga=2.147122167.150058226.1669140286-87315227.1669140286) and should not be used for modelling these two species.

#wrangling eBird data with the auk package requires installation of AWK on windows computers. Please see #https://cornelllabofornithology.github.io/auk/articles/auk.html.

#eBird data has not been zerofilled because there was no species filtering done and we are assuming that all stationary counts have at least 1 bird observed.

#The column "sensor" currently only differentiates between ARU & human point count data types. Future versions should consider differentiating between SM2 ARU data types and other ARU data types due to differences in the perceptibility of these two approaches, either via QPAD or a correction factor.

#The replace TMTTs script will be replaced by a wildRtrax function in the near future.

#The filter by temporal covariates and study area components of the code should be moved from script 4 to this script for the next iteration for efficiency's sake.

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(wildRtrax) #to download data from wildtrax
library(data.table) #for binding lists into dataframes
library(lubridate) #date wrangling
library(auk) #eBird wrangling
library(purrr) #Mapping functions to read in multiple files
library(sf) #read in & handle study area
library(terra) #raster management
library(dggridR) #to make grid for checking for duplicates
library(QPAD) #to load species list for bird data

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_NationalModels/NationalModels5.0/Data"

#A. DOWNLOAD DATA FROM WILDTRAX#######################

#1. Login to WildTrax----
config <- "script/login.R"
source(config)

#1. Get list of projects from WildTrax----
wt_auth()

#sensor = PC gives all ARU and point count projects
project.list <- wt_get_download_summary(sensor_id = 'PC')

#2. Convert to a plain dataframe----
projects <- data.frame(project = as.character(project.list$project),
                       project_id = as.numeric(project.list$project_id),
                       sensorId = as.character(project.list$sensorId),
                       tasks = as.numeric(project.list$tasks),
                       status = as.character(project.list$status))

#3. Loop through projects to download data----
dat.list <- list()
error.log <- data.frame()
for(i in 1:nrow(projects)){
  
  #Do each sensor type separately because the reports have different columns and we need different things for each sensor type
  if(projects$sensorId[i]=="ARU"){
    
    #Get summary report
    report.try <- try(wt_download_report(project_id = projects$project_id[i], sensor_id = projects$sensorId[i], weather_cols = F, report = "summary"))
    
    #Get task report for ARU model
    task.try <- try(wt_download_report(project_id = projects$project_id[i], sensor_id = projects$sensorId[i], weather_cols = F, report = "task"))
    
    if(class(report.try)=="data.frame"){
      dat.try <- report.try %>% 
        left_join(task.try  %>% 
                    dplyr::select("organization", "project_name", "location", "recording_date", "longitude", "latitude", "method", "status", "observer", "observer_id", "equipment_used", "buffer"))
    }
  }
  
  if(projects$sensorId[i]=="PC"){
    
    dat.try <- try(wt_download_report(project_id = projects$project_id[i], sensor_id = projects$sensorId[i], weather_cols = F, report="report"))
    
  }
  
  if(class(dat.try)=="data.frame"){
    dat.list[[i]] <- dat.try
  }
  
  #Log projects that error
  if(class(dat.try)!="data.frame"){
    error.log <- rbind(error.log, 
                       projects[i,])
    
  }
  
  print(paste0("Finished dataset ", projects$project[i], " : ", i, " of ", nrow(projects), " projects"))
  
}

#4. Go download error log projects from wildtrax.ca----

#5. Read in error projects----
error.files <- list.files(file.path(root, "errorFiles"), full.names = TRUE)

raw.error <- data.frame()
for(i in 1:length(error.files)){
  raw.error <- read.csv(error.files[i]) %>% 
    rbind(raw.error)
}

#6. Collapse list----
#standardize column names between sensor types

all.wt <- rbindlist(dat.list, fill=TRUE)  %>% 
  rbind(raw.error, fill=TRUE) %>% 
  mutate(project = ifelse(is.na(project), project_name, project),
         speciesCode = ifelse(is.na(speciesCode), species_code, speciesCode),
         date = ifelse(is.na(date), recording_date, date),
         buffer = ifelse(is.na(buffer), bufferRadius.m., buffer)) %>%
  dplyr::select(-project_name, -recording_date, -species_code, -bufferRadius.m.) %>% 
  left_join(projects %>% 
              dplyr::rename(project_status = status))

#7. Filter out projects that shouldn't be used----
#nothing in BU training & all "DO NOT USE" projects in projectInventory file
#filter out 'NONE' method ARU projects later after this field is parsed out

instructions <- read.csv(file.path(root, "projectInventory", "projectInstructions.csv"))

raw.wt <- all.wt %>% 
  anti_join(instructions) %>%
  dplyr::filter(organization!="BU-TRAINING")

#9. Save date stamped data & project list----
save(raw.wt, projects, error.log, file=paste0(root, "/wildtrax_raw_", Sys.Date(), ".Rdata"))

#B. GET EBIRD DATA##########################

#Note: loop currently has to be run by hand because auk_set_ebd_path() requires restarting the session after running the command

#1. Get list of ebd objects to process----
ebd.files <- list.files(file.path(root, "ebd_raw"), pattern="ebd_*")

#2. Set up loop----
for(i in 1:length(ebd.files)){
  
  #3. Set ebd path----
  auk_set_ebd_path(file.path(root, "ebd_raw", ebd.files[i]), overwrite=TRUE)
  
  #4. Define filters----
  filters <- auk_ebd(paste0(ebd.files[i], ".txt")) %>% 
    auk_protocol("Stationary") %>% 
    auk_duration(c(0, 10)) %>% 
    auk_complete()
  
  #5. Filter data----
  #select columns to keep
  filtered <- auk_filter(filters, file=file.path(root, "ebd_filtered", paste0(ebd.files[i], ".txt")), overwrite=TRUE,
                         keep = c("group identifier", "sampling_event_identifier", "scientific name", "common_name", "observation_count", "latitude", "longitude", "locality_type", "observation_date", "time_observations_started", "observer_id", "duration_minutes"))
  
}

#5. Check list of processed files----
ebd.files.done <- list.files(file.path(root, "ebd_filtered"), pattern="ebd_*")

#C. HARMONIZE###############################

#1. Set desired columns----
colnms <- c("source", "organization", "project", "sensor", "tagMethod", "equipment", "location", "buffer", "lat", "lon", "year", "date", "observer", "duration", "distance", "species", "abundance", "isSeen", "isHeard")

#2. Wrangle wildtrax data-----

load(file.path(root, "wildtrax_raw_2023-01-20.Rdata"))

#2a. A bit of prep----
dat.wt <- raw.wt %>% 
  dplyr::select(-observer) %>% 
  rename(lat = latitude, lon = longitude, species = speciesCode, equipment = equipment_used, observer=observer_id) %>% 
  full_join(projects %>% 
              rename(sensor = sensorId) %>% 
              dplyr::select(project_id, project, sensor)) %>% 
  mutate(source = "WildTrax", 
         date = ymd_hms(date),
         year = year(date)) %>% 
  separate(buffer, into=c("buffer"), sep=" ", remove=TRUE) %>% 
  mutate(buffer = as.numeric(ifelse(!buffer %in% c("50", "10000"), str_sub(buffer, -100, -2), buffer)),
         buffer = ifelse(is.na(buffer), 0, buffer))

#2b. Get the point count data----
#wrangle distance and duration maximums
#remove counts with unknown duration and distance
pc.wt.meth <- dat.wt %>% 
  dplyr::filter(sensor=="PC",
                durationMethod!="UNKNOWN",
                distanceMethod!="UNKNOWN") %>% 
  dplyr::select(durationMethod, distanceMethod) %>% 
  unique() %>% 
  rowwise() %>% 
  mutate(durationMethod2 = ifelse(str_sub(durationMethod, -1, -1)=="+", str_sub(durationMethod, -100, -2), durationMethod),
         chardur = str_locate_all(durationMethod2, "-"),
         chardurmax = max(chardur),
         duration = as.numeric(str_sub(durationMethod2, chardurmax+1, -4)),
         chardis = str_locate_all(distanceMethod, "-"),
         chardismax = max(chardis),
         distance1 = str_sub(distanceMethod, chardismax+1, -2),
         distance = ifelse(distance1 %in% c("AR", "IN"), Inf, as.numeric(distance1))) %>% 
  dplyr::select(distanceMethod, durationMethod, distance, duration)

pc.wt <- dat.wt %>% 
  dplyr::filter(sensor=="PC") %>% 
  mutate(tagMethod = ifelse(distanceMethod=="0m-INF-ARU", "1SPT", "PC"),
         equipment = "Human") %>% 
  left_join(pc.wt.meth) %>% 
  dplyr::select(all_of(colnms)) %>% 
  data.frame()

#2c. Get the aru data----
#filter ARU data to first detection of each individual
#wrangle duration
#remove 'NONE' method
#wrangle equipment type
aru.wt.equip <- dat.wt %>% 
  dplyr::filter(sensor=="ARU",
                !is.na(equipment)) %>% 
  dplyr::select(equipment) %>% 
  unique() %>% 
  rowwise() %>% 
  mutate(equipmentend = str_locate(equipment, "\\("),
         equipmentlong = str_sub(equipment, 1, equipmentend[,1]-1),
         equipmentshort = case_when(equipmentlong=="?+" ~ "unknown",
                                    str_detect(equipmentlong, "SM2")==TRUE ~ "SM2",
                                    str_detect(equipmentlong, "SM3")==TRUE ~ "SM3",
                                    str_detect(equipmentlong, "SM4")==TRUE ~ "SM4",
                                    str_detect(equipmentlong, "mini")==TRUE ~ "mini",
                                    str_detect(equipmentlong, "Mini")==TRUE ~ "mini",
                                    is.na(equipmentlong) ~ "unknown"))

aru.wt <- dat.wt %>% 
  dplyr::filter(sensor=="ARU") %>% 
  separate(method, into=c("duration", "tagMethod"), remove=TRUE) %>% 
  dplyr::filter(tagMethod %in% c("1SPM", "1SPT")) %>% 
  mutate(duration = as.numeric(str_sub(duration, -100, -2))/60,
         distance = Inf) %>% 
  group_by(source, organization, project, sensor, tagMethod, equipment, location, buffer, lat, lon, year, date, observer, duration, distance, species, abundance, isSeen, isHeard, individual_appearance_order) %>%
  mutate(first_tag = min(tag_start_s)) %>%
  ungroup() %>%
  dplyr::filter(tag_start_s == first_tag) %>% 
  left_join(aru.wt.equip) %>% 
  dplyr::mutate(equipment = ifelse(is.na(equipmentshort), "unknown", equipmentshort),
                isSeen = "f",
                isHeard = "t") %>% 
  dplyr::select(all_of(colnms))

#2d. Replace TMTTs with predicted abundance----
tmtt <- read.csv("C:/Users/elly/Documents/ABMI/WildTrax/TMTT/data/tmtt_predictions_mean.csv") %>% 
  rename(species = species_code)

tmtt.wt <- aru.wt %>% 
  dplyr::filter(abundance=="TMTT") %>%
  mutate(species = ifelse(species %in% tmtt$species, species, "species"),
         observer_id = as.integer(ifelse(observer %in% tmtt$observer_id, observer, 0))) %>% 
  data.frame() %>% 
  left_join(tmtt) %>% 
  mutate(abundance = round(pred)) %>% 
  dplyr::select(colnames(aru.wt))

#2e. Put back together----
#summarize abundance
use.wt <- aru.wt %>% 
  dplyr::filter(abundance!="TMTT") %>% 
  rbind(tmtt.wt) %>% 
  rbind(pc.wt)

#3. Wrangle ebird data----
#Note this assumes observations with "X" individuals are 1s
#Remove hotspot data

ebd.files.done <- list.files(file.path(root, "eBird", "ebd_filtered"), pattern="ebd_*", full.names=TRUE)

#Note this next line takes a long time to run (couple hours)
raw.ebd <- purrr::map(.x=ebd.files.done, .f=~read_ebd(.)) %>% 
  rbindlist()

tax.wt <- read.csv(file.path(root, "lu_species.csv")) %>% 
  mutate(scientific_name = paste(species_genus, species_name)) %>% 
  rename(species = species_code) %>% 
  dplyr::select(scientific_name, species)

use.ebd <- raw.ebd %>% 
  dplyr::filter(locality_type!="H") %>% 
  mutate(source = "eBird",
         organization = "eBird",
         project= "eBird",
         sensor="PC",
         equipment="human",
         tagMethod="PC",
         singlesp="n",
         buffer=0,
         date = ymd_hms(paste0(observation_date, time_observations_started)),
         year = year(date),
         distance = Inf,
         abundance = as.numeric(ifelse(observation_count=="X", 1, observation_count)),
         isSeen=NA,
         isHeard=NA) %>% 
  rename(lat = latitude,
         lon = longitude,
         observer = observer_id,
         duration = duration_minutes,
         location = sampling_event_identifier) %>% 
  left_join(tax.wt) %>% 
  dplyr::select(all_of(colnms)) %>% 
  dplyr::filter(!is.na(date))
  
#D. PUT TOGETHER############################

#1. Put everything together----
use <- rbind(use.wt, use.ebd)

#2. Separate out locations only----
loc <- use %>% 
  dplyr::select(source, organization, project, sensor, location, lat, lon) %>% 
  unique() %>% 
  dplyr::filter(!is.na(lat))
  
#3. Clip by study area----

#3a. Read in study area----
sa <- read_sf("G:/Shared drives/BAM_NationalModels/NationalModels5.0/Regions/GEE_BufferedNatMod/GEE_BufferedNatMod.shp")

#3v. Create raster (much faster than from polygon)
r <- rast(ext(sa), resolution=1000)

sa.r <- rasterize(x=vect(sa), y=r, field="subUnit")

#3c. Extract raster locations and remove NAs----
loc.sa <- loc %>% 
  dplyr::select(lat, lon) %>% 
  st_as_sf(coords=c("lon", "lat"), crs=4326) %>% 
  st_transform(crs(sa)) %>% 
  vect() %>%
  terra::extract(x=sa.r) %>% 
  cbind(loc) %>% 
  dplyr::filter(!is.na(subUnit))

#E. REMOVE DUPLICATES----

#1. Set up grid----
grid <- dgconstruct(area=0.2, metric=TRUE)

#2. Identify locations that share a grid cell----
loc.grid <- loc.sa %>% 
  mutate(cell = dgGEO_to_SEQNUM(grid, lon, lat)$seqnum) 

n.grid <- loc.grid %>% 
  group_by(cell) %>% 
  summarize(n=n()) %>%
  ungroup() %>% 
  dplyr::filter(n > 1) %>% 
  left_join(loc.grid)

#3. Link to visit data and identify surveys that share duration, date time, total abundance----
visit.grid <- n.grid %>% 
  left_join(use) %>% 
  mutate(abundance = as.numeric(abundance)) %>% 
  dplyr::filter(!is.na(abundance)) %>% 
  group_by(source, project, sensor, location, lat, lon, cell, year, date, observer, duration, distance) %>%
  summarize(abundance = sum(abundance)) %>% 
  ungroup()

visit.n <- visit.grid %>% 
  group_by(cell, date, duration, abundance) %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  dplyr::filter(n > 1) %>% 
  left_join(visit.grid) %>% 
  unique()

#4. Take out all the ebird surveys in the duplicates----
visit.ebird.remove <- visit.n %>% 
  dplyr::filter(source=="eBird")

#5. Check again----
visit.n2 <- visit.grid %>%
  anti_join(visit.ebird.remove)  %>%
  group_by(cell, date, duration, abundance) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  dplyr::filter(n > 1) %>%
  left_join(visit.grid)

#6. Randomly select one to remove with same lat/lon----
set.seed(1234)
visit.n2.keep <- visit.n2 %>% 
  group_by(date, duration, abundance, lat, lon) %>% 
  sample_n(1) %>% 
  ungroup()

visit.n2.remove <- visit.n2 %>% 
  anti_join(visit.n2.keep)

#7. Check again----
visit.n3 <- visit.grid %>%
  anti_join(visit.ebird.remove)  %>%
  anti_join(visit.n2.remove) %>% 
  group_by(cell, date, duration, abundance) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  dplyr::filter(n > 1) %>%
  left_join(visit.grid)

#5. Put it all back together----
#Remove buffered locations
#Remove locations outside QPAD limits
dat <- use %>% 
  inner_join(loc.sa) %>% 
  anti_join(visit.ebird.remove %>% 
              dplyr::select(-abundance)) %>% 
  anti_join(visit.n2.remove %>% 
              dplyr::select(-abundance)) %>% 
  dplyr::select(all_of(colnms)) %>% 
  dplyr::filter(!buffer %in% c(50, 10000),
                lon >= -164,
                lon <= -52,
                lat >= 39,
                lat <= 69)

#F. SEPARATE INTO VISITS AND OBSERVATIONS####
#1. Identify unique visits----
#add primary key
visit <- dat %>% 
  dplyr::select(-species, -abundance, -isSeen, -isHeard) %>% 
  unique() %>% 
  mutate(id = row_number()) %>% 
  dplyr::filter(!is.na(duration),
                !is.na(distance))

#2. Tidy bird data----

#2a. Get list of species from qpad----
load_BAM_QPAD(4)
spp <- getBAMspecieslist()

#2b. Filter----
#remove unknown abundance
#filter to QPAD V4 species list
#link to primary key and make wide
bird <- dat %>% 
  mutate(abundance = as.numeric(abundance)) %>% 
  dplyr::filter(!is.na(abundance),
                !is.na(duration),
                !is.na(distance),
                species %in% spp) %>% 
  full_join(visit) %>% 
  mutate(species = ifelse(is.na(species), "NONE", species),
         abundance = ifelse(is.na(abundance), 1, abundance)) %>% 
  pivot_wider(id_cols=id, names_from=species, values_from=abundance, values_fn=sum, values_fill=0, names_sort=TRUE)

#G. SAVE!####
save(visit, bird, file=file.path(root, "01_NM5.0_data_clean.R"))
