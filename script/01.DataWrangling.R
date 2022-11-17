# ---
# title: National Models 4.1 - get data
# author: Elly Knight, Anna Drake
# created: Novemeber 17, 2022
# ---

library(tidyverse)
library(wildRtrax)
library(data.table)
library(RODBC)
library(lubridate)

#NOTE: THIS IS INTENDED FOR QPAD DATA WRANGLING; NEEDS UPDATES####

#A. DOWNLOAD DATA FROM WILDTRAX####

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
for(i in 1:nrow(projects)){
  
  if(projects$sensorId[i]=="ARU"){
    dat.try <- try(wt_download_report(project_id = projects$project_id[i], sensor_id = projects$sensorId[i], weather_cols = F, report = "summary"))
  }
  
  if(projects$sensorId[i]=="PC"){
    dat.try <- try(wt_download_report(project_id = projects$project_id[i], sensor_id = projects$sensorId[i], weather_cols = F, report="report"))
  }
  
  if(class(dat.try)=="data.frame"){
    dat.list[[i]] <- dat.try
  }
  
  print(paste0("Finished dataset ", projects$project[i], " : ", i, " of ", nrow(projects), " projects"))
  
}

#4. Collapse list----
raw <- rbindlist(dat.list, fill=TRUE) %>% 
  mutate(project = ifelse(is.na(project), project_name, project),
         speciesCode = ifelse(is.na(speciesCode), species_code, speciesCode)) %>% 
  dplyr::select(-project_name)

#5. Save date stamped data & project list----
save(raw, projects, file=paste0("G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels4.1/Data/wildtrax_raw_", Sys.Date(), ".Rdata"))

#B. COMPARE TO V6 DATABASE####

#WE MIGHT NOT NEED THIS IF THE PATCH IS ADEQUATE####

#1. Load wildtrax dataset----
load("G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels4.1/Data/wildtrax_raw_2022-11-06.Rdata")

#2. Connect to BAM v6 dataset----
v6 <- odbcConnectAccess2007("G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/DataStuff/Avian.Data/BAM-V6/BAM-V6-USE.accdb")

#3. Get list of v6 projects----
v6.tbls <- sqlTables(v6)
v6.projects <- data.frame(PCODE = unique(sqlFetch(v6, "location")$dataset_fk_V4)) %>% 
  mutate(V6 = 1)

#4. Get list of wildtrax projects----
#remove ones that are not complete in wildtrax
load("G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels4.1/Data/BAMpatch.Rdata")

wt.projects <- data.frame(project = unique(raw$project)) %>% 
  dplyr::filter(!project %in% patch_bind$project) %>% 
  mutate(WT=1)

#5. Read in lookup table----
lu.projects <- read.csv("G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/QPAD/Data/VersionDatasetComparison_MHrevised.csv") %>% 
  dplyr::select(project, PCODE)

#6. Evaluate project overlap----
wtv6.projects <- wt.projects %>% 
  full_join(lu.projects) %>% 
  full_join(v6.projects)

#7. Projects not in WT----
missing.projects <- wtv6.projects %>% 
  dplyr::filter(V6==1 & is.na(WT))

#C. BAM V6 DATA####

#1. Get data for missing projects----
v6.loc <- sqlFetch(v6, "location") %>% 
  dplyr::filter(dataset_fk_V4 %in% missing.projects$PCODE)
v6.visit <- sqlFetch(v6, "pc_visit") %>% 
  dplyr::filter(location_name_4 %in% v6.loc$SS_V4)
v6.bird <- sqlFetch(v6, "pc_detection") %>% 
  dplyr::filter(PKEY %in% v6.visit$PKEY_V4)

#2. Format to match wt----
v6.raw <- v6.bird %>% 
  dplyr::select(PKEY, species_code, duration_interval, distance_band, abundance, heard, seen) %>% 
  rename(visit = PKEY, speciesCode = species_code, isHeard = heard, isSeen = seen, distance_band_numid = distance_band) %>% 
  mutate(isHeard = case_when(isHeard=="No" ~ "f",
                             isHeard=="Yes" ~ "t",
                             !is.na(isHeard) ~ "",
                             is.na(isHeard) ~ ""),
         isSeen = case_when(isSeen=="No" ~ "f",
                            isSeen=="Yes" ~ "t",
                            !is.na(isSeen) ~ "",
                            is.na(isSeen) ~ "")) %>% 
  inner_join(v6.visit %>% 
               mutate(date = ymd_hms(paste(as.character(survey_date), str_sub(as.character(survey_time), 12, 19)))) %>% 
               dplyr::select(PKEY_V4, location_name_4, method_duration, method_distance, date, observer_id) %>% 
               rename(visit = PKEY_V4, location = location_name_4, protocol_distance_numid = method_distance, protocol_duration_id = method_duration, observer = observer_id)) %>% 
  inner_join(v6.loc %>% 
               dplyr::select(SS_V4, latitude, longitude, dataset_fk_V4, location_buffer_radius_m) %>% 
               rename(location = SS_V4, project = dataset_fk_V4, buffer = location_buffer_radius_m)) %>% 
  dplyr::select(-visit) %>% 
  mutate(organization=NA)

#3. Replace distance & duration codes with text----
durmeth <- sqlFetch(v6, "lu_pc_protocol_duration") %>% 
  rename(durationMethod = protocol_duration_range)
dismeth <- sqlFetch(v6, "lu_pc_protocol_distance") %>% 
  rename(distanceMethod = protocol_distance_range)
durint <- sqlFetch(v6, "lu_pc_duration_interval") %>% 
  rename(durationInterval = duration_description)
disint <- sqlFetch(v6, "lu_pc_distance_band") %>% 
  rename(distanceBand = distance_band_description)

patch.raw <- v6.raw %>% 
  left_join(durmeth) %>% 
  left_join(dismeth) %>% 
  left_join(durint) %>% 
  left_join(disint) %>% 
  dplyr::select(-protocol_duration_id, -protocol_distance_numid, -duration_interval, -distance_band_numid) %>% 
  mutate(dataset = "BAM",
         method = "PC",
         tag_start_s = NA,
         individual_appearance_order = NA,
         ARUmethod = NA)

#D. PUT DATASETS TOGETHER & TIDY####

#1. Put together----
raw.all <- raw %>% 
  mutate(dataset = "WT",
         date = ymd_hms(ifelse(is.na(date), recording_date, date))) %>% 
  rename(buffer = bufferRadius.m.,
         ARUmethod = method) %>% 
  left_join(projects %>% 
              dplyr::select(project, sensorId) %>% 
              rename(method = sensorId)) %>% 
  dplyr::select(colnames(patch.raw)) %>% 
  rbind(patch.raw)

#2. Remove datasets without permission----
remove <- c("Weyerhaeuser Point Counts Westworth 1995",
            "Weyerhaeuser Point Counts Ursus Ecosystem Management 2001, 2004, 2007",
            "Weyerhaeuser Point Counts STRIX Ecological 2006 to 2020",
            "Weyerhaeuser Point Counts Aspen Ecological 2001, 2004", 
            "WHPC")

use <- raw.all %>% 
  dplyr::filter(!project %in% remove)

save(use, file="G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/NationalModels4.1/Data/NationalModels4.1v4_dat_2022-11-06.Rdata")


#date vs recording_date
#project vs project_name
#speciesCode vs species_code