# ---
# title: National Models 4.1 - get data
# author: Elly Knight, Melina Houle, Anna Drake
# created: Novemeber 17, 2022
# ---

#NOTES################################

#TO DO: GET EBIRD FOR OUTSIDE CANADA####
#TO DO: GET CALLING LAKE DATA#####

#The "BAMProjects_WildTrax.csv" file is a list of all projects currently in WildTrax that BAM can use in the national models. This file should be updated for each iteration of the national models in collaboration with Erin Bayne. Future versions of this spreadsheet can hopefully be derived by a combination of organization and a google form poll for consent from other organizations.

#The "BAMProjects_WildTrax.csv" file also contains information on which ARU projects are processed for a single species or taxa and therefore those visits should only be used for models for the appropriate taxa. This file should be updated for each iteration of the national models in collaboration with Erin Bayne.

#There are a handful of projects that are not downloading properly via wildRtrax. An issue is open on this. These projects are listed in the error.log object.

#BAM patch was compiled by Melina Houle to include datasets that are not yet loaded into WildTrax. Future iterations of the national models should not need this patch, with the exception of the BBS data which is likely too large to ever be uploaded in WT.

#BC/YK patch was provided by Anna Drake to include datasets from these regions that are unavailable in the BAM database. These datasets will be harmonized for inclusion in WildTrax under BAM and so future iterations of the national models should not need this patch.

#raw eBird data is downloaded from the eBird interface at https://ebird.org/data/download/ebd prior to wrangling with the "auk" package and will require a request for access. Use the custom download tool to download only the datasets for Canada and the US instead of the global dataset. Note you will also need the global sampling file to use the auk package for zero filling.

#raw eBird data omits Great Grey Owl & Northern Hawk Owl as sensitive species (https://support.ebird.org/en/support/solutions/articles/48000803210?b_id=1928&_gl=1*xq054u*_ga*ODczMTUyMjcuMTY2OTE0MDI4Ng..*_ga_QR4NVXZ8BM*MTY2OTE0MDI4NS4xLjEuMTY2OTE0MDM3OC4zNS4wLjA.&_ga=2.147122167.150058226.1669140286-87315227.1669140286) and should not be used for modelling these two species.

#wrangling eBird data with the auk package requires installation of AWK on windows computers. Please see #https://cornelllabofornithology.github.io/auk/articles/auk.html.

#eBird data has not been zerofilled because there was no species filtering done and we are assuming that all stationary counts have at least 1 bird observed.

#The column "sensor" currently only differentiates between ARU & human point count data types. Future versions should consider differentiating between SM2 ARU data types and other ARU data types due to differences in the perceptibility of these two approaches, either via QPAD or a correction factor.

#The replace TMTTs script will be replaced by a wildRtrax function in the near future.

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(wildRtrax) #to download data from wildtrax
library(data.table) #for binding lists into dataframes
library(lubridate) #date wrangling
library(auk) #eBird wrangling

#2. Set root path for data on google drive----
root <- "G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/NationalModelsV4.1/PointCount/"

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

#3. Filter by list of projects we can use----
#Fix names to match report
project.names <- data.frame(project_id = c(325, 432, 856, 977, 978),
                            newname = c("Tłı̨chǫ Winter Road CWS Northern Region 2019",
                                        "Ts’udé Nilįné Tuyeta PA CWS Northern Region 2020",
                                        "Ts’udé Nilįné Tuyeta PA CWS Northern Region 2021",
                                        "Nááts'įhch'oh NPR 2018 CWS Northern Region",
                                        "Nááts'įhch'oh NPR 2019 CWS Northern Region"))

use <- read.csv(file.path(root, "BAMProjects_WildTrax.csv"))

projects.wt <- projects %>% 
  dplyr::filter(project_id %in% use$project_id) %>% 
  left_join(project.names) %>% 
  mutate(project = ifelse(is.na(newname), project, newname)) %>% 
  dplyr::select(-newname)

#4. Loop through projects to download data----
dat.list <- list()
error.log <- data.frame()
for(i in 1:nrow(projects.wt)){
  
  #Do each sensor type separately because the reports have different columns and we need different things for each sensor type
  if(projects.wt$sensorId[i]=="ARU"){
    dat.try <- try(wt_download_report(project_id = projects.wt$project_id[i], sensor_id = projects.wt$sensorId[i], weather_cols = F, report = "summary"))
  }
  
  if(projects.wt$sensorId[i]=="PC"){
    dat.try <- try(wt_download_report(project_id = projects.wt$project_id[i], sensor_id = projects.wt$sensorId[i], weather_cols = F, report="report"))
  }
  
  if(class(dat.try)=="data.frame"){
    dat.list[[i]] <- dat.try
  }
  
  #Log projects that error
  if(class(dat.try)!="data.frame"){
    error.log <- rbind(error.log, 
                       projects.wt[i,])
    
  }
  
  print(paste0("Finished dataset ", projects.wt$project[i], " : ", i, " of ", nrow(projects.wt), " projects"))
  
}

#5. Collapse list----
#standardize column names between sensor types
#fix special character project names
report.names <- data.frame(project=c("EdÃ©hzhÃ­e National Wildlife Area CWS Northern Region 2016",
                                        "EdÃ©hzhÃ­e National Wildlife Area CWS Northern Region 2019",
                                        "TÅ‚Ä±Ì¨chÇ« Winter Road CWS Northern Region 2019",
                                        "Tsâ€™udÃ© NilÄ¯nÃ© Tuyeta PA CWS Northern Region 2020",
                                        "ForÃªt Montmorency long-term bird survey 2000",
                                        "RÃ©gularisation du Lac KÃ©nogami EIA 2001",
                                        "Tsâ€™udÃ© NilÄ¯nÃ© Tuyeta PA CWS Northern Region 2021",
                                        "NÃ¡Ã¡ts'Ä¯hch'oh NPR 2018 CWS Northern Region",
                                        "NÃ¡Ã¡ts'Ä¯hch'oh NPR 2019 CWS Northern Region",
                                        "EdÃ©hzhÃ­e National Wildlife Area CWS Northern Region 2021",
                                        "Thaidene NÃ«nÃ© PA CWS Northern Region 2022"),
                            newname=c("Edéhzhíe National Wildlife Area CWS Northern Region 2016",
                                         "Edéhzhíe National Wildlife Area CWS Northern Region 2019",
                                         "Tchicho Winter Road CWS Northern Region 2019",
                                         "Ts’udé Niliné Tuyeta PA CWS Northern Region 2020",
                                         "Forêt Montmorency long-term bird survey 2000",
                                         "Régularisation du Lac Kénogami EIA 2001",
                                         "Ts’udé Niliné Tuyeta PA CWS Northern Region 2021",
                                         "Nááts'ihch'oh NPR 2018 CWS Northern Region",
                                         "Nááts'ihch'oh NPR 2019 CWS Northern Region",
                                         "Edéhzhíe National Wildlife Area CWS Northern Region 2021",
                                         "Thaidene Nëné PA CWS Northern Region 2022"))

raw.wt <- rbindlist(dat.list, fill=TRUE) %>%
  mutate(project = ifelse(is.na(project), project_name, project),
         speciesCode = ifelse(is.na(speciesCode), species_code, speciesCode),
         date = ifelse(is.na(date), recording_date, date)) %>%
  dplyr::select(-project_name, -recording_date, -species_code) %>% 
  left_join(report.names) %>% 
  mutate(project = ifelse(is.na(newname), project, newname)) %>% 
  dplyr::select(-newname)

#7. Save date stamped data & project list----
save(raw.wt, projects.wt, error.log, file=paste0(root, "/wildtrax_raw_", Sys.Date(), ".Rdata"))
load(file.path(root, "wildtrax_raw_2022-11-22.Rdata"))

#B. GET PATCH DATA###############################

#1. BAM patch----
raw.bam <- readRDS(file.path(root, "pc_patch.rds"))

#2. BC/YK patch----
raw.bcyk <- readRDS(file.path(root, "BC_YK_patch.rds"))

#C. GET EBIRD DATA##########################

#Note: loop currently has to be run by hand because auk_set_ebd_path() requires restarting the session after running the command

library(auk)
root <- "G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/NationalModelsV4.1/PointCount/"

#1. Get list of ebd objects to process----
ebd.files <- list.files(file.path(root, "ebd_raw"), pattern="ebd_*")

i <- 24
#2. Set up loop----
for(i in 3:length(ebd.files)){
  
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

#3. Check list of processed files----
ebd.files.done <- list.files(file.path(root, "ebd_filtered"), pattern="ebd_*")


#D. HARMONIZE###############################

#1. Set desired columns----
colnms <- c("source", "project", "sensor", "singlesp", "location", "buffer", "lat", "lon", "year", "date", "observer", "duration", "distance", "species", "abundance")

#2. Wrangle wildtrax data-----

#2a. A bit of prep----
dat.wt <- raw.wt %>% 
  rename(lat = latitude, lon = longitude, species = speciesCode, buffer = bufferRadius.m.) %>% 
  full_join(projects.wt %>% 
              rename(sensor = sensorId) %>% 
              dplyr::select(project_id, project, sensor)) %>% 
  mutate(source = "WildTrax", 
         date = ymd_hms(date),
         year = year(date))

#2b. Get the point count data----
#wrangle distance and duration maximums
#remove counts with unknown duration and distance
pc.wt <- dat.wt %>% 
  dplyr::filter(sensor=="PC",
                durationMethod!="UNKNOWN",
                distanceMethod!="UNKNOWN") %>% 
  rowwise() %>% 
  mutate(durationMethod = ifelse(str_sub(durationMethod, -1, -1)=="+", str_sub(durationMethod, -100, -2), durationMethod),
         chardur = str_locate_all(durationMethod, "-"),
         chardurmax = max(chardur),
         duration = as.numeric(str_sub(durationMethod, chardurmax+1, -4)),
         chardis = str_locate_all(distanceMethod, "-"),
         chardismax = max(chardis),
         distance1 = str_sub(distanceMethod, chardismax+1, -2),
         distance = ifelse(distance1 %in% c("AR", "IN"), Inf, as.numeric(distance1)),
         singlesp = "n") %>% 
  dplyr::select(colnms) %>% 
  data.frame()

#2c. Get the aru data----
#filter ARU data to first detection of each individual
#wrangle duration
#identify datasets that should only be used for single species

ssp <- read.csv(file.path(root, "BAMProjects_WildTrax.csv")) %>% 
  dplyr::filter(single.species=="y")

aru.wt <- dat.wt %>% 
  dplyr::filter(sensor=="ARU") %>% 
  mutate(singlesp = ifelse(project_id %in% ssp$project_id, "y", "n")) %>% 
  separate(method, into=c("duration", "method"), remove=TRUE) %>% 
  mutate(duration = as.numeric(str_sub(duration, -100, -2)),
         distance = Inf) %>% 
  group_by(source, project, sensor, singlesp, location, buffer, lat, lon, year, date, observer, duration, distance, species, abundance, individual_appearance_order) %>%
  mutate(first_tag = min(tag_start_s)) %>%
  ungroup() %>%
  dplyr::filter(tag_start_s == first_tag) %>% 
  dplyr::select(all_of(colnms))

#2d. Replace TMTTs with predicted abundance----
tmtt <- read.csv("C:/Users/Elly Knight/Documents/ABMI/Projects/TMTT/data/tmtt_predictions.csv") %>% 
  rename(species = species_code)
user <- read.csv("C:/Users/Elly Knight/Documents/ABMI/Projects/TMTT/data/app_user.csv") %>% 
  rename(observer = user_name) %>% 
  dplyr::select(observer, user_id)

tmtt.wt <- aru.wt %>% 
  dplyr::filter(abundance=="TMTT") %>% 
  left_join(user) %>% 
  mutate(user_id = ifelse(is.na(user_id), observer, user_id))%>% 
  mutate(boot = round(runif(max(row_number()), 1, 100)),
         species = ifelse(species %in% tmtt$species, species, "species"),
         user_id = as.integer(ifelse(user_id %in% tmtt$user_id, user_id, 0))) %>% 
  data.frame() %>% 
  left_join(tmtt) %>% 
  mutate(abundance = round(pred)) %>% 
  dplyr::select(colnames(aru.wt))

#2e. Put back together----
use.wt <- aru.wt %>% 
  dplyr::filter(abundance!="TMTT") %>% 
  rbind(tmtt.wt)

#3. Wrangle BAM patch data----
#wrangle distance and duration maximums
#remove counts with odd duration method entries
#replace all unknown dates (MN-BBATLAS, NEFBMP2012-19) with June 15 of 2012
use.bam <- raw.bam %>% 
  dplyr::filter(!durationMethod %in% c("", " during the 10 minutes.", " seemed to bring food then brooded; assume nestlings stil", " timeperiod C")) %>% 
  rowwise() %>% 
  mutate(source = "BAM",
         sensor = "PC",
         singlesp = "n",
         surveyDateTime = ifelse(visitDate=="", paste0("2012-06-15", surveyDateTime), surveyDateTime),
         date = ymd_hms(surveyDateTime),
         year = year(date),
         observer = NA,
         durationMethod = ifelse(str_sub(durationMethod, -1, -1)=="+", str_sub(durationMethod, -100, -2), durationMethod),
         chardur = str_locate_all(durationMethod, "-"),
         chardurmax = max(chardur),
         duration = as.numeric(str_sub(durationMethod, chardurmax+1, -4)),
         chardis = str_locate_all(distanceMethod, "-"),
         chardismax = max(chardis),
         distance1 = str_sub(distanceMethod, chardismax+1, -2),
         distance = ifelse(distance1 %in% c("AR", "IN"), Inf, as.numeric(distance1))) %>% 
  rename(buffer = bufferRadiusMeters, lat = latitude, lon = longitude) %>% 
  data.frame() %>% 
  dplyr::select(all_of(colnms))

#4. Wrangle BC/YT patch data----
use.bcyk <- raw.bcyk %>% 
  rename(project = Source, duration = dur, distance = dis, species = spp, abundance = count) %>% 
  mutate(source = "BCYK",
         sensor = ifelse(str_sub(project, -3, -1)=="ARU", "ARU", "PC"),
         singlesp = "n",
         location = NA,
         buffer = NA,
         date = ymd_hm(paste0(dt, tm)),
         year = year(date),
         observer = NA) %>% 
  dplyr::select(all_of(colnms))

#5. Wrangle ebird data----
#Note this assumes observations with "X" individuals are 1s

raw.ebd <- read_ebd("ebd_data_filtered.txt")

use.ebd <- raw.ebd %>% 
  mutate(date = ymd_hms(paste0(observation_date, time_observations_started)),
         abundance = as.numeric(ifelse(observation_count=="X", 1, observation_count)),
         route = NA,
         stop = NA) %>% 
  separate(state_code, into=c("country", "province")) %>% 
  rename(lat = latitude,
         lon = longitude,
         duration = duration_minutes) %>% 
  dplyr::select(year, province, route, stop, date, count, lat, lon, time, day, duration)
  
#E. PUT TOGETHER############################

#1. Put everything together----
use <- rbind(use.wt, use.bam, use.bcyk, use.ebd)

#2. Separate into visit and detection objects----
  
#3. Clip by study area----

#4. Remove duplicates----
#dggridR::
#setup for 200m, flag, then check for identical timestamps


#F. SAVE!####
