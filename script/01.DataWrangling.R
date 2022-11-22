# ---
# title: National Models 4.1 - get data
# author: Elly Knight, Anna Drake
# created: Novemeber 17, 2022
# ---

#PREAMBLE####

#1. Load packages----

library(tidyverse)
library(wildRtrax)
library(data.table)
library(lubridate)

#2. Set root path for data on google drive----
root <- "G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/NationalModelsV4.1/PointCount/"

#A. DOWNLOAD DATA FROM WILDTRAX####

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
save(raw, projects, file=paste0("G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/QPAD/Data/wildtrax_raw_", Sys.Date(), ".Rdata"))

#B. GET BAM PATCH DATA####

#C. GET BC/YK PATCH DATA####
#D. GET EBIRD DATA####

#E. PUT TOGETHER####
