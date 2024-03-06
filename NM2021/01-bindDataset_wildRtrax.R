#####################################################################
# NOTES:
#   The script needs a config file to point to where credentials are kept. 
#
#   WT has duplicates on project name. Reports don't include organization to allow us to 
#   distinguish between the two differents projects.  It is thus not possible to merge summary with reports to recover sensor and status. 
#
#   Bugs on projects using special characters. I raised the issue to the WildTrax team. I fix it manually for now 
#
#   Some projects don't download: Error in x[[1]] : subscript out of bounds. project_id are  c(84, 446, 590, 805, 809, 814, 848, 854, 884, 913, 1007, 1084).
#   WildTrax has been informed. 
#
#  Lots of data uploaded on Wildtrax are workshop material, lab, class, course, demo dataset. 
#  I used only the ones published for now (!=Active). But that makes me loose many
#
#  Weyer are in the list. I think we can't use them for the National Model.
#
#  
#
#####################################################################
install.packages("remotes")
remotes::install_github("ABbiodiversity/wildRtrax")

library(utils)
library(stringr) #str_match
library(plyr) #rbind.fill
library(data.table)
library(readr)
library(tidyverse)
library(wildRtrax)

# Source wildTrax credential
config <- "E:/MelinaStuff/BAM/WildTrax/WT-Integration/config.R"
source(config)

# Authentify
wt_auth()

#Download the project summary
my_projects <- wt_get_download_summary(sensor_id = 'PC')
# Delete  from list status= Active
# Need to find a proper way to target dataset.
my_projects <- my_projects[!(my_projects$status=="Active"),]
nrow(my_projects) # n = 171 dataset

project_lid <-  sort(as.vector(unlist(my_projects$project_id)))


# Problem at download. Issue open on WT GitHub
notworking <- c(84, 446, 809)
# test project to delete from list

#Remove from project_list the ones that doesn't download  / test     
project_lid <- subset(project_lid, !(project_lid %in% notworking))

#update project_list
my_projects <-subset(my_projects,my_projects$project_id %in% project_lid) #n = 159

# Loop using chunk to avoid memory issue
chunk_length <- 40  
lll <- split(project_lid, ceiling(seq_along(project_lid) / chunk_length))
i<-0
mylist<-list()
for (loop in lll) {
  i <- i+1
  lilo <- lapply(loop, function(x) {
    print(x)
    #f_wt <- read.csv(file.path(wdir, wt_data, x), sep=",", header = TRUE, stringsAsFactors = FALSE, encoding = "UTF-8")
    f_wt <- wt_download_report(project_id = x, sensor_id = 'PC', cols_def = F, weather_cols = F)
    return(f_wt)
  })
  m <-do.call(rbind, lilo) 
  mylist[[i]]<-m #add each element to list
}
# Bind data 
wt_data <-do.call(rbind, mylist) 

# delete project = NA
wt_data <- wt_data[!(is.na(wt_data$project)),]


# Fix project name in report to match project name in summary
# *****MOST OF THEM HAVE STATUS = ACTIVE and won't download if we filter using status
wt_data["project"][wt_data["project"] == "EdÃ©hzhÃ­e National Wildlife Area CWS Northern Region 2016"]  <-"Edéhzhíe National Wildlife Area CWS Northern Region 2016"
wt_data["project"][wt_data["project"] == "EdÃ©hzhÃ­e National Wildlife Area CWS Northern Region 2019"] <-"Edéhzhíe National Wildlife Area CWS Northern Region 2019"  
wt_data["project"][wt_data["project"] == "TÅ‚Ä±Ì¨chÇ« Winter Road CWS Northern Region 2019"] <- "Tchicho Winter Road CWS Northern Region 2019"  
wt_data["project"][wt_data["project"] == "Tsâ€™udÃ© NilÄ¯nÃ© Tuyeta PA CWS Northern Region 2020"] <-"Ts’udé Niliné Tuyeta PA CWS Northern Region 2020"
wt_data["project"][wt_data["project"] == "ForÃªt Montmorency long-term bird survey 2000"] <- "Forêt Montmorency long-term bird survey 2000"   
wt_data["project"][wt_data["project"] == "RÃ©gularisation du Lac KÃ©nogami EIA 2001"] <-"Régularisation du Lac Kénogami EIA 2001"
wt_data["project"][wt_data["project"] == "Tsâ€™udÃ© NilÄ¯nÃ© Tuyeta PA CWS Northern Region 2021"]  <-"Ts’udé Niliné Tuyeta PA CWS Northern Region 2021"
wt_data["project"][wt_data["project"] == "NÃ¡Ã¡ts'Ä¯hch'oh NPR 2018 CWS Northern Region"]  <- "Nááts'ihch'oh NPR 2018 CWS Northern Region"           
wt_data["project"][wt_data["project"] == "NÃ¡Ã¡ts'Ä¯hch'oh NPR 2019 CWS Northern Region"]  <- "Nááts'ihch'oh NPR 2019 CWS Northern Region"         
wt_data["project"][wt_data["project"] == "EdÃ©hzhÃ­e National Wildlife Area CWS Northern Region 2021"] <-"Edéhzhíe National Wildlife Area CWS Northern Region 2021" 
wt_data["project"][wt_data["project"] == "Thaidene NÃ«nÃ© PA CWS Northern Region 2022"]<- "Thaidene Nëné PA CWS Northern Region 2022"

# Fix project name in summary to match project name in report
# *****MOST OF THEM HAVE STATUS = ACTIVE and won't download if we filter using status
my_projects["project"][my_projects["project_id"] == 325] <- "Tchicho Winter Road CWS Northern Region 2019"
my_projects["project"][my_projects["project_id"] == 432] <- "Ts’udé Niliné Tuyeta PA CWS Northern Region 2020"
my_projects["project"][my_projects["project_id"] == 856] <- "Ts’udé Niliné Tuyeta PA CWS Northern Region 2021"
my_projects["project"][my_projects["project_id"] == 977] <- "Nááts'ihch'oh NPR 2018 CWS Northern Region"
my_projects["project"][my_projects["project_id"] == 978] <- "Nááts'ihch'oh NPR 2019 CWS Northern Region"


# Should be none
unique(subset(wt_data$project, !(wt_data$project %in% my_projects$project)))

# Recover sensor and status from summary
wt_data <- merge(wt_data, my_projects, by = "project")


#------- STOP merge not working. organization not found in summary
###

# Set local tempdir where patch data are found.
wdir <- "E:/MelinaStuff/BAM/WildTrax/WTdownload/"

# PATCH (i.e., not loaded yet. Issues open on WildTrax)
patch_data <- "patch"
patchList <- list.files(file.path(wdir, patch_data))
patchList <-unique(sub("_[^_]*$", "", patchList))
patch_report <- lapply(patchList, function(x) {
  print (x)
  location_patch <- read.csv(file.path(wdir, patch_data, paste0(x, "_location.csv")), sep=",", header = TRUE, stringsAsFactors = FALSE, encoding = "UTF-8")
  location_patch$sensor <- "PC"
  location_patch$status <- "Active"
  names(location_patch)[names(location_patch) == "comments"] <- "comments_location"
  visit_patch <- read.csv(file.path(wdir, patch_data, paste0(x, "_visit.csv")), sep=",", header = TRUE, stringsAsFactors = FALSE, encoding = "UTF-8")
  names(visit_patch)[names(visit_patch) == "comments"] <- "comments_visit"
  survey_patch <- read.csv(file.path(wdir, patch_data, paste0(x, "_survey.csv")), sep=",", header = TRUE, stringsAsFactors = FALSE, encoding = "UTF-8")
  names(survey_patch)[names(survey_patch) == "comments"] <- "comments_survey"
  extended_patch <- read.csv(file.path(wdir, patch_data, paste0(x, "_extended.csv")), sep=",", header = TRUE, stringsAsFactors = FALSE, encoding = "UTF-8", row.names= NULL)
  lv_patch <- merge(location_patch[ , c("location", "latitude", "longitude", "bufferRadiusMeters", "sensor", "status", "comments_location")], visit_patch[ , c("location", "visitDate", "comments_visit")], by = "location")
  survey_patch$visitDate <- sub("\\ .*", "", survey_patch$surveyDateTime)
  lvs_patch <- merge(lv_patch, survey_patch, by = c("location","visitDate"))
  #extended_field <- c("organization",	"project",	"location",	"surveyDateTime",	"species",	"distanceband",	"durationinterval")
  #full_patch <- merge(extended_patch[ , extended_field], lvs_patch, by = c("location",	"surveyDateTime",	"species",	"distanceband",	"durationinterval"))
  lvs_patch$comments <- paste(lvs_patch$comments_location, lvs_patch$comments_survey, lvs_patch$comments_survey, sep = ";")
  #f_aru <- read_csv(fi, sep=",")
  return(lvs_patch)
})

patch_bind <-do.call(rbind, patch_report) # 

#Fix colnames to fit WTdownload
names(patch_bind)[names(patch_bind) == "bufferRadiusMeters"] <- "bufferRadius.m."
names(patch_bind)[names(patch_bind) == "surveyDateTime"] <- "date"
names(patch_bind)[names(patch_bind) == "species"] <- "speciesCode"
names(patch_bind)[names(patch_bind) == "distanceband"] <- "distanceBand"
names(patch_bind)[names(patch_bind) == "durationinterval"] <- "durationInterval"
patch_bind$comments <- gsub("NA;", "", patch_bind$comments)
patch_bind$comments <- gsub(";NA", "", patch_bind$comments)

# BIND wt_download and patch
wt_bind <-rbind.fill(wt_data, patch_bind) 
col_order <-  c("organization", "project", "sensor", "status", "location", "latitude", "longitude", "date", 
               "distanceMethod", "durationMethod", "bufferRadius.m.", "distanceBand", "durationInterval", 
               "speciesCode", "speciesCommonName", "abundance", "isHeard", "isSeen", "observer", "scientificName", 
               "daily_weather_station_nm", "daily_weather_station_elevation", "daily_weather_station_distance", 
               "daily_min_temp", "daily_max_temp", "daily_mean_temp", "daily_total_rain_mm", "daily_total_snow_cm", 
               "daily_precipitation_mm", "daily_snow_on_ground_cm", "hourly_weather_station_nm", 
               "hourly_weather_station_elevation", "hourly_weather_station_distance", "hourly_temp", 
               "hourly_dew_point", "hourly_rel_humidity", "hourly_precipitation_mm", "hourly_wind_direction", 
               "hourly_wind_speed", "hourly_visibility_km","hourly_station_pressure", "hourly_humidex", 
               "hourly_wind_chill", "hourly_weather_attributes", "comments")       

wt_bind <- wt_bind[, col_order]
# Save
save(wt_bind, file="./01_wtBind.RData")
## END
##############################################################################################################################


#####################
#--  CHECK using  XY mapping 
#####################
library(sp)
library(sf)

wt_geo <-wt_bind[!(is.na(wt_bind$latitude)),]
wt_geo <-wt_geo[!(is.na(wt_geo$longitude)),]

# Projection
DD <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

coordinates(wt_geo) <- c("longitude", "latitude")
proj4string(wt_geo) <- DD

pc <- st_as_sf(wt_geo, coords = c("longitude", "latitude"), crs = DD)


## Map Studyregion
f_studyregion <- "E:/MelinaStuff/BAM/request/AnnaDrake/SampleRegion/ClimateSamplingRegion_AnnaDrake.shp"
studyregion <- st_read(f_studyregion)
plot(st_geometry(studyregion))

## Add point
sensor <- pc$sensor
plot(st_geometry(pc), pch  = 20, col=ifelse(sensor=="ARU", "red",ifelse(sensor== "PC", "black", "blue")), add= TRUE)

## Map canada extent
f_canada <-"E:/MelinaStuff/BAM/NationalModelv4.1/runNationalModel/data/basemaps/CanadaLAEA.shp"
canada<- st_read(f_canada)
canada_DD <-st_transform(canada, DD)
plot(st_geometry(canada_DD), border="grey", add = TRUE)
plot(st_geometry(studyregion), lwd=3, add = TRUE)

pc_shp <- st_as_sf(pc, coords = c("longitude", "latitude"), crs = DD)

st_write(pc_export, "E:/MelinaStuff/BAM/request/AnnaDrake/wt_location1.shp")
st_write(studyregion, "E:/MelinaStuff/BAM/request/AnnaDrake/studyregion.shp")


