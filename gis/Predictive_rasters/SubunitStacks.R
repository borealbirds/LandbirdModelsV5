# ---
# title: National Models 5.0 - stack and crop rasters for subunit prediction
# author: Elly Knight
# created: March 5, 2024
# ---

#NOTES################################

#This script creates separate raster stacks for each of the modelling subunits for prediction.

#First we create use the bcr shapefile and canada/USA shapefiles to make a polygon for each subunit with a 100 km buffer. We also do this for regions that we have mreged.

#Next, we identify the covariates that we use for that region, stack them, and then crop them by the polygon for that region. We save out each stack separately.

#The stacks can then be moved over to compute canada for model prediction.

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #basic data wrangling
library(terra) #basic raster handling
library(sf) #basic shapefile handling
library(nngeo) #fill holes in polygons

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_NationalModels/NationalModels5.0"

#3. Load data packages with covariate lookup table----
load(file.path(root, "Data", "04_NM5.0_data_stratify.R"))

#SUBUNIT POLYGONS#####################

#1. Read in country shapefiles----
can <- read_sf(file.path(root, "Regions", "CAN_adm", "CAN_adm0.shp")) %>% 
  st_transform(crs=5072) 
usa <- read_sf(file.path(root, "Regions", "USA_adm", "USA_adm0.shp")) %>% 
  st_transform(crs=5072) 

#2. Read in BCR shapefile----
#Remove subunit 1
bcr <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel.shp")) %>% 
  dplyr::filter(subUnit!=1) %>% 
  st_transform(crs=5072)

ggplot(bcr) +
  geom_sf(aes(fill=factor(subUnit)))

#3. Identify BCRs for each country----
bcr.ca <- bcr %>% 
  st_intersection(can) %>% 
  mutate(country="can") %>% 
  dplyr::select(subUnit, country)

bcr.usa <- bcr %>% 
  st_intersection(usa) %>% 
  mutate(country="usa") %>% 
  dplyr::select(subUnit, country)

#4. Make merged subunits----
bcr.can4142 <- bcr.ca %>% 
  dplyr::filter(subUnit %in% c(41, 42)) %>% 
  st_union() %>% 
  nngeo::st_remove_holes() %>% 
  st_sf() %>% 
  mutate(country="can",
         subUnit = 4142)

bcr.can41423 <- bcr.ca %>% 
  dplyr::filter(subUnit %in% c(41, 42, 3)) %>% 
  st_union() %>% 
  nngeo::st_remove_holes() %>% 
  st_sf() %>% 
  mutate(country="can",
         subUnit = 41423)

bcr.usa414232 <- bcr.usa %>% 
  dplyr::filter(subUnit %in% c(41, 42, 3, 2)) %>% 
  st_union() %>% 
  nngeo::st_remove_holes() %>% 
  st_sf() %>% 
  mutate(country="usa",
         subUnit = 414232)

bcr.can8182 <- bcr.ca %>% 
  dplyr::filter(subUnit %in% c(81, 82)) %>% 
  st_union() %>% 
  nngeo::st_remove_holes() %>% 
  st_sf() %>% 
  mutate(country="can",
         subUnit = 8182)

#polygons to remove
bcr.remove <- data.frame(country=c("can", "usa", "usa", "can", "usa"),
                         subUnit=c(41, 42, 3, 42, 41))

#5. Put together----
bcr.country <- rbind(bcr.ca, bcr.usa, bcr.can4142, bcr.can41423, bcr.usa414232, bcr.can8182) %>% 
  anti_join(bcr.remove)

#6. Set up loop for BCR buffering----
bcr.out <- data.frame()
for(i in 1:nrow(bcr.country)){
  
  #7. Filter & buffer shapefile----
  bcr.buff <- bcr.country %>% 
    dplyr::filter(row_number()==i) %>% 
    st_buffer(100000)
  
  #8. Crop to international boundary----
  if(bcr.buff$country=="can"){ bcr.i <- st_intersection(bcr.buff, can)}
  if(bcr.buff$country=="usa"){ bcr.i <- st_intersection(bcr.buff, usa)}
  
  #9. Put together----
  bcr.out <- rbind(bcr.out, bcr.i %>% 
                     dplyr::select(country, subUnit))
  
  print(paste0("Finished bcr ", i, " of ", nrow(bcr.country)))
  
}

#10. Save----
write_sf(bcr.out, file.path(root, "Regions", "BAM_BCR_NationalModel_Buffered.shp"))

#GET PREDICTION RASTER LOCATIONS######

#1. Get list of files in the folder----
#Remove Archived files and a couple weirdos, 100m ones
files <- data.frame(path = list.files(file.path(root, "PredictionRasters"), full.names = TRUE, recursive = TRUE, pattern="*.tif"),
                    file = list.files(file.path(root, "PredictionRasters"), full.names = FALSE, recursive = TRUE, pattern="*.tif")) %>% 
  dplyr::filter(!str_sub(file, 1, 7)=="Archive") %>% 
  rowwise() %>% 
  mutate(int = str_locate_all(file, "/"),
         int = max(int)) %>% 
  ungroup() %>% 
  mutate(file = str_sub(file, int+1, 100)) %>% 
  separate(file, into=c("var", "scale", "year", "filetype"), remove=FALSE) %>% 
  mutate(year = as.numeric(ifelse(is.na(filetype), NA, year))) %>% 
  dplyr::filter(!var %in% c("Standardized", "VLCE100"),
                scale!="100m") %>% 
  mutate(cov = paste0(var, "_", scale))

#2. Add the two that are duplicates with temporal mistmatch---
add <- files %>% 
  dplyr::filter(cov %in% c("ERAPPTsm_1km", "ERATavesm_1km")) %>% 
  mutate(cov = gsub(pattern="sm", replacement="smt", x=cov))

files.use <- rbind(files, add) %>% 
  arrange(cov, year)

#3. Check against the covariate list----
missing <- data.frame(list = names(covlist)) %>% 
  dplyr::filter(!list %in% c(files.use$cov, "bcr"))
nrow(missing)

#4. Get the t-1 ones----
meth <- readxl::read_excel(file.path(root, "NationalModels_V5_VariableList.xlsx"), sheet = "ExtractionLookup")

lag <- meth %>% 
  dplyr::filter(YearMatch==-1)

#5. Match to each desired year of prediction----

#Get the covariate list
cov <- unique(files.use$cov)

#Set the years for prediction
years <- seq(1985, 2020, 5)

year.out <- data.frame()
for(i in 1:length(cov)){
  
  #Get the available years for that covariate
  cov.i <- files.use %>% 
    dplyr::filter(cov==cov[[i]])
  
  #Identify if is a cov that doesn't have annual layers---
  if(is.na(cov.i$year)[1]){
    
    year.out <- rbind(year.out,
                      data.frame(cov = cov[i],
                                 year = NA,
                                 predyear = years))
  } else {
    
    #Adjust the time lag as needed
    if(cov[i] %in% lag$Label){ cov.i$year <- cov.i$year + 1}
    
    #Match to nearest year of available data
    dt = data.table::data.table(year=unique(cov.i$year), val = unique(cov.i$year))
    data.table::setattr(dt, "sorted", "year")
    data.table::setkey(dt, year)
    year.i <- dt[J(years), roll = "nearest"]$val
    
    #Put it together----
    year.out <- rbind(year.out,
                      data.frame(cov = cov[i],
                                 year = year.i,
                                 predyear = years))
    
  }
  
}

#6. Put together with file paths----
files.year <- year.out %>% 
  left_join(files.use,
            multiple="all") %>% 
  unique()
  
#MAKE STACKS####

#1. Set up loop for subunits * year----
units <- bcr.out %>% 
  st_drop_geometry() %>% 
  mutate(bcr = paste0(country, subUnit)) %>% 
  unique() %>% 
  expand_grid(year = seq(1985, 2020, 5))

for(i in 1:nrow(units)){
  
  #2. Get subunit----
  bcr.i <- paste0(units$bcr[i])
  year.i <- paste0(units$year[i])
  
  #3. Get covariate list----
  covlist.i <- covlist %>% 
    dplyr::filter(bcr==bcr.i) %>% 
    pivot_longer(ERAMAP_1km:mTPI_1km, names_to="cov", values_to="use") %>% 
    dplyr::filter(use==TRUE)
  
  #4. Get file paths----
  files.i <- files %>% 
    dplyr::filter(cov %in% covlist.i$cov,
                  year.rd == year.i)
  
  #5. Read them in----
  rast(files.i$path)
  
  #6. Get the bcr shp----
  shp.i <- bcr.out %>% 
    dplyr::filter(country==units$country[i],
                  subUnit==units$subUnit[i])
  
  #7. Crop----
  
  #8. Sanity check plot----
  
  #9. Save----
  
  
}