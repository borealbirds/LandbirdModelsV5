# ---
# title: National Models 5.0 - extract covariates
# author: Elly Knight
# created: December 22, 2022
# ---

#NOTES################################

#This script uses CONUS Albers equal-area conic (EPSG: 5072)

#This script extracts the list of model covariates agreed upon by the BAM National Model team. The full list of those covariates can be found at https://docs.google.com/spreadsheets/d/1XATuq8BOYC2KkbJObaturaf4NFD2ufvn/edit?usp=sharing&ouid=104837701987164094932&rtpof=true&sd=true

#Records that return NA values from the layers stored in google drive are locations that are outside the area of the raster (i.e., coastlines)

#This script was written in pieces while assembling the various covariates and is unlikely inefficient. It should be revised for the next iteration of the national models. Potential improvements for next version include:

#This script extracts covariates for every combination of location & year. Efficiency could be improved by ~15% by taking only the unique locations for all rows in the method table with TemporalResolution=="static"

#This script extracts covariates for all locations. Efficiency could be improved by filtering to layer extent prior to extraction.

#This script extracts covariates for all data points. Efficiency could be improved by removing surveys outside of acceptable survey windows and the survey area prior to extraction

#NOTES FOR V6:

#1. One way to improve the reproducibility and adding new layers would be to compile the list of layers to extract by comparing the full list to the column names of existing compiled objects instead of using the "running" and "complete" lookup columns. Similarly for GEE extraction and the downloaded files.

#2. Oh right, and should probably parallelize it all.

#3. Should swap out reduceregions rgee code for ee_extract in next version for faster implementation

#4. Figure out why ALAN has negative values in extraction

#5. Add sanity check: covariate distribution against prediction layers. BEtter sanity checks for categorical analyses that don't rely on extrapolation.

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

#2. Create sf object of just location*year for annual layers----
loc.yr <- visit |> 
  dplyr::select(project, location, lat, lon, year) |> 
  unique() |> 
  mutate(id=paste(project, location, lat, lon, year, sep="_")) |> 
  st_as_sf(coords=c("lon", "lat"), crs=4326, remove=FALSE) |> 
  st_transform(crs=5072)

#EXTRA. Test sample dataset----
# set.seed(1234)
# loc.n <- loc.yr |>
#   sample_n(100000)
loc.n <- loc.yr

#3. Buffer location objects----
#200 m radius for local extent
#2 km radius to landscape extent
loc.buff <- st_buffer(loc.n, 200)
loc.buff2 <- st_buffer(loc.n, 2000)

#B. EXTRACT COVARIATES FROM GOOGLE DRIVE####

#1. Get list of layers to run----
meth.gd <- dplyr::filter(meth, Source=="Google Drive")

#2. Plain dataframe for joining to output----
 # loc.gd <- data.frame(loc.n) |>
 #   dplyr::select(-geometry)

#Read in existing dataframe if not starting from scratch
loc.gd <- read.csv(file=file.path(root, "Data", "Covariates", "03_NM5.0_data_covariates_GD.csv"))

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
      loc.cov <- loc.n |> 
        st_transform(crs(rast.i)) |> 
        terra::extract(x=rast.i, ID=FALSE)
      
      loc.cov$id <- loc.n$id
    }
    
    if(meth.gd.i$Extraction[1]=="radius"){
      
      #Extract point value, filter out NAs, reproject to UTM, buffer, reproject to raster CRS, extract value, rejoin to full df
      
      loc.cov.i <- loc.n |> 
        st_transform(crs(rast.i)) |> 
        terra::extract(x=rast.i, bind=TRUE) |> 
        st_as_sf() |> 
        rename(cov = names(rast.i)[1]) |> 
        dplyr::filter(!is.na(cov)) |> 
        dplyr::select(-cov) |> 
        st_transform(crs=5072) |> 
        st_buffer(meth.gd.i$RadiusExtent[1]) |> 
        st_transform(crs(rast.i))
      
      loc.cov <- cbind(loc.cov.i,
                       exact_extract(x=rast.i,
                                     y=loc.cov.i,
                                     meth.gd.i$RadiusFunction[1],
                                     force_df=TRUE)) |> 
        st_drop_geometry() |> 
        right_join(loc.n |> 
                     st_drop_geometry())
      
    }
    
  }
  
  #8. Determine if temporally matching----
  if(meth.gd.i$TemporalResolution[1]=="match"){
    
    #9. Get list of individual files----
    
    #For Landfire file structure
    if(any(str_detect(meth.gd.i$Link,"LandFire"))){
      files.i <- data.frame(Link=list.files(meth.gd.i$Link, full.names = TRUE, recursive=TRUE, pattern="*.tif"),
                            file=list.files(meth.gd.i$Link, recursive=TRUE, pattern="*.tif")) |> 
        separate(file, into=c("year", "region", "name", "tif"), remove=FALSE) |> 
        mutate(name = case_when(name %in% c("105cbd", "130cbd", "140cbd", "220cbd") ~ "biomass",
                                name %in% c("105cc", "130cc", "140cc", "220cc") ~ "closure",
                                name %in% c("105ch", "130ch", "140ch", "220ch") ~ "height"),
               year = as.numeric(year)) |> 
        arrange(year) |> 
        unique()
      
      #Just height for CV extraction
      if(any(str_detect(meth.gd.i$Label,"cv"))){
        files.i <- files.i |> 
          dplyr::filter(name=="height")
      }
      
    }
    
    #For NLCD file structure
    else if(str_detect(meth.gd.i$Link[1],"NLCD")){
      files.i <- data.frame(Link=list.files(meth.gd.i$Link, full.names = TRUE, pattern="*.img"),
                            file=list.files(meth.gd.i$Link, pattern="*.img")) |> 
        mutate(year = as.numeric(str_match(file, "cd_\\s*(.*?)\\s*_lan")[,2])) |> 
        arrange(year) |> 
        unique() |> 
        dplyr::filter(!is.na(year))
    }
    
    #Everything else
    else {
      files.i <- data.frame(Link=list.files(meth.gd.i$Link, full.names = TRUE, pattern="*tif", recursive=TRUE),
                            file=list.files(meth.gd.i$Link, pattern="*tif", recursive=TRUE)) |> 
        mutate(year = as.numeric(str_extract(file, "[^_]+(?=\\.tif$)"))) |> 
        arrange(year) |> 
        unique() |> 
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
      files.j <- dplyr::filter(files.i, year==yrs[j], name %in% meth.gd.i$RasterName)
      rast.i <- rast(files.j$Link)
      names(rast.i) <- meth.gd.i$Label
      
      #13. NA out -9999s for Landfire----
      if(str_detect(files.j$Link[1],"LandFire")){ rast.i <- subst(x=rast.i, from=-9999, to=NA)}
      
      #14. Extract - determine point or buffer extraction and buffer extent---
      if(meth.gd.i$Extraction[1]=="point"){
        loc.cov <- loc.n.j |> 
          st_transform(crs(rast.i)) |> 
          terra::extract(x=rast.i, ID=FALSE) |> 
          data.table::setnames(meth.gd.i$Label) |> 
          cbind(loc.n.j |> 
                  st_drop_geometry() |> 
                  dplyr::select(id)) |> 
          rbind(loc.cov)
      }
      
      if(meth.gd.i$Extraction[1]=="radius" & meth.gd.i$RadiusExtent[1]==200){
        loc.cov <- loc.buff.j |> 
          st_transform(crs(rast.i)) |> 
          exact_extract(x=rast.i, meth.gd.i$RadiusFunction, force_df=TRUE) |> 
          data.table::setnames(meth.gd.i$Label) |> 
          cbind(loc.n.j |> 
                  st_drop_geometry() |> 
                  dplyr::select(id)) |> 
          rbind(loc.cov)
      }
      
      if(meth.gd.i$Extraction[1]=="radius" & meth.gd.i$RadiusExtent[1]==2000){
        loc.cov <- loc.buff2.j |> 
          st_transform(crs(rast.i)) |> 
          exact_extract(x=rast.i, meth.gd.i$RadiusFunction, force_df=TRUE) |> 
          data.table::setnames(meth.gd.i$Label) |> 
          cbind(loc.n.j |> 
                  st_drop_geometry() |> 
                  dplyr::select(id)) |> 
          rbind(loc.cov) 
      }
      
      print(paste0("Finished year ", j, " of ", length(yrs)))
      
    }
    
  }
  
  #15. Fix column names----
  if(any(str_detect(meth.gd.i$Link,"CONUS"))){
    colnames(loc.cov) <- c(paste0(meth.gd.i$Label, "_conus"), "id")
    nms <- c(colnames(loc.gd), paste0(meth.gd.i$Label, "_conus"))
  } else {nms <- c(colnames(loc.gd), meth.gd.i$Label) }
  
  #16. Add output to main file----
  loc.gd <- left_join(loc.gd, loc.cov)
  colnames(loc.gd) <- nms
  
  #17. Save----
  write.csv(loc.gd, file=file.path(root, "Data", "Covariates", "03_NM5.0_data_covariates_GD.csv"), row.names=FALSE)
  
  #18. Report status----
  print(paste0("Finished stack category ", i, " of ", length(loop)))
  
}

#19. Merge AK & CONUS columns for landfire----

if("LFheigthcv_1km_conus" %in% colnames(loc.gd)){
  loc.gd2 <- loc.gd  |>
    mutate(LFbiomass_1km = ifelse(!is.na(LFbiomass_1km), LFbiomass_1km, LFbiomass_1k_conus),
           LFcrownclosure_1km = ifelse(!is.na(LFcrownclosure_1km), LFcrownclosure_1km, LFcrownclosure_1km_conus),
           LFheigth_1km = ifelse(!is.na(LFheigth_1km), LFheigth_1km, LFheigth_1km_conus),
           LFheigthcv_1km = ifelse(!is.na(LFheigthcv_1km), LFheigthcv_1km, LFheigthcv_1km_conus),
           LFbiomass_5x5 = ifelse(!is.na(LFbiomass_5x5), LFbiomass_5x5, LFbiomass_5x5_conus),
           LFcrownclosure_5x5 = ifelse(!is.na(LFcrownclosure_5x5), LFcrownclosure_5x5, LFcrownclosure_5x5_conus),
           LFheigth_5x5 = ifelse(!is.na(LFheigth_5x5), LFheigth_5x5, LFheigth_5x5_conus),
           LFheigthcv_5x5 = ifelse(!is.na(LFheigthcv_5x5), LFheigthcv_5x5, LFheigthcv_5x5_conus)) |>
    dplyr::select(-LFbiomass_1km_conus, -LFcrownclosure_1km_conus, -LFheigth_1km_conus, -LFheigthcv_1km_conus,
                  -LFbiomass_5x5_conus, -LFcrownclosure_5x5_conus, -LFheigth_5x5_conus, -LFheigthcv_5x5_conus)
  
} else { loc.gd2 <- loc.gd }

write.csv(loc.gd2, file=file.path(root, "Data", "Covariates", "03_NM5.0_data_covariates_GD.csv"), row.names=FALSE)

#C. EXTRACT COVARIATES FROM SCANFI####

#1. Get list of covariates to run----
meth.scanfi <- dplyr::filter(meth, Source=="SCANFI", Running==1)

#2. Get & wrangle list of SCANFI layers----
#Note: landcover layer has a different name in 1985 (VegTypeClass vs nfiLandCover)
files.scanfi <- data.frame(Link=list.files("G:/.shortcut-targets-by-id/11nj6IZyUe3EqrOEVEmDfrnkDCEQyinSi/SCANFI_share", full.names = TRUE,recursive=TRUE, pattern="*.tif"),
                           file=list.files("G:/.shortcut-targets-by-id/11nj6IZyUe3EqrOEVEmDfrnkDCEQyinSi/SCANFI_share", recursive=TRUE, pattern="*.tif")) |> 
  dplyr::filter(str_detect(file, "tif")) |> 
  mutate(file = str_sub(file, 13, 100)) |> 
  separate(file, into=c("scanfi", "covtype", "variable", "S", "year", "v0", "filetype")) |> 
  mutate(year = as.numeric(ifelse(filetype=="v0", v0, year)),
         RasterName = case_when(variable=="VegTypeClass" ~ "scanfilcc",
                          variable=="prcC" ~ "conifer",
                          variable=="prcD" ~ "deciduous",
                          variable=="prcB" ~ "deciduous",
                          variable=="LodepolePine" ~ "lodgepolepine",
                          variable=="nfiLandCover" ~ "scanfilcc",
                          !is.na(variable) ~ tolower(variable))) |> 
  dplyr::filter(scanfi%in%c("SCANFI", "CaNFIR"),
                RasterName %in% meth.scanfi$RasterName) |> 
  dplyr::select(Link, RasterName, year) |> 
  arrange(year, RasterName)
#Check everything is there
table(files.scanfi$RasterName, files.scanfi$year)

#3. Remove points outside of SCANFI coverage----
#do this for SCANFI but not for others because scanfi is so much more time intensive (many high resolution layers)
rast.scanfi <- rast(files.scanfi$Link[1])

loc.scanfi.sa <- loc.n |> 
  st_transform(crs(rast.scanfi)) |> 
  terra::extract(x=rast.scanfi, ID=FALSE)
colnames(loc.scanfi.sa) <- "scanfi"

#4. Rebuffer----
loc.scanfi.buff <- cbind(loc.n, loc.scanfi.sa) |> 
  dplyr::filter(!is.na(scanfi)) |> 
  st_transform(5072) |> 
  st_buffer(200) |> 
  st_transform(crs(rast.scanfi))

loc.scanfi.buff2 <- cbind(loc.n, loc.scanfi.sa) |> 
  dplyr::filter(!is.na(scanfi)) |> 
  st_transform(5072) |> 
  st_buffer(2000) |> 
  st_transform(crs(rast.scanfi))

#5. Match year of data to year of SCANFI----
#Use only buffer object because all extractions are radius method
years.scanfi <- unique(files.scanfi$year)
dt = data.table::data.table(year=years.scanfi, val=years.scanfi)
data.table::setattr(dt, "sorted", "year")
data.table::setkey(dt, year)
loc.scanfi.buff$year.rd <- dt[J(loc.scanfi.buff$year), roll = "nearest"]$val
loc.scanfi.buff2$year.rd <- dt[J(loc.scanfi.buff2$year), roll = "nearest"]$val

#6. Read in/create dataframe----
loc.scanfi <- data.frame()
#loc.scanfi <- read.csv(file.path(root, "Data", "Covariates", "03_NM5.0_data_covariates_SCANFI.csv"))

#7. Set up to loop through years of SCANFI----
years.scanfi <- unique(files.scanfi$year)
for(i in 2:length(years.scanfi)){
  
  loc.buff.yr <- dplyr::filter(loc.scanfi.buff, year.rd==years.scanfi[i]) |> 
    arrange(lat, lon) |> 
    mutate(loop = ceiling(row_number()/1000))
  
  loc.buff2.yr <- dplyr::filter(loc.scanfi.buff2, year.rd==years.scanfi[i]) |> 
    arrange(lat, lon) |> 
    mutate(loop = ceiling(row_number()/1000))
  
  #7. Read in rasters for 200m mean extraction----
  files.i <- meth.scanfi |>
    dplyr::select(-Link) |>
    dplyr::filter(RadiusFunction=="mean", RadiusExtent==200) |>
    left_join(files.scanfi, multiple="all") |>
    dplyr::filter(year==years.scanfi[i])
  rast.i <- rast(files.i$Link)
  names(rast.i) <- files.i$Label

  #8. Read in rasters for 200m cv extraction----
  filesd.i <- meth.scanfi |>
    dplyr::select(-Link) |>
    dplyr::filter(RadiusFunction=="cv", RadiusExtent==200) |>
    left_join(files.scanfi, multiple="all") |>
    dplyr::filter(year==years.scanfi[i])
  rastsd.i <- rast(filesd.i$Link)
  names(rastsd.i) <- filesd.i$Label
  
  #9. Read in rasters for 200m mode extraction----
  filesmd.i <- meth.scanfi |>
    dplyr::select(-Link) |>
    dplyr::filter(RadiusFunction=="mode", RadiusExtent==200) |>
    left_join(files.scanfi, multiple="all") |>
    dplyr::filter(year==years.scanfi[i])
  rastmd.i <- rast(filesmd.i$Link)
  names(rastmd.i) <- filesmd.i$Label
  
  #7. Read in rasters for 1.4km mean extraction----
  files2.i <- meth.scanfi |>
    dplyr::select(-Link) |>
    dplyr::filter(RadiusFunction=="mean", RadiusExtent==2000) |>
    left_join(files.scanfi, multiple="all") |>
    dplyr::filter(year==years.scanfi[i])
  rast2.i <- rast(files2.i$Link)
  names(rast2.i) <- files2.i$Label

  #8. Read in rasters for 1.4km cv extraction----
  filesd2.i <- meth.scanfi |>
    dplyr::select(-Link) |>
    dplyr::filter(RadiusFunction=="cv", RadiusExtent==2000) |>
    left_join(files.scanfi, multiple="all") |>
    dplyr::filter(year==years.scanfi[i])
  rastsd2.i <- rast(filesd2.i$Link)
  names(rastsd2.i) <- filesd2.i$Label
  
  #8. Set up loops----
  loc.scanfi.list <- list()
  for(j in 1:nrow(loc.buff.yr)){
    
    loc.i <- loc.buff.yr[j,]
    loc2.i <- loc.buff2.yr[j,]
    
    #9. Extract----
    loc.scanfi.i <- try(exact_extract(x=rast.i, y=loc.i, "mean", force_df=TRUE))
    loc.scanfi.sd.i <- try(exact_extract(x=rastsd.i, y=loc.i, fun="coefficient_of_variation", force_df=TRUE))
    loc.scanfi.md.i <- try(exact_extract(x=rastmd.i, y=loc.i, fun="mode", force_df=TRUE))
    loc.scanfi2.i <- try(exact_extract(x=rast2.i, y=loc2.i, "mean", force_df=TRUE))
    loc.scanfi.sd2.i <- try(exact_extract(x=rastsd2.i, y=loc2.i, fun="coefficient_of_variation", force_df=TRUE))
    colnames(loc.scanfi.sd2.i) <- "coefficient_of_variation_ls"
    
    #10. Save output if extraction works----
    if(class(loc.scanfi.md.i)=="data.frame"){
      loc.scanfi.list[[j]] <- cbind(loc.i, loc.scanfi.md.i
#                                    loc.scanfi.sd.i, loc.scanfi.md.i, loc.scanfi2.i, loc.scanfi.sd2.i
                                    ) |> 
        data.frame() |> 
        dplyr::select(-geometry, -year.rd, -loop, -scanfi)
    }
    
    #11. Report progress----
    print(paste0("Finished loop ", j, " of ", nrow(loc.buff.yr), " for year ", i, " of ", length(years.scanfi), ": ", years.scanfi[i]))
    
    flush.console()
    
  }
  
  loc.scanfi.bind <- data.table::rbindlist(loc.scanfi.list)
  names(loc.scanfi.bind) <- c(colnames(loc.i)[1:6], names(rast.i), names(rastsd.i), "landcover.2", names(rast2.i), names(rastsd2.i))
  
  loc.scanfi <- rbind(loc.scanfi, loc.scanfi.bind)
  
  #13. Save----
  write.csv(loc.scanfi, file=file.path(root, "Data", "Covariates", "03_NM5.0_data_covariates_SCANFI.csv"), row.names = FALSE)
  
}

#D. EXTRACT COVARIATES FROM GEE - TEMPORALLY STATIC####

#1. Create GCS bucket----
#This only needs to be done once!
# project_id <- ee_get_earthengine_path() |>
#   list.files(., "\\.json$", full.names = TRUE) |>
#   jsonlite::read_json() |>
#   '$'(project_id) # Get the Project ID
# 
# googleCloudStorageR::gcs_create_bucket("national_models", projectId = project_id)

#2. Get list of static layers to run----
meth.gee <- dplyr::filter(meth, Source=="Google Earth Engine", Running==1, TemporalResolution=="static")

#3. Make/get dataframe----
#loc.gee.static <- data.frame()
loc.gee.static <- read.csv(file.path(root, "Data", "Covariates", "03_NM5.0_data_covariates_GEE-static.csv"))

#3. Make method loop----
loop <- meth.gee |> 
  dplyr::select(RadiusFunction, RadiusExtent) |> 
  unique()

for(h in 1:nrow(loop)){
  
  meth.gee.h <- meth.gee |> 
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
    loc.gee <- loc.buff |> 
      mutate(loop = ceiling(row_number()/1000))
  }
  
  if(loop$RadiusExtent[h]==2000){
    loc.gee <- loc.buff2 |> 
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
                                         scale=meth.gee$GEEScale[i])
    }
    
    if(loop$RadiusFunction[h]=="cv"){
      img.red <- img.stack$reduceRegions(reducer=ee$Reducer$stdDev(),
                                         collection=poly,
                                         scale=meth.gee$GEEScale[i])
    }
    
    task.list[[i]] <- ee_table_to_gcs(collection=img.red,
                                      bucket="national_models",
                                      fileFormat = "CSV")
    task.list[[i]]$start()
    
    #10. Monitor----
    ee_monitoring(task.list[[i]], max_attempts=1000)
    
    #11.Download to local----
    ee_gcs_to_local(task = task.list[[i]], dsn=file.path(root, "Data", "Covariates", "GEE", "Static", paste0("03_NM5.0_data_covariates_GEE_static_", loop$RadiusFunction[h], "_", loop$RadiusExtent[h], "_", i, ".csv")))
    
    print(paste0("Finished batch ", i, " of ", max(loc.gee$loop)))
    
  }
  
  print(paste0("Finished loop ", h, " of ", nrow(loop)))
  
}

#12. Read in output files----
files.gee <- data.frame(path = list.files(file.path(root, "Data", "Covariates", "GEE", "Static"), full.names = TRUE)) |>
  separate(path, into=c("j1", "j2", "j3", "j4", "j5", "j6", "j7", "method", "extent", "n"), sep="_", remove=FALSE)

files.gee.cv.200 <- dplyr::filter(files.gee, method=="cv", extent=="200")
loc.gee.cv.200 <- purrr::map(files.gee.cv.200$path, read.csv) |> 
  data.table::rbindlist()

files.gee.cv.2000 <- dplyr::filter(files.gee, method=="cv", extent=="2000")
loc.gee.cv.2000 <- purrr::map(files.gee.cv.2000$path, read.csv) |> 
  data.table::rbindlist()

files.gee.mean.200 <- dplyr::filter(files.gee, method=="mean", extent=="200")
loc.gee.mean.200 <- purrr::map(files.gee.mean.200$path, read.csv) |> 
  data.table::rbindlist()

files.gee.mean.2000 <- dplyr::filter(files.gee, method=="mean", extent=="2000")
loc.gee.mean.2000 <- purrr::map(files.gee.mean.2000$path, read.csv) |> 
  data.table::rbindlist()

#13. Fix column names----
lab.cv.200 <- meth.gee |> 
  dplyr::filter(RadiusFunction=="cv",
                RadiusExtent=="200")

lab.cv.2000 <- meth.gee |> 
  dplyr::filter(RadiusFunction=="cv",
                RadiusExtent=="2000")

lab.mean.200 <- meth.gee |> 
  dplyr::filter(RadiusFunction=="mean",
                RadiusExtent=="200")

lab.mean.2000 <- meth.gee |> 
  dplyr::filter(RadiusFunction=="mean",
                RadiusExtent=="2000")

#14. Fix column names, put together, calculate----
loc.gee <- loc.gee.cv.200 |> 
  data.table::setnames(c(colnames(loc.gee.cv.200)[1:7],
                         lab.cv.200$Label,
                         colnames(loc.gee.cv.200)[8:9])) |> 
  dplyr::select(c(id, lab.cv.200$Label)) |> 
  full_join(loc.gee.cv.2000 |> 
              data.table::setnames(c(colnames(loc.gee.cv.2000)[1:7],
                                     lab.cv.2000$Label,
                                     colnames(loc.gee.cv.2000)[8:9])) |> 
              dplyr::select(c(id, lab.cv.2000$Label))) |> 
  full_join(loc.gee.mean.200 |> 
              data.table::setnames(c(colnames(loc.gee.mean.200)[1],
                                     lab.mean.200$Label[1:2],
                                     colnames(loc.gee.mean.200)[4:8],
                                     lab.mean.200$Label[3],
                                     colnames(loc.gee.mean.200)[10],
                                     lab.mean.200$Label[4:5],
                                     colnames(loc.gee.mean.200)[13:14])) |> 
              dplyr::select(c(id, lab.mean.200$Label))) |> 
  full_join(loc.gee.mean.2000 |> 
              data.table::setnames(c(colnames(loc.gee.mean.2000)[1],
                                     lab.mean.2000$Label[1:2],
                                     colnames(loc.gee.mean.2000)[4:8],
                                     lab.mean.2000$Label[3],
                                     colnames(loc.gee.mean.2000)[10:12])) |> 
              dplyr::select(c(id, lab.mean.2000$Label))) |> 
  mutate(ETHheightcv_1km = ETHheightcv_1km/ETHheight_1km,
         ETHheightcv_5x5 = ETHheightcv_5x5/ETHheight_5x5) |> 
  mutate(ETHheightcv_1km = ifelse(is.infinite(ETHheightcv_1km), 0, ETHheightcv_1km),
         ETHheightcv_5x5 = ifelse(is.infinite(ETHheightcv_5x5), 0, ETHheightcv_5x5))

#15. Zerofill----
zerocols <- meth.gee |> 
  dplyr::filter(Zerofill==1)
zero.gee <- loc.gee |> 
  dplyr::select(zerocols$Label) |> 
  replace_na(list(occurrence=0, recurrence=0, seasonality=0, occurrence_ls=0))
loc.gee.static <- loc.gee |> 
  dplyr::select(-zerocols$Label) |> 
  cbind(zero.gee)

#16. Save----
write.csv(loc.gee.static, file=file.path(root, "Data", "Covariates", "03_NM5.0_data_covariates_GEE-static.csv"), row.names = FALSE)

#E. EXTRACT COVARIATES FROM GEE - TEMPORALLY MATCHED####

#1. Get list of static layers to run----
meth.gee <- dplyr::filter(meth, Source=="Google Earth Engine", TemporalResolution=="match")

#2. Plain dataframe for joining to output----
#loc.gee <- data.frame(loc.n) |> 
#   dplyr::select(-geometry)
loc.gee <- read.csv(file=file.path(root, "Data", "Covariates", "03_NM5.0_data_covariates_GEE-match.csv"))

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
    
    loc.n.yr <- dplyr::filter(loc.n.i, yearrd==years.gee[j]) |> 
      mutate(loop = ceiling(row_number()/1000))
    loc.buff2.yr <- dplyr::filter(loc.buff2.i, yearrd==years.gee[j]) |> 
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
          loc.k[[k]] <- ee_extract(x=img.i, y=point, scale=meth$GEEScale[i])
        }
        
        if(meth.gee$Extraction[i]=="radius"){
          
          if(meth.gee$RadiusFunction[i]=="mean"){
            img.red <- img.i$reduceRegions(reducer=ee$Reducer$mean(),
                                           collection=poly,
                                           scale=meth.gee$GEEScale[i])
          }
          
          if(meth.gee$RadiusFunction[i]=="mode"){
            img.red <- img.i$reduceRegions(reducer=ee$Reducer$mode(),
                                           collection=poly, 
                                           scale=meth.gee$GEEScale[i])
          }
          
          task.list[[k]] <- ee_table_to_gcs(collection=img.red,
                                            bucket="national_models",
                                            fileFormat = "CSV")
          task.list[[k]]$start()
          
          ee_monitoring(task.list[[k]], max_attempts=1000)
          
          ee_gcs_to_local(task = task.list[[k]], dsn=file.path(root, "Data", "Covariates", "GEE", "Match", paste0("03_NM5.0_data_covariates_GEE_match_", meth.gee$Label[i], "_", years.gee[j], "_", k, ".csv")))
          
          loc.k[[k]] <- read.csv(file.path(root, "Data", "Covariates", "GEE", "Match", paste0("03_NM5.0_data_covariates_GEE_match_", meth.gee$Label[i], "_", years.gee[j], "_", k, ".csv")))
          
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
        loc.j[[j]] <- data.table::rbindlist(loc.k, fill=TRUE) |> 
          dplyr::select(id, mode)
        
        #13. Fix column names----
        colnames(loc.j[[j]]) <- c("id", meth.gee$Label[i])
      }
      
      
    }
    
    print(paste0("Finished year ", j, " of ", length(years.gee)))
  
  #14. Collapse data across years----
  loc.i <- data.table::rbindlist(loc.j, fill=TRUE)
  loc.gee <- cbind(loc.gee, loc.i |> 
                     dplyr::select(meth.gee$Label[i]))
  
  write.csv(loc.gee, file=file.path(root, "Data", "Covariates", "03_NM5.0_data_covariates_GEE-match.csv"), row.names=FALSE)
  
  print(paste0("FINISHED LAYER ", i, " of ", nrow(meth.gee)))
  
}

#15. Read in output files for radius extraction----
files.gee <- data.frame(path = list.files(file.path(root, "Data", "Covariates", "GEE", "Match"), full.names = TRUE)) |>
  separate(path, into=c("j1", "j2", "j3", "j4", "j5", "j6", "j7", "covariate", "extent", "year"), sep="_", remove=FALSE)

vars <- files.gee |> 
  dplyr::select(extent, covariate) |> 
  unique()

for(i in 1:nrow(vars)){
  
  files.gee.i <- dplyr::filter(files.gee, covariate==vars$covariate[i],
                               extent==vars$extent[i])
  loc.gee.i <- purrr::map(files.gee.i$path, read.csv) |> 
    data.table::rbindlist()
  
  label <- paste0(vars$covariate[i], "_", vars$extent[i])
  
  loc.gee <- loc.gee.i |> 
    data.table::setnames(c(colnames(loc.gee.i)[1:6],
                           label,
                           colnames(loc.gee.i)[8:11])) |> 
    dplyr::select(c(id, all_of(label))) |> 
    full_join(loc.gee)
  
}

#16. Save again----
write.csv(loc.gee, file=file.path(root, "Data", "Covariates", "03_NM5.0_data_covariates_GEE-match.csv"), row.names=FALSE)

#E. ASSEMBLE####

#1. Load data----
load(file.path(root, "Data", "02_NM5.0_data_offsets.R"))

#2. Load extracted covariates----
loc.gee.match <- read.csv(file.path(root, "Data", "Covariates", "03_NM5.0_data_covariates_GEE-match.csv"))

loc.gee.static <- read.csv(file.path(root, "Data", "Covariates", "03_NM5.0_data_covariates_GEE-static.csv"))

loc.gd <- read.csv(file=file.path(root, "Data", "Covariates", "03_NM5.0_data_covariates_GD.csv"))

loc.scanfi <- read.csv(file.path(root, "Data", "Covariates", "03_NM5.0_data_covariates_SCANFI.csv"))

#3. Attribute to country----
can <- read_sf(file.path(root, "Regions", "CAN_adm", "CAN_adm0.shp")) |> 
  dplyr::select(geometry)

visit.country <- visit |> 
  st_as_sf(coords=c("lon", "lat"), crs=4326, remove=FALSE) |> 
  st_intersection(can) |> 
  mutate(country = "can") |> 
  full_join(visit) |> 
  mutate(country = ifelse(is.na(country), "usa", country)) |> 
  st_drop_geometry()

#3. Add id to visit and join together and wrangle----
meth.use <- meth |> 
  dplyr::filter(Use==1)

#Remove lat lon fields due to rounding errors that cause mismatches
#Remove 0 values for some layers that should be NAs
#format landcover classes as factors
#create eBird method and make method a factor
#zero out heights < 0.1 and NA the height cv values for those
#convert SCANFI species variables to biomass and remove biomass

scanfi_1km <- c("SCANFIBalsamFir_1km", "SCANFIBlackSpruce_1km", "SCANFIprcC_1km", "SCANFIprcD_1km", "SCANFIDouglasFir_1km", "SCANFIJackPine_1km", "SCANFILodgepolePine_1km", "SCANFIPonderosaPine_1km", "SCANFITamarack_1km", "SCANFIWhiteRedPine_1km")

scanfi_5x5 <- c("SCANFIBalsamFir_5x5", "SCANFIBlackSpruce_5x5", "SCANFIprcC_5x5", "SCANFIprcD_5x5", "SCANFIDouglasFir_5x5", "SCANFIJackPine_5x5", "SCANFILodgepolePine_5x5", "SCANFIPonderosaPine_5x5", "SCANFITamarack_5x5", "SCANFIWhiteRedPine_5x5")

visit.covs <- visit.country |> 
  rename(row = id) |> 
  mutate(id=paste(project, location, lat, lon, year, sep="_")) |> 
  left_join(loc.gd |> 
              dplyr::select(-lat, -lon)) |> 
  left_join(loc.scanfi |> 
              dplyr::select(-lat, -lon)) |> 
  left_join(loc.gee.static) |> 
  left_join(loc.gee.match) |> 
  dplyr::select(-id) |> 
  rename(id = row) |> 
  dplyr::select(all_of(colnames(visit.country)), all_of(meth.use$Label)) |> 
  mutate(NLCD_1km = ifelse(NLCD_1km==0, NA, NLCD_1km),
         VLCE_1km = ifelse(VLCE_1km==0, NA, VLCE_1km),
         usroad_1km = ifelse(country=="usa", usroad_1km, NA),
         usroad_5x5 = ifelse(country=="usa", usroad_5x5, NA),
         canroad_1km = ifelse(country=="can", canroad_1km, NA),
         canroad_5x5 = ifelse(country=="can", canroad_5x5, NA),) |> 
  mutate(MODISLCC_1km = factor(MODISLCC_1km),
         MODISLCC_5x5 = factor(MODISLCC_5x5),
         SCANFI_1km = factor(SCANFI_1km),
         NLCD_1km = factor(NLCD_1km),
         ABoVE_1km = factor(ABoVE_1km),
         VLCE_1km = factor(VLCE_1km),
         tagMethod = factor(tagMethod)) |> 
  mutate(method = ifelse(source=="eBird", "eBird", as.character(tagMethod)),
         method = factor(method, levels=c("PC", "eBird", "1SPM", "1SPT"))) |> 
  mutate(SCANFIheight_1km = ifelse(SCANFIheight_1km < 0.1, 0, SCANFIheight_1km),
         SCANFIheight_5x5 = ifelse(SCANFIheight_5x5 < 0.1, 0, SCANFIheight_5x5),
         SCANFIheightcv_1km = ifelse(SCANFIheight_1km < 0.1, NA, SCANFIheightcv_1km),
         SCANFIheightcv_5x5 = ifelse(SCANFIheight_5x5 < 0.1, NA, SCANFIheightcv_5x5),
         ETHheight_1km = ifelse(ETHheight_1km < 0.1, 0, ETHheight_1km),
         ETHheight_5x5 = ifelse(ETHheight_5x5 < 0.1, 0, ETHheight_5x5),
         ETHheightcv_1km = ifelse(ETHheight_1km < 0.1, NA, ETHheightcv_1km),
         ETHheightcv_5x5 = ifelse(ETHheight_5x5 < 0.1, NA, ETHheightcv_5x5),
         LFheigth_1km = ifelse(LFheigth_1km < 0.1, 0, LFheigth_1km),
         LFheigth_5x5 = ifelse(LFheigth_5x5 < 0.1, 0, LFheigth_5x5),
         LFheigthcv_1km = ifelse(LFheigth_1km < 0.1, NA, LFheigthcv_1km),
         LFheigthcv_5x5 = ifelse(LFheigth_5x5 < 0.1, NA, LFheigthcv_5x5)) |>
  mutate_at(all_of(scanfi_1km), ~(. * SCANFIbiomass_1km),
            all_of(scanfi_5x5), ~(. & SCANFIbiomass_5x5)) |> 
  dplyr::select(-SCANFIbiomass_1km, -SCANFIbiomass_5x5)

#5. Sanity checks----

#Take a sample of the visits to look at each cov
set.seed(1234)
visit.sample <- visit.covs |> 
  sample_n(10000)

#Remove biomass vars
meth.plot <- meth.use |> 
  dplyr::filter(!Label %in% c("SCANFIbiomass_1km", "SCANFIbiomass_5x5"))

#Plot to a folder
for(i in 1:nrow(meth.plot)){
  
  visit.i <- visit.sample |> 
    dplyr::select(lat, lon, meth.plot$Label[i]) |> 
    data.table::setnames(c("lat", "lon", "cov")) |> 
    mutate(na = ifelse(is.na(cov), "NA", "VALUE"))
  
  plot.na.i <- ggplot(visit.i) +
    geom_point(aes(x=lon, y=lat, colour=na))
  
  ggsave(plot.na.i, filename=file.path(root, "Data", "Covariates", "Plots", "NA", paste0(meth.plot$Label[i], ".jpeg")), width=10, height=8)
  
  plot.i <- ggplot(visit.i) +
    geom_point(aes(x=lon, y=lat, colour=cov))
  
  ggsave(plot.i, filename=file.path(root, "Data", "Covariates", "Plots", "Cov", paste0(meth.plot$Label[i], ".jpeg")), width=10, height=8)
  
  print(paste0("Finished plot ", i, " of ", nrow(meth.plot)))
  
}

#G. SAVE#####
visit <- visit.covs

save(visit, bird, offsets, file=file.path(root, "Data", "03_NM5.0_data_covariates.R"))
