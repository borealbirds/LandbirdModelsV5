# ---
# title: National Models 5.0 - stratify by BCR*country
# author: Elly Knight
# created: September 29, 2022
# ---

#NOTES################################

#In this script, we stratify & prepare the data for modelling. Steps include:

#1. BCR Attribution: Diving the data into BCRs for separate models. Each BCR is buffered by 100km. We do this so that we can feather predictions from adjacent regions together. The exception is the international boundaries between US and Canada, which we don't buffer because spatial layers for the covariates are different on either side of the border. In this case, we use a shapefiles of country boundaries to intersect with the buffered regions.

#2. Filtering the data. We remove data outside the study area, data before a minimum year, and data outside of dawn surveys. Additional dataset filtering steps can be added here as necessary.

#3. Aggregating subunits based on minimum sample size threshold. Sample sizes are very low in some subunits in Alaska and northwestern BC.

#4. Thinning. First we assign each survey to a grid cell (1.9 km spacing) Then we randomly selecting one survey per cell per year for each bootstrap. Note on line 186-192, we use an alternative to sample_n(), which does provide a truly random sample and biases towards lines of data near the top of the dataframe. This step takes a while to run.

#5. Indicating species that should be modelled for each subunit based on a minimum detection threshold.

#6. Determining which covariates should be used for each subunit. This is done using the same extraction lookup table as the '03.ExtractCovariates.R' script. For each subunit, we choose the version of a covariate with the highest priority that has full coverage for that subunit. We then rename the covariates so that they match the naming conventions of the prediction layers for the prediction step.

#Removal of non-dawn surveys is a secondary line of defense in case temporal patterns of QPAD offsets are unrealistic. In other words, we are restricting the dataset to times of day when p(availability) should be relatively high. This filter should be removed once QPAD offsets have been revisited. It should also be revisited if/when the national models include nocturnal species.

#TODO: BIRD OBJECT HAS NO PKEY

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(terra) #basic raster handling
library(sf) #basic shapefile handling
library(fasterize) #fast rasterization of shapefiles
library(exactextractr) #fast & efficient raster extraction
library(dggridR) #grid for spatial thinning
library(Matrix) #for sparse matrix conversion

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_NationalModels/NationalModels5.0"

#3. Load data packages with offsets and covariates----
load(file.path(root, "Data", "03_NM5.0_data_covariates.R"))

#4. Turn off scientific notation---
options(scipen=99999)

#BCR ATTRIBUTION#####################

#1. Make visit a vector object----
visit.v <- visit %>% 
  dplyr::select(id, lat, lon) %>% 
  st_as_sf(coords=c("lon", "lat"), crs=4326) %>% 
  st_transform(5072) %>% 
  vect()

#2. Read in country shapefiles----
can <- read_sf(file.path(root, "Regions", "CAN_adm", "CAN_adm0.shp")) %>% 
  st_transform(crs=5072) 
usa <- read_sf(file.path(root, "Regions", "USA_adm", "USA_adm0.shp")) %>% 
  st_transform(crs=5072) 

#3. Read in BCR shapefile----
#Remove subunit 1
bcr <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel.shp")) %>% 
  dplyr::filter(subUnit!=1) %>% 
  st_transform(crs=5072)

ggplot(bcr) +
  geom_sf(aes(fill=factor(subUnit)))

#4. Identify BCRs for each country----
bcr.ca <- bcr %>% 
  st_intersection(can) %>% 
  mutate(country="can")

bcr.usa <- bcr %>% 
  st_intersection(usa) %>% 
  mutate(country="usa")

bcr.country <- rbind(bcr.ca, bcr.usa)
  
#5. Set up loop for BCR attribution----
bcr.df <- data.frame(id=visit$id)
for(i in 1:nrow(bcr.country)){
  
  #6. Filter & buffer shapefile----
  bcr.buff <- bcr.country %>% 
    dplyr::filter(row_number()==i) %>% 
    st_buffer(100000)
  
  #7. Crop to international boundary----
  if(bcr.buff$country=="can"){ bcr.i <- st_intersection(bcr.buff, can)}
  if(bcr.buff$country=="usa"){ bcr.i <- st_intersection(bcr.buff, usa)}
  
  #8. Convert to raster for fast extraction----
  r <- rast(ext(bcr.i), resolution=1000, crs=crs(bcr.i))
  bcr.r <- rasterize(x=bcr.i, y=r, field="subUnit")
  
  #9. Extract raster value----
  bcr.out <- data.frame(subUnit=extract(x=bcr.r, y=visit.v)[,2]) %>% 
    mutate(use = ifelse(!is.na(subUnit), TRUE, FALSE))
  bcr.df[,(i+1)] <- bcr.out$use
  
  print(paste0("Finished bcr ", i, " of ", nrow(bcr.country)))
  
}

colnames(bcr.df) <- c("id", paste0(bcr.country$country, bcr.country$subUnit))

#FILTERING######################

#1. Set filters----

#Minimum year of data
minyear <- 1980

#Survey time
mintssr <- -1
maxtssr <- 6

#day of year
minday <- 135 #May 15
maxday <- 196 #July 15

#2. Remove points outside study area
bcr.in <- bcr.df[rowSums(ifelse(bcr.df==TRUE, 1, 0)) > 0,]

#3. Filter visits----
visit.use <- visit %>% 
  dplyr::filter(id %in% bcr.in$id,
                year >= minyear,
                tssr >= mintssr,
                tssr <= maxtssr,
                jday >= minday,
                jday <= maxday)

#4. Filter offset, bird, and bcr objects by remaining visits----
offsets.use <- offsets %>% 
  dplyr::filter(id %in% visit.use$id)

bird.use <- bird %>% 
  dplyr::filter(id %in% visit.use$id)

bcr.use <- bcr.in %>% 
  dplyr::filter(id %in% visit.use$id)

#SAMPLE SIZE AGGREGATION##############

#1. Summarize subunit sample sizes----
bcr.n <- data.frame(bcr = unique(colnames(bcr.use)[-1]),
                    n = colSums(bcr.use[-1]==TRUE)) %>% 
  arrange(n)

write.csv(bcr.n, file.path(root, "Regions", "SubUnitSampleSizes.csv"), row.names=FALSE)

#2. Plot sample sizes----
bcr.sf <- bcr.country %>% 
  mutate(bcr = paste0(country, "_", subUnit)) %>% 
  left_join(bcr.n)

ggplot(bcr.sf) +
  geom_sf(aes(group=bcr, fill=log(n))) +
  scale_fill_viridis_c()

ggsave(filename=file.path(root, "Figures", "SubUnitSampleSizes.jpeg"), width = 8, height = 4)

#3. Aggregate subunits----
#CA: merge 42 with 41
#USA: merge 42 and 3 with 41

#Newfoundland and BCR 2 still have < 1000 but merging may be less useful than retaining

#Also try:
#merging 82 (Newfoundland) with mainland 81 and comparing to only 82 with small sample size
#merging 2 (coastal AK) with 41 and comparing to only 41 with small sample size

bcr.agr <- bcr.use %>% 
  mutate(can4142 = ifelse(can41==TRUE | can42==TRUE, TRUE, FALSE),
         usa41423 = ifelse(usa42==TRUE | usa41==TRUE | usa3==TRUE, TRUE, FALSE),
         usa414232 = ifelse(usa42==TRUE | usa41==TRUE | usa3==TRUE | usa2==TRUE, TRUE, FALSE),
         can8182 = ifelse(can82==TRUE | can81==TRUE, TRUE, FALSE)) %>% 
  dplyr::select(-can41, -usa42, -usa3, -can42, -usa41)

#4. Check sample sizes again
bcr.n2 <- data.frame(bcr = unique(colnames(bcr.agr)[-1]),
                    n = colSums(bcr.agr[-1]==TRUE)) %>% 
  arrange(n)

write.csv(bcr.n2, file.path(root, "Regions", "SubUnitSampleSizes_AfterAggregation.csv"), row.names=FALSE)

#SPATIAL GRID FOR BOOTSTRAP THINNING##############

#1. Get list of BCRs----
bcrs <- sort(unique(colnames(bcr.agr)[-1]))

#2. Create grid for thinning----
grid <- dgconstruct(spacing = 2.5, metric=TRUE)

visit.grid <- visit.use %>% 
  mutate(cell = dgGEO_to_SEQNUM(grid, lon, lat)$seqnum)

length(unique(visit.grid$cell))

#UPDATE BIRD LIST BY BCR##############

#1. Set detection threshold for inclusion----
nmin <- 30

#2. Set up loop----
birdlist <- data.frame(bcr=bcrs)

for(i in 1:length(bcrs)){
  
  #3. Select visits within BCR----
  bcr.i <- bcr.agr[,c("id", bcrs[i])] %>%
    data.table::setnames(c("id", "use")) %>%
    dplyr::filter(use==TRUE)
  
  #4. Filter bird data----
  bird.i <- bird.use %>% 
    dplyr::filter(id %in% bcr.i$id)
  
  #5. Determine whether exceeds threshold----
  birdlist[i,c(2:ncol(bird.use))] <- colSums(bird.i[2:ncol(bird.use)]) > nmin
  
}

#6. Rename columns with bird ID----
colnames(birdlist) <- c("bcr", colnames(bird.use[2:ncol(bird.use)]))

#F. COVARIATE LOOKUP TABLE###############

#1. Get extraction methods lookup table----
meth <- readxl::read_excel(file.path(root, "NationalModels_V5_VariableList.xlsx"), sheet = "ExtractionLookup")

#2. Set up dataframe for variables that are global---
meth.global <- meth %>% 
  dplyr::filter(Extent=="global")

cov.global <- expand.grid(bcr=bcrs, cov = meth.global$Label) %>% 
  mutate(val = TRUE) %>% 
  pivot_wider(names_from=cov, values_from=val)

#3. Set up dataframe for variables that need prioritization----
#collapse AK & CONUS priority fields - only needed for extraction
meth.prior <- meth %>% 
  dplyr::filter(Extent!="global") %>% 
  mutate(Priority = as.numeric(str_sub(Priority, 1, 1))) %>% 
  dplyr::select(Category, Name, Label, Priority) %>% 
  unique() %>% 
  arrange(Category, Name, Priority)

cov.prior <- matrix(ncol=nrow(meth.prior), nrow=length(bcrs),
                    dimnames=list(bcrs, meth.prior$Label)) %>% 
  data.frame()

#4. Loop for subcateogries of variables that need prioritization----
subcat <- unique(meth.prior$Name)

#4. Set up bcr loop----
for(i in 1:length(bcrs)){
  
  #5. Filter data to BCR----
  visit.i <- bcr.agr[,c("id", bcrs[i])] %>%
    data.table::setnames(c("id", "use")) %>%
    dplyr::filter(use==TRUE) %>% 
    dplyr::select(-use) %>% 
    left_join(visit.use, by="id")
  
  #6. Set up Name loop----
  for(j in 1:length(subcat)){
    
    #7. Select covs to compare----
    meth.j <- dplyr::filter(meth.prior, Name %in% subcat[j])
    
    cov.j <- visit.i %>% 
      dplyr::select(all_of(meth.j$Label))
    
    #8. Determine which cov to use----
    #Count non-na and non-zero values per cov, round to nearest 100 to avoid slight differences, and choose highest priority to use for that Name
    na.j <- cov.j %>% 
      mutate(id = row_number()) %>% 
      pivot_longer(cols=all_of(meth.j$Label), names_to="Label", values_to="Value") %>% 
      dplyr::filter(!is.na(Value) & as.numeric(Value) > 0) %>% 
      group_by(Label) %>% 
      summarize(n = round(n(), -2)) %>% 
      ungroup() %>% 
      left_join(meth.prior, by="Label")
    
    use.j <- na.j %>% 
      dplyr::filter(n==max(n, na.rm=TRUE)) %>% 
      dplyr::filter(Priority==min(Priority, na.rm=TRUE)) %>% 
      mutate(use = 1) %>% 
      right_join(meth.j, by = join_by(Label, Category, Name, Priority)) %>% 
      mutate(use = ifelse(is.na(use), 0, use))
    
    #9. Update lookup dataframe----
    for(k in 1:nrow(use.j)){
      
      if(use.j$use[k]==1){
        cov.prior[bcrs[i],use.j$Label[k]] <- TRUE
      } else {
        cov.prior[bcrs[i],use.j$Label[k]] <- FALSE
      }
    }
    
    #10. Remove covariates with only 1 option and are in Canada----
    if(nrow(use.j)==1){
      cov.prior[20:34, use.j$Label[1]] <- FALSE
    }
    
  }
  
  print(paste0("Finished BCR ", i, " of ", length(bcrs)))
  
}

#11. Put them together----
covlist <- cbind(cov.global, cov.prior)

#SAVE#####

#1. Package----
#rename objects
#thin out covariate fields for minimizing RAM
#convert bird object to sparse matrix for minimizing RAM (doesn't save space to do this for other objects because they don't have enough zeros)
cov <- visit.use %>% 
  dplyr::select(-source, -organization, -project, -sensor, -equipment, -location, -buffer, -lat, -lon, -year, -date, -observer, -duration, -distance, -tssr, -jday) %>% 
  mutate(tagMethod = as.factor(tagMethod))

visit <- visit.use %>% 
  dplyr::select(id, source, organization, project, sensor, tagMethod, equipment, location, buffer, lat, lon, year, date, observer, duration, distance, tssr, jday) 

gridlist <- visit.grid %>% 
  dplyr::select(id, year, cell)

bird <- bird.use %>% 
  column_to_rownames("id") %>% 
  as.matrix() %>% 
  as("dgCMatrix")
offsets <- offsets.use
bcrlist <- bcr.agr

#2. Save----
save(visit, cov, bird, offsets, covlist, birdlist, bcrlist, gridlist, file=file.path(root, "Data", "04_NM5.0_data_stratify.R"))
