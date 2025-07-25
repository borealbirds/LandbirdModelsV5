# ---
# title: National Models 5.0 - stratify by BCR*country
# author: Elly Knight
# created: September 29, 2023
# ---

#NOTES################################

#In this script, we stratify & prepare the data for modelling. Steps include:

#1. BCR Attribution: Diving the data into BCRs for separate models. Each BCR is buffered by 100km. We do this so that we can feather predictions from adjacent regions together. The exception is the international boundaries between US and Canada, which we don't buffer because spatial layers for the covariates are different on either side of the border. In this case, we use a shapefiles of country boundaries to intersect with the buffered regions.

#2. Filtering the data. We remove data outside the study area, data before a minimum year, and data outside of dawn surveys. Additional dataset filtering steps can be added here as necessary.

#3. Aggregating subunits based on minimum sample size threshold. Sample sizes are very low in some subunits in Alaska and northwestern BC.

#4. Thinning. First we assign each survey to a grid cell (1.9 km spacing) Then we randomly selecting one survey per cell per year for each bootstrap. Note on line 186-192, we use an alternative to sample_n(), which does provide a truly random sample and biases towards lines of data near the top of the dataframe. This step takes a while to run.

#5. Indicating species that should be modelled for each subunit based on a minimum detection threshold.

#6. Determining which covariates should be used for each subunit. This is done using the same extraction lookup table as the '03.ExtractCovariates.R' script. For each subunit, we choose the version of a covariate with the highest priority that has full coverage for that subunit. We then rename the covariates so that they match the naming conventions of the prediction layers for the prediction step.

#7. Using variance inflation factor (VIF) with a threshold of 10 to remove highly collinear variables to improve efficiency and interpretibility.

#Removal of non-dawn surveys is a secondary line of defense in case temporal patterns of QPAD offsets are unrealistic. In other words, we are restricting the dataset to times of day when p(availability) should be relatively high. This filter should be removed once QPAD offsets have been revisited. It should also be revisited if/when the national models include nocturnal species.

#To consider for V6:
#stratification across years for year-specific validation.
#Explore tidy models

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(terra) #basic raster handling
library(sf) #basic shapefile handling
library(fasterize) #fast rasterization of shapefiles
library(exactextractr) #fast & efficient raster extraction
library(dggridR) #grid for spatial thinning
library(Matrix) #for sparse matrix conversion
library(usdm) #VIF

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_NationalModels5"

#3. Load data packages with offsets and covariates----
load(file.path(root, "data", "03_NM5.0_data_covariates.Rdata"))

#4. Turn off scientific notation---
options(scipen=99999)

#BCR ATTRIBUTION#####################

#1. Make visit a vector object----
visit.v <- visit |> 
  dplyr::select(id, lat, lon) |> 
  st_as_sf(coords=c("lon", "lat"), crs=4326) |> 
  st_transform(5072) |> 
  vect()

#2. Read in buffered & unbuffered shapefiles----
bcr <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel_Unbuffered.shp"))
bcr.buff <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel_Buffered.shp"))
  
#3. Set up loop for BCR attribution----
bcr.df <- data.frame(id=visit$id)
for(i in 1:nrow(bcr.buff)){
  
  #4. Filter bcr shapefile----
  bcr_buff.i <- bcr.buff |> 
    dplyr::filter(row_number()==i)
  bcr.i <- bcr |> 
    dplyr::filter(subUnit==bcr_buff.i$subUnit,
                  country==bcr_buff.i$country)
  
  #5. Convert to raster for fast extraction----
  r_buff <- rast(ext(bcr_buff.i), resolution=1000, crs=crs(bcr_buff.i))
  r <- rast(ext(bcr.i), resolution=1000, crs=crs(bcr.i))
  
  bcr_buff.r <- rasterize(x=bcr_buff.i, y=r_buff, field="subUnit")
  bcr.r <- rasterize(x=bcr.i, y=r, field="subUnit")
  
  #6. Extract raster value----
  bcr_buff.out <- data.frame(subUnit=extract(x=bcr_buff.r, y=visit.v)[,2]) |> 
    mutate(buff = ifelse(!is.na(subUnit), TRUE, FALSE))
  bcr.out <- data.frame(subUnit=extract(x=bcr.r, y=visit.v)[,2]) |> 
    mutate(bcr = ifelse(!is.na(subUnit), TRUE, FALSE))
  
  #7. Wrangle to 3 classes: 0 = exclude from model, 1 = use in model, 2 = use in validation
  bcr.status <- bcr_buff.out |> 
    mutate(bcr = bcr.out$bcr) |> 
    mutate(use = case_when(buff==FALSE & bcr==FALSE ~ 0,
                           bcr==FALSE & buff==TRUE ~ 2,
                           buff==TRUE ~ 1))
  
  bcr.df[,(i+1)] <- bcr.status$use
  
  print(paste0("Finished bcr ", i, " of ", nrow(bcr.buff)))
  
}

colnames(bcr.df) <- c("id", paste0(bcr$country, bcr$subUnit))

#FILTERING######################

#1. Set filters----

#Minimum year of data
minyear <- 1980

#Survey time
mintssr <- -1
maxtssr <- 6

#day of year
minday <- 152 #June 1
maxday <- 196 #July 15

#2. Remove points outside study area
bcr.in <- bcr.df |> 
  dplyr::filter(rowSums(across(-id)) > 0) |> 
  mutate(id = as.integer(id))

#3. Filter visits----
visit.use <- visit |> 
  dplyr::filter(id %in% bcr.in$id,
                year >= minyear,
                tssr >= mintssr,
                tssr <= maxtssr,
                jday >= minday,
                jday <= maxday,
                duration > 0)

#4. Filter offset, bird, and bcr objects by remaining visits----
offsets.use <- offsets |> 
  dplyr::filter(id %in% visit.use$id)

bird.use <- bird |> 
  dplyr::filter(id %in% visit.use$id)

bcr.use <- bcr.in |> 
  dplyr::filter(id %in% visit.use$id)

#SAMPLE SIZE CHECK##############

#1. Summarize subunit sample sizes----
bcr.n <- data.frame(bcr = unique(colnames(bcr.use)[-1]),
                    n = colSums(bcr.use[-1]>0)) |> 
  arrange(n)

write.csv(bcr.n, file.path(root, "Regions", "SubUnitSampleSizes.csv"), row.names=FALSE)

#2. Plot sample sizes----
bcr.sf <- bcr.buff |> 
  mutate(bcr = paste0(country, subUnit)) |> 
  left_join(bcr.n)

ggplot(bcr.sf) +
  geom_sf(aes(group=bcr, fill=log(n))) +
  scale_fill_viridis_c()

ggsave(filename=file.path(root, "Figures", "SubUnitSampleSizes.jpeg"), width = 8, height = 4)

#BOOTSTRAP THINNING##############

#1. Get list of BCRs----
bcrs <- sort(unique(colnames(bcr.use)[-1]))

#2. Create grid for thinning----
grid <- dgconstruct(spacing = 2.5, metric=TRUE)

visit.grid <- visit.use |> 
  mutate(cell = dgGEO_to_SEQNUM(grid, lon, lat)$seqnum)

length(unique(visit.grid$cell))

#3. Set number of bootstraps ----
boots <- 32

#4. Make the bootstrap list ----
#Choose one sample id for each unique grid cell & year ----
for(i in 1:boots){
  
  set.seed(i)
  
  boot.i <- visit.grid |> 
    dplyr::select(id, year, cell) |> 
    group_by(year, cell) |> 
    mutate(rowid = row_number(),
           use = sample(1:max(rowid), 1)) |> 
    ungroup() |> 
    dplyr::filter(rowid==use) |> 
    dplyr::select(cell, year, id) |> 
    arrange(cell, year)
  
  if(i==1){
    bootlist <- boot.i
  } else {
    bootlist[,i+2] <- boot.i$id
  }
  
  cat(i, "    ")

}

#5. Fix column names ----
colnames(bootlist) <- c("cell", "year", paste0("b", seq(1:boots)))

#UPDATE BIRD LIST BY BCR##############

#1. Set detection threshold for inclusion----
nmin <- 30

#2. Set up loop----
birdlist <- data.frame(bcr=bcrs)

for(i in 1:length(bcrs)){
  
  #3. Filter bird data by BCR & bootstrap 1----
  bird.i <- bcr.use[,c("id", bcrs[i])] |>
    data.table::setnames(c("id", "use")) |>
    dplyr::filter(use > 0) |> 
    dplyr::select(-use) |> 
    left_join(bird.use, by="id") |> 
    dplyr::filter(id %in% bootlist$b1) |> 
    dplyr::select(all_of(colnames(offsets)))
  
  #4. Determine whether exceeds threshold----
  birdlist[i,c(2:ncol(bird.i))] <- sapply(bird.i[2:ncol(bird.i)], function(x) sum(x > 0, na.rm=TRUE)) > nmin
  
}

#5. Rename columns with bird ID----
colnames(birdlist) <- c("bcr", colnames(offsets[2:ncol(offsets)]))

#F. COVARIATE LOOKUP TABLE###############

#1. Get extraction methods lookup table----
meth <- readxl::read_excel(file.path(root, "NationalModels_V5_VariableList.xlsx"), sheet = "ExtractionLookup") |> 
  dplyr::filter(Use==1)

#2. Set up dataframe for variables that are global---
meth.global <- meth |> 
  dplyr::filter(Global==1)

cov.global <- expand.grid(bcr=bcrs, cov = meth.global$Label) |> 
  mutate(val = TRUE) |> 
  pivot_wider(names_from=cov, values_from=val)

#3. Set up dataframe for variables that are in Canada----
meth.can <- meth |> 
  dplyr::filter(Canada==1)

cov.can <- expand.grid(bcr=bcrs, cov = meth.can$Label) |> 
  mutate(val = TRUE) |> 
  pivot_wider(names_from=cov, values_from=val)
cov.can[str_detect(bcrs, "usa"), -1] <- FALSE

#4. Set up dataframe for variables that are in US----
meth.us <- meth |> 
  dplyr::filter(US==1)

cov.us <- expand.grid(bcr=bcrs, cov = unique(meth.us$Label)) |> 
  mutate(val = TRUE) |> 
  pivot_wider(names_from=cov, values_from=val)
cov.us[str_detect(bcrs, "can"), -1] <- FALSE

#5. Set up dataframe for variables that are in Alaska----
meth.ak <- meth |> 
  dplyr::filter(US_AK==1)

cov.ak <- expand.grid(bcr=bcrs, cov = meth.ak$Label) |> 
  mutate(val = TRUE) |> 
  pivot_wider(names_from=cov, values_from=val)
cov.ak[str_detect(bcrs, "can"), -1] <- FALSE
cov.ak[!(bcrs %in% c("usa2", "usa41423", "usa5", "usa43", "usa40")),-1] <- FALSE

#6. Set up dataframe for variables that are in CONUS----
meth.conus <- meth |> 
  dplyr::filter(US_CONUS==1)

cov.conus <- expand.grid(bcr=bcrs, cov = meth.conus$Label) |> 
  mutate(val = TRUE) |> 
  pivot_wider(names_from=cov, values_from=val)
cov.conus[str_detect(bcrs, "can"), -1] <- FALSE
cov.conus[bcrs %in% c("usa2", "usa41423", "usa5", "usa43", "usa40"),-1] <- FALSE

#7. Set up VLCE dataframe----
meth.vlce <- meth |> 
  dplyr::filter(Label=="VLCE_1km")

cov.vlce <- expand.grid(bcr=bcrs, cov = meth.vlce$Label) |> 
  mutate(val = TRUE) |> 
  pivot_wider(names_from=cov, values_from=val)
cov.vlce[str_detect(bcrs, "usa"), -1] <- FALSE
cov.vlce[bcrs %in% c("can11"), -1] <- FALSE

#8. Put all together----
covlist.in <- cbind(cov.global, cov.can[,-1], cov.us[,-1], cov.ak[,-1], cov.conus[,-1], cov.vlce[,-1])
rownames(covlist.in) <- covlist.in$bcr

#9. Set up bcr loop----
covlist <- covlist.in
vif <- list()
for(i in 1:length(bcrs)){
  
  #10. Filter data to BCR----
  visit.i <- bcr.use[,c("id", bcrs[i])] |>
    data.table::setnames(c("id", "use")) |>
    dplyr::filter(use > 0) |> 
    dplyr::select(-use) |> 
    left_join(visit.use, by="id")
  
  #11. Get covariate list for this BCR-----
  covlist.i <- covlist.in |> 
    dplyr::filter(bcr==bcrs[i]) |> 
    pivot_longer(-bcr, names_to="cov", values_to="use") |> 
    dplyr::filter(use==TRUE)

  #12. Get the covariates for bootstrap 1----
  #Drop factors, use them always
  cov.i <- visit.i |> 
    dplyr::filter(id %in% bootlist$b1) |> 
    dplyr::select(all_of(covlist.i$cov)) |> 
    select_if(is.numeric) |> 
    data.frame()
  
  #14. Run VIF----
  #Keep height,road, and alan variables for all BCRs as well as % deciduous and % coniferous for Canada (SCANFI) to allow consistent interpretation of models across BCRs
  vif.i <- vifstep(cov.i, th=10, keep=c("SCANFIheight_1km", "SCANFIheight_5x5", "SCANFIprcC_1km", "SCANFIprcC_5x5", "SCANFIprcD_1km", "SCANFIprcD_5x5", "LFheigth_1km", "LFheigth_5x5", "ETHheight_1km", "ETHheight_5x5", "canroad_1km", "canroad_5x5", "usaroad_1km", "usaroad_5x5", "CCNL_1km"))
  
  #15. Update covariate list----
  covlist[bcrs[i], vif.i@excluded] <- FALSE
  
  #16. Save the vif info to look at later----
  vif[[i]] <- vif.i
  save(vif, file=file.path(root, "data", "Covariates", "VIF.Rdata"))
  
  print(paste0("Finished BCR ", i, " of ", length(bcrs)))
  
}

#SAVE#####

#1. Package----
#rename objects
#thin out covariate fields for minimizing RAM
#convert bird object to sparse matrix for minimizing RAM (doesn't save space to do this for other objects because they don't have enough zeros)
cov <- visit.use |> 
  dplyr::select(-source, -organization, -project, -sensor, -equipment, -location, -buffer, -lat, -lon, -year, -date, -observer, -duration, -distance, -tssr, -jday, -tagMethod, -method, -country)

visit <- visit.use |> 
  dplyr::select(id, source, organization, project, sensor, tagMethod, method, equipment, location, country, buffer, lat, lon, year, date, observer, duration, distance, tssr, jday) 

bird <- bird.use |> 
  dplyr::select(all_of(colnames(offsets))) |> 
  column_to_rownames("id") |> 
  as.matrix() |> 
  as("dgCMatrix")

offsets <- offsets.use
bcrlist <- bcr.use

#2. Some checks ----

#Site IDs
sites <- bcrlist |> 
  pivot_longer(-id, names_to="bcr", values_to="use") |> 
  dplyr::filter(use==TRUE) |> 
  left_join(visit |> dplyr::select(id, lat, lon))

ggplot(sites |> sample_n(100000)) +
  geom_point(aes(x=lon, y=lat, colour=bcr))

#3. Save----
save(visit, cov, bird, offsets, covlist, birdlist, bcrlist, bootlist, file=file.path(root, "Data", "04_NM5.0_data_stratify.Rdata"))
