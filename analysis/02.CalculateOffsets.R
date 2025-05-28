# ---
# title: National Models 5.0 - calculate offsets
# author: Elly Knight
# created: December 22, 2022
# ---

#NOTES################################

# This script uses the qpad-offsets package, which requires downloading that R project from github (https://github.com/borealbirds/qpad-offsets).

# This script calculates offsets for ALL species available in QPAD, some of which may not fully satisfy the assumptions of QPAD. Further filtering of the species list should be conducted.

# This script uses the traditional QPAD method for offset calculation. Future versions of the models should consider the joint estimation approach for more accurate density estimation.

#This script uses a custom (slower) version of the qpad-offsets repo code to accommodate locations outside the covarage of the available rasters and areas where there is true sunrise.

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(lubridate) #temporal data wrangling
library(QPAD) #to get offset model estimates
library(terra) #to handle covariate extraction from qpad rasters
library(sf) #coordinate projection
library(lutz) #timezone lookup
library(suncalc) #required for QPAD for time since sunrise calculation
library(data.table) #collapse list to dataframe

library(maptools)
library(intrval)
library(raster)

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_NationalModels5"

#A. LOAD QPAD#####

#1. Estimates----
load_BAM_QPAD(version = 3)

#2. Raster data----
# from https://github.com/borealbirds/qpad-offsets
rlcc <- rast("C:/Users/elly/Documents/BAM/QPAD/qpad-offsets/data/lcc.tif")
rtree <- rast("C:/Users/elly/Documents/BAM/QPAD/qpad-offsets/data/tree.tif")
rd1 <- rast("C:/Users/elly/Documents/BAM/QPAD/qpad-offsets/data/seedgrow.tif")

#3. Source functions----
# adapted from https://github.com/borealbirds/qpad-offsets
# source("analysis/00.QPADfunctions.R")

#B. PREP DATA####

#1. Load data package from script 01----
load(file.path(root, "data", "01_NM5.0_data_clean.Rdata"))

#2. Project visit coordinates ----
set.seed(1234)
visit_v <- visit |> 
#  sample_n(1000) |> 
  st_as_sf(coords=c("lon", "lat"), crs=4326, remove=FALSE) |> 
  st_transform(crs(rtree)) |> 
  vect()

#3. Define lcc values ----
lcclevs <- data.frame(lcc = seq(0, 19, 1),
                      class = c(NA, "Conif", "Conif", NA, NA, "DecidMixed", "DecidMixed", NA, "Open", NA, "Open", "Open", "Open", "Open", "Wet", "Open", "Open", "Open", NA, NA)) |> 
  mutate(lcc2 = case_when(class %in% c("Conif", "DecidMixed") ~ "Forest",
                          class %in% c("Open", "Wet") ~ "OpenWet"))

#4. Extract raster values ----
visit_r <- terra::extract(rlcc, visit_v, bind=TRUE) |> 
  terra::extract(x=rtree, bind=TRUE) |> 
  terra::extract(x=rd1, bind=TRUE) |> 
  st_as_sf() |> 
  st_drop_geometry() |> 
  left_join(lcclevs) |> 
  rename(datetime = date) |> 
  mutate(date = as.Date(datetime),
         tree = tree/100,
         tree = case_when(tree < 0 ~ 0,
                          tree > 1 ~ 0,
                          !is.na(tree) ~ tree))

#5. Get time zones ----
#look up timezones with lutz package
visit_r$tz <- tz_lookup_coords(visit_r$lat, visit_r$lon, method="fast")

#calculate sunrise for each timezone
tzs <- unique(visit_r$tz)
visit_s <- data.frame()
for(i in 1:length(tzs)){
  
  #get just the data for that timezone
  visit_i <- dplyr::filter(visit_r, tz==tzs[i]) |>
    mutate(datetime2 = ymd_hms(datetime, tz=tzs[i]),
           date = as.Date(datetime2, tz=tzs[i]))
  
  #get sunrise
  visit_s <- getSunlightTimes(data=visit_i, keep=c("sunrise", "nadir")) |> 
    cbind(visit_i |> dplyr::select(-c(date, lat, lon))) |> 
    mutate(sr = case_when(is.na(sunrise) ~ nadir,
                          !is.na(sunrise) ~ sunrise),
           tssr = as.numeric(difftime(datetime2, sr, units="hours"))) |> 
    rbind(visit_s)
  
  cat("Finished time zone", i, "of", length(tzs), "\n")
}

#6. Finish formatting -----
visit_x <- visit_s |> 
  rename(TREE = tree, LCC2 = lcc2, LCC4 = class, DSLS = seedgrow) |> 
  mutate(TSSR = tssr/24,
         JDAY = yday(datetime)/366,
         LCC4 = factor(LCC4, levels=c("DecidMixed", "Conif", "Open", "Wet")),
         LCC2 = factor(LCC2, levels=c("DecidMixed", "Conif")),
         MAXDUR = duration,
         MAXDIS = distance/100)

#7. Check tssr distribution ----
#just take a sample to speed it up
set.seed(1234)
ggplot(visit_x |> sample_n(1000)) +
  geom_hex(aes(x=lon, y=tssr)) +
  geom_smooth(aes(x=lon, y=tssr)) +
  facet_wrap(~source) +
  scale_fill_viridis_c()

set.seed(1234)
ggplot(visit_x |> sample_n(1000)) +
  geom_histogram(aes(x=tssr)) +
  geom_vline(aes(xintercept=0), linetype="dashed") +
  facet_wrap(~source)

#D. CALCULATE OFFSETS####

#1. Get list of species----
spp <- getBAMspecieslist()

#2. Set up output & loop----
offsets <- data.frame(id=visit_x$id)

for (i in 1:length(spp)) {
  
  #3. Get best model & coeff----
  mi <- bestmodelBAMspecies(spp[i], type="BIC")
  cfi <- coefBAMspecies(spp[i], mi$sra, mi$edr)
  ## constant for NA cases
  cf0 <- exp(unlist(coefBAMspecies(spp[i], 0, 0)))
  
  #4. Make Design matrices for singing rates (`Xp`) and for EDR (`Xq`) -----
  Xp <- visit_x |>
    mutate('(Intercept)' = 1,
           TSSR2 = TSSR^2,
           JDAY2 = JDAY^2,
           DSLS2 = DSLS^2) |> 
    dplyr::select('(Intercept)', TSSR, JDAY, TSSR2, JDAY2, DSLS, DSLS2) |> 
    as.matrix()
  
  Xq <- visit_x |> 
    mutate('(Intercept)' = 1,
           LCC2OpenWet=ifelse(LCC4 %in% c("Open", "Wet"), 1, 0),
           LCC4Conif=ifelse(LCC4=="Conif", 1, 0),
           LCC4Open=ifelse(LCC4=="Open", 1, 0),
           LCC4Wet=ifelse(LCC4=="Wet", 1, 0)) |> 
    dplyr::select('(Intercept)', TREE, LCC2OpenWet, LCC4Conif, LCC4Open, LCC4Wet) |> 
    as.matrix()
  
  p <- rep(NA, nrow(visit_x))
  A <- q <- p
  
  #5. Update design matrices to match the coefs ----
  Xp2 <- Xp[,names(cfi$sra), drop=FALSE]
  if("DSLS" %in% colnames(Xp2)){OKp <- !is.na(Xp2[,"DSLS"])} else {OKp <- rowSums(is.na(Xp2))}
  Xq2 <- Xq[,names(cfi$edr), drop=FALSE]
  OKq <- rowSums(is.na(Xq2)) == 0
  
  #6. calculate p, q, and A based on constant phi and tau for the respective NAs----
  p[!OKp] <- sra_fun(visit_x$MAXDUR[!OKp], cf0[1])
  unlim <- ifelse(visit_x$MAXDIS[!OKq] == Inf, TRUE, FALSE)
  A[!OKq] <- ifelse(unlim, pi * cf0[2]^2, pi * visit_x$MAXDIS[!OKq]^2)
  q[!OKq] <- ifelse(unlim, 1, edr_fun(visit_x$MAXDIS[!OKq], cf0[2]))
  
  #7. Calculate time/lcc varying phi and tau for non-NA cases----
  phi1 <- exp(drop(Xp2[OKp,,drop=FALSE] %*% cfi$sra))
  tau1 <- exp(drop(Xq2[OKq,,drop=FALSE] %*% cfi$edr))
  p[OKp] <- sra_fun(visit_x$MAXDUR[OKp], phi1)
  unlim <- ifelse(visit_x$MAXDIS[OKq] == Inf, TRUE, FALSE)
  A[OKq] <- ifelse(unlim, pi * tau1^2, pi * visit_x$MAXDIS[OKq]^2)
  q[OKq] <- ifelse(unlim, 1, edr_fun(visit_x$MAXDIS[OKq], tau1))
  
  #8. Fix 0s----
  #log(0) is not a good thing, apply constant instead
  ii <- which(p == 0)
  p[ii] <- sra_fun(visit_x$MAXDUR[ii], cf0[1])
  
  #9. Calculate offset -----
  offsets[,i+1] = round(log(p) + log(A) + log(q), 4)

  cat(spp[i], "\n")
  
}
colnames(offsets) <- c("id", spp)

#10. Fix GRAJ####
offsets <- rename(offsets, CAJA = GRAJ)

#11. Reorder by ID-----
offsets <- offsets[match(visit$id, offsets$id),]

#12. Clean up the visit object ----
visit <- visit |> 
  left_join(visit_x |> 
              dplyr::select(id, tssr)) |> 
  mutate(jday = yday(date))

#E. SAVE####

#1. Save----
save(visit, bird, offsets, file=file.path(root, "data", "02_NM5.0_data_offsets.Rdata"))
