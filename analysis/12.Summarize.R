# ---
# title: National Models 5.0 - summarize validation
# author: Elly Knight
# created: March 5, 2024
# ---

#NOTES################################

# This script summarizes the output from the `10.Validate.R` alongside additional model metadata.

# The output is a .xlsx file similar to the V4 "results" object at file.path(root, "data", "Lookups", "BAMv4-results-2025-05-21.xlsx"))

# The metadata table should be edited mannually as a csv outside of R from the previous version

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #data wrangling
library(data.table) #rbinding lists to dataframe with filling unmatched columns
library(openxlsx) #read and save to excel
library(sf) #shapefiles
library(BAMexploreR) #variable importance & pop est function

#2. Set root path----
root <- "G:/Shared drives/BAM_NationalModels5"

#3. Restrict scientific notation----
options(scipen=9999)

#4. Load dataset ----
load(file.path(root, "data", "04_NM5.0_data_stratify.Rdata"))

#SPECIES###########

#1. Get WT list ----
species.wt <- read.csv(file.path(root, "data", "lookups", "lu_species.csv")) |> 
  rename(id = species_code,
         scientific = species_name,
         english = species_common_name,
         french = species_french_name,
         family = species_family) |>
  dplyr::select(id, scientific, english, french, family)

#2. Get list of packaged species----
packaged <- data.frame(file = list.files(file.path(root, "output", "10_packaged"), pattern="*.tif", recursive=TRUE))  |> 
  separate(file, into=c("sppfolder", "bcrfolder", "spp", "bcr", "year", "filetype"), remove=FALSE) |>  
  mutate(year = as.numeric(year)) |> 
  dplyr::filter(sppfolder!="extrapolation")

#3. Get list of validatd species----
validated <- data.frame(file = list.files(file.path(root, "output", "11_validation"), pattern=".Rdata")) |> 
  separate(file, into=c("spp", "filetype")) |> 
  dplyr::select(spp)

#4. Get list of species ready to summarize----
spp <- packaged |> 
  dplyr::select(spp) |> 
  unique() |> 
  inner_join(validated)

#5. Join to species info from v4 ----
species <- species.wt |> 
  inner_join(spp |> 
               rename(id = spp)) |> 
  dplyr::select(all_of(colnames(species.wt))) |> 
  dplyr::filter(id %in% c("BBWA", "BBWO", "BLPW", "CAWA", "CONW", "LEYE", "OSFL", "OVEN", "RUBL", "SOSA", "TEWA"))

#6. Check against to do list ----
missing <- data.frame(id=colnames(dplyr::select(birdlist, -bcr))) |> 
  anti_join(species)

#REGIONS############

#1. Get BCR shapefile for modelling ----
bcr_mod <- read_sf(file.path(root, "gis", "Subregions_unbuffered.shp"))

#2. Get original BCR shapefile for BCR names ----
bcr_og <- read_sf(file.path(root, "Regions", "archived", "BCR_Terrestrial_master.shp")) |> 
  st_drop_geometry() |> 
  mutate(name = str_to_title(BCRNAME)) |> 
  dplyr::select(BCR, name) |> 
  unique() |> 
  rename(bcr_link = BCR)

#3. Put together ----
bcr <- bcr_mod |> 
  st_drop_geometry() |> 
  rename(region = bcr,
         bcr = subUnit) |> 
  mutate(country = ifelse(country=="can", "Canada", "USA"),
         area_km2 = round(as.numeric(st_area(bcr_mod))/1000000),
         bcr_link = case_when(bcr %in% c(40:42) ~ 4,
                              bcr %in% c(60:62) ~ 6,
                              bcr %in% c(70:72) ~ 7,
                              bcr %in% c(80:82) ~ 8,
                              !is.na(bcr) ~ bcr)) |> 
  left_join(bcr_og) |> 
  mutate(name = case_when(bcr==41423 ~ "Northwestern Interior Forest & Artic Plains And Mountains",
                          !is.na(name) ~ name),
         type = "Single model") |> 
  dplyr::select(region, type, country, bcr, name, area_km2)

#4. Make mosaic list ----
area.ca <- bcr |> 
  dplyr::filter(country=="Canada") |> 
  summarize(area = sum(area_km2))

area.ak <- bcr |> 
  dplyr::filter(region %in% c("usa41423", "usa2", "usa40", "usa43", "usa5")) |> 
  summarize(area = sum(area_km2))

area.48 <- bcr |> 
  dplyr::filter(region %in% c("usa5", "usa9", "usa10", "usa11", "usa13", "usa14", "usa23", "usa28")) |> 
  summarize(area = sum(area_km2))

mosaic <- data.frame(region=c("Canada", "Alaska", "Lower48"),
                     type = rep("Mosaic", 3),
                     country=c("Canada", "USA", "USA"),
                     bcr = rep(NA, 3),
                     name = c("Canada mosaic", "Alaska mosaic", "Lower 48 hemiboreal mosaic"),
                     area_km2 = c(area.ca$area, area.ak$area, area.48$area))

#5. Put together ----
regions <- rbind(bcr, mosaic)

#VARIABLES########

#1. Get it from google drive----
variables <- read.csv(file.path(root, "data", "Lookups", "covariate_metadata.csv"))

#IMPORTANCE#########

#1. Get it from BAMexploreR ----
imp_bamx <- bam_predictor_importance_v5

#2. Wrangle ----
importance <- imp_bamx |> 
  rename(id = spp,
         region = bcr,
         variable = predictor,
         importance_mean = mean_rel_inf,
         importance_sd = sd_rel_inf) |> 
  inner_join(species) |> 
  dplyr::select(id, scientific, english, region, variable, importance_mean, importance_sd) |> 
  arrange(id, region, -importance_mean)
  
#VALIDATION#############################

#1. Set up loop----
eval.out <- list()
for(i in 1:nrow(spp)){
  
  species.i <- spp$spp[i]
  
  #2. Load validation output ----
  load(file.path(root, "output", "11_validation", paste0(species.i, ".Rdata")))
  
  #4. Summarize----
  eval.out[[i]] <- rbindlist(eval, use.names = TRUE, idcol = "region") |> 
#    dplyr::select(-boot) |> 
    mutate(pseudo_R2 = ifelse(pseudo_R2 < 0, 0, pseudo_R2)) |> 
    group_by(region) |> 
    summarize(prevalence = mean(prevalence),
           AUC_init = mean(AUC_init),
           AUC_final = mean(AUC_final),
           pseudo_R2 = mean(pseudo_R2)) |> 
    ungroup() |> 
    left_join(rbindlist(oc, idcol = "region") |> 
                dplyr::select(region, occc, oprec, oaccu)) |> 
    mutate(id = species.i)
  
  cat(i, " ")

}

#5. Final object----
validation <- rbindlist(eval.out, fill=TRUE) |> 
  inner_join(species) |> 
  dplyr::select(id, scientific, english, region, prevalence, AUC_init, AUC_final, pseudo_R2, occc, oprec, oaccu)

#ABUNDANCES##########

#1. Set up the function ----
pop_sum <- function(file, type, region, id, year){
  
  #read in raster stack of bootstraps, remove the buffer (single models only)
  if(type=="Single model"){ 
    bcr <- dplyr::filter(bcr_mod, bcr==region) |> 
      vect()
    r <- rast(file.path(root, "output", "07_predictions", file)) |> 
      crop(bcr, mask=TRUE)
    }
  if(type=="Mosaic"){r <- rast(file.path(root, "output", "08_mosaics", file))}
  
  #get the truncation value from the packaged raster
  p <- rast(file.path(root, "output", "10_packaged", id, region,
                      paste0(id, "_", region, "_", year, ".tif")))
  q99 <- global(p$mean, max, na.rm=TRUE)[,1]
  
  #truncate
  t <- clamp(r, upper = q99, values=TRUE)
  
  #estimate population
  abun <- global(t * 1 * 100, sum, na.rm=TRUE)[,1]
  
  #estimate density
  mn <- global(t, mean, na.rm=TRUE)[,1]
  
  #results
  out <- data.frame(
    id = id,
    region = region,
    year = year,
    abundance_estimate = round(median(abun)/1e6, 4),
    population_lower = round(quantile(abun, 0.05)/1e6, 4),
    population_upper = round(quantile(abun, 0.95)/1e6, 4),
    density_estimate = round(median(mn), 4),
    density_lower = round(quantile(mn, 0.05), 4),
    density_upper = round(quantile(mn, 0.95), 4)
  )
  
  return(out)
  
}

#2. Get BCR shapefile for modelling ----
bcr_mod <- read_sf(file.path(root, "gis", "Subregions_unbuffered.shp")) |> 
  st_transform(crs=5072)

#3. List of models to run ----
pred <- data.frame(file = list.files(file.path(root, "output", "07_predictions"), recursive = TRUE)) |> 
  separate(file, into=c("f1", "id", "region", "year", "filetype"), remove=FALSE) |> 
  dplyr::select(-f1, -filetype) |> 
  mutate(type = "Single model")

mosaic <- data.frame(file = list.files(file.path(root, "output", "08_mosaics"), recursive = TRUE)) |> 
  separate(file, into=c("region", "f1", "id", "year", "filetype"), remove=FALSE) |> 
  dplyr::select(-f1, -filetype) |> 
  mutate(type = "Mosaic")
  
todo <- rbind(pred, mosaic) |> 
  inner_join(species) |> 
  dplyr::select(file, type, region, id, year) |> 
  dplyr::filter(region=="Canada")

#4. Apply function ----
abundance_out <- purrr::pmap(todo, pop_sum)

#5. Final object ----
abundances <- rbindlist(abundance_out) |> 
  left_join(species) |> 
  dplyr::select(id, scientific, english, region, year, population_estimate, population_lower, population_upper, density_estimate, density_lower, density_upper)

#METADATA#########

#1. Get metadata -----
metadata <- read.csv(file.path(root, "data", "Lookups", "BAMv5-results-metadata.csv"))

#PACKAGE#########

#1. Put it together ----
out <- list(species, regions, variables, importance, validation, abundances)
names(out) <- c("species", "regions", "variables", "importance", "validation", "abundances")

write.xlsx(out, file = file.path(root, "output", "12_BAMV5-results_example.xlsx"))
