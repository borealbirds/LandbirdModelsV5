# ---
# title: National Models 5.0 - summarize validation
# author: Elly Knight
# created: March 5, 2024
# ---

#NOTES################################

# This script summarizes the output from the `10.Validate.R` alongside additional model metadata.

# The output is a .xlsx file similar to the V4 "results" object.

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #data wrangling
library(data.table) #rbinding lists to dataframe with filling unmatched columns
library(readxl)

#2. Set root path----
root <- "G:/Shared drives/BAM_NationalModels5"

#3. Restrict scientific notation----
options(scipen=999999)

#SPECIES###########

#1. Get V4 object -----
species.4 <- read_excel(file.path(root, "data", "Lookups", "BAMv4-results-2025-05-21.xlsx"), sheet="species")

#2. Get list of packaged species----
packaged <- data.frame(file = list.files(file.path(root, "output", "10_packaged"), pattern="*.tif", recursive=TRUE))  |> 
  separate(file, into=c("sppfolder", "bcrfolder", "spp", "bcr", "year", "filetype"), remove=FALSE) |> 
  mutate(year = as.numeric(year)) |> 
  dplyr::filter(sppfolder!="extrapolation",
                bcr=="mosaic",
                year==2020) |> 
  dplyr::select(spp)

#3. Get list of validatd species----
validated <- data.frame(file = list.files(file.path(root, "output", "11_validated"), pattern=".Rdata")) |> 
  separate(file, into=c("spp", "filetype")) |> 
  dplyr::select(spp)

#4. Get list of species ready to summarize----
spp <- inner_join(packaged, validated)

#5. Join to species info from v4 ----
species <- species.4 |> 
  inner_join(spp |> 
               rename(id = spp))

#REGIONS############


#VARIABLES########

#1. Get it from google drive----
variables <- read.csv(file.path(root, "data", "Lookups", "covariate_metadata.csv"))

#IMPORTANCE#########

#VALIDATION#############################

#TO DO:UPDATE FOR THE NEW FORMAT#######

#1. Set up loop----
eval.out <- list()
for(i in 1:nrow(spp)){
  
  species.i <- spp$spp[i]
  
  #2. List of Rdata files----
  loop.i <- dplyr::filter(validated, species==species.i)
  
  #3. Read them in----
  eval.list <- list()
  for(j in 1:nrow(loop.i)){
    
    boot.j <- loop.i$boot[j]
    
    load(loop.i$path[j])
    
    eval.list[[j]] <- rbindlist(eval, fill=TRUE) |> 
      mutate(boot = boot.j)
    
  }
  
  #4. Summarize----
  eval.out[[i]] <- rbindlist(eval.list, fill=TRUE) |> 
    dplyr::select(-boot) |> 
    pivot_longer(-c(bcr, spp), names_to="metric", values_to="value") |> 
    group_by(spp, bcr, metric) |> 
    summarize(mean = mean(value),
              sd = sd(value)) |> 
    ungroup()
  
}

#5. Final object----
validation <- rbindlist(eval.out, fill=TRUE)

#ABUNDANCES##########



#DENSITIES#############


#METADATA#########

#1. Get V4 object -----
metadata.4 <- read_excel(file.path(root, "data", "Lookups", "BAMv4-results-2025-05-21.xlsx"), sheet="metadata")

#PACKAGE#########