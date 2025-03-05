# ---
# title: National Models 5.0 - summarize validation
# author: Elly Knight
# created: March 5, 2024
# ---

#NOTES################################

# This script summarizes the output from the `10.Validate.R` script into a single dataframe to be included with model products

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #data wrangling
library(data.table) #rbinding lists to dataframe with filling unmatched columns

#2. Set root path----
root <- "G:/Shared drives/BAM_NationalModels5"

#3. Restrict scientific notation----
options(scipen=999999)

#INVENTORY#############################

#1. List of validated models-----
validated <- data.frame(file = list.files(file.path(root, "output", "validation")),
                        path = list.files(file.path(root, "output", "validation"), full.names = TRUE)) |> 
  separate(file, into=c("species", "boot", "filetype"))

#2. Get todo list----
loop <- validated |> 
  group_by(species) |> 
  summarize(boots = n()) |> 
  ungroup() |> 
  dplyr::filter(boots==10)

#SUMMARIZE###########

#1. Set up loop----
eval.out <- list()
for(i in 1:nrow(todo)){
  
  species.i <- loop$species[i]
  
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
eval <- rbindlist(eval.out, fill=TRUE)

#6. Save ----
write.csv(eval, file.path(root, "output", "packaged", "BAMV5ModelEvaluationSummary.csv"), row.names = FALSE)
