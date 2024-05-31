# ---
# title: National Models 5.0 - covariate interpretation
# author: Mannfred Boehm
# created: May 22, 2024
# ---


#NOTES################################

# This script extracts covariate contributions to model predictions. 

# The outputs will be analysed for identifying covariates of importance across several metrics (e.g. BCR, ecology, etc.)



#PREAMBLE############################

#1. Attach packages----
library(gbm)
library(tidyverse)


#2. Create a list of dataframes containing relative influence per covariate----
#   Every list element represents a bootstrap sample 

# connect to BAM Drive and find bootstrap files 
root <- "G:/Shared drives/BAM_NationalModels/NationalModels5.0/"

gbm_objs <- list.files(file.path(root, "output", "bootstraps"))

gbm_objs_test <- gbm_objs[1:3] # for testing


# import extraction lookup table to obtain covariate classes (`var_class`)
# lookup table is missing "Year" and "Method", so manually adding here
lookup <- readxl::read_xlsx(file.path(root, "NationalModels_V5_VariableList.xlsx"), sheet="ExtractionLookup") |>
  dplyr::select(Category, Label) |> 
  dplyr::rename(var_class = Category, var = Label) |> 
  tibble::add_row(var_class = "Method", var = "method") |> 
  tibble::add_row(var_class = "Year", var = "year")


# create an index containing the species (FLBC), BCR, and bootstrap replicate
sample_id <- 
  gbm_objs_test |> 
  stringr::str_split_fixed(pattern="_", n=3) |> 
  gsub(".R", "", x = _) |>
  dplyr::as_tibble() |> 
  magrittr::set_colnames(c("spp", "bcr", "boot")) |>
  dplyr::arrange(spp, bcr)


# create a list of dataframes containing relative influence per covariate
covs <- list()
for(i in 1:length(gbm_objs_test)){
  
   # loads a `gbm` object named `b.i`
  load(file=file.path(root, "output", "bootstraps", gbm_objs_test[i]))
  
  # `summary(b.i)` is a `data.frame` with columns `var` and `rel.inf`
  covs[[i]] <- b.i |>
    gbm::summary.gbm(plotit = FALSE) |>
    dplyr::left_join(lookup, by="var") |>
    as_tibble() |> 
    dplyr::cross_join(x = _, sample_id[i,])
   
  # print progress
  cat(paste("\riteration", i))
    Sys.sleep(0.5)
}


# flatten list of dataframes
covs_all <- purrr::reduce(covs, full_join)


#WRITE FUNCTION#######################

#' @param species is `"all"` or a `character` with the FLBC denoting the species of interest. 
#' See `unique(sample_id$spp)` for possible species. 
#' 
#' @param method not sure yet...
#'
#' @return a `data.frame` to be passed to a plotting function...tbd
#'
#' @examples ...tbd

# what a user might want:
# rel.inf by species
# rel.inf by bcr


rel_inf_bcr <- function(species=c("all", ...), method){

  # sum relative influence by BCR and variable class
  rel_inf_sum <- 
    covs_all |> 
    group_by(bcr, var_class) |> 
    summarise(sum_influence = sum(rel.inf))
  
  
  # sum rel. influence of variable classes 
  # var_sum <- 
  #   covs_all |> 
  #   group_by(var_class)  |>  
  #   summarise(sum_var_class = sum(rel.inf))
  
  
  # sum rel. influence of BCR
  bcr_sum <-
    covs_all |>  
    group_by(bcr) |>  
    summarise(sum_bcr = sum(rel.inf))
  
  
  bcr_proportion_inf <- 
    rel_inf_sum |> 
    left_join(x = _, bcr_sum, by="bcr") |> 
    mutate(prop = sum_influence/sum_bcr ) 
  
  
  # colour blind palette from: https://colorbrewer2.org/#type=qualitative&scheme=Set3&n=12
  
  cbPalette <- c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5","#ffed6f")
  
  ggplot() +
    geom_bar(aes(x=bcr, y=prop, fill=var_class), data=bcr_proportion_inf, stat = "identity") + 
    scale_fill_manual(values = cbPalette) +
    theme_classic()
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.title =element_blank()) +
    # theme_bw() 
  
  print("plotting proportion of variable influence per BCR")
  
}
                        


# root <- "G:/Shared drives/BAM_NationalModels/NationalModels5.0/output/bootstraps"
# 
# prediction_files <- list.files(root)
# 
# loop <- 
#   prediction_files %>% 
#   stringr::str_split_fixed(pattern="_", n=3) %>% 
#   dplyr::as_tibble() %>% 
#   magrittr::set_colnames(c("spp", "bcr", "boot")) %>% 
#   dplyr::arrange(spp, bcr)
# 
# bcr.i <- loop$bcr[i]
# spp.i <- loop$spp[i]
# boot.i <- loop$boot[i]




