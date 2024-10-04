library(gbm)
library(tidyverse)


# import extraction lookup table to obtain covariate classes (`var_class`)
# lookup table is missing "Year" and "Method", so manually adding here
root1 <- "G:/Shared drives/BAM_NationalModels4/NationalModels4.0/Manuscript/variable importance/"
varclass_lookup <- read_csv(file.path(root1, "BAM_VariableListClasses.csv")) 
  

# import covariate importance data from "01_prepare_data.R"
root2 <- "C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/NationalModelsV5/"
bam_covariate_importance <- 
  readRDS(file.path(root2, "v4", "covariate_interactions", "covariate_importance.rds")) |> 
  dplyr::filter(rel.inf >= 1.00) |> 
  dplyr::group_by(spp, bcr, var) |> 
  dplyr::summarise(
    mean_rel_inf = mean(rel.inf, na.rm = TRUE),  
    sd_rel_inf = sd(rel.inf, na.rm = TRUE),
    n_boots = n()) |> 
  dplyr::filter(n_boots >= 5) |> 
  dplyr::arrange(spp, bcr, desc(mean_rel_inf)) |> 
  dplyr::left_join(varclass_lookup, by="var") 

write.csv(bam_covariate_importance, file=file.path(root1, "BAM_RelativeImportance.csv"))
  
  