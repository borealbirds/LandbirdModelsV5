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
  dplyr::mutate(boot = as.integer(boot)) |> 
  tidyr::complete(spp, bcr, var, boot=1:32, fill=list(rel.inf=0)) |>  # ensure all bootstraps are present, fill missing with 0
  dplyr::group_by(spp, bcr, var) |>  # Group once by spp, bcr, and var
  dplyr::summarise(mean_rel_inf = mean(rel.inf, na.rm = TRUE),  # calculate mean across all bootstraps (including 0s)
            sd_rel_inf = sd(rel.inf, na.rm = TRUE),      # calculate standard deviation across all bootstraps (including 0s)
            n_boots = sum(rel.inf > 0)) |>    # count the number of bootstraps where the variable appeared (non-zero rel.inf)
  dplyr::filter(mean_rel_inf > 0) |> 
  dplyr::arrange(spp, bcr, desc(mean_rel_inf)) |> 
  dplyr::left_join(varclass_lookup, by = "var") |> 
  dplyr::slice_max(mean_rel_inf, n=20, with_ties = FALSE)

write.csv(bam_covariate_importance, file=file.path(root1, "BAM_RelativeImportance.csv"))
  
  