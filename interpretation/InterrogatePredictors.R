# ---
# title: National Models 5.0 - covariate interpretation
# author: Mannfred Boehm
# created: May 22, 2024
# ---


#NOTES################################

# This script extracts covariate contributions to model predictions. 

# The outputs will be analysed for identifying covariates of importance across several metrics (e.g. BCR, ecology, etc.)

# The first part of this script generates the exported data for the R package function. 

# The second part is the function that summarises the v.5.0 models. 






#WRITE FUNCTION#######################
#'
#'@param covariate_data An exported `data.frame` (see: `data(bam_covariate_importance)` where rows are covariates and columns denote the relative influence of for a given bootstrap replicate by species by BCR permutation. 
#'
#'@param species The species to be summarised. Default is `"all"` but can also be a `character` with names (common, scientific, or FLBCs) denoting the species of interest. We follow the taxonomy of the AOS 64th Supplement (July 21, 2023). 
#' Perhaps we can work a list of available species into an exported data.frame, e.g. via `usethis::use_data(unique(sample_id$sci_name))`. 
#' 
#'@param bcr The Bird Conservation Regions to be summarised. Works in the same manner as `species`.
#' 
#'@param traits One of `avonet` (Tobias et al 2022) or `acad`. Can also be a `data.frame` with species as rows (see `species` argument above) and traits as columns (can we handle continuous traits?). 
#' 
#'@param group_by `<tidy-select>` An unquoted expression specifying 1-2 grouping variable(s) by which to summarise the relative influence of model covariates. 
#'E.g. `group_by = bcr` summarises the relative influence of model covariates by Bird Conservation Region, while `species` summarises covariates by taxon. 
#'Expressions must be columns that exist in `bam_covariate_importance` or one of the trait databases (see: `traits` argument).  
#' 
#'@param plot If `TRUE`, creates a stacked bar plot with relative influence.
#' 
#'@param colours A `character` of hex codes specifying the colours if `plot = TRUE`. 
#'
#'@param export If `TRUE`, exports to an object the dataframe underlying any plots created by this function.
#' 
#'@return ...tbd
#'
#'@examples ...tbd

# colour blind palette from: https://colorbrewer2.org/#type=qualitative&scheme=Set3&n=12
cbPalette <- c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5","#ffed6f")
groups <- c("bcr", "var_class")

bam_explore(data=covs_all, groups=c("bcr", "var_class"), traits=NULL)


bam_explore <- function(data = covs_all, groups=NULL, species="all", bcr="all", traits = NULL, plot = FALSE, colours=NULL){

  
  # Filter the dataset by species if specified
  user_species <- species
  
  if (species != "all") {
    covs_all <- covs_all |> dplyr::filter(species %in% user_species)
  }
  
  
  # If user specifies built-in dataset, load it via `data()`, otherwise use user-specified dataset
  # if (traits %in% c("avonet", "acad")){
  #   traits <- data(traits)
  # } else {traits <- trait} 
    
  
  # Check that trait data is the right class
  if (is.data.frame(traits) == FALSE) {
    print("trait data must be a `data.frame`")
  }
  
  # get mean and std. deviation across bootstraps for each species x bcr x var_class permutation
  group_sym <- rlang::syms(groups) # for group_by
  
  boot_variance <-
    covs_all |> 
    group_by(!!!group_sym, sci_name) |> 
    summarise(mean_boot = mean(rel.inf), sd_boot = sd(rel.inf))
  
  
  # summarise covariate importance across user-specified groups
  rel_inf_sum <- 
    covs_all |> 
    group_by(!!!group_sym) |> 
    summarise(sum_influence = sum(rel.inf))

  
  # calculate the sum of rel. influence of group 1
  var_sum1 <-
    rel_inf_sum |>
    group_by(!!group_sym[[1]])  |>
    summarise(sum_group1 = sum(sum_influence))
  
  
  # calculate the sum of rel. influence of group 2
  # var_sum2 <-
  #   rel_inf_sum |>  
  #   group_by(!!group_sym[[2]]) |>  
  #   summarise(sum_group2 = sum(sum_influence))
  
  
  proportion_inf <- 
    rel_inf_sum |> 
    left_join(x = _, var_sum1, by=groups[1]) |> 
    mutate(prop = sum_influence/sum_group1) 
  
  
  ggplot() +
    geom_bar(aes(x=groups[1], y=prop, fill=var_class), data=proportion_inf, stat = "identity") + 
    scale_fill_manual(values = colours) +
    theme_classic()
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.title =element_blank()) +
    # theme_bw() 
  
  print(paste("plotting proportion of variable influence per"), grouping_vars, sep=" ")
  
}
                        
