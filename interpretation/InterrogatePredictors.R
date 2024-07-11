# ---
# title: National Models 5.0 - covariate interpretation
# author: Mannfred Boehm
# created: May 22, 2024
# ---

#'@param covariate_data An exported `data.frame` (see: `data(bam_covariate_importance)` where rows are covariates and columns denote the relative influence of for a given bootstrap replicate by species by BCR permutation. 
#'
#'@param species The species to be summarised. Default is `"all"` but can also be a `character` with names (common, scientific, or FLBCs) denoting the species of interest. We follow the taxonomy of the AOS 64th Supplement (July 21, 2023). 
#' Perhaps we can work a list of available species into an exported data.frame, e.g. via `usethis::use_data(unique(sample_id$sci_name))`. 
#' 
#'@param bcr The Bird Conservation Regions to be summarised. Works in the same manner as `species`.
#' 
#'@param traits One of `avonet` (Tobias et al 2022) or `acad`. Can also be a `data.frame` with species as rows (see `species` argument above) and traits as columns (can we handle continuous traits?). 
#' 
#'@param groups A `character` specifying two grouping variables by which to summarise the relative influence of model covariates. 
#'E.g. `groups = c("bcr", "var_class)"` summarises the relative influence of model covariates by Bird Conservation Region and variable class (e.g. Landcover, Biomass, Climate).
#'Note that the first grouping variable must be discrete and will constitute the x-axis of a stacked barplot, while the relative influence will be estimated for the second grouping variable.  
#'Expressions must be columns that exist in `bam_covariate_importance` or one of the trait databases (see: `traits` argument).  
#' 
#'@param plot If `TRUE`, creates a stacked bar plot with relative influence.
#' 
#'@param colours A `character` of hex codes specifying the colours if `plot = TRUE`. 
#'
#'@param export If `TRUE`, exports to an object the dataframe underlying any plots created by this function.
#' 
#'@return A stacked barchart. The y-axis is the proportion of covariate importance.
#'
#'@examples ...tbd



bamexplorer_stackedbarchart <- function(data = bam_covariate_importance, groups = NULL, species = "all", bcr = "all", traits = NULL, plot = FALSE, colours = NULL){

  
  # Filter the dataset by species if specified
  user_species <- species
  
  if (user_species != "all" & 
      (all(user_species %in% data$common_name)==TRUE | all(user_species %in% data$sci_name)==TRUE)) {
    bam_covariate_importancel <- 
      bam_covariate_importance |> 
      dplyr::filter(species %in% user_species)
  }
  
  
  # If user specifies built-in dataset, load it via `data()`, otherwise use user-specified dataset
  # if (traits %in% c("avonet", "acad")){
  #   traits <- data(traits)
  # } else {traits <- trait} 
    
  
  # Check that trait data is the right class
  if (is.null(traits) == FALSE & is.data.frame(traits) == FALSE) {
    print("argument `traits` must be NULL or a `data.frame`")
  }
  
  # for dplyr::group_by
  group_sym <- rlang::syms(groups) 
  
  
  # sum covariate importance across for every permutation of group1 and group2
  rel_inf_sum <- 
    bam_covariate_importance |> 
    group_by(!!!group_sym) |> 
    summarise(sum_influence = sum(rel.inf), .groups="keep")

  
  # sum of covariate importance for each of group1 (all group2 sums are amalgamated into group1 bins)
  group1_sum <-
    rel_inf_sum |>
    group_by(!!group_sym[[1]])  |>
    summarise(sum_group1 = sum(sum_influence), .groups="keep")
  
  # get the %contribution of group2 covariates to overall covariate importance for a given group1
  proportion_inf <- 
    rel_inf_sum |> 
    left_join(x = _, group1_sum, by=groups[1]) |> 
    mutate(prop = sum_influence/sum_group1) 
  
  
  barplot <- 
    ggplot() +
    geom_bar(aes(x=!!group_sym[[1]], y=prop, fill=!!group_sym[[2]]), data=proportion_inf, stat = "identity") + 
    # scale_fill_manual(values = colours) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
  
  print(barplot)
  print(paste("plotting proportion of variable influence by", groups[1], "and", groups[2], sep=" "))
  
}

# test function
# colour blind palette from: https://colorbrewer2.org/#type=qualitative&scheme=Set3&n=12
# colours <- c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5","#ffed6f")
bamexplorer_stackedbarchart(groups=c("common_name", "var_class"))



  











# side-by-side boxplots displaying bootstrap variation by variable class
#'@param covariate_data An exported `data.frame` (see: `data(bam_covariate_importance)` where rows are covariates and columns denote the relative influence of for a given bootstrap replicate by species by BCR permutation. 
#'
#'
#'@return ...tbd
#'
#'@examples ...tbd


bamexplorer_boxplots <- function(data = bam_covariate_importance, group = NULL, species = "all", bcr = "all", traits = NULL, plot = FALSE, colours = NULL){
  
  
  # for dplyr::group_by
  group_sym <- rlang::syms(c(group, "var_class", "boot"))
  
  # need to be able to specify what BCRs (or species or bird group, etc) to plot by
  # sum rel. influence of the grouped variable (e.g. species) per `var_class` and `boot` replicate
  rel_inf_sum <- 
    bam_covariate_importance |> 
    group_by(!!!group_sym) |> 
    summarise(sum_influence = sum(rel.inf), .groups="keep")
  
  
  # sum of covariate importance for each of group1 (all var_class sums are amalgamated into group1 bins)
  group1_sum <-
    rel_inf_sum |>
    group_by(!!group_sym[[1]], !!group_sym[[2]])  |>
    summarise(sum_group1 = sum(sum_influence), .groups="keep")
  
  
  proportion_inf <- 
    rel_inf_sum |> 
    left_join(x = _, group1_sum, by=c(group, "var_class")) |> 
    mutate(prop = sum_influence/sum_group1) 
  
  
  ggplot(proportion_inf, aes(x = var_class, y = prop, fill = !!group_sym[[1]])) +
    geom_boxplot(position = position_dodge(width = 0.75), alpha=0.05) +
    geom_point(aes(colour=factor(!!group_sym[[1]])),  position = position_dodge(width = 0.75), alpha=0.7, size=2.5) +
    labs(x = "Variable Class", y = "Relative Importance (%)", 
         title = paste("Covariate importance by", group, sep=" ")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
}

# test function
bamexplorer_boxplots(group="common_name")










# partial dependence plot
#'@param covariate_data An exported `data.frame` (see: `data(bam_covariate_importance)` where rows are covariates and columns denote the relative influence of for a given bootstrap replicate by species by BCR permutation. 
#'
#'@param ... This allows users to pass additional arguments to `gbm::plot.gbm()`
#'
#'@return ...tbd
#'
#'@examples ...tbd


# this loop creates a list of lists of `data.frame`s. 
# each top level element is a bootstrap replicate, each second-level element represents a covariate for that bootstrap.
# each covariate has a corresponding `data.frame`, with column 1 being the covariate domain and column 2 being the predicted range.
boot_pts <- list()
for (q in 1:length(gbm_objs)){
  
  # load a bootstrap replicate 
  load(file.path(root, "output", "bootstraps", gbm_objs[q]))
  
  # find evaluation points (`data.frame`) for every covariate (indexed by `i`)  
  pts <- list()
  for (i in 1:length(b.i$var.names)){
  pts[[i]] <- plot.gbm(x=b.i, return.grid = TRUE, i.var = i)
  }
  
  boot_pts[[q]] <- pts
  
  # print progress
  cat(paste("\riteration", q))
  Sys.sleep(0.001)
}

# now, for every species x bcr x var permutation, i want to find the average and s.d. x and y grid

# group `b.i` objects by bcr x species x var (or var_class?)
# group_keys() finds the names of the group permutations
# 1617 permutations of  bcr x common_name x var 
boot_group_keys <- 
  bam_covariate_importance |> 
  dplyr::group_by(bcr, common_name, var) |> 
  dplyr::group_keys()


# go into boot_group_rows[i] and find row numbers for the ith species x bcr x var permutation
# the name of the sbv permutation is boot_group_keys[i,]
# use boot_group_keys[i,]$var to get the variable name; use that name to search the column names of boot_pts[[i]]
# when the column 1 name matches boot_group_keys[i], then assign it to newlist[[i]]


# output: a list of lists where the top level elements are species x bcr x var permutations, 
# and the second-level elements are dataframes of the associated bootstrapped model predictions

boot_pts_sorted <- list()
for (z in 1:length(boot_group_rows)){
  
  key_z <- boot_group_keys[z,]
  
  # `sample_id` and `boot_pts` have the same length and order
  # so we can annotate elements in `boot_pts` using info from `sample_id`
  boot_pts_index <- which(sample_id$bcr == key_z$bcr & sample_id$common_name == key_z$common_name)
  
  # a list of bootstrap predictions for zth species x bcr permutation 
  spp_bcr_list <- boot_pts[boot_pts_index]
  
  # search for all bootstrap predictions for the zth covariate
  spp_bcr_var <- list()
  for (w in 1:length(spp_bcr_list)) {
    
    # get covariate names for the current spp x bcr permutation 
    var_names <- 
      lapply(spp_bcr_list[[w]], colnames) |> 
      lapply(X=_, `[[`, 1) |> #don't need the name of the y variable
      purrr::flatten_chr()
    
    # find which elements match the current covariate of interest
    spp_bcr_var[[w]] <- spp_bcr_list[[w]][which(var_names %in% key_z)] 
  }
  
  boot_pts_sorted[[z]] <- spp_bcr_var
  names(boot_pts_sorted[[z]]) <- toString(boot_group_keys[z,])
  
  # print progress
  cat(paste("\riteration", z))
  Sys.sleep(0.001)
}

 




