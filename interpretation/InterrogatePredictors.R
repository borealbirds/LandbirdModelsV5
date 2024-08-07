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
  
  ifelse(bcr == "all", bcr_to_filter <- unique(bam_covariate_importance$bcr), bcr_to_filter <- bcr)
  
  # for dplyr::group_by
  group_sym <- rlang::syms(c(group, "var_class", "boot"))
  
  # need to be able to specify what BCRs (or species or bird group, etc) to plot by
  # sum rel. influence of the grouped variable (e.g. species) per `var_class` and `boot` replicate
  rel_inf_sum <- 
    bam_covariate_importance |> 
    group_by(!!!group_sym) |> 
    dplyr::filter(bcr %in% bcr_to_filter) |> 
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
#'@param data An exported `list` of `data.frame`s (see: `data(boot_pts_sorted)` where list elements are xy coordinates of covariates and predicted responses for a given species x bcr x bootstrap permutation.
#'
#'@param ... This allows users to pass additional arguments to `gbm::plot.gbm()`
#'
#'@return ...tbd
#'
#'@examples ...tbd

# now, for every species x bcr x var permutation, i want to find the average and s.d. x and y grid
library(splines)

bamexplorer_partial_dependence <- function(data = boot_pts_sorted, bcr, common_name, covariate) {
  
  # construct the key for accessing the desired data frame
  key <- paste(bcr, common_name, covariate, sep = "_")
  
  # check if the key exists in the data
  if (!key %in% names(data)) {
    stop("The specified combination of species, bcr, and covariate does not exist.")
  }
  
  # combine all the relevant into a single data frame with a column`replicate` with its bootstrap ID
  combined_df <- bind_rows(data[[key]], .id = "replicate")
  
  # convert the covariate to a symbol for dynamic grouping
  covariate_sym <- rlang::sym(covariate)
  
  # create a vector of x values (covering the domain of the combined bootstraps) for prediction using `smooth.spline()`
  x_grid <- seq(min(combined_df[[covariate]]), max(combined_df[[covariate]]), length.out = 1000)
  
  # fit a smoothing function to each bootstrap replicate and predict over the domain of x values
  predictions <- 
    lapply(data[[key]], function(df) {
    fit <- smooth.spline(df[[covariate]], df$y)
    predict(fit, x_grid)$y 
  })
  
  # combine predictions into a data frame
  prediction_df <- do.call(cbind, predictions)
  
  # simplify column names
  colnames(prediction_df) <- paste0("replicate_", seq_along(predictions))
  prediction_df <- as_tibble(prediction_df)
  prediction_df[[covariate]] <- x_grid
  
  # calculate summary statistics (mean and error bounds) for each x value
  summary_df <- 
    prediction_df |>
    pivot_longer(data = _, cols = starts_with("replicate_"), names_to = "replicate", values_to = "predicted_response") |>
    group_by(.data = _, !!covariate_sym) |>
    summarise(
      mean_response = mean(predicted_response),
      lower_bound = quantile(predicted_response, 0.025),
      upper_bound = quantile(predicted_response, 0.975)
    )
  
  # create a partial dependence plot with error envelope
  ggplot(summary_df, aes(x = !!covariate_sym, y = mean_response)) +
    geom_line(color = "blue") +
    geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), alpha = 0.2, fill = "blue") +
    labs(title = paste("Partial Dependence Plot for", common_name, "in BCR", bcr, "and Covariate", covariate, "across 10 bootstraps"),
         x = covariate, y = "singing males per Ha") +
    theme_minimal()
  
  
}







# 2-way interactions
#'@param data An exported nested `list` (see: `data(boot_pts_sorted_i2)`) where the top level elements are species x bcr x 2-way interaction tuples, 
# and the second-level elements are matrices of the associated average and SD of the bootstrapped prediction spaces
#'
#'@param bcr A `character` specifying the Bird Conservation Regions (BCRs) of interest.
#'
#'@param common_name A `character` specifying the species of interest.
#'
#'@param covariates By default `all` where every covariate interaction (for a given BCR and species) is searched. 
#' Otherwise, the search can be narrowed by specifying a `character` specifying a single covariate (for `n_interactions=2` or `n_interactions=3`) 
#' or two covariates (for `n_interactions=3`) to narrow the search to interactions of interest.
#'
#'@param n_interactions A `numeric` specifying the degree of interactions. Can be `n=2` for two-way interactions or `n=3` for three-way interactions. 
#'
#'@param threshold A `numeric` specifying the number of most influential interactions to output. E.g. `threshold=10` will print the top 5 most influential interactions.
#'
#'@return ...tbd
#'
#'@examples ...tbd
#'
#'



bamexplorer_interactions <- function(data = boot_pts_sorted_i2, bcr, common_name){
  
  
  # construct the keys for accessing the desired gbm objects
  # note that `paste()` can handle vectors so (e.g.) `bcr` or `common_name` may have length > 1.
  keys <- as.character(outer(bcr, common_name, FUN=paste, sep="."))

  # extract bcr x spp info from boot_pts_sorted_i2
  boot_bcr_spp <- 
    paste(stringr::str_split_i(names(boot_pts_sorted_i2), "\\.", 1),
          stringr::str_split_i(names(boot_pts_sorted_i2), "\\.", 2),
          sep=".")
  
  # subset boot_pts_sorted_i2 to the bcr x spp queried
  queried_bcr_spp <- boot_pts_sorted_i2[which(boot_bcr_spp %in% keys)]
  
  # for every element of `queried_bcr_spp` (a single 2-way interaction)
  # 1. get mean of the mean +/- sd 2-way interaction space
  # 2. enter as an element into a list
  # 3. flatten list into dataframe
  # 4. sort by highest y
  
  queried_means_list <- list()
  for (v in 1:length(queried_bcr_spp)){
    
    if (purrr::is_empty(queried_bcr_spp[[v]]) == FALSE){
    
    # bring bcr, spp, var info along
    info_v <- stringr::str_split(names(queried_bcr_spp)[v], "\\.", simplify=TRUE)
    
    queried_means_list[[v]] <- 
      purrr::list_c(queried_bcr_spp[[v]]) |> 
      tibble::as_tibble() |> 
      dplyr::rename(y_mean = V1, y_sd = V2) |> 
      dplyr::summarise(mean_y_mean = mean(y_mean), mean_y_sd = mean(y_sd)) |> 
      dplyr::mutate(bcr = info_v[1], common_name = info_v[2], var_1 = info_v[3], var_2 = info_v[4])
        
        
    # print progress
    cat(paste("\rsearching", v, "of", length(queried_bcr_spp), "interactions (some are NULL)"))
    Sys.sleep(0.000001)
    
    } # close if()
  } # close loop
  
  # gather non-NULL entries
  nulls <- sapply(queried_means_list, is.null)
  
  # gather all interaction space means into a single table and sort
  # `round()` is used for tie-breakers where the means are functionally identical but there is large differences in SD
  queried_means <- 
    queried_means_list[!nulls] |> 
    purrr::list_rbind() |> 
    mutate(mean_y_mean = round(mean_y_mean, 3)) |> 
    dplyr::arrange(desc(mean_y_mean), mean_y_sd)
  
  return(queried_means)
} 
  
# sanity check by plotting
boot_pts_sorted_i2 <- readRDS(file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boot_pts_sorted_i2.rds")
test <- bamexplorer_interactions(bcr = c("can10", "can11"), common_name = "Alder Flycatcher")
load(file=file.path(root, "output", "bootstraps", gbm_objs[1]))


# lowest y_mean (0.016 +/- 0.006)
plot.gbm(x=b.i, return.grid = FALSE, i.var = c("SCANFIclosure_1km", "CanHF_1km"), type="response")


# highest y_mean (0.062 +/- 0.0194)
plot.gbm(x=b.i, return.grid = FALSE, i.var = c("SCANFIBlackSpruce_5x5", "Peatland_1km"), type="response")
 



  
  




#################################################



