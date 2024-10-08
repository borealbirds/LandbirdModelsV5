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
  
  # determine if the variable is a factor or continuous
  if (is.factor(combined_df[[covariate]])) {
    # For factors, use the levels as the x_grid
    x_grid <- levels(combined_df[[covariate]])
    
    # calculate the mean response for each level of the factor
    summary_df <- combined_df %>%
      group_by(!!rlang::sym(covariate)) %>%
      summarise(
        mean_response = mean(y, na.rm = TRUE),
        lower_bound = quantile(y, 0.025, na.rm = TRUE),
        upper_bound = quantile(y, 0.975, na.rm = TRUE)
      )
    
    # create the partial dependence plot for factor variables
    ggplot(summary_df, aes_string(x = covariate, y = "mean_response")) +
      geom_point(size = 3, color = "blue") +
      geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), width = 0.2, color = "blue") +
      labs(title = paste("Partial Dependence Plot for", common_name, "in BCR", bcr, "and Covariate", covariate, "across 10 bootstraps"),
           x = covariate, y = "Singing Males per Ha") +
      theme_minimal()
    
  } else {
    
    # for continuous variables, create a vector of x values (covering the domain of the combined bootstraps) for 
    # generating predictions via `smooth.spline()`
    x_grid <- seq(min(combined_df[[covariate]], na.rm = TRUE), 
                  max(combined_df[[covariate]], na.rm = TRUE), 
                  length.out = 1000)
    
    # fit a smoothing function to each bootstrap replicate and predict over the domain of x values
    predictions <- lapply(data[[key]], function(df) {
      fit <- smooth.spline(df[[covariate]], df$y)
      predict(fit, x_grid)$y
    })
    
    # combine predictions into a data frame
    prediction_df <- as_tibble(do.call(cbind, predictions))
    colnames(prediction_df) <- paste0("replicate_", seq_along(predictions))
    prediction_df[[covariate]] <- x_grid
    
    # Calculate summary statistics (mean and error bounds) for each x value
    summary_df <- prediction_df %>%
      pivot_longer(cols = starts_with("replicate_"), names_to = "replicate", values_to = "predicted_response") %>%
      group_by(!!rlang::sym(covariate)) %>%
      summarise(
        mean_response = mean(predicted_response, na.rm = TRUE),
        lower_bound = quantile(predicted_response, 0.025, na.rm = TRUE),
        upper_bound = quantile(predicted_response, 0.975, na.rm = TRUE)
      )
    
    # Create the partial dependence plot for continuous variables
    ggplot(summary_df, aes(x = !!rlang::sym(covariate), y = mean_response)) +
      geom_line(color = "blue") +
      geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), alpha = 0.2, fill = "blue") +
      labs(title = paste("Partial Dependence Plot for", common_name, "in BCR", bcr, "and Covariate", covariate, "across 10 bootstraps"),
           x = covariate, y = "Singing Males per Ha") +
      theme_minimal()
  }# close ifelse
} # close function
  
  




# faceted partial dependence plot
#'@param data An exported `list` of `data.frame`s (see: `data(boot_pts_sorted)` where list elements are xy coordinates of covariates and predicted responses for a given species x bcr x bootstrap permutation.
#'
#'@param ... This allows users to pass additional arguments to `gbm::plot.gbm()`
#'
#'@return ...tbd
#'
#'@examples ...tbd


bamexplorer_facet_pd <- function(data = boot_pts_sorted, var_type=NULL, bcr, common_name, n = 9) {
  
  # ensure the user specifies the variable type
  if (is.null(var_type) || !var_type %in% c("continuous", "factor")) {
    stop("`var_type` must be 'continuous' or 'factor'.")
  }
  
  # filter the top `n` covariates based on relative importance
  # `bam_covariate_importance` needs to be called from backend
  # NOTE: data not yet filtered for `var_type`
  ranked_covariates <- 
    bam_covariate_importance |> 
    dplyr::filter(common_name == common_name & bcr == bcr) |> 
    dplyr::group_by(var) |>  # group by variable name
    dplyr::summarise(mean_rel_inf = mean(rel.inf, na.rm = TRUE)) |>  # calculate the mean relative importance
    dplyr::arrange(desc(mean_rel_inf))  # sort by mean relative importance
  
  
  # iterate over each row in covariate_importance to filter for variable type
  specified_covariates <- c()
  for (k in seq_len(nrow(ranked_covariates))) {
    var <- ranked_covariates$var[k]
    key <- paste(bcr, common_name, var, sep = "_")
    
    # check if the key exists in the data
    if (!key %in% names(data)) {
      next  # skip if the key does not exist
    } #close if()
    
    # combine all bootstrap replicates into one data frame
    combined_df <- bind_rows(data[[key]], .id = "replicate")
    
    # if var[i] matches the user-specified type, add it to the top `n` list
    if (var_type == "factor" && is.factor(combined_df[[var]]) || 
         (var_type == "continuous" && !is.factor(combined_df[[var]]))) {
    
        specified_covariates <- c(specified_covariates, var)
    } #close if()
    
    # stop when we have `n` covariates
    if (length(specified_covariates) == n) {
      break
    } #close if()
    
  } #close for()
 
  
  # give a warning if fewer than `n` covariates of the specified type are available
  if (length(specified_covariates) < n) {
    warning(paste("Only", length(specified_covariates), "matching covariates found for the specified type."))
  }
  

  # this loop creates a dataframe of density predictions for every covariate 
  # in `specified_covarites`. The result is a list with as many elements as
  # `length(specified_covariates)`
  specified_covariates_dfs <- list()
  
  for (v in 1:length(specified_covariates)) {
    
    var <- specified_covariates[v]
    
    # construct the key to access the data frame
    key <- paste(bcr, common_name, var, sep = "_")
    
    # combine all bootstrap replicates of a specific spp x bcr x var tuple
    # into one data frame with a `replicate` column
    combined_df <- bind_rows(data[[key]], .id = "replicate")
    covariate_sym <- rlang::sym(var)
    
    if (var_type == "factor") {
      # For factors, calculate mean response for each level
      summary_df <- combined_df %>%
        group_by(!!covariate_sym) %>%
        summarise(
          mean_response = mean(y, na.rm = TRUE),
          lower_bound = quantile(y, 0.025, na.rm = TRUE),
          upper_bound = quantile(y, 0.975, na.rm = TRUE)
        ) %>%
        mutate(covariate_value = !!covariate_sym,
               covariate_name = var) %>%
        select(covariate_value, mean_y, lower_y, upper_y, covariate_name)  # add a covariate identifier for faceting
      
    } else if (var_type == "continuous") {
      
      # for continuous variables, define the prediction domain
      x_grid <- seq(min(combined_df[[var]], na.rm = TRUE), 
                    max(combined_df[[var]], na.rm = TRUE), 
                    length.out = 1000)
      
      # fit smoothing splines and predict for each bootstrap replicate
      f <- function(df) {fit <- smooth.spline(df[[var]], df$y); predict(fit, x_grid)$y}
      predictions_v <- lapply(data[[key]], f)
      
      # combine predictions into a single data frame
      predictions_v_merged <- as_tibble(do.call(cbind, predictions_v))
      colnames(predictions_v_merged) <- paste0("replicate_", seq_along(predictions_v_merged))
      predictions_v_merged[[var]] <- 
        predictions_v_merged |> 
        dplyr::mutate(covariate_value = x_grid)
      
      # estimate mean and error bounds for each predictor value
      summary_df <- 
        predictions_v_merged[[var]] %>%
        pivot_longer(cols = starts_with("replicate_"), names_to = "replicate", values_to = "predicted_response") %>%
        group_by(covariate_value) %>%
        summarise(
          mean_y = mean(predicted_response, na.rm = TRUE),
          lower_y = quantile(predicted_response, 0.025, na.rm = TRUE),
          upper_y = quantile(predicted_response, 0.975, na.rm = TRUE)
        ) %>%
        mutate(covariate_name = var)   # add a covariate identifier for faceting
    } #close else()
    
    # create a list element for the target covariate
    specified_covariates_dfs[[var]] <- summary_df
    
  } # close for()

  # merge all data frames of predictions 
  merged_df <- dplyr::bind_rows(specified_covariates_dfs)
  
  
  # generate a faceted plot based on the specified variable type
  plot <- 
    ggplot(merged_df, aes(covariate_value, y = mean_y)) +
    labs(title = paste("Partial Dependence Plots for", common_name, "in BCR", bcr),
         x = "Covariate Value", y = "Singing Males per Ha") +
    facet_wrap(~ covariate_name, scales = "free_x") +
    theme_minimal()
  
  if (var_type == "continuous") {
    plot <- plot +
      geom_line(color = "purple") +
      geom_ribbon(aes(ymin = lower_y, ymax = upper_y), alpha = 0.2, fill = "purple")
  } else if (var_type == "factor") {
    plot <- plot +
      geom_point(size = 3, color = "purple") +
      geom_errorbar(aes(ymin = lower_y, ymax = upper_y), width = 0.2, color = "purple")
  } # close if()
  
  print(plot)
} #close function
  





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



bamexplorer_interactions <- function(data = boot_pts_reduced_i2, bcr, common_name){
  
  # construct the keys for accessing the desired gbm objects
  # note that `paste()` can handle vectors so (e.g.) `bcr` or `common_name` may have length > 1.
  keys <- as.character(outer(bcr, common_name, FUN=paste, sep="."))

  # extract bcr x spp info from boot_pts_sorted_i2
  boot_bcr_spp <- 
    paste(stringr::str_split_i(names(boot_pts_reduced_i2), "\\.", 1),
          stringr::str_split_i(names(boot_pts_reduced_i2), "\\.", 2),
          sep=".")
  
  # subset boot_pts_sorted_i2 to the bcr x spp queried
  queried_bcr_spp <- boot_pts_reduced_i2[which(boot_bcr_spp %in% keys)]
  
  # for every element of `queried_bcr_spp` (a single 2-way interaction)
  # 1. get mean of the mean +/- sd 2-way interaction space
  # 2. enter as an element into a list
  # 3. flatten list into dataframe
  # 4. sort by highest y
  
  queried_means_list <- list()
  for (v in 1:length(queried_bcr_spp)){
    
    # bring bcr, spp, var info along
    info_v <- stringr::str_split(names(queried_bcr_spp)[v], "\\.", simplify=TRUE)
    
    queried_means_list[[v]] <- 
      queried_bcr_spp[[v]] |> 
      tibble::as_tibble() |> 
      dplyr::mutate(bcr = info_v[1], common_name = info_v[2], var_1 = info_v[3], var_2 = info_v[4])
        
        
    # print progress
    cat(paste("\rsearching", v, "of", length(queried_bcr_spp), "2-way interactions"))
    Sys.sleep(0.000001)
    
  }
  
  # gather all interaction space means into a single table and sort
  # `round()` is used for tie-breakers where the means are functionally identical but there is large differences in SD
  queried_means <- 
    queried_means_list |> 
    purrr::list_rbind() |> 
    mutate(mean_y_mean = round(mean_y_mean, 3)) |> 
    dplyr::arrange(desc(mean_y_mean), mean_y_sd)
  
  return(queried_means)
} 
  
# sanity check by plotting
boot_pts_reduced_i2 <- readRDS(file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boot_pts_reduced_i2.rds")
test <- bamexplorer_interactions(bcr = c("can12", "can11"), common_name = "Alder Flycatcher")
load(file=file.path(root, "output", "bootstraps", gbm_objs[1]))


# lowest y_mean (0.016 +/- 0.006)
plot.gbm(x=b.i, return.grid = FALSE, i.var = c("SCANFITamarack_5x5", "year"), type="response")


# highest y_mean (0.062 +/- 0.0194)
plot.gbm(x=b.i, return.grid = FALSE, i.var = c("SCANFIJackPine_5x5", "SCANFITamarack_1km"), type="response")
 



# PROBLEM: if threshold is 0.10, what do we do when bootstraps vary between e.g. 0.08-0.011? The lower bootstraps will be lost




#################################################



