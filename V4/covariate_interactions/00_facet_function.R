# this function is called `bamexplorer_facet_pd` in BAMexploreR
# as of Oct 7, 2024

# need to indicate scales= "free", fixed", free_x", "free_y" (see: https://ggplot2.tidyverse.org/reference/facet_grid.html)
# covariate_importance_zeroed <- readRDS("rds_files/02_covariate_importance_zeroed.rds")
# boot_pts_sorted <- readRDS("rds_files/02_boot_pts_sorted.rds")

bamexplorer_facet_pd <- function(data = boot_pts_sorted, var_type=NULL, bcr, spp, n = 9, scales=NULL) {
  
  # ensure the user specifies the variable type
  if (is.null(var_type) || !var_type %in% c("continuous", "factor")) {
    stop("`var_type` must be 'continuous' or 'factor'.")
  }
  
  # filter the top `n` covariates based on relative importance
  # `covariate_importance_zeroed` needs to be called from backend
  # NOTE: data not yet filtered for `var_type`
  ranked_covariates <- 
    covariate_importance_zeroed |> 
    dplyr::filter(spp == !!spp & bcr == !!bcr) |>  # in case `spp` or `bcr` are objects in the local env.
    slice_max(mean_rel_inf, n = n, with_ties = FALSE)
  
  
  # iterate over each row in covariate_importance to filter for variable type
  specified_covariates <- c()
  for (k in seq_len(nrow(ranked_covariates))) {
    var <- ranked_covariates$var[k]
    key <- paste(bcr, spp, var, sep = "_")
    
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
    key <- paste(bcr, spp, var, sep = "_")
    
    # combine all bootstrap replicates of a specific spp x bcr x var tuple
    # into one data frame with a `replicate` column
    combined_df <- dplyr::bind_rows(data[[key]], .id = "replicate")
    covariate_sym <- rlang::sym(var)
    
    if (var_type == "factor") {
      # For factors, calculate mean response for each level
      summary_df <- combined_df |> 
        dplyr::group_by(!!covariate_sym) |> 
        dplyr::summarise(
          mean_response = mean(y, na.rm = TRUE),
          lower_bound = quantile(y, 0.025, na.rm = TRUE),
          upper_bound = quantile(y, 0.975, na.rm = TRUE)
        ) |> 
        dplyr::mutate(covariate_value = !!covariate_sym,
               covariate_name = var) |> 
        dplyr::select(covariate_value, mean_y, lower_y, upper_y, covariate_name)  # add a covariate identifier for faceting
      
    } else if (var_type == "continuous") {
      
      # for continuous variables, define the prediction domain
      x_grid <- seq(min(combined_df[[var]], na.rm = TRUE), 
                    max(combined_df[[var]], na.rm = TRUE), 
                    length.out = 1000)
      
      # fit smoothing splines and predict for each bootstrap replicate
      f <- function(df) {fit <- smooth.spline(df[[var]], df$y); predict(fit, x_grid)$y}
      predictions_v <- lapply(data[[key]], f)
      
      # combine predictions into a single data frame
      predictions_v_merged <- tibble::as_tibble(do.call(cbind, predictions_v))
      colnames(predictions_v_merged) <- paste0("replicate_", seq_along(predictions_v_merged))
      predictions_v_merged[[var]] <- 
        predictions_v_merged |> 
        dplyr::mutate(covariate_value = x_grid)
      
      # estimate mean and error bounds for each predictor value
      summary_df <- 
        predictions_v_merged[[var]] |> 
        tidyr::pivot_longer(cols = starts_with("replicate_"), names_to = "replicate", values_to = "predicted_response") |> 
        dplyr::group_by(covariate_value) |> 
        dplyr::summarise(
          mean_y = mean(predicted_response, na.rm = TRUE),
          lower_y = quantile(predicted_response, 0.025, na.rm = TRUE),
          upper_y = quantile(predicted_response, 0.975, na.rm = TRUE)
        ) |> 
        dplyr::mutate(covariate_name = var)   # add a covariate identifier for faceting
    } #close else()
    
    # create a list element for the target covariate
    specified_covariates_dfs[[var]] <- summary_df
    
  } # close for()
  
  
  # extract the order of covariates from ranked_covariates
  # this helps `ggplot` to create facets from most influential to least
  covariate_order <- 
    ranked_covariates |> 
    dplyr::filter(var %in% specified_covariates)  |> 
    dplyr::pull(var) |> 
    unique()
  
  # merge all data frames of predictions 
  # treat `covariate_name` as a factor to prevent alphabetizing by `ggplot`
   merged_df <- 
     specified_covariates_dfs |> 
     dplyr::bind_rows() |> 
     dplyr::mutate(covariate_name = factor(covariate_name, levels = covariate_order)) 
  
  
  # generate a faceted plot based on the specified variable type
  plot <- 
    ggplot(merged_df, aes(covariate_value, y = mean_y)) +
    labs(title = paste("Partial Dependence Plots for", spp, "in", bcr),
         x = "Covariate Value", y = "Singing Males per Ha") +
    facet_wrap(~ covariate_name, scales = scales) +
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
