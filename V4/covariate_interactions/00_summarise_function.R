# this function is called `bamexplorer_interactions` in BAMexploreR
# as of Oct 7, 2024 


summarise_interactions <- function(data = boot_pts_reduced_i2, bcr, common_name){
  
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
