ExDet<-function (ref, tg, xp) 
{
  tg <- as.matrix(tg)
  ref <- as.matrix(ref)
  a <- apply(ref, 2, min, na.rm = TRUE)
  b <- apply(ref, 2, max, na.rm = TRUE)
  minref <- matrix(a, nrow = nrow(tg), ncol = ncol(tg), byrow = TRUE)
  maxref <- matrix(b, nrow = nrow(tg), ncol = ncol(tg), byrow = TRUE)
  nt1.df <- data.frame(apply(array(data = c(tg - minref, maxref - 
                                              tg, rep(0, nrow(tg) * ncol(tg))), dim = c(dim(tg), 3)), 
                             c(1, 2), min)/(maxref - minref))
  names(nt1.df) <- xp
  nt1 <- rowSums(nt1.df)
  mic_nt1 <- apply(nt1.df, 1, FUN = function(x) base::which.min(x))
  univ.rge <- which(nt1 == 0)
  mic_nt1[univ.rge] <- NA
  tg.univ <- matrix(tg[univ.rge, ], ncol = ncol(tg))
  aa <- apply(ref, 2, mean, na.rm = TRUE)
  bb <- stats::var(ref, na.rm = TRUE)
  mah.ref <- stats::mahalanobis(x = ref, center = aa, cov = bb, tol=1e-20) # increase tolerance
  mah.pro <- stats::mahalanobis(x = tg.univ, center = aa, cov = bb, tol=1e-20)# increase tolerance
  mah.max <- max(mah.ref[is.finite(mah.ref)])
  nt2 <- mah.pro/mah.max
  nt1[univ.rge] <- nt2
  if (length(xp) == 1) 
    cov.combs <- matrix(1)
  if (length(xp) > 1) 
    cov.combs <- utils::combn(x = 1:ncol(tg.univ), m = length(xp) - 
                                1)
  cov.combs <- as.list(data.frame(cov.combs))
  if (length(xp) == 1) {
    cov.aa <- cov.combs %>% purrr::map(., ~apply(as.matrix(ref[, 
                                                               .]), 2, mean))
    cov.bb <- cov.combs %>% purrr::map(., ~var(as.matrix(ref[, 
                                                             .])))
  }
  else {
    cov.aa <- cov.combs %>% purrr::map(., ~apply(as.matrix(ref[, 
                                                               .]), 2, mean, na.rm = TRUE))
    cov.bb <- cov.combs %>% purrr::map(., ~var(as.matrix(ref[, 
                                                             .]), na.rm = TRUE))
  }
  if (nrow(tg.univ) < 2) {
    warning("Only one prediction point within analogue conditions. Mahalanobis distances cannot be calculated.")
    mah_nt2 <- vector(mode = "list", length = length(cov.combs))
  }
  else {
    mah_nt2 <- purrr::pmap(.l = list(cov.combs, cov.aa, cov.bb), 
                           .f = function(a, b, c) stats::mahalanobis(x = as.matrix(tg.univ[, 
                                                                                           a]), center = b, cov = c, tol=1e-20))# increase tolerance
  }
  mah_nt2 <- mah_nt2 %>% purrr::set_names(., xp)
  mah_nt2 <- mah_nt2 %>% purrr::map_df(., cbind)
  mah_nt2 <- as.matrix(mah_nt2)
  mic_nt2 <- 100 * (mah.pro - mah_nt2)/mah_nt2
  mic_nt2 <- apply(mic_nt2, 1, FUN = function(x) base::which.max(x))
  results <- tibble::tibble(ExDet = nt1, mic_univariate = mic_nt1, 
                            mic_combinatorial = NA)
  if (nrow(tg.univ) > 1) 
    results$mic_combinatorial[univ.rge] <- mic_nt2
  results <- results %>% dplyr::mutate(mic_combinatorial = ifelse(ExDet >= 
                                                                    0 & ExDet <= 1, NA, mic_combinatorial))
  results <- results %>% dplyr::mutate(mic = rowSums(.[2:3], 
                                                     na.rm = TRUE))
  return(results)
}


compute_extrapolation <- function(samples,
                                  segments,
                                  covariate.names,
                                  prediction.grid,
                                  coordinate.system,
                                  resolution = NULL,
                                  verbose = TRUE){
  
  #---------------------------------------------
  # Perform function checks
  #---------------------------------------------
  
  calls <- names(sapply(match.call(), deparse))[-1]
  
  if(any("segments" %in% calls)) {
    if(verbose) warning("The 'segments' argument is deprecated, please use 'samples' instead.")
    samples <- segments
  }
  
  if(!"data.frame"%in%class(prediction.grid)) stop("pred.grid must be of class data.frame")
  if(!"data.frame"%in%class(samples)) stop("samples must be of class data.frame")
  
  if(!"x"%in%names(prediction.grid) | !"y"%in%names(prediction.grid)) stop("pred.grid must contain x and y coordinates")
  
  if(!all(covariate.names%in%names(samples))) stop("Missing/unrecognised covariates in the sample data")
  if(!all(covariate.names%in%names(prediction.grid))) stop("Missing/unrecognised covariates in the prediction grid")
  
  coordinate.system <- coordinate.system
  
  samples <- na.omit(samples) # Cannot have NA values
  prediction.grid <- na.omit(prediction.grid)
  
  #---------------------------------------------
  # Check if prediction grid is regular
  #---------------------------------------------
  
  check.grid <- prediction.grid %>%
    dplyr::select(x, y) %>%
    dplyr::mutate(z = 1)
  
  grid.regular <- try(raster::rasterFromXYZ(check.grid), silent = TRUE)
  
  # grid.regular <- try(raster::rasterFromXYZ(check.grid,
  #                     res = ifelse(is.null(resolution), c(NA,NA), resolution)),
  #                     silent = TRUE)
  
  #---------------------------------------------
  # If grid is irregular, rasterise prediction.grid based on specified resolution
  #---------------------------------------------
  
  if (class(grid.regular) == "try-error") {
    if (is.null(resolution)) stop("Prediction grid cells are not regularly spaced.\nA target raster resolution must be specified. See package documentation for details.")
    
    if (verbose) warning("Prediction grid cells are not regularly spaced.\nData will be rasterised and covariate values averaged. See package documentation for details.")
    
    RasteriseGrid <- TRUE
  } else if (class(grid.regular) == "RasterLayer" & !is.null(resolution)) {
    if (verbose) warning("New resolution specified.\nData will be rasterised and covariate values averaged. See package documentation for details.")
    
    RasteriseGrid <- TRUE
    
  } else {
    RasteriseGrid <- FALSE
  }
  
  if(RasteriseGrid){
    
    check.grid$z <- NULL
    sp::coordinates(check.grid) <- ~x+y
    sp::proj4string(check.grid) <- coordinate.system
    
    # Create empty raster with desired resolution
    
    ras <- raster::raster(raster::extent(check.grid), res = resolution)
    raster::crs(ras) <- coordinate.system
    
    # Create individual rasters for each covariate
    
    ras.list <- purrr::map(.x = covariate.names,
                           .f = ~raster::rasterize(as.data.frame(check.grid), ras,
                                                   prediction.grid[,.x], fun = mean_ras)) %>%
      purrr::set_names(., covariate.names)
    
    # Combine all rasters
    
    ras.list <- raster::stack(ras.list)
    
    # Update prediction grid
    
    prediction.grid <- raster::as.data.frame(ras.list, xy = TRUE, na.rm = TRUE)
    
    
  } # End if
  
  if(verbose) message("Computing ...")
  
  #---------------------------------------------
  # Define reference and target systems
  #---------------------------------------------
  
  reference <- samples[, covariate.names]
  target <- prediction.grid[, covariate.names]
  
  #---------------------------------------------
  # Run the exdet tool from Mesgaran et al. (2014)
  #---------------------------------------------
  
  mesgaran <- ExDet(ref = reference,
                    tg = target,
                    xp = covariate.names)
  
  #---------------------------------------------
  # Add coordinates
  #---------------------------------------------
  
  mesgaran <- prediction.grid %>%
    dplyr::select(x,y) %>%
    cbind(., mesgaran)
  
  #---------------------------------------------
  # Return a list with univariate, combinatorial, and analog conditions as separate elements
  #---------------------------------------------
  
  reslist <- list(data=NULL, rasters=NULL)
  
  reslist$data$all <- mesgaran
  
  reslist$data$univariate <- mesgaran %>%
    dplyr::filter(., ExDet < 0)
  
  reslist$data$combinatorial <- mesgaran %>%
    dplyr::filter(., ExDet > 1)
  
  reslist$data$analogue <- mesgaran %>%
    dplyr::filter(., ExDet >= 0 & ExDet <= 1)
  
  #---------------------------------------------
  # Create rasters from extrapolation/MIC values
  #---------------------------------------------
  
  reslist$rasters$ExDet <- reslist$data %>%
    purrr::map(., ~ dplyr::select(., x, y, ExDet)) %>%
                 #safe_raster(.))%>%
    purrr::map(., "result")
  
  reslist$rasters$mic <- reslist$data %>%
    purrr::map(., ~ dplyr::select(., x, y, mic)) %>%
                 #safe_raster(.)) %>%
    purrr::map(., "result")
  
  #---------------------------------------------
  # Check that rasters have been produced for each extrapolation type
  #---------------------------------------------
  
  null.check <- purrr::map_lgl(.x = reslist$rasters$ExDet, .f = ~is.null(.x))
  
  ms <- names(null.check[null.check])
  ms <- purrr::map_dbl(.x = reslist$data[ms], .f = ~nrow(.x))
  # ms <- names(ms[ms>0])
  
  if(length(ms)>0){
    
    for(i in 1:length(ms)){
      
      if(ms[i]>0){
        
        # Extract data
        
        ds <- reslist$data[names(ms[i])]
        
        # Build raster
        
        predr <- prediction.grid %>%
          dplyr::select(x,y) %>%
          dplyr::mutate("ID" = 1) %>%
          raster::rasterFromXYZ(xyz = ., crs = coordinate.system)
        
        rs <- ps <- purrr::map(.x = ds,
                               .f= ~raster::rasterize(x = .x[,c("x", "y")], y = predr))
        
        # Reassign values
        
        for(i in 1:length(rs)){
          rs[[i]][!is.na(rs[[i]])] <- ds[[i]]$ExDet
          ps[[i]][!is.na(ps[[i]])] <- ds[[i]]$mic
        }
        
        reslist$rasters$ExDet <- append(reslist$rasters$ExDet, rs) %>%
          purrr::discard(is.null)
        
        reslist$rasters$mic <- append(reslist$rasters$mic, ps) %>%
          purrr::discard(is.null)
        
        
      } # End if ms[i] > 0
      
    } # End for loop length(ms)
  } # End if(length(ms)>0)
  
  #---------------------------------------------
  # Project rasters
  #---------------------------------------------
  
  for(r in 1:length(reslist$rasters$ExDet)){
    if(!is.null(reslist$rasters$ExDet[[r]]))raster::projection(reslist$rasters$ExDet[[r]]) <- coordinate.system}
  
  for(r in 1:length(reslist$rasters$mic)){
    if(!is.null(reslist$rasters$mic[[r]]))raster::projection(reslist$rasters$mic[[r]]) <- coordinate.system}
  
  #  #---------------------------------------------
  #  # Print/save summary
  #  #---------------------------------------------
  
  sumres <- summarise_extrapolation(extrapolation.object = reslist,
                                    covariate.names = covariate.names,
                                    extrapolation = TRUE,
                                    mic = TRUE)
  
  class(sumres) <- c("extrapolation_results_summary", class(sumres))
  reslist <- append(x = reslist, values = list(summary = sumres))
  
  # Add function inputs to obviate need to specify them in map()
  reslist <- append(x = reslist, values = list(
    covariate.names = covariate.names,
    samples = samples,
    prediction.grid = prediction.grid,
    coordinate.system = coordinate.system))
  
  reslist <- append(list(type = c("extrapolation", "mic")), reslist)
  
  # Keep it classy
  class(reslist) <- c("extrapolation_results", class(reslist))
  
  if(verbose) message("Done!")
  return(reslist)
  
}


map_extrapolation <- function(map.type = NULL,
                              extrapolation.object = NULL,
                              base.layer = "ocean",
                              sightings = NULL,
                              tracks = NULL,
                              verbose = TRUE){
  
  #---------------------------------------------
  # Perform function checks
  #---------------------------------------------
  
  if(is.null(map.type)) stop("Argument 'maptype' must be specified.")
  if(!map.type%in%c("extrapolation", "nearby", "mic")) stop("Unknown map type requested.")
  if(is.null(extrapolation.object) & map.type == "extrapolation") stop("Argument 'extrapolation.values' cannot be NULL when maptype is set to 'extrapolation'.")
  if(is.null(extrapolation.object) & map.type == "mic") stop("Argument 'extrapolation.values' cannot be NULL when maptype is set to 'mic'.")
  if(!map.type %in% extrapolation.object$type) stop("Map type undefined for the input extrapolation.object")
  
  #---------------------------------------------
  # Extract data
  #---------------------------------------------
  
  covariate.names <- extrapolation.object$covariate.names
  prediction.grid <- extrapolation.object$prediction.grid
  coordinate.system <- extrapolation.object$coordinate.system
  if(map.type == "nearby") gower.values <- extrapolation.object$raster
  if(map.type %in% c("extrapolation", "mic")) extrapolation.values <- extrapolation.object
  
  
  if(!is.null(sightings)){
    
    if(!class(sightings)%in%c("data.frame", "matrix", "SpatialPoints", "SpatialPointsDataFrame")) stop("Sightings must be of class matrix, data.frame, SpatialPoints or SpatialPointsDataFrame")
    
    if(is.data.frame(sightings) | is.matrix(sightings)){
      
      if(sum(c("x","y")%in%names(sightings))<2) {
        if(verbose) message("Missing x,y coordinates: sightings not shown")
        sightings = NULL
        latlon.sightings = NULL}
      
      if(!is.null(sightings)) latlon.sightings <- sum(all(range(sightings$x)>=-180 & range(sightings$x)<=180))
      
    }else{coords.sightings <- sp::coordinates(sightings)
    latlon.sightings <- sum(all(as.vector(apply(coords.sightings,2, range))>=-180 & as.vector(apply(coords.sightings,2, range))<=180))
    }
  }
  
  if(!is.null(tracks)){
    
    if(!class(tracks)%in%c("data.frame", "matrix", "SpatialLine", "SpatialLines", "SpatialLinesDataFrame")) stop("Tracks must be of class matrix, data.frame, SpatialLines or SpatialLinesDataFrame")
    
    if(is.data.frame(tracks) | is.matrix(tracks)){
      
      if(sum(c("x","y")%in%names(tracks))<2) {
        if(verbose) message("Missing x,y coordinates: tracks not shown")
        tracks = NULL
        latlon.tracks = NULL}
      
      if(!is.null(tracks)) latlon.tracks <- sum(all(range(tracks$x)>=-180 & range(tracks$x)<=180))
      
    }else{
      
      coords.tracks <- sp::coordinates(tracks)
      
      if(is.list(coords.tracks)){
        coords.tracks <- purrr::flatten(coords.tracks)
        coords.tracks <- do.call("rbind", coords.tracks)}
      
      latlon.tracks <- sum(all(as.vector(apply(coords.tracks, 2, range))>=-180 & as.vector(apply(coords.tracks, 2, range))<=180))
    }
  }
  
  baselyr <- switch(base.layer,
                    "ocean" = "Esri.OceanBasemap",
                    "world" = "Esri.WorldImagery",
                    "gray" = "Esri.WorldGrayCanvas")
  
  coordinate.system <- check_crs(coordinate.system = coordinate.system)
  
  #---------------------------------------------
  # Define coordinate systems
  #---------------------------------------------
  
  suppressWarnings(crs.ll <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  # EPSG:3857 (Web Mercator) is needed by leaflet for plotting
  
  suppressWarnings(crs.webmerc <- sp::CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"))
  
  #---------------------------------------------
  # Check which type of extrapolation occurred
  #---------------------------------------------
  
  if(map.type %in% c("mic", "extrapolation")){
    types <- purrr::map_dbl(.x = extrapolation.values$data[2:4],
                            .f = ~nrow(.))
    types <- c("Univariate", "Combinatorial", "Analogue")[types>0]}
  
  #---------------------------------------------
  # Defines survey extent (in lat/lon coords)
  #---------------------------------------------
  
  survey_extent <- methods::as(raster::extent(range(prediction.grid$x),
                                              range(prediction.grid$y)),
                               "SpatialPolygons")
  
  sp::proj4string(survey_extent) <- coordinate.system
  survey_extent <- sp::spTransform(survey_extent, CRSobj = crs.ll)
  survey_extent <- raster::extent(survey_extent)
  
  #---------------------------------------------
  # Converts segments to SpatialLines
  #---------------------------------------------
  
  if(!is.null(tracks)){
    
    if(latlon.tracks==1){crstracks <- crs.ll
    }else{crstracks <- coordinate.system}
    
    if(is.data.frame(tracks) | is.matrix(tracks)){
      
      if("TransectID"%in%names(tracks)){
        
        tracks <- tracks %>%
          split(.$ID) %>%
          purrr::map(.x = ., .f = ~dplyr::select(.x, x, y) %>%
                       as.matrix(.) %>%
                       raster::spLines(., crs = crstracks)) %>%
          do.call("rbind", .) %>%
          sp::SpatialLinesDataFrame(., data = tracks) %>%
          sp::spTransform(., CRSobj = crs.ll)
        
      }else{
        
        if(verbose) message("No transect ID detected: Linking all segments")
        
        tracks <- tracks %>%
          dplyr::select(x, y) %>%
          as.matrix(.) %>%
          raster::spLines(., crs = crstracks) %>%
          sp::SpatialLinesDataFrame(., data = tracks) %>%
          sp::spTransform(., CRSobj = crs.ll)
        
      }
      
      
    }else{tracks <- sp::spTransform(tracks, CRSobj = crs.ll)}
  }
  
  #---------------------------------------------
  # Sightings
  #---------------------------------------------
  
  if(!is.null(sightings)){
    
    if(latlon.sightings==1){crssightings <- crs.ll
    }else{crssightings <- coordinate.system}
    
    if(is.data.frame(sightings) | is.matrix(sightings)){
      
      sp::coordinates(sightings) <- ~ x + y
      sp::proj4string(sightings) <- crssightings
      
    }else{sightings <- sp::spTransform(x = sightings, CRSobj = crs.ll)}
  }
  
  #---------------------------------------------
  # Basemap
  #---------------------------------------------
  
  exleaf <- leaflet::leaflet(sightings) %>%
    leaflet::setView(lng = mean(c(survey_extent[1], survey_extent[2])),
                     lat = mean(c(survey_extent[3],survey_extent[4])),
                     zoom = 5) %>%
    leaflet::fitBounds(survey_extent[1],
                       survey_extent[3],
                       survey_extent[2],
                       survey_extent[4]) %>%
    
    leaflet::addProviderTiles(provider = baselyr)
  # leaflet::addProviderTiles(provider = providers[which(names(providers)==baselyr)][[1]])
  
  
  if(map.type == "extrapolation"){
    
    #---------------------------------------------
    # Project rasters to lat/lon for plotting
    #---------------------------------------------
    
    projr <- proj_rasters(ll = extrapolation.values$rasters$ExDet, coordinate.system = coordinate.system)
    
    #---------------------------------------------
    # Adds rasters and colour palettes
    #---------------------------------------------
    
    if("univariate"%in%names(projr)){
      
      pal.univariate <- leaflet::colorNumeric(c("#7f2704", "#fd8d3c", "#fee6ce"),
                                              raster::getValues(projr$univariate),na.color = "transparent")
      suppressWarnings(
        exleaf <- exleaf %>%
          leaflet::addRasterImage(map = .,
                                  x = projr$univariate,
                                  colors = pal.univariate,
                                  group = "Univariate",
                                  project = FALSE,
                                  opacity = 1) %>%
          
          #---------------------------------------------
        # legend
        #---------------------------------------------
        
        addLegend_decreasing(map = .,
                             pal = pal.univariate,
                             opacity = 1,
                             r.values = raster::getValues(projr$univariate),
                             decreasing = TRUE,
                             group = "Univariate",
                             title = "Univariate")
      )}
    
    if("combinatorial"%in%names(projr)){
      
      
      if(base.layer=="ocean"){
        pal.combinatorial <- leaflet::colorNumeric(c("#eadef7", "#a56bd6", "#47026d"),
                                                   raster::getValues(projr$combinatorial), na.color = "transparent")
      }else{
        pal.combinatorial <- leaflet::colorNumeric(c("#deebf7", "#6baed6", "#08519c"),
                                                   raster::getValues(projr$combinatorial), na.color = "transparent")}
      
      suppressWarnings(
        exleaf <- exleaf %>%
          leaflet::addRasterImage(map = .,
                                  x = projr$combinatorial,
                                  colors = pal.combinatorial,
                                  group = "Combinatorial",
                                  project = FALSE,
                                  opacity = 1) %>%
          
          #---------------------------------------------
        # legend
        #---------------------------------------------
        
        addLegend_decreasing(map = .,
                             pal = pal.combinatorial,
                             opacity = 1,
                             decreasing = TRUE,
                             group="Combinatorial",
                             r.values = raster::getValues(projr$combinatorial),
                             title = "Combinatorial")
      )}
    
    if("analogue"%in%names(projr)){
      
      pal.analogue <- leaflet::colorNumeric(c("#e5f5e0", "#74c476", "#00441b"),
                                            raster::getValues(projr$analogue),
                                            na.color = "transparent")
      
      suppressWarnings(
        exleaf <- exleaf %>%
          leaflet::addRasterImage(map = .,
                                  x = projr$analogue,
                                  colors = pal.analogue,
                                  group ="Analogue",
                                  project = FALSE,
                                  opacity = 1) %>%
          
          #---------------------------------------------
        # legend
        #---------------------------------------------
        
        addLegend_decreasing(map = .,
                             pal = pal.analogue,
                             opacity = 1,
                             decreasing = TRUE,
                             group="Analogue",
                             r.values = raster::getValues(projr$analogue),
                             title = "Analogue")
      )}
    
  }else if(map.type == "mic"){
    
    #---------------------------------------------
    # Adds rasters and colour palettes
    #---------------------------------------------
    
    pal <- leaflet::colorFactor(palette = c("#B2FF8C", dichromat::colorschemes$Categorical.12[c(1:4, 7:12)]),
                                domain = extrapolation.values$data$all$mic,
                                na.color = "transparent")
    
    micvars <- extrapolation.values$data$all %>%
      dplyr::count(mic)
    
    allvars <- tibble::tibble(vars = c("None", covariate.names), ID = seq(0, length(covariate.names)))
    
    mapvars <- dplyr::left_join(x = micvars, y = allvars, by = c("mic" = "ID"))
    
    suppressWarnings(
      exleaf <- exleaf %>% leaflet::addRasterImage(map = .,
                                                   x = raster::as.factor(extrapolation.values$rasters$mic$all),
                                                   colors = pal,
                                                   group = "MIC",
                                                   opacity = 1) %>%
        
        #---------------------------------------------
      # legend
      #---------------------------------------------
      
      addLegend_decreasing(map = .,
                           pal = pal,
                           labFormat = labelFormat(
                             prefix = "",
                             suffix = paste0(" - ", mapvars$vars)),
                           opacity = 1,
                           decreasing = TRUE,
                           group = "MIC",
                           r.values = raster::getValues(extrapolation.values$rasters$mic$all),
                           title = "MIC"))
    
    
  }else if(map.type == "nearby"){
    
    #---------------------------------------------
    # Adds rasters and colour palettes
    #---------------------------------------------
    
    pal.near <- leaflet::colorNumeric(pals::parula(100),
                                      raster::getValues(gower.values),
                                      na.color = "transparent")
    suppressWarnings(
      exleaf <- exleaf %>%
        leaflet::addRasterImage(map = .,
                                colors = pal.near,
                                x = gower.values,
                                group="% nearby",
                                opacity = 1) %>%
        
        #---------------------------------------------
      # legend
      #---------------------------------------------
      
      addLegend_decreasing(map = .,
                           pal = pal.near,
                           opacity = 1,
                           decreasing = TRUE,
                           group ="% nearby",
                           r.values = raster::getValues(gower.values),
                           title = "% nearby"))
    
    
  }
  
  
  #---------------------------------------------
  # Adds survey tracks and sightings
  #---------------------------------------------
  
  if(is.null(tracks)==FALSE)  exleaf <- exleaf %>%
    leaflet::addPolylines(data = tracks,
                          col = "black",
                          weight = 1,
                          group = "Tracks")
  
  if(is.null(sightings)==FALSE){
    
    if("size"%in%names(sightings)){
      
      # Create html labels for group sizes
      
      sizelist <- mapply(function(x, y) {
        htmltools::HTML(sprintf("%s: %s",
                                htmltools::htmlEscape(x),
                                htmltools::htmlEscape(y)))},
        "Group size", sightings$size, SIMPLIFY = FALSE)
      sizelist <- purrr::set_names(sizelist, NULL)
      
      exleaf <- exleaf %>%
        leaflet::addCircleMarkers(lng = sp::coordinates(sightings)[,1],
                                  lat = sp::coordinates(sightings)[,2],
                                  radius=~(log(sightings$size)+1)*3,
                                  label = sizelist,
                                  labelOptions = lapply(1:nrow(sightings),
                                                        function(x) {
                                                          labelOptions(direction='auto')}),
                                  stroke = FALSE,
                                  color="black",
                                  fillOpacity = 0.5,
                                  group = "Sightings")
    }else{
      
      if(verbose) message("No 'size' column detected: Group size not shown")
      
      exleaf <- exleaf %>%
        leaflet::addCircleMarkers(lng = sp::coordinates(sightings)[,1],
                                  lat = sp::coordinates(sightings)[,2],
                                  stroke = FALSE,
                                  radius = 4,
                                  color="black",
                                  fillOpacity = 0.5,
                                  group = "Sightings")
      
    }
  }
  
  #---------------------------------------------
  # Layer controls
  #---------------------------------------------
  
  toggles <- paste0(switch(map.type,
                           "extrapolation" = "1",
                           "mic" = "2",
                           "nearby" = "3"),
                    paste0(as.character(as.numeric(c(!is.null(tracks),
                                                     !is.null(sightings)))), collapse=""))
  
  lyr.controls <- switch(toggles,
                         "001" = c("Sightings"),
                         "010" = c("Tracks"),
                         "011" = c("Tracks", "Sightings"),
                         "100" = c(types),
                         "101" = c(types, "Sightings"),
                         "110" = c(types, "Tracks"),
                         "111" = c(types, "Tracks", "Sightings"),
                         "200" = c("MIC"),
                         "201" = c("MIC", "Sightings"),
                         "210" = c("MIC", "Tracks"),
                         "211" = c("MIC", "Tracks", "Sightings"),
                         "300" = c("% nearby"),
                         "301" = c("% nearby", "Sightings"),
                         "310" = c("% nearby", "Tracks"),
                         "311" = c("% nearby", "Tracks", "Sightings"))
  
  exleaf <- exleaf %>%
    addLayersControl(position = "topleft",
                     overlayGroups = lyr.controls,
                     options = layersControlOptions(collapsed = FALSE))
  
  if(verbose) warning('map_extrapolation relies on the leaflet package, which is built around a Web Mercator projection (EPSG:3857), and therefore requires rasters to be reprojected for plotting. As a result, minor discrepancies may  occur between the interactive maps shown in the viewer, and the underlying raw data. The latter can be accessed directly from extrapolation object returned by <compute_extrapolation> and visualised using alternative packages such as ggplot2.')
  
  return(exleaf)
}
  
summarise_extrapolation <- function(extrapolation.object,
                                    covariate.names = NULL,
                                    extrapolation = TRUE,
                                    mic = TRUE){
  #.............................................
  # Extract extrapolation values
  #---------------------------------------------
  
  ex.data <- extrapolation.object$data$all$ExDet
  
  #---------------------------------------------
  # Extract output values
  #---------------------------------------------
  
  res <- n_and_p(ex.data)
  
  #---------------------------------------------
  # Format as nice-looking table
  #---------------------------------------------
  
  # Which type(s) extrapolation did not occur?
  
  zeroes <- purrr::map_lgl(.x = res, .f = ~.x==0)
  zeroes.names <- names(zeroes[!zeroes])
  
  tb.names <- gsub(pattern = ".n", replacement = "", x = zeroes.names, fixed = TRUE) %>%
    gsub(pattern = ".p", replacement = "", x = ., fixed = TRUE) %>%
    unique(.) %>%
    tools::toTitleCase(.)
  
  #---------------------------------------------
  # Filter accordingly
  #---------------------------------------------
  
  res <- res[zeroes.names]
  
  #---------------------------------------------
  # Most influential covariates - by extrapolation type
  #---------------------------------------------
  
  if(mic){
    
    mic_data_univariate <- covariate.names[extrapolation.object$data$all$mic_univariate]
    mic_data_combinatorial <- covariate.names[extrapolation.object$data$all$mic_combinatorial]
    
    #---------------------------------------------
    # Tabulate the results
    #---------------------------------------------
    
    mic_data <- list()
    
    #---------------------------------------------
    # Univariate extrapolation
    #---------------------------------------------
    
    if(all(is.na(mic_data_univariate))){
      
      mic_data_univariate <- list()
      
    }else{
      
      mic_data_univariate <- mic_data_univariate %>%
        table(.) %>%
        as.data.frame(.) %>%
        dplyr::mutate(type = "Univariate") %>%
        dplyr::mutate(perc = 100*Freq/length(extrapolation.object$data$all$ExDet))
      
    }
    
    #---------------------------------------------
    # Combinatorial extrapolation
    #---------------------------------------------
    
    if(all(is.na(mic_data_combinatorial))){
      
      mic_data_combinatorial <- character(0)
      
    }else{
      
      mic_data_combinatorial <- mic_data_combinatorial %>%
        table(.) %>%
        as.data.frame(.) %>%
        dplyr::mutate(type = "Combinatorial") %>%
        dplyr::mutate(perc = 100*Freq/length(extrapolation.object$data$all$ExDet))
    }
    
    #---------------------------------------------
    # Combine into single list
    #---------------------------------------------
    
    if(!purrr::is_empty(mic_data_univariate)) mic_data <- append(mic_data, list(mic_data_univariate))
    if(!purrr::is_empty(mic_data_combinatorial)) mic_data <- append(mic_data, list(mic_data_combinatorial))
    
    #---------------------------------------------
    # Rename columns
    #---------------------------------------------
    
    if(!purrr::is_empty(mic_data)){
      
      mic_data <- purrr::map(.x = mic_data,
                             .f = ~purrr::set_names(.x, c("covariate", "freq", "Type", "perc")))
      
      #---------------------------------------------
      # Format list data
      #---------------------------------------------
      
      mic_res <- purrr::map(.x = mic_data,
                            .f = ~strsplit(paste(.x$freq, .x$perc), " ")) %>%
        purrr::map2(.x = ., .y = mic_data, .f = ~set_names(.x, sort(.y$covariate))) %>%
        purrr::map(.x = ., .f = ~purrr::map(.x = ., .f = ~as.numeric(.x))) %>%
        purrr::map(.x = ., .f = ~purrr::map(.x = ., .f = ~list(.n = .x[1], .p = .x[2]))) %>%
        purrr::flatten()
    }
    
  }
  #---------------------------------------------
  # Return output
  #---------------------------------------------
  
  if(mic){
    
    invisible(list(extrapolation = res,
                   mic = mic_data))
  }else{
    
    invisible(list(extrapolation = res))
  }
  
  
}
  

#---------------------------------------------
# Function to tally the number/% of cells subject to extrapolation
#---------------------------------------------

n_and_p <- function(x){
  
  exl <- list(univariate.n = length(x[x < 0]),
              univariate.p = 100 * length(x[x < 0])/length(x),
              combinatorial.n = length(x[x > 1]),
              combinatorial.p = 100 * length(x[x > 1])/length(x),
              analogue.n = length(x[x >= 0 & x <=1]),
              analogue.p = 100 * length(x[x >= 0 & x <=1])/length(x))
  return(exl)
}

check_crs <- function(coordinate.system){
  
  if(!class(coordinate.system)=="CRS"){
    
    coord.err <- tryCatch(expr = sp::CRS(coordinate.system),
                          error = function(e) return(NA))
    
    if(is.na(coord.err)){stop('Unrecognised coordinate system')
    }else{suppressWarnings(coordinate.system <- sp::CRS(coordinate.system))}
  }
  
  return(coordinate.system)
  
}     

addLegend_decreasing <- function (map, position = c("topright", "bottomright", "bottomleft", "topleft"), pal, r.values, na.label = "NA", bins = 7, colors, opacity = 0.5, labels = NULL, labFormat = labelFormat(), title = NULL, className = "info legend", layerId = NULL, group = NULL, data = getMapData(map), decreasing = FALSE) {
  
  position <- match.arg(position)
  type <- "unknown"
  na.color <- NULL
  extra <- NULL
  if (!missing(pal)) {
    if (!missing(colors))
      stop("You must provide either 'pal' or 'colors' (not both)")
    if (missing(title) && inherits(r.values, "formula"))
      title <- deparse(r.values[[2]])
    r.values <- leaflet::evalFormula(r.values, data)
    type <- attr(pal, "colorType", exact = TRUE)
    args <- attr(pal, "colorArgs", exact = TRUE)
    na.color <- args$na.color
    if (!is.null(na.color) && grDevices::col2rgb(na.color, alpha = TRUE)[[4]] ==
        0) {
      na.color <- NULL
    }
    if (type != "numeric" && !missing(bins))
      warning("'bins' is ignored because the palette type is not numeric")
    if (type == "numeric") {
      cuts <- if (length(bins) == 1)
        pretty(r.values, bins)
      else bins
      
      if (length(bins) > 2)
        if (!all(abs(diff(bins, differences = 2)) <=
                 sqrt(.Machine$double.eps)))
          stop("The vector of breaks 'bins' must be equally spaced")
      n <- length(cuts)
      r <- range(r.values, na.rm = TRUE)
      cuts <- cuts[cuts >= r[1] & cuts <= r[2]]
      n <- length(cuts)
      p <- (cuts - r[1])/(r[2] - r[1])
      extra <- list(p_1 = p[1], p_n = p[n])
      p <- c("", paste0(100 * p, "%"), "")
      if (decreasing == TRUE){
        colors <- pal(rev(c(r[1], cuts, r[2])))
        labels <- rev(labFormat(type = "numeric", cuts))
      }else{
        colors <- pal(c(r[1], cuts, r[2]))
        labels <- rev(labFormat(type = "numeric", cuts))
      }
      colors <- paste(colors, p, sep = " ", collapse = ", ")
      
    }
    else if (type == "bin") {
      cuts <- args$bins
      n <- length(cuts)
      mids <- (cuts[-1] + cuts[-n])/2
      if (decreasing == TRUE){
        colors <- pal(rev(mids))
        labels <- rev(labFormat(type = "bin", cuts))
      }else{
        colors <- pal(mids)
        labels <- labFormat(type = "bin", cuts)
      }
      
    }
    else if (type == "quantile") {
      p <- args$probs
      n <- length(p)
      cuts <- stats::quantile(r.values, probs = p, na.rm = TRUE)
      mids <- stats::quantile(r.values, probs = (p[-1] + p[-n])/2, na.rm = TRUE)
      if (decreasing == TRUE){
        colors <- pal(rev(mids))
        labels <- rev(labFormat(type = "quantile", cuts, p))
      }else{
        colors <- pal(mids)
        labels <- labFormat(type = "quantile", cuts, p)
      }
    }
    else if (type == "factor") {
      v <- sort(unique(stats::na.omit(r.values)))
      colors <- pal(v)
      labels <- labFormat(type = "factor", v)
      if (decreasing == TRUE){
        colors <- pal(rev(v))
        labels <- rev(labFormat(type = "factor", v))
      }else{
        colors <- pal(v)
        labels <- labFormat(type = "factor", v)
      }
    }
    else stop("Palette function not supported")
    if (!any(is.na(r.values)))
      na.color <- NULL
  }
  else {
    if (length(colors) != length(labels))
      stop("'colors' and 'labels' must be of the same length")
  }
  legend <- list(colors = I(unname(colors)), labels = I(unname(labels)),
                 na_color = na.color, na_label = na.label, opacity = opacity,
                 position = position, type = type, title = title, extra = extra,
                 layerId = layerId, className = className, group = group)
  invokeMethod(map, data, "addLegend", legend)
}


proj_rasters <- function(ll, coordinate.system){
  
  suppressWarnings(crs.webmerc <- sp::CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"))
  
  llr <- ll # Copy list
  
  # Univariate extrapolation is negative by definition
  # When only a small number of cells are subject to UE, the resampling
  # may result in the loss of some of them.
  # By recording the indices of UE cells, we can perform a simplistic
  # correction to make sure they show up on the map.
  
  analogue.xy <- raster::as.data.frame(llr$analogue, xy = TRUE) %>% stats::na.omit(.)
  analogue.xy <- sp::SpatialPointsDataFrame(coords = analogue.xy[, c("x", "y")],
                                            data = analogue.xy,
                                            proj4string = coordinate.system)
  analogue.xy <- sp::spTransform(analogue.xy, CRSobj = crs.webmerc)
  
  univariate.ind <- raster::Which(llr$univariate < 0, cells = TRUE)
  univariate.values <- llr$univariate[univariate.ind]
  
  univariate.xy <- raster::as.data.frame(llr$univariate, xy = TRUE) %>% stats::na.omit(.)
  univariate.xy <- sp::SpatialPointsDataFrame(coords = univariate.xy[, c("x", "y")],
                                              data = univariate.xy,
                                              proj4string = coordinate.system)
  univariate.xy <- sp::spTransform(univariate.xy, CRSobj = crs.webmerc)
  
  llr$all <- NULL
  llr <- purrr::discard(.x = llr, is.null)
  
  suppressWarnings(
    llr <- purrr::map(.x = llr, # Same extent as the full raster, allows correct alignment
                      .f = ~raster::projectRaster(from = .x,
                                                  to = ll$all,
                                                  method = 'ngb')) %>%
      purrr::map(.x = ., # CRS used by leaflet
                 .f = ~raster::projectRaster(from = .,
                                             crs = crs.webmerc,
                                             method = 'ngb')))
  
  llr.univariate.ind <- raster::cellFromXY(object = llr$univariate,
                                           xy = sp::coordinates(univariate.xy))
  
  llr$univariate[llr.univariate.ind[which(is.na(llr$univariate[llr.univariate.ind]))]] <-    univariate.values[which(is.na(llr$univariate[llr.univariate.ind]))]
  
  r1 <- raster::as.data.frame(llr$univariate, xy = TRUE)
  r2 <- raster::as.data.frame(llr$analogue, xy = TRUE)
  names(r1) <- names(r2) <- c("x", "y", "ExDet")
  
  duplicate.cells <- rbind(r1, r2) %>%
    stats::na.omit(.) %>%
    dplyr::select(., x, y) %>%
    .[duplicated(.),]
  
  llr.analogue.ind <- raster::cellFromXY(object = llr$analogue,
                                         xy = duplicate.cells)
  
  
  llr$analogue[llr.analogue.ind] <- NA
  
  
  return(llr)}

compute_nearby <- function (samples,
                            covariate.names,
                            prediction.grid,
                            coordinate.system,
                            nearby,
                            max.size = 1e7,
                            no.partitions = 10,
                            resolution = NULL,
                            verbose = TRUE) {
  
  #---------------------------------------------
  # Perform function checks
  #---------------------------------------------
  
  if(nearby<=0) stop("nearby must be strictly positive")
  if(!is.numeric(nearby)) stop("Non-numeric input to argument: nearby")
  if(max.size<=0) stop("max.size must be strictly positive")
  if(!is.numeric(max.size)) stop("Non-numeric input to argument: max.size")
  if(no.partitions>nrow(prediction.grid)) stop("Number of partitions too large")
  
  coordinate.system <- check_crs(coordinate.system = coordinate.system)
  
  #---------------------------------------------
  # Check if prediction grid is regular
  #---------------------------------------------
  
  check.grid <- prediction.grid %>%
    dplyr::select(x, y) %>%
    dplyr::mutate(z = 1)
  
  grid.regular <- try(raster::rasterFromXYZ(check.grid), silent = TRUE)
  
  #---------------------------------------------
  # If grid is irregular, rasterise prediction.grid based on specified resolution
  #---------------------------------------------
  
  if(class(grid.regular)=="try-error"){
    
    if(is.null(resolution)) stop('Prediction grid cells are not regularly spaced.\nA target raster resolution must be specified.')
    
    warning('Prediction grid cells are not regularly spaced.\nData will be rasterised and covariate values averaged.')
    
    check.grid$z <- NULL
    sp::coordinates(check.grid) <- ~x+y
    sp::proj4string(check.grid) <- coordinate.system
    
    # Create empty raster with desired resolution
    
    ras <- raster::raster(raster::extent(check.grid), res = resolution)
    raster::crs(ras) <- coordinate.system
    
    # Create individual rasters for each covariate
    
    ras.list <- purrr::map(.x = covariate.names,
                           .f = ~raster::rasterize(as.data.frame(check.grid), ras,
                                                   prediction.grid[,.x], fun = mean_ras)) %>%
      purrr::set_names(., covariate.names)
    
    # Combine all rasters
    
    ras.list <- raster::stack(ras.list)
    
    # Update prediction grid
    
    prediction.grid <- raster::as.data.frame(ras.list, xy = TRUE, na.rm = TRUE)
    
    # warning('New prediction grid (pred.grid) saved to global environment.')
    # assign(x = 'pred.grid', prediction.grid, envir = .GlobalEnv)
    
    
  } # End if class(grid.regular)
  
  #---------------------------------------------
  # Check size of input datasets
  #---------------------------------------------
  
  big.data <- ifelse(prod(nrow(samples), nrow(prediction.grid)) > max.size, TRUE, FALSE)
  
  #---------------------------------------------
  # Compute counterfactuals
  #---------------------------------------------
  
  if(big.data){
    
    counterfact <- whatif.opt(formula = NULL,
                              data = make_X(calibration_data = samples,
                                            test_data = samples,
                                            var_name = covariate.names),
                              cfact = make_X(calibration_data = samples,
                                             test_data = prediction.grid,
                                             var_name = covariate.names),
                              nearby = nearby,
                              no.partitions = no.partitions,
                              verbose = verbose)
    
  }else{
    
    counterfact <- whatif(formula = NULL,
                          data = make_X(calibration_data = samples,
                                        test_data = samples,
                                        covariate.names),
                          cfact = make_X(calibration_data = samples,
                                         test_data = prediction.grid,
                                         covariate.names),
                          nearby = nearby,
                          choice = "distance",
                          verbose = verbose)
    
  }
  
  #---------------------------------------------
  # Convert to raster
  #---------------------------------------------
  
  rgow <- cbind(prediction.grid[, c("x", "y")],
                100*counterfact$sum.stat) # in percent
  names(rgow)[3] <- 'perc_nearby'
  
  rgow <- raster::rasterFromXYZ(xyz = rgow,
                                crs = coordinate.system)
  
  
  reslist <- list(type = "nearby",
                  raster = rgow,
                  covariate.names = covariate.names,
                  samples = samples,
                  prediction.grid = prediction.grid,
                  coordinate.system = coordinate.system)
  
  
  if(verbose) message('Done!')
  return(reslist)}

whatif <- function (formula = NULL, data, cfact, range = NULL, freq = NULL,
                    nearby = 1, distance = "gower", miss = "list", choice = "both",
                    return.inputs = FALSE, return.distance = FALSE, mc.cores = 2, verbose = TRUE,
                    ...)
{
  
  if (mc.cores <= 0)
    stop("mc.cores must be an integer greater than 0.", call. = FALSE)
  if(verbose) message("Preprocessing data ...")
  
  # if (grepl("Zelig*", class(data)) & missing(cfact))
  #   cfact <- zelig_setx_to_df(data)
  # if (grepl("Zelig*", class(data)) & !missing(cfact)) {
  #   formula <- formula(delete.response(terms(data$formula)))
  #   data <- data$zelig.out$z.out[[1]]$model
  # }
  
  if (!((is.character(cfact) && is.vector(cfact) && length(cfact) ==
         1) || is.data.frame(cfact) || (is.matrix(cfact) && !is.character(cfact)))) {
    stop("'cfact' must be either a string, a R data frame, or a R non-character matrix")
  }
  if (is.character(cfact)) {
    cfact <- utils::read.table(cfact)
  }
  if (dim(cfact)[1] == 0) {
    stop("no counterfactuals supplied: 'cfact' contains zero rows")
  }
  if (!any(complete.cases(cfact))) {
    stop("there are no cases in 'cfact' without missing values")
  }
  if ("(Intercept)" %in% dimnames(cfact)[[2]]) {
    cfact <- cfact[, -(which(dimnames(cfact)[[2]] == "(Intercept)"))]
  }
  if (is.list(data) && !(is.data.frame(data))) {
    if (!((("formula" %in% names(data)) || ("terms" %in%
                                            names(data))) && (("data" %in% names(data)) || ("model" %in%
                                                                                            names(data))))) {
      stop("the list supplied to 'data' is not a valid output object")
    }
    tt <- terms(data)
    attr(tt, "intercept") <- rep(0, length(attr(tt, "intercept")))
    if ("data" %in% names(data)) {
      if (is.data.frame(data$data)) {
        data <- model.matrix(tt, model.frame(tt, data = data$data,
                                             na.action = NULL))
      }
      else {
        data <- model.matrix(tt, model.frame(tt, data = eval(data$data,
                                                             envir = .GlobalEnv), na.action = NULL))
      }
    }
    else {
      data <- model.matrix(tt, data = data$model)
    }
    if (!(is.matrix(data))) {
      stop("observed covariate data could not be extracted from output object")
    }
    rm(tt)
  }
  else {
    if (!((is.character(data) && is.vector(data) && length(data) ==
           1) || is.data.frame(data) || (is.matrix(data) &&
                                         !is.character(data)))) {
      stop("'data' must be either a string, a R data frame, a R non-character matrix, or an output object")
    }
    if (is.character(data)) {
      data <- read.table(data)
    }
  }
  if (dim(data)[1] == 0) {
    stop("no observed covariate data supplied: 'data' contains zero rows")
  }
  if (!any(complete.cases(data))) {
    stop("there are no cases in 'data' without missing values")
  }
  if (!(is.null(formula))) {
    if (identical(class(formula), "formula")) {
      if (!(is.data.frame(as.data.frame(data)))) {
        stop("'data' must be coercable to a data frame in order to use 'formula'")
      }
      if (!(is.data.frame(as.data.frame(cfact)))) {
        stop("'cfact' must be coercable to a data frame in order to use 'formula'")
      }
      formula <- update.formula(formula, ~. - 1)
      ttvar <- all.vars(formula)
      for (i in 1:length(ttvar)) {
        if (!(ttvar[i] %in% dimnames(data)[[2]])) {
          stop("variables in 'formula' either unlabeled or not present in 'data'")
        }
        if (!(ttvar[i] %in% dimnames(cfact)[[2]])) {
          stop("variable(s) in 'formula' either unlabeled or not present in 'cfact'")
        }
      }
      rm(ttvar)
      data <- model.matrix(formula, data = model.frame(formula,
                                                       as.data.frame(data), na.action = NULL))
      cfact <- model.matrix(formula, data = model.frame(formula,
                                                        as.data.frame(cfact), na.action = NULL))
    }
    else {
      stop("'formula' must be of class 'formula'")
    }
  }
  if (!(identical(complete.cases(cfact), rep(TRUE, dim(cfact)[1])))) {
    cfact <- na.omit(cfact)
    if(verbose) message("Note:  counterfactuals with missing values eliminated from cfact")
  }
  if (is.data.frame(data)) {
    if (is.character(as.matrix(data))) {
      stop("observed covariate data not coercable to numeric matrix due to character column(s)")
    }
    data <- suppressWarnings(data.matrix(data))
  }
  else {
    data <- data.matrix(as.data.frame(data))
  }
  if (is.data.frame(cfact)) {
    if (is.character(as.matrix(cfact))) {
      stop("counterfactual data not coercable to numeric matrix due to character column(s)")
    }
    cfact <- suppressWarnings(data.matrix(cfact))
  }
  else {
    cfact <- data.matrix(as.data.frame(cfact))
  }
  if (!(is.matrix(data) && is.numeric(data))) {
    stop("observed covariate data not coercable to numeric matrix")
  }
  if (!(is.matrix(cfact) && is.numeric(cfact))) {
    stop("counterfactual data not coercable to numeric matrix")
  }
  na.fail(cfact)
  if (!identical(ncol(cfact), ncol(data))) {
    stop("number of columns of 'cfact' and 'data' are not equal")
  }
  if (!(is.null(range))) {
    if (!(is.vector(range) && is.numeric(range))) {
      stop("'range' must be a numeric vector")
    }
    if (!identical(length(range), ncol(data))) {
      stop("length of 'range' does not equal number of columns of 'data'")
    }
  }
  if (!(is.null(freq))) {
    if (!(is.vector(freq) && is.numeric(freq))) {
      stop("'freq' must be a numeric vector")
    }
    na.fail(freq)
  }
  if (!(is.null(nearby))) {
    if (!(is.numeric(nearby) && is.vector(nearby) && length(nearby) ==
          1 && nearby >= 0)) {
      stop("'nearby' must be numeric, greater than or equal to 0, and a scalar")
    }
  }
  if (!(identical(miss, "list") || identical(miss, "case"))) {
    stop("'miss' must be either ''case'' or ''list''")
  }
  if (!(identical(distance, "gower") || identical(distance,
                                                  "euclidian"))) {
    stop("'distance' must be either ''gower'' or ''euclidian''")
  }
  if (!(identical(choice, "both") || identical(choice, "hull") ||
        identical(choice, "distance"))) {
    stop("'choice' must be either ''both'', ''hull'', or ''distance''")
  }
  if (!(is.logical(return.inputs))) {
    stop("'return.inputs' must be logical, i.e. either TRUE or FALSE")
  }
  if (!(is.logical(return.distance))) {
    stop("'return.distance' must be logical, i.e. either TRUE or FALSE")
  }
  n = nrow(data)
  convex.hull.test <- function(x, z, mc.cores = mc.cores) {
    one_core_pb <- mc.cores == 1
    n <- nrow(x)
    k <- ncol(x)
    m <- nrow(z)
    if (one_core_pb && m == 1)
      one_core_pb <- FALSE
    if (one_core_pb)
      pb <- txtProgressBar(min = 1, max = m, style = 3)
    A <- rbind(t(x), rep(1, n))
    C <- c(rep(0, n))
    D <- c(rep("=", k + 1))
    in_ch <- function(i, one_core_pb = FALSE) {
      B <- c(z[i, ], 1)
      lp.result <- lpSolve::lp(objective.in = C, const.mat = A,
                               const.dir = D, const.rhs = B)
      if (one_core_pb)
        setTxtProgressBar(pb, i)
      if (lp.result$status == 0)
        return(TRUE)
      else return(FALSE)
    }
    if (one_core_pb) {
      hull <- sapply(1:m, in_ch, one_core_pb = one_core_pb)
    }
    else {
      if (.Platform$OS.type == "windows")
        hull <- parallel::mclapply(1:m, in_ch, mc.cores = mc.cores)
      else hull <- pbmclapply(1:m, in_ch, mc.cores = mc.cores)
      hull <- unlist(hull)
    }
    if (one_core_pb)
      close(pb)
    return(hull)
  }
  calc.gd <- function(dat, cf, range) {
    n <- nrow(dat)
    m <- nrow(cf)
    dat = t(dat)
    dist = matrix(0, m, n, dimnames = list(1:m, 1:n))
    for (i in 1:m) {
      temp <- abs(dat - cf[i, ])/range
      if (any(range == 0)) {
        temp[is.nan(temp)] <- 0
        temp[temp == Inf] <- NA
      }
      dist[i, ] <- colMeans(temp, na.rm = T)
    }
    return(t(dist))
  }
  calc.ed <- function(dat, cf) {
    n <- nrow(dat)
    m <- nrow(cf)
    dat <- t(dat)
    dist = matrix(0, m, n, dimnames = list(1:m, 1:n))
    for (i in 1:m) {
      temp <- (dat - cf[i, ])^2
      dist[i, ] <- (colSums(temp))
    }
    return(t(dist))
  }
  geom.var <- function(dat, rang) {
    n <- nrow(dat)
    dat <- t(dat)
    ff <- function(x) {
      temp <- abs(dat - x)/rang
      if (any(rang == 0)) {
        temp[is.nan(temp)] <- 0
        temp[temp == Inf] <- NA
      }
      tmp <- sum(colMeans(temp, na.rm = TRUE))
      return(tmp)
    }
    sum.gd.x <- sum(apply(dat, 2, ff), na.rm = TRUE)
    gv.x <- (0.5 * sum.gd.x)/(n^2)
    return(gv.x)
  }
  calc.cumfreq <- function(freq, dist) {
    m <- length(freq)
    n <- ncol(dist)
    res <- matrix(0, n, m)
    for (i in 1:m) res[, i] <- (colSums(dist <= freq[i]))/nrow(dist)
    return(res)
  }
  if (identical(miss, "list")) {
    data <- na.omit(data)
    n <- nrow(data)
  }
  if ((choice == "both") | (choice == "hull")) {
    if(verbose) message("Performing convex hull test ...")
    test.result <- convex.hull.test(x = na.omit(data), z = cfact,
                                    mc.cores = mc.cores)
  }
  if ((choice == "both") | (choice == "distance")) {
    if(verbose) message("Calculating distances ....")
    if (identical(distance, "gower")) {
      samp.range <- apply(data, 2, max, na.rm = TRUE) -
        apply(data, 2, min, na.rm = TRUE)
      if (!is.null(range)) {
        w <- which(!is.na(range))
        samp.range[w] <- range[w]
      }
      if (identical(TRUE, any(samp.range == 0))) {
        if(verbose) message("Note:  range of at least one variable equals zero")
      }
      dist <- calc.gd(dat = data, cf = cfact, range = samp.range)
    }
    else {
      dist <- calc.ed(dat = na.omit(data), cf = cfact)
    }
    if(verbose) message("Calculating the geometric variance...")
    if (identical(distance, "gower")) {
      gv.x <- geom.var(dat = data, rang = samp.range)
    }
    else {
      gv.x <- 0.5 * mean(calc.ed(dat = na.omit(data), cf = na.omit(data)))
    }
    if (identical(miss, "case") && identical(distance, "euclidian")) {
      summary <- colSums(dist <= nearby * gv.x) * (1/nrow(na.omit(data)))
    }
    else {
      summary <- colSums(dist <= nearby * gv.x) * (1/n)
    }
    if(verbose) message("Calculating cumulative frequencies ...")
    if (is.null(freq)) {
      if (identical(distance, "gower")) {
        freqdist <- seq(0, 1, by = 0.05)
      }
      else {
        min.ed <- min(dist)
        max.ed <- max(dist)
        freqdist <- round(seq(min.ed, max.ed, by = (max.ed -
                                                      min.ed)/20), 2)
      }
    }
    else {
      freqdist <- freq
    }
    cumfreq <- calc.cumfreq(freq = freqdist, dist = dist)
    dimnames(cumfreq) <- list(seq(1, nrow(cfact), by = 1),
                              freqdist)
  }
  if(verbose) message("Finishing up ...")
  if (return.inputs) {
    if (choice == "both") {
      if (return.distance) {
        out <- list(call = match.call(), inputs = list(data = data,
                                                       cfact = cfact), in.hull = test.result, dist = t(dist),
                    geom.var = gv.x, sum.stat = summary, cum.freq = cumfreq)
      }
      else {
        out <- list(call = match.call(), inputs = list(data = data,
                                                       cfact = cfact), in.hull = test.result, geom.var = gv.x,
                    sum.stat = summary, cum.freq = cumfreq)
      }
    }
    if (choice == "distance") {
      if (return.distance) {
        out <- list(call = match.call(), inputs = list(data = data,
                                                       cfact = cfact), dist = t(dist), geom.var = gv.x,
                    sum.stat = summary, cum.freq = cumfreq)
      }
      else {
        out <- list(call = match.call(), inputs = list(data = data,
                                                       cfact = cfact), geom.var = gv.x, sum.stat = summary,
                    cum.freq = cumfreq)
      }
    }
    if (choice == "hull") {
      out <- list(call = match.call(), inputs = list(data = data,
                                                     cfact = cfact), in.hull = test.result)
    }
  }
  else {
    if (choice == "both") {
      if (return.distance) {
        out <- list(call = match.call(), in.hull = test.result,
                    dist = t(dist), geom.var = gv.x, sum.stat = summary,
                    cum.freq = cumfreq)
      }
      else {
        out <- list(call = match.call(), in.hull = test.result,
                    geom.var = gv.x, sum.stat = summary, cum.freq = cumfreq)
      }
    }
    if (choice == "distance") {
      if (return.distance) {
        out <- list(call = match.call(), dist = t(dist),
                    geom.var = gv.x, sum.stat = summary, cum.freq = cumfreq)
      }
      else {
        out <- list(call = match.call(), geom.var = gv.x,
                    sum.stat = summary, cum.freq = cumfreq)
      }
    }
    if (choice == "hull") {
      out <- list(call = match.call(), in.hull = test.result)
    }
  }
  class(out) <- "whatif"
  return(invisible(out))
}



whatif.opt <- function (formula = NULL,
                        data, cfact,
                        nearby = 1,
                        miss = "list",
                        no.partitions,
                        verbose = TRUE)
{
  
  if(verbose) message("Preprocessing data ...")
  
  #---------------------------------------------
  # Perform function checks
  #---------------------------------------------
  
  if (!((is.character(cfact) && is.vector(cfact) && length(cfact) ==
         1) || is.data.frame(cfact) || (is.matrix(cfact) && !is.character(cfact)))) {
    stop("'cfact' must be either a string, a R data frame, or a R non-character matrix")
  }
  
  if (is.character(cfact)) {
    cfact <- utils::read.table(cfact)
  }
  
  if (dim(cfact)[1] == 0) {
    stop("no counterfactuals supplied: 'cfact' contains zero rows")
  }
  
  if (!any(stats::complete.cases(cfact))) {
    stop("there are no cases in 'cfact' without missing values")
  }
  
  if ("(Intercept)" %in% dimnames(cfact)[[2]]) {
    cfact <- cfact[, -(which(dimnames(cfact)[[2]] == "(Intercept)"))]
  }
  
  if (is.list(data) && !(is.data.frame(data))) {
    if (!((("formula" %in% names(data)) || ("terms" %in%
                                            names(data))) && (("data" %in% names(data)) || ("model" %in%
                                                                                            names(data))))) {
      stop("the list supplied to 'data' is not a valid output object")
    }
    
    tt <- terms(data)
    attr(tt, "intercept") <- rep(0, length(attr(tt, "intercept")))
    if ("data" %in% names(data)) {
      if (is.data.frame(data$data)) {
        data <- model.matrix(tt, model.frame(tt, data = data$data,
                                             na.action = NULL))
      }else {
        data <- model.matrix(tt, model.frame(tt, data = eval(data$data,
                                                             envir = .GlobalEnv), na.action = NULL))
      }
    }else {
      data <- model.matrix(tt, data = data$model)
    }
    if (!(is.matrix(data))) {
      stop("observed covariate data could not be extracted from output object")
    }
    rm(tt)
  }else {
    if (!((is.character(data) && is.vector(data) && length(data) ==
           1) || is.data.frame(data) || (is.matrix(data) &&
                                         !is.character(data)))) {
      stop("'data' must be either a string, a R data frame, a R non-character matrix, or an output object")
    }
    if (is.character(data)) {
      data <- utils::read.table(data)
    }
  }
  
  if (dim(data)[1] == 0) {
    stop("no observed covariate data supplied: 'data' contains zero rows")
  }
  
  if (!any(stats::complete.cases(data))) {
    stop("there are no cases in 'data' without missing values")
  }
  
  if (!(is.null(formula))) {
    if (identical(class(formula), "formula")) {
      if (!(is.data.frame(as.data.frame(data)))) {
        stop("'data' must be coercable to a data frame in order to use 'formula'")
      }
      if (!(is.data.frame(as.data.frame(cfact)))) {
        stop("'cfact' must be coercable to a data frame in order to use 'formula'")
      }
      formula <- update.formula(formula, ~. - 1)
      ttvar <- all.vars(formula)
      for (i in 1:length(ttvar)) {
        if (!(ttvar[i] %in% dimnames(data)[[2]])) {
          stop("variables in 'formula' either unlabeled or not present in 'data'")
        }
        if (!(ttvar[i] %in% dimnames(cfact)[[2]])) {
          stop("variable(s) in 'formula' either unlabeled or not present in 'cfact'")
        }
      }
      rm(ttvar)
      data <- model.matrix(formula, data = model.frame(formula,
                                                       as.data.frame(data), na.action = NULL))
      cfact <- model.matrix(formula, data = model.frame(formula,
                                                        as.data.frame(cfact), na.action = NULL))
    }else {
      stop("'formula' must be of class 'formula'")
    }
  }
  
  if (!(identical(stats::complete.cases(cfact), rep(TRUE, dim(cfact)[1])))) {
    cfact <- na.omit(cfact)
    if(verbose) message("Note:  counterfactuals with missing values eliminated from cfact")
  }
  
  if (is.data.frame(data)) {
    if (is.character(as.matrix(data))) {
      stop("observed covariate data not coercable to numeric matrix due to character column(s)")
    }
    data <- suppressWarnings(data.matrix(data))
  }else {
    data <- data.matrix(as.data.frame(data))
  }
  
  if (is.data.frame(cfact)) {
    if (is.character(as.matrix(cfact))) {
      stop("counterfactual data not coercable to numeric matrix due to character column(s)")
    }
    cfact <- suppressWarnings(data.matrix(cfact))
  }else{
    cfact <- data.matrix(as.data.frame(cfact))
  }
  
  if (!(is.matrix(data) && is.numeric(data))) {
    stop("observed covariate data not coercable to numeric matrix")
  }
  
  if (!(is.matrix(cfact) && is.numeric(cfact))) {
    stop("counterfactual data not coercable to numeric matrix")
  }
  na.fail(cfact)
  
  if (!identical(ncol(cfact), ncol(data))) {
    stop("number of columns of 'cfact' and 'data' are not equal")
  }
  
  
  if (!(is.null(nearby))) {
    if (!(is.numeric(nearby) && is.vector(nearby) && length(nearby) ==
          1 && nearby >= 0)) {
      stop("'nearby' must be numeric, greater than or equal to 0, and a scalar")
    }
  }
  
  if (!(identical(miss, "list") || identical(miss, "case"))) {
    stop("'miss' must be either ''case'' or ''list''")
  }
  
  n = nrow(data)
  
  #---------------------------------------------
  # Define functions
  #---------------------------------------------
  
  # Original functions
  
  calc.gd <- function(dat, cf, range) {
    n <- nrow(dat)
    m <- nrow(cf)
    dat = t(dat)
    dist = matrix(0, m, n, dimnames = list(1:m, 1:n))
    for (i in 1:m) {
      temp <- abs(dat - cf[i, ])/range
      if (any(range == 0)) {
        temp[is.nan(temp)] <- 0
        temp[temp == Inf] <- NA
      }
      dist[i, ] <- colMeans(temp, na.rm = T)
    }
    return(t(dist))
  }
  
  geom.var <- function(dat, rang) {
    n <- nrow(dat)
    dat <- t(dat)
    ff <- function(x) {
      temp <- abs(dat - x)/rang
      if (any(rang == 0)) {
        temp[is.nan(temp)] <- 0
        temp[temp == Inf] <- NA
      }
      tmp <- sum(colMeans(temp, na.rm = TRUE))
      return(tmp)
    }
    sum.gd.x <- sum(apply(dat, 2, ff), na.rm = TRUE)
    gv.x <- (0.5 * sum.gd.x)/(n^2)
    return(gv.x)
  }
  
  
  calcgd <- function(dat, cf, range, split.factor = no.partitions) {
    
    # Split matrices into smaller chunks
    
    nlist <- split(1:nrow(dat),
                   cut(seq_along(1:nrow(dat)),
                       split.factor, labels = FALSE))
    mlist <- split(1:nrow(cf),
                   cut(seq_along(1:nrow(cf)),
                       split.factor, labels = FALSE))
    
    chunkdat <- purrr::map(.x = nlist, .f = ~dat[.x,])
    chunkcf <- purrr::map(.x = mlist, .f = ~cf[.x,])
    
    # split samples then rbind
    # split predgrid then cbind
    
    pb <- dplyr::progress_estimated(split.factor, 0) # Progress bar
    
    chunk.results <- purrr::map(.x = chunkdat,
                                .f = function(x) {
                                  pb$tick()$print()
                                  purrr::map(.x = chunkcf,
                                             function(y) calc.gd(dat = x, cf = y, range = range)) %>%
                                    do.call(cbind, .)})
    
    chunk.results <- do.call(rbind, chunk.results)
    return(chunk.results)
    
  } # End calc.gd
  
  geomvar <- function(dat, rang) {
    
    n <- nrow(dat)
    dat <- t(dat)
    
    pbb <- dplyr::progress_estimated(ncol(dat))
    # assign(x = 'pbb', value = dplyr::progress_estimated(ncol(dat)), envir = .GlobalEnv)
    
    fff <- function(x, dat, rang){
      # pbb$tick()$print()
      return(colMeans(abs(dat - dat[,x])/rang))
    }
    
    temp <- purrr::map(1:ncol(dat),
                       ~{pbb$tick()$print()
                         fff(x = .x, dat = dat, rang = rang) %>%
                           sum(.)})
    
    temp <- Reduce('+', temp)
    gv.x <- (0.5 * temp)/(n^2)
    
  }
  
  if (identical(miss, "list")) {
    data <- na.omit(data)
    n <- nrow(data)
  }
  
  #---------------------------------------------
  # Perform calculations
  #---------------------------------------------
  
  if(verbose) message("Calculating distances ....")
  
  samp.range <- apply(data, 2, max, na.rm = TRUE) - apply(data, 2, min, na.rm = TRUE)
  
  
  if (identical(TRUE, any(samp.range == 0))) {
    if(verbose) message("Note:  range of at least one variable equals zero")
  }
  
  dist <- calcgd(dat = data,
                 cf = cfact,
                 range = samp.range,
                 split.factor = no.partitions)
  
  gc()
  
  if(verbose) message("\n")
  if(verbose) message("Calculating the geometric variance ...")
  
  gv.x <- geomvar(dat = data, rang = samp.range)
  
  gc()
  
  summary <- colSums(dist <= nearby * gv.x) * (1/n)
  
  #---------------------------------------------
  # Wrap up
  #---------------------------------------------
  
  if(verbose) message("\n")
  if(verbose) message("Finishing up ...")
  
  out <- list(call = match.call(), geom.var = gv.x,
              sum.stat = summary)
  
  
  class(out) <- "whatif"
  return(invisible(out))
}


make_X <- function (calibration_data,
                    test_data,
                    var_name){
  # Changed from: rescale_cov(ynew = test_data[, k], y = calibration_data[, k]) -- June 1st, 2021
  X <- sapply(var_name, function (k) {rescale_cov(ynew = test_data[[k]], y = calibration_data[[k]])})
  X <- as.data.frame (X)
  names (X) <- var_name
  return (X)}

rescale_cov <- function (ynew, y) { return ((ynew - mean(y, na.rm = TRUE)) / (stats::sd(y, na.rm = TRUE))) }
