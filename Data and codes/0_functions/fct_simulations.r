### Authors: Dr Charlotte Lambert & Dr Matthieu Authier

### simuler une distribution d'animaux dans l'espace
simul_spat <- function(predgrid,
                       density_map,
                       N = 10000, # taille estimee de la pop
                       mean_group_size = 6.04,
                       overdispersion = NULL,
                       projection_espg = 2154, # lambert93
                       verbose = FALSE,
                       soap = TRUE,
                       seed = NULL,
                       is.raster = FALSE # predgrid & density_map doivent être deja projetees
                       ) {
  if(is.raster == FALSE){
    # predict on a regular grid
    pts <- predgrid %>%
      st_centroid() %>%
      st_transform(crs = st_crs(4326)) %>%
      st_coordinates() %>%
      as.data.frame()
    names(pts) <- c("lon", "lat")
    
    # response
    study_area <- density_map %>%
      filter(!is.na(density)) %>%
      st_area() %>%
      sum() %>%
      units::drop_units()
    
    average_density_per_m = N / (mean_group_size * study_area)
    
    y <- density_map %>%
      mutate(z = log(density)) %>%
      select(lon, lat, z) %>%
      st_drop_geometry() %>%
      as.data.frame()
    
    ## load soap?
    if(!soap) {
      soap_output <- prepare_soap(
        data = y,
        # contour = NEA_sf,
        polygon_pred = study,
        ratio_simplify = 0.01,
        N = 15
      )
      
      bnd <- soap_output$bnd
      knots <- soap_output$knots
      y <- soap_output$df_cropped_data
    } else {
      load(file = "data/soap_preparation.RData")
    }
    
    mod <- mgcv::gam(z ~ 1 + s(lon, lat, bs = 'so', xt = list(bnd = bnd)),
                     knots = knots,
                     data = y
    )
    
    grid <- predgrid %>%
      st_centroid() %>%
      st_coordinates() %>%
      as.data.frame() %>%
      mutate(density = exp(predict(mod, newdata = pts)),
             density = average_density_per_m * density / mean(density, na.rm = TRUE)
      )
    
    # convert to grid
    coordinates(grid) <- ~ X + Y
    proj4string(grid) <- CRS(st_crs(projection_espg)$proj4string)
    gridded(grid) <- TRUE
    
    # to im object, must use maptools to convert an sf object to a correct im object
    # https://rdrr.io/cran/maptools/man/as.ppp.html
    X <- maptools::as.im.SpatialGridDataFrame(grid) # spatstat::integral(X)
    # use spatsat.geom::as.im() for class data.frame (sur data.frame(raster))
    
    ### generate observations according to an Inhomogeneous Poisson Point Process with spatstat IPP
    if(is.null(seed)) {
      seed <- sample.int(1e6, 1)
    }
    set.seed(seed)
    y_obs <- spatstat.core::rpoispp(lambda = X, drop = TRUE)
    
    pts <- obs <- data.frame(x = y_obs$x,
                             y = y_obs$y
    )
    
    pts <- pts %>%
      st_as_sf(coords = c("x", "y"), crs = st_crs(projection_espg)) %>%
      st_transform(crs = st_crs(4326)) %>%
      st_coordinates()
    obs$lon <- pts[, 1]
    obs$lat <- pts[, 2]
    n_obs <- nrow(obs)
    ### generate a mark: group size
    if(!is.null(overdispersion)) {
      if(!is.numeric(overdispersion)) {
        if(verbose) {
          message("Overdispersion parameter is not numeric\n\tSimulating from a Poisson")
        }
        omega <- rep(1, n_obs)
      } else {
        if(overdispersion <= 1) {
          if(verbose) {
            message("Overdispersion parameter must be > 1\n\tSimulating from a Poisson")
          }
          omega <- rep(1, n_obs)
        } else {
          omega <- rgamma(n_obs, shape = overdispersion, rate = overdispersion)
        }
      }
    } else {
      omega <- rep(1, n_obs)
    }
    obs$size <- rpois(n_obs, lambda = (mean_group_size - 1) * omega) + 1
    return(obs)
  }
  if(is.raster == TRUE){
    # predict on a regular grid
    pts <- as.data.frame(predgrid, xy = T)
    
    # response
    pixel_surface <- res(density_map)[1] * res(density_map)[2] # comme le raster est projeté, on utilise sa résolution (en m)
    study_area <- (pixel_surface * ncell(density_map)) #/ 1000000

    # average_density_per_m = N / (mean_group_size * study_area)

    # density_map <- average_density_per_m * density_map / mean(density_map, na.rm = TRUE) # rescale
    
    # to im object
    # https://rdrr.io/cran/maptools/man/as.ppp.html
    density_map <- density_map / 1000000 #(pour passer en m²)
    df_map <- as.data.frame(density_map, xy = T)
    df_map <- df_map %>% drop_na()
    
    X <- spatstat.geom::as.im(df_map)
    # use spatstat.geom::as.im() for class data.frame (sur data.frame(raster))
    
    ### generate observations according to an Inhomogeneous Poisson Point Process with spatstat IPP
    if(is.null(seed)) {
      seed <- sample.int(1e6, 1)
    }
    set.seed(seed)
    y_obs <- spatstat.core::rpoispp(lambda = X, drop = TRUE)
    
    pts <- obs <- data.frame(x = y_obs$x,
                             y = y_obs$y
    )
    
    pts <- pts %>%
      st_as_sf(coords = c("x", "y"), crs = st_crs(projection_espg)) %>%
      st_transform(crs = st_crs(4326)) %>%
      st_coordinates()
    obs$lon <- pts[, 1]
    obs$lat <- pts[, 2]
    n_obs <- nrow(obs)
    
    ### generate a mark: group size
    if(!is.null(overdispersion)) {
      if(!is.numeric(overdispersion)) {
        if(verbose) {
          message("Overdispersion parameter is not numeric\n\tSimulating from a Poisson")
        }
        omega <- rep(1, n_obs)
      } else {
        if(overdispersion <= 1) {
          if(verbose) {
            message("Overdispersion parameter must be > 1\n\tSimulating from a Poisson")
          }
          omega <- rep(1, n_obs)
        } else {
          omega <- rgamma(n_obs, shape = overdispersion, rate = overdispersion)
        }
      }
    } else {
      omega <- rep(1, n_obs)
    }
    obs$size <- rpois(n_obs, lambda = (mean_group_size - 1) * omega) + 1
    return(obs)
    
  }
 }

# distance to line transects and detection process
detection_process <- function(pts, # spatial point (en sf) projeté
                              transects, # projeté comme pts
                              sigma = NULL, 
                              verbose = FALSE) {
  ## compute distances: with geosphere (use longitude and latitude)
  # d1 <- geosphere::dist2Line(p = pts[, c("lon", "lat")],
  #                           line = transects,
  #                           distfun = geosphere::distGeo
  #                           ) %>%
  #   as.data.frame() %>%
  #   mutate(distance_km = distance / 1e3,
  #          size = pts$size
  #          ) %>%
  #   as.data.frame()
  
  ## compute distance with st_distance (use projected layers) FAAAAAR MORE FASTER
  dist_mat <- st_distance(
    x = pts,
    y = transects) 
  pts$distance <- apply(dist_mat, 1, min)
  pts$distance_km <- pts$distance / 1e3
  pts$closest.seg <- apply(dist_mat, 1, function(x) min(which(x == min(x, na.rm = TRUE))))
  pts$segID <- transects$segID[pts$closest.seg]
  
  ## detection process
  if(is.null(sigma)) {
    if(verbose) {
      message("Strip-transect: perfect detection within the first 200 meters")
    }
    esw <- 0.2
    pts$detected <- ifelse(pts$distance_km <= esw, 1, 0)
  } else {
    if(verbose) {
      message("Line-transect: assuming a half-normal key for detection function")
    }
    esw <- (pnorm(+Inf, 0, sigma) - 0.5) / dnorm(0, 0, sigma)
    proba <- exp(- (pts$distance_km)^2 / (2 * sigma * sigma))
    pts$detected <- rbinom(nrow(pts), size = 1, prob = proba)
  }
  pts <- pts %>% st_drop_geometry()

  theme_set(theme_bw(base_size = 12))
  g <- ggplot() +
    geom_sf(data = transects, color = "black") +
    geom_point(data = pts[pts$detected == 0, ],
               aes(x = x, y = y), alpha = 0.3, shape = 20
               ) +
    geom_point(data = pts[pts$detected == 1, ],
               aes(x = x, y = y), shape = 21, fill = "midnightblue"
               ) +
    annotation_scale(location = "tl", width_hint = 0.5) +
    annotation_north_arrow(location = "tr",
                           which_north = "true",
                           pad_x = unit(0.2, "cm"),
                           pad_y = unit(0.1, "cm"),
                           style = north_arrow_fancy_orienteering
                           ) +
    labs(title = "Detection process",
         caption = paste("Sightings = ", sum(pts$detected), sep = " ")
         ) +
    theme(legend.position = "bottom",
          legend.key.width = unit(0.5, "cm"),
          legend.text = element_text(size = 6),
          panel.grid = element_line(colour = "transparent"), #element_blank(),
          plot.title = element_text(lineheight = 0.8, face = "bold"),
          axis.text = element_text(size = 10),
          strip.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill = "azure"),
          panel.border = element_rect(fill = NA)
          )

  return(list(data = pts,
              fig = g,
              esw = esw
              ))
}


# # cbindlist (from mvmeta)
# cbindlist <- function(list) {
#   n <- length(list)
#   res <- NULL
#   for (i in seq(n)) res <- cbind(res, list[[i]])
#   return(res)
# }

### simuler une distribution d'animaux uniforme dans l'espace
simul_spat_unif <- function(N = 10000, # taille estimee de la pop
                       mean_group_size = 6.04,
                       overdispersion = NULL,
                       projection_espg = 2154, # lambert93
                       verbose = FALSE,
                       seed = NULL) {
    ### generate observations according to an Inhomogeneous Poisson Point Process with spatstat IPP
    if(is.null(seed)) {
      seed <- sample.int(1e6, 1)
    }
    set.seed(seed)
    y_obs <- spatstat.core::runifpoint(win = spatstat.geom::as.owin(list(xrange = c(78695.73, 510417.1),
                                                                      yrange = c(6358857, 6695834))),
                                        n = N)
    
    pts <- obs <- data.frame(x = y_obs$x,
                             y = y_obs$y  )
    
    pts <- pts %>%
      st_as_sf(coords = c("x", "y"), crs = st_crs(projection_espg)) %>%
      st_transform(crs = st_crs(4326)) %>%
      st_coordinates()
    obs$lon <- pts[, 1]
    obs$lat <- pts[, 2]
    n_obs <- nrow(obs)
    
    ### generate a mark: group size
    if(!is.null(overdispersion)) {
      if(!is.numeric(overdispersion)) {
        if(verbose) {
          message("Overdispersion parameter is not numeric\n\tSimulating from a Poisson")
        }
        omega <- rep(1, n_obs)
      } else {
        if(overdispersion <= 1) {
          if(verbose) {
            message("Overdispersion parameter must be > 1\n\tSimulating from a Poisson")
          }
          omega <- rep(1, n_obs)
        } else {
          omega <- rgamma(n_obs, shape = overdispersion, rate = overdispersion)
        }
      }
    } else {
      omega <- rep(1, n_obs)
    }
    obs$size <- rpois(n_obs, lambda = (mean_group_size - 1) * omega) + 1
    return(obs)
}

###
simul <- function(predgrid,
                       truth,
                       is.raster = TRUE, # predgrid & density_map doivent être deja projetees
                       N_obs,
                       mean_group_size,
                       projection_espg = 2154, # L93
                       design_transects, # projeté aussi
                       strip = FALSE,
                       esw = 0.163,
                       soap = TRUE,
                       seed_id = NULL,
                       n_sim = 100,
                       distance.max = 700 # in meters
                       ) {
  
  simulations <- lapply(1:n_sim, function(i) {
    if(i%%10 == 0) {
      message(paste0("\t\t simulations  = ", i))
    }
    ### order transects to ensure consistency between indexes and ID
    design_transects <- design_transects %>%
      st_as_sf(design_transects, crs = st_crs(projection_espg)) %>%
      arrange(seg_id)
    
    
    ### simule point process
    simultpp <- simul_spat(predgrid,
                           density_map = truth,
                           N = N_obs,
                           mean_group_size = mean_group_size,
                           overdispersion = NULL,
                           is.raster = T
    )
    
    
    ### exclude points not around sampled transects
    # option 1 : buffer around transects, then spatial join to keep points inside
    # # buffer around design
    # design <- design_transects %>%
    #   st_as_sf %>%
    #   st_buffer(dist = 500,
    #             endCapStyle = "FLAT", joinStyle = "MITRE") %>%
    #   as_Spatial()
    # 
    # # clip point process on buffered design
    # which.over <- secr::pointsInPolygon(simultpp[,c("x", "y")], design)
    # pts <- simultpp[which(which.over == T),]
    # # plot(design)
    # # points(simultpp[, c("x", "y")], pch = 19, cex = 0.5)
    # # points(pts[, c("x", "y")], pch = 19, cex = 0.5, col = "red")

    # option 2 : test if points are within a given distance from transects
    # set points as sf points
    pts_sf <- st_as_sf(simultpp, coords = c("x", "y"), crs = st_crs(projection_espg)) %>% # proj lambert93
      dplyr::mutate(x = sf::st_coordinates(.)[,1],
                    y = sf::st_coordinates(.)[,2])
    # check which are within distance from transects
    sparse_mat <- st_is_within_distance(x = pts_sf, 
                                        y = design_transects,
                                        dist = distance.max, 
                                        sparse = F) # get a matrix with rows corresponding to x (pts), cols to lines (design)
    which.true <- apply(sparse_mat, 1, any)
    pts <- pts_sf[which(which.true),]
    # plot(design)
    # points(simultpp[, c("x", "y")], pch = 19, cex = 0.5)
    # points(pts[, c("x", "y")], pch = 19, cex = 0.5, col = "red")
    
    
    ### detection process
    if(strip) {
      sigma_esw <- NULL
    } else {
      sigma_esw <- pelaStan::scale_hn(esw = esw)
    }
    detected <- detection_process(pts = pts_sf,
                                  transects = design_transects,
                                  sigma = sigma_esw
    )
    
    
    ### distdata
    dd <- detected$data %>%
      filter(detected == 1) %>%
      mutate(ID = closest.seg)
    
    
    ### legdata
    legdata <- NULL
    legdata <- design_transects %>%
      mutate(ID = 1:n()
      ) %>%
      st_drop_geometry() %>%
      left_join(dd %>%
                  group_by(ID) %>%
                  summarize(n_detected = sum(detected),
                            n_ind = sum(size)
                  ),
                by = "ID"
      ) %>%
      mutate(n_detected = ifelse(is.na(n_detected), 0, n_detected),
             n_ind = ifelse(is.na(n_ind), 0, n_ind)
      ) %>%
      dplyr::select(survey, transect, sea_state, subjective, lon, lat, seg_id, seg_length, day, n_detected, n_ind)

    ### output
    # return(list(N_detec = legdata %>% dplyr::select(!n_ind),
    #             N_ind = legdata %>% dplyr::select(!n_detected)))
    return(legdata)
  })
  

  ### wrap-up
  Detections <- simulations %>% 
    reduce(left_join,
           by = c("survey", "transect", "sea_state", "subjective", "lon", "lat", "seg_id", "seg_length", "day")) %>%
    dplyr::select(-contains("n_ind")) %>%
    mutate(Sample.Label = str_c(survey, seg_id, sep = "_"))
  
  Individus <- simulations %>% 
                    reduce(left_join,
                           by = c("survey", "transect", "sea_state", "subjective", "lon", "lat", "seg_id", "seg_length", "day")) %>%
    dplyr::select(-contains("n_detected")) %>%
    mutate(Sample.Label = str_c(survey, seg_id, sep = "_"))
  
  out <- list(Detections = Detections,
              Individus = Individus  )
  return(out)
}


simul_unif <- function(N_obs,
                  mean_group_size,
                  projection_espg = 2154, # L93
                  design_transects, # projeté aussi
                  strip = FALSE,
                  esw = 0.163,
                  soap = TRUE,
                  seed_id = NULL,
                  n_sim = 100,
                  distance.max = 700 # in meters
) {
  
  simulations <- lapply(1:n_sim, function(i) {
    if(i%%10 == 0) {
      message(paste0("\t\t simulations  = ", i))
    }
    ### order transects to ensure consistency between indexes and ID
    design_transects <- design_transects %>%
      st_as_sf(design_transects, crs = st_crs(projection_espg)) %>%
      arrange(seg_id)
    
    
    ### simule point process
    simultpp <- simul_spat_unif(N = N_obs,
                           mean_group_size = mean_group_size,
                           overdispersion = NULL
    )
    
     
    # test if points are within a given distance from transects
    # set points as sf points
    pts_sf <- st_as_sf(simultpp, coords = c("x", "y"), crs = st_crs(projection_espg)) %>% # proj lambert93
      dplyr::mutate(x = sf::st_coordinates(.)[,1],
                    y = sf::st_coordinates(.)[,2])
    # check which are within distance from transects
    sparse_mat <- st_is_within_distance(x = pts_sf, 
                                        y = design_transects,
                                        dist = distance.max, 
                                        sparse = F) # get a matrix with rows corresponding to x (pts), cols to lines (design)
    which.true <- apply(sparse_mat, 1, any)
    pts <- pts_sf[which(which.true),]
    # plot(design)
    # points(simultpp[, c("x", "y")], pch = 19, cex = 0.5)
    # points(pts[, c("x", "y")], pch = 19, cex = 0.5, col = "red")
    
    
    ### detection process
    if(strip) {
      sigma_esw <- NULL
    } else {
      sigma_esw <- pelaStan::scale_hn(esw = esw)
    }
    detected <- detection_process(pts = pts_sf,
                                  transects = design_transects,
                                  sigma = sigma_esw
    )
    
    
    ### distdata
    dd <- detected$data %>%
      filter(detected == 1) %>%
      mutate(ID = closest.seg)
    
    
    ### legdata
    legdata <- NULL
    legdata <- design_transects %>%
      mutate(ID = 1:n()
      ) %>%
      st_drop_geometry() %>%
      left_join(dd %>%
                  group_by(ID) %>%
                  summarize(n_detected = sum(detected),
                            n_ind = sum(size)
                  ),
                by = "ID"
      ) %>%
      mutate(n_detected = ifelse(is.na(n_detected), 0, n_detected),
             n_ind = ifelse(is.na(n_ind), 0, n_ind)
      ) %>%
      dplyr::select(survey, transect, sea_state, subjective, lon, lat, seg_id, seg_length, day, n_detected, n_ind)
    
    ### output
    # return(list(N_detec = legdata %>% dplyr::select(!n_ind),
    #             N_ind = legdata %>% dplyr::select(!n_detected)))
    return(legdata)
  })
  
  
  ### wrap-up
  Detections <- simulations %>% 
    reduce(left_join,
           by = c("survey", "transect", "sea_state", "subjective", "lon", "lat", "seg_id", "seg_length", "day")) %>%
    dplyr::select(-contains("n_ind")) %>%
    mutate(Sample.Label = str_c(survey, seg_id, sep = "_"))
  
  Individus <- simulations %>% 
    reduce(left_join,
           by = c("survey", "transect", "sea_state", "subjective", "lon", "lat", "seg_id", "seg_length", "day")) %>%
    dplyr::select(-contains("n_detected")) %>%
    mutate(Sample.Label = str_c(survey, seg_id, sep = "_"))
  
  out <- list(Detections = Detections,
              Individus = Individus  )
  return(out)
}

