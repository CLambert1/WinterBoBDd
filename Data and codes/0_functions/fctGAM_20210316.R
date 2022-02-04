
fit_all_gam <- function(envdata, outcome, predictors, rescale = TRUE,
                        family = "tweedie", esw = NULL, 
                        weight_by_g_naught = 1, splines_by = NULL,
                        max_cor = 0.5, nb_max_pred = 4, nb_min_pred = 1, complexity = 4
) {
  ## outcome = nom de la variable réponse
  ## esw = nom de la colonne avec l'a surface échantillonnées'esw ; data must have a "Effort" column with length of effort chunks
  ## complexity = nombre de degre de liberte (sans compter l'offset)
  ## family must be one of "negative binomial", "poisson" or "tweedie"
  ## default is "negative binomial"
  
  
  ## design matrix
  X <- envdata[, predictors]
  
  ## standardize
  if(rescale){  envdata[, predictors] <- apply(X, 2, function(x) { (x - mean(x))/sd(x) } )  }
  
  ## prepare smooth terms
  writeLines("using gam() with thin-plate splines")
  intercept <- "~ 1"
  # smoothers <- paste("s(", predictors, ", k = ", complexity, ", bs = 'tp'", ")", sep = "")
  if (is.null(splines_by)) {
    smoothers <- paste("s(", predictors, ", k = ", complexity, 
                       ", bs = 'tp')", sep = "")
  } else {
    if (!is.character(splines_by)) {
      stop("Arg. 'splines_by' must be character")
    }
    if (!any(names(envdata) == splines_by)) {
      stop(paste("No column named '", splines_by, "' in dataframe 'envdata'", 
                 sep = ""))
    }
    smoothers <- paste("te(", predictors, ",", splines_by, ", k = ", complexity, 
                       ", bs = 'tp')", sep = "")
  }
  
  ## all combinations among nb_max_pred
  writeLines("\t* Create list of models, please wait")
  all_x <- pbapply::pblapply(nb_min_pred:nb_max_pred, combn, x = length(predictors))
  
  ## check whether cross-correlation needs to be evaluated
  if(nb_max_pred == 1) {
    rm_combn <- c(rep(0, length(predictors) + 1))
  } else {
    ## identify which combination is to be removed
    rm_combn <- lapply(all_x[-1], function(mat) {
      sapply(1:ncol(mat), function(i) {
        rho <- cor(X[, mat[, i]]) ; diag(rho) <- 0
        return(max(abs(as.numeric(rho))))
      })
    })
    rm_combn <- c(c(rep(0, length(predictors) + 1)), unlist(rm_combn))
  }
  
  ## Create list of models
  mlist <- function(n, y, predictors) {
    paste(y, 
          apply(X = combn(predictors, n), MARGIN = 2, paste, collapse = " + "),
          sep = paste(intercept, "+", sep = " ")
    )
  }
  all_mods <- c(paste(outcome, intercept, sep = " "),
                unlist(lapply(nb_min_pred:nb_max_pred, mlist, y = outcome, predictors = smoothers))
  )
  
  ## remove combinations of variables that are too correlated
  all_mods2 <- all_mods[which(rm_combn < max_cor)]
  
  # suppress warnings
  options(warn=-1)
  
  ## fit the models
  if(is.null(esw)) { 
    stop("Must Provide a value for esw") 
    } else {
    if(length(esw) == 1 | length(esw) == nrow(envdata)) {
      if(any(names(envdata) == outcome) == FALSE) { 
        stop("No response variable in envdata") 
      } else {
        if(length(weight_by_g_naught) == 1 | length(weight_by_g_naught) == nrow(envdata)) {
          # rescale weights
          if(length(weight_by_g_naught) == 1) {
            w <- rep(1, nrow(envdata))
          } else{
            w <- weight_by_g_naught/mean(weight_by_g_naught)
          }
          my_gam_fct <- function(x) {
            if(family == "negative binomial") {
              model <- mgcv::gam(as.formula(x), data = envdata, offset = log(2*esw*envdata$Effort), family = nb(), method = "REML", weights = w)
            }
            if(family == "poisson") {
              model <- mgcv::gam(as.formula(x), data = envdata, offset = log(2*esw*envdata$Effort), family = poisson(), method = "REML", weights = w)
            }
            if(family == "tweedie") {
              model <- mgcv::gam(as.formula(x), data = envdata, offset = log(2*esw*envdata$Effort), family = tw(), method = "REML", weights = w)
            }
            ### store some results in a data frame
            data.frame(model = x,
                       Convergence = ifelse(model$converged, 1, 0),
                       AIC = model$aic,
                       # GCV = model$gcv.ubre,
                       ResDev = model$deviance,
                       NulDev = model$null.deviance,
                       ExpDev = 100*round(1 - model$deviance/model$null.deviance, 3),
                       RMSE = qpcR::RMSE(model)
            )
          }
          writeLines("\t* Fitting all possible models, please wait")
          all_fits <- pbapply::pblapply(all_mods2, my_gam_fct)
          all_fits <- do.call("rbind", all_fits)
          all_fits <- all_fits %>% arrange(AIC)
          ## Collapse to a data frame
          return(all_fits)
        } else { 
          stop("Please check weight for g(0): provide either a single value or a vector with same length as envdata") }
      }
    } else {
      stop("Please check esw: provide either a single value or a vector with same length as envdata")
    }
  }
}

