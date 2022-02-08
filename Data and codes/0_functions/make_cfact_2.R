# Author: Mathieu Genu & Dr Matthieu Authier
# These functions are now implemented within the pelaDSM package

#' Assess neighbourhood and extra/intepolations
#'
#' @param calibration_data a dataframe. The original data from which a design matrix is to be constructed
#' @param test_data a dataframe. Data to compare to the original data
#' @param var_name a vector of character. Names of the columns to construct design matrices.
#' @param howmany a numeric scalar: radius of the neighboorhood (as units of the geometric mean)
#' @param eps an integer: number of digits to keep for computations
#' @param rm.dup.test logical: should duplicated rows in design matrices be removed?
#' Default to FALSE.
#' @param percent logical: should the percent of neighbours/extrapolation be returned?
#' Default to TRUE.
#' @param near_by logical: should the neighbourhood be assessed instead of extrapolation?
#' Default to FALSE.
#'
#' @importFrom assertthat assert_that is.number is.count
#' @importFrom WhatIf whatif
#'
#' @return a scalar or a vector of neighbours or extrapolation flags
#' @export
#'
#' @examples
#' set.seed(123)
#' n_obs <- 1e2
#' n_test <- 1e3
#' rho <- 0.5
#' Xcal <- replicate(2, rnorm(n_obs)) %*% chol(matrix(c(1, rho, rho, 1), nrow = 2, byrow = TRUE))
#' Xcal <- as.data.frame(Xcal)
#' names(Xcal) <- paste0("X", 1:ncol(Xcal))
#' plot(Xcal,
#'      pch = 22, bg = "steelblue",
#'      type = 'p', bty = 'n', las = 1,
#'      xlab = quote(X[1]), ylab = quote(X[2])
#'      )
#'
#' Xtest <- replicate(2, rnorm(n_test))
#' Xtest <- as.data.frame(Xtest)
#' names(Xtest) <- paste0("X", 1:ncol(Xtest))
#' plot(Xtest,
#'      pch = 21, bg = "lightgreen",
#'      type = 'p', bty = 'n', las = 1,
#'      xlab = quote(X[1]), ylab = quote(X[2])
#'      )
#'
#' ### percentage of extrapolation
#' make_cfact_2(calibration_data = Xcal,
#'              test_data = Xtest
#'              )
#'
#' ### find which points are extrapolations compared to Xcal
#' Xtest$extra <- make_cfact_2(calibration_data = Xcal,
#'                             test_data = Xtest,
#'                             percent = FALSE
#'                             )
#' plot(Xtest$X1, Xtest$X2,
#'      pch = 21, col = "black",
#'      bg = ifelse(Xtest$extra == 1, "tomato", "midnightblue"),
#'      type = 'p', bty = 'n', las = 1,
#'      xlab = quote(X[1]), ylab = quote(X[2])
#'      )
#' ## add the convex hull of the calibration data
#' points(Xcal[c(chull(Xcal), chull(Xcal)[1]), ],
#'        col = "steelblue", lwd = 2,
#'        type = 'l'
#'        )
#'
#' ### mean percentage of data in Xcal to inform a prediction in Xtest
#' make_cfact_2(calibration_data = Xcal,
#'              test_data = Xtest,
#'              near_by = TRUE
#'              )
#'
#' Xtest$nearby <- make_cfact_2(calibration_data = Xcal,
#'                              test_data = Xtest,
#'                              percent = FALSE,
#'                              near_by = TRUE
#'                              )
#' plot(Xtest$X1, Xtest$X2,
#'      cex = Xtest$nearby / mean(Xtest$nearby),
#'      pch = 21, bg = "lightgreen",
#'      type = 'p', bty = 'n', las = 1,
#'      xlab = quote(X[1]), ylab = quote(X[2])
#'      )
#' points(Xcal$X1, Xcal$X2, pch = 22, bg = "steelblue", col = "white")
make_cfact_2 <- function(calibration_data,
                         test_data,
                         var_name = NULL,
                         howmany = 1,
                         eps = 6,
                         rm.dup.test = FALSE,
                         percent = TRUE,
                         near_by = FALSE
) {
  ### sanity checks
  assert_that(is.data.frame(calibration_data))
  assert_that(is.data.frame(test_data))
  assert_that(is.null(var_name) || all(sapply(var_name, is.character)))
  assert_that(is.number(howmany) && howmany > 0)
  assert_that(is.count(eps))
  assert_that(is.logical(rm.dup.test))
  assert_that(is.logical(percent))
  assert_that(is.logical(near_by))
  
  ### custom code
  if(is.null(var_name)) { var_name = names(calibration_data) }
  
  ## standardize new data to predict from
  # this simplifies computation A LOT!
  make_X <- function(calibration_data, test_data, var_name){
    X <- sapply(var_name,
                function(k) {
                  rescale2(ynew = test_data[, k],
                           y = calibration_data[, k]
                  )}
    )
    X <- as.data.frame(X)
    names(X) <- var_name
    return(X)
  }
  ### standardize
  Xcal = make_X(calibration_data = calibration_data, test_data = calibration_data, var_name)
  Xtest = make_X(calibration_data = calibration_data, test_data = test_data, var_name)
  
  # Round the standardized values
  Xcal <- round(Xcal, eps)
  Xtest <- round(Xtest, eps)
  
  # Remove duplicates
  dup <- duplicated(Xcal[, var_name])
  Xcal <- Xcal[dup == FALSE, ]
  rm(dup)
  if(rm.dup.test) {
    dup <- duplicated(Xtest[, var_name])
    Xtest <- Xtest[dup == FALSE, ]
    rm(dup)
  }
  
  # rename rows
  if(is.null(dim(Xcal))) {
    Xcal <- as.data.frame(Xcal)
    names(Xcal) <- var_name
  }
  if(is.null(dim(Xtest))) {
    Xtest <- as.data.frame(Xtest)
    names(Xtest) <- var_name
  }
  row.names(Xcal) <- 1:nrow(Xcal)
  row.names(Xtest) <- 1:nrow(Xtest)
  
  # compute counterfactuals
  if(near_by) {
    cf <- whatif(formula = NULL, data = Xcal, cfact = Xtest, choice = "distance", nearby = howmany, mc.cores = 1)
    if(percent) {
      out <- round(mean(cf$geom.var), 3)
    } else {
      out <- as.numeric(cf$sum.stat)
      out <- out / mean(out)
    }
  } else {
    cf <- whatif(formula = NULL, data = Xcal, cfact = Xtest, choice = "hull", mc.cores = 1)
    cf <- ifelse(cf$in.hull, 0, 1)
    if(percent) {
      out <- round(mean(cf), 3)
    } else {
      out <- cf
    }
  }
  return(out)
}
