#' @export
cds_detection <- function(x, ...){
  UseMethod("cds_detection")
  NextMethod()
}


#' Distance sampling analysis.
#'
#' Adjust a detection function on distance observation for one taxonomic group.
#' It gives the ESW (effective-strip witdth), its standard deviation,
#' a distance object from \link[Distance]{ds} function,
#' and it gives the detection plot of the detection function adjusted.
#'
#' @param distdata data.frame previously built with  \code{\link{prepare_obs_dsm}}.
#' @param bin \code{vector} of breaks of detection plot
#' @param key Choice bewtween two key function :
#'            \itemize{
#'              \item "halfnorm" : Half-Normal.
#'              \item "hazard" : hazard-rate.
#'            }
#' @param upper maximum value of the distance range of detection estimation. Default = \code{NULL}
#' @param is_strip_transect Is it a strip-transect analysis ?
#' If \code{TRUE} it is assumed to be in the 0.2km bandwidth.
#' @return This function return a list with :
#'         \itemize{
#'           \item{esw}{effective strip width. with confidence interval at 95 percent}
#'           \item{esw_cv}{coefficient of variation of esw in percent}
#'           \item{graph}{barplot of detections depending of distance, with confindence
#'           estimated with mcmc chains}
#'           \item{distfit}{distance object from \link[Distance]{ds} function}
#'         }
#' @import ggplot2 units glue assertthat
#' @importFrom coda HPDinterval as.mcmc
#' @importFrom Distance ds
#' @importFrom mvtnorm rmvnorm
#' @importFrom usethis ui_info
#' @export
#' @examples
#' library(pelaDSM)
#'
#' # JELLY : strip-transect
#' distdata_JELLY <- SPEE_data_dsm$list_prepare_obs_by_sp$JELLY_obs_output$distdata
#'
#' cds_JELLY <- cds_detection4dsm(
#'   distdata = distdata_JELLY,
#'   bin = NULL,           # not needed for strip-transect
#'   key = NULL,
#'   upper = NULL,
#'   is_strip_transect = TRUE
#' )
#' # get all info (abundance, density, ...) on cds adjustment
#' summary(cds_JELLY$distFit)
#'
#' # MOLMOL : "standard" distance sampling analysis
#' distdata_MOLMOL <- SPEE_data_dsm$list_prepare_obs_by_sp$MOLMOL_obs_output$distdata
#' MOLMOL_truncation <- SPEE_data_dsm$list_prepare_obs_by_sp$MOLMOL_obs_output$trunc
#'
#' cds_MOLMOL <- cds_detection4dsm(
#'   distdata = distdata_MOLMOL,
#'   bin = seq(0, as.numeric(MOLMOL_truncation), 0.05),
#'   key = "halfnorm",
#'   upper = MOLMOL_truncation,
#'   is_strip_transect = FALSE
#' )
#' # detection graph with esw in red
#' cds_MOLMOL$graph
#' # get all info (abundance, density, ...) on cds adjustment
#' summary(cds_MOLMOL$distFit)
cds_detection.distdata_dsm <- function(distdata,
                                       bin = NULL,
                                       key = NULL,
                                       upper = NULL,
                                       is_strip_transect,
                                       quiet = FALSE) {

  # modify class of distdata because in ds it must be nothing than a data.frame
  class(distdata) <- "data.frame"

  if(nrow(distdata) < 2) {
    stop(paste0("distdata should have at least 2 rows"))
  }

  # drop units and make sure everything in the same format
  distdata <- distdata %>%
    units::drop_units()


  # Seabird case
  if (is_strip_transect == TRUE) {

    if(quiet == FALSE) {
      ui_info("Strip-transect method used with a 0.2km stripwidth each side")
    }

    temp <- distdata
    temp$distance <- sample(c(0.0, 0.05, 0.1, 0.15, 0.2),
                            size = nrow(temp), replace = TRUE) * ifelse(is.na(distdata$distance), NA, 1)
    fit <- Distance::ds(temp,
                        key = "unif",
                        adjustment = "poly",
                        formula = ~1,
                        truncation = 0.2,
                        cutpoints = c(0.0, 0.05, 0.1, 0.15, 0.2),
                        quiet = TRUE)

    return(list(distFit = fit))

  # Non-seabird case
  } else {

    if(quiet == FALSE) {
      ui_info("conventionnal method used")
    }

    optional_args <- c(is.null(bin), is.null(key))
    if(any(optional_args)) {
      stop(
        glue("You must provide a value for : {glue_collapse(c('bin','key')[optional_args], ', ', last = ' and ')}")
      )
    }

    upper <- as.numeric(upper)

    fit = ds(distdata,
             key = ifelse(key == "halfnorm", "hn", "hr"),
             adjustment = NULL,
             truncation = upper,
             quiet = TRUE
    )
    ### Compute the effective strip width of the transect
    hn <- function(x, sigma) { exp(-(x^2) / (2 * sigma^2)) }
    hz <- function(x, sigma, nu) { 1 - exp(-(x/sigma)^(-nu)) }

    get_summary <- function(x, alpha = 0.05) {
      c(mean(x), HPDinterval(as.mcmc(x), prob = 1 - alpha))
    }

    x <- seq(min(bin), max(bin), length.out = 1000)

    ### approximation to the posterior distribution
    my_coef <- as.numeric(fit$ddf$par)
    my_mat <- as.matrix(solve(fit$ddf$hessian))

    if(key == "halfnorm") {
      simpleMC <- exp(rnorm(1000, my_coef, sqrt(as.numeric(my_mat))))
      esw <- (pnorm(upper, 0, simpleMC) - 0.5) / dnorm(0, 0, simpleMC)
      y <- sapply(x, function(i) {
        sapply(simpleMC, hn, x = i)
      })
    }

    if(key == "hazard") {
      simpleMC <- exp(rmvnorm(1000, my_coef, my_mat))
      esw <- sapply(1:1000,function(k){
        integrate(hz, 0, upper, sigma = simpleMC[k, 2], nu = simpleMC[k, 1])$value
      })
      y <- sapply(x, function(i) {
        sapply(1:1000, function(k) {
          hz(x = i, sigma = simpleMC[k, 2], nu = simpleMC[k, 1])
        })
      })
    }

    y <- as.data.frame(cbind(x, t(apply(y, 2, get_summary))))
    names(y) <- c("distance", "moyenne", "binf", "bsup")

    # histogramme des distances
    distdata <- subset(distdata, distance <= upper) # ajout pour assurer que le graphique prenne en compte la troncation
    dd <- hist(distdata$distance, breaks = bin, plot = FALSE)
    step <- diff(bin)
    area_hist <- sum(step*dd$density)

    drect <- data.frame(x1 = numeric(length(step)), x2 = numeric(length(step)),
                        y1 = numeric(length(step)), y2 = numeric(length(step)))

    for(i in 2:length(dd$breaks)){
      drect[i-1, ] <- c(dd$breaks[i-1], dd$breaks[i], 0, dd$density[i-1]/area_hist*mean(esw))
    }

    g_plot_detection <- ggplot() +
      geom_rect(data = drect, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
                colour = "black", fill = "white") +
      geom_vline(xintercept = mean(esw), col ="red", linetype = "dotted", size = 1.5) +
      geom_line(data = y, mapping = aes(x = distance, y = moyenne), size = 1, col = "midnightblue") +
      geom_line(data = y, mapping = aes(x = distance, y = binf), linetype = "dashed", col = "midnightblue",
                size = 1, alpha = 0.6) +
      geom_line(data = y, mapping = aes(x = distance, y = bsup), linetype = "dashed", col = "midnightblue",
                size = 1, alpha = 0.6) +
      scale_x_continuous(breaks = bin, limits = range(bin)) + xlab("Perpendicular Distance") +
      scale_y_continuous(breaks = seq(0, 1, 0.2)) + ylab("Pr(Detection)") +
      theme(plot.title = element_text(lineheight = 0.8, face = "bold"),
            axis.text = element_text(size = 10))

    esw = round(get_summary(esw), 3)
    esw_cv = 100 * round(sd(esw)/mean(esw), 3)

    return(list(graph = g_plot_detection,
                esw = esw,
                esw_cv = esw_cv,
                distFit = fit))

  }

}
