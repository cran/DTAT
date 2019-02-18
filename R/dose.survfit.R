#' Calculate a dose-survival curve from a dose titration study, adding a
#' confidence band
#' 
#' The 'dose-survival curve' is nothing other than an empirical cumulative
#' distribution for MTDi in the sampled population. The term 'survival' is
#' suggested in part by our application of the Kaplan-Meier estimator to
#' interval-censored toxicity information.
#' 
#' TODO: Describe details of degeneracy avoidance, once these have stabilized.
#' 
#' @param de A dose titration experiment like the \code{data} slot of class
#'   \code{\link[=DE-class]{DE}}
#' @param method The method to be used by \code{\link[km.ci]{km.ci}} when
#'  calculating CI
#' @param avoid.degeneracy When TRUE, this parameter directs the function to
#' introduce artificial events into the dose titration experiment, to avoid
#' degeneracies at the lower and upper ends of the dose-survival curve.
#' @param conf.level Confidence level for KM confidence band.
#' @return An object of class \code{survfit}.
#' @author David C. Norris
#' @seealso \code{\link{dose.survival}}, \code{\link[km.ci]{km.ci}}
#' @keywords survival
#' @examples
#' CV <- 0.7; mean_mtd <- 1.0
#' shape <- CV^-2; scale <- mean_mtd/shape
#' trial <- new("DE", doses=0.25 * 1.4^(0:6),
#'              MTDi=rgamma(24, shape=shape, scale=scale),
#'              units="mg")
#' trial <- titration(trial, periods=10)
#' sf <- dose.survfit(trial@data)
#' summary(sf)
#' 
#' @importFrom km.ci km.ci
#' @export
dose.survfit <- function(de, method="rothman", avoid.degeneracy=TRUE, conf.level=0.8){
  # To avoid degeneracy, plant an artificial DLT at lowest dose (unless already present!)
  # and also do the converse (artificial 'O') at highest dose when it shows only DLTs.
  artif.x <- 0.5
  artif.o <- 0.25 # weight applied to create 'fractional' artificial individual
  weights <- rep(1, length(unique(de$id)))
  if(avoid.degeneracy){
    if(with(de, !sum(dlt[dose==1]))){
      de <- rbind(data.frame(id=0, period=0, dose=1, dlt=TRUE), de)
      weights <- c(artif.x, weights)
    }
    if(with(de[de$dose==max(de$dose),], all(dlt))){
      de <- rbind(de, data.frame(id=Inf, period=max(de$period), dose=max(de$dose), dlt=FALSE))
      weights <- c(weights, artif.o)
    }
  }
  S <- dose.survival(de)
  fit <- survfit(S ~ 1, weights=weights)
  stopifnot(max(fit$time) == max(de$dose) || !avoid.degeneracy) # assert degeneracy avoided
  fit <- km.ci(fit, method=method, conf.level=conf.level)
}
