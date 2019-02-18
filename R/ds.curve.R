# Examine the dose-survival curve (with confidence band)


#' Extract the dose-survival curve, with its upper and lower confidence band
#' limits
#' 
#' This utility function simply makes the results of \code{dose.survfit}
#' available in the convenient form of a list.
#' 
#' 
#' @param de A data frame describing a dose-titration study.
#' @param \dots Passed through to function \code{dose.survfit}
#' @return A list with components \code{surv}, \code{upper} and \code{lower},
#' each containing a vector that can be indexed by dose level.
#' @author David C. Norris
#' @seealso \code{\link{dose.survfit}}
#' @keywords survival
#' @examples 
#' CV <- 0.7; mean_mtd <- 1.0
#' shape <- CV^-2; scale <- mean_mtd/shape
#' trial <- new("DE", doses=0.25 * 1.4^(0:6),
#'              MTDi=rgamma(24, shape=shape, scale=scale),
#'              units="mg")
#' trial <- titration(trial, periods=10)
#' ds.curve(trial@data)
#' 
#' @importFrom stats approx
#' @export
ds.curve <- function(de, ...){
  # 2. Get survfit
  fit <- dose.survfit(de, ...)
  surv  <- approx(fit$time, fit$surv , method="constant", xout=1:max(fit$time))$y
  lower <- approx(fit$time, fit$lower, method="constant", xout=1:max(fit$time))$y
  upper <- approx(fit$time, fit$upper, method="constant", xout=1:max(fit$time))$y
  list(surv=surv, lower=lower, upper=upper)
}
