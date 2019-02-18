#' Extract interval-censored dose tolerance data from a dose titration study
#' 
#' Constructs a \code{\link[survival]{Surv}} object from a dose-escalation
#' experiment, using interval-censoring constructs of \code{type='interval2'}.
#' 
#' @param de A data frame describing a dose-titration study
#' @return A \code{Surv} object of type='interval2'
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
#' dose.survival(trial@data)
#' 
#' @importFrom stats aggregate
#' @importFrom dplyr full_join
#' @export
dose.survival <- function(de){
  suppressMessages({
    L <- stats::aggregate(dose ~ id, data=de, FUN=max, subset=!de$dlt)
    names(L)[2] <- 'doseL'
    suppressWarnings({ # expect min() to yield Inf often below
      R <- stats::aggregate(dose ~ id, data=de, FUN=min, subset=de$dlt)
      names(R)[2] <- 'doseR'
    })
    ds <- dplyr::full_join(L, R)
  })
  ds$doseL[is.na(ds$doseL)] <- 0
  ds$doseR[is.na(ds$doseR)] <- Inf
  ds <- ds[order(ds$id),]
  # Compute the Surv object from ds data frame with cols (id, doseL, doseR)
  S <- with(ds, Surv(time=doseL, time2=doseR, type='interval2'))
  S
}
