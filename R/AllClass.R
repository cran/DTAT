library(methods)
# AllClass.R
#
# To achieve greater clarity through this refactoring, I must first
# abstract the object(s) that I now pass around under the unfortunate
# name 'de'.
# Let us allow -doses- and -MTDi- to be specified in absolute terms,
# with the conversion to dose-number scale being done internally.
# Note that the -de- list I build in 'step' is hardly necessary,
# except for the odd design decision to use attr(.,'stop.esc')
# for remembering when certain decisions got made. A proper class
# will use specially designated slots to hold such info!

#' An S4 class for simulating dose-titration study designs
#'
#' @slot doses A numeric vector of prospectively-determined discrete doses to
#'   trial.
#' @slot units A string indicating dose units, e.g. `"mg/kg"`.
#' @slot MTDi A numeric vector of optimal doses for simulated study
#'   participants. Optionally a call to an `r<distribution>(...)` function which
#'   may be parsed to calculate the `mtd_quantiles` slot.
#' @slot mtd_quantiles A numeric vector of quantiles of the distribution from
#'   which the MTDi slot was simulated. Intended mainly to support visualization
#'   of this distribution, e.g. as an transparent overlay on the dose-survival
#'   plot. NULL in case `MTDi` is provided verbatim.
#' @slot fractol A numeric vector of probabilities for the simulated MTDi slot.
#'   Intended mainly to support visualization, e.g. plotting of 'MTD pointers'
#'   on the interactive dose-survival plot.
#' @slot data A data.frame with columns:
#' \itemize{
#'   \item \code{id} Participant identifier
#'   \item \code{period} DLT assessment period, numbered consecutively from 1
#'   \item \code{dose} Dose level, numbered consecutively starting from 1
#'   \item \code{dlt} A logical indicator: did this this participant experience
#'     a DLT during this period?
#' }
#' @slot stop_esc integer Period in which `stop rule' was triggered
#' @slot ds_conf_level numeric Confidence level for confidence band around
#'  Kaplan-Meier estimate of the dose-survival curve.
#' @slot dose_drop_threshold numeric Threshold for triggering the `bypass rule'.
#' @slot stop_esc_under numeric Threshold for triggering the `stop rule'.
#' @slot undo_esc_under numeric Threshold for triggering the `rollback rule'.
#'
#' @export
setClass("DE"
        , slots = c(doses = "numeric"
                   ,units = "character"
                   ,MTDi = "numeric"
                   ,mtd_quantiles = "numeric"
                   ,fractol = "numeric"
                   ,data = "data.frame"  # like the former 'de[[n]]'
                   ,stop_esc = "integer" # period(s?) when escalation stops
                   ,ds_conf_level = "numeric"
                   ,dose_drop_threshold = "numeric"
                   ,stop_esc_under = "numeric"
                   ,undo_esc_under = "numeric"
                   )
        )

setMethod("initialize", "DE",
    function(.Object, doses, MTDi, ...){
      .Object <- callNextMethod(.Object, ...) # invoke the default method
      .Object@doses <- doses
      .Object@MTDi <- MTDi
      # Generate an mtd_quantiles slot if possible
      if(is.call(Q <- substitute(MTDi)) && substr(Q[[1]],1,1)=='r'){
        distname <- substring(as.character(Q[[1]]), 2) # drop the 'r'
        Q1 <- as.character(Q[[1]])
        Q[[1]] <- as.name(paste0('q',distname))
        Q[[2]] <- as.name(".p_")
        env <- parent.frame(n = 3) # caller is 3 frames back!
        assign(".p_", 1:49/50, envir = env)
        .Object@mtd_quantiles <- eval(Q, envir = env) # invoke q<dist>(p,...)
        Q[[1]] <- as.name(paste0('p',distname))
        Q[[2]] <- as.name(".q_")
        assign(".q_", MTDi, envir = env)
        .Object@fractol <- 1 - eval(Q, envir = env) # invoke p<dist>(q,...)
      } else {
        .Object@mtd_quantiles <- numeric(0)
        .Object@fractol <- numeric(0)
      }
      .Object@data <- data.frame(id=integer(0),
                                 period=integer(0),
                                 dose=integer(0),
                                 dlt=logical(0))
      .Object@data[1:3, 'id'] <- 1:3
      .Object@data[1:3, 'period'] <- 1
      .Object@data[1:3, 'dose'] <- 1
      .Object@data[1:3, 'dlt'] <- .Object@MTDi[1:3] < .Object@doses[1]
      .Object@stop_esc <- as.integer(NA)
      .Object@ds_conf_level = 0.8
      .Object@dose_drop_threshold = 0.8
      .Object@stop_esc_under = 1/3
      .Object@undo_esc_under = 1/4
      # Attach selected ... arguments as slots
      dots <- list(...)
      if(!is.null(dots$units)){
        .Object@units = dots$units
      }
      .Object
    })

setGeneric("step_time", def=function(x, ...) NULL)

setMethod("step_time", "DE",
    function(x, verbose=is.null(sys.call(-1))){
      # 0. Obtain the dose-survival curve & decide (stop|undo).esc
      dsc <- ds.curve(x@data)
      domax <- max(x@data$dose) # TODO: Move this inside the IF-block below?
      if(verbose){ cat("dose-survival curve:\n"); print(as.data.frame(dsc)); cat("domax =", domax, "\n") }
      # NB: I'll use dotted local variables 'stop.esc' and 'undo.esc',
      #     to preserve distinction against underscored @stop_esc slot.
      undo.esc <- FALSE
      if(is.na(x@stop_esc)){ # Catch the transition to the (terminal) stop.esc state
        stop.esc <- (x@stop_esc_under > dsc$upper[domax])
        if(stop.esc){ # Perhaps we must undo.esc as well; let's decide...
          undo.esc <- (x@undo_esc_under > dsc$upper[domax])
          if(verbose && undo.esc) cat("Whoa! Backing away from a too-high dose", domax, "\n")
        }
      }
      else {
        stop.esc <- x@stop_esc
      }
      if(verbose) cat("stop.esc =", stop.esc, "\n")
      # 1. Obtain the *last* period in de, 'sufficient' for intra-individual escalation decisions
      permax <- max(x@data$period)
      last <- x@data[x@data$period == permax,]
      # 2. Determine whether there will be a dose escalation
      top <- last[last$dose == max(last$dose),]
      if(!stop.esc && sum(!top$dlt) >= 3) # simple condition: must have 3+ IDs to titrate up
      {
        domax <- domax + 1
        if(domax > length(x@doses)){
          domax <- domax - 1 # "Nope."
          # Treat a pre-set max dose WLOG as though 'found' during trial:
          stop.esc <- TRUE
        }
      }
      if(undo.esc)
        domax <- domax - 1
      if(verbose) cat("domax =", domax, "\n")
      # 3. Drop individuals who have crossed their MTDi going in either direction,
      #    or who cannot tolerate the lowest dose on trial. (Note that these IDs
      #    will still be counted in the dose-survival curve; the point is simply
      #    that we do not need to carry them forward any further for purposes of
      #    visualization or escalation/dose-dropping decisions.)
      if(permax == 1){ # trivial case - drop 1st-period DLTs
        follow <- last[!last$dlt,]
        reduce.ids <- NULL
      } else { # we have 2 periods to look back on
        last2 <- x@data[x@data$period >= permax - 1 & x@data$id <= 3*(permax-1),]
        # For IDs who have crossed their MTDi's, we will have sum(dlt) = T+F (or F+T) = 1.
        # Otherwise, we'll have sum(dlt) = F+F = 0 or T+T = 2.
        crossings <- aggregate(dlt ~ id, data=last2, FUN=function(dlt) sum(dlt)==1)
        cross.ids <- crossings[crossings$dlt,]$id
        follow <- last[!(last$id %in% cross.ids) & !(last$dose==1 & last$dlt),]
        all.dlt <- aggregate(dlt ~ id, data=x@data, FUN=all)
        min.dose <- aggregate(dose ~ id, data=x@data, FUN=min)
        reduce.ids <- intersect(all.dlt[all.dlt$dlt,]$id, min.dose[min.dose$dose>1,]$id)
        # Verify that all reduce.ids are retained in 'follow':
        stopifnot(all(reduce.ids %in% follow$id))
        # Omit top-dose finalizers if we have stopped escalation
        if(stop.esc){
          maxout.ids <- x@data[x@data$period == permax & x@data$dose >= domax
                               & !x@data$dlt,]$id
          follow <- follow[!(follow$id %in% maxout.ids),]
        }
      }
      # 4. Carry forward subjects at (same | escalated | reduced) dose
      follow$period <- follow$period + 1
      follow$dose <- with(follow,
                          ifelse(id %in% reduce.ids, dose - 1, pmin(dose + 1, domax)))
      follow$dlt <- x@doses[follow$dose] > x@MTDi[follow$id]
      # 5. Abandon the lowest dose if the lower limit of 80% confidence band
      #    (i.e., 10% quantile) of computed dose-survival curve lies above 80%.
      # TODO: Move this logic up into step 1 above, right after getting 'last'
      enrolling.dose <- last$dose[which.max(last$id)] # i.e., initial dose of last-enrolled cohort
      if(x@dose_drop_threshold < dsc$lower[enrolling.dose]){
        if(verbose){
          cat(paste(dsc$lower[enrolling.dose], "<", x@dose_drop_threshold,
                    "in period", permax, "-- dropping dose", enrolling.dose, "\n"))
        }
        enrolling.dose <- enrolling.dose + 1
      }
      # 6. Add any new subjects at lowest (remaining) dose
      n <- length(unique(x@data$id))
      if(length(x@MTDi) > n){
        enroll.ids <- (n+1):min(n+3, length(x@MTDi))
        # Extract info from df on whether DLTs occur at enrolling.dose
        enroll <- data.frame(id=enroll.ids,
                             period=permax+1, #follow$period[1],
                             dose=enrolling.dose,
                             dlt=x@MTDi[enroll.ids] < x@doses[enrolling.dose])
        x@data <- rbind(x@data, enroll)
      }
      # 7. Return new value of de
      if(stop.esc && is.na(x@stop_esc)){
        if(verbose) cat("Setting 'stop.esc' attribute <-", permax, "\n")
        x@stop_esc <- as.integer(permax)
      }
      x@data <- rbind(x@data, follow)
      x
    })

#' Simulate a \sQuote{3+3/PC} dose-titration trial
#'
#' @docType methods
#' @param x An object of S4 class \code{\link[=DE-class]{DE}}
#' @param periods The number of DLT assessment periods to titrate over.
#'   Should be a positive integer.
#' @param ... May be used to pass \code{verbatim = 'TRUE'} to internal
#'   \code{step_time} method.
#'
#' @references Norris DC. Precautionary Coherence Unravels Dose Escalation Designs.
#'    \emph{bioRxiv}. December 2017:240846. doi:10.1101/240846.
#'    \url{https://www.biorxiv.org/content/10.1101/240846v1}
#'
#' @export
setGeneric("titration", def=function(x, periods, ...) NULL)

#' @rdname titration
setMethod("titration", c("DE", "numeric"),
    function(x, periods, ...){
      # To avoid pedantry, we accept 'numeric' periods:
      periods <- as.integer(round(periods))
      t <- max(x@data$period)
      while(t < periods){
        t <- t + 1
        x <- step_time(x, ...)
      }
      x
    })

#' Convert a DE object to JSON
#'
#' @param x An object of class \code{DE}
#' @param ... Unused.
#'
#' @docType methods
#' @importFrom stats lm
#' @importFrom data.table rbindlist
#' @export
setMethod("as_d3_data", "DE",
    function(x, ...){
      # Assemble a data list suitable for passing in r2d3(data=).
      #
      # Utility function for converting 'actual' doses (expressed in mg/kg, say)
      # to corresponding 'ordinal' doses on a logarithmically spaced scale.
      # Log-linear interpolation is used within the range of the @doses slot,
      # and a log-linear regression is employed to extrapolate beyond.
      # TODO: Consider converting this to a DE method,
      #       or even attaching it as a slot.
      rel_dose <- function(act_dose){
        rel <- approx(x=log(x@doses)
                     ,y=seq(length(x@doses))
                     ,xout=log(act_dose)
                     ,method="linear")$y
        dose.mult <- exp(stats::lm(log(x@doses) ~ seq(along=x@doses))$coef[2])
        extrapolated <- 1 + log(act_dose/x@doses[1]) / log(dose.mult)
        where_na <- is.na(rel)
        rel[is.na(rel)] = extrapolated[is.na(rel)]
        rel
      }
      data <- list(mtd = data.frame(id = seq_along(x@MTDi)
                                    ,mtd = x@MTDi
                                    ,doscale = rel_dose(x@MTDi)
                                    ,fractol = x@fractol
                                    )
                   ,doses = x@doses
                   ,dunit = x@units
                   ,trial = x@data
                   ,stop_esc = x@stop_esc
                   ,mtd_quantiles = rel_dose(x@mtd_quantiles)
                   ,dose_drop_threshold = x@dose_drop_threshold
                   ,stop_esc_under = x@stop_esc_under
                   ,ds = vector("list", max(x@data$period))
                   )
      # Fill out the $ds component
      for(period in 1:length(data$ds)){
        dsc <- as.data.frame(ds.curve(x@data[x@data$period <= period,]))
        dsc$dose <- seq(nrow(dsc))
        # TODO: Delegate the following data 'tweak'
        #       to the visualization code, since the
        #       need for it arises purely from
        #       visualization considerations:
        ###dsc <- dsc[c(1,1:nrow(dsc)),] # duplicate dose=1 row
        # Hmm.. What if I did the same with the top end?
        ##dsc <- dsc[c(1,1:nrow(dsc),nrow(dsc)),] # duplicate dose=1 row
        dsc <- dsc[c(1:nrow(dsc),nrow(dsc)),] # duplicate top-dose row
        dsc$dose[length(dsc$dose)] <- max(dsc$dose) + 0.5
        dsc$dose <- dsc$dose - 0.5 # shift DS curve 'steps' to *between* doses
        dsc$dose[1] <- 0.8 # avoid spilling over
        data$ds[[period]] <- dsc
      }

      data$ds <- rbindlist(data$ds, idcol="period")
      # NB: We use a *generic* transformation to JSON,
      #     and not r2d3's special default approach
      #     designed to reduce size of data transmitted:
      jsonlite::toJSON(data)
    })

#' Plot a DE object as an interactive htmlwidget
#'
#' @param x An object of class \code{DE}
#' @param y Unused; included for S4 generic consistency
#' @param ... Passed to \code{\link[r2d3]{r2d3}}, enabling caller to (e.g.) the
#'   override the default \code{viewer = "internal"}.
#' @param devtree Logical indicator used to select local package dir
#'
#' @docType methods
#' @importFrom r2d3 r2d3
#' @export
setMethod("plot", c("DE","missing"),
    function(x, y, ..., devtree=FALSE){
      script <- "htmlwidgets/lib/main.js"
      dependencies <- file.path("htmlwidgets/lib",
                                c("margins.js", "exx.js",
                                  "ox-plot.js", "ds-plot.js",
                                  "swim-plot.js", "th-plot.js"))
      if (devtree) { # refer to local package under development
        script <- file.path("DTAT/inst", script)
        dependencies <- file.path("DTAT/inst", dependencies)
        message("script = ", script)
      } else { # refer to installed package
        script <- system.file(script, package="DTAT")
        dependencies <- system.file(dependencies, package="DTAT")
      }
      r2d3::r2d3(data=as_d3_data(x),
                 script = script,
                 d3_version = 4, container = "div",
                 dependencies = dependencies,
                 ...)
    })
