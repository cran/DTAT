#' Perform neutrophil-guided dose titration of a chemotherapy drug.
#'
#' This is included in package DTAT mainly for archival purposes, with the aim
#' to document a reproduction of Figure 5 from the 2017 \emph{F1000Research}
#' paper (referenced below), using a clearer and more general software design
#' than is found in the online code supplement available at https://osf.io/vwnqz/.
#'
#' @param draw.days Integer days on which ANC is to be measured
#' @param Ncycles Number of chemo cycles through which to simulate titration
#' @param doserange Range of doses to consider
#' @param dta A Dose Titration Algorithm (DTA) to drive the titration
#' @return A list with 2 components:
#'  \item{course}{A data frame containing cycle-wise measures
#'   of each id's titration course}
#'  \item{anc.ts}{A data frame detailing high-frequency ANC measures for each id}
#' @references Norris DC. Dose Titration Algorithm Tuning (DTAT) should
#' supersede \sQuote{the} Maximum Tolerated Dose (MTD) in oncology dose-finding
#' trials. \emph{F1000Research}. 2017;6:112. doi:10.12688/f1000research.10624.3.
#' \url{https://f1000research.com/articles/6-112/v3}
#' @author David C. Norris
#' @importFrom Hmisc label<- upData
#' @importFrom stats time
#'
#' @examples
#' if(interactive()){
#' # Reproduce Figure 5 from the F1000Research paper (run time > 10 s).
#' # 1. Set up sim$pop & sim$pkpd by running the repro for Figures 1 & 3:
#' example(topic="Onoue.Friberg", package="DTAT", ask=FALSE)
#' # 2. Do the neutrophil-nadir-guided dose titration:
#' chemo <- titrate(doserange = c(50, 3000),
#'                  dta=newton.raphson(dose1 = 50,
#'                                     omega = 0.75,
#'                                     slope1 = -2.0,
#'                                     slopeU = -0.2)
#'                  )
#' library(latticeExtra)
#' newton <- chemo$course
#' new.ts <- chemo$anc.ts
#' anc.tics <- c(200,500,1500,4000,10000)
#' right <- xYplot(ANC ~ time | id, data=new.ts
#'                 , as.table=TRUE, type="l"
#'                 , layout=c(5,5)
#'                 , scales=list(y=list(log=TRUE, lim=c(100,1.5e4)
#'                                      , at=anc.tics
#'                                      , lab=as.character(anc.tics)),
#'                               x=list(at=seq(0,30,3)))
#' )
#' newton <- upData(newton
#'                  , time = 3*(cycle-1)
#'                  , labels = c(time="Time")
#'                  , units = c(time="weeks")
#'                  , print = FALSE)
#' dose.tics <- c(50, 200, 600, 1500, 3000)
#' left <- xYplot(scaled.dose ~ time | id, data=newton
#'                , as.table=TRUE, type='p', pch='+', cex=1.5
#'                , layout=c(5,5)
#'                , scales=list(y=list(lim=DTAT:::scaled(c(30,3200))
#'                                     , at=DTAT:::scaled(dose.tics)
#'                                     , lab=as.character(dose.tics)),
#'                              x=list(lim=c(-1,31)
#'                                     , at=seq(0,30,3)
#'                                     , lab=c('0','','6','','12','','18','','24','','30')))
#' )
#' update(doubleYScale(left, right, add.ylab2=TRUE)
#'        , par.settings = simpleTheme(col=brewer.pal(4,"PRGn")[c(4,1)])
#' )
#' }
#'
#' @export
titrate <-
function(draw.days=NULL, Ncycles=10,
         doserange=c(50,500), dta=NULL) {
  stopifnot(sim$N <= nrow(sim$pop))
  # Find the ANC nadirs of all 25 IDs, checking ANCs on (integer-vector) draw.days
  # We will accumulate data about each course of treatment into this data frame.
  anc.ts <- data.frame() # This will be used to collect an hourly 'Circ' time series
  course <- expand.grid(cycle=1:Ncycles, id=1:sim$N, Cc=0.0, Cp=0.0
                        , Prol=NA, Tx.1=NA, Tx.2=NA, Tx.3=NA, Circ=NA
                        , dose=NA
                        , CircMin=NA, tNadir=NA
                        , scaled.dose=NA
  )
  for (day in draw.days) {
    newcolumn <- paste("ANC", day, sep="_d")
    course[,newcolumn] <- NA
    units(course[,newcolumn]) <- "cells/mm^3"
    label(course[,newcolumn]) <- paste("Day-",day," ANC", sep="")
  }
  course$dose <- seq(scaled, from=min(doserange), to=max(doserange), length.out=max(course$cycle), digits=0)[course$cycle]
  statevector <- c('Cc','Cp','Prol','Tx.1','Tx.2','Tx.3','Circ')
  course[,statevector[-(1:2)]] <- sim$pop$Circ0[course$id] # Prol(0)=Tx.1(0)=Tx.2(0)=Tx.3(0)=Circ(0):=Circ0
  paramset <-
    function(id, states=NULL, Tinfusion=1.0, dose1=50){
      id <- as.integer(id)
      params <- unlist(sim$pop[id,c('Circ0','gamma','Emax','EC50','CL','Q','Vc','Vp','kTR')])
      params['sigma'] <- 0.15
      params['duration'] <- Tinfusion
      if (is.null(states)) {
        params[c('Cc.0','Cp.0')] <- 0.0
      } else {
        statenames <- c('Cc','Cp','Prol','Tx.1','Tx.2','Tx.3','Circ','dose')
        stopifnot(setequal(names(states), statenames))
        params[paste(statenames,'0',sep='.')] <- states[statenames]
      }
      params['dose'] <- dose1
      unlist(params)
    }
  for (id in 1:sim$N) { # outer loop over IDs permits state cycling
    params <- paramset(id)
    Circ0 <- unname(params['Circ0'])
    recycle.state <- c(Cc = 0.0
                      ,Cp = 0.0
                      ,Prol = Circ0
                      ,Tx.1 = Circ0
                      ,Tx.2 = Circ0
                      ,Tx.3 = Circ0
                      ,Circ = Circ0
    )

    for (cycle in 1:max(course$cycle)) {
      idx <- which(course$cycle==cycle & course$id==id)
      if (!is.null(dta)) { # Override preconfigured dose
        course$dose[idx] <- dta(id, cycle, course)
        if (cycle>1)
          recycle.state <- unlist(traj[nrow(traj),statevector[1:7]]) # set components of 'real state'
      }
      params['dose'] <- course$dose[idx]
      pkpd <- pomp(sim$pkpd, rinit = function(...) recycle.state)
      traj <- trajectory(pkpd, params=params, format="data.frame")
      to.add <- data.frame(id=rep(id,length(traj$time))
                           , time=traj$time + (cycle-1)*max(pkpd@times)
                           , ANC=traj$Circ)
      anc.ts <- rbind(anc.ts, to.add)
      course[idx,statevector] <- traj[which.max(traj$time),statevector]
      # When a function y(x) is sampled around a minimum at equispaced (x-dx, x, x+dx)
      # with corresponding values (y-dy1, y, y+dy2) such that dy1 < 0 < dy2 (i.e., with
      # the middle value y being lowest value), then a quadratic interpolation yields
      # estimated minimum at (x_, y_) given by:
      #   x_ = x - dx/2 * (dy1 + dy2)/(dy2 - dy1)
      #   y_ = y - 1/8 [(dy1 + dy2)/(dy2 - dy1)]^2.
      nadirIdx <- which.min(traj$Circ)
      Dy1Dy2 <- diff(traj$Circ[nadirIdx + (-1:1)])
      SdyDdy <- sum(Dy1Dy2)/diff(Dy1Dy2)
      course[idx,'CircMin'] <- traj$Circ[nadirIdx] - (1/8)*SdyDdy^2
      course[idx,'tNadir'] <- traj$time[nadirIdx] - 0.5*diff(traj$time[nadirIdx+(0:1)])*SdyDdy
      for (day in draw.days) {
        day.idx <- which(traj$time==day*24)
        course[idx,paste("ANC", day, sep="_d")] <- traj[day.idx,'Circ']
      }
    }
  }

  course$id <- ordered(paste("id",course$id,sep=""), levels=paste("id",1:sim$N,sep=""))
  course$tNadir <- course$tNadir/24
  course$scaled.dose <- scaled(course$dose)
  course <- upData(course
                   , labels=c(CircMin="ANC nadir"
                              ,tNadir="Time of ANC nadir"
                              ,dose="Drug Dose"
                              ,scaled.dose="Drug Dose (rescaled)")
                   , units=c(CircMin="cells/mm^3"
                             ,tNadir="days"
                             ,dose="mg"
                             ,scaled.dose="mg")
                   , print=FALSE
  )

  anc.ts <- upData(anc.ts
                   , id = ordered(paste("id",id,sep="")
                                  ,levels=paste("id",1:sim$N,sep=""))
                   , time = time/(24*7)
                   , labels=c(ANC="ANC")
                   , units=c(ANC="cells/mm^3"
                             ,time="weeks")
                   , print=FALSE
  )

  list(course=course, anc.ts=anc.ts)
}
