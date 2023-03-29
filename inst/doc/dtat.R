## ----setup, include=FALSE---------------------------------------------------------------
.First <- function() {
    library(knitr)
    library(rms)
    library(tidyr)
    library(latticeExtra)
    library(RColorBrewer)
    library(pomp)
}
.First()
# set global chunk options
options(replace.assign=TRUE,width=90)
options(datadist="dd")
options(contrasts=c("contr.treatment","contr.treatment"))
options(tinytex.engine_args = '-shell-escape')
opts_knit$set(eval.after="fig.cap")
opts_chunk$set(include=TRUE, tidy=FALSE, tidy.opts=list(width.cutoff=70), cache=FALSE, autodep=TRUE)
set.seed(2016)

## ----Population, results='asis'---------------------------------------------------------
N <- 25
dtx.mm <- 0.808 # molar mass (mg/µM) of docetaxel
pop <- data.frame(id=1:N
                  # From Friberg et al 2002 (Table 4, row 1), taking sdlog ~= CV
                  ,Circ0=rlnorm(N, meanlog=log(5050), sdlog=0.42) # units=cells/mm^3
                  ,MTT=rlnorm(N, meanlog=log(89.3), sdlog=0.16)  # mean transit time
                  ,gamma=rlnorm(N, meanlog=log(0.163), sdlog=0.039) # feedback factor
                  ,Emax=rlnorm(N, meanlog=log(83.9), sdlog=0.33)
                  ,EC50=rlnorm(N, meanlog=log(7.17*dtx.mm), sdlog=0.50)
                  # PK params from 2-compartment docetaxel model of Onoue et al (2016)
                  ,CL=rlnorm(N, meanlog=log(32.6), sdlog=0.295)
                  ,Q =rlnorm(N, meanlog=log(5.34), sdlog=0.551)
                  ,Vc=rlnorm(N, meanlog=log(5.77), sdlog=0.1) # Onoue gives no CV% for V1
                  ,Vp=rlnorm(N, meanlog=log(11.0), sdlog=0.598) # Called 'V2' in Onoue
)

pop <- upData(pop
             ,kTR = 4/MTT
             ,units = c(Circ0="cells/mm^3"
                       ,MTT="hours"
                       ,kTR="1/hour"
                       ,CL="L/h"
                       ,Q="L/h"
                       ,Vc="L"
                       ,Vp="L"
                       )
             ,print=FALSE
             )

latex(describe(pop), file="")

## ----pomp-PKPD--------------------------------------------------------------------------
pkpd.skel <- "
  double c2p = Q*( Cc - Cp ); // central-to-peripheral flux
  DCc = (dose/duration)*(t < duration ? 1.0 : 0.0)/Vc - (CL/Vc)*Cc - c2p/Vc;
  DCp = c2p/Vp;
  // Myelosuppression model (Emax model, then dynamics per eqs 3-7 from Friberg et al 2002
  double Edrug = Emax*Cc/(Cc + EC50); // classic 'Emax model'
  DProl = (1-Edrug) * Prol * kTR * pow((Circ0 / Circ), gamma)  -  kTR * Prol;
  DTx_1 = kTR * (Prol - Tx_1);
  DTx_2 = kTR * (Tx_1 - Tx_2);
  DTx_3 = kTR * (Tx_2 - Tx_3);
  DCirc = kTR * (Tx_3 - Circ);
"
pkpd.txform <- "
  T_Circ0 = log(Circ0);
  T_kTR = log(kTR);
  T_Emax = log(Emax);
  T_EC50 = log(EC50);
  T_CL = log(CL);
  T_Q = log(Q);
  T_Vc = log(Vc);
  T_Vp = log(Vp);
  T_sigma = log(sigma);
  T_dose = log(dose);
  T_duration = log(duration);
"
pkpd.txback <- "
  Circ0 = exp(T_Circ0);
  kTR = exp(T_kTR);
  Emax = exp(T_Emax);
  EC50 = exp(T_EC50);
  CL = exp(T_CL);
  Q = exp(T_Q);
  Vc = exp(T_Vc);
  Vp = exp(T_Vp);
  sigma = exp(T_sigma);
  dose = exp(T_dose);
  duration = exp(T_duration);
"

Tmax <- 21*24 # solve for full 21 days of 3-week cycle
df <- data.frame(time=c(seq(0.0, 1.95, 0.05) # q3min for 2h, 
                       ,seq(2.0, Tmax, 1.0)) # then hourly until Tmax
                ,y=NA)
pkpd <- pomp(data = df
            , times="time", t0=0
            , skeleton=vectorfield(Csnippet(pkpd.skel))
            , statenames=c("Cc", "Cp", "Prol", "Tx.1", "Tx.2", "Tx.3", "Circ")
            , paramnames=c("Circ0","kTR","gamma","Emax","EC50","CL","Q","Vc","Vp"
                          ,"sigma","dose","duration"
                          ,"Cc.0","Cp.0","Prol.0","Tx.1.0","Tx.2.0","Tx.3.0","Circ.0")
           , partrans=parameter_trans(
               toEst = Csnippet(pkpd.txform)
              ,fromEst = Csnippet(pkpd.txback)
             )
)

## ----JustForTesting, echo=FALSE, eval=FALSE---------------------------------------------
#  id <- 7
#  traj <- trajectory(pkpd,
#                     params=c(pop[pop$id==id, -which(names(pop) %in% c('id','MTT'))]
#                              , sigma=0.05
#                              , dose=50
#                              , duration=1
#                              , Cc.0 = 0.0
#                              , Cp.0 = 0.0
#                              , Prol.0 = Circ0
#                              , Tx.1.0 = Circ0
#                              , Tx.2.0 = Circ0
#                              , Tx.3.0 = Circ0
#                              , Circ.0 = Circ0)) |> as.data.frame()

## ----FirstDose, tidy=FALSE--------------------------------------------------------------
# initial conditions and output times
dose1 <- 100 # mg
Tinfusion <- 1 # 1-hour infusion
allout <- data.frame() # accumulator for nrow(pop) individual ODE solutions
parms <- function(id, dose, duration){
  parms <- unlist(pop[id,c('Circ0','kTR','gamma','Emax','EC50','CL','Q','Vc','Vp')])
  parms['sigma'] <- 0.05
  parms['dose'] <- dose
  parms['duration'] <- duration
  parms['Cc.0'] <- parms$Cp.0 <- 0.0
  parms[c('Prol.0','Tx.1.0','Tx.2.0','Tx.3.0','Circ.0')] <- parms$Circ0
  parms
}
for (id in 1:nrow(pop)) {
  Circ0 <- pop$Circ0[id]
  out <- trajectory(pkpd,
                    params=c(pop[pop$id==id, -which(names(pop) %in% c('id','MTT'))]
                             , sigma=0.05
                             , dose=dose1
                             , duration=Tinfusion
                             , Cc.0 = 0.0
                             , Cp.0 = 0.0
                             , Prol.0 = Circ0
                             , Tx.1.0 = Circ0
                             , Tx.2.0 = Circ0
                             , Tx.3.0 = Circ0
                             , Circ.0 = Circ0)) |> as.data.frame()
  out <- out[,c('time','Cc','Cp','Prol','Tx.1','Tx.2','Tx.3','Circ')]
  out$id <- paste("id",id,sep="")
  allout <- rbind(allout, out)
}

library(data.table)
allout <- as.data.table(allout)

## When a function y(x) is sampled around a minimum at equispaced (x-dx, x, x+dx)
## with corresponding values (y-dy1, y, y+dy2) such that dy1 < 0 < dy2 (i.e., with
## the middle value y being lowest value), then a quadratic interpolation yields
## estimated minimum at (x_, y_) given by:
##   x_ = x - dx/2 * (dy1 + dy2)/(dy2 - dy1)
##   y_ = y - 1/8 [(dy1 + dy2)/(dy2 - dy1)]^2.

for(.id in unique(allout$id)) {
  with(subset(allout, id == .id), {
    nadirIdx <- which.min(Circ)
    Dy1Dy2 <- diff(Circ[nadirIdx + (-1:1)])
    SdyDdy <- sum(Dy1Dy2)/diff(Dy1Dy2)
    allout[id == .id, CircMin := Circ[nadirIdx] - (1/8)*SdyDdy^2]
    allout[id == .id, tNadir := time[nadirIdx] - 0.5*diff(time[nadirIdx+(0:1)])*SdyDdy]
  })
}

allout <- upData(allout
      ,id = ordered(id, levels=paste("id",1:N,sep="")) # Note that we may
      ,units=c(Prol="cells/mm^3" # treat the non-circulating compartments
              ,Tx.1="cells/mm^3" # on a 'circulating-equivalent basis';
              ,Tx.2="cells/mm^3" # thus we attach 'cells/mm^3' units to
              ,Tx.3="cells/mm^3" # these quantities as well.
              ,Circ="cells/mm^3"
              ,Cc="ng/mL"
              ,Cp="ng/mL"
              ,time="hours")
      ,print=FALSE
       )

## ----PKplot, fig.cap="Two-compartment pharmacokinetics of the drug."--------------------
cout <- gather(allout, key="Series", value="Concentration"
               , Cc, Cp
               , factor_key = TRUE)

label(cout$Concentration) <- "Drug Concentration"

xYplot(Concentration ~ time | id, group=Series
       , data=cout, subset=time<6
       , layout=c(5,5)
       , type='l', as.table=TRUE
       , label.curves=FALSE
       , par.settings=list(superpose.line=list(lwd=2,col=brewer.pal(4,"PRGn")[c(1,4)]))
       , scales=list(y=list(log=TRUE, lim=c(10^-3,10^1)))
       , auto.key=list(points=FALSE, lines=TRUE, columns=2))

## ----MyelosuppressionPlot, fig.cap=paste("Myelosuppression in the 5-compartment model of Friberg *et al.*, with an infused dose of", dose1, "mg. **Prol:** proliferative compartment; **Tx.n:** transit compartments; **Circ:** circulating cells.")----
mout <- gather(allout, key="Series", value="ANC"
               , Prol, Tx.1, Tx.2, Tx.3, Circ
               , factor_key = TRUE)

mout <- upData(mout
               , time = time/24
               , units = c(time="days")
               , print = FALSE)

xYplot(ANC ~ time | id, group=Series
       , data=mout
       , layout=c(5,5)
       , type='l', as.table=TRUE
       , label.curves=FALSE
       , par.settings=list(superpose.line=list(lwd=2,col=brewer.pal(11,"RdYlBu")[c(1,3,4,8,10)]))
       , scales=list(y=list(log=TRUE, lim=c(100,15000), at=c(200, 500, 1000, 2000, 5000, 10000)))
       , auto.key=list(points=FALSE, lines=TRUE, columns=5))

## ----defseq, echo=FALSE-----------------------------------------------------------------
# Define a seq.function method supporting custom-scaled plot axes.
seq.function <- function(scalefun, from, to, length.out, digits=NULL){
  x <- seq(from=scalefun(from), to=scalefun(to), length.out=length.out)
  y <- numeric(length.out)
  for(i in seq(length(y)))
    y[i] <- uniroot(function(y) scaled(y)-x[i], c(from, to))$root
  if(!is.null(digits))
    y <- round(y, digits=digits)
  y
}

## ----ANCnadirs--------------------------------------------------------------------------
scaled <- function(dose, a=4.0) dose^(1/a)
dose.scale <- list(lim=scaled(c(40,550))
                   , at=scaled(c(50, 100, 250, 500))
                   , lab=c('50','100','250','500')) # TODO: Redo limits & labels when right dosing is found
doChemo <- function(draw.days=NULL, Tcyc=3*7*24, adapt.dosing=c('Newton'), omega=0.75) {
  # Find the ANC nadirs of all 20 IDs, checking ANCs on (integer-vector) draw.days
  # We will accumulate data about each course of treatment into this data frame.
  # TODO: Consider doing away entirely with columns Cc..Circ, since these inits are available in 'out'
  hourly <- which(abs(time(pkpd) - round(time(pkpd))) < .Machine$double.eps^0.5)
  anc.ts <- data.frame() # This will be used to collect an hourly 'Circ' time series
  course <- expand.grid(cycle=1:10, id=1:nrow(pop), Cc=0.0, Cp=0.0
                        , Prol=NA, Tx.1=NA, Tx.2=NA, Tx.3=NA, Circ=NA
                        , dose=NA, ANC.nadir=NA, time.nadir=NA, scaled.dose=NA
                        , ANC.nadir.est=NA, time.nadir.est=NA)
  for (day in draw.days) {
    newcolumn <- paste("ANC", day, sep="_d")
    course[,newcolumn] <- NA
    units(course[,newcolumn]) <- "cells/mm^3"
    label(course[,newcolumn]) <- paste("Day-",day," ANC", sep="")
  }
  course$dose <- seq(scaled, from=50, to=500, length.out=max(course$cycle), digits=0)[course$cycle]
  statevector <- c('Cc','Cp','Prol','Tx.1','Tx.2','Tx.3','Circ')
  course[,statevector[-(1:2)]] <- pop$Circ0[course$id] # Prol(0)=Tx.1(0)=Tx.2(0)=Tx.3(0)=Circ(0):=Circ0
  for (id in 1:nrow(pop)) { # outer loop over IDs permits state cycling
    params <- unlist(pop[id,c('Circ0','gamma','Emax','EC50','CL','Q','Vc','Vp','kTR')])
    params['sigma'] <- 0.05
    params['duration'] <- Tinfusion
    params[c('Cc.0','Cp.0')] <- 0.0
    params[c('Prol.0','Tx.1.0','Tx.2.0','Tx.3.0','Circ.0')] <- params['Circ0']
    for (cycle in 1:max(course$cycle)) {
      idx <- which(course$cycle==cycle & course$id==id)
      if (cycle>1) {
        lag_1 <- which(course$cycle==(cycle-1) & course$id==id)
        if (adapt.dosing=='Newton') { # Override preconfigured dose
          params[paste(statevector,'0',sep='.')] <- traj[nrow(traj),statevector] # recycle end-state
          if (cycle==2) {
            slope <- -2.0
          } else { # cycle >= 3 so we also have lag -2 to look back at
            lag_2 <- which(course$cycle==(cycle-2) & course$id==id)
            dY <- log(course[lag_1,'ANC.nadir'] / course[lag_2,'ANC.nadir'])
            dX <- scaled(course[lag_1,'dose']) - scaled(course[lag_2,'dose'])
            slope <- dY/dX
            if (slope >= 0) slope <- NA # to avoid instability from dY/dX>=0 due to hysteresis
          }
          delta.scaleddose <- ifelse(is.na(slope), 0, log(500 / course[lag_1,'ANC.nadir']) / slope)
          # For safety's sake, we (asymmetrically) apply relaxation factor 'omega' to any proposed dose increase:
          delta.safer <- ifelse(delta.scaleddose > 0
                                , omega*delta.scaleddose
                                , delta.scaleddose)
          new.scaleddose <- scaled(course[lag_1,'dose']) + delta.safer
          course$dose[idx] <- uniroot(function(y) scaled(y)-new.scaleddose, c(0,100000))$root
        } 
      }
      params['dose'] <- course$dose[idx]
      traj <- trajectory(pkpd, params=params) |> as.data.frame()
      to.add <- data.frame(id=rep(id,length(hourly))
                           , time=traj$time[hourly]+(cycle-1)*Tmax
                           , ANC=traj$Circ[hourly])
      anc.ts <- rbind(anc.ts, to.add)
      course[idx,c('ANC.nadir','time.nadir')] <- traj[which.min(traj$Circ),c('Circ','time')]
      course[idx,statevector] <- traj[which.max(traj$time),statevector]
      for (day in draw.days) {
        day.idx <- which(traj$time==day*24)
        course[idx,paste("ANC", day, sep="_d")] <- traj[day.idx,'Circ']
      }
      if (length(draw.days)) {
        # TODO: Here, estimate nadir level and timing based on fitted spline
        ys <- course[idx,paste("ANC", draw.days, sep="_d")]
        fit <- spline(draw.days, log(ys), n=200)
        course[idx,'ANC.nadir.est'] <- exp(min(fit$y))
        course[idx,'time.nadir.est'] <- fit$x[which.min(fit$y)]
      }
    }
  }
  
  course <- upData(course[order(course$cycle),]
                   , id = ordered(paste("id",id,sep="")
                          ,levels=paste("id",1:N,sep=""))
                   , time.nadir = time.nadir/24
                   , scaled.dose = scaled(dose)
                   , labels=c(ANC.nadir="ANC nadir"
                              ,time.nadir="Time of ANC nadir"
                              ,dose="Drug Dose"
                              ,scaled.dose="Drug Dose (rescaled)")
                   , units=c(ANC.nadir="cells/mm^3"
                             ,time.nadir="days"
                             ,dose="mg"
                             ,scaled.dose="mg")
                   , print=FALSE
  )

  anc.ts <- upData(anc.ts
                   , id = ordered(paste("id",id,sep="")
                          ,levels=paste("id",1:N,sep=""))
                   , time = time/(24*7)
                   , labels=c(ANC="ANC")
                   , units=c(ANC="cells/mm^3"
                             ,time="weeks")
                   , print=FALSE
  )
  
  list(course=course, anc.ts=anc.ts)
} # end of function

## ----Linearize--------------------------------------------------------------------------
chemo <- doChemo(draw.days=c(5,6,7,8,11), adapt.dosing=FALSE)
course <- chemo$course
anc.ts <- chemo$anc.ts

## ----ANCnadirsPlot, fig.cap="ANC nadir vs cytotoxic dose. Note the dose axis is scaled to achieve near-linearity."----
xYplot(ANC.nadir ~ scaled.dose | id, data=course, as.table=TRUE
      , layout=c(5,5)
      , scales=list(y=list(log=TRUE, lim=c(100,10000), at=c(200, 500, 1000, 2000, 5000)),
                    x=dose.scale)
)

## ----DoseSlopes, fig.cap="Slopes of $\\log(\\textit{ANC}_{nadir})$ vs $\\sqrt[4]{\\textit{dose}}$ for 25 simulated study subjects."----
slopeForID <- function(x){
  fit <- lm(log(ANC.nadir) ~ scaled.dose, data=subset(course, id==x))
  slope <- fit$coefficients['scaled.dose']
  unname(slope)
}

slope <- sapply(levels(course$id), slopeForID)
densityplot(~slope, plot.points="rug")

## ----NadirTime, fig.cap="Timing of ANC nadir vs initial cytotoxic dose. Note the dose axis is scaled to achieve near-linearity."----
xYplot(time.nadir ~ scaled.dose | id, data=course, as.table=TRUE
      , layout=c(5,5)
      , scales=list(x=dose.scale)
)

## ----NadirTimePaths, fig.cap="Timing and level of ANC nadir vs initial cytotoxic dose. The doses are equally spaced on a power-law scale that yields a roughly linear parametrization of the nadir coordinate paths."----
##, subset=(dose %in% c(500,1000,1500,2500,4000))
xYplot(ANC.nadir ~ time.nadir | id, group=dose, data=course
      , as.table=TRUE
      , layout=c(5,5)
      , auto.key=list(columns=5, title="Dose (mg)", lines=FALSE)
      , par.settings=list(superpose.symbol=list(col=gray.colors(10, start=0.7, end=0.0)))
      , scales=list(y=list(log=TRUE, lim=c(10,6000), equispaced.log=FALSE))
       )


## ----Nadir_vs_Day7, fig.cap="\\label{fig:Nadir_vs_Day7}Day-7 ANC does not serve as suitable proxy for ANC nadir across the whole modeled population. Note that for some individuals (e.g., id3, id11, id20), the day-7 ANC may dangerously underestimate the actual nadir."----
xYplot(ANC.nadir ~ ANC_d7, data=course, group=id, type='l', as.table=TRUE
       , aspect=1
       , label.curves=list(method="offset")
       , scales=list(x=list(log=TRUE, lim=c(40, 6000), equispaced.log=FALSE),
                     y=list(log=TRUE, lim=c(40, 6000), equispaced.log=FALSE))
       , par.settings=list(superpose.line=list(col="black"))
)

## ----Newtonize, fig.cap="Dose titration by Newton's method, to a goal ANC nadir of 500 cells/mm$^3$."----
chemo <- doChemo(adapt.dosing='Newton')
newton <- chemo$course
new.ts <- chemo$anc.ts
anc.tics <- c(200,500,1500,4000,10000)
right <- xYplot(ANC ~ time | id, data=new.ts
                , as.table=TRUE, type="l"
                , layout=c(5,5)
                , scales=list(y=list(log=TRUE, lim=c(100,1.5e4)
                                    , at=anc.tics
                                    , lab=as.character(anc.tics)),
                              x=list(at=seq(0,30,3)))
                )
newton <- upData(newton
                 , time = 3*(cycle-1)
                 , labels = c(time="Time")
                 , units = c(time="weeks")
                 , print = FALSE)
dose.tics <- c(50, 200, 600, 1500, 3000)
left <- xYplot(scaled.dose ~ time | id, data=newton
               , as.table=TRUE, type='p', pch='+', cex=1.5
               , layout=c(5,5)
               , scales=list(y=list(lim=scaled(c(30,3200))
                                    , at=scaled(dose.tics)
                                    , lab=as.character(dose.tics)),
                             x=list(lim=c(-1,31)
                                    , at=seq(0,30,3)
                                    , lab=c('0','','6','','12','','18','','24','','30')))
               )
update(doubleYScale(left, right, add.ylab2=TRUE)
       , par.settings = simpleTheme(col=brewer.pal(4,"PRGn")[c(4,1)])
       )

## ----SessionInfo, echo=TRUE-------------------------------------------------------------
sessionInfo()

