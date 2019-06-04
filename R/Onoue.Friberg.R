#' POMP PK/PD model for docetaxel, combining Onoue et al (2016) with Friberg et
#' al (2002)
#' 
#' This function produces a POMP model combining docetaxel pharmacokinetics
#' (PK) drawn from Table 2 of Onoue et al (2016) with myelosuppression dynamics
#' drawn from Friberg et al (2002). This model enables simulation of
#' neutrophil-guided dose titration of docetaxel, as done in Norris (2017).
#' 
#' 
#' @param N Size of simulated population.
#' @param cycle.length.days Duration (in days) of chemotherapy cycle to be
#' simulated.
#' @param data Passed through as the \code{data} argument of the \code{pomp}
#' constructor.
#' @param delta.t Time-step (in hours) of pomp's \code{euler} plug-in.
#' @return No value is returned; rather, the function sets global variables
#' in package environment \code{DTAT::sim}.
#' 
#' @references
#' 1. Onoue H, Yano I, Tanaka A, Itohara K, Hanai A, Ishiguro H, et al.
#'    Significant effect of age on docetaxel pharmacokinetics in Japanese
#'    female breast cancer patients by using the population modeling approach.
#'    \emph{Eur J Clin Pharmacol}. 2016 Jun;72(6):703-10.
#'    doi:10.1007/s00228-016-2031-3.
#' 
#' 2. Friberg LE, Henningsson A, Maas H, Nguyen L, Karlsson MO. Model of
#'    chemotherapy-induced myelosuppression with parameter consistency across
#'    drugs. \emph{J Clin Oncol}. 2002 Dec 15;20(24):4713-21.
#'    doi:10.1200/JCO.2002.02.140.
#' 
#' 3. Norris DC. Dose Titration Algorithm Tuning (DTAT) should supersede
#'    \sQuote{the} Maximum Tolerated Dose (MTD) in oncology dose-finding trials.
#'    \emph{F1000Research}. 2017;6:112. doi:10.12688/f1000research.10624.3.
#'    \url{https://f1000research.com/articles/6-112/v3}
#' 
#' @author David C. Norris
#' @seealso \code{\link{pomp}}, \code{\link{sim}}
#' @examples
#' # Reproduce the sim$pkpd model and sim$pop population from reference #3:
#' library(pomp)
#' Onoue.Friberg(N=25)
#' sim$pop # NB: this differs from pop of original paper...
#' 
#' # Whereas the present version of Onoue.Friberg() draws simulated populations
#' # using pomp::rprior(), to reproduce the original F1000Research paper [3] we
#' # re-draw sim$pop as originally & prosaically done (see https://osf.io/vwnqz/):
#' set.seed(2016)
#' N <- 25
#' dtx.mm <- 0.808 # molar mass (mg/µM) of docetaxel
#' sim$pop$Circ0 <- rlnorm(N, meanlog=log(5050), sdlog=0.42) # units=cells/mm^3
#' sim$pop$MTT <- rlnorm(N, meanlog=log(89.3), sdlog=0.16)  # mean transit time
#' sim$pop$gamma <- rlnorm(N, meanlog=log(0.163), sdlog=0.039) # feedback factor
#' sim$pop$Emax <- rlnorm(N, meanlog=log(83.9), sdlog=0.33)
#' sim$pop$EC50 <- rlnorm(N, meanlog=log(7.17*dtx.mm), sdlog=0.50)
#' # PK params from 2-compartment docetaxel model of Onoue et al (2016)
#' sim$pop$CL <- rlnorm(N, meanlog=log(32.6), sdlog=0.295)
#' sim$pop$Q  <- rlnorm(N, meanlog=log(5.34), sdlog=0.551)
#' sim$pop$Vc <- rlnorm(N, meanlog=log(5.77), sdlog=0.1) # Onoue gives no CV% for V1
#' sim$pop$Vp <- rlnorm(N, meanlog=log(11.0), sdlog=0.598) # Called 'V2' in Onoue
#' sim$pop$kTR=4/sim$pop$MTT
#' 
#' # Now we run the sim$pkpd model, separately for each of N simultated individuals:
#' allout <- data.frame() # accumulator for N individual ODE solutions
#' for (id in 1:sim$N) {
#'   out <- trajectory(sim$pkpd,
#'                     params=c(sim$pop[sim$pop$id==id, -which(names(sim$pop) %in% c('id','MTT'))]
#'                              , sigma=0.05, dose=100, duration=1),
#'                     format="data.frame")
#'   # drop 'traj' and shift 'time' to first column
#'   out <- out[,c('time',setdiff(colnames(out),c('time','traj')))]
#'   out$id <- paste("id",id,sep="")
#'   allout <- rbind(allout, out)
#' }
#' 
#' library(Hmisc)
#' allout <- upData(allout
#'                  , id = ordered(id, levels=paste("id",1:sim$N,sep=""))
#'                  , units=c(Prol="cells/mm^3", Tx.1="cells/mm^3",
#'                            Tx.2="cells/mm^3", Tx.3="cells/mm^3",
#'                            Circ="cells/mm^3", CircMin="cells/mm^3",
#'                            tNadir="hours", Cc="ng/mL", Cp="ng/mL",
#'                            time="hours"), print=FALSE)
#' 
#' library(tidyr)
#' cout <- gather(allout, key="Series", value="Concentration"
#' , Cc, Cp
#' , factor_key = TRUE)
#' 
#' label(cout$Concentration) <- "Drug Concentration"
#' 
#' # Figure 1 from reference [3]:
#' library(RColorBrewer)
#' xYplot(Concentration ~ time | id, group=Series
#'        , data=cout, subset=time<6
#'        , layout=c(5,NA)
#'        , type='l', as.table=TRUE
#'        , label.curves=FALSE
#'        , par.settings=list(superpose.line=list(lwd=2,col=brewer.pal(4,"PRGn")[c(1,4)]))
#'        , scales=list(y=list(log=TRUE, lim=c(10^-3,10^1)))
#'        , auto.key=list(points=FALSE, lines=TRUE, columns=2))
#'
#' mout <- gather(allout, key="Series", value="ANC"
#' , Prol, Tx.1, Tx.2, Tx.3, Circ
#' , factor_key = TRUE)
#' 
#' mout <- upData(mout
#'                , time = time/24
#'                , units = c(time="days")
#'                , print = FALSE)
#' 
#' # Figure 3 from citation [3]:
#' xYplot(ANC ~ time | id, group=Series
#'        , data=mout
#'        , layout=c(5,5)
#'        , type='l', as.table=TRUE
#'        , label.curves=FALSE
#'        , par.settings=list(superpose.line=list(lwd=2,col=brewer.pal(11,"RdYlBu")[c(1,3,4,8,10)]))
#'        , scales=list(y=list(log=TRUE, lim=c(100,15000), at=c(200, 500, 1000, 2000, 5000, 10000)))
#'        , auto.key=list(points=FALSE, lines=TRUE, columns=5))
#'
#' @export
Onoue.Friberg <-
function(N, cycle.length.days=21,
         data=data.frame(time=c(seq(0.0, 1.95, 0.05), # q3min for 2h, 
                                seq(2.0, cycle.length.days*24, 1.0)), # then hourly until Tmax
                         y=NA),
         delta.t=0.1){
  # Implement a lognormal measurement model
  pkpd.rmeas <- "
  y = rlnorm(log(Circ), sigma);
"
  pkpd.dmeas <- "
  lik = dlnorm(y, log(Circ), sigma, give_log);
"
  pkpd.rprior <- "
  // From Friberg et al 2002 (Table 4, row 1), taking sdlog ~= CV
  Circ0 = rlnorm(log(5050), 0.42);
  double MTT = rlnorm(log(89.3), 0.16);
  kTR = 4/MTT;
  gamma = rlnorm(log(0.163), 0.039);
  Emax  = rlnorm(log(83.9), 0.33);
  double dtx_mm = 0.808; // molar mass (g/mM) of docetaxel
  EC50  = rlnorm(log(7.17*dtx_mm), 0.50);
  // PK params from 2-compartment docetaxel model of Onoue et al (2016)
  CL = rlnorm(log(32.6), 0.295);
  Q  = rlnorm(log(5.34), 0.551);
  Vc = rlnorm(log(5.77), 0.1);
  Vp = rlnorm(log(11.0), 0.598);
"
  pkpd.dprior <- "
  // Parameter #4 setting 'log=1' returns log-density
  lik  = dlnorm(Circ0, log(5050), 0.42, 1);
  double MTT = 4/kTR;
  lik += dlnorm(MTT, log(89.3), 0.16, 1);
  lik += dlnorm(gamma, log(0.163), 0.039, 1);
  lik += dlnorm(Emax, log(83.9), 0.33, 1);
  double dtx_mm = 0.808; // molar mass (g/mM) of docetaxel
  lik += dlnorm(EC50, log(7.17*dtx_mm), 0.50, 1);
  lik += dlnorm(CL, log(32.6), 0.295, 1);
  lik += dlnorm(Q,  log(5.34), 0.551, 1);
  lik += dlnorm(Vc, log(5.77), 0.1, 1);
  lik += dlnorm(Vp, log(11.0), 0.598, 1);
  if (give_log != 1) lik = exp(lik);
"
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
  // We implement nadir-finding by integrating CircMin and tNadir into the state:
  int initialHump = CircMin == Circ_0; // Circ may overshoot Circ.0 initially
  int dropToNadir = DCirc < 0.0 && Circ < Circ_0 && (t-tNadir) < 1.0; // falling segment Circ.0-->nadir
  DCirc_0 = 0.0; // this pseudo-state merely remembers initial value
  DCircMin = dropToNadir ? DCirc : 0.0;
  DtNadir  = initialHump || dropToNadir ? 1.0 : 0.0;
"
  pkpd.step <- "
  double c2p = Q*( Cc - Cp ); // central-to-peripheral flux
  Cc += dt*( (dose/duration)*(t < duration ? 1.0 : 0.0)/Vc - (CL/Vc)*Cc - c2p/Vc );
  Cp += dt * c2p/Vp;
  // Myelosuppression model (Emax model, then dynamics per eqs 3-7 from Friberg et al 2002
  double Edrug = Emax*Cc/(Cc + EC50); // classic 'Emax model'
  Prol += dt*( (1-Edrug) * Prol * kTR * pow((Circ0 / Circ), gamma)  -  kTR * Prol );
  Tx_1 += dt * kTR * (Prol - Tx_1);
  Tx_2 += dt * kTR * (Tx_1 - Tx_2);
  Tx_3 += dt * kTR * (Tx_2 - Tx_3);
  double _DCirc = kTR * (Tx_3 - Circ);
  // We implement nadir-finding by integrating CircMin and tNadir into the state:
  int initialHump = CircMin == Circ_0; // Circ may overshoot Circ.0 initially
  int dropToNadir = _DCirc < 0.0 && Circ < Circ_0 && (t-tNadir) < 1.0; // falling segment Circ.0-->nadir
  Circ += dt * _DCirc;
  CircMin += dt*( dropToNadir ? _DCirc : 0.0 );
  tNadir  += dt*( initialHump || dropToNadir ? 1.0 : 0.0 );
"
  #Tmax <- cycle.length.days*24 # solve for full 21 days of 3-week cycle
  #df <- data.frame(time=c(seq(0.0, 1.95, 0.05), # q3min for 2h, 
  #                        seq(2.0, Tmax, 1.0)), # then hourly until Tmax
  #                 y=NA)
  # This 'rinit' function implements the rinit_spec design new in pomp v2:
  rinit_fun <- function(Circ0, t0=0, ...){
    state <- c(Cc = 0.0
               ,Cp = 0.0
               ,Prol = Circ0
               ,Tx.1 = Circ0
               ,Tx.2 = Circ0
               ,Tx.3 = Circ0
               ,Circ = Circ0
               ,Circ.0 = Circ0
               ,CircMin = Circ0
               ,tNadir = t0
    )
    state
  }
  # NULL defs to avoid R CMD check NOTEs 'no visible global function def':
  vectorfield <- euler <- parameter_trans <- NULL
  pkpd <- pomp(data = data
               , times="time", t0=0
               , skeleton=vectorfield(Csnippet(pkpd.skel))
               , rprocess = euler(step.fun = Csnippet(pkpd.step), delta.t = delta.t)
               , rmeasure = Csnippet(pkpd.rmeas)
               , dmeasure = Csnippet(pkpd.dmeas)
               , rinit = rinit_fun
               , rprior = Csnippet(pkpd.rprior)
               , dprior = Csnippet(pkpd.dprior)
               , statenames=c("Cc", "Cp", "Prol", "Tx.1", "Tx.2", "Tx.3", "Circ","Circ.0","CircMin","tNadir")
               , paramnames=c("Circ0","kTR","gamma","Emax","EC50","CL","Q","Vc","Vp"
                              ,"sigma","dose","duration"
               )
               , partrans = parameter_trans(log=c("Circ0","kTR","Emax","EC50",
                                                  "CL","Q","Vc","Vp","sigma",
                                                  "dose","duration"))
  )
  params.default <- c(Circ0=5050, kTR=4/89.3, gamma=0.163, Emax=83.9, EC50=7.17*0.808,
                      CL=32.6, Q=5.34, Vc=5.77, Vp=11.0, sigma=0.05, dose=50, duration=1.0
  )
  pop <- data.frame(id = integer()
                    ,Circ0= numeric()
                    ,gamma= numeric()
                    ,Emax = numeric()
                    ,EC50 = numeric()
                    ,CL = numeric()
                    ,Q  = numeric()
                    ,Vc = numeric()
                    ,Vp = numeric()
                    ,kTR = numeric()
  )
  for(id in 1:N){
    pop[id,'id'] <- id
    pop[id,-1] <- rprior(pkpd, params=params.default)[colnames(pop)[-1],]
  }
  
  pop$MTT <- 4/pop$kTR
  pop <- upData(pop
               #,MTT = 4/kTR
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
  sim$pkpd <- pkpd
  sim$params.default <- params.default
  sim$pop <- pop
  sim$N <- N
}
