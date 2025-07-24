#' Simulated \sQuote{3+3/PC} dose-titration study from bioRxiv paper no. 240846
#'
#' This is a length-10 list of data frames, summarizing the simulated trial
#' from this paper, at the end of periods 1, 2, ..., 10. This structure reflects
#' an awkward S3 implementation that package DTAT v0.3 reimplemented using S4.
#' This data set is retained to support regression tests.
#'
#'
#' @name de.bioRxiv.240846
#' @docType data
#' @format A length-10 list of data frames, each with the following columns:
#'
#' \describe{
#'   \item{id}{Participant identifier}
#'   \item{period}{DLT assessment period, numbered consecutively from 1}
#'   \item{dose}{Dose level, numbered consecutively starting from 1}
#'   \item{dlt}{A logical indicator: did this this participant experience
#'     a DLT during this period?}
#' }
#' @details A \code{stop.esc} attribute is attached to data frames in this list,
#' indicating when escalation stopped during the simulated trial.
#'
#' @references Norris DC. Precautionary Coherence Unravels Dose Escalation
#' Designs. \emph{bioRxiv}. December 2017:240846. \doi{10.1101/240846}.
#' \url{https://www.biorxiv.org/content/10.1101/240846v1}
#' @keywords datasets
#' @examples
#'
#' data(de.bioRxiv.240846)
#' # Demonstrate that the new S4 3+3/PC implementation reproduces the
#' # simulated trial from the paper:
#' set.seed(2017)
#' CV <- 0.7; mean_mtd <- 1.0
#' shape <- CV^-2; scale <- mean_mtd/shape
#' trial <- new("DE", doses=0.25 * 1.4^(0:6),
#'              MTDi=rgamma(24, shape=shape, scale=scale),
#'              units="mg")
#' trial <- titration(trial, periods=10)
#' stopifnot(all(trial@data == de.bioRxiv.240846[[10]]))
#' stopifnot(trial@stop_esc == attr(de.bioRxiv.240846[[10]],'stop.esc'))
#'
NULL





#' @name DTAT-package
#' @title Dose Titration Algorithm Tuning: a Framework for Dose Individualization
#'  in Drug Development
#' @description Dose Titration Algorithm Tuning (DTAT) is a methodologic framework
#' allowing dose individualization to be conceived as a continuous learning process
#' that begins in early-phase clinical trials and continues throughout drug development,
#' on into clinical practice.
#' This package includes code that researchers may use to reproduce or extend key results
#' of the DTAT research programme, plus tools for trialists to design and simulate a
#' '3+3/PC' dose-finding study. Please see Norris (2017a) \doi{10.12688/f1000research.10624.3}
#' and Norris (2017c) \doi{10.1101/240846}.
#' @aliases DTAT-package DTAT
#' @author David C. Norris
#'
#' @references
#' 1. Norris DC. Dose Titration Algorithm Tuning (DTAT) should supersede
#'    \sQuote{the} Maximum Tolerated Dose (MTD) in oncology dose-finding trials.
#'    \emph{F1000Research}. 2017;6:112. \doi{10.12688/f1000research.10624.3}.
#'    \url{https://f1000research.com/articles/6-112/v3}
#'
#' 2. Norris DC. Costing \sQuote{the} MTD. \emph{bioRxiv}. August 2017:150821.
#'    \doi{10.1101/150821}.
#'    \url{https://www.biorxiv.org/content/10.1101/150821v3}
#'
#' 3. Norris DC. Precautionary Coherence Unravels Dose Escalation Designs.
#'    \emph{bioRxiv}. December 2017:240846. \doi{10.1101/240846}.
#'    \url{https://www.biorxiv.org/content/10.1101/240846v1}
#'
#' 4. Norris DC. One-size-fits-all dosing in oncology wastes money, innovation
#'    and lives. \emph{Drug Discov Today}. 2018;23(1):4-6.
#'    \doi{10.1016/j.drudis.2017.11.008}.
#'    \url{https://precisionmethods.guru/DTAT/Norris%20(2018)%20One-size-fits-all%20dosing%20in%20oncology%20wastes%20money,%20innovation%20and%20lives.pdf}
#'
#' 5. Norris DC. Costing \sQuote{the} MTD ... in 2-D. \emph{bioRxiv}. July 2018:370817.
#'    \doi{10.1101/370817}.
#'    \url{https://www.biorxiv.org/content/10.1101/370817v1}
#'
#' @importFrom pomp pomp rprior trajectory Csnippet
#' @import survival
#' @import methods
#' @keywords internal
"_PACKAGE"





#' Precomputed neutrophil-guided chemotherapy dose titration for 1000 simulated
#' subjects.
#'
#' This dataset is provided to support fast reproduction of a forthcoming
#' pharmacoeconomic paper that includes examination of the empirical
#' distribution of MTDi in N=1000 simulated subjects.
#'
#' Running the examples interactively, you can verify the reproducibility of
#' this dataset. (That demo is included in a \code{donttest} block to spare the
#' CRAN servers.)
#'
#' @name dtat1000
#' @docType data
#' @format A data frame showing end-of-cycle state of neutrophil-guided dose
#' titration for 1000 simulated subjects, across 10 cycles of chemotherapy.
#' \describe{
#'   \item{cycle}{Cycle number 1..10}
#'   \item{id}{Subject identifiers; an ordered factor with levels
#' \code{id1} < \dots{} < \code{id1000}}
#'   \item{Cc}{Central-compartment drug concentration}
#'   \item{Cp}{Peripheral-compartment drug concentration}
#'   \item{Prol}{Progenitor cells in proliferating compartment of
#'               Friberg et al. (2002) model}
#'   \item{Tx.1}{Transit compartment 1}
#'   \item{Tx.2}{Transit compartment 1}
#'   \item{Tx.3}{Transit compartment 1}
#'   \item{Circ}{Concentration (cells/mm^3) of circulating neutrophils}
#'   \item{dose}{Dose of 1-hour infusion administered this cycle}
#'   \item{CircMin}{Neutrophil nadir (cells/mm^3)}
#'   \item{tNadir}{Time (days) of neutrophil nadir}
#'   \item{scaled.dose}{Fourth root of dose}
#'   \item{time}{Time (weeks) of dose administration}
#' }
#' @references 1. Norris DC. Dose Titration Algorithm Tuning (DTAT) should
#' supersede \sQuote{the} Maximum Tolerated Dose (MTD) in oncology dose-finding
#' trials. \emph{F1000Research}. 2017;6:112. \doi{10.12688/f1000research.10624.3}.
#' \url{https://f1000research.com/articles/6-112/v3}
#'
#' 2. Norris DC. Costing \sQuote{the} MTD. \emph{bioRxiv}. August 2017:150821.
#' \doi{10.1101/150821}.
#' \url{https://www.biorxiv.org/content/10.1101/150821v3}
#' @keywords datasets
#' @examples
#'
#' data(dtat1000)
#' # 1. Extract the N final doses, assuming convergence by the tenth course
#' MTD_i <- with(dtat1000, dose[time==27])
#' MTD_i <- MTD_i[MTD_i < 5000] # Exclude few outliers
#' # 2. Do a kernel density plot
#' library(Hmisc)
#' library(latticeExtra)
#' hist <- histogram(~MTD_i, breaks=c(0,100,200,300,400,600,900,1500,2500,4000,5000)
#'                   , xlab=expression(MTD[i]))
#' approx <- data.frame(mtd_i=seq(0, 5000, 10))
#' approx <- upData(approx,
#'                  gamma = dgamma(mtd_i, shape=1.75, scale=200))
#' dist <- xyplot(gamma ~ mtd_i, data=approx, type='l', col='black', lwd=2)
#' library(grid)
#' hist + dist
#' grid.text(expression(MTD[i] %~%
#'                      paste("Gamma(", alpha==1.75, ", ", beta==1/200,")"))
#'          , x=unit(0.5,"npc")
#'          , y=unit(0.75,"npc")
#'          )
#' ## A very long repro, which a user of this package may well wish to verify
#' ## by running the examples interactively, although it takes many minutes
#' ## to compute.  (Enclosed in a dontest block to avoid overburdening CRAN.)
#' \donttest{
#' # Demonstrate close reproduction of original titration (the titration takes many minutes!)
#' set.seed(2016)
#' library(pomp)
#' Onoue.Friberg(N=1000)
#' # This titration may take an hour to run ...
#' chemo <- titrate(doserange = c(50, 3000),
#'                  dta=newton.raphson(dose1 = 100,
#'                                     omega = 0.75,
#'                                     slope1 = -2.0,
#'                                     slopeU = -0.2)
#' )
#'
#' dtat1k <- upData(chemo$course
#'                 , time = 3*(cycle-1)
#'                 , labels = c(time="Time")
#'                 , units = c(time="weeks")
#'                 , print = FALSE)
#'
#' c10dose1k <- subset(dtat1k, cycle==10)$scaled.dose
#' c10dose1000 <- subset(dtat1000, cycle==10)$scaled.dose
#' stopifnot(0.999 < cor(c10dose1k, c10dose1000))
#' }
#'
NULL





#' Environment for simulation global variables.
#'
#' To simplify the code of package DTAT, as well as client tasks, this exported
#' environment contains a handful of global variables useful for the
#' simulations.
#'
#' Global variables maintained within environment \code{sim} are:
#' \enumerate{
#'   \item \code{pkpd}: The population PK/PD model to be simulated.
#'   \item \code{pop}: A sample drawn from the population model.
#'   \item \code{N}: Restricts simulation to first \code{N} subjects in
#'     \code{pop}.
#'   \item \code{params.default}: Default parameters.
#' }
#'
#' @export
#' @name sim
#' @docType data
#' @keywords datasets
#' @examples
#'
#' # Even when nrow(pop) is large, one may easily restrict
#' # time-consuming simulations to pop[1:N,], as follows:
#' sim$N <- 25
#' # Now perform simulation work
#' \dontrun{
#' titrate(...)
#' }
#'
NULL



