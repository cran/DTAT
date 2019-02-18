#' A seq method supporting custom-scaled plot axes.
#' 
#' This provides a \code{seq} method for class \code{function}, supporting a
#' natural axis scaling idiom.
#' 
#' 
#' @param scalefun A numeric function that will be invoked componentwise, and
#' so need not be vectorized)
#' @param from,to The starting and ending values of the sequence returned
#' @param length.out Desired length of the sequence
#' @param digits If non-NULL, returned value is rounded accordingly
#' @param \dots Unused; included for S3 generic/method consistency.
#' @return A numeric vector that (not considering the effect of any rounding
#' applied), becomes an arithmetic sequence after application of
#' \code{scalefun} to it. The initial and final elements of that vector are
#' \code{from} and \code{to}.
#' @author David C. Norris
#' @examples
#' 
#' # Provide evenly-spaced length-6 sequence from 100 to 1000,
#' # evenly spaced on a fourth-root scale:
#' seq(function(dose, a=4.0) dose^(1/a), from=100, to=1000, length.out=6, digits=0)
#' 
#' @importFrom stats uniroot
#' @export
seq.function <-
function(scalefun, from, to, length.out, digits=NULL, ...){
  x <- seq(from=scalefun(from), to=scalefun(to), length.out=length.out)
  y <- numeric(length.out)
  # Handle limiting case where from==to, to avert a uniroot error
  if(from == to){
    y1 <- uniroot(function(y) scaled(y)-x[1], c(from-1,from+1))$root
    return(rep(y1, length.out))
  }
  for(i in seq(length(y)))
    y[i] <- uniroot(function(y) scaled(y)-x[i], sort(c(from, to)))$root
  if(!is.null(digits))
    y <- round(y, digits=digits)
  y
}
