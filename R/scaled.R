#' Power-law scaling for doses
#' 
#' Implement an inverse power-law scaling for drug dose.
#' 
#' 
#' @param dose A numeric vector of doses
#' @param a A numeric exponent for power-law rescaling
#' @return A rescaled vector of doses
#' @author David C. Norris
scaled <-
function(dose, a=4.0) dose^(1/a)
