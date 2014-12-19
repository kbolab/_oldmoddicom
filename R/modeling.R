#' Function for generating a series of outcomes related to a \code{dvhmatrix} class object
#' @param dvh a \code{dvhmatrix} class object
#' @param a factor for calculating EUD value according Niemierko formula
#' @param TD50 Dose that gives the 50\% probability of outcome
#' @param gamma50 Slope of dose-response curve at TD50 of administered dose
#' @export
DR.generate<-function(dvh, a = 1, TD50 = 45, gamma50 = 1.5) {
  
}

## Dose/Response according Lyman 1985, Kutcher and Burman 1989, Deasy 2000 (probit) ##
DR.Lyman <- function (TD50=45, gamma50=1.5, aa=1, diffdvh=NULL, dose=NULL) {
  # check the single choice between dvh matrix or dose series
  if (!is.null(dose) && !is.null(diffdvh)) stop("Select either a DVH or a point dose to calculate NTCP")
  # check the single choice between dvh matrix or dose series 
  if (!is.null(dose)) p <- pnorm(q=((dose - TD50)*gamma50*sqrt(2*pi))/TD50)
  if (!is.null(diffdvh)) p <- pnorm(q=((EUD(dvh.matrix=diffdvh, a=aa) - TD50)*gamma50*sqrt(2*pi))/TD50)
  return(p)
}

## Dose/Response according Goitein 1979 (logprobit) ##
DR.Goitein <- function (TD50=45, gamma50=1.5, aa=1, diffdvh=NULL, dose=NULL) {
  # check the single choice between dvh matrix or dose series
  if (!is.null(dose) && !is.null(diffdvh)) stop("Select either a DVH or a point dose to calculate NTCP")
  # check the single choice between dvh matrix or dose series 
  if (!is.null(dose)) p <- pnorm(q=(log(dose/TD50)*gamma50*sqrt(2*pi)))
  if (!is.null(diffdvh)) p <- pnorm(q=(log(EUD(dvh.matrix=diffdvh, a=aa)/TD50)*gamma50*sqrt(2*pi)))
  return(p)
}

## Dose/Response according Niemierko 1991 (loglogit) ##
DR.Niemierko <- function (TD50=45, gamma50=1.5, aa=1, diffdvh=NULL, dose=NULL) {
  # check the single choice between dvh matrix or dose series
  if (!is.null(dose) && !is.null(diffdvh)) stop("Select either a DVH or a point dose to calculate NTCP")
  # check the single choice between dvh matrix or dose series  
  if (is.vector(dose) && is.null(diffdvh)) p <- 1/(1+(TD50/dose)^(4*gamma50))
  if (is.null(dose) && is.matrix(diffdvh)) p <- 1/(1+(TD50/EUD(dvh.matrix=diffdvh, a=aa))^(4*gamma50))
  return(p)
}

## Dose/Response according Munro, Gilbert, Kallman 1992 (Poisson approximation) ##
DR.Munro <- function (TD50=45, gamma50=1.5, aa=1, diffdvh=NULL, dose=NULL) {
  # check the single choice between dvh matrix or dose series
  if (!is.null(dose) && !is.null(diffdvh)) stop("Select either a DVH or a point dose to calculate NTCP")
  # check the single choice between dvh matrix or dose series  
  if (is.vector(dose) && is.null(diffdvh)) p <- 2^(-(exp(exp(1)*gamma50*(1-dose/TD50))))
  if (is.null(dose) && is.matrix(diffdvh)) p <- 2^(-(exp(exp(1)*gamma50*(1-EUD(dvh.matrix=diffdvh, a=aa)/TD50))))
  return(p) 
}

## Dose/Response according Okunieff 1995 (logit) ##
DR.Okunieff <- function (TD50=45, gamma50=1.5, aa=1, diffdvh=NULL, dose=NULL) {
  # check the single choice between dvh matrix or dose series
  if (!is.null(dose) && !is.null(diffdvh)) stop("Select either a DVH or a point dose to calculate NTCP")
  # check the single choice between dvh matrix or dose series  
  if (is.vector(dose) && is.null(diffdvh)) p <- 1/(1+exp(4*gamma50*(1-(dose/TD50))))
  if (is.null(dose) && is.matrix(diffdvh)) p <- 1/(1+exp(4*gamma50*(1-(EUD(dvh.matrix=diffdvh, a=aa)/TD50))))
  return(p) 
}

## Dose/Response according Warkentin 2004 (Poisson) ##
DR.Warkentin <- function (TD50=45, gamma50=1.5, aa=1, diffdvh=NULL, dose=NULL) {
  # check the single choice between dvh matrix or dose series
  if (!is.null(dose) && !is.null(diffdvh)) stop("Select either a DVH or a point dose to calculate NTCP")
  # check the single choice between dvh matrix or dose series  
  if (is.vector(dose) && is.null(diffdvh)) p <- 0.5^(exp(2*gamma50/log(2)*(1-dose/TD50)))
  if (is.null(dose) && is.matrix(diffdvh)) p <- 0.5^(exp(2*gamma50/log(2)*(1-EUD(dvh.matrix=diffdvh, a=aa)/TD50)))
  return(p) 
}

## Dose/Response according Bentzen 1997 (log Poisson) ##
DR.Bentzen <- function (TD50=45, gamma50=1.5, aa=1, diffdvh=NULL, dose=NULL) {
  # check the single choice between dvh matrix or dose series
  if (!is.null(dose) && !is.null(diffdvh)) stop("Select either a DVH or a point dose to calculate NTCP")
  # check the single choice between dvh matrix or dose series  
  if (is.vector(dose) && is.null(diffdvh)) p <- 0.5^(TD50/dose)^(2*gamma50/log(2))
  if (is.null(dose) && is.matrix(diffdvh)) p <- 0.5^(TD50/EUD(dvh.matrix=diffdvh, a=aa))^(2*gamma50/log(2))
  return(p) 
}
