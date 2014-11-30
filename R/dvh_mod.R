#' Function that creates DVH matrix, 1st column is the dose bins, other columns are volume steps
#' 
#' @param dvh.number The number of DVHs to be generated
#' @param type The type of DVH to be generated: \code{"random"} generates all possible variation in dose distribution, 
#'        \code{"convex"} for DVH similar to PTV structures with good dose homogeneity to high levels,
#'        \code{"concave"} for DVH with low level of dose delivered to the volume as a spared structure,
#'        \code{"mix"} for DVH where dose distribution is averaged in the whole volume spectrum.
#' @param dvh.type The type of volume distribution to be created: \code{"differential"} or \code{"cumulative"}.
#' @param vol.distr Defines if the volume bins have to be divided by the total volume of the structure (\code{"relative"}) or not (\code{"absolute"}).
#' @param max.dose The upper bound of the dose distribution to be simulated.
#' @param dose.bin The dose bin in Gy to compute the value of volume parts in the final DVH.
#' @param volbin.side The value of the side of each cube (in mm) that builds the final volume of the simulated structure.
#' @param min.vol The minimum volume of the range to be simulated.
#' @param max.vol The maximum volume of the range to be simulated.
#' @description  Creates an object of class \code{dvhmatrix}.
#' @details Dose Volume Histograms (DVH) are the basic representation of how the radiation doses distribute inside structures.
#'          Usually they are shown in two forms: \code{"differential"} or \code{"cumulative"} being related each other because
#'          \code{"cumultive"} DVHs are the integral form of \code{"differential"} ones. This function provides a siumlated series
#'          of structures with the features defined by function parameters. The simulated series are useful for producing simulated
#'          Dose/Response models using the modeling functions defined in this package.
#'          class objects.
#' @references Van den Heuvela F. \emph{Decomposition analysis of differential dose volume histograms.} Med Phys. 2006 Feb;33(2):297-307. PubMed PMID: 16532934.
#' @return An object of \code{"dvhmatrix"} class.
#' @export
HI.generate<-function(dvh.number, type=c("random","convex","concave","mix"), 
                      dvh.type=c("differential", "cumulative"), vol.distr=c("relative", "absolute"),
                      max.dose = 75, dose.bin = 0.5, volbin.side = 2.5, min.vol=180, max.vol=220) {
  if (min.vol<2) {
    warning("Minimum volumes under 2 cc aren't allowed, min.vol set to 2.")
    min.vol<-2
  }
  # create the vector of volumes
  volumes<-runif(n = dvh.number, min = min.vol, max = max.vol)
  volbin.num<-round(volumes/((volbin.side/10)^3)) # number of bins for each volume
  convex.dvh<-function(n) {
    return(rnorm(n = n, mean = runif(n = 1, min = max.dose - (max.dose/10), max = max.dose), sd = 2.5))
  }
  return(sapply(X = volbin.num, FUN = convex.dvh))
}