#' Function that creates DVH matrix, 1st column is the dose bins, other columns are volume steps
#' 
#' @param dvhnumber The number of DVHs to be generated
#' @param type The type of DVH to be generated: \code{"random"} generates all possible variation in dose distribution, 
#'        \code{"convex"} for DVH similar to PTV structures with good dose homogeneity to high levels,
#'        \code{"concave"} for DVH with low level of dose delivered to the volume as a spared structure,
#'        \code{"mix"} for DVH where dose distribution is averaged in the whole volume spectrum.
#' @description  Creates an object of class \code{dvhmatrix}.
#' @export
HI.generate<-function(dvhnumber, type=c("random","convex","concave","mix"), 
                      dvh.type=c("differential", "cumulative"), vol.distr=c("relative", "absolute"),
                      mindose = 0, maxdose = 75, dosebin = 0.5, mean.vol=200) {
  
}