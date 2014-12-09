#' Function that creates DVH matrix, 1st column is the dose bins, other columns are volume steps
#' 
#' @param dvh.number The number of DVHs to be generated
#' @param type The type of DVH to be generated: \code{random} generates all possible variation in dose distribution, 
#'        \code{convex} for DVH similar to PTV structures with good dose homogeneity to high levels,
#'        \code{concave} for DVH with low level of dose delivered to the volume as a spared structure,
#'        \code{mix} for DVH where dose distribution is averaged in the whole volume spectrum.
#' @param dvh.type The type of volume distribution to be created: \code{differential} or \code{cumulative}.
#' @param vol.distr Defines if the volume bins have to be divided by the total volume of the structure (\code{relative}) or not (\code{absolute}).
#' @param max.dose The upper bound of the dose distribution to be simulated.
#' @param dose.bin The dose bin in Gy to compute the value of volume parts in the final DVH.
#' @param volbin.side The value of the side of each cube (in mm) that builds the final volume of the simulated structure.
#' @param min.vol The minimum volume of the range to be simulated.
#' @param max.vol The maximum volume of the range to be simulated.
#' @description  Creates an object of class \code{dvhmatrix}.
#' @details Dose Volume Histograms (DVH) are the basic representation of how the radiation doses distribute inside structures.
#'          Usually they are shown in two forms: \code{differential} or \code{cumulative} being related each other because
#'          \code{cumultive} DVHs are the integral form of \code{"differential"} ones. This function provides a siumlated series
#'          of structures with the features defined by function parameters. The simulated series are useful for producing simulated
#'          Dose/Response models using the modeling functions defined in this package.
#'          class objects.
#' @references Van den Heuvela F. \emph{Decomposition analysis of differential dose volume histograms.} Med Phys. 2006 Feb;33(2):297-307. PubMed PMID: 16532934.
#' @return An object of \code{dvhmatrix} class.
#' @export
DVH.generate<-function(dvh.number, type=c("random","convex","concave","mix"), 
                      dvh.type=c("differential", "cumulative"), vol.distr=c("relative", "absolute"),
                      max.dose = 75, dose.bin = 0.5, volbin.side = 2.5, min.vol=180, max.vol=220) {
  if (min.vol<2) {
    warning("Minimum volumes under 2 cc aren't allowed, min.vol set to 2.")
    min.vol <- 2
  }
  if (max.dose <= 10) stop("max.dose must by higher than 10 Gy")
  # create the vector of volumes
  volumes <- runif(n = dvh.number, min = min.vol, max = max.vol)
  volbin.num <- round(volumes/((volbin.side/10)^3)) # number of bins for each volume
  # function for generating convex DVHs voxels series
  convex.dvh <- function(n) {
    mean.dose <- runif(n = 1, min = max.dose - (max.dose/4), max = max.dose - max.dose/8)
    sd.dose <- (max.dose - mean.dose)/2.5
    return(rnorm(n = n, mean = mean.dose, sd = sd.dose))
  }
  # function for generating concave DVHs voxels series
  concave.dvh <- function(n) {
    mean.dose <- runif(n = 1, min = 2, max = max(10, max.dose/4))
    sd.dose <- mean.dose/3
    result<-rnorm(n = n, mean = mean.dose, sd = sd.dose)
    # takes only the voxels with dose>=0
    result<-result[which(result>=0)]
    # adds more voxels to reach the expected number with random uniform distribution
    return(c(result, runif(n = n - length(result), min = 0, max = 2*mean.dose)))
  }
  # function for creating mix DVHs voxels series
  mix.dvh <- function(n) {
    contrib<-c(runif(n = 1, min = .2, max = .6), runif(n = 1, min = 0, max = .3))
    # proportions in contributions to final DVH
    contrib<-c(contrib[1], 1 - sum(contrib), contrib[2])
    part1<-concave.dvh(n = round(contrib[1] * n))
    part3<-convex.dvh(n = round(contrib[3] * n))
    part2<-rnorm(n = n - (length(part1) + length(part3)), 
                 mean = runif(n = 1, min = max.dose/12, max = max.dose - max.dose/5), 
                 sd = runif(n = 1, min = max.dose/15, max = max.dose/9))
    result<-c(part1, part2, part3)
    result<-result[which(result>=0)]    
    # compensate the negative values when available
    return(c(result, runif(n = n - length(result), min = 0, max = max.dose)))    
  }
  # function that creates random DVHs according the previous three given functions
  random.dvh <- function(n) {
    FUN <- sample(x = c(convex.dvh, concave.dvh, mix.dvh), size = 1, replace = T)
    return(FUN[[1]](n))
  }
  
  # voxels creation
  type <- match.arg(type)
  if (type=="convex") dose.voxels<-sapply(X = volbin.num, FUN = convex.dvh)
  if (type=="concave") dose.voxels<-sapply(X = volbin.num, FUN = concave.dvh)
  if (type=="mix") dose.voxels<-sapply(X = volbin.num, FUN = mix.dvh)
  if (type=="random") dose.voxels<-sapply(X = volbin.num, FUN = random.dvh)
  VolBin<-volbin.side^3/1000 # Volume Bin in cc
  # creates the vector of structures volumes
  volume<-unlist(lapply(X = dose.voxels, FUN = function(x) VolBin * length(x))) 
  result<-new("dvhmatrix")
  # creates the dvhmatrix object
  dvh.type<-match.arg(arg = dvh.type)
  vol.distr<-match.arg(arg = vol.distr)
  result@dvh.type<-dvh.type
  result@vol.distr<-vol.distr
  result@volume<-volume
  #browser()
  # creates the list of differential histograms
  hlist<-lapply(X = dose.voxels, FUN = hist, 
                breaks = seq(from = 0, to = max(c(unlist(lapply(X = dose.voxels, FUN = max))) + dose.bin * 4, 
                             max.dose), by = dose.bin), plot = FALSE)
  if (dvh.type=="differential") {
    result@dvh<-cbind(hlist[[1]]$mids, sapply(X = hlist, FUN = function(x) cbind(x$counts * VolBin), simplify = TRUE, USE.NAMES = FALSE))
    if (vol.distr=="relative") for (n in 2:ncol(result@dvh)) result@dvh[,n]<-result@dvh[,n]/result@volume[n-1]
  }
  if (dvh.type=="cumulative") {
    result@dvh<-cbind(hlist[[1]]$mids, sapply(X = hlist, FUN = function(x) cbind(x$counts * VolBin), simplify = TRUE, USE.NAMES = FALSE))
    if (vol.distr=="relative") result<-DVH.diff.to.cum(dvh = result, relative = TRUE) else result<-DVH.diff.to.cum(dvh = result, relative = FALSE)
  }
  return(result)
}

#' Function that converts differential DVHs into cumulative ones
#' 
#' @param dvh Either an object of class \code{dvhmatrix} or \code{matrix} type.
#' @param relative If \code{TRUE} the structure volume bins are divided by the total volume.
#' @description Function that converts an object of class \code{dvhmatrix} or a simple \code{matrix} that 
#'              represents a differential DVH into a cumulative one.
#' @return Either an object of class \code{dvhmatrix} or a \code{matrix} according the argument \code{dvh}.
#' @export
DVH.diff.to.cum <- function(dvh, relative=TRUE) {
  if ((!is.matrix(dvh))&&(class(dvh)!="dvhmatrix")) stop("dvh MUST be either an bobject of class dvhmatrix or a matrix")
  if (class(dvh)=="dvhmatrix") dvh.matrix<-dvh@dvh else dvh.matrix<-dvh  
  dvh.size <- dim(dvh.matrix)
  DVHList <- matrix(nrow=dvh.size[1] + 1, ncol=dvh.size[2])   # create the matrix of cumulative DVHs     
  for (m in 2:dvh.size[2]) {                                  # loop for columns (volumes) 
    total.volume <- sum(dvh.matrix[,m])                       # calculate the total volume
    if (relative==TRUE) {
      for (n in 1:dvh.size[1]) {                                    # loop for rows
        DVHList[n+1, m] <- (total.volume - sum(dvh.matrix[c(1:n),m]))/total.volume # elements of the matrix as relative volume
      }
      DVHList[1,m] <- 1 # first element is 1 by default if relative==TRUE
    } else {
      for (n in 1:dvh.size[1]) {                                    # loop for rows
        DVHList[n+1, m] <- total.volume - sum(dvh.matrix[c(1:n),m]) # elements of the matrix as relative volume
      }
      DVHList[1,m] <- total.volume  # first element is total volume by default
    }                                       
  }
  DVHList[1,1]<-0
  #  DVHList[2:nrow(DVHList),1]<-dvh.matrix[,1]
  for (n in 2:nrow(DVHList)) DVHList[n,1]<-2*dvh.matrix[n-1,1]-DVHList[n-1,1]
  if (class(dvh)=="dvhmatrix") {
    dvh@dvh<-DVHList
    dvh@dvh.type<-"cumulative"
    if (relative==TRUE) dvh@vol.distr<-"relative" else dvh@vol.distr<-"absolute"
    return(dvh)
  } else return(DVHList)
}



extractDVH<-function(x, maxDose=NULL, stepDose=.25, dvh.type=c("differential","cumulative"), 
                     vol.distr=c("relative","absolute"), createObj=FALSE, VolBin=0.015625) {
  # default VolBin is given in cm3
  TotalVol<-length(x)*VolBin
  dvh.type=match.arg(dvh.type)
  vol.distr=match.arg(vol.distr)
  if (is.null(x=maxDose)) maxDose<-max(x)+stepDose*4
  h<-hist(x=x, breaks=seq(from=0, to=maxDose,by=stepDose), plot=FALSE)
  diff<-cbind(h$mids, h$density/sum(h$density)*TotalVol)
  # return matrix without dvhmatrix class structure
  if ((dvh.type=="differential") && (vol.distr=="absolute"))  final.matrix<-diff
  if ((dvh.type=="differential") && (vol.distr=="relative"))  final.matrix<-rel.diff.dvh(diff)
  if ((dvh.type=="cumulative")   && (vol.distr=="absolute"))  final.matrix<-cum.dvh(dvh.matrix=diff, relative=FALSE)
  if ((dvh.type=="cumulative")   && (vol.distr=="relative"))  final.matrix<-cum.dvh(dvh.matrix=diff, relative=TRUE)
  if (createObj=="FALSE") return(final.matrix)
  # return matrix within a dvhmatrix object
  if (createObj==TRUE) return(new("dvhmatrix", dvh=final.matrix, dvh.type=dvh.type, vol.distr=vol.distr, volume=TotalVol))  
}