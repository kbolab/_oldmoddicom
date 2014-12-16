#' Function that creates a \code{dvhmatrix} class object
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
#' @examples # creates a dvhmatrix object with 200 differential - relative histograms
#' a<-DVH.generate(dvh.number=200, dvh.type="differential", vol.distr="relative")
#' 
#' # creates a dvhmatrix object containing 150 cumulative - absolute histograms
#' # maximum dose in simulated series is 60 Gy
#' b<-DVH.generate(dvh.number=150, dvh.type="cumulative", vol.distr="absolute", max.dose=60)
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
    mean.dose <- runif(n = 1, min = max.dose - (max.dose/3), max = max.dose - max.dose/6)
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
    contrib<-c(runif(n = 1, min = .05, max = .4), runif(n = 1, min = .05, max = .4))
    # proportions in contributions to final DVH
    contrib<-c(contrib[1], 1 - sum(contrib), contrib[2])
    part1<-concave.dvh(n = round(contrib[1] * n))
    part3<-convex.dvh(n = round(contrib[3] * n))
    part2<-rnorm(n = n - (length(part1) + length(part3)), 
                 mean = runif(n = 1, min = max.dose/15, max = max.dose - max.dose/8), 
                 sd = runif(n = 1, min = max.dose/15, max = max.dose/10))
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
  result@dvh.type<-"differential" # default value, corrected by DVH.diff.to.cum if dvh.type = "cumulative"
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
  if ((!is.matrix(dvh))&&(class(dvh)!="dvhmatrix")) stop("dvh MUST be either an object of class dvhmatrix or a matrix")
  if (class(dvh)=="dvhmatrix") dvh.matrix<-dvh@dvh else dvh.matrix<-dvh  
  if (class(dvh)=="dvhmatrix") if (dvh@dvh.type=="cumulative") {
    if (relative==TRUE) dvh<-DVH.relative(dvh = dvh)
    return(dvh)
  }
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


#' Function that converts cumulative DVHs into differential ones
#' 
#' @param dvh Either an object of class \code{dvhmatrix} or \code{matrix} type.
#' @param relative If \code{TRUE} the structure volume bins are divided by the total volume.
#' @description Function that converts an object of class \code{dvhmatrix} or a simple \code{matrix} that 
#'              represents a cumulative DVH into a differential one.
#' @return Either an object of class \code{dvhmatrix} or a \code{matrix} according the argument \code{dvh}.
#' @export
DVH.cum.to.diff <- function(dvh, relative=TRUE) {
  if ((!is.matrix(dvh))&&(class(dvh)!="dvhmatrix")) stop("dvh MUST be either an object of class dvhmatrix or a matrix")
  if (class(dvh)=="dvhmatrix") dvh.matrix<-dvh@dvh else dvh.matrix<-dvh
  if (class(dvh)=="dvhmatrix") if (dvh@dvh.type=="differential") {
    if (relative==TRUE) dvh<-DVH.relative(dvh)
    return(dvh)
  }
  dvh.size <- dim(dvh.matrix)
  DVHList <- matrix(nrow=dvh.size[1] - 1, ncol=dvh.size[2])   # create the matrix of differential DVHs
  for (m in 2:dvh.size[2]) {                                  # loop for columns (volumes)
    total.volume<-dvh.matrix[1,m]                             # set total volume to 1st element of the cumulative DVH
    for (n in 2:dvh.size[1]) {                                # loop for rows
      if (relative==TRUE) DVHList[n-1,m] <- (dvh.matrix[n-1,m]-dvh.matrix[n,m])/total.volume # calculate differential volume
      else DVHList[n-1,m] <- (dvh.matrix[n-1,m]-dvh.matrix[n,m])   # differential volume for absolute distribution
    }    
  }       
  for (m in 2:dvh.size[1]) {                                   # loop for column of dose values
    DVHList[m-1,1]<-dvh.matrix[m-1,1]+(dvh.matrix[m,1]-dvh.matrix[m-1,1])/2 # sets the dose values
  }
  if (class(dvh)=="dvhmatrix") {
    dvh@dvh<-DVHList
    dvh@dvh.type<-"differential"
    if (relative==TRUE) dvh@vol.distr<-"relative" else dvh@vol.distr<-"absolute"
    return(dvh)
  } else return(DVHList)
}

#' Extracts a DVH from a vector of dose bins
#' @param x vector of dose bins.
#' @param max.dose Upper dose bound for limiting the DVH computation.
#' @param dose.bin The dose bin for DVH computation in Gy.
#' @param dvh.type The type of volume distribution to be created: \code{differential} or \code{cumulative}.
#' @param vol.distr Defines if the volume bins have to be divided by the total volume of the structure (\code{relative}) or not (\code{absolute}).
#' @param createObj if \code{TRUE} returns a \code{dvhmatrix} class object.
#' @param volbin.side The value of the side of each cube (in mm) that builds the final volume of the structure.
#' @description Function that given a vector of dose bins (either a vector got from sampling a 3D mesh or a vector of
#' simple values of dose) extracts the DVH that summarizes that vector.
#' @return Either an object of class \code{dvhmatrix} or a \code{matrix} according the argument \code{dvh} type.
#' @examples # simulate a vector of dose bins
#' doses <- c(rnorm(n = 10000, mean = 45, sd = 3), rnorm(n = 7000, mean = 65, sd = 2.5))
#' 
#' # creates a dvhmatrix class object
#' DVH<-DVH.extract(x = doses)
#' @export
DVH.extract<-function(x, max.dose=NULL, dose.bin=.25, dvh.type=c("differential","cumulative"), 
                     vol.distr=c("relative","absolute"), createObj=TRUE, volbin.side=2.5) {
  # default VolBin is given in cm3
  VolBin<-(volbin.side/10)^3 
  TotalVol<-length(x)*VolBin
  dvh.type=match.arg(dvh.type)
  vol.distr=match.arg(vol.distr)
  if (is.null(x=max.dose)) max.dose<-max(x)+dose.bin*4
  h<-hist(x=x, breaks=seq(from=0, to=max.dose,by=dose.bin), plot=FALSE)
  diff<-cbind(h$mids, h$density/sum(h$density)*TotalVol)
  # return matrix without dvhmatrix class structure
  if ((dvh.type=="differential") && (vol.distr=="absolute"))  final.matrix<-diff
  if ((dvh.type=="differential") && (vol.distr=="relative"))  final.matrix<-rel.diff.dvh(diff)
  if ((dvh.type=="cumulative")   && (vol.distr=="absolute"))  final.matrix<-cum.dvh(dvh.matrix=diff, relative=FALSE)
  if ((dvh.type=="cumulative")   && (vol.distr=="relative"))  final.matrix<-cum.dvh(dvh.matrix=diff, relative=TRUE)
  if (createObj==FALSE) return(final.matrix)
  # return matrix within a dvhmatrix object
  if (createObj==TRUE) return(new("dvhmatrix", dvh=final.matrix, dvh.type=dvh.type, vol.distr=vol.distr, volume=TotalVol))  
}

#' mean of DVHs
#' @param dvh A \code{dvhmatrix} object
#' @description Function that gives the value of the mean dose of a \code{dvhmatrix} class object
#' @return a vector with the means of doses in the \code{dvhmatrix} object
#' @export
#' @examples # generate a dataset of DVHs
#' a<-DVH.generate(dvh.number = 100)
#' m<-DVH.mean(a)
DVH.mean<-function(dvh)  {
  if (class(dvh)!="dvhmatrix") stop("dvh MUST be a dvhmatrix class object")
  if (dvh@dvh.type=="cumulative") dvh<-DVH.cum.to.diff(dvh = dvh)
  return(DVH.eud(dvh = dvh, a = 1)) # use a = 1 that corresponds to mean dose
}

#' Converts absolute \code{dvhmatrix} class objects to relative
#' @param dvh A \code{dvhmatrix} class object
#' @description Function that converts an object of \code{dvhmatrix} class where DVH are stored in absoulte mode into
#' another \code{dvhmatrix} class object where DVH are stored in relative mode.
#' @return A \code{dvhmatrix} class object.
#' @export
#' @examples # generate a dataset of absolute DVHs
#' a<-DVH.generate(dvh.number = 10, vol.distr = "absolute")
#' b<-DVH.relative(dvh = a)
DVH.relative<-function(dvh) {
  if (dvh@vol.distr=="absolute")
    for (n in 1:(ncol(dvh@dvh) - 1)) dvh@dvh[,n+1]<-dvh@dvh[,n+1]/dvh@volume[n]
  dvh@vol.distr<-"relative"
  return(dvh)
}

#' Converts relative \code{dvhmatrix} class objects to absolute
#' @param dvh A \code{dvhmatrix} class object
#' @description Function that converts an object of \code{dvhmatrix} class where DVH are stored in relative mode into
#' another \code{dvhmatrix} class object where DVH are stored in absolute mode.
#' @return A \code{dvhmatrix} class object.
#' @export
#' @examples # generate a dataset of relative DVHs
#' a<-DVH.generate(dvh.number = 10, vol.distr = "relative")
#' b<-DVH.absolute(dvh = a)
DVH.absolute<-function(dvh) {
  if (dvh@vol.distr=="relative")
    for (n in 1:(ncol(dvh@dvh) - 1)) dvh@dvh[,n+1]<-dvh@dvh[,n+1]*dvh@volume[n]
  dvh@vol.distr<-"absolute"
  return(dvh)
}

#' Calculates Equivalent Uniform Dose for a \code{dvhmatrix} object
#' @param dvh A \code{dvhmatrix} class object
#' @param a Factor for parallel-serial correlation in radiobiological response
#' @description Function that calculates the value of Equivalent Uniform Dose (EUD) for a \code{dvhmatrix} object.
#' @return A vector containing the values of EUD(s) for the given DVH(s)
#' @export
#' @useDynLib moddicom
DVH.eud<-function(dvh, a = 1) {
  dvh<-DVH.cum.to.diff(dvh = dvh, relative = TRUE)
  Ncol<-ncol(dvh@dvh) - 1
  Nrow<-nrow(dvh@dvh)
  ceud<-rep.int(x = 0, times = Ncol) 
  dosebin<-dvh@dvh[,1]
  volumebin<-dvh@dvh[,2:(Ncol + 1)]
  result<-.C("cEUD", as.double(dosebin), as.double(volumebin), as.double(a), as.integer(Nrow), 
             as.integer(Ncol), as.double(ceud))
  return(result[[6]])
}

#' Merge two different \code{dvhmatrix} class objects into one
#' @param receiver The \code{dvhmatrix} object that will receive the \code{addendum} object.
#' @param addendum The \code{dvhmatrix} object to be merged with \code{addendum}.
#' @description This function can e used to create a \code{dvhmatrix} class object from two different
#'              ones. The slots of the resulting object will be the same in the \code{receiver} one,
#'              so some convertions can be automatically realized to create the homogeneous final result.
#' @return A \code{dvhmatrix} class object.
#' @export
#' @examples # creates two different dvhmatrx objects
#' a<-DVH.generate(dvh.number = 100, dvh.type="differential", vol.distr = "relative")
#' b<-DVH.generate(dvh.number = 100, dvh.type="cumulative", vol.distr = "absolute")
#' ab<-DVH.merge(receiver = a, addendum = b)
DVH.merge<-function(receiver=NULL, addendum=NULL) {
  # check dvhmatrix class
  if ((class(receiver)!="dvhmatrix") || (class(addendum)!="dvhmatrix")) 
    stop("BOTH dvh.to.add AND destination MUST be dvhmatrix class objects")
  # creates list of dose bins
  dbin<-list(receiver=receiver@dvh[3,1]-receiver@dvh[2,1], addendum=addendum@dvh[3,1]-addendum@dvh[2,1])
  # detect the volume distribution  
  if (receiver@vol.distr=="relative") rel<-TRUE else rel<-FALSE  # relative state of receiver
  
  # STEP (1): identify conversion function if DVH types are different and convert addendum according receiver type
  if (receiver@dvh.type!=addendum@dvh.type) {
    if (receiver@dvh.type=="cumulative") conv.funct.name<-"DVH.diff.to.cum"
    if (receiver@dvh.type=="differential") conv.funct.name<-"DVH.cum.to.diff"
    conv.funct<-match.fun(FUN=conv.funct.name)  
    # correct the addendum if receiver has absolute distribution
    if ((rel==FALSE) && (addendum@vol.distr=="relative")) 
      for (n in 2:ncol(addendum@dvh)) addendum@dvh[,n]<-addendum@dvh[,n]*addendum@volume[n-1]
    addendum@dvh<-conv.funct(dvh=addendum@dvh, relative=rel)    
  }
  
  # STEP (2): check the dose-bin for interpolating addendum or receiver according the lowest dbin
  if (dbin$receiver!=dbin$addendum) {
    warning("Different dose bins between receiver and addendum: linear interpolation performed.")
    # create temp cumulative dvh, needed for avoiding interpolation artifacts
    if (receiver@dvh.type=="differential") {
      temp.receiver<-DVH.diff.to.cum(dvh=receiver@dvh, relative=rel)
      temp.addendum<-DVH.diff.to.cum(dvh=addendum@dvh, relative=rel)
    } else {
      temp.receiver<-receiver@dvh
      temp.addendum<-addendum@dvh
    }    
    if (dbin$receiver<dbin$addendum) {
      # if receiver's dbin is smaller then interpolates the addendum according the Dose column of receiver
      Dout<-seq(from=0, to=max(temp.receiver[,1], temp.addendum[,1]), by=dbin$receiver)
      # temp matrix
      temp.out<-matrix(nrow=length(Dout), ncol=ncol(temp.addendum))
      temp.out[,1]<-Dout
      for (n in 2:ncol(temp.addendum)) temp.out[,n]<-approx(x=temp.addendum[,1], y=temp.addendum[,n], xout=Dout, yright=0)$y
      if (receiver@dvh.type=="differential") {
        temp.addendum<-DVH.cum.to.diff(dvh=temp.out, relative=rel)
        temp.receiver<-receiver@dvh
      }
      else temp.addendum<-temp.out
    }  
    if (dbin$receiver>dbin$addendum) {
      # if addendum's dbin is smaller then interpolates the receiver according the Dose column of addendum
      Dout<-seq(from=0, to=max(temp.receiver[,1], temp.addendum[,1]), by=dbin$addendum)
      # temp matrix
      temp.out<-matrix(nrow=length(Dout), ncol=ncol(temp.receiver))
      temp.out[,1]<-Dout
      for (n in 2:ncol(temp.receiver)) temp.out[,n]<-approx(x=temp.receiver[,1], y=temp.receiver[,n], xout=Dout, yright=0)$y
      if (receiver@dvh.type=="differential") {
        temp.receiver<-DVH.cum.to.diff(dvh=temp.out, relative=rel)
        temp.addendum<-addendum@dvh
      }
      else temp.receiver<-temp.out
    } 
  } else {
    temp.receiver<-receiver@dvh
    temp.addendum<-addendum@dvh
  }  
  
  # STEP (3): check if the two dvhs have equal length
  if (nrow(temp.receiver)>nrow(temp.addendum)) {
    for (n in 1:(nrow(temp.receiver)-nrow(temp.addendum))) temp.addendum<-rbind(temp.addendum, 0)
    temp.addendum[,1]<-temp.receiver[,1]
  } 
  if (nrow(temp.receiver)<nrow(temp.addendum)) {
    for (n in 1:(nrow(temp.addendum)-nrow(temp.receiver))) temp.receiver<-rbind(temp.receiver, 0)
    temp.receiver[,1]<-temp.addendum[,1]
  }  
  
  # STEP (4): joins the two dvhs
  result.DVH<-cbind(temp.receiver, temp.addendum[,2:ncol(temp.addendum)])
  return(new("dvhmatrix", dvh=result.DVH, dvh.type=receiver@dvh.type, vol.distr=receiver@vol.distr, 
             volume=c(receiver@volume, addendum@volume)))
}


# method for class dvhmatrix
# enhance: is a vector of number of DVHs to be enhanced in the plot
# mean: plots the mean DVH overlapped to the background
# elements: if equal to a string "all" (default value) plots all the DVHs, if it is a vector plots only the
#   chosen DVHs
#' method for plotting \code{dvhmatrix} class objects
#' @description This method overrides the default method \code{plot} and handles the plots of \code{dvhmatrix}
#'              class objects.
#' @param x The \code{dvhmatrix} object to be plotted.
#' @param y missing argument.
#' @param elements A vector containing the indeces of the DVHs to be plotted.
#' @param enhance A vector containing the indeces of the DVHs curves to be enhanced in the plot by changing de default color.
#' @param mean.dvh A \code{logical} value, if \code{TRUE} the mean dvh is plotted.
#' @param median.dvh A \code{logical} value, if \code{TRUE} the median dvh is plotted.
#' @param mean.median.alone A \code{logical} value, if \code{TRUE} only the mean and/or median DVH(s) is/are plotted.
#'        Mean or median DVHs are plotted according the values in \code{mean.dvh} and \code{median.dvh}.
#' @param el.color The color for plotting DVHs. The default value is \code{"black"}.
#' @param en.color The color for plotting the enhanced DVHs. The default value is \code{"red"}.
#' @param mean.color The color for the mean dvh plot. The default value is \code{"blue"}.
#' @param median.color The color for the median dvh plot. The default value is \code{"red"}.
#' @param lwd The size of the lines in the plot. Default value is 1.
#' @param C.I.dvh Plot the confidence interval of mean and/or median DVH(s) according the values in \code{mean.dvh} and
#'        \code{median.dvh}.
#' @param C.I.dvh.width The width of confidence interval of mean, median DVH(s) and overall dvh series.
#' @param C.I.dvh.range Plot the confidence interval of the whole DVHs series according the value in \code{C.I.dvh.width}.
#' @param n.boot The number of bootstrap repetitions for calculating the confidence interval of mean and median DVHs.
#' @param ... Other parameters to be passed to \code{plot} function.
#' @details The dvh are stored in the slot \code{dvh} of a \code{dvhmatrix} class object. Of course it is possible for7
#'          users to use these matrices for plotting the dvhs considering that the structure of the matrix is referenced
#'          under \code{\link{dvhmatrix-class}}. In order to simplify the plotting operations this method has been coded
#'          allowing users to automatically obtain some useful features in the graphs, as the mean and median dvh curves,
#'          and the confidence intervals for mean, median and overall dvh matrix. The confidence intervals of mean and
#'          median are calculated by bootstrapping the whole \code{dvhmatrix} and calculating a given number of mean and
#'          median dvh values that finally are reduced by using a quantile function to create the final plot.
#' @exportMethod plot
#' @examples # create a dvhmatrix object
#' b <- DVH.generate(dvh.number = 100, dvh.type = "cumulative", vol.distr = "absolute")
#' plot(x = b, mean.dvh = TRUE, median.dvh = TRUE, mean.median.alone = TRUE, 
#'      C.I.dvh = TRUE, C.I.dvh.range = TRUE, C.I.dvh.width = .67, lwd = 2)
setMethod("plot", signature(x="dvhmatrix", y="missing"), 
          function(x, y, elements=NULL, enhance=NULL, mean.dvh=FALSE, median.dvh=FALSE, 
                   mean.median.alone=FALSE, el.color="black", en.color="red",  mean.color="blue", median.color="red",
                   lwd=1, C.I.dvh=FALSE, C.I.dvh.width=.95, C.I.dvh.range=FALSE, n.boot=2000, ...){
            dvh <- slot(x,"dvh")
            dvh.type <- slot(x, "dvh.type")
            vol.distr <- slot(x, "vol.distr")
            # define the lables for the y axis
            if ((dvh.type=="cumulative") && (vol.distr=="relative")) ylab<-"Volume [*100%]"
            if ((dvh.type=="cumulative") && (vol.distr=="absolute")) ylab<-"Volume [cc]"
            if ((dvh.type=="differential") && (vol.distr=="relative")) ylab<-"dVolume/dDose [*100%/Gy]"
            if ((dvh.type=="differential") && (vol.distr=="absolute")) ylab<-"dVolume/dDose [cc/Gy]"
            # max value of y axes
            ymax<-max(dvh[,2:ncol(dvh)])
            # direct plot for a single DVH
            if (ncol(dvh)==2) plot(x=dvh[,1], y=dvh[,2], type="l", ylab=ylab, xlab="Dose [Gy]", lwd=lwd)
            # plot for more than one DVH
            if ((ncol(dvh)>2) && (mean.median.alone==FALSE)) {
              # creates the vector of elements to be plotted and enhanced
              if (is.null(elements)) elements<-c(2:ncol(dvh)) else elements<-elements+1
              if (!is.null(enhance)) for (n in 1:length(enhance)) elements<-elements[elements!=(enhance[n]+1)]              
              # open the plot panel and raws the 1st plot
              if (length(elements)>0) plot(x=dvh[,1], y=dvh[,elements[1]], ylim=c(0, ymax), col=el.color, 
                                           lwd=lwd, type="l", xlab="Dose [Gy]", ylab=ylab)
              else
                if (length(elements)>0) lines(x=dvh[,1], y=dvh[,elements[1]], ylim=c(0, ymax), col=el.color, 
                                              lwd=lwd)              
              # plots the remaining DVHs into "elements"
              if (length(elements)==0) {
                start.enhance<-2
                plot(x=dvh[,1], y=dvh[,enhance[1]+1], col=en.color, lwd=lwd, type="l")
              } else start.enhance<-1
              if (length(elements)>1)
                for (n in 2:length(elements)) lines(x=dvh[,1], y=dvh[,elements[n]], col=el.color, lwd=lwd)
              # plots DVHs to be enhanced
              if (length(enhance)>0)
                for (n in start.enhance:length(enhance)) lines(x=dvh[,1], y=dvh[,enhance[n]+1], col=en.color, lwd=lwd)
            }
            # plot the mean DVH
            if (ncol(dvh)>2) {
              mean.DVH<-apply(X=dvh[,2:ncol(dvh)], MARGIN=1, FUN=mean)
              median.DVH<-apply(X=dvh[,2:ncol(dvh)], MARGIN=1, FUN=median)
            }
            if (mean.median.alone==FALSE) { 
              if ((mean.dvh==TRUE) && (ncol(dvh)>2)) 
                lines(x=dvh[,1], y=mean.DVH, col=mean.color, lwd=lwd)
              if ((median.dvh==TRUE) && (ncol(dvh)>2)) 
                lines(x=dvh[,1], y=median.DVH , col=median.color, lwd=lwd)
            }       
            
            if (mean.median.alone==TRUE) {
              # plot both mean and median
              if ((mean.dvh==TRUE) && (median.dvh==TRUE)) {
                plot(x=dvh[,1], y=mean.DVH, ylim=c(0, ymax), col=mean.color, 
                     lwd=lwd, type="l", xlab="Dose [Gy]", ylab=ylab)
                lines(x=dvh[,1], y=median.DVH , col=median.color, lwd=lwd)
              }
              # plot both only mean
              if ((mean.dvh==TRUE) && (median.dvh==FALSE)) {
                plot(x=dvh[,1], y=mean.DVH, ylim=c(0, ymax), col=mean.color, 
                     lwd=lwd, type="l", xlab="Dose [Gy]", ylab=ylab)
              }
              # plot only median
              if ((mean.dvh==FALSE) && (median.dvh==TRUE)) {
                plot(x=dvh[,1], y=median.DVH, ylim=c(0, ymax), col=median.color, 
                     lwd=lwd, type="l", xlab="Dose [Gy]", ylab=ylab)
              }
            }
            # plot of C.I. area, both 95% C.I. of the mean dose and 95% with quantiles are calculated
            if ((C.I.dvh==TRUE) && (ncol(dvh)>2)) {
              # load parallel computing packages if system is Linux
              if (Sys.info()["sysname"]=="Linux") {                     # check OS for using parallel computing under Linux
                require(foreach)
                require(doMC)
                require(parallel)
                registerDoMC(cores=detectCores())
              }
              # definition of function for sampling DVHs and calculating mean DVHs matrix
              sample.mean.dvh<-function(d) {
                # sample temp index of DVHs
                index<-sample(x=c(2:ncol(d)),size=ncol(d)-1,replace=T)
                # creates temp dvh
                temp.dvh<-cbind(d[,index])
                mean.v  <-apply(X=temp.dvh, MARGIN=1, FUN=mean)
                return (mean.v)
              }
              # definition of function for sampling DVHs and calculating median DVHs matrix
              sample.median.dvh<-function(d) {
                # sample temp index of DVHs
                index<-sample(x=c(2:ncol(d)),size=ncol(d)-1,replace=T)
                # creates temp dvh
                temp.dvh<-cbind(d[,index])
                median.v<-apply(X=temp.dvh, MARGIN=1, FUN=median)
                return (median.v)
              }
              # definition of function for sampling DVHs and calculating mean median DVHs matrix
              sample.mean.median.dvh<-function(d) {
                # sample temp index of DVHs
                index<-sample(x=c(2:ncol(d)),size=ncol(d)-1,replace=T)
                # creates temp dvh
                temp.dvh<-cbind(d[,index])
                median.v<-apply(X=temp.dvh, MARGIN=1, FUN=median)
                mean.v  <-apply(X=temp.dvh, MARGIN=1, FUN=mean)
                return(list(mean.dvh=mean.v, median.dvh=median.v))
              }
              if ((mean.dvh==TRUE)&&(median.dvh==FALSE)){
                # empty object for temp mean DVHs
                mean.v<-c()
                # create matrix of mean of resampled DVHs
                if (Sys.info()["sysname"]=="Linux") {
                  mean.v<-foreach(n=1:n.boot) %dopar% sample.mean.dvh(dvh)
                  mean.v<-unlist(mean.v)
                  mean.v<-matrix(data=mean.v, nrow=nrow(dvh))
                } else for (n in 1:n.boot)
                  mean.v  <-cbind(mean.v, sample.mean.dvh(dvh))             
                lowerCI.mean <- apply(X=mean.v, MARGIN=1, FUN=quantile, probs=(1-C.I.dvh.width)/2)
                higherCI.mean<- apply(X=mean.v, MARGIN=1, FUN=quantile, probs=(1+C.I.dvh.width)/2)
                xpoly.mean<-c(dvh[,1], rev(dvh[,1]))             # vector of x coordinates of polygon of CI
                ypoly.mean<-c(lowerCI.mean, rev(higherCI.mean))  # vector of y coordinates of polygon of CI of mean
                # draws the polygon of mean CI
                polygon(x=xpoly.mean, y=ypoly.mean, col=adjustcolor(col=mean.color, alpha.f=.25))
              }
              
              if ((mean.dvh==FALSE)&&(median.dvh==TRUE)){
                # empty object for temp median DVHs
                median.v<-c()                # create matrix of mean of resampled DVHs
                if (Sys.info()["sysname"]=="Linux") {
                  median.v<-foreach(n=1:n.boot) %dopar% sample.median.dvh(dvh)
                  median.v<-unlist(median.v)
                  median.v<-matrix(data=median.v, nrow=nrow(dvh))
                } else for (n in 1:n.boot)
                  median.v<-cbind(median.v, sample.median.dvh(dvh))                
                lowerCI.median <- apply(X=median.v, MARGIN=1, FUN=quantile, probs=(1-C.I.dvh.width)/2)
                higherCI.median<- apply(X=median.v, MARGIN=1, FUN=quantile, probs=(1+C.I.dvh.width)/2)
                xpoly.median<-c(dvh[,1], rev(dvh[,1]))                # vector of x coordinates of polygon of CI
                ypoly.median<-c(lowerCI.median, rev(higherCI.median)) # vector of y coordinates of polygon of CI of median
                # draws the polygon of median CI
                polygon(x=xpoly.median, y=ypoly.median, col=adjustcolor(col=median.color, alpha.f=.25))
              }
              
              if ((median.dvh==TRUE)&&(mean.dvh==TRUE)){
                # empty object for temp median DVHs
                median.v<-c()
                mean.v<-c()                
                # create matrix of mean of resampled DVHs
                if (Sys.info()["sysname"]=="Linux") {
                  synth.dvh<-foreach(n=1:n.boot) %dopar% sample.mean.median.dvh(dvh)
                  # big vector with all the DVHS for mean and median
                  synth.dvh<-unlist(synth.dvh)
                  nrow.dvh<-nrow(dvh)-1 # number of rows of dvh decreased by 1
                  build.index<-function(x) c(x:(x+nrow.dvh)) # function creating indeces numbers
                  start.mean<-which(names(synth.dvh)=="mean.dvh1")     # indeces of start mean dvh values
                  start.median<-which(names(synth.dvh)=="median.dvh1") # indeces of start median dvh values
                  # build mean dvh matrix
                  mean.v<-matrix(data=synth.dvh[sapply(X=start.mean, FUN=build.index, simplify=TRUE)], nrow=nrow(dvh))
                  # build median dvh matrix
                  median.v<-matrix(data=synth.dvh[sapply(X=start.median, FUN=build.index, simplify=TRUE)], nrow=nrow(dvh))
                  
                } else for (n in 1:n.boot){
                  median.v  <-cbind(median.v, sample.mean.median.dvh(dvh)$median.dvh)
                  mean.v  <-cbind(mean.v, sample.mean.median.dvh(dvh)$mean.dvh)
                }   
                lowerCI.mean <- apply(X=mean.v, MARGIN=1, FUN=quantile, probs=(1-C.I.dvh.width)/2)
                higherCI.mean<- apply(X=mean.v, MARGIN=1, FUN=quantile, probs=(1+C.I.dvh.width)/2)
                xpoly.mean<-c(dvh[,1], rev(dvh[,1]))             # vector of x coordinates of polygon of CI
                ypoly.mean<-c(lowerCI.mean, rev(higherCI.mean))  # vector of y coordinates of polygon of CI of mean
                lowerCI.median <- apply(X=median.v, MARGIN=1, FUN=quantile, probs=(1-C.I.dvh.width)/2)
                higherCI.median<- apply(X=median.v, MARGIN=1, FUN=quantile, probs=(1+C.I.dvh.width)/2)
                xpoly.median<-c(dvh[,1], rev(dvh[,1]))                # vector of x coordinates of polygon of CI
                ypoly.median<-c(lowerCI.median, rev(higherCI.median)) # vector of y coordinates of polygon of CI of median
                # draws the polygon of mean CI
                polygon(x=xpoly.mean, y=ypoly.mean, col=adjustcolor(col=mean.color, alpha.f=.25))
                # draws the polygon of median CI
                polygon(x=xpoly.median, y=ypoly.median, col=adjustcolor(col=median.color, alpha.f=.25))
              }
            }
            
            # plot the C.I. of data range of DVHs
            if (C.I.dvh.range==TRUE) {
              # quantiles for DVHs data range
              lowerCI <- c()
              higherCI <- c()
              for (n in 1:nrow(dvh)) {
                lowerCI <- c(lowerCI,  quantile(x=dvh[n,2:ncol(dvh)], probs=((1-C.I.dvh.width)/2)))
                higherCI<- c(higherCI, quantile(x=dvh[n,2:ncol(dvh)], probs=((1+C.I.dvh.width)/2)))
              }
              xpoly<-c(dvh[,1], rev(dvh[,1]))      # vector of x coordinates of polygon of CI
              ypoly<-c(lowerCI, rev(higherCI))     # vector of y coordinates of polygon of CI              
              polygon(x=xpoly, y=ypoly, col=adjustcolor(col="grey", alpha.f=.25))  # draws the polygon
            }            
          }
)

#' Function for extracting the V-Dose from cumulative DVH(s)
#' @description Function for calculating the value of V-Dose of a \code{dvhmatrix} object
#' @param dvh A \code{dvhmatrix} class object
#' @param Dose The dose for calculating the V-Dose value(s)
#' @details V-Dose is, in a given Dose Volume Histogram, the value of the volume of the structure receiving
#'          \strong{at least} the chosen level of dose.
#' @return A vector containing the value(s) of V-Dose(s).
#' @export
#' @examples # create a dvhmatrix class object
#' a<-DVH.generate(dvh.number = 100)
#' DVH.Vdose(dvh = a, Dose = 50)
DVH.Vdose <- function(dvh, Dose) {
  dvh<-DVH.diff.to.cum(dvh = dvh, relative = dvh@vol.distr)
  dvh.matrix<-dvh@dvh
  # try to find exact correspondence between a value in the dose column and the input Dose
  index<-which(x=dvh.matrix[,1]==Dose)
  # check the Dose
  if (Dose >= max(dvh.matrix[,1])) stop("Dose can not be >= maximum dose in dvh matrix")
  if (Dose <= 0) stop("Dose can not be <= 0")
  if (length(x=index)>0) return(as.numeric(dvh.matrix[index,2:ncol(dvh.matrix)])) else {
    # find the index of the Dose value in DVH closest to the input Dose
    index<-which.min(abs(dvh.matrix[,1]-Dose))        
    # store the number of columns in dvh.matrix    
    c<-ncol(dvh.matrix)
    if (dvh.matrix[index,1] < Dose) {
      # calculate the dose bin  
      Dbin <- dvh.matrix[index + 1,1] - dvh.matrix[index,1]
      # calculate the volume bin
      Vbin <- dvh.matrix[index,2:c] - dvh.matrix[index + 1,2:c]
      # calculate the step dose
      incD <- dvh.matrix[index + 1,1] - Dose
      # calucalte the step volume
      incV <- incD * Vbin / Dbin
      # calculation of Vdose
      return(as.numeric(incV + dvh.matrix[index + 1,2:c]))
    } else {
      # calculate the dose bin  
      Dbin <- dvh.matrix[index,1] - dvh.matrix[index - 1,1]
      # calculate the volume bin
      Vbin <- dvh.matrix[index - 1,2:c] - dvh.matrix[index,2:c]
      # calculate the step dose
      incD <- dvh.matrix[index,1] - Dose
      # calculation of Vdose
      incV <- Vbin / Dbin * incD
      return(as.numeric(incV + dvh.matrix[index,2:c]))
    }
  }
}

#' Function for extracting the D-Volume from cumulative DVH(s)
#' @description Function for calculating the value of D-Volume of a \code{dvhmatrix} object
#' @param dvh A \code{dvhmatrix} class object
#' @param The volume for calculating the D-Volume value(s)
#' @details D-Volume is, in a given Dose Volume Histogram, the value of the dose received
#'          by a given volume of the structure(s) of interest.
#' @return A vector containing the value(s) of D-Volume(s).
#' @export
#' @examples # create a dvhmatrix class object
#' a<-DVH.generate(dvh.number = 100)
#' DVH.Dvolume(dvh = a, Volume = 0.5)
DVH.Dvolume <- function(dvh,  Volume=0.001) {
  dvh.matrix<-DVH.diff.to.cum(dvh = dvh, relative = dvh@vol.distr)@dvh
  # create the vector wit indeces of Doses corresponding to
  # minimum difference between Volume threshold and given volume
  Dv <- c()
  Dbin <- dvh.matrix[2,1] - dvh.matrix[1,1] # dose bin
  for (n in 2:ncol(dvh.matrix)) {
    DvIndex <- which.min(abs(dvh.matrix[,n] - Volume))
    if (DvIndex==nrow(dvh.matrix)) {      # if Dvolume is lower than the minimum volume in DVH
      if (Volume < dvh.matrix[DvIndex, n]) {
        Dv <- c(Dv, NA)
        warning("One Dvolume is lower than minimum volume in DVH matrix")
        next
      }
    }
    DvSign <- sign(Volume - dvh.matrix[DvIndex, n])
    if (DvSign==0) { 
      Dv<-c(Dv, dvh.matrix[DvIndex, 1]) # the Dvolume is equal to a number in the matrix
      next
    }
    if (DvSign==-1) {   # Dvolume is greater than dose at DvIndex
      incV <- dvh.matrix[DvIndex, n] - Volume
      Vbin <- dvh.matrix[DvIndex, n] - dvh.matrix[DvIndex + 1, n]      
      incD <- Dbin *incV/Vbin
      Dv   <- c(Dv, incD + dvh.matrix[DvIndex, 1])
      next
    }    
    if (DvSign==1) {    # Dvolume is lower than dose at DvIndex
      incV <- Volume - dvh.matrix[DvIndex - 1, n]      
      Vbin <- dvh.matrix[DvIndex, n] - dvh.matrix[DvIndex - 1, n]      
      incD <- Dbin *incV/Vbin
      Dv   <- c(Dv, incD + dvh.matrix[DvIndex - 1, 1])
      next
    }
  }
  return(Dv)
}

#' Function for calculating the mean dvh with confidence interval
#' @description This function calculates the mean dvh from a \code{dvhmatrix} class object. The mean dvh is
#'              calculated with its confidence interval that is given by a bootstrapped dvh series from
#'              the dvh given in the \code{dvh} object.
#' @param dvh A \code{dvhmatrix} class object
#' @param C.I.width The width of confidence interval
#' @param n.boot The number of bootstrapped dvhs for computing the quantile in C.I. calculation
#' @return A \code{dvhmatrix} object where the column n. 2 is the mean dvh and columns 3 and 4 are low and high 
#'        boundaries of the confidence interval
#' @export
DVH.mean.dvh<-function(dvh, C.I.width = .95, n.boot = 2000) {
  # calculate mean dvh
  mean.dvh<-apply(X = dvh@dvh[,2:ncol(dvh@dvh)], MARGIN = 1, FUN = mean)
  Vdvh<-as.vector(dvh@dvh[,2:ncol(dvh@dvh)])
  # create the mean dvh vector
  meanV<-rep.int(x = 0, times = nrow(dvh@dvh) * n.boot)
  return(mean.dvh)
}