################################################################################
## script containing the function for generating and converting simulated DVH ##
## calculation of gEUD, mean dose, Vdose and Dvolume over DVH series          ##
################################################################################

# function for creating convex DVHs
convex.dvh <- function(convex.mind = 0, convex.maxd = 75, convex.dsbin = 0.5)
{
  nbins <- convex.maxd/convex.dsbin                                  # calculation of bin numbers  
  convex.maxd <- convex.maxd - 5                                     # set max dose  
  max_min <- convex.maxd - convex.mind                               # set new max
  newmin <- convex.mind + 0.8 * runif(1,0,max_min)                   # set new minimum dose
  meandose <- runif(1, newmin, convex.maxd)                          # set mean dose
  sddose <- (convex.maxd - meandose) * runif(1) + 1                  # set standard deviation
  temp.dvh <- dnorm(c(1:nbins)*convex.dsbin,meandose,sddose)         # output DVH
}

# convex.dvh <-function(convex.mind = 0, convex.maxd = 75, volume = 100, convex.dsbin = 0.5, volumeBin = .1) {
#   convex.maxd <- convex.maxd - 5                                     # set max dose  
#   max_min <- convex.maxd - convex.mind                               # set new max
#   newmin <- convex.mind + 0.8 * runif(1,0,max_min)                   # set new minimum dose
#   meandose <- runif(1, newmin, convex.maxd)                          # set mean dose
#   sddose <- (convex.maxd - meandose) * runif(1) + 1                  # set standard deviation
#   n <- 0
#   temp.d<-c()
#   while (n < volume) {
#     n <- n + volumeBin
#     temp.d <- c(temp.d, rnorm(n = 1, mean = meandose, sd = sddose))    
#   }
#   temp.dvh<-hist(x = temp.d, breaks = seq(from = 0, to = ceiling(max(temp.d)), by = convex.dsbin), plot = FALSE)$counts
# }

# function for creating concave DVHs
# concave.dvh <- function(concave.mind = 0, concave.maxd = 75, concave.dsbin = 0.5, sdadd = 2.5, meanpos=.1)
# {
#   nbins <- concave.maxd/concave.dsbin                                # calculation of bin numbers
#   concave.maxd <- (concave.maxd - concave.mind) * meanpos            # decreases the maxdose using mean position
#   meandose <- runif(1, concave.mind, concave.maxd)                   # mean dose calculation
#   sddose <- (concave.maxd - meandose) * runif(1) + sdadd             # standard deviation plus addendum to smooth DVH
#   temp.dvh <- dnorm(c(1:nbins)*concave.dsbin,meandose,sddose)        # temporary DVH
# }

# function for creating DVH with 2 mixed series
mix2.dvh <- function(mind = 0, maxd = 75, dsbin = 0.5)
{
  coeff1 <- runif(1,0,1)                                             # coefficient 1st temp DVH
  outputDVH <- convex.dvh(convex.mind=mind, convex.maxd=maxd, convex.dsbin=dsbin) * coeff1
  coeff2 <- runif(1,0,1)                                             # coefficient 2nd temp DVH
  outputDVH <- outputDVH + concave.dvh(concave.mind=mind, concave.maxd=maxd, concave.dsbin=dsbin) * coeff2                                  
}

# function for creating DVH with 3 mixed series
mix3.dvh <- function(mind = 0, maxd = 75, dsbin = 0.5)
{
  coeff1 <- runif(1,0,1)                                             # coefficient 1st temp DVH
  outputDVH <- convex.dvh(convex.mind=mind, convex.maxd=maxd, convex.dsbin=dsbin) * coeff1
  coeff2 <- runif(1,0,1)                                             # coefficient 2nd temp DVH
  outputDVH <- outputDVH + concave.dvh(concave.mind=mind, concave.maxd=maxd,
                                       concave.dsbin=dsbin,sdadd=1,meanpos=.2) * coeff2
  coeff3 <- runif(1,0,1)
  outputDVH <- outputDVH + concave.dvh(concave.mind=mind, concave.maxd=maxd,
                                       concave.dsbin=dsbin,sdadd=2.5,meanpos=.5) * coeff3
}


# function for creating different simulated differential DVHs
# type of DVH:
#   1: convex
#   2: concave
#   3: mix2
#   4: mix3
#   5: random
#   relative 'TRUE' option sets the outcome as relative DVHs

gen.dvh <- function(dvhnumber, type=c("random","convex","concave","mix2","mix3"), 
                    dvh.type=c("differential", "cumulative"), vol.distr=c("relative", "absolute"),
                    mindose = 0, maxdose = 75, dosebin = 0.5, mean.vol=200) {
  # mathces the type of DVH defined by the user
  dvh.type=match.arg(dvh.type)
  vol.distr=match.arg(vol.distr)
  # start with DVH generation
  nbins <- maxdose/dosebin                                          # number of DVH bins
  DVHList <- matrix(nrow=nbins, ncol=dvhnumber + 1)                 # create the matrix of DVHs
  doseseries <- seq(mindose, maxdose, dosebin)                      # series of bin of doses (x axes)
  type<-match.arg(type)                                             # default type is random
  if (type=="convex") {                                             # convex DVH generator
    for (m in 1:dvhnumber) {
      DVHList[,m + 1] <- convex.dvh(convex.mind=mindose, convex.maxd=maxdose,
                                    convex.dsbin=dosebin)           # append temp DVH in the matrix      
    }    
  }
  if (type=="concave") {                                            # concave DVH generator
    for (m in 1:dvhnumber) {
      DVHList[,m + 1] <- concave.dvh(concave.mind=mindose, 
                                     concave.maxd=maxdose, 
                                     concave.dsbin=dosebin)         # append temp DVH in the matrix
    }
  }
  if (type=="mix2") {                                               # mixed type DVH generator sum of 2 distributions
    for (m in 1:dvhnumber) {
      DVHList[,m + 1] <- mix2.dvh(mind=mindose,maxd=maxdose,dsbin=dosebin)
    }
  }  
  if (type=="mix3") {                                               # mixed type DVH generator sum of 3 distributions
    for (m in 1:dvhnumber) {
      DVHList[,m + 1] <- mix3.dvh(mind=mindose,maxd=maxdose,dsbin=dosebin)
    }
  }
  if (type=="random") {                                             # generate multiple type DVHs
    for (m in 1:dvhnumber) {
      s <- trunc(runif(1,1,5))                                      # seed for random DVH generation
      if (s==5) {s <- 4}                                            # exceptional result of randomization decreased to 4
      if (s==1) {                                                   # with seed=1 generate convex DVH 
        DVHList[,m + 1] <- convex.dvh(convex.mind=mindose, convex.maxd=maxdose,
                                      convex.dsbin=dosebin) 
      }
      if (s==2) {                                                   # with seed=2 generate concave DVH 
        DVHList[,m + 1] <- concave.dvh(concave.mind=mindose, 
                                       concave.maxd=maxdose, 
                                       concave.dsbin=dosebin, 
                                       sdadd=2.5, meanpos=.4)    
      }
      if (s==3) {                                                   # with seed=3 generate mix2 DVH
        DVHList[,m + 1] <- mix2.dvh(mind=mindose,maxd=maxdose,dsbin=dosebin)        
      }
      if (s==4) {                                                   # with seed=3 generate mix3 DVH
        DVHList[,m + 1] <- mix3.dvh(mind=mindose,maxd=maxdose,dsbin=dosebin)        
      }
    }
  }  
  DVHList[1,1]<-dosebin/2
  for (m in 2:nbins) DVHList[m, 1] <- DVHList[m-1,1] + dosebin      # add value of doses in 1st column of matrix
  
  # function for creating volumes according normal distribution
  multiply<-function(X) {abs(rnorm(n=1, mean=mean.vol, sd=mean.vol/10) * X)}  
  # generate return for the function
  if (dvh.type=="differential")
    if (vol.distr=="relative")  
      return(new("dvhmatrix", dvh=rel.diff.dvh(DVHList), vol.distr=vol.distr, 
                 dvh.type=dvh.type, volume=rnorm(n=(ncol(DVHList)-1), mean=mean.vol, sd=mean.vol/10)))
      else {
        DVHList<-rel.diff.dvh(DVHList)
        # generate the absolute value of the volumes in the DVH list
        abs.vol.dvh<-cbind(DVHList[,1], apply(X=DVHList[,2:(dvhnumber + 1)], MARGIN=2, FUN=multiply))
        return(new("dvhmatrix", dvh=abs.vol.dvh, vol.distr=vol.distr, 
                   dvh.type=dvh.type, volume=apply(X=abs.vol.dvh[,2:(dvhnumber + 1)], MARGIN=2, FUN=sum)))
      } 
    else
  if (dvh.type=="cumulative")
    if (vol.distr=="relative") 
      return(new("dvhmatrix", dvh=cum.dvh(dvh.matrix=DVHList, relative=TRUE), vol.distr=vol.distr, 
                 dvh.type=dvh.type, volume=rnorm(n=dvhnumber, mean=mean.vol, sd=mean.vol/10)))
      else {
        DVHList<-rel.diff.dvh(DVHList)
        # generate the absolute value of the volumes in the DVH list
        abs.vol.dvh<-cum.dvh(cbind(DVHList[,1], apply(X=DVHList[,2:(dvhnumber + 1)], MARGIN=2, FUN=multiply)), relative=FALSE)
      return(new("dvhmatrix", dvh=abs.vol.dvh, vol.distr=vol.distr, 
                 dvh.type=dvh.type, volume=abs.vol.dvh[1,2:(dvhnumber + 1)]))
      }
}

# function for converting differential DVHs into
# relative differential DVHS
rel.diff.dvh <- function(dvh.matrix) {
  dvh.size <- dim(dvh.matrix)
  DVHList <- matrix(nrow=dvh.size[1], ncol=dvh.size[2])       # create the matrix of  DVHs
  for (m in 2:dvh.size[2]) {                                  # loop for columns (volumes) 
    total.volume <- sum(dvh.matrix[,m])                       # calculate the total volume
    for (n in 1:dvh.size[1]) {                                # loop for rows
      DVHList[n, m] <- dvh.matrix[n, m]/total.volume          # elements of the matrix as relative volume
    }    
  }
  DVHList[,1]<-dvh.matrix[,1]                                 # sets the dose column  
  return(DVHList)
}

# function for converting cum DVHs into
# relative cum DVHS
rel.cum.dvh <- function(dvh.matrix) {
  dvh.size <- dim(dvh.matrix)
  DVHList <- matrix(nrow=dvh.size[1], ncol=dvh.size[2])       # create the matrix of  DVHs
  for (m in 2:dvh.size[2]) {                                  # loop for columns (volumes) 
    total.volume <- dvh.matrix[1,m]                           # calculate the total volume
    for (n in 1:dvh.size[1]) {                                # loop for rows
      DVHList[n, m] <- dvh.matrix[n, m]/total.volume          # elements of the matrix as relative volume
    }    
  }
  DVHList[,1]<-dvh.matrix[,1]                                 # sets the dose column  
  return(DVHList)
}

# function for converting a matrix of differential DVH 
# into matrix of cumulative DVH (relative volume)  alias=DVH.diff.to.cum
cum.dvh <- function(dvh.matrix, relative=TRUE)
{
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
  return(DVHList)
}

# function for converting a matrix of cumulative DVH into matrix of differential DVH
# alias DVH.cum.to.diff
diff.dvh <- function(dvh.matrix, relative=TRUE) {
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
  return(DVHList)
}

# function for extracting the Vdose from CUMULATIVE DVHs
Vdose <- function(dvh.matrix, Dose) {
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

# function for extracting the Dvolume from a CUMULATIVE DVH
# with maximum dose delivered to a given fraction (Volume) of volume
# default value 0.001 is for detecting MAXIMUM dose
Dvolume <- function(dvh.matrix,  Volume=0.001) {
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

# optimized function for EUD calculation
# works ONLY with relative differential DVHs
EUD <- function(dvh.matrix, a) {  
  dvh.size<-dim(dvh.matrix)                                         # size of DVH
  doseV<-dvh.matrix[,1]^a                                           # vector of doses
  if (dvh.size[2]==2) {
    addenda.EUD<-dvh.matrix[,2]*doseV
    eud<-(sum(addenda.EUD))^(1/a)
  } else {
    addenda.EUD<-dvh.matrix[,2:dvh.size[2]]*doseV                   # matrix of addenda EUD
    eud<-(apply(X=addenda.EUD,MARGIN=2,FUN=sum))^(1/a)
  }
  eud 
}

# optimized function for Mean dose calculation
MeanDose<-function(dvh.matrix) {
  dvh.size<-dim(dvh.matrix)
  doseV<-dvh.matrix[,1]
  if (dvh.size[2]==2) {
    addenda.mean<-dvh.matrix[,2]*doseV
    MeanD<-sum(addenda.mean)
  } else {
    addenda.mean<-dvh.matrix[,2:dvh.size[2]]*doseV
    MeanD<-apply(X=addenda.mean,MARGIN=2,FUN=sum)
  }
  MeanD
}

## function for correcting dose vector in matrices by Linear Quadratic model
correctLQ<-function(doses,FractionsNum,RefFractDose=2,alfabeta=3) {
  # calculate the vector of LQ corrected doses
  BNDoses<-((doses/FractionsNum)+alfabeta)/(RefFractDose+alfabeta)*doses
  return(BNDoses)
}

## function for extracting DVH from they array of dose points achieved by DICOM object
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


# function for joining two dvhmatrix objects
join.dvh<-function(receiver=NULL, addendum=NULL) {
  # check dvhmatrix class
  if ((class(receiver)!="dvhmatrix") || (class(addendum)!="dvhmatrix")) 
    stop("BOTH dvh.to.add AND destination MUST be dvhmatrix class objects")
  # creates list of dose bins
  dbin<-list(receiver=receiver@dvh[3,1]-receiver@dvh[2,1], addendum=addendum@dvh[3,1]-addendum@dvh[2,1])
  # detect the volume distribution  
  if (receiver@vol.distr=="relative") rel<-TRUE else rel<-FALSE  # relative state of receiver
 
  # STEP (1): identify conversion function if DVH types are different and convert addendum according receiver type
  if (receiver@dvh.type!=addendum@dvh.type) {
    if (receiver@dvh.type=="cumulative") conv.funct.name<-"cum.dvh"
    if (receiver@dvh.type=="differential") conv.funct.name<-"diff.dvh"
    conv.funct<-match.fun(FUN=conv.funct.name)  
    # correct the addendum if receiver has absolute distribution
    if ((rel==FALSE) && (addendum@vol.distr=="relative")) 
      for (n in 2:ncol(addendum@dvh)) addendum@dvh[,n]<-addendum@dvh[,n]*addendum@volume[n-1]
    addendum@dvh<-conv.funct(dvh.matrix=addendum@dvh, relative=rel)    
  }
  
  # STEP (2): check the dose-bin for interpolating addendum or receiver according the lowest dbin
  if (dbin$receiver!=dbin$addendum) {
    warning("Different dose bins between receiver and addendum: linear interpolation performed.")
    # create temp cumulative dvh, needed for avoiding interpolation artifacts
    if (receiver@dvh.type=="differential") {
      temp.receiver<-cum.dvh(dvh.matrix=receiver@dvh, relative=rel)
      temp.addendum<-cum.dvh(dvh.matrix=addendum@dvh, relative=rel)
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
        temp.addendum<-diff.dvh(dvh.matrix=temp.out, relative=rel)
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
        temp.receiver<-diff.dvh(dvh.matrix=temp.out, relative=rel)
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
# setMethod("plot", signature(x="dvhmatrix", y="missing"),
#           function(x, y, elements=NULL, enhance=NULL, mean.dvh=FALSE, median.dvh=FALSE, 
#                    mean.median.alone=FALSE, el.color="black", en.color="red",  mean.color="blue", median.color="red",
#                    lwd=1, C.I.dvh=FALSE, C.I.dvh.width=.95, min.sample.num=1000,...){
#             dvh <- slot(x,"dvh")
#             dvh.type <- slot(x, "dvh.type")
#             vol.distr <- slot(x, "vol.distr")
#             # define the lables for the y axis
#             if ((dvh.type=="cumulative") && (vol.distr=="relative")) ylab<-"Volume [*100%]"
#             if ((dvh.type=="cumulative") && (vol.distr=="absolute")) ylab<-"Volume [cc]"
#             if ((dvh.type=="differential") && (vol.distr=="relative")) ylab<-"dVolume/dDose [*100%/Gy]"
#             if ((dvh.type=="differential") && (vol.distr=="absolute")) ylab<-"dVolume/dDose [cc/Gy]"
#             # max value of y axes
#             ymax<-max(dvh[,2:ncol(dvh)])
#             # direct plot for a single DVH
#             if (ncol(dvh)==2) plot(x=dvh[,1], y=dvh[,2], type="l", ylab=ylab, xlab="Dose [Gy]", lwd=lwd)
#             # plot for more than one DVH
#             if ((ncol(dvh)>2) && (mean.median.alone==FALSE)) {
#               # creates the vector of elements to be plotted and enhanced
#               if (is.null(elements)) elements<-c(2:ncol(dvh)) else elements<-elements+1
#               if (!is.null(enhance)) for (n in 1:length(enhance)) elements<-elements[elements!=(enhance[n]+1)]              
#               # open the plot panel and raws the 1st plot
#               if (length(elements)>0) plot(x=dvh[,1], y=dvh[,elements[1]], ylim=c(0, ymax), col=el.color, 
#                                              lwd=lwd, type="l", xlab="Dose [Gy]", ylab=ylab)
#                 else
#               if (length(elements)>0) lines(x=dvh[,1], y=dvh[,elements[1]], ylim=c(0, ymax), col=el.color, 
#                                              lwd=lwd)              
#               # plots the remaining DVHs into "elements"
#               if (length(elements)==0) {
#                 start.enhance<-2
#                 plot(x=dvh[,1], y=dvh[,enhance[1]+1], col=en.color, lwd=lwd, type="l")
#               } else start.enhance<-1
#               if (length(elements)>1)
#                 for (n in 2:length(elements)) lines(x=dvh[,1], y=dvh[,elements[n]], col=el.color, lwd=lwd)
#               # plots DVHs to be enhanced
#               if (length(enhance)>0)
#                 for (n in start.enhance:length(enhance)) lines(x=dvh[,1], y=dvh[,enhance[n]+1], col=en.color, lwd=lwd)
#             }
#             # plot the mean DVH
#             mean.DVH<-apply(X=dvh[,2:ncol(dvh)], MARGIN=1, FUN=mean)
#             median.DVH<-apply(X=dvh[,2:ncol(dvh)], MARGIN=1, FUN=median)
#             if (mean.median.alone==FALSE) { 
#               if ((mean.dvh==TRUE) && (ncol(dvh)>2)) 
#                 lines(x=dvh[,1], y=mean.DVH, col=mean.color, lwd=lwd)
#               if ((median.dvh==TRUE) && (ncol(dvh)>2)) 
#                 lines(x=dvh[,1], y=median.DVH , col=median.color, lwd=lwd)
#             }       
#             
#             if (mean.median.alone==TRUE) {
#               # plot both mean and median
#               if ((mean.dvh==TRUE) && (median.dvh==TRUE)) {
#                 plot(x=dvh[,1], y=mean.DVH, ylim=c(0, ymax), col=mean.color, 
#                      lwd=lwd, type="l", xlab="Dose [Gy]", ylab=ylab)
#                 lines(x=dvh[,1], y=median.DVH , col=median.color, lwd=lwd)
#               }
#               # plot both only mean
#               if ((mean.dvh==TRUE) && (median.dvh==FALSE)) {
#                 plot(x=dvh[,1], y=mean.DVH, ylim=c(0, ymax), col=mean.color, 
#                      lwd=lwd, type="l", xlab="Dose [Gy]", ylab=ylab)
#               }
#               # plot only median
#               if ((mean.dvh==FALSE) && (median.dvh==TRUE)) {
#                 plot(x=dvh[,1], y=median.DVH, ylim=c(0, ymax), col=median.color, 
#                      lwd=lwd, type="l", xlab="Dose [Gy]", ylab=ylab)
#               }
#             }
#             
#             # plot of C.I. area, both 95% C.I. of the mean dose and 95% with quantiles are calculated
#             if ((C.I.dvh==TRUE) && (ncol(dvh)>2)) {
#               # calculate the bootstrap C.I. of DVHs
#               ndvh<-ncol(dvh) # use bootstrap if the number of dvh is lower than min.sample.num
#               lowerCI.mean<-c()    # vector of lower bound of C.I. of the mean
#               higherCI.mean<-c()   # vector of higher bound of C.I. of the mean
#               lowerCI.median<-c()  # vector of lower bound of C.I. of the median
#               higherCI.median<-c() # vector of higher bound of C.I. of the median
#               
#               # STEP (1): C.I. for the mean 
#               for (n in 1:nrow(dvh)) {
#                 bstrap.mean<-c()     # bootstrap vector of the mean
#                 bstrap.median<-c()   # bootstrap vector of the median
#                 for (i in 1:min.sample.num) {
#                   bsample<-sample(dvh[n,2:ndvh], (ndvh-1), replace=T)
#                   bestimate.mean<-mean(bsample)
#                   bstrap.mean<-c(bstrap.mean, bestimate.mean)
#                   bestimate.median<-median(bsample)
#                   bstrap.median<-c(bstrap.median, bestimate.median)
#                 }
#                 lowerCI.mean <- c(lowerCI.mean,  quantile(x=bstrap.mean, probs=((1-C.I.dvh.width)/2)))
#                 higherCI.mean<- c(higherCI.mean, quantile(x=bstrap.mean, probs=((1+C.I.dvh.width)/2)))
#                 lowerCI.median <- c(lowerCI.median,  quantile(x=bstrap.median, probs=((1-C.I.dvh.width)/2)))
#                 higherCI.median<- c(higherCI.median, quantile(x=bstrap.median, probs=((1+C.I.dvh.width)/2)))
#               } 
#               # draws mean C.I.
#               lowerCI.mean<-smooth.spline(x=dvh[,1], y=lowerCI.mean)$y
#               higherCI.mean<-smooth.spline(x=dvh[,1], y=higherCI.mean)$y  # reverse the vector of higher confidence interval         
#               xpoly.mean<-c(dvh[,1], rev(dvh[,1])) # vector of x coordinates of polygon of CI
#               ypoly.mean<-c(lowerCI.mean, rev(higherCI.mean))     # vector of y coordinates of polygon of CI of mean
#               if (mean.dvh==TRUE)
#                 polygon(x=xpoly.mean, y=ypoly.mean, col=adjustcolor(col=mean.color, alpha.f=.25))  # draws the polygon
#               
#               
#               if (median.dvh==TRUE){
#                 if (dvh.type=="cumulative") {                
#                   require(zoo)
#                   # function for smoothing high bound median CI curves
#                   smooth.y.high<-function(y) {
#                     n<-1
#                     succ<-c()
#                     while (n<length(y)) {
#                       if (length(which(y[(n+1):length(y)]>y[n]))>0) # devi identificare i valori per la rampa bassa
#                         succ<-c(succ, NA) else succ<-c(succ, y[n])
#                       n<-n+1
#                     }  
#                     if (length(is.na(succ))>0)
#                       succ<-c(na.approx(object=succ), y[length(y)]) else succ<-y
#                     return(succ)
#                   }
#                   # function for smoothing low bound median CI curves
#                   smooth.y.low<-function(y) {
#                     n<-2
#                     succ<-c()
#                     while (n<=length(y)) {
#                       if (length(which(y[1:(n-1)]<y[n]))>0) # devi identificare i valori per la rampa bassa
#                         succ<-c(succ, NA) else succ<-c(succ, y[n])
#                       n<-n+1
#                     }  
#                     if (length(is.na(succ))>0)
#                       succ<-c(y[1], na.approx(object=succ)) else succ<-y
#                     return(succ)
#                   }
#                   # smooth the median curves
#                   higherCI.median<-smooth.y.high(y=higherCI.median)                
#                   lowerCI.median <- smooth.y.low(y=lowerCI.median)
#                 }
#                 
#                 # further spline smoothing
#                 higherCI.median<-smooth.spline(x=dvh[,1], y=higherCI.median)$y
#                 lowerCI.median <- smooth.spline(x=dvh[,1], y=lowerCI.median)$y
#                 xpoly.median<-c(dvh[,1], rev(dvh[,1])) # vector of x coordinates of polygon of CI
#                 ypoly.median<-c(lowerCI.median, rev(higherCI.median))
#                 # draws the polygon
#                 polygon(x=xpoly.median, y=ypoly.median, col=adjustcolor(col=median.color, alpha.f=.25))
#               } 
#               # STEP (2): quantiles for DVHs data range
#               lowerCI <- c()
#               higherCI <- c()
#               for (n in 1:nrow(dvh)) {
#                 lowerCI <- c(lowerCI,  quantile(x=dvh[n,2:ndvh], probs=((1-C.I.dvh.width)/2)))
#                 higherCI<- c(higherCI, quantile(x=dvh[n,2:ndvh], probs=((1+C.I.dvh.width)/2)))
#               }
#               xpoly<-c(dvh[,1], rev(dvh[,1]))      # vector of x coordinates of polygon of CI
#               ypoly<-c(lowerCI, rev(higherCI))     # vector of y coordinates of polygon of CI              
#               polygon(x=xpoly, y=ypoly, col=adjustcolor(col="grey", alpha.f=.25))  # draws the polygon
#             }            
#           }
#   )



# method for class dvhmatrix
# enhance: is a vector of number of DVHs to be enhanced in the plot
# mean: plots the mean DVH overlapped to the background
# elements: if equal to a string "all" (default value) plots all the DVHs, if it is a vector plots only the
#   chosen DVHs
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
              # resample matrix of DVHs
              index<-sample(x=c(2:(ncol(dvh))), size=(ncol(dvh)-1)*n.boot, replace=T)
              # array of resampled DVHs
              resample.array<-array(data=dvh[,index], dim=c(nrow(dvh), ncol(dvh)-1, n.boot))
              
              if (mean.dvh==TRUE)
                # create matrix of mean of resampled DVHs
                mean.v  <- apply(X=resample.array, MARGIN=c(1,3), FUN=mean)
                lowerCI.mean <- apply(X=mean.v, MARGIN=1, FUN=quantile, probs=(1-C.I.dvh.width)/2)
                higherCI.mean<- apply(X=mean.v, MARGIN=1, FUN=quantile, probs=(1+C.I.dvh.width)/2)
                xpoly.mean<-c(dvh[,1], rev(dvh[,1]))             # vector of x coordinates of polygon of CI
                ypoly.mean<-c(lowerCI.mean, rev(higherCI.mean))  # vector of y coordinates of polygon of CI of mean
                # draws the polygon of means
                polygon(x=xpoly.mean, y=ypoly.mean, col=adjustcolor(col=mean.color, alpha.f=.25))  # draws the polygon
              
              
              if (median.dvh==TRUE){
                # create matrix of median of resampled DVHs
                median.v<- apply(X=resample.array, MARGIN=c(1,3), FUN=median)
                lowerCI.median <- apply(X=median.v, MARGIN=1, FUN=quantile, probs=(1-C.I.dvh.width)/2)
                higherCI.median<- apply(X=median.v, MARGIN=1, FUN=quantile, probs=(1+C.I.dvh.width)/2)
                xpoly.median<-c(dvh[,1], rev(dvh[,1]))                # vector of x coordinates of polygon of CI
                ypoly.median<-c(lowerCI.median, rev(higherCI.median)) # vector of y coordinates of polygon of CI of median
                # draws the polygon
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



# method for class dvhmatrix
# enhance: is a vector of number of DVHs to be enhanced in the plot
# mean: plots the mean DVH overlapped to the background
# elements: if equal to a string "all" (default value) plots all the DVHs, if it is a vector plots only the
#   chosen DVHs
setMethod("plot", signature(x="dvhmatrix", y="missing"), 
          function(x, y, elements=NULL, enhance=NULL, mean.dvh=FALSE, median.dvh=FALSE, 
                   mean.median.alone=FALSE, el.color="black", en.color="red",  mean.color="blue", median.color="red",
                   lwd=1, C.I.dvh=FALSE, C.I.dvh.width=.95, C.I.dvh.range=FALSE, n.boot=2000,...){
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

