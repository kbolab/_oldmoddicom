#' function to calculate the first order features
#' 
#' @description  calculates Shannon entropy, kursosis and Skewness of a given list of arrays of voxels
#' @param inputData a list where each element is an array of the voxel of the image. Each element of the list normally refers to a patient.             
#' @return three lists: the first list contains the entropies, the second the kurtosis and the third the skewness
#' @export
#' @examples \dontrun{
#' # Create an instante of new.mmButo and load some cases
#' obj<-new.mmButo()
#' obj$loadCollection(Path = '/progetti/immagini/urinaEasy')
#' 
#' # get the three ROIs
#' Retto<-obj$getROIVoxel(ROIName="Retto")  
#' 
#' # get the possible biopsy
#' aa<-RAD.firstOrderFeatureImage(inputData = Retto )
#' aa$entropy
#' }#' #' 
#' @import entropy moments 
RAD.firstOrderFeatureImage <- function ( inputData )
{
  # set some variables;
  numPatient<-length(inputData)
  obj.mButo<-new.mmButo()
  ImageEntropy <- array(data = c(0), dim = c(numPatient))
  ImageKurtosis <- array(data = c(0), dim = c(numPatient))
  ImageSkewness <- array(data = c(0), dim = c(numPatient))
  histSamples<-500
  
  voxel.stats<-obj.mButo$getROIVoxelStats( inputData )
  maxVoxelValue<-max(voxel.stats$summary$max)
  minVoxelValue<-min(voxel.stats$summary$min)
  histSamples.array<-seq( from = minVoxelValue, to=maxVoxelValue, by = (maxVoxelValue-minVoxelValue)/histSamples   )
  
  # loop on each patient
  for (i in 1:numPatient)  {
    istogr <- c();    freq <- c()

    voxelCube.values<-unlist(inputData[[i]]$masked.images$voxelCube)
    voxelCube.values<-voxelCube.values[ voxelCube.values!=0  ] 
    # Calcola l'istogramma dei grigi
    istogr <- hist(voxelCube.values, breaks = histSamples,  plot=FALSE)
    # Calcola le frequenze da dover utilizzare nel calcolo dell'entropia
    freq <- freqs(y = istogr$counts)
    # Calcola l'entropia di Shannon per ogni paziente
    ImageEntropy[i] <- entropy.plugin(freqs = freq, unit = c("log2"))
    # Calcola la Kurtosis per ogni paziente (Kurtosis=0 distribuzione normale, Kurtosis > 0 distribuzione leptocurtica
    # cioè stretta, Kurtosis < 0 distribuzione platicurtica cioè larga)
    ImageKurtosis[i] <- kurtosis (x = voxelCube.values)
    # Calcola la Skewness per ogni paziente (Skewness = 0 simmetria perfetta, Skewness > 0 asimmetrica verso destra
    # Skewness < 0 asimmetrica verso sinistra)
    ImageSkewness[i] <- skewness(x = voxelCube.values)
  }
  #Restituisce una lista di array; ciascun valore corrisponde al singolo paziente ()
  return(list ("entropy"=ImageEntropy, "kurtosis"=ImageKurtosis, "skewness"=ImageSkewness) ) 
}
#' function to calculate Are/Volume and related measures
#' 
#' @description  calculates Are, Volume, Area/Volume Ratio and equivolumetric Spherical Area Ratio
#' @param listaROIVoxels an output of a \code{obj$getROIVoxel()} method
#' @return a list containing, for each patient the indicated measures
#' @export
#' @examples \dontrun{
#' # Create an instante of new.mmButo and load some cases
#' obj<-new.mmButo()
#' obj$loadCollection(Path = '/progetti/immagini/urinaEasy')
#' 
#' # get the three ROIs
#' Retto<-obj$getROIVoxel(ROIName="Retto")  
#' 
#' # get the possible biopsy
#' uu<-RAD.areaVolume(listaROIVoxels = Retto)
#' }#' #' 
RAD.areaVolume<-function( listaROIVoxels ) {
  objS<-services();
  obj.mmButo<-new.mmButo();
  arrayAV<-list()
  for ( i in names(listaROIVoxels) ) {
    geometry<-listaROIVoxels[[ i ]]$geometricalInformationOfImages;
    pSX<-geometry$pixelSpacing[1]
    pSY<-geometry$pixelSpacing[2]
    pSZ<-as.numeric(geometry$SliceThickness  )
    # expand the cropped voxelCube
    voxelCube <- obj.mmButo$mmButoLittleCube.expand(   listaROIVoxels[[i]] )
    arrayAV[[ i ]]$Area<-objS$SV.rawSurface(voxelMatrix = voxelCube, pSX = pSX, pSY=pSY,pSZ=pSZ)    
    if ( arrayAV[[ i ]]$Area == -1 ) {
      arrayAV[[ i ]]$Volume<- -1
      arrayAV[[ i ]]$equivolumetricSphericAreaRatio<- -1
    }
    else {
      arrayAV[[ i ]]$Volume<-length(which(voxelCube!=0))*pSX*pSY*pSZ
      arrayAV[[ i ]]$equivolumetricSphericAreaRatio<- ( 4*pi* (   (3/(4*pi))*arrayAV[[ i ]]$Volume   )^(2/3) ) / arrayAV[[ i ]]$Area
    }    
  }
  return(arrayAV)
}

#' Extract the possible Biopsies for a given ROIVoxelData
#' 
#' @description  Extract the possible Biopsies for a given ROIVoxelData. It can extract all the possible biopsies considering different volumes within a specified volume range along the axes. Please consider that if you give order to extract voxelCubes from dx.max=2,dy.max=2,dy.max=0 it will consider as biggest voxelCubes cube of 2*2+1,2*2+1,1 because with max.dx,max,dy,max,dz you specify only the distance from the centroid (the centroids itself is the '+1' in each direction)
#' @param ROIVoxelData as got from of a \code{obj$getROIVoxel()} method. It considers the cropped version and provide internally to explode it
#' @param dx.min the minumim number of voxels to be considered along the x dimension: default is equal to 2.
#' @param dy.min the minumim number of voxels to be considered along the y dimension: default is equal to 2
#' @param dz.min the minumim number of voxels to be considered along the z dimension: default is equal to 0' 
#' @param dx.max the maximum number of voxel to be considered for biopsy along x axes;
#' @param dy.max the maximum number of voxel to be considered for biopsy along y axes;
#' @param dz.max the maximum number of voxel to be considered for biopsy along z axes; 
#' @return a list containing a lot of things
#' @examples \dontrun{
#' # Create an instante of new.mmButo and load some cases
#' obj<-new.mmButo()
#' obj$loadCollection(Path = '/progetti/immagini/urinaEasy')
#' 
#' # get the three ROIs
#' Retto<-obj$getROIVoxel(ROIName="Retto")  
#' 
#' # get the possible biopsy
#' biopsy<-RAD.VirtualBiopsy(dx.max = 3,dy.max = 3,dz.max = 2,ROIVoxelData = Retto,dx.min = 2,dy.min = 2, dz.min = 1)
#' }#' 
#' @export
#' @useDynLib moddicom
RAD.VirtualBiopsy <- function ( ROIVoxelData, dx.min=2, dy.min=2, dz.min=0, dx.max, dy.max, dz.max)  {

  # instance the object justo to use the static methods
  # (i.e. to explode the voxel cubes)
  obj.mmButo<-new.mmButo()
  # set and initialize general variables
  h <- 1;  lista <- list();  pazienti <- list();
  NumPatients<-length( ROIVoxelData )
  # loop for each patient
  for (i in seq(1,NumPatients) )  {
    print(  paste(c("Processing: ",names(ROIVoxelData)[i]),collapse='' )  )
    # set some variables
    exam <- c();    carot <- c();
    # estraggo i dati del paziente i-esimo con le relative dimensioni. I need to expand the because
    # they are in the cropped format
    exam <- obj.mmButo$mmButoLittleCube.expand(   ROIVoxelData[[i]] )
    # get the dimension of the big voxel cube and set examArray
    dim1 <- dim(exam)[1];    dim2 <- dim(exam)[2];    dim3 <- dim(exam)[3]
    examArray <- array(data = exam, dim = (dim1*dim2*dim3))
    # build the expand.grid for each <dx,dy,dz> possible values 
    stepCheck <- expand.grid(
      indiceX = seq (dx.min, dx.max), 
      indiceY = seq (dy.min, dy.max), 
      indiceZ = seq (dz.min, dz.max))
    
    NumCheck <- nrow(stepCheck)
    stepCheckX <- array(data = stepCheck[,1], dim = c(NumCheck))
    stepCheckY <- array(data = stepCheck[,2], dim = c(NumCheck))
    stepCheckZ <- array(data = stepCheck[,3], dim = c(NumCheck))
    
    # setting some variables;
    lungh<-c();    expand<-c();    p<-1;    controlli<-c();    r<-c();  iteraz<-0;
    # loop for each possible configuration of the interested <dx,dy,dz>
    for ( s in seq(1,NumCheck) )   {
      if (stepCheckX[s]!=0 && stepCheckY[s]!=0) {
        iteraz <- iteraz+1
        # calculates how many controls are needed for each expand.grid (? ask Carlotta)
        expand <- expand.grid (indiceX=seq(-stepCheck[s,1],stepCheck[s,1]), 
                               indiceY=seq(-stepCheck[s,2],stepCheck[s,2]),
                               indiceZ=seq(-stepCheck[s,3],stepCheck[s,3]))        
        lungh[p] <- nrow(expand)
        controlli[p] <- s
        p <- p+1
      }
    }
    
    # beginning from the first <dx,dy,dz>
    for ( p in seq(1,iteraz) )    {
      # set some variables
      stepX <- stepCheckX[(controlli[p])];  stepY <- stepCheckY[(controlli[p])];  stepZ <- stepCheckZ[(controlli[p])]
      lunghezza <- lungh[p]; uni<-0;
      esameCarot <- array(data = c(0), dim = dim1*dim2*dim3)
      # call the .C procedure
      prova.carot <- .C( "new_virtualBiopsy", as.double (examArray), as.integer (dim1), as.integer (dim2), 
                         as.integer (dim3), as.integer(stepX), as.integer(stepY), 
                         as.integer(stepZ), as.integer (lunghezza), as.integer(esameCarot) ,
                         as.integer(uni) )
      carot <- array(data = prova.carot[[9]], dim = c(dim1, dim2, dim3))
      
      # calculates the coords <x,y,z> of the centroids
      carot.index <- which(carot==1, arr.ind = T )
      # save the desider output (hte <dx,dy,dz of the single analysis, the number of possible 
      # samples and the coords of the centroids)
      realX<-as.numeric(ROIVoxelData[[i]]$geometricalInformationOfImages$pixelSpacing[1])*stepX
      realX<-realX*2+realX
      realY<-as.numeric(ROIVoxelData[[i]]$geometricalInformationOfImages$pixelSpacing[2])*stepY
      realY<-realY*2+realY
      realZ<-as.numeric(ROIVoxelData[[i]]$geometricalInformationOfImages$SliceThickness)*stepY
      realZ<-realZ*2+realZ
      listLabel<-paste(c(stepX,"_",stepY,"_",stepZ),collapse='')
      lista[[ listLabel ]] <- list(
        "dx_dy_dz"=c(stepX*2+1,stepY*2+1, stepY*2+1), 
        "volume"=(stepX*2+1) * (stepY*2+1) *(stepY*2+1),
        "real_dx_dy_dz"=c(realX,realY,realZ), 
        "real_volume"=realX * realY *realZ,
        "NumCarotaggi"=prova.carot[[10]], 
        "IndexBiopsy"=carot.index)
      
    }
    pazienti[[ names(ROIVoxelData)[i] ]] <- lista
  }
  # define the returning object as member of the class 'virtualBiopsy'
  class(pazienti)<-'virtualBiopsy'
  return(pazienti)
}

