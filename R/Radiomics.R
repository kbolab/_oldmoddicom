#' function to calculate the first order features
#' 
#' @description  calculates Shannon entropy, kursosis and Skewness of a given list of arrays of voxels
#' @param inputData a list where each element is an array of the voxel of the image. Each element of the list normally refers to a patient.             
#' @param histSamples in order to avoid distorsion in the calculus of the entropy the histogram can be sampled to a wished number of bins. The defaul is 150, any other integer value can be set using this parameter.
#' @return three lists: the first list contains the entropies, the second the kurtosis and the third the skewness
#' @export
#' @import entropy moments 
RAD.firstOrderFeatureImage <- function ( inputData, histSamples = 150)
{
  numPatient<-length(inputData)
  ImageEntropy <- array(data = c(0), dim = c(numPatient))
  ImageKurtosis <- array(data = c(0), dim = c(numPatient))
  ImageSkewness <- array(data = c(0), dim = c(numPatient))
  for (i in 1:numPatient)
  {
    paziente <- c()
    istogr <- c()
    freq <- c()
    # Carica i dati dell'iesimo paziente
    paziente <- unlist(as.array(x = inputData[[i]]))
    # Calcola l'istogramma dei grigi
    istogr <- hist(paziente, breaks = histSamples, freq = FALSE, main = paste(c("Histogram Breaks=",histSamples),collapse='') )
    # Calcola le frequenze da dover utilizzare nel calcolo dell'entropia
    freq <- freqs(y = istogr$counts)
    # Calcola l'entropia di Shannon per ogni paziente
    ImageEntropy[i] <- entropy.plugin(freqs = freq, unit = c("log2"))
    # Calcola la Kurtosis per ogni paziente (Kurtosis=0 distribuzione normale, Kurtosis > 0 distribuzione leptocurtica
    # cioè stretta, Kurtosis < 0 distribuzione platicurtica cioè larga)
    ImageKurtosis[i] <- kurtosis (x = paziente)
    # Calcola la Skewness per ogni paziente (Skewness = 0 simmetria perfetta, Skewness > 0 asimmetrica verso destra
    # Skewness < 0 asimmetrica verso sinistra)
    ImageSkewness[i] <- skewness(x = paziente)
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
RAD.areaVolume<-function( listaROIVoxels ) {
  objS<-services()
  arrayAV<-list()
  for ( i in names(listaROIVoxels) ) {
    geometry<-listaROIVoxels[[ i ]]$geometricalInformationOfImages;
    pSX<-geometry$pixelSpacing[1]
    pSY<-geometry$pixelSpacing[2]
    pSZ<-as.numeric(geometry$SliceThickness  )
    voxelCube<-listaROIVoxels[[ i ]]$masked.images;
    arrayAV[[ i ]]$Area<-objS$SV.rawSurface(voxelMatrix = voxelCube, pSX = pSX, pSY=pSY,pSZ=pSZ)    
    if ( arrayAV[[ i ]]$Area == -1 ) {
      arrayAV[[ i ]]$Volume<- -1
      arrayAV[[ i ]]$equivolumetricSphericAreaRadio<- -1
    }
    else {
      arrayAV[[ i ]]$Volume<-length(which(voxelCube!=0))*pSX*pSY*pSZ
      arrayAV[[ i ]]$equivolumetricSphericAreaRadio<- ( 4*pi* (   (3/(4*pi))*arrayAV[[ i ]]$Volume   )^(2/3) ) / arrayAV[[ i ]]$Area
    }    
  }
  return(arrayAV)
}



