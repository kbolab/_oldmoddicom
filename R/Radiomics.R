#' function to calculate the first order features
#' 
#' @description  Calculates Shannon entropy, kursosis, Skewness, mean, standard deviation and energy of a given list of arrays of voxels
#' @param inputData a list where each element is an array of the voxel of the image. Each element of the list normally refers to a patient.             
#' @return six lists: the first list contains the entropies, the second the kurtosis, the third the skewness, the fourth the mean, the fifth the standard deviation and the sixth the energy
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
RAD.firstOrderFeatureImage <- function ( inputData ) { 
  # NON USATA, DA RANZARE?
  if(class(inputData) == "geoLetStructureVoxelList") 
    return(RAD.firstOrderFeatureImage.geoLet( inputData))
  if(class(inputData) == "mmButoStructureVoxelList")  
    return(RAD.firstOrderFeatureImage.mmButo( inputData))
}
RAD.firstOrderFeatureImage.geoLet<-function(inputData) {
  # NON USATA, DA RANZARE?
  numPatient<-1
  ImageEntropy <- array(data = c(0), dim = c(numPatient))
  ImageKurtosis <- array(data = c(0), dim = c(numPatient))
  ImageSkewness <- array(data = c(0), dim = c(numPatient))
  ImageMean <- array(data = c(0), dim = c(numPatient))
  ImageStandDeviat <- array(data = c(0), dim = c(numPatient))
  ImageEnergy <- array(data = c(0), dim = c(numPatient))
  
  histSamples<-500
  maxVoxelValue<-max(inputData$masked.images$voxelCube);  minVoxelValue<-min(inputData$masked.images$voxelCube);
  histSamples.array<-seq( from = minVoxelValue, to=maxVoxelValue, by = (maxVoxelValue-minVoxelValue)/histSamples   )

  istogr <- c();    freq <- c();  i<-1;
  if(is.list(inputData) == TRUE  ) {
    #voxelCube.values<-inputData$masked.images$voxelCube[ inputData$masked.images$voxelCube!=0  ]
    voxelCube.values<-inputData$masked.images$voxelCube[ !is.na(inputData$masked.images$voxelCube)  ] 
    
    # Calcola l'istogramma dei grigi
    #istogr <- hist(voxelCube.values, breaks = histSamples,  plot=FALSE)
    istogr <- discretize(x = voxelCube.values[which(!is.na(voxelCube.values),arr.ind = TRUE ) ], numBins = histSamples)        
    # Calcola le frequenze da dover utilizzare nel calcolo dell'entropia
    #freq <- freqs(y = istogr$counts)
    # Calcola l'entropia di Shannon per ogni paziente
    #ImageEntropy[i] <- entropy.plugin(freqs = freq, unit = c("log2"))
    ImageEntropy[i] <- entropy.plugin(freqs = istogr, unit = c("log2"))
    # Calcola la Kurtosis per ogni paziente (Kurtosis=0 distribuzione normale, Kurtosis > 0 distribuzione leptocurtica
    # cioè stretta, Kurtosis < 0 distribuzione platicurtica cioè larga)
    ImageKurtosis[i] <- kurtosis (x = voxelCube.values[which(!is.na(voxelCube.values),arr.ind = TRUE ) ])
    # Calcola la Skewness per ogni paziente (Skewness = 0 simmetria perfetta, Skewness > 0 asimmetrica verso destra
    # Skewness < 0 asimmetrica verso sinistra)
    ImageSkewness[i] <- skewness(x = voxelCube.values[which(!is.na(voxelCube.values),arr.ind = TRUE ) ])
    #Calcola media dei grigi dei singoli pazienti
    ImageMean[i] <- mean(x = voxelCube.values[which(!is.na(voxelCube.values),arr.ind = TRUE ) ])
    #Calcola deviazione standard
    ImageStandDeviat[i] <- sd (x = voxelCube.values[which(!is.na(voxelCube.values),arr.ind = TRUE ) ])
    
    #Calcola la Enery, permette di valutare l'uniformità dell'immagine. Può assumere un valore tra 0 e 1:
    # 0 = non uniforme; 1 = uniforme;
    Energy <- 0
    for (j in 1:length(istogr$density))
    {
      Energy <- Energy+(istogr$density[j])^2
    }
    ImageEnergy[i] <- Energy
  } else {
    ImageEntropy[i] <-NA;
    ImageKurtosis[i]<-NA;
    ImageSkewness[i]<-NA;
    ImageMean[i]<-NA;
    ImageStandDeviat[i]<-NA;
    ImageEnergy[i]<-NA;
  }  
  #Restituisce una lista di array; ciascun valore corrisponde al singolo paziente ()
  return(list ("entropy"=ImageEntropy, "kurtosis"=ImageKurtosis, "skewness"=ImageSkewness, "mean"=ImageMean, 
               "standardDeviation"=ImageStandDeviat, "energy"=ImageEnergy))   
}
RAD.firstOrderFeatureImage.mmButo <- function ( inputData ) {
  
  # set some variables;
  numPatient<-length(inputData)
  obj.mButo<-new.mmButo()
  ImageEntropy <- array(data = c(0), dim = c(numPatient))
  ImageKurtosis <- array(data = c(0), dim = c(numPatient))
  ImageSkewness <- array(data = c(0), dim = c(numPatient))
  ImageMean <- array(data = c(0), dim = c(numPatient))
  ImageStandDeviat <- array(data = c(0), dim = c(numPatient))
  ImageEnergy <- array(data = c(0), dim = c(numPatient))
  
  histSamples<-500
  
  voxel.stats<-obj.mButo$getROIVoxelStats( inputData )
  maxVoxelValue<-max(voxel.stats$summary$max)
  minVoxelValue<-min(voxel.stats$summary$min)
  histSamples.array<-seq( from = minVoxelValue, to=maxVoxelValue, by = (maxVoxelValue-minVoxelValue)/histSamples   )
  
  # loop on each patient
  for (i in 1:numPatient)  {
    istogr <- c();    freq <- c()
    if((TRUE %in% is.na(inputData[[i]])) == FALSE  ) {
        voxelCube.values<-unlist(inputData[[i]]$masked.images$voxelCube)
        voxelCube.values<-voxelCube.values[ voxelCube.values!=0  ] 
        # Calcola l'istogramma dei grigi
        #istogr <- hist(voxelCube.values, breaks = histSamples,  plot=FALSE)
        
        #istogr <- discretize(x = voxelCube.values, numBins = histSamples)
        istogr <- discretize(x = voxelCube.values[which(!is.na(voxelCube.values),arr.ind = TRUE ) ], numBins = histSamples)        
        # Calcola le frequenze da dover utilizzare nel calcolo dell'entropia
        #freq <- freqs(y = istogr$counts)
        # Calcola l'entropia di Shannon per ogni paziente
        #ImageEntropy[i] <- entropy.plugin(freqs = freq, unit = c("log2"))
        ImageEntropy[i] <- entropy.plugin(freqs = istogr, unit = c("log2"))
        # Calcola la Kurtosis per ogni paziente (Kurtosis=0 distribuzione normale, Kurtosis > 0 distribuzione leptocurtica
        # cioè stretta, Kurtosis < 0 distribuzione platicurtica cioè larga)
        ImageKurtosis[i] <- kurtosis (x = voxelCube.values[which(!is.na(voxelCube.values),arr.ind = TRUE ) ])
        # Calcola la Skewness per ogni paziente (Skewness = 0 simmetria perfetta, Skewness > 0 asimmetrica verso destra
        # Skewness < 0 asimmetrica verso sinistra)
        ImageSkewness[i] <- skewness(x = voxelCube.values[which(!is.na(voxelCube.values),arr.ind = TRUE ) ])
        #Calcola media dei grigi dei singoli pazienti
        ImageMean[i] <- mean(x = voxelCube.values[which(!is.na(voxelCube.values),arr.ind = TRUE ) ])
        #Calcola deviazione standard
        ImageStandDeviat[i] <- sd (x = voxelCube.values[which(!is.na(voxelCube.values),arr.ind = TRUE ) ])
        
        #Calcola la Enery, permette di valutare l'uniformità dell'immagine. Può assumere un valore tra 0 e 1:
        # 0 = non uniforme; 1 = uniforme;
        Energy <- 0
        for (j in 1:length(istogr))
        {
          Energy <- Energy+(freqs(istogr)[j])^2
        }
        
        ImageEnergy[i] <- Energy
    } else {
      ImageEntropy[i] <-NA;
      ImageKurtosis[i]<-NA;
      ImageSkewness[i]<-NA;
      ImageMean[i]<-NA;
      ImageStandDeviat[i]<-NA;
      ImageEnergy[i]<-NA;
    }
  }
  #Restituisce una lista di array; ciascun valore corrisponde al singolo paziente ()
  return(list ("entropy"=ImageEntropy, "kurtosis"=ImageKurtosis, "skewness"=ImageSkewness, "mean"=ImageMean, 
               "standardDeviation"=ImageStandDeviat, "energy"=ImageEnergy)) 
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
RAD.areaVolume <- function ( listaROIVoxels ) {
  if(class(listaROIVoxels) == "geoLetStructureVoxelList") 
    return(RAD.areaVolume.geoLet( listaROIVoxels ))
  if(class(listaROIVoxels) == "mmButoStructureVoxelList")  
    return(RAD.areaVolume.mmButo( listaROIVoxels ))
}
RAD.areaVolume.geoLet<-function( listaROIVoxels ) {
  objS<-services();  arrayAV<-list(); obj.mmButo<-new.mmButo(); i<-1;
  if(is.list(listaROIVoxels) == TRUE  ) {
    arrayAV[[ i ]]<-list();
    geometry<-listaROIVoxels$geometricalInformationOfImages;
    pSX<-geometry$pixelSpacing[1]
    pSY<-geometry$pixelSpacing[2]
    pSZ<-as.numeric(geometry$SliceThickness  )
    # expand the cropped voxelCube
    voxelCube <- listaROIVoxels$masked.images$voxelCube
    arrayAV[[ i ]]$Area<-objS$SV.rawSurface(voxelMatrix = voxelCube, pSX = pSX, pSY=pSY,pSZ=pSZ)    
    if ( arrayAV[[ i ]]$Area == -1 ) {
      arrayAV[[ i ]]$Volume<- -1
      arrayAV[[ i ]]$equivolumetricSphericAreaRatio<- -1
    }
    else {
      arrayAV[[ i ]]$Volume<-length(which(voxelCube!=0))*pSX*pSY*pSZ
      arrayAV[[ i ]]$equivolumetricSphericAreaRatio<- ( 4*pi* (   (3/(4*pi))*arrayAV[[ i ]]$Volume   )^(2/3) ) / arrayAV[[ i ]]$Area
    }
  } else {
    arrayAV[[ i ]]$Volume<-NA
    arrayAV[[ i ]]$equivolumetricSphericAreaRatio<-NA
    arrayAV[[ i ]]$Area<-NA
  }

  return(arrayAV[[i]])  
}
RAD.areaVolume.mmButo<-function( listaROIVoxels ) {
  objS<-services();
  obj.mmButo<-new.mmButo();
  arrayAV<-list()
  
  # progression bar  
  iterazione<-0
  pb = txtProgressBar(min = 0, max = length(listaROIVoxels), initial = 0)
  
  for ( i in names(listaROIVoxels) ) {
    if((TRUE %in% is.na(listaROIVoxels[[i]])) == FALSE  ) {
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
    } else {
      arrayAV[[ i ]]$Volume<-NA
      arrayAV[[ i ]]$equivolumetricSphericAreaRatio<-NA
      arrayAV[[ i ]]$Area<-NA
    }
    setTxtProgressBar(pb,iterazione)
    iterazione<-iterazione+1
    
  }
  close(pb)
  
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
#' @param sampleResultAt indicates the maximum number of biopsy for each case. Often having 500 or 10000000 biopsy is the same with the exception of the memory allocated: by this parameter the user can forse a \code{sample} in order to reduce the occupation in memory. Default is \code{Inf}
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
RAD.VirtualBiopsy <- function ( ROIVoxelData, dx.min=2, dy.min=2, dz.min=0, dx.max, dy.max, dz.max, sampleResultAt = Inf)  {
  # instance the object justo to use the static methods
  # (i.e. to explode the voxel cubes)
  obj.mmButo<-new.mmButo()
  # set and initialize general variables
  h <- 1;  lista <- list();  pazienti <- list();
  NumPatients<-length( ROIVoxelData )
  # loop for each patient
  for (i in seq(1,NumPatients) )  {
    print(  paste(c("Processing: ",names(ROIVoxelData)[i]),collapse='' )  )
    
    if((TRUE %in% is.na(ROIVoxelData[[i]])) == FALSE  ) {
      
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
        
        minValue<- min(examArray[which(!is.na(examArray),arr.ind = T)]) - 10000;
        examArray[which(is.na(examArray),arr.ind = T)]<- minValue
        
        # call the .C procedure
        prova.carot <- .C( "new_virtualBiopsy", as.double (examArray), as.integer (dim1), as.integer (dim2), 
                           as.integer (dim3), as.integer(stepX), as.integer(stepY), 
                           as.integer(stepZ), as.integer (lunghezza), as.integer(esameCarot) ,
                           as.integer(uni), as.double(minValue) )
        carot <- array(data = prova.carot[[9]], dim = c(dim1, dim2, dim3))
        
        carot[which(examArray == minValue,arr.ind = T)]<- NA
        examArray[which(examArray == minValue,arr.ind = T)]<- NA

        # calculates the coords <x,y,z> of the centroids
        carot.index <- which(carot==1, arr.ind = T )
        # save the desider output (hte <dx,dy,dz of the single analysis, the number of possible 
        # samples and the coords of the centroids)
        realX<-as.numeric(ROIVoxelData[[i]]$geometricalInformationOfImages$pixelSpacing[1])*stepX
        realX<-realX*2+realX
        realY<-as.numeric(ROIVoxelData[[i]]$geometricalInformationOfImages$pixelSpacing[2])*stepY
        realY<-realY*2+realY
        realZ<-as.numeric(ROIVoxelData[[i]]$geometricalInformationOfImages$SliceThickness)*stepZ
        realZ<-realZ*2+realZ
        listLabel<-paste(c(stepX,"_",stepY,"_",stepZ),collapse='')
        # fai un resample se non interessa avere indietro TUTTI i punti carotabili
        if( sampleResultAt != Inf  && dim(carot.index)[1]>sampleResultAt ) {
          listaPossibiliSamples<-seq(1,dim(carot.index)[1])
          listaPossibiliSamples<-sample(listaPossibiliSamples, size = sampleResultAt, replace = FALSE)
          carot.index.resampled<-carot.index[ listaPossibiliSamples,  ]
          prova.carot[[10]]<-listaPossibiliSamples
        }
        else carot.index.resampled<-carot.index
        
        lista[[ listLabel ]] <- list(
          "dx_dy_dz"=c(stepX*2+1,stepY*2+1, stepZ*2+1), 
          "volume"=(stepX*2+1) * (stepY*2+1) *(stepZ*2+1),
          "real_dx_dy_dz"=c(realX,realY,realZ), 
          "real_volume"=realX * realY *realZ,
          "NumCarotaggi"=prova.carot[[10]], 
          "IndexBiopsy"=carot.index.resampled)
        
      }
      pazienti[[ names(ROIVoxelData)[i] ]] <- lista
    }
    else  {
      pazienti[[ names(ROIVoxelData)[i] ]] <-NA;
    }
  }
  # define the returning object as member of the class 'virtualBiopsy'
  class(pazienti)<-'virtualBiopsy'
  return(pazienti)
}
#' Apply erosion to a set of voxelCubes
#' 
#' @description  Apply the erosion to a set of ROIs internal voxels
#' @param ROIVoxelData as got from of a \code{obj$getROIVoxel()} method. It considers the cropped version and provide internally to explode it
#' @param margin.x the erosion margin along the x axes: default is 2.
#' @param margin.y the erosion margin along the y axes: default is 2.
#' @param margin.z the erosion margin along the z axes: default is 1.
#' @return a list containing the eroded voxelCubes and some stats
#' @examples \dontrun{
#' # Create an instante of new.mmButo and load some cases
#' obj<-new.mmButo()
#' obj$loadCollection(Path = '/progetti/immagini/urinaEasy')
#' 
#' # get the three ROIs
#' Retto<-obj$getROIVoxel(ROIName="Retto")  
#' 
#' # get the possible biopsy
#' erodedCubes<-RAD.applyErosion( ROIVoxelData = Retto )
#' }#' 
#' @export
#' @useDynLib moddicom
RAD.applyErosion<-function(  ROIVoxelData, margin.x=2, margin.y=2, margin.z=1 ) {
  res<-list()
  print("Beginnig erosion....");
  
  for ( i in names(ROIVoxelData) ) {
    if((TRUE %in% is.na(ROIVoxelData[[i]])) == FALSE  ) {
      
      print( paste( c("Eroding: ",i)   ,collapse = '')    );
      # declare the lists
      res[[i]]<-list();    res[[i]]$stat<-list();

      # get the voxel cube and prepare the erosion
      erodedVoxelCube<-ROIVoxelData[[i]]$masked.images$voxelCube;
      minValue<- min(erodedVoxelCube[which(!is.na(erodedVoxelCube),arr.ind = T)]) - 10000;
      erodedVoxelCube[which(is.na(erodedVoxelCube),arr.ind = T)]<- minValue
      # get the dimensions and set the desired margins
      nX<-dim(erodedVoxelCube)[1];    nY<-dim(erodedVoxelCube)[2];    nZ<-dim(erodedVoxelCube)[3];
      mx<-margin.x; my<-margin.y; mz<-margin.z;
      iterator<-0; # this is just to avoid infinite loops...
      # erode it!
      
      aa<-.C("erosion",as.double(erodedVoxelCube), as.integer(nX), as.integer(nY), 
             as.integer(nZ),as.integer(margin.x),as.integer(margin.y), 
             as.integer(margin.z), as.integer(iterator), as.double(minValue)) 
      
      erodedVoxelCube<-array(aa[[1]], dim=c(nX,nY,nZ))
      erodedVoxelCube[which(erodedVoxelCube == minValue,arr.ind = T)]<- NA
      
      ROIVoxelData[[i]]$masked.images$voxelCube<-erodedVoxelCube;
    }
    else {
      print( paste( c("Eroding: ",i)   ,collapse = '')    );
      ROIVoxelData[[i]]<-NA
    }
  }
  return( ROIVoxelData )
}


#' Apply erosion to a set of voxelCubes
#' 
#' @description  Apply the erosion to a set of ROIs internal voxels
#' @export
RAD.compareSignals<-function( ROIVoxelData, bootStrappedTimes = 1000 ) {
  interpolatedDensity<-list();
  listaDensity<-list()
  valoreMassimoX<-0;valoreMassimoY<-0;
  addingMatrix<-c()
  bootstrappedAddingMatrix<-c();
  
  for(i in names( ROIVoxelData )) {
    if((TRUE %in% is.na(ROIVoxelData[[i]])) == FALSE  ) {
      valori<-ROIVoxelData[[i]]$masked.images$voxelCube[which(ROIVoxelData[[i]]$masked.images$voxelCube!=0)]
      listaDensity[[i]]<-density(valori)
      valoreMassimoX<-max(valoreMassimoX,listaDensity[[i]]$x)
      valoreMassimoY<-max(valoreMassimoY,listaDensity[[i]]$y)
    }
  }
  
  for(i in names( ROIVoxelData )) {
    if((TRUE %in% is.na(ROIVoxelData[[i]])) == FALSE  ) {
      interpolatedDensity[[i]]<-approx(listaDensity[[i]]$x,listaDensity[[i]]$y,n=valoreMassimoX, xout=seq( from=0 , to=max(valoreMassimoX) )) 
      # azzera gli NA
      interpolatedDensity[[i]]$y[which(is.na(interpolatedDensity[[i]]$y))]<-0   
      addingMatrix<-rbind(addingMatrix,interpolatedDensity[[i]]$y);
    } 
  }

  bootstrappedAddingMatrix<-apply(addingMatrix,2,sample,size=bootStrappedTimes,replace=T)
  quantileMatrix<-apply(bootstrappedAddingMatrix, 2, quantile, probs = c(.025, .975), na.rm = TRUE)   # versione bootstrappato
  
  x<-interpolatedDensity[[i]]$x
  x[which(is.na(x))]<-0
  
  quantileMatrixSPU<-smooth.spline( x = x, y = quantileMatrix[1,])
  quantileMatrixSPL<-smooth.spline( x = x, y = quantileMatrix[2,])
  
  meanMatrix<-colMeans(addingMatrix, na.rm = TRUE)      # media normale
  
  # plot it
  plot(x = 0, y=0, xlim=c(0,valoreMassimoX),ylim=c(0,valoreMassimoY))
  polygon(c(x,rev(x)),c(meanMatrix,rev(quantileMatrixSPU$y)), col=rgb(.7, .7, .7, 0.2), lty = c("dashed"))
  polygon(c(x,rev(x)),c(meanMatrix,rev(quantileMatrixSPL$y)), col=rgb(.7, .7, .7, 0.2), lty = c("dashed"))
  lines( x = x , y=meanMatrix )   

}

#' return a simple list of voxelcubes
#' 
#' @description  convertes the output of a \code{obj$getROIVoxel()} to an easy-to-handle list of voxelCubes
#' @param ROIVoxelData as got from of a \code{obj$getROIVoxel()} method. Cropped or expanded, it's the same
#' @return a list containing the eroded voxelCubes presented in a 'easy to handle' shape
#' @examples \dontrun{
#' # Create an instante of new.mmButo and load some cases
#' obj<-new.mmButo()
#' obj$loadCollection(Path = '/progetti/immagini/urinaEasy')
#' 
#' # get the three ROIs
#' GTV<-obj$getROIVoxel(ROIName="GTV")  
#' 
#' # get the possible biopsy
#' res<-RAD.easyROI( GTV )
#' }#' 
#' @export
RAD.easyROI<-function( ROIVoxelData ) {
  res<-list();
  for(i in names(ROIVoxelData)) {
    if((TRUE %in% is.na(ROIVoxelData[[i]])) == FALSE  ) {
      res[[i]]<-ROIVoxelData[[i]]$masked.images$voxelCube
    } 
    else {
      res[[i]]<-NA
    }
  }
  return(res);
}




RAD.getNeighbourhood<-function( ROIVoxelData, margin.x = 1, margin.y = 1, margin.z = 0, sampleSize = 1000) {
  
  wErode<-RAD.applyErosion(  ROIVoxelData, margin.x=margin.x+1, margin.y=margin.y+1, margin.z=margin.z );
  
  totalMatrix<-list();
  submatrix<-c(); id<-1;
  for(patientName in names(wErode)) {
    coords<-which( !is.na(wErode[[ patientName ]]$masked.images$voxelCube),arr.ind = TRUE  )
    cat("\n extracting ",patientName,'......');
    sliceNumber<-as.integer(dim(wErode[[ patientName ]]$masked.images$voxelCube)[3]/2)
    riga<-sliceNumber;
    matriceCentroidi<-coords[which(coords[,3]==sliceNumber),]
    matriceCentroidi<-matriceCentroidi[ sample( seq(1,nrow(matriceCentroidi)),size = min(sampleSize, nrow(matriceCentroidi)),replace = FALSE)   ,]
    coords<-matriceCentroidi

    matriceAddendum<-c();
    for(riga in seq(1,dim(matriceCentroidi)[1] )) {

      from.x<-(coords[riga,1]-margin.x);  to.x<-(coords[riga,1]+margin.x);
      from.y<-(coords[riga,2]-margin.y);  to.y<-(coords[riga,2]+margin.y);
      from.z<-(coords[riga,3]-margin.z);  to.z<-(coords[riga,3]+margin.z);
      matriceCoordsAround<-expand.grid( seq(from.x,to.x),  seq(from.y,to.y), seq(from.z,to.z) );

      submatrix<-GTV[[patientName]]$masked.images$voxelCube[ (coords[riga,1]-margin.x):(coords[riga,1]+margin.x),
                                                              (coords[riga,2]-margin.y):(coords[riga,2]+margin.y),
                                                              (coords[riga,3]-margin.z):(coords[riga,3]+margin.z)   ]
      matriceCoordsAround<-cbind(
                                rep(id,nrow(matriceCoordsAround)),
                                rep(coords[riga,1],nrow(matriceCoordsAround)),
                                rep(coords[riga,2],nrow(matriceCoordsAround)),
                                rep(coords[riga,3],nrow(matriceCoordsAround)),
                                matriceCoordsAround, array(submatrix))
#      matriceAddendum<-rbind(matriceAddendum,matriceCoordsAround)
#      totalMatrix[[patientName]]<-rbind(totalMatrix[[patientName]],matriceAddendum)
#      matriceAddendum<-rbind(matriceAddendum,matriceCoordsAround)
      totalMatrix[[patientName]]<-rbind(totalMatrix[[patientName]],matriceCoordsAround)
      
      id<-id+1
    }
    
    colnames(totalMatrix[[patientName]])<-c("id","x.centr","y.centr","z.centr","x.neigh","y.neigh","z.neigh","value")
  }
  cat("\n");
  return(totalMatrix)
}


#' return the biopsy 
#' 
#' @description  get the list of possible centroids and return the list of biopsy
#' @param possBio as got from of a \code{RAD.VirtualBiopsy()} method.
#' @param ROIVoxelData is the voxelData (cropped) got from of a \code{obj$getROIVoxel()} method.
#' @param x the x dim of the biopsy (as passed to the \code{RAD.VirtualBiopsy()} method: the real cube will have (2*x+1) along x axes)
#' @param y the x dim of the biopsy (as passed to the \code{RAD.VirtualBiopsy()} method: the real cube will have (2*y+1) along y axes) 
#' @param z the x dim of the biopsy (as passed to the \code{RAD.VirtualBiopsy()} method: the real cube will have (2*z+1) along z axes)  
#' @return a list (quite big) containing all the biopsy
#' @examples \dontrun{
#' # Create an instante of new.mmButo and load some cases
#' obj<-new.mmButo()
#' obj$loadCollection(Path = '/progetti/immagini/urinaEasy')
#' 
#' # get the three ROIs
#' GTV<-obj$getROIVoxel(ROIName="GTV")  
#' 
#' # get the possible centroids
#' a<-RAD.VirtualBiopsy(ROIVoxelData = GTV,dx.min = 3,dy.min = 4,dz.min = 1,dx.max = 4,dy.max = 4,dz.max = 1)
#' 
#' # get the biopsy
#' listaBiopsie<-RAD.getBiopsy(possBio = a,ROIVoxelData = GTV,x = 4,y = 4,z = 1)
#' }#' 
#' @export
RAD.getBiopsy<-function(possBio, ROIVoxelData, x = 4, y = 4, z = 1) {
  obj.mmButo<-new.mmButo()
  submatrix<-list()
  stringaDim<-paste(c(x,y,z),collapse='_');
  #scorri tutti i pazienti
  for(paziente in names(possBio)) {
    print(paziente);
    if((TRUE %in% is.na(possBio[[paziente]])) == FALSE  ) {
      submatrix[[paziente]]<-list()
      # prendi la lista dei centroidi per fare la biopsia
      listaPunti<-possBio[[ paziente ]][[ stringaDim ]]$IndexBiopsy
      # espandi il voxelCube
      voxelCube <- obj.mmButo$mmButoLittleCube.expand(   ROIVoxelData[[ paziente ]] )
      # per ogni centroide estrai la corrispondente biopsia
      for( riga in seq(1,dim(listaPunti)[1])) {
        # coordinate del centroide
        xC<-listaPunti[riga,1]; yC<-listaPunti[riga,2]; zC<-listaPunti[riga,3];
        # carota
        submatrix[[paziente]][[riga]]<-voxelCube[  (xC-x):(xC+x),(yC-y):(yC+y),(zC-z):(zC+z)   ]
        # se c'è anche solo uno '0', 'epic fail!'
        if( length(which(submatrix[[paziente]][[riga]]==0)) !=0 ) { 
          print("Epic fail");
          browser();
        }
      }
    }
    else {
      submatrix[[paziente]]<-NA
    }
  }
  return(submatrix);
}
#' Calculate entropy and stdev maps
#' 
#' @description  Calculate entropy maps and stdev maps of the images, eventually normalizing the image
#' @param obj.mmButo an object of class \code{new.mmButo}
#' @param ROIName the name of the interested ROI
#' @param margin.x for the biopsy; default is 3 (the box has 7 voxels along x direction)
#' @param margin.y for the biopsy; default is 3 (the box has 7 voxels along y direction)
#' @param margin.z for the biopsy; default is 1 (the box has 3 voxels along z direction)
#' @param erosion.x for the erosion; default is 3 
#' @param erosion.y for the erosion; default is 3 
#' @param erosion.z for the erosion; default is 1 
#' @param collection the interested collection; the dafault collection is \code{default}
#' @param normalizationROIName the name of the ROI used to normaliza the signal: default is \code{NA}
#' @param kindOfNormalization the kind of normalization: \code{'linear'} or \code{'log'}. Default is \code{'linear'}
#' @param kindOfOutput is a string which can be '\code{normal}' or '\code{extended}'. In the first case it returna a complex list, cropped, with all the information needed to extend the images; in the second it provides the extended images avoiding to return all the ancillary (geometry) data
#' @return a list containing the entropy maps and the stdev maps.
#' @examples \dontrun{
#' # Create an instante of new.mmButo and load some cases
#' obj<-new.mmButo()
#' obj$loadCollection(Path = '/progetti/immagini/urinaEasy')
#' 
#' # get the three ROIs
#' Retto<-obj$getROIVoxel(ROIName="Retto")  
#' 
#' # get the possible biopsy
#' erodedCubes<-RAD.applyErosion( ROIVoxelData = Retto )
#' }#' 
#' @import entropy  
#' @export
RAD.borderTextureMap<-function(obj.mmButo, ROIVoxelData, margin.x=3,margin.y=3,margin.z=1,collection="default", ROINameForNormalization = NA, typeOfCorrection='linear', erosion.x = 3,erosion.y = 3,erosion.z = 1, kindOfOutput="normal", usingCache = FALSE) {
  objS<-services();
#   # prendi le matrici dei voxel completi della ROI di interesse (croppati)
#   ROIVoxelData<-obj.mmButo$getROIVoxel(ROIName = ROIName)
#   if(!is.na(ROINameForNormalization)) {
#     ROIForCorrection<-obj.mmButo$getROIVoxel(ROIName = ROINameForNormalization)
#     ROIVoxelData<-obj.mmButo$getCorrectedROIVoxel(inputROIVoxel = ROIVoxelData,correctionROIVoxel = ROIForCorrection, typeOfCorrection = typeOfCorrection)
#   }
  # calcola l'erosione
  eroded<-RAD.applyErosion(ROIVoxelData = ROIVoxelData,margin.x = erosion.x,margin.y = erosion.y,margin.z = erosion.z)

  # prendi la lista di oggetti geoLet
  list_geoLet<-obj.mmButo$getAttribute("list_geoLet");
  entropyMap<-list();
  standardDev<-list();
  ct<-1
  for(patient in names(ROIVoxelData)) {
    
    print( paste("ENTROPY - Now processing:",patient),collapse='' )
    
    # calcola il bordo (come la differenza)
    bordo<-ROIVoxelData[[patient]]$masked.images$voxelCube - eroded[[patient]]$masked.images$voxelCube
    # prendi la struttura di VoxelData
    voxelData.bordo<-ROIVoxelData[[patient]]
    # ed ad essa sostituisci l'immagine del bordo
    voxelData.bordo$masked.images$voxelCube<-bordo
    # così' sono pronto per estendere
    bordoEspanso <- obj.mmButo$mmButoLittleCube.expand(   voxelData.bordo )
    # i punti non a zero di tale struttura sono quindi quelli in cui posso "carotare".
    centr<-which(bordoEspanso!=0,arr.ind = T)
    if ( usingCache == TRUE) list_geoLet[[collection]][[patient]]$cacheLoad()
    originalMR<-list_geoLet[[collection]][[patient]]$getImageVoxelCube()

    entropyMap[[patient]]<-array(0,dim=c(  dim(originalMR)[1],dim(originalMR)[2],dim(originalMR)[3]   ))
    standardDev[[patient]]<-array(0,dim=c(  dim(originalMR)[1],dim(originalMR)[2],dim(originalMR)[3]   ))
    
    for ( i in seq(1,dim(centr)[1]) ) {
      
      if( 
        (centr[i,][3]+margin.z)<=dim(originalMR)[3] && (centr[i,][3]-margin.z)>=1 && 
        (centr[i,][2]+margin.y)<=dim(originalMR)[2] && (centr[i,][2]-margin.y)>=1 &&
        (centr[i,][1]+margin.x)<=dim(originalMR)[1] && (centr[i,][1]-margin.x)>=1
      ) {      
        subCube<-originalMR[  
          (centr[i,][1]-margin.x):(centr[i,][1]+margin.x), 
          (centr[i,][2]-margin.y):(centr[i,][2]+margin.y), 
          (centr[i,][3]-margin.z):(centr[i,][3]+margin.z) ]
        
        istogr <- hist(subCube,  plot=FALSE)    
        freq <- freqs(y = istogr$counts)
        ImageEntropy <- entropy.plugin(freqs = freq, unit = c("log2"))
        entropyMap[[patient]][ centr[i,][1] , centr[i,][2], centr[i,][3] ]<-ImageEntropy
        standardDev[[patient]][ centr[i,][1] , centr[i,][2], centr[i,][3] ]<-sd(subCube)
      }
    }
    # ora Croppali, o facciamo notte
    entropyMap[[patient]]<-objS$cropCube( entropyMap[[patient]] )
    standardDev[[patient]]<-objS$cropCube( standardDev[[patient]]  )
    ct<-ct+1
  }
  if( kindOfOutput == "normal") {
    return( list("entropyMap"=entropyMap,"stdMap"=standardDev)   )
  }
  if( kindOfOutput == "extended") {
    entropyExtended<-ROIVoxelData
    standardDevExtended<-ROIVoxelData
    ct<-1
    listaEntropy<-list();
    listaSD<-list();
    for( patient in names(ROIVoxelData)) {
      entropyExtended[[patient]]$masked.images$voxelCube<-entropyMap[[patient]]$voxelCube
      standardDevExtended[[patient]]$masked.images$voxelCube<-standardDev[[patient]]$voxelCube
      listaEntropy[[patient]]<-obj.mmButo$mmButoLittleCube.expand(  entropyExtended[[patient]] )
      listaSD[[patient]]<-obj.mmButo$mmButoLittleCube.expand(  standardDevExtended[[patient]] )
      ct<-ct+1
#      if(ct==5) break;   # just for debug      
    }
    return( list("entropyMap"=listaEntropy,"stdMap"=listaSD)   )
  }  
  
}
#' Apply a user defined function over a defined region of mmButo voxel maps
#' 
#' @description  It applies a user defined \code{function} over a voxel maps of in a \code{new.mmButo} object. The maps can be the normal maps or maps obtained by erosion and/or biopsy.
#' @param obj.mmButo an object of class \code{new.mmButo}. Ignored if \code{ROIName} and \code{normalizationROIName} are not character but \code{new.mmButo$getROIVoxel()} object types
#' @param ROIVoxelData the ROIVoxelData as obtained from a \code{getROIVoxel()} or a \code{getCorrectedROIVoxel}
#' @param collection the interested collection; the dafault collection is '\code{default}'
#' @param biopsyDim.xyz is an integer numeric vector with the \code{c(x,y,z)} dimensions in voxels of the cube for biopsy. The default is \code{NA} which means that no biopsy is perfomed
#' @param erosion.xyz is an integer vector with the \code{c(x,y,z)} dimensions of the erosion which has to be applied to the voxelCubes respect the \code{ROIName}
#' @param applyToBorder a logical value, if \code{TRUE} the function \code{FUN} is applied to the border, if \code{FALSE} the FUNCTION is applied to all the voxel within the \code{ROIName} (eventually eroded)
#' @param FUN is the \code{function} to be be applied to any biopsy (if a \code{biopsyDim.xyz}) is given, otherwise the function is applied to all voxels
#' @param cropResult a logical value. If \code{TRUE} the result will be cropped and will occupy much less RAM, otherwise it will have the original size (default value is \code{TRUE})
#' @param ... Other parameters to be passed to the function \code{FUN}
#' @return a list containing the output of the function \code{FUN} applied over biopsies or over the whole ROI voxels series#'
#' @examples \dontrun{
#' # Create an instante of new.mmButo and load some cases
#' obj<-new.mmButo()
#' obj$loadCollection(Path = '/progetti/immagini/urinaEasy')
#' 
#' Urina<-obj$getROIVoxel(ROIName = "Urina")
#' GTV<-obj$getROIVoxel(ROIName = "GTV")
#' GTV.C <- obj$getCorrectedROIVoxel(inputROIVoxel = GTV,correctionROIVoxel = Urina,typeOfCorrection = "log")
#' 
#' calcolaSD<-function(a, listOfFUNParameters=list()) {	
#'  b<-sd(a);
#'  if(listOfFUNParameters$noise==TRUE)  b<-b+b*runif(1)
#'  return(b)
#' }  
#' 
#' # no noise, not BORDER, on biopsy of <2,2,0> and erosion of <2,2,0>
#' a<-RAD.imagesApplyFUN(obj.mmButo = obj,ROIVoxelData = GTV.C,FUN = calcolaSD,erosion.xyz = c(5,5,3),applyToBorder = FALSE,listOfFUNParameters = list("noise"=FALSE),biopsyDim.xyz = c(2,2,0))
#' 
#' }#' 
#' @export
RAD.imagesApplyFUN<-function(obj.mmButo, ROIVoxelData, collection="default",
                             biopsyDim.xyz=NA,erosion.xyz = NA,applyToBorder = FALSE, FUN = NA, cropResult = TRUE, ... ) {
  if(applyToBorder==TRUE && is.na(erosion.xyz) ) stop("Provide a valid 'erosion.xyz' value to apply a function on the border voxels")
  if( !is.function(FUN)) stop("No function provided")  
  objS<-services();

  # calcola l'erosione (se necessario)
  eroded<-NA
  if(!is.na(erosion.xyz[1]) && !is.na(erosion.xyz[2]) && !is.na(erosion.xyz[3])) { 
    eroded<-RAD.applyErosion(ROIVoxelData = ROIVoxelData,margin.x = erosion.xyz[1],margin.y = erosion.xyz[2],margin.z = erosion.xyz[3])
  }
  
  # prendi la lista di oggetti geoLet
  list_geoLet<-obj.mmButo$getAttribute("list_geoLet");
  ct<-1
  FUNMap<-list();
  for(patient in names(ROIVoxelData)) {
  
    if((TRUE %in% is.na(ROIVoxelData[[patient]])) == FALSE  ) {
    
      print( paste("FUN - Now processing:",patient),collapse='' )
 
      voxelData.ready<-ROIVoxelData[[patient]]
      
      # calcola il bordo (come la differenza)
      if(applyToBorder==TRUE) {
        bordo<-ROIVoxelData[[patient]]$masked.images$voxelCube - eroded[[patient]]$masked.images$voxelCube
        # prendi la struttura di VoxelData
        # ed ad essa sostituisci l'immagine del bordo
        voxelData.ready$masked.images$voxelCube<-bordo      
      }
      if(applyToBorder==FALSE && !is.na(eroded)  ) {
        voxelData.ready$masked.images$voxelCube<-eroded[[patient]]$masked.images$voxelCube
      }
      
      # ora sono pronto per estendere
      voxelData.ready.espanso <- obj.mmButo$mmButoLittleCube.expand(   voxelData.ready   )
      # ora vediamo se c'è da carotare
      # i punti non a zero di tale struttura sono quindi quelli in cui posso "carotare".
      centr<-which(voxelData.ready.espanso!=0,arr.ind = T)
      originalMR<-list_geoLet[[collection]][[patient]]$getImageVoxelCube()
      
      # prepara la FUNMap
      FUNMap[[patient]]<-array(0,dim=c(  dim(originalMR)[1],dim(originalMR)[2],dim(originalMR)[3]   ))
      for ( i in seq(1,dim(centr)[1]) ) {
        if(!is.na(biopsyDim.xyz[1]) && !is.na(biopsyDim.xyz[2]) && !is.na(biopsyDim.xyz[3]) ) {
          margin.x<-biopsyDim.xyz[1]; margin.y<-biopsyDim.xyz[2]; margin.z<-biopsyDim.xyz[3]
        }
        else {margin.x=0; margin.y=0; margin.z=0;}
        if( 
          (centr[i,][3]+margin.z)<=dim(originalMR)[3] && (centr[i,][3]-margin.z)>=1 && 
          (centr[i,][2]+margin.y)<=dim(originalMR)[2] && (centr[i,][2]-margin.y)>=1 &&
          (centr[i,][1]+margin.x)<=dim(originalMR)[1] && (centr[i,][1]-margin.x)>=1
        ) {      
          subCube<-originalMR[  
            (centr[i,][1]-margin.x):(centr[i,][1]+margin.x), 
            (centr[i,][2]-margin.y):(centr[i,][2]+margin.y), 
            (centr[i,][3]-margin.z):(centr[i,][3]+margin.z) ]
          FUNMap[[patient]][ centr[i,][1] , centr[i,][2], centr[i,][3] ]<-FUN(subCube , ...)
        }
      }
     
      # se richiesto, CROPPA!
      if ( cropResult == TRUE )  FUNMap[[patient]]<-objS$cropCube( FUNMap[[patient]] )    
    }
    else {
      FUNMap[[patient]]<-NA
    }
  }
  return(FUNMap)
}