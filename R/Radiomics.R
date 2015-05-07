#' Calculates image voxels internal respect a given ROI
#' @description This function can be used to calculate internal voxels, for each ROIs, in a \code{dataStorage} data structure taken from a \code{geoLet} object
#' @param dataStorage is the structure returned from a \code{obj$getAttribute("dataStructure")} where \code{obj} is an instance of a \code{geoLet}
#' @param Structure is the ROI name od the ROI that should be extracted
#' @param SeriesInstanceUID Is the interested series Instance UID.
#' @return A list with unknown meaning
#' @useDynLib moddicom
#' @export
RAD.NewMultiPIPOblique<-function(dataStorage, Structure, SeriesInstanceUID) {
  objService<-services()
  
  # define some variables to make more clear the code
  numberOfRows<-as.numeric(dataStorage$info[[SeriesInstanceUID]][[1]]$Rows);
  numberOfColumns<-as.numeric(dataStorage$info[[SeriesInstanceUID]][[1]]$Columns);
  numberOfSlices<-length(dataStorage$img[[SeriesInstanceUID]]);
  
  # initialize the image array with the right dimension
  image.arr<-array( data = -1, dim = c(numberOfRows, numberOfColumns, numberOfSlices ) )
  # index values listed as characters: creates empty array of DICOM orientation matrices
  index<-as.character(sort(as.numeric(names( dataStorage$img[[SeriesInstanceUID]]) )))  
  
  # create and fill the vector DOM, of DICOM orientation matrices
  DOM<-c();  nn<-0
  for (n in index) {
    DOM<-c(DOM, dataStorage$info[[SeriesInstanceUID]][[n]]$orientationMatrix[c(1:3,5:7,13:15)])
    nn<-nn+1
    image.arr[,,nn]<-dataStorage$img[[SeriesInstanceUID]][[n]]    
  }  
  
  # fills the vectors of X and Y coordinates 
  # and other Vectors 'associatedInstanceNumberVect' and 'arrayInstanceNumberWithROI'
  TotalX<- -10000;  TotalY<- -10000;  arrayAssociationROIandSlice<- -10000;
  OriginX<- -10000;   OriginY<- -10000
  associatedInstanceNumberVect<- -10000

  contatoreROI<-1; indiceDOM<-1;
  # for each instance number
  for (n in index) {
    # check if there is a ROI for such slice
    for (m in which(dataStorage$info[[SeriesInstanceUID]][[n]]$ROIList[,1]==Structure)) {
      
      # find the slice and gets the key for accessing at coordinates vectors
      key<-dataStorage$info[[SeriesInstanceUID]][[n]]$ROIList[m,2]
      
      # calculate how many ROIs are co-planar
      numeroROIComplanari<-length(dataStorage$structures[[Structure]][[key]])
      
      # for each one of them concat the array
      for(indiceROI in seq(1,numeroROIComplanari)) {          
        TotalX<-c(TotalX, dataStorage$structures[[Structure]][[key]][[indiceROI]][,1])
        TotalY<-c(TotalY, dataStorage$structures[[Structure]][[key]][[indiceROI]][,2])
        
        # calculate how many points compose the ROI
        numeroPunti<-length(dataStorage$structures[[Structure]][[key]][[indiceROI]][,1])
        
        # for each point write which is the related Slice in the cube-matrix
        arrayAssociationROIandSlice<-c(arrayAssociationROIandSlice,rep(   which( index == n ) -1  , numeroPunti ))
        
        # Usa OriginX and OriginY as terminator
        TotalX<-c(TotalX, OriginX)
        TotalY<-c(TotalY, OriginY)      
        arrayAssociationROIandSlice<-c(arrayAssociationROIandSlice,OriginX)
        
        contatoreROI<-contatoreROI+1      
      } 
    }      
    indiceDOM<-indiceDOM+1;
  }

  # ok, call the Wrapper!
  final.array<-RAD.NewMultiPointInPolyObl(
    # array of DICOM Orientation Matrices
    DICOMOrientationVector = DOM, 
    # X and Y vector Points
    totalX = TotalX, totalY = TotalY, 
    # association between ROIs and Slices in the 3D Matrix
    arrayAssociationROIandSlice = arrayAssociationROIandSlice,
    # matrices dimensions (rows and columns)
    nX = numberOfColumns, 
    nY = numberOfRows,
    nZ = numberOfSlices
  )
  
  final.array<-array(data = final.array, dim = c(   numberOfColumns, numberOfRows, numberOfSlices )   )
    # In Example: 
    #
    # > TotalX[0:70]
    # [1] -10000.00     -5.16     -3.28     -1.41      0.47      2.34      4.22      5.12      6.09      7.61      7.97      9.48      9.84     10.96     11.72
    # [16]     12.14     12.93     13.59     13.61     14.19     14.70     14.85     14.23     13.59     13.45     12.68     11.72     11.11      9.84      8.78
    # [31]      7.97      6.09      4.22      2.84      2.34      0.85      0.47      0.14     -1.41     -3.28     -3.92     -5.16     -5.73     -6.56     -7.03
    # [46]     -8.91    -10.26    -10.78    -11.25    -11.86    -12.01    -11.83    -11.29    -10.78    -10.64     -9.88     -8.91     -8.81     -7.51     -7.03
    # [61]     -5.46     -5.16 -10000.00     -5.16     -3.28     -1.41      0.47      2.34      4.22      5.38
    # 
    # > arrayAssociationROIandSlice[0:70]
    # [1] -10000      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8
    # [22]      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8
    # [43]      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      9 -10000
    # [64]      9      9      9      9      9      9      9
  
    # ROTATE THE MATRIX
    for ( i in seq(1,dim(image.arr)[3] )) {
      image.arr[,,i]<-objService$SV.rotateMatrix(image.arr[,,i])
      final.array[,,i]<-t(objService$SV.rotateMatrix(final.array[,,i],rotations=3))
    }
  
#    return(list(TotalX=TotalX, TotalY=TotalY, FullZ=FullZ, Offset=Offset, 
#                DOM=array(DOM, dim = c(3,3,length(index))), final.array=final.array, masked.images=final.array*image.arr))
    return(list("DOM"=array(DOM, dim = c(3,3,length(index))), "final.array"=final.array, "masked.images"=final.array*image.arr))
}
# ==========================================================================================================================
#' Wrapper for C function
#' @useDynLib moddicom
#RAD.NewMultiPointInPolyObl<-function(totalX, totalY, nX, nY, associatedInstanceNumberVect,
                                     #NumSlices, NumSlicesWithROI, 
                                     #arrayInstanceNumber, arrayInstanceNumberWithROI, arrayPosizioneInstanceNumberWithROI,
                                     #DICOMOrientationVectorWithROI, DICOMOrientationVector,minX,maxX,minY,maxY ) {
RAD.NewMultiPointInPolyObl<-function(DICOMOrientationVector,totalX,totalY,arrayAssociationROIandSlice,nX,nY,nZ ) {  
  
  maxX<-max(totalX)
  minX<-min(totalX[which(totalX>-10000)])
  maxY<-max(totalY)
  minY<-min(totalY[which(totalY>-10000)])
  
  # creates the PIPvector
  PIPvector<-rep.int(x = 0, times = nX * nY * nZ)  
  numberOfPoints<-length(totalX);
  result<-.C("NewMultiPIPObl", 
             as.integer(PIPvector), as.double(totalX), as.double(totalY), as.integer(numberOfPoints), 
             as.integer(nX), as.integer(nY), as.integer(nZ),             
             as.integer(arrayAssociationROIandSlice), 
             as.double(DICOMOrientationVector),as.double(minX),as.double(maxX),as.double(minY),as.double(maxY))  
  
  return(result[[1]])
}


#' Load and handle a tree
#' @description this is the old (and deprecated) version of mmButo. It is able to look into a folder serching all DICOM studies and return a structure containint all the voxel cubes.
#'               The available methods are:
#'               \itemize{
#'               \item \code{openTreeMultiROIs(Path, structureList)} 
#'               is a method used to open a chosen folder. This method loads all the DICOM objects into
#'               the indicated folder (without recursion) as attribute of the object.
#'               Information can be retrieved using \code{getAttribute} method or, simply, getting the result of this method.#'               
#'               }
#' @import MASS colorRamps 
#' @export
RAD.mmButo<-function() {
  
  dataStructure<-list()
  attributeList<-list()
  logObj<-list()
  arrayAR<-list()           # arrayAlgorithmResult

  # ========================================================================================
  # openTreeMultiROIs
  # it load all the subfolder of a given folder searching for all the available studyes
  # ========================================================================================  
  openTreeMultiROIs<-function(Path, structureList) {
    cubeVoxelList<-list()
    listaFolders<-list.dirs( Path )
    
    counter<-1
    for( folderName in listaFolders[2:length(listaFolders)] ) {
      #      setTxtProgressBar(pb, counter)
      # Instantiate the object
      obj<-geoLet()
      obj$openDICOMFolder( folderName )      
      listaROI<-obj$getROIList();
      
      cubeVoxelList[[ folderName ]]<-list()
      cubeVoxelList[[ folderName ]][["voxelCubes"]]<-list()
      cubeVoxelList[[ folderName ]][["info"]]<-list()
      cubeVoxelList[[ folderName ]][["ROIPointList"]]<-list()
      
      ct<-1;
      for( ROIName in obj$getROIList()[2,] ) {
        if ((missingArg(structureList)) ||  (ROIName %in% structureList)  ) {
          ds<-obj$getAttribute("dataStorage")
          ROIPointList<-obj$getROIPointList(ROIName)
          SS<-names(ds$img);
          
          result<-RAD.NewMultiPIPOblique(dataStorage = ds, Structure = ROIName, SeriesInstanceUID = SS)
          
          cubeVoxelList[[ folderName ]][["voxelCubes"]][[ ROIName ]]<-result$masked.images    
          cubeVoxelList[[ folderName ]][["info"]]<-result$DOM
          
          cubeVoxelList[[ folderName ]][["geometricData"]]<-list();
          cubeVoxelList[[ folderName ]][["geometricData"]]$pixelSpacing<-ds$info[[1]][[1]]$pixelSpacing
          cubeVoxelList[[ folderName ]][["geometricData"]]$SliceThickness<-as.numeric(ds$info[[1]][[1]]$SliceThickness)
          cubeVoxelList[[ folderName ]][["ROIPointList"]][[ ROIName ]]<-ROIPointList
                  
          info_struct<-ds$info[[SS]]                            # structure which contains information
          img_struct<-ds$img[[SS]]                              # structure which contains images
          nnRows<-as.numeric(info_struct[[1]]$Rows)                                     # Rows
          nnColumns<-as.numeric(info_struct[[1]]$Columns)                               # Columns
          nnOfSlices<-length(img_struct)                                                # Number of slices
          instanceNumberList<-as.character( sort( as.numeric( names(img_struct) ) ) )   # List of the instance number (key for img_struct)
          
          image.arr<-array(data=-1, dim = c( nnRows, nnColumns, nnOfSlices ) )
          nn<-0
          for (n in instanceNumberList) {
            nn<-nn+1
            image.arr[,,nn]<-img_struct[[n]]    
          }
          image.arr=aperm(a = image.arr, perm = c(2,1,3))[,nnColumns:1,]
          cubeVoxelList[[ folderName ]][["image.arr"]]<-image.arr
        }
      }  
      
      counter<-counter+1;
    }
    
    dataStructure<<-cubeVoxelList      
  }
  # ========================================================================================
  # ROIStats
  # return stats of a ROI
  # ========================================================================================
  ROIStats<-function(ROIName) {
    res<-list();
    res$details<-list();
    res$total<-list();
    for(i in names(dataStructure)) {
      res$details$stdev[[i]]<-sd(dataStructure[[i]]$voxelCubes$Urina[which(dataStructure[[i]]$voxelCubes$Urina!=0)])
      res$details$mean[[i]]<-mean(dataStructure[[i]]$voxelCubes$Urina[which(dataStructure[[i]]$voxelCubes$Urina!=0)])
    }
    res$total$mean<-mean(res$details$mean)
    res$total$min<-min(res$details$mean)
    res$total$max<-max(res$details$mean)
    return(res);    
  }
  # ========================================================================================
  # execAlgorithm
  # execute an algorithm on the loaded DICOM studies.
  # ========================================================================================
  execAlgorithm<-function( algorithm, ROIName , grayTuniningValue, ROIName4Tuning , additionalParams=c()) {
    
    XmaxVal<-c() 
    YmaxVal<-c() 
    array4VoxelCube<-list();
    normalizedVoxelCube<-list();

    # BAVA
    if( algorithm == "BAVA" ) {      
      interpolatedDensity<-list();
      if(length(grayTuniningValue)==0) logObj$sendLog("'grayTuniningValue' is missing", NMI = TRUE)
      if(length(ROIName4Tuning)==0) logObj$sendLog("'ROIName4Tuning' is missing", NMI = TRUE)
      
      # massimo valore di grigio della ROI di tutti i pazienti 
      # (lo uso per sapere la X massima per il resampling)
      maxGreyOfROIName<-ROIStats(ROIName=ROIName)
      maxGreyOfROIName<-maxGreyOfROIName$total$max+1; 
      valoreMassimoX<-c()
      valoreMassimoY<-c()
      # loop for each Series Instance UID
      for( i in names(dataStructure)) {   

        # Prendi il valore medio della vescita del paziente in esame
        maxSpecificROI4Tuning<-mean(dataStructure[[i]]$voxelCubes[[ROIName4Tuning]][which(dataStructure[[i]]$voxelCubes[[ROIName4Tuning]]!=0)])
        # Normalizza
        normalizedVoxelCube[[i]] <- dataStructure[[i]]$voxelCubes[[ROIName]] * ( grayTuniningValue / maxSpecificROI4Tuning)
        
        #voxelCube <- dataStructure[[i]]$voxelCubes[[ROIName]]
        voxelCube <- normalizedVoxelCube[[i]]
        array4VoxelCube[[i]] <- voxelCube[which(  voxelCube  !=0 )];
        # passa alla KDF così ho la versione continua
        array4VoxelCube[[i]]<-density(    array4VoxelCube[[i]]   )        
        # memorizza il valore massimo della x (per poter poi fare il resampling sulla stessa scala)
        valoreMassimoX<-c(valoreMassimoX,array4VoxelCube[[i]]$x)
        valoreMassimoY<-c(valoreMassimoY,array4VoxelCube[[i]]$y)
      }       
      # ora fai i resampling
      for( i in names(dataStructure)) { 
        interpolatedDensity[[i]]<-approx(array4VoxelCube[[i]]$x,array4VoxelCube[[i]]$y,n=valoreMassimoX, xout=seq( from=0 , to=max(valoreMassimoX) )) 
        # azzera gli NA
        interpolatedDensity[[i]]$y[which(is.na(interpolatedDensity[[i]]$y))]<-0        
      }
      # write the results in the array
      arrayAR$BAVA<<-list();
      arrayAR$BAVA$details<<-list();
#      arrayAR$BAVA$details$density<<-interpolatedDensity     
      arrayAR$BAVA$details$interpolatedD<<-interpolatedDensity     
      arrayAR$BAVA$summary<<-list();
      arrayAR$BAVA$summary$XmaxVal<<-max( valoreMassimoX )
      arrayAR$BAVA$summary$YmaxVal<<-max( valoreMassimoY )      
      return();
    }

    # Kernel Density Function
    if( algorithm == "KDF" ) {      
      arrayAR[["KDF"]]<<-list();      
      arrayAR[["KDF"]]$details<<-list();
      arrayAR[["KDF"]]$summary<<-list(); 
      interpolatedDensity<-list();
      # loop for each Series Instance UID
      for( SeriesInstanceUID in names(dataStructure)) {          
        XmaxVal <- c(XmaxVal,max(dataStructure[[SeriesInstanceUID]]$image.arr));      
        # get only the 'not zero' voxel value
        voxelCube <- dataStructure[[SeriesInstanceUID]]$voxelCubes[[ROIName]]
        array4VoxelCube[[SeriesInstanceUID]] <- voxelCube[which(  voxelCube  !=0 )];
      } 
      # NORMALIZATION (for upper bound only)
      for( SeriesInstanceUID in names(dataStructure)) {
       tmpArr<-as.array(array4VoxelCube[[SeriesInstanceUID]])
       maxRMN<-max(  dataStructure[[SeriesInstanceUID]]$image.arr  )
       array4VoxelCube[[SeriesInstanceUID]]<-density(    tmpArr /  max(tmpArr)     )
       interpolatedDensity[[SeriesInstanceUID]]<-approx(array4VoxelCube[[SeriesInstanceUID]]$x,array4VoxelCube[[SeriesInstanceUID]]$y,n=100,xout=seq(from=0,to=1,by = .01))
       YmaxVal <- c(YmaxVal,max(array4VoxelCube[[SeriesInstanceUID]]$y));    # sospetto errore.   
      }
      # write the results in the array
      arrayAR$KDF$details$density<<-array4VoxelCube
      arrayAR$KDF$details$interpolatedD<<-interpolatedDensity

      arrayAR$KDF$summary<<-list();
      arrayAR$KDF$summary$XmaxVal<<-max( XmaxVal )
      arrayAR$KDF$summary$YmaxVal<<-max( YmaxVal )
      return();
    }
    #  CARLOTTAGGIO
    if( algorithm == "virtualBiopsy" ) { 
      arrayAR$Biopsy<<-list();  
      arrayAR$Biopsy$results<<-list()
      arrayAR$Biopsy$results<<-allPopulationVirtualBiopsy( nx=2,ny=2,nz=0, ROIName4Normalization=ROIName4Tuning, normalization=TRUE, grayTuniningValue = grayTuniningValue)
      return();
    }
    # AREA/VOLUME
    if( algorithm == "rawAreaVolume" ) { 
      objS<-services()
      arrayAR$AreaVolume<<-list()
      for ( i in names(dataStructure)) {
        pSX<-dataStructure[[ i ]]$geometricData$pixelSpacing[[1]]
        pSY<-dataStructure[[ i ]]$geometricData$pixelSpacing[[2]]
        pSZ<-dataStructure[[ i ]]$geometricData$SliceThickness
        arrayAR$AreaVolume[[i]]$Area<<-objS$SV.rawSurface(voxelMatrix = dataStructure[[ i ]]$voxelCubes[[ROIName]], pSX = pSX, pSY=pSY,pSZ=pSZ)
        arrayAR$AreaVolume[[i]]$Volume<<-length(which(dataStructure[[ i ]]$voxelCubes[[ ROIName ]]!=0))*pSX*pSY*pSZ
      }
      return();
    }
    # Not yet implemented NMI error
    logObj$sendLog(message = "Not yet implemented", NMI = TRUE);
  }
  # ========================================================================================
  # plotResults
  # It plots the results of a chosen algorithm
  # ========================================================================================
  plotResults<-function( singleLines = TRUE, meanLine=TRUE, algorithm , ylim = c(), xlim = c(), add=FALSE, colMean = "blue", color = "red", xlab=c() , main=c(), meanDensityLine=2) {      
    
    if( algorithm == "KDF" | algorithm == "BAVA")  {
      ct<-1;
      addingMatrix<-c()
      if( algorithm == "KDF" ) {
        if(length(xlab)==0) xlab<-"Normalized greylevel Histogram"
        if(length(main)==0) main="Kernel Density Function"
      }
      if( algorithm == "BAVA") {
        if(length(xlab)==0) xlab<-"Normalized greylevel Histogram"
        if(length(main)==0) main="BAVA comparison"
      }
        
      for(pathName in names(arrayAR[[algorithm]]$details$interpolatedD) ) {        
        if(length(ylim)==0) ylim <- c(0,arrayAR[[algorithm]]$summary$YmaxVal)

        if( singleLines == TRUE) {
          if( ct == 1  & add==FALSE) {      
            if ( length(xlim) == 0 ) 
              plot( arrayAR[[algorithm]]$details$interpolatedD[[pathName]], ylim = ylim, col = color , xlab = xlab, main = main, type='l' ) 
            else 
              plot( arrayAR[[algorithm]]$details$interpolatedD[[pathName]], ylim = ylim, col = color , xlab = xlab, main = main, xlim = xlim, type='l' ) 
          }
          else 
            lines( arrayAR[[algorithm]]$details$interpolatedD[[pathName]] , col = color ) 
        }
        
        addingMatrix<-rbind(addingMatrix,arrayAR[[algorithm]]$details$interpolatedD[[pathName]]$y)        
        ct<-ct+1
      }
      
      if( meanLine == TRUE ) {
        
        if( singleLines == FALSE & add==FALSE) {
            xlim<-c( min(arrayAR[[algorithm]]$details$interpolatedD[[pathName]]$x) ,max(arrayAR[[algorithm]]$details$interpolatedD[[pathName]]$x)  )
            if ( length(xlim) == 0 ) 
              plot( x=c(), y=c(), ylim = ylim, col = color , xlab = xlab, main = main, type='l',lwd=meanDensityLine ) 
            else 
              plot(  x=c(), y=c(), ylim = ylim, col = color , xlab = xlab, main = main, xlim = xlim, type='l',lwd=meanDensityLine ) 
        }        
        
        # calculate the mean
#        stDevMatrix<-apply(addingMatrix, 2, sd)

        x<-arrayAR[[algorithm]]$details$interpolatedD[[pathName]]$x
        x[which(is.na(x))]<-0
        addingMatrix[which(is.na(addingMatrix))]<-0

        bootstrappedAddingMatrix<-apply(addingMatrix,2,sample,size=1000,replace=T)
        #quantileMatrix<-apply(addingMatrix, 2, quantile, probs = c(.025, .975), na.rm = TRUE)   # versione NON bootstrappata
        quantileMatrix<-apply(bootstrappedAddingMatrix, 2, quantile, probs = c(.025, .975), na.rm = TRUE)   # versione bootstrappato

        quantileMatrixSPU<-smooth.spline( x = x, y = quantileMatrix[1,])
        quantileMatrixSPL<-smooth.spline( x = x, y = quantileMatrix[2,])

        meanMatrix<-colMeans(addingMatrix, na.rm = TRUE)      # media normale
        
        
        # plot it
        polygon(c(x,rev(x)),c(meanMatrix,rev(quantileMatrixSPU$y)), col=rgb(.7, .7, .7, 0.2), lty = c("dashed"))
        polygon(c(x,rev(x)),c(meanMatrix,rev(quantileMatrixSPL$y)), col=rgb(.7, .7, .7, 0.2), lty = c("dashed"))
        lines( x = x , y=meanMatrix , col = colMean ) 
      }
    }
  } 
  bootstrapColMatrix<-function( x ){
    d<-sample(x=x,size=10000,replace=T)
    return(mean(d))
  }
  # ========================================================================================
  # virtualBiopsy
  # Calculates the position of elements which is possible to do virtual Biopsy
  # ======================================================================================== 
  #' Calculates the position of elements which is possible to do virtual Biopsy
  #' description: This function can be used to calculate the index of elements for virtual Biopsy along a given distance along x,y,z
  #' param: voxelCubes is the voxel space along x
  #' param: nx is the voxel space along x
  #' param: ny is the voxel space along y
  #' param: nz is the voxel space along z
  #' return: a list 

  virtualBiopsy <- function (voxelCubes,nx,ny,nz) {
    
    i<-0;    j<-0;    k<-0
    #voxelCubes <- ds.positive[[1]]$voxelCubes[[1]]
    exam <- array(voxelCubes,dim = dim(voxelCubes)[1]*dim(voxelCubes)[2]*dim(voxelCubes)[3])
    size1 <- (dim(voxelCubes)[1]);    size2 <- (dim(voxelCubes)[2]);    size3 <- (dim(voxelCubes)[3])
    cmp <- as.matrix(expand.grid(indiceX=seq(i-nx,i+nx),indiceY=seq(j-ny,j+ny),indiceZ=seq(k-nz,k+nz)))
    expand <- array(cmp, dim = dim(voxelCubes)[1]*dim(voxelCubes)[2])
    control <- (dim(cmp)[1])
    lungh <- length(exam)
    carotaggioVolume <- array(data = c(0), dim = dim(voxelCubes)[1]*dim(voxelCubes)[2]*dim(voxelCubes)[3])
    
    #load virtual biopsy C function to complete matrix with 1 and 0
    biopsy <- .C(   "virtualBiopsy", as.double (exam), as.integer (size1), as.integer (size2),
                    as.integer (size3), as.integer(nx), as.integer(ny), as.integer(nz), as.integer(expand),
                    as.integer (control), as.integer(lungh), as.integer(carotaggioVolume)   )
    
    
    #organize output into matrix structure
    virtual.biopsy <- array (biopsy[11][[1]], dim=c(size1,size2,size3))
    
    # create list of Grey of Virtual Biopsy
    index.carot <- which(virtual.biopsy == 1, arr.ind=T)   #indici di punti carotabili
    
    #input C grey function
    #initialitation
    indX <- index.carot[,1];    indY <- index.carot[,2];    indZ <- index.carot[,3]
    numCentroidi <- nrow(index.carot)
    lista.grigi2 <- data.frame()
    outputGrigi <- array(data = c(0), dim = control*numCentroidi)
    media <- c();    deviance <- c()
    
    # load C Grey function
    grigi <- .C(   "greyMatrix", as.integer(indX), as.integer(indY), as.integer(indZ), as.integer(numCentroidi),
                   as.integer(nx), as.integer(ny), as.integer(nz), as.integer(outputGrigi), as.double(exam),
                   as.integer (size1), as.integer (size2), as.integer (size3), as.integer (control)
    )
    #Grey array into matrix
    grey.matrix <- array (grigi[[8]], dim=c(control,numCentroidi))  
    
    # Evaluate mean and sd
    k <- 1
    for(j in 1:numCentroidi) {
      media[k] <- mean(grey.matrix[,j])
      deviance[k] <- sd (grey.matrix[,j])
      k <- k+1
    }
    
    #Grey matrix into a list. Each list's element is a virtual biopsy
    lista.grigi2 <- split(grey.matrix, rep(1:ncol(grey.matrix), each = nrow(grey.matrix)))
    
    # print the means, sd and grey list
    return(list("medie"=media, "devianze"= deviance, "lista.grigi"=lista.grigi2))
  }

  # ========================================================================================
  # allPopulationVirtualBiopsy
  # set an attribute of the class
  # ========================================================================================
  allPopulationVirtualBiopsy<-function( nx=2,ny=2,nz=0, ROIName4Normalization, grayTuniningValue, normalization=TRUE) {
    
    # if requested, get the higher value in order to "normalize" the greylevel
    if( normalization == TRUE) {
      UpperBoundDiNormalizzazione <- grayTuniningValue
#      UpperBoundDiNormalizzazione <- ROIStats(ROIName4Normalization)$total$max;
    }
    else {UpperBoundDiNormalizzazione=1;}
    
    # get the number of patient to treat
    total<-length(names(dataStructure));
    
    medie<-list();
    devianze<-list();
    
    for(i in names(dataStructure) )  { 
      
      # do you need a normalization ROI?
      if( normalization == TRUE)
        { piscioPaziente4Tuning<-mean(dataStructure[[ i ]]$voxelCubes[[ROIName4Normalization]][which(dataStructure[[ i ]]$voxelCubes[[ROIName4Normalization]]!=0)])}
      else {piscioPaziente4Tuning=1;}
      # get the biopsy for the given patient
      a<-virtualBiopsy(   (dataStructure[[ i ]]$voxelCubes$GTV)*(UpperBoundDiNormalizzazione/piscioPaziente4Tuning)    ,nx,ny,nz)
      # add the result in terms of mean and std.deviation
      medie[[i]]<-a$medie
      devianze[[i]]<-a$devianze  
    }
    return( list( "medie"=medie, "devianze"=devianze  ) )    
  } 
  # ========================================================================================
  # setAttribute
  # set an attribute of the class
  # ========================================================================================
  setAttribute<-function(attribute, value) {
    if(attribute=="verbose") {
      if(!is.list(value)) return;
      for(i in names(value)) {
        attributeList$verbose[[ i ]] <- value[[i]]
      }
      logObj$setOutput( attributeList$verbose )
      return;  
    }
    #attributeList[[ attribute ]]<<-value     
  } 
  # ========================================================================================
  # setAttribute
  # get an attribute of the class
  # ========================================================================================
  getAttribute<-function(attribute) {
    if(attribute == "dataStorage") return(dataStructure)
    if(attribute == "results") return(arrayAR)
  }
  # ========================================================================================
  # constructor
  # ========================================================================================
  constructor<-function() {
    dataStructure<<-list()
    attributeList<<-list()
    arrayAR<<-list()
    # initialize LOG attributes
    attributeList$verbose<<-list("lv1"=TRUE,"lv2"=TRUE,"lv3"=FALSE,"onScreen"=TRUE,"onFile"=FALSE)    
    # define LOG object handler
    logObj<<-logHandler()
    # setting log defaults
    logObj$setOutput( list("onScreen" = attributeList$verbose$onScreen,   "onFile" = attributeList$verbose$onFile )  )
  }
  constructor()
  return(list(
      openTreeMultiROIs=openTreeMultiROIs,
      setAttribute=setAttribute,
      getAttribute=getAttribute,
      execAlgorithm=execAlgorithm,
      plotResults=plotResults,
      
      ROIStats=ROIStats))
}

# example -- mmButo stile Toronto
#rm(obj.positive)
#obj.positive<-RAD.mmButo()
#obj.positive$openTreeMultiROIs("/progetti/immagini/CONTOURED/Positive/easy", structureList=c("GTV","Retto"))
#b<-obj.positive$getAttribute("dataStorage")
#image(b[[1]]$voxelCubes$GTV[,,10])
#obj.positive$execAlgorithm(algorithm="KDF",ROIName="GTV");

# a<-obj$openTreeMultiROIs("/progetti/immagini/CONTOURED/Positive/easy", structureList=c("GTV","Retto"))

# example -- ossa stile Sarhos
#obj<-geoLet();
#obj$openDICOMFolder("/progetti/immagini/SarhoshTest")
#a<-obj$getAttribute("dataStorage")
#image(a$dose[[1]][,,10], col = gray.colors(2^8))
#ROIName<-"Ossa, NOS"
#ds<-obj$getAttribute("dataStorage")
#ROIPointList<-obj$getROIPointList(ROIName)
#SS<-names(ds$img);
#SeriesInstanceUID<-"1.2.840.113619.2.278.3.17485314.39.1392360071.244"
#Structure<-ROIName
#dataStorage<-a
#result<-RAD.NewMultiPIPOblique(dataStorage = ds, Structure = ROIName, SeriesInstanceUID = SS)
## dyn.load("./src/PointInPolygon.so"); 


