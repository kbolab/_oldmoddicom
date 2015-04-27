#' Calculates image voxels internal respect a given ROI
#' @description This function can be used to calculate internal voxels, for each ROIs, in a \code{dataStorage} data structure taken from a \code{geoLet} object
#' @param dataStorage is the structure returned from a \code{obj$getAttribute("dataStructure")} where \code{obj} is an instance of a \code{geoLet}
#' @param Structure is the ROI name od the ROI that should be extracted
#' @param SeriesInstanceUID Is the interested series Instance UID.
#' @export
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
# 
# > arrayAssociationROIandSlice[0:70]
# [1] -10000      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8
# [22]      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8
# [43]      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      9 -10000
# [64]      9      9      9      9      9      9      9
#     
  
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
#' @export
RAD.mmButo<-function() {
  
  dataStructure<-list()
  attributeList<-list()
  logObj<-list()
  arrayAR<-list()           # arrayAlgorithmResult
  
  openTreeMultiROIs<-function(Path, structureList) {
    cubeVoxelList<-list()
    listaFolders<-list.dirs( Path )
    
    # a bit of froceries
    #    pb <- txtProgressBar(min = 0, max = length(listaFolders)-1, style = 3)
    counter<-1
    for( folderName in listaFolders[2:length(listaFolders)] ) {
      #      setTxtProgressBar(pb, counter)
      # Instantiate the object
      obj<-geoLet()
      #obj$setAttribute(attribute="verbose",value=FALSE)  
#      print( folderName )
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
          
          #result<-MultiPIPOblique(dataStorage = ds, Structure = ROIName, SeriesInstanceUID = SS, output=c("masked.images","DICOMInformationMatrix","resampled.array"))
          #result<-RAD.MultiPIPOblique(dataStorage = ds, Structure = ROIName, SeriesInstanceUID = SS)
          result<-RAD.NewMultiPIPOblique(dataStorage = ds, Structure = ROIName, SeriesInstanceUID = SS)
          
          cubeVoxelList[[ folderName ]][["voxelCubes"]][[ ROIName ]]<-result$masked.images    
          #cubeVoxelList[[ folderName ]][["info"]]<-result$DICOMInformationMatrix
          cubeVoxelList[[ folderName ]][["info"]]<-result$DOM
          
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
    #    close(pb) 
    dataStructure<<-cubeVoxelList      
  }


  execAlgorithm<-function( algorithm, ROIName ) {
    
    XmaxVal<-c() 
    YmaxVal<-c() 
    array4VoxelCube<<-list();

    # Kernel Density Function
    if( algorithm == "KDF" ) {      
      arrayAR[["KDF"]]<<-list();      
      # loop for each Series Instance UID
      for( SeriesInstanceUID in names(dataStructure)) {          
        XmaxVal <- c(XmaxVal,max(dataStructure[[SeriesInstanceUID]]$image.arr));      
        # get only the 'not zero' voxel value
        voxelCube <- dataStructure[[SeriesInstanceUID]]$voxelCubes[[ROIName]]
        array4VoxelCube[[SeriesInstanceUID]] <- voxelCube[which(  voxelCube  !=0 )];
      } 
      # NORMALIZATION (for upper bound only)
      for( SeriesInstanceUID in names(dataStructure)) {
       array4VoxelCube[[SeriesInstanceUID]]<-density(array4VoxelCube[[SeriesInstanceUID]]/max( XmaxVal ))
       YmaxVal <- c(YmaxVal,max(array4VoxelCube[[SeriesInstanceUID]]$y));      
      }
      # write the results in the array
      arrayAR$KDF$details<<-array4VoxelCube
      arrayAR$KDF$summary<<-list();
      arrayAR$KDF$summary$XmaxVal<<-max( XmaxVal )
      arrayAR$KDF$summary$YmaxVal<<-max( YmaxVal )
      return();
    }
    # Not yet implemented NMI error
    logObj$sendLog(message = "Not yet implemented", NMI = TRUE);
  }

  plotResults<-function( algorithm="all" ) {        
    if( algorithm == "all" | algorithm == "KDF")  {
      ct<-1;
      for(pathName in arrayAR$KDF$details ) {        
        plot( arrayAR$KDF$details[[1]], ylim=c(0,arrayAR$KDF$summary$YmaxVal) ) 
        lines( arrayAR$KDF$details[[2]] ) 
        lines( arrayAR$KDF$details[[3]] )        
      }      
    }
  } 

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
  getAttribute<-function(attribute) {
    if(attribute == "dataStorage") return(dataStructure)
    if(attribute == "algorithmsResult") return(arrayAR)
  }
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
  return(list(openTreeMultiROIs=openTreeMultiROIs,setAttribute=setAttribute,getAttribute=getAttribute,execAlgorithm=execAlgorithm,plotResults=plotResults))
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



