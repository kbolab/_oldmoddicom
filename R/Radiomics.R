#' Calculates image voxels internal respect a given ROI
#' @description This function can be used to calculate internal voxels, for each ROIs, in a \code{dataStorage} data structure taken from a \code{geoLet} object
#' @param dataStorage is the structure returned from a \code{obj$getAttribute("dataStructure")} where \code{obj} is an instance of a \code{geoLet}
#' @param Structure is the ROI name od the ROI that should be extracted
#' @param SeriesInstanceUID Is the interested series Instance UID.
#' @export
#' @return A list with unknown meaning
RAD.MultiPIPOblique<-function(dataStorage, Structure, SeriesInstanceUID) {
  objService<-services()
  # initialize the image array with the right dimension
  image.arr<-array(data=-1, dim = c(dataStorage$info[[SeriesInstanceUID]][[1]]$Rows, 
                                    dataStorage$info[[SeriesInstanceUID]][[1]]$Columns,
                                    length(dataStorage$img[[SeriesInstanceUID]])))
  # index values listed as characters. 
  index<-as.character(sort(as.numeric(names(dataStorage$img[[SeriesInstanceUID]]))))  
  # creates the array of images ordered according referenceInstanceUID
  nn<-0
  # creates empty array of DICOM orientation matrices
  DOM<-c()
  for (n in index) {
    # fills the vector of DICOM orientation matrices
    DOM<-c(DOM, dataStorage$info[[SeriesInstanceUID]][[n]]$orientationMatrix[c(1:3,5:7,13:15)])
    nn<-nn+1
    image.arr[,,nn]<-dataStorage$img[[SeriesInstanceUID]][[n]]    
  }  
  
  # fills the vectors of X and Y coordinates, offset and empty or full slices
  TotalX<- -10000
  TotalY<- -10000
  Offset<-0
  FullZ <-c()
  OriginX<- -10000
  OriginY<- -10000
  for (n in index) {
    if (length(which(dataStorage$info[[SeriesInstanceUID]][[n]]$ROIList[,1]==Structure))==0) {
      FullZ<-c(FullZ, 0)      
      tempOffset<-1
      TotalX<-c(TotalX, OriginX)
      TotalY<-c(TotalY, OriginY)
      Offset<-c(Offset, Offset[length(Offset)] + 1)
    } else {
      # fills full slices
      FullZ<-c(FullZ, 1)
      for (m in which(dataStorage$info[[SeriesInstanceUID]][[n]]$ROIList[,1]==Structure)) {
        # find the slice and gets the key for accessing at coordinates vectors
        key<-dataStorage$info[[SeriesInstanceUID]][[n]]$ROIList[m,2]
        TotalX<-c(TotalX, dataStorage$structures[[Structure]][[key]][[1]][,1])
        TotalY<-c(TotalY, dataStorage$structures[[Structure]][[key]][[1]][,2])
        TotalX<-c(TotalX, OriginX)
        TotalY<-c(TotalY, OriginY)
        # set the tempOffset
        tempOffset<-length(dataStorage$structures[[Structure]][[key]][[1]][,1]) + 1
        Offset<-c(Offset, Offset[length(Offset)] + tempOffset)
      }      
    }       
  }
  
  final.array<-RAD.MultiPointInPolyObl(DICOMOrientationVector = DOM, totalX = TotalX, totalY = TotalY, NumSlices = length(FullZ),  
                                   Offset = Offset,  nX = as.numeric(dataStorage$info[[SeriesInstanceUID]][[1]]$Columns), nY = as.numeric(dataStorage$info[[SeriesInstanceUID]][[1]]$Rows), FullZ = FullZ)
  final.array<-array(data = final.array, dim = c(dataStorage$info[[SeriesInstanceUID]][[1]]$Columns, dataStorage$info[[SeriesInstanceUID]][[1]]$Rows, length(FullZ)))
  
  for ( i in seq(1,dim(image.arr)[3] )) {
    image.arr[,,i]<-objService$rotateMatrix(image.arr[,,i])
    final.array[,,i]<-t(objService$rotateMatrix(final.array[,,i],rotations=3))
  }
  #image(rotateMatrix(image.arr[,,2]) * t(rotateMatrix(final.array[,,2],rotations=3)),col = grey(seq(0, 1, length = 256)))
  
  return(list(TotalX=TotalX, TotalY=TotalY, FullZ=FullZ, Offset=Offset, 
              DOM=array(DOM, dim = c(3,3,length(index))), final.array=final.array, masked.images=final.array*image.arr))
}
#' Wrapper for C function
#' @useDynLib moddicom
RAD.MultiPointInPolyObl<-function(DICOMOrientationVector, totalX, totalY, NumSlices, Offset, FullZ, nX, nY) {
  objService<-services()
#  dyn.load(objService$SV.LoadAccordingOSType(library.name = "PointInPolygon"))    
  # creates the PIPvector
  PIPvector<-rep.int(x = 0, times = nX * nY * NumSlices)  
  result<-.C("MultiPIPObl", as.double(totalX), as.double(totalY), as.integer(nX), as.integer(nY), 
             as.integer(NumSlices), as.integer(Offset), as.integer(PIPvector), as.integer(FullZ), as.double(DICOMOrientationVector))  
  return(result[[7]])
}
#' Load and handle a tree
#' @description This function can be used to load, in one shot, many DICOM studies referring different patients. Each study is stored in a separate element of a list and basics analysis can be done. Even if it is quite deprecated, it can still be useful for 'quick & dirty' image analysis.
#' @export
RAD.mmButo<-function() {
  
  dataStructure<-list()
  attributeList<-list()
  
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
      print( folderName )
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
          result<-RAD.MultiPIPOblique(dataStorage = ds, Structure = ROIName, SeriesInstanceUID = SS)
          
          cubeVoxelList[[ folderName ]][["voxelCubes"]][[ ROIName ]]<-result$masked.images    
          #cubeVoxelList[[ folderName ]][["info"]]<-result$DICOMInformationMatrix
          cubeVoxelList[[ folderName ]][["info"]]<-result$DOM
          
          cubeVoxelList[[ folderName ]][["ROIPointList"]][[ ROIName ]]<-ROIPointList
          
          # -im
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
    return( cubeVoxelList );
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
  constructor<-function() {
    dataStructure<<-list()
    attributeList<<-list()
    attributeList$verbose<<-list("lv1"=TRUE,"lv2"=TRUE,"lv3"=FALSE,"onScreen"=TRUE,"onFile"=FALSE)    
  }
  constructor()
  return(list(openTreeMultiROIs=openTreeMultiROIs,setAttribute=setAttribute))
}
# example
# obj<-RAD.mmButo()
# a<-obj$openTreeMultiROIs("/progetti/immagini/CONTOURED/Positive/easy", structureList=c("GTV","Retto"))
