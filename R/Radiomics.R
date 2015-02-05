#' Calculates image voxels internal respect a given ROI
#' @description This function can be used to calculate internal voxels, for each ROIs, in a \code{dataStorage} data structure taken from a \code{geoLet} object
#' @param dataStorage is the structure returned from a \code{obj$getAttribute("dataStructure")} where \code{obj} is an instance of a \code{geoLet}
#' @param Structure is the ROI name od the ROI that should be extracted
#' @param SeriesInstanceUID Is the interested series Instance UID.
#' @export
#' @return A list with unknown meaning
RAD.MultiPIPOblique<-function(dataStorage, Structure, SeriesInstanceUID) {
  
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
  
  final.array<-MultiPointInPolyObl(DICOMOrientationVector = DOM, totalX = TotalX, totalY = TotalY, NumSlices = length(FullZ),  
                                   Offset = Offset,  nX = as.numeric(dataStorage$info[[SeriesInstanceUID]][[1]]$Columns), nY = as.numeric(dataStorage$info[[SeriesInstanceUID]][[1]]$Rows), FullZ = FullZ)
  final.array<-array(data = final.array, dim = c(dataStorage$info[[SeriesInstanceUID]][[1]]$Columns, dataStorage$info[[SeriesInstanceUID]][[1]]$Rows, length(FullZ)))
  
  for ( i in seq(1,dim(image.arr)[3] )) {
    image.arr[,,i]<-rotateMatrix(image.arr[,,i])
    final.array[,,i]<-t(rotateMatrix(final.array[,,i],rotations=3))
  }
  #image(rotateMatrix(image.arr[,,2]) * t(rotateMatrix(final.array[,,2],rotations=3)),col = grey(seq(0, 1, length = 256)))
  
  return(list(TotalX=TotalX, TotalY=TotalY, FullZ=FullZ, Offset=Offset, 
              DOM=array(DOM, dim = c(3,3,length(index))), final.array=final.array, masked.images=final.array*image.arr))
}