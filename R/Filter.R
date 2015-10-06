# ================================================================
#' Apply a filter to a 2d or a 3d array
#' 
#' @description  It applies filter to a 2d or 3d array. '...' is banned, so shit cannot flies ( ;) ).
#' @param arr2App is the array containing the image that should be filtered
#' @param kernel.type is the kernel you want to use: \code{gaussian}, \code{laplacian}, \code{LoG},\code{shrpen}, \code{emboss}, \code{sobel}
#' @param sigma is the \code{sigma} value for gaussian filter
#' @return the filteres array (same geometry of the one in input)
#' @examples \dontrun{
#' # Create an instante of new.mmButo and load some cases
#' a<-array( 0, dim=c(10,10,10))
#' a[5,5,5]<-177; a[5,6,5]<-198; a[6,6,5]<-45; a[6,5,5]<-12;
#' b<-FIL.applyFilter( a , kernel.type="gaussian", sigma=1.4)
#' }#' 
#' @export
#' @import spatialfil
FIL.applyFilter<-function( arr2App, kernel.type, sigma = 1.4 ) {
  if( is.null(sigma) ) sigma = 1.4;
  kern2Apply<-convKernel(  sigma = sigma, k = kernel.type);
  ret <- applyFilter(x =arr2App, kernel = kern2Apply);
}

#' Apply a filter to mmButo object
#' 
#' @description  It applies filter to an mButo object to simplify RADIOMICS issues
#' @param obj.mmButo is the \code{new.mmButo} object
#' @param ROINameForNormalization if specified indicates which ROI should be used for signal normalization
#' @param ROIName is the name of the ROI that should be extracted and filtered
#' @param filter.pipeline is a \code{list} where is indicated the pipeline of filtering to be applied (in sequence)
#' @param collection is the collection: the default is '\code{default}'
#' @param cropResult is a boolean (\code{TRUE} or \code{FALSE}) which indicates if the result should be cropped or not, in order to save memory. Default is \code{TRUE}
#' @return a list containing the filtered images (cropped or not)
#' @export
#' @examples \dontrun{
#' # Create an instante of new.mmButo and load some cases
#' obj<-new.mmButo()
#' obj$loadCollection("/progetti/immagini/urinaEasy")
#' 
#' GTV<-obj$getROIVoxel(ROIName = "GTV")
#' Urina<-obj$getROIVoxel(ROIName = "Urina")
#' GTV.C<-obj$getCorrectedROIVoxel(inputROIVoxel = GTV,correctionROIVoxel = Urina,typeOfCorrection = "log")
#'
#' # build the filtering pipeline
#' filterPipeline<- list()
#' filterPipeline[[1]]<-list("kernel.type"="gaussian")
#' filterPipeline[[2]]<-list("kernel.type"="emboss", "sigma"=1.5)
#'
#' # filter the images (normalizing the GTV signal with Urin) cropping the result
#' a<-FIL.applyFilterToStudy(obj.mmButo = obj,ROIVoxelData = GTV.C, filter.pipeline = filterPipeline )
#'
#' }#' 
FIL.applyFilterToStudy<-function(obj.mmButo, ROIVoxelData, filter.pipeline ,collection="default",cropResult=TRUE) {
  objS<-services();
  FUNMap<-list();

  # prendi la lista di oggetti geoLet
  list_geoLet<-obj.mmButo$getAttribute("list_geoLet");
  for(patient in names(ROIVoxelData)) {
    
    if((TRUE %in% is.na(ROIVoxelData[[patient]])) == FALSE  ) {
      print( paste("FUN - Now filtering:",patient),collapse='' )
      
      # prendi i voxelData
      voxelData.ready<-ROIVoxelData[[patient]]
      # ora sono pronto per estendere
      voxelData.ready.espanso <- obj.mmButo$mmButoLittleCube.expand(   voxelData.ready   ) 
      # e creare la relativa maschera
      voxelData.ready.espanso[voxelData.ready.espanso!=0]<-1
      
      # prendi l'original voxelCute non tagliato dalla ROI
      # è su questo ceh dovrò applicare il filtro
      originalMR<-list_geoLet[[collection]][[patient]]$getImageVoxelCube()
      
      for(i in seq(1,length(filter.pipeline))) {
        print( paste("     => applying",filter.pipeline[[i]]$kernel.type),collapse='' )
        
        originalMR<-FIL.applyFilter( originalMR, 
                                     kernel.type = filter.pipeline[[i]]$kernel.type,
                                     sigma = filter.pipeline[[i]]$sigma)
      }
      # l'output deve essere mascherato
      FUNMap[[patient]]<-originalMR * voxelData.ready.espanso;
      
      # se richiesto, CROPPA!
      if ( cropResult == TRUE )  FUNMap[[patient]]<-objS$cropCube( FUNMap[[patient]] )         
    }
    else {
      FUNMap[[patient]]<-NA
    }
  }
  return(FUNMap)
}



# 
# 
# # ================================================================
# #' Apply a filter to a 2d or a 3d array
# #' 
# #' @description  It applies filter to a 2d or 3d array. '...' is banned, so shit cannot flies ( ;) ).
# #' @param arr2App is the array containing the image that should be filtered
# #' @param kernel.type is the kernel you want to use: \code{gaussian}, \code{laplacian} or \code{unsharp}
# #' @param sigma is the \code{sigma} value for gaussian filter
# #' @param nx is the \code{nx} value for gaussian filter
# #' @param ny is the \code{ny} value for gaussian filter
# #' @param alpha is the \code{alpha} value for laplacian and unsharp filters
# #' @return the filteres array (same geometry of the one in input)
# #' @examples \dontrun{
# #' # Create an instante of new.mmButo and load some cases
# #' a<-array( 0, dim=c(10,10,10))
# #' a[5,5,5]<-177; a[5,6,5]<-198; a[6,6,5]<-45; a[6,5,5]<-12;
# #' b<-FIL.applyFilter( a , kernel.type="unsharp", alpha=1)
# #' }#' 
# old.FIL.applyFilter<-function( arr2App, kernel.type , sigma=NA , alpha=NA, nx=NA, ny=NA ) {
#   if(kernel.type=="gaussian" && (is.na(nx) || is.na(ny) || is.na(sigma))  )  stop("for gaussian filter nx,ny,sigma have to be specified (i.e.: 15,15,1)")
#   if(kernel.type=="LoG" && (is.na(nx) || is.na(ny) || is.na(sigma))  )  stop("for LoG filter nx,ny,sigma have to be specified (i.e.: 15,15,1)")
#   if(kernel.type=="laplacian" && (is.na(alpha) )  )  stop("for laplacian filter alpha has to be specified (i.e.: 0)")
#   if(kernel.type=="unsharp" && (is.na(alpha) )  )  stop("for unsharp filter alpha has to be specified (i.e.: 0)")  
#   # se è in 2D passalo in 3D (una sola slice in 3D)  così la restante
#   # trattazione è unica per entrambi i casi
#   dimensioni<-3;
#   if(length(dim(arr2App))==2) {
#     arr2App<-array(arr2App,dim = c(dim(arr2App)[1],dim(arr2App)[2],1))
#     dimensioni<-2;
#   }
#   # ora cicla per ogni slice
#   for( z in seq(1,dim(arr2App)[3])  ) {
#     
#     if(kernel.type=="gaussian") newIm<-kernel2dsmooth( arr2App[,,z] , kernel.type="gauss", nx=nx, ny=ny, sigma=sigma)
#     if(kernel.type=="LoG") newIm<-kernel2dsmooth( arr2App[,,z] , kernel.type="LoG", nx=nx, ny=ny, sigma=sigma)
#     if(kernel.type=="laplacian")  newIm<-kernel2dsmooth( arr2App[,,z], K=kernel2dmeitsjer(type = "laplacian",alpha=alpha))
#     if(kernel.type=="unsharp") newIm<-kernel2dsmooth( arr2App[,,z], kernel.type="unsharp", alpha=alpha)
#     arr2App[,,z]<-newIm
#   }
#   if(dimensioni==2) return(arr2App[,,1])
#   if(dimensioni==3) return(arr2App)
# }
# 
# #' Apply a filter to mmButo object
# #' 
# #' @description  It applies filter to an mButo object to simplify RADIOMICS issues
# #' @param obj.mmButo is the \code{new.mmButo} object
# #' @param ROINameForNormalization if specified indicates which ROI should be used for signal normalization
# #' @param ROIName is the name of the ROI that should be extracted and filtered
# #' @param filter.pipeline is a \code{list} where is indicated the pipeline of filtering to be applied (in sequence)
# #' @param collection is the collection: the default is '\code{default}'
# #' @param cropResult is a boolean (\code{TRUE} or \code{FALSE}) which indicates if the result should be cropped or not, in order to save memory. Default is \code{TRUE}
# #' @return a list containing the filtered images (cropped or not)
# #' @examples \dontrun{
# #' # Create an instante of new.mmButo and load some cases
# #' obj<-new.mmButo()
# #' obj$loadCollection("/progetti/immagini/urinaEasy")
# #' 
# #' # build the filtering pipeline
# #' filterPipeline<- list()
# #' filterPipeline[[1]]<-list("kernel.type"="gaussian", "nx"=10, "ny"=10, "sigma"=1)
# #' filterPipeline[[2]]<-list("kernel.type"="laplacian", "alpha"=1.5)
# #'
# #' # filter the images (normalizing the GTV signal with Urin) cropping the result
# #' a<-FIL.applyFilterToStudy(obj.mmButo = obj,ROIName = "GTV",ROINameForNormalization='Urina', filter.pipeline = filterPipeline )
# #'
# #' }#' 
# old.FIL.applyFilterToStudy<-function(obj.mmButo, ROINameForNormalization=NA, ROIName, filter.pipeline ,collection="default",cropResult=TRUE) {
#   objS<-services();
#   FUNMap<-list();
#   # prendi la ROI
#   ROIVoxelData<-obj.mmButo$getROIVoxel(ROIName = ROIName)
#   # correggi se è il caso di farlo
#   if(!is.na(ROINameForNormalization)) {
#     if (!is.list(ROINameForNormalization))
#       ROIForCorrection <- obj.mmButo$getROIVoxel(ROIName = ROINameForNormalization)
#     else ROIForCorrection <- ROINameForNormalization
#     ROIVoxelData<-obj.mmButo$getCorrectedROIVoxel(inputROIVoxel = ROIVoxelData, correctionROIVoxel = ROIForCorrection)
#   }  
#   
#   # prendi la lista di oggetti geoLet
#   list_geoLet<-obj.mmButo$getAttribute("list_geoLet");
#   for(patient in names(ROIVoxelData)) {
#     
#     if((TRUE %in% is.na(ROIVoxelData[[patient]])) == FALSE  ) {
#       print( paste("FUN - Now filtering:",patient),collapse='' )
#       
#       # prendi i voxelData
#       voxelData.ready<-ROIVoxelData[[patient]]
#       # ora sono pronto per estendere
#       voxelData.ready.espanso <- obj.mmButo$mmButoLittleCube.expand(   voxelData.ready   ) 
#       # e creare la relativa maschera
#       voxelData.ready.espanso[voxelData.ready.espanso!=0]<-1
#       
#       # prendi l'original voxelCute non tagliato dalla ROI
#       # è su questo ceh dovrò applicare il filtro
#       originalMR<-list_geoLet[[collection]][[patient]]$getImageVoxelCube()
# 
#       for(i in seq(1,length(filter.pipeline))) {
#         print( paste("     => applying",filter.pipeline[[i]]$kernel.type),collapse='' )
#         originalMR<-FIL.applyFilter( originalMR, 
#                                      kernel.type = filter.pipeline[[i]]$kernel.type, 
#                                      sigma = filter.pipeline[[i]]$sigma, 
#                                      alpha = filter.pipeline[[i]]$alpha, 
#                                      nx = filter.pipeline[[i]]$nx, 
#                                      ny = filter.pipeline[[i]]$ny )
#       }
#       # l'output deve essere mascherato
#       FUNMap[[patient]]<-originalMR * voxelData.ready.espanso;
#       
#       # se richiesto, CROPPA!
#       if ( cropResult == TRUE )  FUNMap[[patient]]<-objS$cropCube( FUNMap[[patient]] )         
#     }
#     else {
#       FUNMap[[patient]]<-NA
#     }
#   }
#   return(FUNMap)
# }