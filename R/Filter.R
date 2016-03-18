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
#' @param ROINameForNormalization the name of the ROI used to normalize the signal
#' @param valueForNormalization the highest value used to normalize the 'ROINameForNormalization' voxel values
#' @param ROIName the name of the ROI to extract
#' @param filter.pipeline is a \code{list} where is indicated the pipeline of filtering to be applied (in sequence)
#' @param collection is the collection: the default is '\code{default}'
#' @param cropResult is a boolean (\code{TRUE} or \code{FALSE}) which indicates if the result should be cropped or not, in order to save memory. Default is \code{TRUE}
#' @param scaleFactor can be 'voxel' or 'space'. If 'voxel' (the default) it consider the sam sigma values (if specified) for all the geoLet object, if 'space' it normalize the sigma according to the different pixelSpacing.
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
#' a<-FIL.applyFilterToStudy(obj.mmButo = obj,ROIName = "GTV",ROINameForNormalization='Urina', valueForNormalization=10000,filter.pipeline = filterPipeline )
#'
#' }#' 
FIL.applyFilterToStudy<-function(obj.mmButo, ROINameForNormalization=NA, valueForNormalization = NA, ROIName, filter.pipeline ,collection="default",cropResult=TRUE, scaleFactor="voxel") {
  objS<-services();
  FUNMap<-list();
  if(scaleFactor != "voxel" & scaleFactor != "space") stop("\n Error: 'scaleFactor' can be only 'voxel' or 'space'.");
  
  # prendi i pixelSpacing
  pixelSpacingArr<-obj.mmButo$getAttribute(attributeName = "list_PixelSpacing")
  
  # prendi la ROI
  ROIVoxelData<-obj.mmButo$getROIVoxel(ROIName = ROIName)
  
#   # ERRORE se non è indicata una ROI per la normalizzazione
#   if(is.na(ROINameForNormalization)) {
#     stop("ERROR: 'ROINameforNormalization' must be specified");
#   }  

  if(!is.na(ROINameForNormalization)) {
    # prendi le ROI per la normalizzazione e calcolane le statistiche
    ROIVoxelDataForNormalization<-obj.mmButo$getROIVoxel(ROIName = ROINameForNormalization)
    statsForNormalization<-obj.mmButo$getROIVoxelStats(ROIVoxelList = ROIVoxelDataForNormalization)

    # se non è stato passato un valore di normalizzazione prendi il massimo delle medie della ROI
    # utilizzata per normalizzare
    if(is.na(valueForNormalization)) {
      arrayNormalizationROIVoxelValue<-c();
      for ( patient in names(statsForNormalization$details) ) {
        arrayNormalizationROIVoxelValue<-c(arrayNormalizationROIVoxelValue,statsForNormalization$details[[ patient ]]$mean)
      }  
      valueForNormalization<-max(arrayNormalizationROIVoxelValue)
    }
    # calcola l'array dei fattori moltiplicativi
    multiplier<-list();
    for ( patient in names(statsForNormalization$details) ) {
      multiplier[[ patient ]] <- valueForNormalization / statsForNormalization$details[[ patient ]]$mean 
    }  
  } 
  else {
    # calcola l'array dei fattori moltiplicativi nel caso in cui non si voglia normalizzazione 
    # (piazzali tutti a uno)
    multiplier<-list();
    for ( patient in names(ROIVoxelData)) {
      multiplier[[ patient ]] <- 1
    }      
  } 

  # prendi la lista di oggetti geoLet
  list_geoLet<-obj.mmButo$getAttribute("list_geoLet");
  
  # inizia a ciclare su tutti i pazienti
  for(patient in names(ROIVoxelData)) {
    
    if((TRUE %in% is.na(ROIVoxelData[[patient]])) == FALSE  ) {
      print( paste("FUN - Now filtering:",patient),collapse='' )
      # prendi i voxelData ed estendili
      voxelData.ready<-ROIVoxelData[[patient]]
      voxelData.ready.espanso <- obj.mmButo$mmButoLittleCube.expand(   voxelData.ready   ) 
      voxelDataDaRiapplicare <- voxelData.ready.espanso;
      
      # e crea la relativa maschera
      #voxelData.ready.espanso[voxelData.ready.espanso!=0]<-1
      voxelData.ready.espanso[!is.na(voxelData.ready.espanso)]<-1
      
      # prendi l'original voxelCute non tagliato dalla ROI
      # è su questo che dovrò applicare il filtro
      originalMR <- list_geoLet[[collection]][[patient]]$getImageVoxelCube()
      # ora moltiplica l'MR originale per il valore di normalizzazione
      originalMR <- originalMR * multiplier[[ patient ]] 
      
      for(i in seq(1,length(filter.pipeline))) {
        print( paste("     => applying",filter.pipeline[[i]]$kernel.type),collapse='' )
        if (scaleFactor == "space") {
          normalizedSigma<-sqrt((filter.pipeline[[i]]$sigma^2)/(pixelSpacingArr[[patient]][1]*pixelSpacingArr[[patient]][2]))
        } else {
          normalizedSigma<-filter.pipeline[[i]]$sigma
        }
        originalMR<-FIL.applyFilter( originalMR, 
                                     kernel.type = filter.pipeline[[i]]$kernel.type,
                                     sigma = normalizedSigma)
                                     #sigma = filter.pipeline[[i]]$sigma)
      }
      # l'output deve essere mascherato
      u<-voxelData.ready.espanso
      
      #u[  u[,,]!=0 ]<-1
      u[which(!is.na(u),arr.ind = TRUE)]<-1
      u[which(is.na(u),arr.ind = TRUE)]<-0
      
      #FUNMap[[patient]]<-originalMR * voxelData.ready.espanso;
      # Applica la maschera

      FUNMap[[patient]]<-originalMR * u;
      FUNMap[[patient]][which(is.na(voxelData.ready.espanso),arr.ind = TRUE)]<-NA
      
      # se richiesto, CROPPA!
      if ( cropResult == TRUE )  FUNMap[[patient]]<-objS$cropCube( FUNMap[[patient]] )   
    }
    else {
      FUNMap[[patient]]<-NA
    }
  }
  return(FUNMap)
}
