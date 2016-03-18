#' class for viewing
#' 
#' @description  allows to view what's happening
#' @useDynLib moddicom
#' @import Rvcg rgl 
#' @export
viewer<-function() {
  
  isosurface<-function( matrice , lower=lower , metodo= "wired", percent = 0.5, pixelSpacing, ... ) {
    parametri<-list(...)
    x<-seq(0, dim(matrice)[1] )
    y<-seq(0, dim(matrice)[2] )
    z<-seq(0, dim(matrice)[3] )
    g<-expand.grid(x,y,z)
    v<-as.array( matrice )
    storage.mode(v) <- "integer"
    mesh <- vcgIsosurface( v, threshold = lower, spacing = pixelSpacing)
    decimface <- vcgQEdecim( mesh, percent = percent)
    if(metodo == "noPlot" ) {
      return(mesh)
    }
    if(metodo == "wired" ) {
      if (is.null (parametri$col ) ) parametri$col = 8
  wire3d( mesh , col = parametri$col )
    }
    if(metodo == "smothed" ) {
      if (is.null (parametri$mu ) ) parametri$mu = -0.53 
      if (is.null (parametri$lambda ) ) parametri$lambda = 0.5
      if (is.null (parametri$col ) ) parametri$col = 8
      #shade3d(vcgSmooth(mesh = mesh, lambda = parametri$lambda, mu = parametri$mu) , col = parametri$col)
      
      shade3d(vcgSmooth(mesh = mesh, lambda = parametri$lambda, mu = parametri$mu) , color = "grey")   
      return();
      return(list(mesh=mesh, lambda=parametri$lambda, mu = parametri$mu))
    } 
  }
  
  becca<-function( obj.geoLet, SeriesInstanceUID = 1, ROIName, pixelSpacing  ) {
    SeriesInstanceUID<-1
    ROIName<-"Ossa, NOS"
    pixelSpacing<-c(0.76,0.76,2.5)
    
    matrice<-obj.geoLet$getROIVoxels("Ossa, NOS", 1)$masked.images
    x<-seq(0, dim(matrice)[1] )
    y<-seq(0, dim(matrice)[2] )
    z<-seq(0, dim(matrice)[3] )
    g<-expand.grid(x,y,z)
    v<-as.array( matrice )
    storage.mode(v) <- "integer"
    mesh <- vcgIsosurface( v, lower = 1, spacing = pixelSpacing)
    mesh <- vcgQEdecim( mesh, percent = .5)
    
    box.x<-c(  min(mesh$vb[1,]  ),max(  mesh$vb[1,])  )
    box.y<-c(  min(mesh$vb[2,]  ),max(  mesh$vb[2,])  )
    box.z<-c(  min(mesh$vb[3,]  ),max(  mesh$vb[3,])  )
    dotCoords<-expand.grid(box.x,box.y,box.z)
    
    wire3d( mesh , col = "gray" )
    points3d(x = dotCoords[,1],y=dotCoords[,2],z=dotCoords[,3], size=8, col="red" )
    planes3d(.6, -.9, .3, 0, alpha=.5,col="red")
  } 

  plotROIs<-function(obj , arrayROINames = NA , sliceNumber = 1, ps.x = NA, ps.y = NA, ps.z = NA) {
    # controlli formali
    if( length(arrayROINames)>20 )  stop("Ahahahahahahah!!!! Ciccio, al massimo 20 ROI, forse meno.....");
    # definisce una palette di circa 20 color (bo'... son 20?)
    paletteColor<-c('#38D747','#F366FD','#C9360E','#ADCDE7','#391C1F','#F3BC6E','#058557','#5E601A','#FCACA9','#F7556B','#E7EDBE','#759D16','#551E50','#A76A19','#E994F2','#2DC3CB','#B3EE5A','#FB8F5B','#40B85C','#04C1A9')
    # verifica che ce ne sia almeno una
    if(length(arrayROINames)==0) return;
    align<-list();  cubeDim<-c();
    # frulla sull'array di ROINames per allinearle una per una
    for( ROIName in seq(1,length(arrayROINames))) {
      tmpAlign <- obj$getAlignedStructureAndVoxelCube(ROIName = arrayROINames[ROIName] , ps.x = ps.x, ps.y = ps.y,ps.z = ps.z)
      # prendi i punti ROI allineati
      align[[ arrayROINames[ ROIName] ]]<-tmpAlign$ROI
      # prendi il voxelCube (anche se riscrive non è importante, tanto è sempre lo stesso)
      vc<-tmpAlign$voxelCube
      # prendi le dimensioni del voxelCube (sovrascrive? No problem, vedi sopra)
      #cubeDim<-tmpAlign$cubeDim
      cubeDim<-dim(tmpAlign$voxelCube)
    }
    
    # plotta l'immagine  
    image(vc[,,sliceNumber],col = grey.colors(255))
    # prendi in numero di slice
    whichROILine<-sliceNumber;
    colorIncrement<-1
    
    # Cicla per plottare ogni ROI nell'array delle ROI
    for( ROIName in seq(1,length(arrayROINames))) {
      # scorri nella lista delle ROIPointList per cercare la fetta corretta
      for( internalROIPlanes in seq(1,length( align[[   arrayROINames[ ROIName]    ]] ))    ) {
        # 
        if(align[[   arrayROINames[ ROIName]    ]][[internalROIPlanes]][[1]][1,3] == whichROILine  ) {
          
          for(subROIOnTheSamePlane in seq(1, length(align[[   arrayROINames[ ROIName]    ]][[internalROIPlanes]])    )){
            ROI<-align[[   arrayROINames[ ROIName]    ]][[internalROIPlanes]][[subROIOnTheSamePlane]]
            
            ROI[,1]<-ROI[,1] / cubeDim[1]
            ROI[,2]<-1-ROI[,2] / cubeDim[2]
            points(x = ROI[,1],y = ROI[,2],type='l', lwd = 3, col = paletteColor[colorIncrement])  
          }
          colorIncrement<-colorIncrement+1;
        }
      }
    }
  }
  plot.DVH<-function(obj, whichDVH=1) {
    a<-obj$getAttribute(attribute = "DVHsFromFile");
    volume<-a[[as.character(whichDVH)]]$DVHData.volume
    dose<-a[[as.character(whichDVH)]]$DVHData.dose
    mainTitle<-paste( c("DVH of '",as.character(whichDVH),"'"),collapse='');
    plot(x=dose,y=volume,type='l', main=mainTitle)

  }
  
  return( list( isosurface = isosurface, 
                plotROIs = plotROIs,
                plot.DVH = plot.DVH)  )
  
}



checkDifferences<-function(obj_geoLet, ROIName="CTV1", plotIt = TRUE, newPixelSpacing = NA,
                           verbose=FALSE, forceReCalculus=FALSE, fastEngine = TRUE,
                           decimation=FALSE, decimation.percentage=0.8, 
                           smoothing=FALSE, smoothing.iterations = 10) {
  obj<-obj_geoLet;
  ds<-obj$getAttribute(attribute = "dataStorage")
  enci<-obj$calculateDVH(ROIName = ROIName, newPixelSpacing = newPixelSpacing,
                         verbose=verbose, forceReCalculus=forceReCalculus, fastEngine = fastEngine,
                         decimation=decimation, decimation.percentage=decimation.percentage, 
                         smoothing=smoothing, smoothing.iterations = smoothing.iterations)
  onci<-ds$info$DVHs[[1]]$DVHFromFile[[ROIName]]$DVHObj
  
  enci@volume<-enci@dvh[1,2]
  onci@volume<-onci@dvh[1,2]
  
  enci.cum<-enci
  onci.cum<-onci
  
  enci.diff<-DVH.cum.to.diff(enci)
  onci.diff<-DVH.cum.to.diff(onci)
  
  both.cum<-DVH.merge(receiver = enci.cum,addendum = onci.cum)
  both.diff<-DVH.merge(receiver = enci.diff,addendum = onci.diff)
  
  # error measures
  total.volume<-both.cum@dvh[1,2]
  volumeratio<-min(both.diff@volume[1]/both.diff@volume[2], both.diff@volume[2]/both.diff@volume[1])
  
  delta.cum.array<-((both.cum@dvh[,2]-both.cum@dvh[,3]))
  delta.diff.array<-((both.diff@dvh[,2]-both.diff@dvh[,3]))
  
#    delta.cum.array[which(is.na(delta.cum.array) | is.infinite(delta.cum.array),arr.ind = T)]<-0
#    delta.diff.array[which(is.na(delta.diff.array) | is.infinite(delta.diff.array),arr.ind = T)]<-0
  max.cum.tot<-max(abs(delta.cum.array))
  max.diff.tot<-max(abs(delta.diff.array))
  
  if(plotIt==TRUE) {
    par(mfrow=c(2,1))
    #cumulatives
    max.y.value <- max(both.cum@dvh[,2], both.cum@dvh[,3])
    par(mar=c(5,4,4,8)+.1)
    
    plot(x = both.cum@dvh[,1],y = both.cum@dvh[,2], type='l',axes=T, main=paste(c('cumulative DVHs, % error (',ROIName,')'),collapse = ''), ylab='Volume', xlab='Dose', col='darkred', ylim=c(0,max.y.value))

    lines(x = both.cum@dvh[,1],y = both.cum@dvh[,3], type='l',col='blue' )
    axis(2, ylim=c(0,max(both.cum@dvh[,3])),col="black",lwd=2)

    par(new=TRUE)
    plot(x = both.cum@dvh[,1], axes=F,y = delta.cum.array, type='l',col='darkgreen', ylim=c(min(delta.cum.array),max(delta.cum.array)), ylab='' , xlab='', lty=3, lwd=2)
    axis(4, ylim=c(min(delta.cum.array),max(delta.cum.array)),lwd=2,line=0)
    abline(h = 0,lty=3)
    mtext("abs error",side=4,line=3)
    
    #  differential
    max.y.value <- max(both.diff@dvh[,2], both.diff@dvh[,3])
    par(mar=c(5,4,4,8)+.1)
    plot(x = both.diff@dvh[,1],y = both.diff@dvh[,2], type='l',axes=T, main=paste(c('differential DVHs, % error (',ROIName,')')), ylab='Volume', xlab='Dose', col='darkred', ylim=c(0,max.y.value))
    lines(x = both.diff@dvh[,1],y = both.diff@dvh[,3], type='l',col='blue' )
    axis(2, ylim=c(0,max(both.diff@dvh[,3])),col="black",lwd=2)
    
    par(new=TRUE)
    plot(x = both.diff@dvh[,1], axes=F,y = delta.diff.array, type='l',col='darkgreen', ylim=c(min(delta.diff.array),max(delta.diff.array)), ylab='' , xlab='', lty=3, lwd=2)
    axis(4, ylim=c(min(delta.diff.array),max(delta.diff.array)),lwd=2,line=0)
    mtext("abs error",side=4,line=3)
    abline(h = 0,lty=3)
    
    }
  
  return( list( "volumeratio"= volumeratio,
                "delta.cum.array" = delta.cum.array,
                "delta.diff.array" = delta.diff.array,
                "max.cum.tot"= max.cum.tot,
                "max.diff.tot"= max.diff.tot,
                "total.volume"=total.volume
               )
            )
}



checkDifferencesForAllROIs<-function(obj_geoLet, ROINameArray, plotIt = TRUE, newPixelSpacing = NA,
                                     verbose=FALSE, forceReCalculus=FALSE, fastEngine = TRUE,
                                     decimation=FALSE, decimation.percentage=0.8, 
                                     smoothing=FALSE, smoothing.iterations = 10) {
  obj<-obj_geoLet

  e.list<-list();
  
  for( nomeROI in ROINameArray) {
    if(length(obj$getROIPointList(ROINumber = nomeROI))!=0) {
      a<-checkDifferences(obj_geoLet = obj,ROIName = nomeROI,plotIt = plotIt, newPixelSpacing = newPixelSpacing,
                          verbose=verbose, forceReCalculus=forceReCalculus, fastEngine = fastEngine,
                          decimation=decimation, decimation.percentage=decimation.percentage, 
                          smoothing=smoothing, smoothing.iterations = smoothing.iterations                          
                          )
      e.list[[nomeROI]]<-list();
      e.list[[nomeROI]]$volumeratio<-a$volumeratio
      e.list[[nomeROI]]$max.cum.tot<-a$max.cum.tot
      e.list[[nomeROI]]$max.diff.tot<-a$max.diff.tot
      e.list[[nomeROI]]$total.volume<-a$total.volume

    }
  }
  return(e.list);
}
test.smoothing<-function( obj_geoLet, ROINameArray, plotIt = TRUE, newPixelSpacing = NA,
                          verbose=FALSE, forceReCalculus=FALSE, fastEngine = TRUE,
                          decimation=FALSE, 
                          smoothing.iterations.array, decimation.percentage.array) {
  
  tabella<-list();
  grid.espansa<-expand.grid(smoothing.iterations.array, decimation.percentage.array )
  
#  for(index in seq(1,length(smoothing.iterations.array))  ) {
  for(index in seq(1,dim(grid.espansa)[1])  ) {

    smoothIterations<-grid.espansa[index,1]
    decimation.percentage<-grid.espansa[index,2]
    uppa<-system.time(a<-checkDifferencesForAllROIs(obj_geoLet = obj,ROINameArray = ROINameArray,
                               plotIt = plotIt, newPixelSpacing = newPixelSpacing,
                               verbose=verbose, forceReCalculus=TRUE, fastEngine = fastEngine,
                               decimation=decimation, decimation.percentage=decimation.percentage, 
                               smoothing = TRUE,smoothing.iterations = smoothIterations))

    for(ROIName in ROINameArray) {
      if( length(tabella[[ROIName]])==0 ) tabella[[ROIName]]<-c();
      riga<-c( as.numeric(uppa[3]),smoothIterations, decimation.percentage, a[[ROIName]]$total.volume,a[[ROIName]]$volumeratio, a[[ROIName]]$max.cum.tot, a[[ROIName]]$max.diff.tot )
      tabella[[ROIName]]<-rbind(tabella[[ROIName]],riga)
    }
  }
  return(  tabella  );
}


# pointLocation<-function( mesh, x, y, z ) {
# 
#   result<-0
#   res<-.C("pointIn3DPolygon",   as.double(mesh[1,]), as.double(mesh[2,]), as.double(mesh[3,]), as.integer(dim(mesh)[2] ),as.double(x), as.double(y), as.double(z), as.integer(result) ) ;
#   
#   
# }