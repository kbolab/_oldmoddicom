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
      cubeDim<-tmpAlign$cubeDim
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
          ROI<-align[[   arrayROINames[ ROIName]    ]][[internalROIPlanes]][[1]]
          ROI[,1]<-ROI[,1] / cubeDim[1]
          ROI[,2]<-1-ROI[,2] / cubeDim[2]
          points(x = ROI[,1],y = ROI[,2],type='l', lwd = 2, col = paletteColor[colorIncrement])  
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




# pointLocation<-function( mesh, x, y, z ) {
# 
#   result<-0
#   res<-.C("pointIn3DPolygon",   as.double(mesh[1,]), as.double(mesh[2,]), as.double(mesh[3,]), as.integer(dim(mesh)[2] ),as.double(x), as.double(y), as.double(z), as.integer(result) ) ;
#   
#   
# }