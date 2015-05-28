#' class for viewing
#' 
#' @description  allows to view what's happening
#' @useDynLib moddicom
#' @import Rvcg rgl 
#' @export
viewer<-function() {
  
  isosurface<-function( matrice , lower , metodo= "wired", percent = 0.5, pixelSpacing, ... ) {
    parametri<-list(...)
    x<-seq(0, dim(matrice)[1] )
    y<-seq(0, dim(matrice)[2] )
    z<-seq(0, dim(matrice)[3] )
    g<-expand.grid(x,y,z)
    v<-as.array( matrice )
    storage.mode(v) <- "integer"
    mesh <- vcgIsosurface( v, lower = lower, spacing = pixelSpacing)
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
      shade3d(vcgSmooth(mesh = mesh, lambda = parametri$lambda, mu = parametri$mu) , col = parametri$col)      
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
  
  
  return( list( isosurface = isosurface )  )
  
}




# pointLocation<-function( mesh, x, y, z ) {
# 
#   result<-0
#   res<-.C("pointIn3DPolygon",   as.double(mesh[1,]), as.double(mesh[2,]), as.double(mesh[3,]), as.integer(dim(mesh)[2] ),as.double(x), as.double(y), as.double(z), as.integer(result) ) ;
#   
#   
# }