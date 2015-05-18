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
    mesh <- vcgIsosurface( v, lower = 1, spacing = pixelSpacing)
    decimface <- vcgQEdecim( mesh, percent = percent)
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

  return( list( isosurface = isosurface )  )
  
}
