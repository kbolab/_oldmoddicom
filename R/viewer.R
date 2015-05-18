#' class for viewing
#' 
#' @description  allows to view what's happening
#' @useDynLib moddicom
#' @import Rvcg rgl 
#' @export
viewer<-function() {
  
  isosurface<-function( matrice , lower , metodo= "wired", percent = 0.5 ) {
    x<-seq(0, dim(matrice)[1] )
    y<-seq(0, dim(matrice)[2] )
    z<-seq(0, dim(matrice)[3] )
    g<-expand.grid(x,y,z)
    v<-as.array( matrice )
    storage.mode(v) <- "integer"
    mesh <- vcgIsosurface(v,lower=1)
    decimface <- vcgQEdecim( mesh, percent = percent)
    if(metodo == "wired" ) wire3d(mesh,col="red")        
  }

  return( list( isosurface = isosurface )  )
  
}