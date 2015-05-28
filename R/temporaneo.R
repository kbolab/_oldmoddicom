# 
# rm(obj);rm(dati);rm(a)
# obj<-geoLet()
# #obj$openDICOMFolder("/progetti/immagini/SarhoshTest")
# obj$openDICOMFolder("/progetti/immagini/CONTOURED/Negative/PRE/Alimenti")
# dati <- obj$getAttribute("dataStorage")
# #ROIName<-"Ossa, NOS"
# ROIName<-"GTV"
# orientationMatrix<-dataStorage$info[[1]][[1]]$orientationMatrix
# 
# estrappolaVoxel<-function( obj,   ROIName ,  Nx, Ny, Nz, orientationMatrix, nPoints.x, nPoints.z, nPoints.y  ) {
#   
#   matrice<-c()
#   dataStorage<-obj$getAttribute("dataStorage")
#   ROIPointList<-obj$getROIPointList( ROIName );
# 
#   a11<-orientationMatrix[1,1]; a21<-orientationMatrix[2,1];
#   a31<-orientationMatrix[3,1]; a12<-orientationMatrix[1,2];
#   a22<-orientationMatrix[2,2]; a32<-orientationMatrix[3,2];
#   Ox<orientationMatrix[1,4]; Oy<-orientationMatrix[2,4];  Oy<-orientationMatrix[3,4];
#   
#   coords.x<-seq(from = 0, to = (nPoints.x-1))
#   coords.y<-seq(from = 0, to = (nPoints.y-1))
#   coords.z<-seq(from = 0, to = (nPoints.z-1))
#   coors.grid <- expand.grid( coords.x , coords.y , coords.z  ) 
#   
#   
#   
#   
#   for( i in names(ROIPointList)) {    
#     for( t in seq(1,length( ROIPointList[[i]] ))) {
#       matrice<-rbind(matrice,ROIPointList[[i]][[t]] )
#     }    
#   }
#   
#   box.minX<-min(matrice[,1]); box.maxX<-max(matrice[,1])
#   box.minY<-min(matrice[,2]); box.maxY<-max(matrice[,2])
#   box.minZ<-min(matrice[,3]); box.maxZ<-max(matrice[,3])
#   
#     
#   
#   
#   
#   
#   samples.x<-seq( from=box.minX,  to=box.maxX,  by =  (box.maxX-box.minX)/Nx )
#   samples.y<-seq( from=box.minY,  to=box.maxY,  by =  (box.maxY-box.minY)/Ny )
#   samples.z<-seq( from=box.minZ,  to=box.maxZ,  by =  (box.maxZ-box.minZ)/Nz )
#   Grid<-expand.grid(x = samples.x, y = samples.y , z = samples.z)
#   
# }