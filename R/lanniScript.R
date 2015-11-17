# 
# 
# path.pre<-"/progetti/immagini/lanni/pre/22872220";
# path.post<-"/progetti/immagini/lanni/post/22872220" 
# ROIName<-"PET+"
# 
# # istanzia gli oggetti
# obj.pre <- geoLet();
# obj.post <- geoLet();
# obj.pet <- geoLet();
# 
# # carica il post e prendi il dataStorage: mi servità per prendere il nome della serie
# obj.post$openDICOMFolder( path.post );
# ds.post<-obj.post$getAttribute("dataStorage")
# serieName<-names(ds.post$info)
# 
# # carica il PRE utilizzando il nome della serie del post
# obj.pre$openDICOMFolder(pathToOpen = path.pre,setValidCTRMNSeriesInstanceUID = serieName)
# 
# # fai le verifiche per vedere che non ci siano troppe differenze nelle geometrie
# vox.post<-obj.post$getImageVoxelCube();
# 
# # CHIODO per la PET (dato che ne hanno inserite due serie)
# petSeries<-"1.3.12.2.1107.5.1.4.11088.30000014021409165065600002670";
# 
# obj.pet$openDICOMFolder(pathToOpen = path.pre,setValidCTRMNSeriesInstanceUID = petSeries)
# ds.pet<-obj.pet$getAttribute("dataStorage")
# 
# 
# ss<-services();
# voxelVolume.pet<-obj.pet$getImageVoxelCube()
# dx.pet<-ds.pet$info[[1]][[1]]$pixelSpacing[1];
# dy.pet<-ds.pet$info[[1]][[1]]$pixelSpacing[2];
# dz.pet<-as.numeric(abs(ds.pet$info[[1]][[1]]$ImagePositionPatient[3]-ds.pet$info[[1]][[2]]$ImagePositionPatient[3]))
# 
# 
# dx.post<-ds.post$info[[1]][[1]]$pixelSpacing[1]
# dy.post<-ds.post$info[[1]][[1]]$pixelSpacing[2]
# dz.post<-as.numeric(abs(ds.post$info[[1]][[1]]$ImagePositionPatient[3]-ds.post$info[[1]][[2]]$ImagePositionPatient[3]))
# 
# coordsTop.post<-ds.post$info[[1]][[1]]$ImagePositionPatient
# coordsTop.pet<-ds.pet$info[[1]][[1]]$ImagePositionPatient
# top.x.post<-coordsTop.post[1];  top.y.post<-coordsTop.post[2];  top.z.post<-coordsTop.post[3]
# top.x.pet<-coordsTop.pet[1];  top.y.pet<-coordsTop.pet[2];  top.z.pet<-coordsTop.pet[3]
# 
# top.x<-top.x.post - top.x.pet
# top.y<-top.y.post - top.y.pet
# top.z<-top.z.post - top.z.pet
# seq.x<-seq(from=top.x,to = dim(vox.post)[1]*dx.post+top.x,by = dx.post)
# seq.y<-seq(from=top.y,to = dim(vox.post)[2]*dy.post+top.y,by = dy.post)
# seq.z<-seq(from=0,to = dim(vox.post)[3]*dz.post,by = dz.post)
# 
# 
# VoxelCubePointPos.post<-expand.grid(seq.x,seq.y,seq.z)
# 
# aa<-ss$new.SV.trilinearInterpolator.onGivenPoints(
#   voxelCube = voxelVolume.pet,  pixelSpacing.old = c(dx.pet,dy.pet,dz.pet),
#   newPointCoords.x = seq.x,newPointCoords.y = seq.y,newPointCoords.z = seq.z) 
# 
# bbb<-aa
# for(i in seq(1,dim(bbb)[3])) {
#   bbb[,,i]<-aa[,,dim(bbb)[3]-i+1]
# }
# aa<-bbb
# rm(bbb)
# 
