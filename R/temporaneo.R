# # 
# obj<-geoLet()
# #obj$openDICOMFolder("/progetti/immagini/SaroshTest")
# obj$openDICOMFolder("/progetti/immagini/doloreAcuto")
# #obj$openDICOMFolder("/progetti/immagini/doloreOttuso")
#  objS<-services();
#  objV<-viewer();
# # 
#  ds<-obj$getAttribute("dataStorage")
# # 
# GTV.voxel<-obj$getROIVoxels(Structure = "Ossa1")
# GTV.rotated<-obj$rotateToAlign(ROIName = "Ossa1")
# GTV.normal<-obj$getROIPointList(ROINumber = "Ossa1")
# 
# CT<-objS$expandCube(littleCube = GTV.voxel$masked.images$voxelCube,x.start = GTV.voxel$masked.images$location$min.x,y.start = GTV.voxel$masked.images$location$min.y,z.start = GTV.voxel$masked.images$location$min.z,fe = GTV.voxel$masked.images$location$fe,se = GTV.voxel$masked.images$location$se,te = GTV.voxel$masked.images$location$te)
# #CT<-obj$getDoseVoxelCube()
# 
# step.x<-obj$getPixelSpacing()[1]
# step.y<-obj$getPixelSpacing()[2]
# step.z<-as.numeric(obj$getAttribute("SliceThickness"))
# x<-seq(from = 0, to = (dim(CT)[1]-1)*step.x, by = step.x)
# y<-seq(from = 0, to = (dim(CT)[2]-1)*step.y, by = step.y)
# z<-seq(from = 0, to = (dim(CT)[3]-1)*step.z, by = step.z)
# contour3d(f = CT,level = 1,x = x,y = y,z = z)
