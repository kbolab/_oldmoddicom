LoadAccordingOSType<-function(library.name) {
  if (Sys.info()["sysname"]=="Windows") 
    return(paste(library.name, ".dll", sep=""))
  if (Sys.info()["sysname"]=="Linux")
    return(paste(library.name, ".so", sep=""))
}

# function for PIP detection using a C algorithm
# returns a vector of points: 0 - outside, 1 - inside the polygons
MultiPointInPoly<-function(X, Y, totalX, totalY, NumSlices, Offset, FullZ) {
  dyn.load(LoadAccordingOSType(library.name = "PointInPolygon"))  
  nX<-length(X)
  nY<-length(Y)
  PIPvector<-rep.int(x = 0, times = length(X)*length(Y)*NumSlices)
  nVert=Offset[length(Offset)]
  result<-.C("MultiPIP", as.double(X), as.double(Y), as.double(totalX), as.double(totalY), as.integer(nX), as.integer(nY), 
             as.integer(NumSlices), as.integer(Offset), as.integer(PIPvector), as.integer(FullZ))  
  return(result[[9]])
}

MultiPointInPolyDist<-function(X, Y, totalX, totalY, NumSlices, Offset, FullZ, Origin, delta.z) {
  dyn.load(LoadAccordingOSType(library.name = "PointInPolygon"))
  nX<-length(X)
  nY<-length(Y)
  PIPvector<-rep.int(x = 0, times = length(X)*length(Y)*NumSlices)
  DISTvector<-rep.int(x = -delta.z, times = length(X)*length(Y)*NumSlices)
  nVert=Offset[length(Offset)]
  result<-.C("MultiPIPDist", as.double(X), as.double(Y), as.double(totalX), as.double(totalY), as.integer(nX), as.integer(nY), 
             as.integer(NumSlices), as.integer(Offset), as.integer(PIPvector), as.integer(FullZ), as.double(DISTvector), as.double(Origin))  
  return(list(PIP=result[[9]], distance=result[[11]]))
}

point.distance<-function(aX, aY, bX, bY) {
  dyn.load(LoadAccordingOSType(library.name = "PointInPolygon"))
  dista<-0
  return(.C("dist", as.double(aX), as.double(aY), as.double(bX), as.double(bY), as.double(dista))[[5]])
}

ComplexEuclNorm<-function(Real, Imaginary) {
  dyn.load(LoadAccordingOSType(library.name = "PointInPolygon"))
  result<-0
  return(.C("QNorm", as.double(Real), as.double(Imaginary), as.double(result))[[3]])
}

# function for calculating the Volume of a Mesh
StructureVolume<-function(mesh, measure.unit=c("cm3", "mm3")) {
  if (class(x = mesh)!="mesh3d") stop("mesh isn't a mesh3d object")
  dyn.load(LoadAccordingOSType(library.name = "PointInPolygon"))
  measure.unit<-match.arg(arg = measure.unit)
  if (measure.unit=="mm3") mu<-1
  if (measure.unit=="cm3") mu<-1000
  Volume=0
  # MeshVolume(double *X, double *Y, double *Z, int *numT, int *V1, int *V2, int *V3, double *Volume)
  return(.C("MeshVolume", as.double(mesh$vb[1,]), as.double(mesh$vb[2,]), as.double(mesh$vb[3,]),
            as.integer(ncol(mesh$it)), as.integer(mesh$it[1,]-1), as.integer(mesh$it[2,]-1),
            as.integer(mesh$it[3,]-1), as.double(Volume))[[8]]/mu)  
}

# function for calculating the surface of a mesh
StructureSurface<-function(mesh, measure.unit=c("cm2", "mm2")) {
  if (class(x = mesh)!="mesh3d") stop("mesh isn't a mesh3d object")
  dyn.load(LoadAccordingOSType(library.name = "PointInPolygon"))
  measure.unit<-match.arg(arg = measure.unit)
  if (measure.unit=="mm2") mu<-1
  if (measure.unit=="cm2") mu<-100
  Surface=0
  # MeshVolume(double *X, double *Y, double *Z, int *numT, int *V1, int *V2, int *V3, double *Volume)
  return(.C("MeshSurface", as.double(mesh$vb[1,]), as.double(mesh$vb[2,]), as.double(mesh$vb[3,]),
            as.integer(ncol(mesh$it)), as.integer(mesh$it[1,]-1), as.integer(mesh$it[2,]-1), 
            as.integer(mesh$it[3,]-1), as.double(Surface))[[8]]/mu)
}

resampleCubeC<-function (xNVoxel,yNVoxel,zNVoxel,xDim,yDim,zDim,newXNVoxel,newYNVoxel,newZNVoxel,values,returnMatrix){
  dyn.load(LoadAccordingOSType(library.name = "moddicom"))
  result<-.C("trilinearInterpolator",as.integer(xNVoxel),as.integer(yNVoxel),as.integer(zNVoxel),as.double(xDim),as.double(yDim),as.double(zDim),
             as.integer(newXNVoxel),as.integer(newYNVoxel),as.integer(newZNVoxel),as.double(values),returnMatrix);result[11][[1]] 
}

# function that calculates the distance 
# from a Point P (defined as a vector of length 2 (x, y)) to a Polygon (defined as matrix)
# without doubling last vertex
Dist.from.Polygon<-function(Polygon, num=400) {
  dyn.load(LoadAccordingOSType(library.name = "PointInPolygon"))
  Pl<-nrow(Polygon) # length of Polygon
  PointV<-rep.int(x = 0, times = num*num)
  totalX<-Polygon[,1] # vector of X coordinates
  totalY<-Polygon[,2] # vector of Y coordinates
  X<-seq(from = min(totalX) - 5, to = max(totalX) + 5, length.out = num)
  Y<-seq(from = min(totalY) - 5, to = max(totalY) + 5, length.out = num)
  result<-.C("DistPolygon", as.double(totalX), as.double(totalY), as.double(X), 
             as.double(Y), as.double(PointV), as.integer(num), as.integer(Pl))[[5]]
  return(matrix(data = result, nrow = num))
}

# function for converting triangle3d to mesh
triangle2mesh <- function(x) {
  v <- list()
  n <- nrow(x$v1)
  nit <- 1:n
  v$vb <- t(cbind(rbind(x$v1,x$v2,x$v3),1))
  v$it <- rbind(nit,nit+n,nit+2*n)
  class(v) <- "mesh3d"
  return(v)
}

# ROIPointList:       list of structure point as achieved by moddicom getROIPointList function
# num:                number of point for sampling Point inP olygon function
# threshold:          threshold for applying Marching Cube algorithm
# smoothing:          smooth the final mesh
# geoLet:             an object of class geoLet
# Structure.Index:    index of structure of interest in geoLet object
# Volume.Resolution:  resolution for calculating the volume of a structure in cm3
RTStruct.Mesh<-function(geoLet, Structure.Index, num=200, threshold=1, smoothing=TRUE, iterations=20, 
                            algorithm=c("PIP", "distance") , ...) {
  algorithm<-match.arg(arg = algorithm)
  message("Preprocessing...")
  ROIPointList<-geoLet$getROIPointList_v2(Structure.Index)
  CTseqZ<-sort(as.numeric(names(geoLet$getAttribute("CT")))) # vector of CT images z coordinates
  # find minima  
  minX<-NULL
  for (n in 1:length(ROIPointList)) {
    minX<-min(minX, min(sapply(X = ROIPointList[[n]], FUN = function(x) min(x[,1]))))
  }
  minY<-NULL
  for (n in 1:length(ROIPointList)) {
    minY<-min(minY, min(sapply(X = ROIPointList[[n]], FUN = function(x) min(x[,2]))))
  }
  minZ<-min(CTseqZ)
    
  # find maxima
  maxX<-NULL
  for (n in 1:length(ROIPointList)) {
    maxX<-max(maxX, max(sapply(X = ROIPointList[[n]], FUN = function(x) max(x[,1]))))
  }
  maxY<-NULL
  for (n in 1:length(ROIPointList)) {
    maxY<-max(maxY, max(sapply(X = ROIPointList[[n]], FUN = function(x) max(x[,2]))))
  }
  maxZ<-max(CTseqZ)

  # creation of a proportioned cube 
  # find maximum width between X and Y
  MAL<-max(maxX-minX, maxY-minY)  
  # axial vector of coordinates for PIP increased by 10 voxels
  axial.coor<-cbind(seq(from = minX-5, to = minX+MAL+5, length.out = num), 
                    seq(from = minY-5, to = minY+MAL+5, length.out = num))
  # creates matrix of points to be checked in polygon
  point.v<-NULL
  for (n in 1:num) 
    point.v<-rbind(point.v, cbind(axial.coor[,1], rep.int(x = axial.coor[n,2], times = num)))
  # vector with results of test of point in polygon 2D for all Z slices

  # creates point in polygon array
  require(Rvcg)
  require(rgl)

  #########################################################################
  # build the Offset vector and totalX, totalY, totalZ vectors version #2 #
  #########################################################################
  
  StructSeqZ<-c()  # sequence of Z coordinates having a contour delineated
  for (n in 1:length(ROIPointList)) {
    StructSeqZ<-c(StructSeqZ, ROIPointList[[n]][[1]][1,3])
  }
  FullZ<-as.numeric(!is.na(match(CTseqZ, StructSeqZ)))  # vector of full and empty slices as 1 (FULL) and 0 (EMPTY)
  Offset<-0
  OriginX<- minX - 10 # origin X point for closing the polygons
  OriginY<- minY - 10 # origin Y point for closing the polygons
  totalX<-minX - 10   
  totalY<-minY - 10
  p = 1  # Slices counter
  n = 1  # contours counter
  while (p <= length(CTseqZ)) {
    if (FullZ[p]==0) {
      Offset<-c(Offset, Offset[length(Offset)] + 1)  # increase the length of vectors by 1
      totalX<-c(totalX, OriginX)
      totalY<-c(totalY, OriginY)
    }
    if (FullZ[p]==1) {
      tempOffset<-0
      for (m in 1:length(ROIPointList[[n]])) {
        tempOffset<-tempOffset + nrow(ROIPointList[[n]][[m]]) + 1 # Offset increases every slice
        totalX<-c(totalX, ROIPointList[[n]][[m]][,1])             # contours increase even more than once per slice  
        totalY<-c(totalY, ROIPointList[[n]][[m]][,2])
        totalX<-c(totalX, OriginX)                                # closes the Polygon by adding one origin
        totalY<-c(totalY, OriginY)      
      }
      Offset<-c(Offset, Offset[length(Offset)] + tempOffset)
      n = n + 1
    }
    p = p + 1
  }
  #browser()
  message("\nSampling points in axial contours...")  
  # use C function for sampling 3D PIP along the slices
  if (algorithm=="PIP")
    pip.v<-MultiPointInPoly(X = axial.coor[,1], Y = axial.coor[,2], totalX = totalX, totalY = totalY, 
                            NumSlices = length(CTseqZ), Offset = Offset, FullZ = FullZ)
  # calculates the difference between z coordinate in two adjacent layers
  delta.z<-abs(CTseqZ[2]-CTseqZ[1])
  if (algorithm=="distance") {
    #X, Y, totalX, totalY, NumSlices, Offset, FullZ, Origin
    pip.v.dist<-MultiPointInPolyDist(X = axial.coor[,1], Y = axial.coor[,2], totalX = totalX, totalY = totalY, 
                                 NumSlices = length(CTseqZ), Offset = Offset, FullZ = FullZ, Origin = OriginX, delta.z=delta.z)
    pip.v<-pip.v.dist[[1]]
    distances<-pip.v.dist[[2]]
  }
  
  message("\nBuilding 3D array of ROI...")
  # at the end of point sampling top and down blank layers are added and the final array is computed
  if (algorithm=="PIP") {
    pip.arr <-array(data = c(rep.int(x = 0, times = 2*num^2), pip.v, rep.int(x = 0, times = 2*num^2)), 
                    dim = c(num, num, length(CTseqZ)+4))  # array of point in polygon with two empty slices added
    pip.arr<-pip.arr*100  # needed for vcgIsosurface function    
    message("\nCreating 3D mesh...")
    # creates evaluation points coordinates
  #  mesh<-vcgIsosurface(vol = pip.arr, lower = threshold, origin = c(minX - 5, minY - 5, minZ - (delta.z * 2)),
  #                      spacing = c(abs(axial.coor[2,1]-axial.coor[1,1]), 
  #                                  abs(axial.coor[2,2]-axial.coor[1,2]),
  #                                  delta.z))
    require("misc3d")
    mesh<-contour3d(f = pip.arr, level = 50, x = axial.coor[,1], y = axial.coor[,2], 
                    z = c(minZ-delta.z*2, minZ-delta.z, CTseqZ, maxZ+delta.z, maxZ+delta.z*2), 
                    engine = "none")
    message("\nArranging mesh3d object...")
    mesh<-triangle2mesh(x = mesh) # correction for the number of polygons
    message("\nCleaning the mesh from supernumerary vertices...")
    mesh<-vcgClean(mesh = mesh, sel = c(0,0,1,1,2,2,3,3,4,4,5,5,6,6,0,0))
  }
  #browser()
  if (algorithm=="distance") {
    distances<-distances+50
    #distances[which(distances<0)]<-0
    pip.arr <-array(data = c(rep.int(x = -delta.z, times = 2*num^2), distances, rep.int(x = -delta.z, times = 2*num^2)), 
                    dim = c(num, num, length(CTseqZ)+4))  # array of point in polygon with two empty slices added
    #browser()
    message("\nCreating 3D mesh...")
    # creates evaluation points coordinates
    mesh<-vcgIsosurface(vol = pip.arr, lower = 50, origin = c(minX - 5, minY - 5, minZ - (delta.z * 2)),
                        spacing = c(abs(axial.coor[2,1]-axial.coor[1,1]), 
                                    abs(axial.coor[2,2]-axial.coor[1,2]),
                                    delta.z))
  }
  

  mesh.orig<-mesh
  message("\nDrawing OpenGL view...")
  if (smoothing==FALSE) {
    shade3d(x = mesh, col="grey", ...)
  }  
  else {
    mesh<-vcgSmooth(mesh = mesh, type = "HClaplace", iteration = iterations)
    shade3d(x = mesh, col="grey", ...)
  }

  return(list(minX=minX,minY=minY,minZ=minZ,maxX=maxX,maxY=maxY,maxZ=maxZ,
              deltaX=abs(axial.coor[2,1]-axial.coor[1,1]), deltaY=abs(axial.coor[2,2]-axial.coor[1,2]),
              deltaZ=delta.z, Max.Axial.Length=MAL,axial.coor=axial.coor, point.v=point.v,pip.arr=pip.arr, 
              mesh=mesh, mesh.orig=mesh.orig, CTseqZ=CTseqZ))
}

RTStruct.Lines<-function(geoLet, Structure.Index, add=TRUE, ...) {
  ROIPointList<-geoLet$getROIPointList_v2(Structure.Index)
  lines3d(x = ROIPointList[[1]][[1]][,1], y = ROIPointList[[1]][[1]][,2], z = ROIPointList[[1]][[1]][,3], add=add, ...)
  for (n in 2:length(ROIPointList)) {
    for (m in 1:length(ROIPointList[[n]]))
      lines3d(x = ROIPointList[[n]][[m]][,1], y = ROIPointList[[n]][[m]][,2], z = ROIPointList[[n]][[m]][,3], add=add, ...)
  }
}


# function that returns the correct value of x, y, z coordinates
# and the permuted Dose matrix of a geoLet object
AlignedDICOMGeometry<-function(geoLet) {
  if (is.list(geoLet$getAttribute("RD"))) 
    DoseGrid<-geoLet$getAttribute("RD")$img
  else stop("\nNo RT Dose image in geoLet object.")
  # starting point of the Dose grid
  x0<- geoLet$getAttribute("geometry")$ImagePositionPatient_RD[1] 
  y0<- geoLet$getAttribute("geometry")$ImagePositionPatient_RD[2] 
  # y is the sequence of x for trnasposed matrices
  y<-seq(from = (dim(DoseGrid)[1] - 1)*geoLet$getAttribute("geometry")$PixelSpacing_RD[2] + y0, to = y0, 
         by =  -geoLet$getAttribute("geometry")$PixelSpacing_RD[2])
  
  x<-seq(from = x0, to = (dim(DoseGrid)[2] - 1)*geoLet$getAttribute("geometry")$PixelSpacing_RD[1] + x0, 
         by =  geoLet$getAttribute("geometry")$PixelSpacing_RD[1])

  z<-sort(geoLet$getAttribute("geometry")$GridFrameOffsetVector_RD)
  # permutation o f3D dose matrix and calculation of correct Dose values according $DoseGridScaling_RD
  Dose3D<-aperm(a = geoLet$getAttribute("RD")$img * geoLet$getAttribute("geometry")$DoseGridScaling_RD, perm = c(2,1,3))  
  return(list(x=x, y=y, z=z, Dose3D=Dose3D))
}
# 
# # function for displaying the 3D dose in different modalities
# # geoLet:           an object of class "geoLet"
# # dose.mode:        how to consider the number in "dose.levels", if "absolute" they are considered as absolute Dose levels
# #                   in measured in Gy (e.g. c(10,20) is considered 10 Gy, 20 Gy); if "relative" they are considered 
# #                   level of percentage (e.g. c(10, 20) is considered as .1 (10%) and .2 (20%) respectively)
# # dose.levels:      a vector containing the values of doses to compute the isosurfaces with dose.mode="isosurface"  OR
# #                   a vector of percent levels for computing the boundaries of doses to be showed in dose.mode="pointcloud".
# #                   In this last case only the higher and lower levels are chosen to define the boundaries.
# #                   If a single number if given in "pointcloud" type it is considered as lower boundary of cloud
# # oversampl.factor: number for interpolating the dose matrix
# # density:          density of points, used to reduce the probability to show a dose point in "pointcloud" type
# # color.scale:      color scale string for colorRamp plotting
# ShowDose3D<-function (geoLet, type=c("pointcloud", "isosurface", "isodoses"), dose.mode=c("absolute", "relative"), 
#                       color.scale=c("black", "darkblue","blue","deepskyblue","green","yellow","orange","red","darkred","white"),
#                       dose.levels, oversampl.factor=2, density, ...) {
#   # match variables
#   type<-match.arg(arg = type)
#   dose.mode<-match.arg(arg = dose.mode)
#   if (missingArg(symbol = geoLet)) stop("\nargument \"geoLet\" is missing, with no default")
#   # check for negative or zero isodoses
#   if (length(which(dose.levels<=0))>0) stop("\nnegative or zero dose levels not allowed")
#   # align Dose3D object
#   Dose3D<-AlignedDICOMGeometry(geoLet = geoLet)  
#   if (nrow(geoLet$getAttribute("geometry")$PlanData_RP)>1) {
#     # check the number of prescription points
#     print(geoLet$getAttribute("geometry")$PlanData_RP)
#     cat("\nThere is more than one Reference Dose item for prescription.",
#         "\nPlease enter the item number for referring 100% of prescription dose:")
#     # request an input from the user if the number of items is higher than 1
#     RefItem<-readline(prompt = "")
#     RefItem<-as.integer(x = RefItem)
#     if ((RefItem<1)||(RefItem>nrow(geoLet$getAttribute("geometry")$PlanData_RP))) stop("\nInvalid item number.")
#     # calculate prescribed dose
#     PrescrDose<-as.numeric(geoLet$getAttribute("geometry")$PlanData_RP[RefItem,5])
#     # calculate maximum Dose
#     MaxDose<-max(Dose3D$Dose3D)
#     # check that the prescribed dose is a numeric value
#     if ((!is.numeric(x = PrescrDose))||is.na(x = PrescrDose)) stop("\nInvalid prescription dose")    
#   }    
#   #################################################
#   #            PointCloud type plot               #
#   #################################################
#   if (type=="pointcloud") {
#     # default boundaries of pointcloud from 10% to maximum
#     if (missingArg(dose.levels)) dose.levels<- c(10, MaxDose)
#     if (length(dose.levels)>1) {
#       # if length of dose.levels = 2 then dose boundaries are from minimum to maximum of dose.levels
#       if (length(dose.levels)>2) message("\nDose levels length > 2: only min and max of dose.levels considered for boudaries of pointcloud color")
#       if (dose.mode=="absolute") dose.levels<- c(min(dose.levels)/PrescrDose, max(dose.levels)/PrescrDose)
#       if (dose.mode=="relative") dose.levels<- c(min(dose.levels)/100, max(dose.levels)/100)
#     }    
#     # if length of dose.levels = 1 then dose.levels is considered as minimum dose for pointcloud
#     if (length(dose.levels)==1) {
#       if (dose.mode=="absolute") dose.levels<- c(dose.levels/PrescrDose, MaxDose/PrescrDose)
#       if (dose.mode=="relative") dose.levels<- c(dose.levels/100, MaxDose/PrescrDose)
#     }
#     # steps of Dose matrix
#     delta.x<-(Dose3D$x[2]-Dose3D$x[1])/oversampl.factor
#     delta.y<-(Dose3D$y[2]-Dose3D$y[1])/oversampl.factor
#     delta.z<-(Dose3D$z[2]-Dose3D$z[1])/oversampl.factor
#     newx<-seq(from = Dose3D$x[1], to = Dose3D$x[length(Dose3D$x)], by = delta.x)
#     newy<-seq(from = Dose3D$y[1], to = Dose3D$y[length(Dose3D$y)], by = delta.y)
#     newz<-seq(from = Dose3D$z[1], to = Dose3D$z[length(Dose3D$z)], by = delta.z)
#     # increment always positive
#     delta.x<-abs(delta.x)
#     delta.y<-abs(delta.y)
#     delta.z<-abs(delta.z)
#     newDimX<-length(newx)
#     newDimY<-length(newy)
#     newDimZ<-length(newz)
#     # calculating the interpolated dose matrix
#     matrix2Return<-array(0,c(newDimX,newDimY,newDimZ))
#     message("\nInterpolating dose matrix...")
#     matrix2Return<-resampleCubeC(dim(Dose3D$Dose3D)[1], dim(Dose3D$Dose3D)[2], dim(Dose3D$Dose3D)[3],
#                                  delta.x, delta.y, delta.z, newDimX, newDimY, newDimZ, 
#                                  Dose3D$Dose3D, matrix2Return)
#     message("\nSampling dose points...")
#     # dose vconverted to relative
#     matrix2Return<-matrix2Return/PrescrDose
#     seeds<-runif(n = length(matrix2Return))
#     # creates plot probability vector
#     if (missingArg(density)) density<-.05
#     PlotV<-(matrix2Return-dose.levels[1])/(dose.levels[2]-dose.levels[1])*density
#     PlotV[which(matrix2Return<dose.levels[1])]<-0
#     # PlotV[which(matrix2Return>dose.levels[2])]<-density  # shows points over higher limit
#     PlotV[which(matrix2Return>dose.levels[2])]<-0   # hides points over higher limit
#     # creating vector of colors
#     message("\nCreating vector of points colors...")
#     colV<-rgb(colorRamp(color.scale)
#       (PlotV/density)/255)
#     message("\nBuilding plot points vectors...")
#     PlotV<-as.numeric(PlotV>seeds)
#     # creating vector of points coordinates
#     xyz<-expand.grid(x=newx, y=newy, z=newz, KEEP.OUT.ATTRS = FALSE)
#     # creating vector of final coordinates and colors    
#     xV<-xyz$x[which(PlotV==1)]
#     yV<-xyz$y[which(PlotV==1)]
#     zV<-xyz$z[which(PlotV==1)]
#     colV<-colV[which(PlotV==1)]
#     # remove memory consuming elements before plotting
#     rm(PlotV, xyz, seeds, matrix2Return)
#     require("rgl")
#     points3d(x = xV, y = yV, z = zV, color=colV, size=5, alpha=.75, add=TRUE)
#     # computing the color legend
#     LegendDoses<-seq(from = dose.levels[1], to = dose.levels[2], by = (dose.levels[2]-dose.levels[1])/2000)
#     # creating legend adding values to maximum if needed 
#     colValue<-(LegendDoses-dose.levels[1])/(dose.levels[2]-dose.levels[1])
#     LegendColor<-rgb(colorRamp(color.scale)(colValue)/255)
#     par(mar=c(5,4,4,5)+.1, bg="black")
#     plot(x = c(-1,1), y = c(dose.levels[1], dose.levels[1]), type="l", col=LegendColor[1], axes = F,
#          ylim=c(dose.levels[1], dose.levels[2]), xlab="", ylab="", xlim=c(-1,1))
#     
#     # axis for Absolute Dose
#     axis(side = 2, at = pretty(x = LegendDoses), labels = pretty(x = LegendDoses) * PrescrDose, col="white", col.axis="white")
#     mtext(text = "Dose [Gy]", side = 2, line = 3, col="white")
#     # Axis for relative Dose
#     axis(side = 4, at = pretty(x = LegendDoses), labels = pretty(x = LegendDoses*100), line = 0, col="white", col.axis="white")
#     mtext(text = "Dose [% of Prescription]", side = 4, line=3, col="white")
#     #lapply(X = LegendDoses, FUN = lines, x = c(-1,1), y = c(LegendDoses, LegendDoses), col=LegendColor)
#     for (n in 2:length(LegendDoses)) lines(x = c(-1,1), y = c(LegendDoses[n], LegendDoses[n]), col = LegendColor[n])
#   }
#   
#   #################################################
#   #            Isosurface type plot               #
#   #################################################
#   if (type=="isosurface") {
#     dose.levels<-sort(dose.levels)
#     if (dose.mode=="absolute") dose.levels<-dose.levels/PrescrDose
#     if (dose.mode=="relative") dose.levels<-dose.levels/100
#     # sort the dose levels
#     dose.levels<-sort(dose.levels)
#     # check for extra levels to be deleted
#     if (max(dose.levels) > MaxDose/PrescrDose) 
#       dose.levels<-dose.levels[1:(max(which(dose.levels<=MaxDose/PrescrDose)))]
#     # normalize doses to maximum and sets isodoses colors
#     colV<-rgb(colorRamp(color.scale)(dose.levels/(MaxDose/PrescrDose))/255)
#     # set alpha level of isosurfaces    
#     if (missingArg(density)) density<-.1
#     require("misc3d")
#     require("Rvcg")
#     for (n in 1:length(dose.levels)) {
#       mesh<-contour3d(f = Dose3D$Dose3D/PrescrDose, level = dose.levels[n], x = Dose3D$x, y = Dose3D$y, z = Dose3D$z, engine = "none")
#       mesh<-triangle2mesh(x = mesh) # correction for the number of polygons
#       mesh<-vcgClean(mesh = mesh, sel = c(0,0,1,1,2,2,3,3,4,4,5,5,6,6,0,0))  # clean the mesh
#       shade3d(x = mesh, color = colV[n], alpha=density)
#     }
#     # plot the legenda
#     par(mar=c(5,4,4,5)+.1, bg="black")    
#     plot(x = c(-1,1), y = c(dose.levels[1]*PrescrDose, dose.levels[1]*PrescrDose), 
#          col=colV[1], type="l", lwd=3, ylim=c(min(dose.levels)*PrescrDose , max(dose.levels)*PrescrDose), xlab="", ylab="", axes=F)
#     # axis for absolute Dose
#     axis(side = 2, at = dose.levels * PrescrDose, labels = dose.levels * PrescrDose, col="white", col.axis="white")
#     mtext(text = "Dose [Gy]", side = 2, line = 3, col="white")
#     # axis for relative Dose
#     axis(side = 4, at = dose.levels * PrescrDose, labels = dose.levels*100, line = 0, col="white", col.axis="white")
#     mtext(text = "Dose [% of Prescription]", side = 4, line=3, col="white")
#     if (length(dose.levels>1)) for (n in 2:length(dose.levels)) lines(x = c(-1,1), 
#                                                                       y = c(dose.levels[n]*PrescrDose, dose.levels[n]*PrescrDose),
#                                                                       col=colV[n], lwd=4)
#   }
#   #################################################
#   #             Isodoses type plot                #
#   #################################################
#   if (type=="isodoses") {
#     dose.levels<-sort(dose.levels)
#     if (dose.mode=="absolute") {
#       dose.levels<-dose.levels/PrescrDose
#       Dose3D$Dose3D<-Dose3D$Dose3D/PrescrDose
#     }
#     if (dose.mode=="relative") {
#       dose.levels<-dose.levels/100
#       Dose3D$Dose3D<-Dose3D$Dose3D/PrescrDose
#     }
#     # sort the dose levels
#     dose.levels<-sort(dose.levels)
#     # check for extra levels to be deleted
#     if (max(dose.levels) > MaxDose/PrescrDose) 
#       dose.levels<-dose.levels[1:(max(which(dose.levels<=MaxDose/PrescrDose)))]
#     # normalize doses to maximum and sets isodoses colors
#     colV<-rgb(colorRamp(color.scale)(dose.levels/(MaxDose/PrescrDose))/255)
#     require("misc3d")
#     require("Rvcg")
#     delta.x<-(Dose3D$x[2]-Dose3D$x[1])
#     delta.y<-(Dose3D$y[2]-Dose3D$y[1])
#     # check the sign of delta.x and delta.y, needed to apply contourLines function
#     if (delta.x<0) x<- -Dose3D$x else x<- Dose3D$x
#     if (delta.y<0) y<- -Dose3D$y else y<- Dose3D$x
#     # calculate the isodoses in the plot
#     isodoses<-apply(X = Dose3D$Dose3D, MARGIN = 3, FUN = contourLines, x = x, y = y, levels = dose.levels)
#     # correct X values
#     if (delta.x<0)
#       for (n in 1:length(isodoses)) {
#         if (length(isodoses[[n]])>0)         
#           for (m in 1:length(isodoses[[n]])) isodoses[[n]][[m]]$x<- - isodoses[[n]][[m]]$x
#       }
#     # correct Y values
#     if (delta.y<0)
#       for (n in 1:length(isodoses)) {
#         if (length(isodoses[[n]])>0)         
#           for (m in 1:length(isodoses[[n]])) isodoses[[n]][[m]]$y<- - isodoses[[n]][[m]]$y    
#       }
#     # compute colors
#     # sort the dose levels
#     dose.levels<-sort(dose.levels)
#     # check for extra levels to be deleted
#     if (max(dose.levels) > MaxDose/PrescrDose) 
#       dose.levels<-dose.levels[1:(max(which(dose.levels<=MaxDose/PrescrDose)))]
#     # normalize doses to maximum and sets isodoses colors
#     colV<-rgb(colorRamp(color.scale)(dose.levels/(MaxDose/PrescrDose))/255)
#     # plot isodose lines
#     require("rgl")
#     for (n in 1:length(isodoses)) {
#       if (length(isodoses[[n]])>0) {
#         for (m in 1:length(isodoses[[n]])) {
#           # set the z level for plotting isodoses
#           z<-rep.int(x = Dose3D$z[n], times = length(isodoses[[n]][[m]]$x))
#           lines3d(x = isodoses[[n]][[m]]$x, y = isodoses[[n]][[m]]$y, z = z, 
#                   color=colV[which(dose.levels==isodoses[[n]][[m]]$level)], add=TRUE)
#         }
#       }
#     }
#     # plot the legenda
#     par(mar=c(5,4,4,5)+.1, bg="black")    
#     plot(x = c(-1,1), y = c(dose.levels[1]*PrescrDose, dose.levels[1]*PrescrDose), 
#          col=colV[1], type="l", lwd=3, ylim=c(min(dose.levels)*PrescrDose , max(dose.levels)*PrescrDose), xlab="", ylab="", axes=F)
#     # axis for absolute Dose
#     axis(side = 2, at = dose.levels * PrescrDose, labels = dose.levels * PrescrDose, col="white", col.axis="white")
#     mtext(text = "Dose [Gy]", side = 2, line = 3, col="white")
#     # axis for relative Dose
#     axis(side = 4, at = dose.levels * PrescrDose, labels = dose.levels*100, line = 0, col="white", col.axis="white")
#     mtext(text = "Dose [% of Prescription]", side = 4, line=3, col="white")
#     if (length(dose.levels>1)) for (n in 2:length(dose.levels)) lines(x = c(-1,1), 
#                                                                       y = c(dose.levels[n]*PrescrDose, dose.levels[n]*PrescrDose),
#                                                                       col=colV[n], lwd=4)
#   }
#   
#   # final values to return in the console
#   cat("\nMaximum Dose =    ", MaxDose, " Gy")
#   cat("\nPrescribed Dose = ", PrescrDose, " Gy")  
# }


# function for displaying the 3D dose in different modalities
# geoLet:           an object of class "geoLet"
# dose.mode:        how to consider the number in "dose.levels", if "absolute" they are considered as absolute Dose levels
#                   in measured in Gy (e.g. c(10,20) is considered 10 Gy, 20 Gy); if "relative" they are considered 
#                   level of percentage (e.g. c(10, 20) is considered as .1 (10%) and .2 (20%) respectively)
# dose.levels:      a vector containing the values of doses to compute the isosurfaces with dose.mode="isosurface"  OR
#                   a vector of percent levels for computing the boundaries of doses to be showed in dose.mode="pointcloud".
#                   In this last case only the higher and lower levels are chosen to define the boundaries.
#                   If a single number if given in "pointcloud" type it is considered as lower boundary of cloud
# density:          mean density of points in space, in pointcloud is equal to mean(number of points)/mm3
# color.scale:      color scale string for colorRamp plotting
ShowDose3D<-function (geoLet, type=c("pointcloud", "isosurface", "isodoses"), dose.mode=c("absolute", "relative"), 
                      color.scale=c("black", "darkblue","blue","deepskyblue","green","yellow","orange","red","darkred","white"),
                      dose.levels, density=1, ...) {
  # match variables
  type<-match.arg(arg = type)
  dose.mode<-match.arg(arg = dose.mode)
  if (missingArg(symbol = geoLet)) stop("\nargument \"geoLet\" is missing, with no default")
  # check for negative or zero isodoses
  if (length(which(dose.levels<=0))>0) stop("\nnegative or zero dose levels not allowed")
  # align Dose3D object
  Dose3D<-AlignedDICOMGeometry(geoLet = geoLet)  
  if (nrow(geoLet$getAttribute("geometry")$PlanData_RP)>1) {
    # check the number of prescription points
    print(geoLet$getAttribute("geometry")$PlanData_RP)
    cat("\nThere is more than one Reference Dose item for prescription.",
        "\nPlease enter the item number for referring 100% of prescription dose:")
    # request an input from the user if the number of items is higher than 1
    RefItem<-readline(prompt = "")
    RefItem<-as.integer(x = RefItem)
    if ((RefItem<1)||(RefItem>nrow(geoLet$getAttribute("geometry")$PlanData_RP))) stop("\nInvalid item number.")
    # calculate prescribed dose
    PrescrDose<-as.numeric(geoLet$getAttribute("geometry")$PlanData_RP[RefItem,5])
    # calculate maximum Dose
    MaxDose<-max(Dose3D$Dose3D)
    # check that the prescribed dose is a numeric value
    if ((!is.numeric(x = PrescrDose))||is.na(x = PrescrDose)) stop("\nInvalid prescription dose")    
  }    
  #################################################
  #            PointCloud type plot               #
  #################################################
  if (type=="pointcloud") {
    # default boundaries of pointcloud from 10% to maximum
    if (missingArg(dose.levels)) dose.levels<- c(10, MaxDose)
    if (length(dose.levels)>1) {
      # if length of dose.levels = 2 then dose boundaries are from minimum to maximum of dose.levels
      if (length(dose.levels)>2) message("\nDose levels length > 2: only min and max of dose.levels considered for boudaries of pointcloud color")
      if (dose.mode=="absolute") dose.levels<- c(min(dose.levels)/PrescrDose, max(dose.levels)/PrescrDose)
      if (dose.mode=="relative") dose.levels<- c(min(dose.levels)/100, max(dose.levels)/100)
    }    
    # if length of dose.levels = 1 then dose.levels is considered as minimum dose for pointcloud
    if (length(dose.levels)==1) {
      if (dose.mode=="absolute") dose.levels<- c(dose.levels/PrescrDose, MaxDose/PrescrDose)
      if (dose.mode=="relative") dose.levels<- c(dose.levels/100, MaxDose/PrescrDose)
    }
    # steps of Dose matrix
    delta.x<-(Dose3D$x[2]-Dose3D$x[1])
    delta.y<-(Dose3D$y[2]-Dose3D$y[1])
    delta.z<-(Dose3D$z[2]-Dose3D$z[1])    
    # create boundaries of matrix to compute the random point
    x0<-min(which(Dose3D$Dose/PrescrDose>=(dose.levels[1]-.075), arr.ind = TRUE)[,1])
    x1<-max(which(Dose3D$Dose/PrescrDose>=(dose.levels[1]-.075), arr.ind = TRUE)[,1])
    y0<-min(which(Dose3D$Dose/PrescrDose>=(dose.levels[1]-.075), arr.ind = TRUE)[,2])
    y1<-max(which(Dose3D$Dose/PrescrDose>=(dose.levels[1]-.075), arr.ind = TRUE)[,2])
    z0<-min(which(Dose3D$Dose/PrescrDose>=(dose.levels[1]-.075), arr.ind = TRUE)[,3])
    z1<-max(which(Dose3D$Dose/PrescrDose>=(dose.levels[1]-.075), arr.ind = TRUE)[,3])
    numPoints<-round((x1-x0)*(y1-y0)*(z1-z0)*density)  # calculates the number of points to be computed
    # create coordinates of points to calculate the dose
    xCoor<-runif(n = numPoints, min = min(Dose3D$x[x0], Dose3D$x[x1]) , max = max(Dose3D$x[x0], Dose3D$x[x1]))
    yCoor<-runif(n = numPoints, min = min(Dose3D$y[y0], Dose3D$y[y1]) , max = max(Dose3D$y[y0], Dose3D$y[y1]))
    zCoor<-runif(n = numPoints, min = min(Dose3D$z[z0], Dose3D$z[z1]) , max = max(Dose3D$z[z0], Dose3D$z[z1]))
    # calculating the interpolated dose matrix
    message("\nInterpolating dose matrix...")
    matrix2Return<-triLinear(A = Dose3D$Dose3D, origin = c(Dose3D$x[1], Dose3D$y[1], Dose3D$z[1]),
                             xVett = xCoor, yVett = yCoor, zVett = zCoor, 
                             dx = delta.x, dy = delta.y, dz = delta.z)
    message("\nSampling dose points...")
    # dose vconverted to relative
    values<-matrix2Return/PrescrDose
    #seeds<-runif(n = length(matrix2Return), min = dose.levels[1], max = dose.levels[2])
    # creates plot probability vector
    PlotV<-(values-dose.levels[1])/(dose.levels[2]-dose.levels[1])
    PlotV[which(values<dose.levels[1])]<-0
    # PlotV[which(matrix2Return>dose.levels[2])]<-density  # shows points over higher limit
    PlotV[which(values>dose.levels[2])]<-0   # hides points over higher limit
    # creating vector of colors
    message("\nCreating vector of points colors...")
    colV<-rgb(colorRamp(color.scale)(PlotV)/255)
    message("\nBuilding plot points vectors...")
    PlotV<-as.numeric(PlotV>0)
    # creating vector of final coordinates and colors    
    xV<-xCoor[which(PlotV==1)]
    yV<-yCoor[which(PlotV==1)]
    zV<-zCoor[which(PlotV==1)]
    colV<-colV[which(PlotV==1)]
    # remove memory consuming elements before plotting
    rm(PlotV, matrix2Return)
    require("rgl")
    points3d(x =xV, y = yV, z = zV, color=colV, size=5, add=TRUE, ...)
    # computing the color legend
    LegendDoses<-seq(from = dose.levels[1], to = dose.levels[2], by = (dose.levels[2]-dose.levels[1])/2000)
    # creating legend adding values to maximum if needed 
    colValue<-(LegendDoses-dose.levels[1])/(dose.levels[2]-dose.levels[1])
    LegendColor<-rgb(colorRamp(color.scale)(colValue)/255)
    par(mar=c(5,4,4,5)+.1, bg="grey30")
    plot(x = c(-1,1), y = c(dose.levels[1], dose.levels[1]), type="l", col=LegendColor[1], axes = F,
         ylim=c(dose.levels[1], dose.levels[2]), xlab="", ylab="", xlim=c(-1,1))
    
    # axis for Absolute Dose
    axis(side = 2, at = pretty(x = LegendDoses), labels = pretty(x = LegendDoses) * PrescrDose, col="white", col.axis="white")
    mtext(text = "Dose [Gy]", side = 2, line = 3, col="white")
    # Axis for relative Dose
    axis(side = 4, at = pretty(x = LegendDoses), labels = pretty(x = LegendDoses*100), line = 0, col="white", col.axis="white")
    mtext(text = "Dose [% of Prescription]", side = 4, line=3, col="white")
    #lapply(X = LegendDoses, FUN = lines, x = c(-1,1), y = c(LegendDoses, LegendDoses), col=LegendColor)
    for (n in 2:length(LegendDoses)) lines(x = c(-1,1), y = c(LegendDoses[n], LegendDoses[n]), col = LegendColor[n])
  }
  
  #################################################
  #            Isosurface type plot               #
  #################################################
  if (type=="isosurface") {
    dose.levels<-sort(dose.levels)
    if (dose.mode=="absolute") dose.levels<-dose.levels/PrescrDose
    if (dose.mode=="relative") dose.levels<-dose.levels/100
    # sort the dose levels
    dose.levels<-sort(dose.levels)
    # check for extra levels to be deleted
    if (max(dose.levels) > MaxDose/PrescrDose) 
      dose.levels<-dose.levels[1:(max(which(dose.levels<=MaxDose/PrescrDose)))]
    # normalize doses to maximum and sets isodoses colors
    colV<-rgb(colorRamp(color.scale)(dose.levels/(MaxDose/PrescrDose))/255)
    # set alpha level of isosurfaces    
    if (missingArg(density)) density<-.1
    require("misc3d")
    require("Rvcg")
    for (n in 1:length(dose.levels)) {
      mesh<-contour3d(f = Dose3D$Dose3D/PrescrDose, level = dose.levels[n], x = Dose3D$x, y = Dose3D$y, z = Dose3D$z, engine = "none")
      mesh<-triangle2mesh(x = mesh) # correction for the number of polygons
      mesh<-vcgClean(mesh = mesh, sel = c(0,0,1,1,2,2,3,3,4,4,5,5,6,6,0,0))  # clean the mesh
      shade3d(x = mesh, color = colV[n], ...)
    }
    # plot the legenda
    par(mar=c(5,4,4,5)+.1, bg="grey30")    
    plot(x = c(-1,1), y = c(dose.levels[1]*PrescrDose, dose.levels[1]*PrescrDose), 
         col=colV[1], type="l", lwd=3, ylim=c(min(dose.levels)*PrescrDose , max(dose.levels)*PrescrDose), xlab="", ylab="", axes=F)
    # axis for absolute Dose
    axis(side = 2, at = dose.levels * PrescrDose, labels = dose.levels * PrescrDose, col="white", col.axis="white")
    mtext(text = "Dose [Gy]", side = 2, line = 3, col="white")
    # axis for relative Dose
    axis(side = 4, at = dose.levels * PrescrDose, labels = dose.levels*100, line = 0, col="white", col.axis="white")
    mtext(text = "Dose [% of Prescription]", side = 4, line=3, col="white")
    if (length(dose.levels>1)) for (n in 2:length(dose.levels)) lines(x = c(-1,1), 
                                                                      y = c(dose.levels[n]*PrescrDose, dose.levels[n]*PrescrDose),
                                                                      col=colV[n], lwd=4)
  }
  #################################################
  #             Isodoses type plot                #
  #################################################
  if (type=="isodoses") {
    # sort the dose levels
    dose.levels<-sort(dose.levels)
    if (dose.mode=="absolute") {
      dose.levels<-dose.levels/PrescrDose
      Dose3D$Dose3D<-Dose3D$Dose3D/PrescrDose
    }
    if (dose.mode=="relative") {
      dose.levels<-dose.levels/100
      Dose3D$Dose3D<-Dose3D$Dose3D/PrescrDose
    }
    # check for extra levels to be deleted
    if (max(dose.levels) > MaxDose/PrescrDose) 
      dose.levels<-dose.levels[1:(max(which(dose.levels<=MaxDose/PrescrDose)))]
    # normalize doses to maximum and sets isodoses colors
    colV<-rgb(colorRamp(color.scale)(dose.levels/(MaxDose/PrescrDose))/255)
    require("misc3d")
    require("Rvcg")
    delta.x<-(Dose3D$x[2]-Dose3D$x[1])
    delta.y<-(Dose3D$y[2]-Dose3D$y[1])
    # check the sign of delta.x and delta.y, needed to apply contourLines function
    if (delta.x<0) x<- -Dose3D$x else x<- Dose3D$x
    if (delta.y<0) y<- -Dose3D$y else y<- Dose3D$x
    # calculate the isodoses in the plot
    isodoses<-apply(X = Dose3D$Dose3D, MARGIN = 3, FUN = contourLines, x = x, y = y, levels = dose.levels)
    # correct X values
    if (delta.x<0)
      for (n in 1:length(isodoses)) {
        if (length(isodoses[[n]])>0)         
          for (m in 1:length(isodoses[[n]])) isodoses[[n]][[m]]$x<- - isodoses[[n]][[m]]$x
      }
    # correct Y values
    if (delta.y<0)
      for (n in 1:length(isodoses)) {
        if (length(isodoses[[n]])>0)         
          for (m in 1:length(isodoses[[n]])) isodoses[[n]][[m]]$y<- - isodoses[[n]][[m]]$y    
      }
    # compute colors
    # sort the dose levels
    dose.levels<-sort(dose.levels)
    # check for extra levels to be deleted
    if (max(dose.levels) > MaxDose/PrescrDose) 
      dose.levels<-dose.levels[1:(max(which(dose.levels<=MaxDose/PrescrDose)))]
    # normalize doses to maximum and sets isodoses colors
    colV<-rgb(colorRamp(color.scale)(dose.levels/(MaxDose/PrescrDose))/255)
    # plot isodose lines
    require("rgl")
    for (n in 1:length(isodoses)) {
      if (length(isodoses[[n]])>0) {
        for (m in 1:length(isodoses[[n]])) {
          # set the z level for plotting isodoses
          z<-rep.int(x = Dose3D$z[n], times = length(isodoses[[n]][[m]]$x))
          lines3d(x = isodoses[[n]][[m]]$x, y = isodoses[[n]][[m]]$y, z = z, 
                  color=colV[which(dose.levels==isodoses[[n]][[m]]$level)], add=TRUE)
        }
      }
    }
    # plot the legenda
    par(mar=c(5,4,4,5)+.1, bg="grey30")    
    plot(x = c(-1,1), y = c(dose.levels[1]*PrescrDose, dose.levels[1]*PrescrDose), 
         col=colV[1], type="l", lwd=3, ylim=c(min(dose.levels)*PrescrDose , max(dose.levels)*PrescrDose), xlab="", ylab="", axes=F)
    # axis for absolute Dose
    axis(side = 2, at = dose.levels * PrescrDose, labels = dose.levels * PrescrDose, col="white", col.axis="white")
    mtext(text = "Dose [Gy]", side = 2, line = 3, col="white")
    # axis for relative Dose
    axis(side = 4, at = dose.levels * PrescrDose, labels = dose.levels*100, line = 0, col="white", col.axis="white")
    mtext(text = "Dose [% of Prescription]", side = 4, line=3, col="white")
    if (length(dose.levels>1)) for (n in 2:length(dose.levels)) lines(x = c(-1,1), 
                                                                      y = c(dose.levels[n]*PrescrDose, dose.levels[n]*PrescrDose),
                                                                      col=colV[n], lwd=4)
  }
  
  # final values to return in the console
  cat("\nMaximum Dose =    ", MaxDose, " Gy")
  cat("\nPrescribed Dose = ", PrescrDose, " Gy")  
}
