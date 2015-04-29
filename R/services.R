#' class for handling logs/warnings/errorss
#' 
#' @description  It handles messages from script to a chosen output (screen, file, etc.)
#' @export
services<-function() {
    # ------------------------------------------------
  # SV.getPointPlaneDistance
  # ------------------------------------------------
  SV.getPointPlaneDistance<-function(Punto,Piano) {
    return(abs(Piano[1]*Punto[1]+Piano[2]*Punto[2]+Piano[3]*Punto[3]+Piano[4])/sqrt(Piano[1]^2+Piano[2]^2+Piano[3]^2))
  }
  # ------------------------------------------------
  # SV.get3DPosFromNxNy
  # ------------------------------------------------
  SV.get3DPosFromNxNy<-function(Nx,Ny,oM) {
    return(xy<-t(oM%*%c(Nx,Ny,0,1)))
  }
  # ------------------------------------------------
  # SV.getPlaneEquationBetween3Points
  # ------------------------------------------------
  SV.getPlaneEquationBetween3Points<-function(Pa,Pb,Pc) {
    ac<-(Pb[2]-Pa[2])*(Pc[3]-Pa[3])-(Pc[2]-Pa[2])*(Pb[3]-Pa[3])
    bc<-(Pb[3]-Pa[3])*(Pc[1]-Pa[1])-(Pc[3]-Pa[3])*(Pb[1]-Pa[1])
    cc<-(Pb[1]-Pa[1])*(Pc[2]-Pa[2])-(Pc[1]-Pa[1])*(Pb[2]-Pa[2])
    dc<--(ac*Pa[1]+bc*Pa[2]+cc*Pa[3])
    return(c(ac,bc,cc,dc))
  }
  SV.rotateMatrix<-function( m , rotations = 1) {
    if(rotations == 1 ) m<-t(m[nrow(m):1,])
    if(rotations == 2 ) m<-m[nrow(m):1,ncol(m):1]
    if(rotations == 3 ) m<-t(m)[ncol(m):1,]
    return( m )
  }
  SV.LoadAccordingOSType<-function(library.name) {
    if (Sys.info()["sysname"]=="Windows") 
      return(paste(library.name, ".dll", sep=""))
    if (Sys.info()["sysname"]=="Linux")
      return(paste(library.name, ".so", sep=""))
  }
  # ========================================================================================
  # virtualBiopsy
  # Calculates the position of elements which is possible to do virtual Biopsy
  # ======================================================================================== 
  #' Calculates the position of elements which is possible to do virtual Biopsy
  #' description: This function can be used to calculate the index of elements for virtual Biopsy along a given distance along x,y,z
  #' param: voxelCubes is the voxel space along x
  #' param: nx is the voxel space along x
  #' param: ny is the voxel space along y
  #' param: nz is the voxel space along z
  #' return: A matrix with 1 and 0 which indicates the index for virtual Biopsy
  virtualBiopsyCentroids<-function (voxelCubes,nx,ny,nz){ 
    
    carotaggio.volume <- array(0, dim = dim(voxelCubes))
    # legge ogni singolo elemento della matrice dei voxelCubes ad una distanza dai bordi pari a nx,ny,nz
    for (i in (nx+1):(dim(voxelCubes)[1]-nx))
    {
      for (j in (ny+1):(dim(voxelCubes)[2]-ny))
      {
        for(k in (nz+1):(dim(voxelCubes)[3]-nz))
        {
          if (voxelCubes[i,j,k]!=0)
          {
            # expand grid delle possibili combinazioni dell'intorno, centrate in i,j,k
            combinazioni.poss <- expand.grid(indiceX=seq(i-nx,i+nx),indiceY=seq(j-ny,j+ny),indiceZ=seq(k-nz,k+nz))
            indici.tumore <- as.matrix(combinazioni.poss)
            somma <- 0
            # check degli elementi intorno a i,j,k ed incrementa la variabile somma se diverso da zero
            for(ct in (1:nrow(indici.tumore))) {
              if(indici.tumore[[ct,1]]!=i | indici.tumore[[ct,2]]!=j | indici.tumore[[ct,3]]!=k) {
                if(voxelCubes[indici.tumore[ct,1],indici.tumore[ct,2],indici.tumore[ct,3]]>0) 
                  somma<-somma+1
              }
            }
            # sovrascrive 1 nella posizione i,j,k nella matrice di output se variabile somma è pari al numero di
            # combinazioni calcolate dall'expand grid -1
            if (somma==((nrow(indici.tumore))-1)){
              carotaggio.volume[i,j,k] <- 1
            } 
          }
        }
      }
    }
    return(carotaggio.volume)
  }  
  return(list(SV.getPointPlaneDistance = SV.getPointPlaneDistance,
              SV.get3DPosFromNxNy = SV.get3DPosFromNxNy,
              SV.getPlaneEquationBetween3Points = SV.getPlaneEquationBetween3Points,
              SV.rotateMatrix = SV.rotateMatrix,
              SV.LoadAccordingOSType = SV.LoadAccordingOSType,
              SV.rotateMatrix = SV.rotateMatrix,
              virtualBiopsyCentroids = virtualBiopsyCentroids
              ))  
}



