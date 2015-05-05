#' class for handling logs/warnings/errorss
#' 
#' @description  It handles messages from script to a chosen output (screen, file, etc.)
#' @useDynLib moddicom
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
#   virtualBiopsyCentroids <- function (voxelCubes,nx,ny,nz) {
#     
#     #voxelCubes <- ds.positive[[1]]$voxelCubes[[1]]
#     exam <- array(voxelCubes, dim = dim(voxelCubes)[1]*dim(voxelCubes)[2]*dim(voxelCubes)[3])
#     size1 <- (dim(voxelCubes)[1])
#     size2 <- (dim(voxelCubes)[2])
#     size3 <- (dim(voxelCubes)[3])
#     cmp <- as.matrix(expand.grid(indiceX=seq(i-nx,i+nx),indiceY=seq(j-ny,j+ny),indiceZ=seq(k-nz,k+nz)))
#     expand <- array(cmp, dim = dim(voxelCubes)[1]*dim(voxelCubes)[2])
#     control <- (dim(cmp)[1])
#     lungh <- length(exam)
#     carotaggioVolume <- array(data = c(0), dim = dim(voxelCubes)[1]*dim(voxelCubes)[2]*dim(voxelCubes)[3])
#     
#     #load virtual biopsy C function
#     biopsy <- .C(   "virtualBiopsy", as.integer (exam), as.integer (size1), as.integer (size2), 
#                     as.integer (size3), as.integer(nx), as.integer(ny), as.integer(nz), as.integer(expand),
#                     as.integer (control), as.integer(lungh), as.integer(carotaggioVolume)   )
#     
#     
#     virtual.biopsy <- array (biopsy[11][[1]], dim=c(size1,size2,size3))
#     return(virtual.biopsy)
#     
#   }

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
#' return: a list 
  virtualBiopsy <- function (voxelCubes,nx,ny,nz) {
    
    i<-0;    j<-0;    k<-0

    exam <- array(voxelCubes,dim = dim(voxelCubes)[1]*dim(voxelCubes)[2]*dim(voxelCubes)[3])
    size1 <- (dim(voxelCubes)[1])
    size2 <- (dim(voxelCubes)[2])
    size3 <- (dim(voxelCubes)[3])
    cmp <- as.matrix(expand.grid(indiceX=seq(i-nx,i+nx),indiceY=seq(j-ny,j+ny),indiceZ=seq(k-nz,k+nz)))
    expand <- array(cmp, dim = dim(voxelCubes)[1]*dim(voxelCubes)[2])
    control <- (dim(cmp)[1])
    lungh <- length(exam)
    carotaggioVolume <- array(data = c(0), dim = dim(voxelCubes)[1]*dim(voxelCubes)[2]*dim(voxelCubes)[3])
    
    #load virtual biopsy C function
    biopsy <- .C(   "virtualBiopsy", as.integer (exam), as.integer (size1), as.integer (size2), 
                    as.integer (size3), as.integer(nx), as.integer(ny), as.integer(nz), as.integer(expand),
                    as.integer (control), as.integer(lungh), as.integer(carotaggioVolume)   )
        
    virtual.biopsy <- array (biopsy[11][[1]], dim=c(size1,size2,size3))
    
    # list of of Virtual Biopsy
    indici.carot <- which(virtual.biopsy == 1, arr.ind=T)
    array.devianze2 <- (     rep(  c(0), times=nrow(indici.carot==1)  )  )
    medie2 <- (     rep(  c(0), times=nrow(indici.carot==1)  )  )
    p <- 1;     q <- 1
    lista.grigi <- list()
    if(nrow(indici.carot)<2) return(list("fottiti"="si"))
    
    for(  i in (1:(nrow(indici.carot)))   )
    {      
      comb.poss2 <- expand.grid(    indiceX=   seq(  indici.carot[i,1]-nx,indici.carot[i,1]+nx  ),
                                    indiceY=   seq(  indici.carot[i,2]-ny,indici.carot[i,2]+ny  ),
                                    indiceZ=   seq(  indici.carot[i,3]-nz,indici.carot[i,3]+nz  )         )
      scala.grigi2 <- c()      
      for(ct in (1:nrow(comb.poss2)))
      {
        scala.grigi2[ct] <- voxelCubes [comb.poss2[ct,1],comb.poss2[ct,2],comb.poss2[ct,3]]
      }
      
      lista.grigi[[i]] <- scala.grigi2
      array.devianze2[p] <- sd(scala.grigi2)
      medie2[q] <- mean(scala.grigi2)
      p <- p+1;      q <- q+1
    }
    
    return(list("lista.grigi"=lista.grigi,"devianze"=array.devianze2,"medie"=medie2,"fottiti"="no"))
  }  

  elaboraCarlottaggio<-function( ds.n ,nx=2,ny=2,nz=0) {
    
    UpperBoundDiNormalizzazione<-max(c(obj.n$ROIStats("Urina")$total$max,obj.p$ROIStats("Urina")$total$max))
    
    total<-length(names(ds.n));
    
    
    #medieNorm<-list();
    #devianzeNorm<-list();
    #centroide<-list();
    
    medie<-list();
    devianze<-list();
    
    for(i in names(ds.n) )  {    
      
      print( i )
      piscioPaziente4Tuning<-mean(ds.n[[i]]$voxelCubes[["Urina"]][which(ds.n[[i]]$voxelCubes[["Urina"]]!=0)])
      
      a<-Biopsy(   (ds.n[[ i ]]$voxelCubes$GTV)*(UpperBoundDiNormalizzazione/piscioPaziente4Tuning)    ,nx,ny,nz)
      
      if(a$fottiti!="si")
        medie[[i]]<-a$medie
      devianze[[i]]<-a$devianze  
      
      #     medie[[i]]<-density(a$medie)
      #     devianze[[i]]<-density(a$devianze)
      
      #medieNorm[[i]]<-approx(medie$x,medie$y,n=UpperBoundDiNormalizzazione, xout=seq( from=0 , to=max(UpperBoundDiNormalizzazione) ))
      #devianzeNorm[[i]]<-approx(devianze$x,devianze$y,n=UpperBoundDiNormalizzazione, xout=seq( from=0 , to=max(UpperBoundDiNormalizzazione) )) 
      
      #medieNorm[[i]]$y[which(is.na(medieNorm[[i]]$y))]<-0
      #devianzeNorm[[i]]$y[which(is.na(devianzeNorm[[i]]$y))]<-0
      
      #centro<-COGravity(x=a$medie,y=a$devianze);
      #centroide[[i]]<-c(   centro[1] , centro[3] )
      
    }
    #return( list( "medie"=a$medie, "devianze"=a$devianze,"medieNorm"=medieNorm, "devianzeNorm"=devianzeNorm, "centroide"=centroide  ) )
    return( list( "medie"=medie, "devianze"=devianze  ) )
    
  }  
  lanciaEsempio<-function()
  # list of the available methods of the class
  return(list(SV.getPointPlaneDistance = SV.getPointPlaneDistance,
              SV.get3DPosFromNxNy = SV.get3DPosFromNxNy,
              SV.getPlaneEquationBetween3Points = SV.getPlaneEquationBetween3Points,
              SV.rotateMatrix = SV.rotateMatrix,
              SV.LoadAccordingOSType = SV.LoadAccordingOSType,
              SV.rotateMatrix = SV.rotateMatrix,
              virtualBiopsy  = virtualBiopsy 
              ))  
}



