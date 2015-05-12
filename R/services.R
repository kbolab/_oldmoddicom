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
  SV.rawSurface<-function(voxelMatrix, pSX, pSY, pSZ) {
    
    if( pSX != pSY ) warning("\n X and Y must have the same pixelSpacing for this implementation");
    
    nX<-dim(voxelMatrix)[1]
    nY<-dim(voxelMatrix)[2]
    nZ<-dim(voxelMatrix)[3]
    arr<-array(voxelMatrix)
    superficie<-0
    res<-.C("rawSurface",as.double(arr),as.integer(nX), as.integer(nY), as.integer(nZ), as.double(pSX), as.double(pSY), as.double(pSZ), as.double(superficie) );
    return(res[[8]]);    
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
#lanciaEsempio<-function() 
  # list of the available methods of the class
  return(list(SV.getPointPlaneDistance = SV.getPointPlaneDistance,
              SV.get3DPosFromNxNy = SV.get3DPosFromNxNy,
              SV.getPlaneEquationBetween3Points = SV.getPlaneEquationBetween3Points,
              SV.rotateMatrix = SV.rotateMatrix,
              SV.LoadAccordingOSType = SV.LoadAccordingOSType,
              SV.rotateMatrix = SV.rotateMatrix,
              SV.rawSurface = SV.rawSurface
              ))  
}




