#' class for handling logs/warnings/errorss
#' 
#' @description  It handles messages from script to a chosen output (screen, file, etc.)
#' @useDynLib moddicom
#' @export
#' @import stringr XML 
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
  SV.getNxNyFrom3D<-function(Px,Py,Pz,oM) {
    Ny<-(oM[2,1]*(Px-oM[1,4])+oM[1,1]*(oM[2,4]-pY)) / (oM[1,2]*oM[2,1]-oM[2,2]*oM[1,1] )
    Nx<-(oM[2,2]*(Px-oM[1,4])+oM[1,2]*(oM[2,4]-Py)) / (oM[1,1]*oM[2,2]-oM[2,1]*oM[1,2] )
    return( list("Nx"=Nx, "Ny"=Ny)  )
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
  SV.rotateMatrix3D<-function( matrice , rotations = 1, inverti=FALSE) {
    for(corsa in seq(1,dim(matrice)[3])) {
      m<-matrice[,,corsa]
      if(rotations == 1 ) m<-t(m[nrow(m):1,])
      if(rotations == 2 ) m<-m[nrow(m):1,ncol(m):1]
      if(rotations == 3 ) m<-t(m)[ncol(m):1,]
      if(inverti==FALSE) nuovaPosizione<-corsa
      else nuovaPosizione<-dim(m)[3]-corsa;
      matrice[,,nuovaPosizione]<-m
    }
    return( matrice )
  }  
  SV.LoadAccordingOSType<-function(library.name) {
    if (Sys.info()["sysname"]=="Windows") 
      return(paste(library.name, ".dll", sep=""))
    if (Sys.info()["sysname"]=="Linux")
      return(paste(library.name, ".so", sep=""))
  }
  SV.rawSurface<-function(voxelMatrix, pSX, pSY, pSZ) {
    nX<-dim(voxelMatrix)[1]
    nY<-dim(voxelMatrix)[2]
    nZ<-dim(voxelMatrix)[3]
    arr<-array(voxelMatrix)
    superficie<-0
    res<-.C("rawSurface",as.double(arr),as.integer(nX), as.integer(nY), as.integer(nZ), as.double(pSX), as.double(pSY), as.double(pSZ), as.double(superficie) );
    return(res[[8]]);    
  } 
  SV.list2XML<-function(  lista, livello=0 , type='1' , outputFileName='') {
    if(!is.list( lista  ))   return (lista);
    objS<-services();
    stringa<-''
    for( i in names(lista) ) {
      if(!is.list(lista[[i]])) {lastLF<-'';} else {lastLF<-'\n';}
      if(type=='1') {
        roiName<-i
        roiName<-str_replace_all(string = roiName,pattern = " ",replacement = "_")
        roiName<-str_replace_all(string = roiName,pattern = ",",replacement = "_")
        stringa<- paste( c(stringa,"\n",rep(" ",livello),"<",roiName,">",objS$SV.list2XML(lista[[i]], livello+1 ),lastLF,rep(" ",livello),"</",roiName,">"    ),collapse = '' )
      }
    }
    if(outputFileName!='' && livello==0) {
      fileConn<-file(outputFileName)
      writeLines(   paste(c("<xml>",stringa,"</xml>"),collapse='')   , fileConn)
      close(fileConn)
      return;
    }
    else return (stringa)
  }
#   SV.trilinearInterpolator<-function(campo,pixelSpacing,newNxNyNzDim) {
#     
#     Nx<-dim(campo)[1];	Ny<-dim(campo)[2];	Nz<-dim(campo)[3]
#     xDim<-pixelSpacing[1];	yDim<-pixelSpacing[2];	zDim<-pixelSpacing[3]
#     newNx<-newNxNyNzDim[1];	newNy<-newNxNyNzDim[2];	newNz<-newNxNyNzDim[3]
#     
#     result<-array(rep(0,newNx*newNx*newNx))	
#     
#     res<-.C("trilinearInterpolator", as.integer(Nx),as.integer(Ny),as.integer(Nz), as.double(xDim),as.double(yDim),as.double(zDim),as.integer(newNx),as.integer(newNy),as.integer(newNz),as.double(campo),as.double(result) );
#     
#     return( array( res[[11]] , dim=c(newNx,newNy,newNz) ) )    
#   }  
  new.SV.trilinearInterpolator<-function(voxelCube = voxelCube,pixelSpacing.new = pixelSpacing.new,pixelSpacing.old = pixelSpacing.old ) {
    
    Nx.old<-dim(voxelCube)[1];	Ny.old<-dim(voxelCube)[2];	Nz.old<-dim(voxelCube)[3]
    xDim.old<-pixelSpacing.old[1];	yDim.old<-pixelSpacing.old[2];	zDim.old<-pixelSpacing.old[3]
    xDim.new<-pixelSpacing.new[1];	yDim.new<-pixelSpacing.new[2];	zDim.new<-pixelSpacing.new[3]
    
    fattoreDiScalaX<-pixelSpacing.old[1]/pixelSpacing.new[1];
    fattoreDiScalaY<-pixelSpacing.old[2]/pixelSpacing.new[2];
    fattoreDiScalaZ<-pixelSpacing.old[3]/pixelSpacing.new[3];
    
    Nx.new<-ceiling(Nx.old * fattoreDiScalaX)
    Ny.new<-ceiling(Ny.old * fattoreDiScalaY)
    Nz.new<-ceiling(Nz.old * fattoreDiScalaZ)
    
    result<-array(rep( 0 , Nx.new * Ny.new * Nz.new ))	
    
    res<-.C("newnewtrilinearInterpolator", 
            as.integer(Nx.old),as.integer(Ny.old),as.integer(Nz.old),
            as.integer(Nx.new),as.integer(Ny.new),as.integer(Nz.new), 
            as.double(pixelSpacing.old[1]),as.double(pixelSpacing.old[2]),as.double(pixelSpacing.old[3]),
            as.double(pixelSpacing.new[1]),as.double(pixelSpacing.new[2]),as.double(pixelSpacing.new[3]),
            as.double(voxelCube),as.double(result) );
    result<-array( res[[14]] , dim=c(Nx.new,Ny.new,Nz.new) )
    return( result )      
  }  
  new.SV.trilinearInterpolator.onGivenPoints<-function(voxelCube ,pixelSpacing.old , newPointCoords.x,newPointCoords.y,newPointCoords.z  ) {
    
    Nx.old<-dim(voxelCube)[1];	Ny.old<-dim(voxelCube)[2];	Nz.old<-dim(voxelCube)[3]
    xDim.old<-pixelSpacing.old[1];	yDim.old<-pixelSpacing.old[2];	zDim.old<-pixelSpacing.old[3]

#     Nx.new<-ceiling(Nx.old * fattoreDiScalaX)
#     Ny.new<-ceiling(Ny.old * fattoreDiScalaY)
#     Nz.new<-ceiling(Nz.old * fattoreDiScalaZ)
    
    punti.asseX<-length(newPointCoords.x);
    punti.asseY<-length(newPointCoords.y);
    punti.asseZ<-length(newPointCoords.z);
    
    result<-array(rep( 0 , punti.asseX * punti.asseY * punti.asseZ ))	

    res<-.C("newnewtrilinearInterpolator_onGivenPoints", 
            as.integer(Nx.old),as.integer(Ny.old),as.integer(Nz.old),
            as.double(pixelSpacing.old[1]),as.double(pixelSpacing.old[2]),as.double(pixelSpacing.old[3]),
            as.double(newPointCoords.x),as.double(newPointCoords.y),as.double(newPointCoords.z),as.integer(punti.asseX),as.integer(punti.asseY),as.integer(punti.asseZ),
            as.double(voxelCube),as.double(result) );
    # result<- res[[14]] 
    result<-array( res[[14]] , dim=c(punti.asseX,punti.asseY,punti.asseZ) )
    return( result )      
  }  
  old_new.SV.trilinearInterpolator.onGivenPoints<-function(voxelCube ,pixelSpacing.old , newPointCoords.x,newPointCoords.y,newPointCoords.z  ) {
    
    Nx.old<-dim(voxelCube)[1];	Ny.old<-dim(voxelCube)[2];	Nz.old<-dim(voxelCube)[3]
    xDim.old<-pixelSpacing.old[1];	yDim.old<-pixelSpacing.old[2];	zDim.old<-pixelSpacing.old[3]
    
    #     Nx.new<-ceiling(Nx.old * fattoreDiScalaX)
    #     Ny.new<-ceiling(Ny.old * fattoreDiScalaY)
    #     Nz.new<-ceiling(Nz.old * fattoreDiScalaZ)
    
    quantiNuoviPunti<-length(newPointCoords.x);
    
    result<-array(rep( 0 , quantiNuoviPunti ))	
    
    res<-.C("newnewtrilinearInterpolator_onGivenPoints", 
            as.integer(Nx.old),as.integer(Ny.old),as.integer(Nz.old),
            as.double(pixelSpacing.old[1]),as.double(pixelSpacing.old[2]),as.double(pixelSpacing.old[3]),
            as.double(newPointCoords.x),as.double(newPointCoords.y),as.double(newPointCoords.z),as.integer(quantiNuoviPunti),
            as.double(voxelCube),as.double(result) );
    result<- res[[12]] 
    return( result )      
  }   
  # ========================================================================================
  # cropCube: crop a voxel cube in order to limit its dimension to the needs
  # ========================================================================================   
  cropCube<-function( bigCube ) {

    matPos<-which(bigCube!=0,arr.ind = T)
    min.x<-min(matPos[,1]);     max.x<-max(matPos[,1])
    min.y<-min(matPos[,2]);     max.y<-max(matPos[,2])
    min.z<-min(matPos[,3]);     max.z<-max(matPos[,3])
    newCube<-bigCube[ min.x:max.x, min.y:max.y , min.z:max.z]
    location<-list( "min.x"=min.x, "max.x"=max.x, "min.y"=min.y, "max.y"=max.y, "min.z"=min.z, "max.z"=max.z  )
    return( list ( "voxelCube"=newCube, "location"=location) )
  }    
  adjCommandLinePar<-function(stringa) {
    if ( Sys.info()["sysname"] == "Windows") {
      nuovaStringa<-str_replace_all(string = stringa,pattern = "'",replacement = "")
    }
    else nuovaStringa<-stringa;
    
    return(nuovaStringa);
  }  
  # ========================================================================================
  # expandCube: expand a cropped voxel cube
  # ========================================================================================     
  expandCube<-function( littleCube,  x.start, y.start, z.start, fe, se, te) {
    
    if(length(dim(littleCube))==3) {
      bigCube<-array(0,dim=c(fe,se,te) )
      for(z in seq(1,dim(littleCube)[3] ) ) {
        for(y in seq(1,dim(littleCube)[2] ) ) {
          for(x in seq(1,dim(littleCube)[1] ) ) {
            bigCube[ x+x.start-1 , y+y.start-1, z+z.start-1  ]<-littleCube[x,y,z]
          }
        }
      }
    }
    if(length(dim(littleCube))==2) {
      bigCube<-array(0,dim=c(fe,se,te) )
      for(y in seq(1,dim(littleCube)[2] ) ) {
        for(x in seq(1,dim(littleCube)[1] ) ) {
          bigCube[ x+x.start-1 , y+y.start-1, 1+z.start-1  ]<-littleCube[x,y]
        }
      }
    }
    if(length(dim(littleCube))==1) stop("FY!");
    return( bigCube )    
  }   
  triangle2mesh <- function(x) {
    v <- list()
    n <- nrow(x$v1)
    nit <- 1:n
    v$vb <- t(cbind(rbind(x$v1,x$v2,x$v3),1))
    v$it <- rbind(nit,nit+n,nit+2*n)
    class(v) <- "mesh3d"
    return(v)
  }  
  getXMLStructureFromDICOMFile<-function(fileName, folderCleanUp = FALSE) {
    fileNameXML<-paste(fileName,".xml")    
    fileNameXML<-str_replace_all(string = fileNameXML , pattern = " .xml",replacement = ".xml")
    pathToStore<-substr(fileName,1,tail(which(strsplit(fileName, '')[[1]]=='/'),1)-1)
    if(!file.exists( fileNameXML ) | folderCleanUp==TRUE) {
      # now check the problem of viewRAY is the 0008,0005 present?
      stringa1<-"dcmdump";
      stringa2<-paste(" ",fileName," | grep '0008,0005'",collapse='')
      options(warn=-1)
      aTMPa<-system2( stringa1,stringa2,stdout=TRUE )
      options(warn=0)
      if(length(aTMPa)==0) {
        stringa1<-"dcmodify";
        stringa2<-paste(" -i '(0008,0005)=ISO_IR 100'  ",fileName,collapse='')
        options(warn=-1)
        system2(stringa1,stringa2,stdout=NULL)
        options(warn=0)        
      }
      
      stringa1<-"dcm2xml";
      stringa2<-paste(" +M  ",fileName,fileNameXML,collapse='')
      options(warn=-1)
      system2(stringa1,stringa2,stdout=NULL)
      options(warn=0)
    }
    # Load the XML file
    doc = xmlInternalTreeParse(fileNameXML)
    return(doc);
  }
#   elaboraCarlottaggio<-function( ds.n ,nx=2,ny=2,nz=0) {
#     
#     UpperBoundDiNormalizzazione<-max(c(obj.n$ROIStats("Urina")$total$max,obj.p$ROIStats("Urina")$total$max))
#     
#     total<-length(names(ds.n));
#     
#     
#     #medieNorm<-list();
#     #devianzeNorm<-list();
#     #centroide<-list();
#     
#     medie<-list();
#     devianze<-list();
#     
#     for(i in names(ds.n) )  {    
#       
#       print( i )
#       piscioPaziente4Tuning<-mean(ds.n[[i]]$voxelCubes[["Urina"]][which(ds.n[[i]]$voxelCubes[["Urina"]]!=0)])
#       
#       a<-Biopsy(   (ds.n[[ i ]]$voxelCubes$GTV)*(UpperBoundDiNormalizzazione/piscioPaziente4Tuning)    ,nx,ny,nz)
#       
#       if(a$fottiti!="si")
#         medie[[i]]<-a$medie
#       devianze[[i]]<-a$devianze  
#       
#       #     medie[[i]]<-density(a$medie)
#       #     devianze[[i]]<-density(a$devianze)
#       
#       #medieNorm[[i]]<-approx(medie$x,medie$y,n=UpperBoundDiNormalizzazione, xout=seq( from=0 , to=max(UpperBoundDiNormalizzazione) ))
#       #devianzeNorm[[i]]<-approx(devianze$x,devianze$y,n=UpperBoundDiNormalizzazione, xout=seq( from=0 , to=max(UpperBoundDiNormalizzazione) )) 
#       
#       #medieNorm[[i]]$y[which(is.na(medieNorm[[i]]$y))]<-0
#       #devianzeNorm[[i]]$y[which(is.na(devianzeNorm[[i]]$y))]<-0
#       
#       #centro<-COGravity(x=a$medie,y=a$devianze);
#       #centroide[[i]]<-c(   centro[1] , centro[3] )
#       
#     }
#     #return( list( "medie"=a$medie, "devianze"=a$devianze,"medieNorm"=medieNorm, "devianzeNorm"=devianzeNorm, "centroide"=centroide  ) )
#     return( list( "medie"=medie, "devianze"=devianze  ) )
#     
#   }  
#lanciaEsempio<-function() 
  # list of the available methods of the class
  return(list(SV.getPointPlaneDistance = SV.getPointPlaneDistance,
              SV.get3DPosFromNxNy = SV.get3DPosFromNxNy,
              SV.getPlaneEquationBetween3Points = SV.getPlaneEquationBetween3Points,
              SV.rotateMatrix = SV.rotateMatrix,
              SV.LoadAccordingOSType = SV.LoadAccordingOSType,
              SV.rotateMatrix = SV.rotateMatrix,
              SV.rotateMatrix3D = SV.rotateMatrix3D,
              SV.rawSurface = SV.rawSurface,
              triangle2mesh = triangle2mesh,
              cropCube = cropCube,
              expandCube = expandCube,
              new.SV.trilinearInterpolator = new.SV.trilinearInterpolator,
              new.SV.trilinearInterpolator.onGivenPoints= new.SV.trilinearInterpolator.onGivenPoints,
              SV.list2XML = SV.list2XML,
              adjCommandLinePar = adjCommandLinePar,
              getXMLStructureFromDICOMFile = getXMLStructureFromDICOMFile
              ))  
}




