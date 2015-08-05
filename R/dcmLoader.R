
# obj<-geoLet();
# obj$openDICOMFolder("/progetti/immagini/urinaEasy/Alimenti")


dcmLoader<-function() {
  
  dataStorage<-list();          # Attribute with ALL the data
  dataChache<-list();           # Cache (internal use)
  SOPClassUIDList<-list();      # SOPClassUIDList
  attributeList<-list()
  logObj<-list()               # log handler
  objServ<-list();
  
  #=================================================================================
  # openDICOMFolder
  # Loads a Folder containing one or more DICOM Studies
  #=================================================================================
  #' Open a folder and load the content
  openDICOMFolder<-function(pathToOpen) {
    
    if( attributeList$verbose$lv1 == TRUE ) logObj$sendLog(pathToOpen)
    # get the dcm file type
    SOPClassUIDList<<-getFolderContent(pathToOpen);
    # Load CT/RMN Scans
    if( attributeList$verbose$lv2 == TRUE ) logObj$sendLog("---------------------------------------")
    if( attributeList$verbose$lv2 == TRUE ) logObj$sendLog("Load Images")
    if( attributeList$verbose$lv2 == TRUE ) logObj$sendLog("---------------------------------------")
    imageStructure<-loadCTRMNRDScans(SOPClassUIDList);
    # put images into dataStorage
    dataStorage<<-imageStructure;
    # Load RTStruct Files
    if( attributeList$verbose$lv2 == TRUE ) logObj$sendLog("---------------------------------------")
    if( attributeList$verbose$lv2 == TRUE ) logObj$sendLog("Load RTStruct")
    if( attributeList$verbose$lv2 == TRUE ) logObj$sendLog("---------------------------------------")
    dataStorage[["structures"]]<<-loadRTStructFiles(SOPClassUIDList);   
    # Associate ROI and Images
    associateROIandImageSlices();
    # set the internal attribute indicating the path
    attributeList[["path"]]<<-pathToOpen
  }
  #=================================================================================
  # loadRTStructFiles
  # Loads a DICOM RT Struct (one x folder)
  #=================================================================================
  loadRTStructFiles<-function(SOPClassUIDList) {    
    imageSerie<-list()
    listaPuntiROI<-list()
    
    # loop over the list    
    # even if the assumption is that only one RTStruct is admitted for a CT scan serie
    for(i in names(SOPClassUIDList)) {
      if(  SOPClassUIDList[[i]]$kind=="RTStructureSetStorage") {
        if( attributeList$verbose$lv2 == TRUE ) logObj$sendLog(i) 
        TMP<-getStructuresFromXML( i );
      }
    }
    
    # now let me use some more easy to handle variable names
    matrice2<-TMP$IDROINameAssociation;
    matrice3<-TMP$tableROIPointList;
    
    listaROI<-list()
    
    # for each ROI
    for(i in matrice2[2,]) {
      # get the points
      subMatrix<-matrice3[which(matrice3[,2]==i,arr.ind = TRUE),]
      # if some points exist
      quantiElementiTrovati<--1
      
      if(is.list(subMatrix) & !is.array(subMatrix)) quantiElementiTrovati<-1
      if(is.matrix(subMatrix) & is.array(subMatrix)) quantiElementiTrovati<-dim(subMatrix)[1]
      if(quantiElementiTrovati==-1) stop("Errore inatteso nel caricamento delle slides. Error code:#0001");
      
      if( quantiElementiTrovati >0 ) {
        #      if( length(subMatrix) >0 ) {
        listaROI[[i]]<-list()
        
        # add properly the points to the 'listaROI' structure
        #        for(contatore in seq(1,dim(subMatrix)[1]) ) {
        for(contatore in seq(1,quantiElementiTrovati) ) {
          #          if(i == "Urina" | i =="Vescica") browser();
          
          if( quantiElementiTrovati == 1) {
            ROIPointStringList<-subMatrix[[3]]
            SOPInstance<-subMatrix[[4]]
          }
          else {
            ROIPointStringList<-subMatrix[contatore,3][[1]]
            SOPInstance<-subMatrix[contatore,4][[1]]
          }
          listaCoords<-strsplit(ROIPointStringList,"\\\\");
          listaCoords<-as.numeric(listaCoords[[1]])
          # if a ROI already exists for the slice, well append it to the list
          if( !( SOPInstance  %in% names(listaROI[[i]])  ) )  listaROI[[i]][[   SOPInstance  ]]<-list()          
          listaROI[[i]][[   SOPInstance  ]][[ length(listaROI[[i]][[   SOPInstance  ]])+1  ]]<-matrix(listaCoords,ncol=3,byrow=T)
          # Add the first one as last (close the loop)
          listaROI[[i]][[   SOPInstance  ]][[ length(listaROI[[i]][[   SOPInstance  ]]) ]]<-rbind(listaROI[[i]][[   SOPInstance  ]][[length(listaROI[[i]][[   SOPInstance  ]])]],listaROI[[i]][[   SOPInstance  ]][[length(listaROI[[i]][[   SOPInstance  ]])]][1,])
        }   
      }
    }
    
    return(listaROI); 
  }  
  #=================================================================================
  # associateROIandImageSlices
  # Create the association between ROI and images
  #=================================================================================
  getAttrFromSOPiUID<-function( serUID, sopUID , attrName) {
    return( dataStorage$info[[serUID]][[sopUID]][[attrName]] )
  }
  #=================================================================================
  # associateROIandImageSlices
  # Create the association between ROI and images
  #=================================================================================
  associateROIandImageSlices<-function() {
    # ok now try to find out which is the RMN plane
    lista<-dataStorage$info[[1]]
    planes<-list()
    assignedSCANSandROI<-matrix(NA,ncol=3,nrow=1)
    
    for(i in names(lista)) {
      planes[[i]]<-dataStorage$info[[1]][[i]]$planeEquation
      dataStorage$info[[1]][[i]][["ROIList"]]<<-matrix("",ncol=2);
      
      for(k in names(dataStorage$structures)) {
        for(t in names(dataStorage$structures[[k]])) {
          punto<-dataStorage$structures[[k]][[t]][[1]][1,]
          distanza<-objServ$SV.getPointPlaneDistance(Punto=punto, Piano=planes[[i]])
          
          if( abs(distanza)<0.1 ) {
            if(dataStorage$info[[1]][[i]][["ROIList"]][1]==''  ) {
              dataStorage$info[[1]][[i]][["ROIList"]]<<-matrix(c(k,t),ncol=2)
            }
            else {
              dataStorage$info[[1]][[i]][["ROIList"]]<<-rbind(dataStorage$info[[1]][[i]][["ROIList"]],c(k,t))   #   <<
            }
            assignedSCANSandROI<-rbind(assignedSCANSandROI,c(i,k,t))
          }
        }
      }      
    }  
    
    # now check if it has assigned all the ROIs
    for( i in names(dataStorage$structures) )  {
      lunghezzaAssegnate<-length(which(assignedSCANSandROI[,2]==i))
      lunghezzaCaricate<-length(dataStorage$structures[[i]])
      if(lunghezzaAssegnate!=  lunghezzaCaricate) { 
        
        cat("\n ROI NON CARICATE CORRETTAMENTE: ",i,"\n"); 
      }
    }
  }
  #=================================================================================
  # loadCTRMRDNScans
  # Loads a DICOM CT/MR Scans
  #=================================================================================  
  loadCTRMNRDScans<-function(SOPClassUIDList) {   
    imageSerie<-list()
    objServ<-services()
    
    # loop over the list    
    for(i in names(SOPClassUIDList)) {
      #      if(SOPClassUIDList[[i]]$kind=="RTDoseStorage" | 
      if(  SOPClassUIDList[[i]]$kind=="CTImageStorage" |
             SOPClassUIDList[[i]]$kind=="MRImageStorage" ) {
        
        
        # get the Series number
        seriesInstanceUID<-getDICOMTag(i,"0020,000e")
        # get the Instance number (number of CT in the serie)
        instanceNumber<-getDICOMTag(i,"0020,0013")              
        imageSerie[["info"]][[seriesInstanceUID]][[instanceNumber]]<-list()
        # get the Patient Position
        imageSerie[["info"]][[seriesInstanceUID]][[instanceNumber]][["PatientPosition"]]<-getAttribute(attribute<-"PatientPosition",fileName=i)
        # get the image data
        immagine<-getDICOMTag(i,"7fe0,0010");
        
        # Do I have to rotate the image?
        imageToBeRotated <- -1
        if( imageSerie[["info"]][[seriesInstanceUID]][[instanceNumber]][["PatientPosition"]] == "HFS" ) imageToBeRotated<-1
        if( imageSerie[["info"]][[seriesInstanceUID]][[instanceNumber]][["PatientPosition"]] == "FFS" ) imageToBeRotated<-1
        if ( imageToBeRotated != -1 ) {
          immagine <- objServ$SV.rotateMatrix( immagine, imageToBeRotated )
        }
        
        # now update the structure in memory
        imageSerie[["img"]][[seriesInstanceUID]][[instanceNumber]]<-immagine              
        imageSerie[["info"]][[seriesInstanceUID]][[instanceNumber]][["SOPClassUID"]]<-SOPClassUIDList[[i]]$kind
        imageSerie[["info"]][[seriesInstanceUID]][[instanceNumber]][["SOPInstanceUID"]]<-getDICOMTag(i,"0008,0018")
        imageSerie[["info"]][[seriesInstanceUID]][[instanceNumber]][["FrameOfReferenceUID"]]<-getDICOMTag(i,"0020,0052")
        imageSerie[["info"]][[seriesInstanceUID]][[instanceNumber]][["fileName"]]<-i
        pixelSpacing<-getAttribute(attribute<-"PixelSpacing",fileName=i)
        #pixelSpacing<-getDICOMTag(i,"0028,0030")
        oM<-getAttribute(attribute<-"orientationMatrix",fileName=i)              
        imageSerie[["info"]][[seriesInstanceUID]][[instanceNumber]][["ImagePositionPatient"]]<-getAttribute(attribute<-"ImagePositionPatient",fileName=i)
        imageSerie[["info"]][[seriesInstanceUID]][[instanceNumber]][["ImageOrientationPatient"]]<-getAttribute(attribute<-"ImageOrientationPatient",fileName=i)
        oM[1,1]<-oM[1,1]*pixelSpacing[1]
        oM[2,1]<-oM[2,1]*pixelSpacing[1]
        oM[3,1]<-oM[3,1]*pixelSpacing[1]
        oM[1,2]<-oM[1,2]*pixelSpacing[2]
        oM[2,2]<-oM[2,2]*pixelSpacing[2]
        oM[3,2]<-oM[3,2]*pixelSpacing[2]
        imageSerie[["info"]][[seriesInstanceUID]][[instanceNumber]][["orientationMatrix"]]<-oM
        imageSerie[["info"]][[seriesInstanceUID]][[instanceNumber]][["pixelSpacing"]]<-pixelSpacing
        imageSerie[["info"]][[seriesInstanceUID]][[instanceNumber]][["Rows"]]<-getDICOMTag(i,"0028,0010")
        imageSerie[["info"]][[seriesInstanceUID]][[instanceNumber]][["Columns"]]<-getDICOMTag(i,"0028,0011")
        imageSerie[["info"]][[seriesInstanceUID]][[instanceNumber]][["SliceThickness"]]<-getDICOMTag(i,"0018,0050")
        instanceNumber<-getDICOMTag(i,"0020,0013")
        
        # three points to find out plane equation
        Pa<-c(oM[1,4],oM[2,4],oM[3,4])  
        Pb<-objServ$SV.get3DPosFromNxNy(1000,0,oM)
        Pc<-objServ$SV.get3DPosFromNxNy(0,1000,oM)
        
        abcd<-objServ$SV.getPlaneEquationBetween3Points(Pa,Pb,Pc) 
        
        piano<-matrix(abcd,nrow=1)
        colnames(piano)<-c("a","b","c","d")
        imageSerie[["info"]][[seriesInstanceUID]][[instanceNumber]][["planeEquation"]]<-piano
        if( attributeList$verbose$lv2 == TRUE ) logObj$sendLog(i)
      }
      
      if(  SOPClassUIDList[[i]]$kind=="RTDoseStorage") {
        # get the Series number
        seriesInstanceUID<-getDICOMTag(i,"0020,000e")              
        studyInstanceUID<-getDICOMTag(i,"0020,000d")              
        imagePositionPatient<-getAttribute("ImagePositionPatient",fileName=i)
        ImageOrientationPatient<-getAttribute("ImageOrientationPatient",fileName=i)
        GridFrameOffsetVector<-getAttribute("GridFrameOffsetVector",fileName=i)
        PatientPosition<-getAttribute(attribute<-"PatientPosition",fileName=i)
        Rows<-getDICOMTag(i,"0028,0010");
        Columns<-getDICOMTag(i,"0028,0011")
        
        pixelSpacing<-getAttribute(attribute<-"PixelSpacing",fileName=i)
        imageSerie[["info"]][[seriesInstanceUID]][["1"]][["imagePositionPatient"]]<-imagePositionPatient
        imageSerie[["info"]][[seriesInstanceUID]][["1"]][["ImageOrientationPatient"]]<-ImageOrientationPatient
        imageSerie[["info"]][[seriesInstanceUID]][["1"]][["PatientPosition"]]<-PatientPosition
        imageSerie[["info"]][[seriesInstanceUID]][["1"]][["Rows"]]<-Rows
        imageSerie[["info"]][[seriesInstanceUID]][["1"]][["Columns"]]<-Columns
        imageSerie[["info"]][[seriesInstanceUID]][["1"]][["pixelSpacing"]]<-pixelSpacing
        imageSerie[["info"]][[seriesInstanceUID]][["1"]][["doseType"]]<-getDICOMTag(i,"0020,000d") 
        imageSerie[["info"]][[seriesInstanceUID]][["1"]][["GridFrameOffsetVector"]]<-GridFrameOffsetVector
        imageSerie[["info"]][[seriesInstanceUID]][["1"]][["DoseGridScaling"]]<-getDICOMTag(i,"3004,000e") 
        imageSerie[["info"]][[seriesInstanceUID]][["1"]][["SOPClassUID"]]<-SOPClassUIDList[[i]]$kind
        imageSerie[["info"]][[seriesInstanceUID]][["1"]][["FrameOfReferenceUID"]]<-getDICOMTag(i,"0020,0052")
        immagine<-getDICOMTag(i,"7fe0,0010");
        imageSerie[["dose"]][[seriesInstanceUID]]<-immagine * as.numeric( imageSerie[["info"]][[seriesInstanceUID]][["1"]][["DoseGridScaling"]] )
      }      
    }
    return(imageSerie);
  }
  #=================================================================================
  # getImageFromRAW
  # build a row data from a DICOM file stored on filesystem and load it 
  # into memory (using DCMTK)
  #================================================================================= 
  getStructuresFromXML<-function(fileName) {    
    
    
    fileNameXML<-paste(fileName,".xml")    
    fileNameXML<-str_replace_all(string = fileNameXML , pattern = " .xml",replacement = ".xml")
    pathToStore<-substr(fileName,1,tail(which(strsplit(fileName, '')[[1]]=='/'),1)-1)
    
    # dcmodify -i "(0008,0005)=ISO_IR 100" ./RTSTRUCT999.61977.8337.20150525122026787.dcm
    
    if(!file.exists( fileNameXML )) {
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
    
    # SEQUENCES: the one with the attribute  tag="3006,0020"  and name="StructureSetROISequence" is the one with association NAME<->ID
    n2XML<-getNodeSet(doc,'/file-format/data-set/sequence[@tag="3006,0020" and @name="StructureSetROISequence"]/item')
    # SEQUENCES: now get the true coords
    n3XML<-getNodeSet(doc,'/file-format/data-set/sequence[@tag="3006,0039" and @name="ROIContourSequence"]/item')
    
    # ROI Names
    matrice2<-c()
    for(i in n2XML) {
      ROINumber<-xpathApply(xmlDoc(i),'/item/element[@tag="3006,0022"]',xmlValue)[[1]]
      ROIName<-xpathApply(xmlDoc(i),'/item/element[@tag="3006,0026"]',xmlValue)[[1]]
      matrice2<-rbind(matrice2,c(ROINumber,ROIName))
    }
    matrice2<-t(matrice2)
    
    # ROI Point list
    matrice3<-c()
    for(i in n3XML) {
      ROINumber<-xpathApply(xmlDoc(i),'/item/element[@tag="3006,0084"]',xmlValue)[[1]]
      ROIName<-matrice2[2,which(matrice2[1,]==ROINumber)]
      listaPuntiDaRavanare<-getNodeSet(xmlDoc(i),'/item/sequence/item')
      for(i2 in listaPuntiDaRavanare)   {
#         print("intervenire qui")
#         stop()
        ReferencedSOPInstanceUID<-xpathApply(xmlDoc(i2),'//element[@tag="0008,1155"]',xmlValue)[[1]]
        ROIPointList<-xpathApply(xmlDoc(i2),'/item/element[@tag="3006,0050"]',xmlValue)
        matrice3<-rbind(matrice3,c(ROINumber,ROIName,ROIPointList,ReferencedSOPInstanceUID))
      }
    }
    
    return(list("IDROINameAssociation"=matrice2,"tableROIPointList"=matrice3))
  }
  #=================================================================================
  # getImageFromRAW
  # build a row data from a DICOM file stored on filesystem and load it 
  # into memory (using DCMTK)
  #=================================================================================
  getImageFromRAW<-function(fileName) {    
    objSV<-services()
    fileNameRAW<-paste(fileName,".0.raw")    
    fileNameRAW<-str_replace_all(string = fileNameRAW , pattern = " .0.raw",replacement = ".0.raw")
    pathToStore<-substr(fileName,1,tail(which(strsplit(fileName, '')[[1]]=='/'),1)-1)
    if(!file.exists( fileNameRAW )) {
      stringa1<-"dcmdump";
      stringa2<-paste(" +W  ",pathToStore,fileName,collapse='')
      options(warn=-1)
      system2(stringa1,stringa2,stdout=NULL)
      options(warn=0)
    }
    rowsDICOM<-as.numeric(getDICOMTag(fileName,'0028,0010'))
    columnsDICOM<-as.numeric(getDICOMTag(fileName,'0028,0011'))
    
    bitsAllocated<-as.numeric(getDICOMTag(fileName,'0028,0100'))
    
    if(SOPClassUIDList[[fileName]]$kind!="RTDoseStorage"){
      if(bitsAllocated!=16) stop("16bit pixel are allowed only for non-RTDoseStorage")
      rn<-readBin(con = fileNameRAW, what="integer", size=2, endian="little",n=rowsDICOM*columnsDICOM)    
      rn<-matrix(rn,ncol=columnsDICOM, byrow = TRUE)
    }
    if(SOPClassUIDList[[fileName]]$kind=="RTDoseStorage"){
      if(bitsAllocated==32) {
        if(SOPClassUIDList[[fileName]]$kind!="RTDoseStorage") stop("32bit pixel are allowed only for RTDoseStorage")
        numberOfFrames<-as.numeric(getDICOMTag(fileName,'0028,0008'))
        rn<-readBin(con = fileNameRAW, what="integer", size=4, endian="little",n=rowsDICOM*columnsDICOM*numberOfFrames) 
        # per ora va via come ciclo FOR, poi ci ragioniamo....        
        matRN<-array(0,c(rowsDICOM,columnsDICOM,numberOfFrames))
        ct<-1
        for( z in seq(1,numberOfFrames)) {
          for(x in seq(1,rowsDICOM)) {
            for(y in seq(1,columnsDICOM)) {
              matRN[x,columnsDICOM-y,z]<-rn[ct]
              ct<-ct+1 
            }
          }
        }        
        new_atRN<-array(0,c(columnsDICOM,rowsDICOM,numberOfFrames))
        for(ct in seq(1:dim(matRN)[3]  )) {
          new_atRN[,,ct]<-t(objSV$SV.rotateMatrix( matRN[,,ct], rotations=2 ))
        }
        rn<-new_atRN
      } else  {
        numberOfFrames<-as.numeric(getDICOMTag(fileName,'0028,0008'))
        rn<-readBin(con = fileNameRAW, what="integer", size=2, endian="little",n=rowsDICOM*columnsDICOM*numberOfFrames)
        matRN<-array(0,c(rowsDICOM,columnsDICOM,numberOfFrames))
        ct<-1
        for( z in seq(1,numberOfFrames)) {
          for(x in seq(1,rowsDICOM)) {
            for(y in seq(1,columnsDICOM)) {
              matRN[x,columnsDICOM-y,z]<-rn[ct]
              ct<-ct+1 
            }
          }
        }        
        new_atRN<-array(0,c(columnsDICOM,rowsDICOM,numberOfFrames))
        for(ct in seq(1:dim(matRN)[3]  )) {
          new_atRN[,,ct]<-t(objSV$SV.rotateMatrix( matRN[,,ct], rotations=2 ))
        }
        rn<-new_atRN        
      }
      #stop("RTDose must have 32 bit words!")
    }    
    return(rn)
  }
  getROIList<-function() {
    mat2Ret<-matrix( c(seq(1,length(names(dataStorage$structures))),names(dataStorage$structures)),nrow=2 ,byrow=T )
    return(mat2Ret)
  }
  getROIPointList<-function(ROINumber) {    
    return(dataStorage$structures[ROINumber][[names(dataStorage$structures[ROINumber])]])
  }
  getAssociationTable<-function( tipoTabella, ROIName ) {
    
    if(tipoTabella=="SOPInstance_vs_SliceLocation") {
      matrice<-c()
      involvedCT<-names(dataStorage$structures[[ROIName]]);
      
      for(index in names(dataStorage$info[[1]]) ) {
        matrice<-rbind(matrice,cbind(dataStorage$info[[1]][[index]]$ROIList,index) );
      }
      matrice<-matrice[which(matrice[,1]==ROIName),]
      return(matrice)
    }
  }  
  rotateToAlign<-function(ROIName) {
    tabella1<-getAssociationTable( "SOPInstance_vs_SliceLocation" , ROIName )
    pointList<-getROIPointList( ROIName )
    newPointList<-pointList
    iterazione<-1
    
    for(index in names( pointList )) {
      for(internalIndex in seq(1,length(pointList[[index]]))) {
        IMGsliceInfo<-dataStorage$info[[1]][[ tabella1[ which(tabella1[,2]==index)  ,3  ]  ]]
        sliceLocation<-as.numeric(tabella1[which(tabella1[,2]==index),3])
        DOM<-IMGsliceInfo$orientationMatrix
        
        m<-pointList[[index]][[internalIndex]];

        numeroRighe<-dim(m)[1]
        for(i in seq(1,numeroRighe)) {
          valori<-calcolaNX(m[i,],DOM )
          newPointList[[index]][[internalIndex]][i,1]<-valori$Nx
          newPointList[[index]][[internalIndex]][i,2]<-valori$Ny
          #newPointList[[index]][[internalIndex]][i,3]<-sliceLocation
          newPointList[[index]][[internalIndex]][i,3]<-valori$Nz
        }
      }
      iterazione<-iterazione+1
    }
    return( list("pointList"=newPointList) ) 
  }
  grid.matrixForExtraction<-function( ROIName , SeriesInstanceUIDToBeMapped) {
    
  }
  
  calcolaNX<-function(riga,DOM) {
    Px<-riga[1];    Py<-riga[2];
    a11<-DOM[1,1];    a21<-DOM[2,1];    a31<-DOM[3,1];    a12<-DOM[1,2];
    a22<-DOM[2,2];    a32<-DOM[3,2];    Sx<-DOM[1,4];     Sy<-DOM[2,4];     Sz<-DOM[3,4]; 
    Nx<-(a22*Px-a12*Py-a22*Sx+a12*Sy)/(a11*a22-a21*a12);
    Ny<-(a11*Py-a21*Px+a21*Sx-a11*Sy)/(a22*a11-a21*a12);
    Nz<-Sz;
    return(list("Nx"=Nx,"Ny"=Ny,"Nz"=Nz))
  }
  
  #=================================================================================
  # getAttribute
  # getAttribute: what else?
  #=================================================================================
  
  getAttribute<-function(attribute,seriesInstanceUID="",fileName="") {
    
    if(fileName == "" ) fileName<-getDefaultFileName(seriesInstanceUID)    
    
    # internal DATA
    if(attribute=="dataStorage")  return(dataStorage)  
    # DICOM TAGS
    if(attribute=="PatientName" | attribute=="(0010,0010)")  return(getDICOMTag(fileName,"0010,0010"))  
    if(attribute=="ROIPointList")  return(dataStorage$structures)
    if(attribute=="SOPClassUIDList") return(SOPClassUIDList) 
    if(attribute=="PatientID" | attribute=="(0010,0020)")  return(getDICOMTag(fileName,"0010,0020"))  
    if(attribute=="Rows" | attribute=="(0028,0010)")  return(getDICOMTag(fileName,"0028,0010"))  
    if(attribute=="Columns" | attribute=="(0028,0011)")  return(getDICOMTag(fileName,"0028,0011"))  
    if(attribute=="StudyDate" | attribute=="(0008,0020)")  return(getDICOMTag(fileName,"0008,0020"))  
    if(attribute=="Modality" | attribute=="(0008,0060)")  return(getDICOMTag(fileName,"0008,0060"))  
    if(attribute=="PatientSex" | attribute=="(0010,0040)")  return(getDICOMTag(fileName,"0010,0040"))  
    if(attribute=="SliceThickness" | attribute=="(0018,0050)")  return(getDICOMTag(fileName,"0018,0050"))      
    if(attribute=="SeriesInstanceUID" | attribute=="(0020,000e)")  return(getDICOMTag(fileName,"0020,000e"))      
    if(attribute=="ImagePositionPatient" | attribute=="(0020,0032)")  return(  splittaTAG(getDICOMTag(fileName,"0020,0032"))   )      
    if(attribute=="ImageOrientationPatient" | attribute=="(0020,0037)")  return(   splittaTAG(getDICOMTag(fileName,"0020,0037"))  )
    if(attribute=="PatientPosition" | attribute=="(0018,5100)")  return(  getDICOMTag(fileName,"0018,5100")   )      
    if(attribute=="PixelSpacing" | attribute=="(0028,0030)")  return(splittaTAG(getDICOMTag(fileName,"0028,0030")))
    if(attribute=="PixelSpacing" | attribute=="(0028,0030)")  return(splittaTAG(getDICOMTag(fileName,"0028,0030")))
    if(attribute=="GridFrameOffsetVector" | attribute=="(3004,000c)")  return(splittaTAG(getDICOMTag(fileName,"3004,000c")))
    
    if(attribute=="orientationMatrix")  return( buildOrientationMatrix(fileName)  )
  }
  #=================================================================================
  # splittaTAG
  # internal and stupid function useful to kill artifacts from 
  # a text taken from DCMTK files
  #=================================================================================
  splittaTAG<-function(stringa) {
    return( as.numeric(strsplit(stringa,split = "\\\\")[[1]])   )  
  }
  #=================================================================================
  # buildOrientationMatrix
  # internal and stupid function to build Orientation Matrix
  #=================================================================================
  buildOrientationMatrix<-function(fileName) {
    iPP<-getAttribute(attribute="ImagePositionPatient",fileName=fileName)
    iOP<-getAttribute(attribute="ImageOrientationPatient",fileName=fileName)
    matrice<-matrix(c(iOP[1],iOP[2],iOP[3],0,iOP[4],iOP[5],iOP[6],0,0,0,0,0,iPP[1],iPP[2],iPP[3],1),ncol=4); 
    return(matrice)
  }
  #=================================================================================
  # getDefaultFileName
  # give back the filename of the first CT Scan of the seriesInstanceUID 
  # (if specified)
  #=================================================================================
  getDefaultFileName<-function(seriesInstanceUID="") {
    if(seriesInstanceUID=="") seriesInstanceUID<-names(dataStorage$info)[1]
    primoIndice<-names(dataStorage$info[[seriesInstanceUID]])[1]
    return(dataStorage$info[[seriesInstanceUID]][[primoIndice]]$fileName)
  }
  #=================================================================================
  # NAME: getDICOMTag
  # it queries a DICOM file for a specific tag (short content, no image data)
  #=================================================================================
  getDICOMTag<-function(fileName="",tag=tag) {
    stringa1<-"dcmdump"
    
    # if he want an image, grab it by a raw dump
    if(tag == "7fe0,0010") return( getImageFromRAW(fileName) );
    
    # otherwise go on
    stringa2<-paste(" +L -Un +P '",tag,"'  ",fileName,collapse='')
    options(warn=-1)
    
    a<-as.character(system2(stringa1,stringa2,stdout=TRUE))
    if( length(a) == 0) return ("");
    options(warn=0)
    valore<-"unknown";
    # parse the output according with the different possible fashions
    if( substr(a[1],16,16) == "[" ) {
      valore<-substr(a,which(strsplit(a, '')[[1]]=='[')+1,which(strsplit(a, '')[[1]]==']')-1)
    } 
    if( substr(a[1],16,16) == "(" ) {
      valore<-substr(a,which(strsplit(a, '')[[1]]=='(')[2]+1,which(strsplit(a, '')[[1]]==')')[2]-1)
    }
    if(valore=="unknown") {
      a[length(a)]<-str_trim(substr(a[length(a)],16,which(strsplit(a[length(a)], '')[[1]]=='#')-1))
      if(length(a) == 1) {valore = a;}
      else {valore<-paste(a,collapse="")}
    }
    return(valore)
  }
  #=================================================================================
  # NAME: getFolderContent
  # check a folder and find out the content in terms of DICOM objects
  #=================================================================================
  getFolderContent<-function(pathToOpen="") {
    
    # if no path is give, use the set one
    if(pathToOpen=="") pathToOpen<-attributeList[["path"]];
    
    DCMFilenameArray<-list.files(pathToOpen,"*.dcm")    
    SOPClassUIDList<-list()
    for(i in 1:length(DCMFilenameArray) ) {
      fileNameWithPath<-paste(pathToOpen,"/",DCMFilenameArray[i] , sep="");
      if( substr(fileNameWithPath,nchar(fileNameWithPath)-3,nchar(fileNameWithPath))=='.dcm' ) {
        valore<-getDICOMTag(fileNameWithPath,"0008,0016")
        # do the system call
        SOPClassUIDList[[fileNameWithPath]]<-list();
        SOPClassUIDList[[fileNameWithPath]]$tag<-valore
        SOPClassUIDList[[fileNameWithPath]]$kind<-"Unknown"
        if( SOPClassUIDList[[fileNameWithPath]]$tag == "1.2.840.10008.5.1.4.1.1.2" ) SOPClassUIDList[[fileNameWithPath]]$kind<-"CTImageStorage"
        if( SOPClassUIDList[[fileNameWithPath]]$tag == "1.2.840.10008.5.1.4.1.1.481.2" ) SOPClassUIDList[[fileNameWithPath]]$kind<-"RTDoseStorage"
        if( SOPClassUIDList[[fileNameWithPath]]$tag == "1.2.840.10008.5.1.4.1.1.481.3" ) SOPClassUIDList[[fileNameWithPath]]$kind<-"RTStructureSetStorage"
        if( SOPClassUIDList[[fileNameWithPath]]$tag == "1.2.840.10008.5.1.4.1.1.481.5" ) SOPClassUIDList[[fileNameWithPath]]$kind<-"RTPlanStorage"
        if( SOPClassUIDList[[fileNameWithPath]]$tag == "1.2.840.10008.5.1.4.1.1.4" ) SOPClassUIDList[[fileNameWithPath]]$kind<-"MRImageStorage"
      }
    } 
    return(SOPClassUIDList);
  }  
  setAttribute<-function(attribute, value) {
    if(attribute=="verbose") {
      if(!is.list(value)) return;
      for(i in names(value)) {
        attributeList$verbose[[ i ]] <- value[[i]]
      }
      logObj$setOutput( attributeList$verbose )
      return;  
    }
    attributeList[[ attribute ]]<<-value    
  }
  
  getROIVoxels<-function( Structure = Structure, SeriesInstanceUID = SeriesInstanceUID) {
    SOPClassUID<- dataStorage$info[[SeriesInstanceUID]][[1]]$SOPClassUID
    
    if (SOPClassUID == "CTImageStorage" | SOPClassUID == "MRImageStorage") {
      return (getROIVoxelsFromCTRMN( Structure = Structure, SeriesInstanceUID = SeriesInstanceUID))
    }
    if (SOPClassUID == "RTDoseStorage" ) {
      return (getROIVoxelsFromRTDose( Structure = Structure, SeriesInstanceUID = SeriesInstanceUID))
    }
    print("SOPClassUID not yet prevista");
    stop();
  }
  getROIVoxelsFromRTDose<-function( Structure = Structure, SeriesInstanceUID = SeriesInstanceUID) {
    # define some variables to make more clear the code
    numberOfRows<-as.numeric(dataStorage$info[[SeriesInstanceUID]][[1]]$Rows);
    #numberOfRows<-as.numeric(dataStorage$info[[1]][[1]]$Rows);
    numberOfColumns<-as.numeric(dataStorage$info[[SeriesInstanceUID]][[1]]$Columns);
    #numberOfRows<-as.numeric(dataStorage$info[[1]][[1]]$Columns);
    numberOfSlices<-dim(dataStorage$dose[[1]])[3]
    # initialize the image array with the right dimension
    image.arr<-array( data = -1, dim = c(numberOfRows, numberOfColumns, numberOfSlices ) )    
    
    stop("now I'm here");    
  }
  #   getROIVoxelsFromCTRMN<-function( Structure = Structure, SeriesInstanceUID = SeriesInstanceUID) {
  #     # define some variables to make more clear the code
  #     numberOfRows<-as.numeric(dataStorage$info[[SeriesInstanceUID]][[1]]$Rows);
  #     #numberOfRows<-as.numeric(dataStorage$info[[1]][[1]]$Rows);
  #     numberOfColumns<-as.numeric(dataStorage$info[[SeriesInstanceUID]][[1]]$Columns);    
  #   }
  getROIVoxelsFromCTRMN<-function( Structure = Structure, SeriesInstanceUID = SeriesInstanceUID) {
    objService<-services()    
    # define some variables to make more clear the code
    numberOfRows<-as.numeric(dataStorage$info[[SeriesInstanceUID]][[1]]$Rows);
    #numberOfRows<-as.numeric(dataStorage$info[[1]][[1]]$Rows);
    numberOfColumns<-as.numeric(dataStorage$info[[SeriesInstanceUID]][[1]]$Columns);
    #numberOfRows<-as.numeric(dataStorage$info[[1]][[1]]$Columns);
    numberOfSlices<-length(dataStorage$img[[SeriesInstanceUID]]);
    
    # initialize the image array with the right dimension
    image.arr<-array( data = -1, dim = c(numberOfRows, numberOfColumns, numberOfSlices ) )
    # index values listed as characters: creates empty array of DICOM orientation matrices
    index<-as.character(sort(as.numeric(names( dataStorage$img[[SeriesInstanceUID]]) )))  
    
    # create and fill the vector DOM, of DICOM orientation matrices
    DOM<-c();  nn<-0
    for (n in index) {
      DOM<-c(DOM, dataStorage$info[[SeriesInstanceUID]][[n]]$orientationMatrix[c(1:3,5:7,13:15)])
      nn<-nn+1
      image.arr[,,nn]<-dataStorage$img[[SeriesInstanceUID]][[n]]    
    }  
    # fills the vectors of X and Y coordinates 
    # and other Vectors 'associatedInstanceNumberVect' and 'arrayInstanceNumberWithROI'
    TotalX<- -10000;  TotalY<- -10000;  arrayAssociationROIandSlice<- -10000;
    OriginX<- -10000;   OriginY<- -10000
    associatedInstanceNumberVect<- -10000
    contatoreROI<-1; indiceDOM<-1;
    # for each instance number
    for (n in index) {
      # check if there is a ROI for such slice
      for (m in which(dataStorage$info[[SeriesInstanceUID]][[n]]$ROIList[,1]==Structure)) {
        
        # find the slice and gets the key for accessing at coordinates vectors
        key<-dataStorage$info[[SeriesInstanceUID]][[n]]$ROIList[m,2]
        
        # calculate how many ROIs are co-planar
        numeroROIComplanari<-length(dataStorage$structures[[Structure]][[key]])
        # for each one of them concat the array
        for(indiceROI in seq(1,numeroROIComplanari)) {          
          TotalX<-c(TotalX, dataStorage$structures[[Structure]][[key]][[indiceROI]][,1])
          TotalY<-c(TotalY, dataStorage$structures[[Structure]][[key]][[indiceROI]][,2])
          
          # calculate how many points compose the ROI
          numeroPunti<-length(dataStorage$structures[[Structure]][[key]][[indiceROI]][,1])
          
          # for each point write which is the related Slice in the cube-matrix
          arrayAssociationROIandSlice<-c(arrayAssociationROIandSlice,rep(   which( index == n ) -1  , numeroPunti ))
          
          # Usa OriginX and OriginY as terminator
          TotalX<-c(TotalX, OriginX)
          TotalY<-c(TotalY, OriginY)      
          arrayAssociationROIandSlice<-c(arrayAssociationROIandSlice,OriginX)
          
          contatoreROI<-contatoreROI+1      
        } 
      }      
      indiceDOM<-indiceDOM+1;
    }
    #    browser();
    indiceDOM<-indiceDOM+1;
    indiceDOM<-indiceDOM-1;
    #    browser();
    # ok, call the Wrapper!
    final.array<-NewMultiPointInPolyObl(
      # array of DICOM Orientation Matrices
      DICOMOrientationVector = DOM, 
      # X and Y vector Points
      totalX = TotalX, totalY = TotalY, 
      # association between ROIs and Slices in the 3D Matrix
      arrayAssociationROIandSlice = arrayAssociationROIandSlice,
      # matrices dimensions (rows and columns)
      nX = numberOfColumns, 
      nY = numberOfRows,
      nZ = numberOfSlices
    )
    final.array<-array(data = final.array, dim = c(   numberOfColumns, numberOfRows, numberOfSlices )   )
    # In Example: 
    #
    # > TotalX[0:70]
    # [1] -10000.00     -5.16     -3.28     -1.41      0.47      2.34      4.22      5.12      6.09      7.61      7.97      9.48      9.84     10.96     11.72
    # [16]     12.14     12.93     13.59     13.61     14.19     14.70     14.85     14.23     13.59     13.45     12.68     11.72     11.11      9.84      8.78
    # [31]      7.97      6.09      4.22      2.84      2.34      0.85      0.47      0.14     -1.41     -3.28     -3.92     -5.16     -5.73     -6.56     -7.03
    # [46]     -8.91    -10.26    -10.78    -11.25    -11.86    -12.01    -11.83    -11.29    -10.78    -10.64     -9.88     -8.91     -8.81     -7.51     -7.03
    # [61]     -5.46     -5.16 -10000.00     -5.16     -3.28     -1.41      0.47      2.34      4.22      5.38
    # 
    # > arrayAssociationROIandSlice[0:70]
    # [1] -10000      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8
    # [22]      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8
    # [43]      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      9 -10000
    # [64]      9      9      9      9      9      9      9
    
    # ROTATE THE MATRIX
    
    
    for ( i in seq(1,dim(image.arr)[3] )) {
      #      image.arr[,,i]<-objService$SV.rotateMatrix(image.arr[,,i])
      final.array[,,i]<-t(objService$SV.rotateMatrix(final.array[,,i],rotations=3))
    }
    
    
    #    return(list(TotalX=TotalX, TotalY=TotalY, FullZ=FullZ, Offset=Offset, 
    #                DOM=array(DOM, dim = c(3,3,length(index))), final.array=final.array, masked.images=final.array*image.arr))
    return(list("DOM"=array(DOM, dim = c(3,3,length(index))), "final.array"=final.array, "masked.images"=final.array*image.arr))
  }
  NewMultiPointInPolyObl<-function(DICOMOrientationVector,totalX,totalY,arrayAssociationROIandSlice,nX,nY,nZ ) {  
    
    maxX<-max(totalX)
    minX<-min(totalX[which(totalX>-10000)])
    maxY<-max(totalY)
    minY<-min(totalY[which(totalY>-10000)])
    
    # creates the PIPvector
    PIPvector<-rep.int(x = 0, times = nX * nY * nZ)  
    numberOfPoints<-length(totalX);
    result<-.C("NewMultiPIPObl", 
               as.integer(PIPvector), as.double(totalX), as.double(totalY), as.integer(numberOfPoints), 
               as.integer(nX), as.integer(nY), as.integer(nZ),             
               as.integer(arrayAssociationROIandSlice), 
               as.double(DICOMOrientationVector),as.double(minX),as.double(maxX),as.double(minY),as.double(maxY))  
    
    return(result[[1]])
  }
  #=================================================================================
  # Constructor
  #=================================================================================
  constructor<-function() {
    dataStorage <<- list();   
    dataChache <<- list();
    attributeList<<-list()
    attributeList$verbose<<-list("lv1"=TRUE,"lv2"=TRUE,"lv3"=FALSE,"onScreen"=TRUE,"onFile"=FALSE)
    logObj<<-logHandler()
    #logObj$setOutput( onScreen = attributeList$verbose$onScreen,   onFile = attributeList$verbose$onFile   )
    logObj$setOutput( list("onScreen" = attributeList$verbose$onScreen,   "onFile" = attributeList$verbose$onFile )  )
    logObj$do("clearOutputFile")
    objServ<<-services()
  }
  constructor()
  return(list(openDICOMFolder=openDICOMFolder,getAttribute=getAttribute,
              getDICOMTag=getDICOMTag,getROIList=getROIList,getROIPointList=getROIPointList,
              setAttribute=setAttribute,getFolderContent=getFolderContent,getROIVoxels=getROIVoxels,
              getAssociationTable=getAssociationTable,rotateToAlign=rotateToAlign))
}
