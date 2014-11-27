#' class for handling DICOM files (load, etc.)
#' 
#' @param 
#' @description   Instantiate an object of the class \code{geoLet}.
#' @import stringr XML
#' @export
geoLet<-function() {
  
  dataStorage<-list();          # Attribute with ALL the data
  dataChache<-list();           # Cache (internal use)
  SOPClassUIDList<-list();      # SOPClassUIDList
  attributeList<-list()
  logObj<-list()               # log handler
  
  # ------------------------------------------------
  # openDICOMFolder
  # Loads a Folder containing one or more DICOM Studies
  # ------------------------------------------------
  openDICOMFolder<-function(pathToOpen) {
    if( attributeList$verbose$lv1 == TRUE ) logObj$sendLog(pathToOpen)
    # get the dcm file type
    SOPClassUIDList<<-getFolderContent(pathToOpen);
    # Load CT/RMN Scans
    if( attributeList$verbose$lv2 == TRUE ) logObj$sendLog("---------------------------------------")
    if( attributeList$verbose$lv2 == TRUE ) logObj$sendLog("Load Images")
    if( attributeList$verbose$lv2 == TRUE ) logObj$sendLog("---------------------------------------")
    imageStructure<-loadCTRMNScans(SOPClassUIDList);
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
  # ------------------------------------------------
  # loadRTStructFiles
  # Loads a DICOM RT Struct (one x folder)
  # ------------------------------------------------
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
      if( dim(subMatrix)[1] >0 ) {
        listaROI[[i]]<-list()
        # add properly the points to the 'listaROI' structure
        for(contatore in seq(1,dim(subMatrix)[1]) ) {
          ROIPointStringList<-subMatrix[contatore,3][[1]]
          listaCoords<-strsplit(ROIPointStringList,"\\\\");
          listaCoords<-as.numeric(listaCoords[[1]])
          # if a ROI already exists for the slice, well append it to the list
          if( !( subMatrix[contatore,4][[1]]  %in% names(listaROI[[i]])  ) )  listaROI[[i]][[   subMatrix[contatore,4][[1]]  ]]<-list()          
          listaROI[[i]][[   subMatrix[contatore,4][[1]]  ]][[ length(listaROI[[i]][[   subMatrix[contatore,4][[1]]  ]])+1  ]]<-matrix(listaCoords,ncol=3,byrow=T)
          # Add the first one as last (close the loop)
          listaROI[[i]][[   subMatrix[contatore,4][[1]]  ]][[ length(listaROI[[i]][[   subMatrix[contatore,4][[1]]  ]]) ]]<-rbind(listaROI[[i]][[   subMatrix[contatore,4][[1]]  ]][[length(listaROI[[i]][[   subMatrix[contatore,4][[1]]  ]])]],listaROI[[i]][[   subMatrix[contatore,4][[1]]  ]][[length(listaROI[[i]][[   subMatrix[contatore,4][[1]]  ]])]][1,])
        }   
      }
    }

    return(listaROI); 
  }  
  # ------------------------------------------------
  # getPlanningSeries
  # Returns the image series used to make the plan (in order: RS, RD, RP)
  # ------------------------------------------------  
  getPlanningSeries<-function() {
    return(list("result"="Not yet implemented"))
  }
  # ------------------------------------------------
  # associateROIandImageSlices
  # Create the association between ROI and images
  # ------------------------------------------------  
  associateROIandImageSlices<-function() {
    # ok now try to find out which is the RMN plane
    lista<-dataStorage$info[[1]]
    planes<-list()
    assignedSCANSandROI<-matrix(NA,ncol=3,nrow=1)
    
    for(i in names(lista)) {
      planes[[i]]<-dataStorage$info[[1]][[i]]$planeEquation
      dataStorage$info[[1]][[i]][["ROIList"]]<<-matrix("",ncol=2);    # <<

      for(k in names(dataStorage$structures)) {
        for(t in names(dataStorage$structures[[k]])) {
          punto<-dataStorage$structures[[k]][[t]][[1]][1,]
          distanza<-SV.getPointPlaneDistance(Punto=punto, Piano=planes[[i]])
          
          if( abs(distanza)<0.1 ) {
            if(dataStorage$info[[1]][[i]][["ROIList"]][1]==''  ) {
              dataStorage$info[[1]][[i]][["ROIList"]]<<-matrix(c(k,t),ncol=2)    #   <<
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
  # ------------------------------------------------
  # loadCTRMNScans
  # Loads a DICOM CT/MR Scans
  # ------------------------------------------------  
  loadCTRMNScans<-function(SOPClassUIDList) {   imageSerie<-list()
    # loop over the list    
    for(i in names(SOPClassUIDList)) {
#      if(SOPClassUIDList[[i]]$kind=="RTDoseStorage" | 
       if(  SOPClassUIDList[[i]]$kind=="CTImageStorage" |
           SOPClassUIDList[[i]]$kind=="MRImageStorage" ) {
              # get the Series number
              seriesInstanceUID<-getDICOMTag(i,"0020,000e")
              # get the Instance number (number of CT in the serie)
              instanceNumber<-getDICOMTag(i,"0020,0013")
              # get the image data
              immagine<-getDICOMTag(i,"7fe0,0010");
              # now update the structure in memory
              imageSerie[["img"]][[seriesInstanceUID]][[instanceNumber]]<-immagine
              imageSerie[["info"]][[seriesInstanceUID]][[instanceNumber]]<-list()
              imageSerie[["info"]][[seriesInstanceUID]][[instanceNumber]][["instanceNumber"]]<-SOPClassUIDList[[i]]$kind
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
              
              immagine<-getDICOMTag(i,"7fe0,0010");
              # three points to find out plane equation
              Pa<-c(oM[1,4],oM[2,4],oM[3,4])  
              Pb<-SV.get3DPosFromNxNy(1000,0,oM)
              Pc<-SV.get3DPosFromNxNy(0,1000,oM)

              abcd<-SV.getPlaneEquationBetween3Points(Pa,Pb,Pc) 
                
              piano<-matrix(abcd,nrow=1)
              colnames(piano)<-c("a","b","c","d")
              imageSerie[["info"]][[seriesInstanceUID]][[instanceNumber]][["planeEquation"]]<-piano
              if( attributeList$verbose$lv2 == TRUE ) logObj$sendLog(i)
      }
    }
    return(imageSerie);
  }
  # ------------------------------------------------
  # getImageFromRAW
  # build a row data from a DICOM file stored on filesystem and load it 
  # into memory (using DCMTK)
  # ------------------------------------------------ 
  getStructuresFromXML<-function(fileName) {    
    fileNameXML<-paste(fileName,".xml")    
    fileNameXML<-str_replace_all(string = fileNameXML , pattern = " .xml",replacement = ".xml")
    pathToStore<-substr(fileName,1,tail(which(strsplit(fileName, '')[[1]]=='/'),1)-1)
    
    if(!file.exists( fileNameXML )) {
      stringa1<-"dcm2xml";
      stringa2<-paste(" +M  ",fileName,fileNameXML,collapse='')
      options(warn=-1)
      cat(stringa1,stringa2);
      system2(stringa1,stringa2,stdout=NULL)
      options(warn=0)
    }
  
    # Load the XML file
    doc = xmlInternalTreeParse(fileNameXML)
    
    #n1<-getNodeSet(doc,"/file-format/data-set/element")    
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
        ReferencedSOPInstanceUID<-xpathApply(xmlDoc(i2),'//element[@tag="0008,1155"]',xmlValue)[[1]]
        ROIPointList<-xpathApply(xmlDoc(i2),'/item/element[@tag="3006,0050"]',xmlValue)
        matrice3<-rbind(matrice3,c(ROINumber,ROIName,ROIPointList,ReferencedSOPInstanceUID))
      }
    }
    
    return(list("IDROINameAssociation"=matrice2,"tableROIPointList"=matrice3))
  }
  # ------------------------------------------------
  # getImageFromRAW
  # build a row data from a DICOM file stored on filesystem and load it 
  # into memory (using DCMTK)
  # ------------------------------------------------ 
  getImageFromRAW<-function(fileName) {    
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
    if(bitsAllocated!=16) stop("Bits not allocated in 16 bit word!")
    rn<-readBin(con = fileNameRAW, what="integer", size=2, endian="little",n=rowsDICOM*columnsDICOM)    
    rn<-matrix(rn,ncol=columnsDICOM, byrow = TRUE)
    return(rn)
  }
  getROIList<-function() {
    mat2Ret<-matrix( c(seq(1,length(names(dataStorage$structures))),names(dataStorage$structures)),nrow=2 ,byrow=T )
    return(mat2Ret)
  }
  getROIPointList<-function(ROINumber) {    
    return(dataStorage$structures[ROINumber][[names(dataStorage$structures[ROINumber])]])
  }
# ------------------------------------------------
# getAttribute
# getAttribute: what else?
# ------------------------------------------------ 
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
    if(attribute=="PixelSpacing" | attribute=="(0028,0030)")  return(splittaTAG(getDICOMTag(fileName,"0028,0030")))
    if(attribute=="orientationMatrix")  return( buildOrientationMatrix(fileName)  )
  }
# ------------------------------------------------
# splittaTAG
# internal and stupid function useful to kill artifacts from 
# a text taken from DCMTK files
# ------------------------------------------------ 
  splittaTAG<-function(stringa) {
    return( as.numeric(strsplit(stringa,split = "\\\\")[[1]])   )  
  }
# ------------------------------------------------
# buildOrientationMatrix
# internal and stupid function to build Orientation Matrix
# ------------------------------------------------ 
  buildOrientationMatrix<-function(fileName) {
    iPP<-getAttribute(attribute="ImagePositionPatient",fileName=fileName)
    iOP<-getAttribute(attribute="ImageOrientationPatient",fileName=fileName)
    matrice<-matrix(c(iOP[1],iOP[2],iOP[3],0,iOP[4],iOP[5],iOP[6],0,0,0,0,0,iPP[1],iPP[2],iPP[3],1),ncol=4); 
    return(matrice)
  }
# ------------------------------------------------
# getDefaultFileName
# give back the filename of the first CT Scan of the seriesInstanceUID 
# (if specified)
# ------------------------------------------------ 
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
      return;  
    }
    attributeList[[ attribute ]]<<-value
  }
  constructor<-function() {
    dataStorage <<-list();   
    dataChache<<-list();
    
    attributeList<<-list()
    attributeList$verbose<<-list("lv1"=TRUE,"lv2"=TRUE,"lv3"=FALSE,"onScreen"=TRUE,"onFile"=FALSE)
    logObj<<-logHandler()
    logObj$setOutput( onScreen = attributeList$verbose$onScreen,   onFile = attributeList$verbose$onFile   )
  }
  constructor()
  return(list(openDICOMFolder=openDICOMFolder,getAttribute=getAttribute,getPlanningSeries=getPlanningSeries,
              getDICOMTag=getDICOMTag,getROIList=getROIList,getROIPointList=getROIPointList,setAttribute=setAttribute,
              getFolderContent=getFolderContent))
}

# --------------------------------------------------------------------------------------------------------------------
#  Wrappers
# --------------------------------------------------------------------------------------------------------------------
GL.openDICOMFolder<-function(obj,path) {
  return(obj$openDICOMFolder(path))
}
GL.getAttribute<-function(obj, attrib=c("PatientName","PatientID","dataStorage"),seriesInstanceUID="",fileName="") {
  obj$getAttribute(attribute=attrib,seriesInstanceUID=seriesInstanceUID,fileName=fileName);
}
GL.getROIList<-function(obj) {
  return(obj$getROIList());
}
GL.getROIPointList<-function(obj,ROINameOrID) {
  return( obj$getROIPointList(ROINameOrID) );
}
GL.getPatientName<-function(obj,seriesInstanceUID,fileName="") {
  return( obj$getAttribute(attribute="PatientName",seriesInstanceUID=seriesInstanceUID,fileName=fileName) );
}
GL.getPatientID<-function(obj,seriesInstanceUID="",fileName="") {
  return( obj$getAttribute(attribute="PatientID",seriesInstanceUID=seriesInstanceUID,fileName=fileName) );
}
GL.setAttribute<-function(obj,attributeName,attributeValue) {
  return( obj$setAttribute(attributeName,attributeValue) );
}
GL.getDICOMTag<-function(obj,fileName,tag) {
  return( obj$getDICOMTag(fileName,tag) )  
}
GL.summary<-function() {
  return();
}
GL.plot<-function() {
  return();
}

# --------------------------------------------------------------------------------------------------------------------
#  Examples
# --------------------------------------------------------------------------------------------------------------------

#obj<-geoLet()
#obj$setAttribute("verbose",list("lv1"=TRUE,"lv2"=TRUE))
#obj$openDICOMFolder("/progetti/immagini/Positive/POST/DAgostini")
