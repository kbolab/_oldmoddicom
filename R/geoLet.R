#' class for loading and presenting DICOM data
#' 
#' @description  Instantiate an object of the class \code{geoLet}.This represents just the classname, 
#'               for each instantiated object many methods are available, not in canonical S3 or S4.
#'               
#'               The available methods are:
#'               \itemize{
#'               \item \code{void openDICOMFolder(string pathToOpen)} 
#'               is a method used to open a chosen folder. This method loads all the DICOM objects into
#'               the indicated folder (without recursion) as attribute of the object.
#'               Information can be retrieved using \code{getAttribute} method or 
#'               more specific methods.
#'               \item \code{list/string/array getAttribute(string attribute,string seriesInstanceUID="",string fileName="")} 
#'               Is a method used to get internal attributes. Allowed values fot \code{attributes} are:
#'               
#'               \itemize{
#'                 \item \code{dataStorage} : return a structured list with all the information retrieved from the 
#'                 DICOM object in the folder. More detailed information about the structur of the returned list
#'                 are available at <inserire link>
#'                 \item \code{ROIPointList} : return the ROI point List for all the stored ROIs;
#'                 \item \code{PatientName} : return the content of the (0010,0010) DICOM tag(*);
#'                 \item \code{PatientID} : return the content of the (0010,0020) DICOM tag(*);
#'                 \item \code{Rows} : return the content of the (0028,0010) DICOM tag(*);
#'                 \item \code{Columns} : return the content of the (0028,0011) DICOM tag(*);
#'                 \item \code{StudyDate} : return the content of the (0028,0011) DICOM tag(*);
#'                 \item \code{Modality} : return the content of the (0008,0060) DICOM tag(*);
#'                 \item \code{PatientSex} : return the content of the (0010,0040) DICOM tag(*);
#'                 \item \code{SeriesInstanceUID} : return the content of the (0020,000e) DICOM tag(*);
#'                 \item \code{SliceThickness} : return the content of the (0018,0050) DICOM tag(*);
#'                 \item \code{ImagePositionPatient} : return the content of the (0020,0032) DICOM tag(*);
#'                 \item \code{ImageOrientationPatient} : return the content of the (0020,0037) DICOM tag(*);
#'                 \item \code{PixelSpacing} : return the content of the (0028,0030) DICOM tag(*);
#'               }
#'               
#'               (*) If no \code{seriesInstanceUID} or \code{filename} are provided, it returns the tag found in 
#'               what seems to be the referencing CT or RMN scan. Pay attention to this point: in case of doubt about the
#'               content of the folder this can lead to errors.
#'               \item \code{string getDICOMTag(string fileName="",string tag) }
#'               it allows to retrieve the value of a specific DICOM tag on a specified DICOM filename. If the name is not
#'               indicated the method gets the first DICOM file he can find. 
#'               \item  \code{matrix getROIList() }
#'               returns the list of the available ROIs
#'               \item \code{list getROIPointList(string ROINumber) }
#'               returns the list of the points which define the indicated ROI. 
#'               \item \code{ void setAttribute(string attributeName,string attributeValue) }
#'               it sets the value of the specified attribute. The list of the availeable attributes is the following:
#'                  \itemize{
#'                    \item \code{verbose} : it define the policy adopted to print warnings and logs during computation. It is handled using an \code{errorHandler} object and works with three different levels. A list indicating the behaviour for each level should be provided.
#'                  }
#'               \item \code{list getFolderContent(string pathToOpen="") }
#'               explores the content of the given folder and returns informarion about the stored DICOM objects
#'               \item \code{getPixelSpacing()} it give back the pixelSpaxing along x,y,z of the image scan. Such information can also be obtained gy \code{getAttribute} but it is so common that this explicit method is provided.
#'               \item \code{cacheSave} : it forces to save che most expensive internal memory structure and drop by memory. Until a next \code{cacheLoad} will not be available. It is useful in order to use the memory 'on demand' instead of having everything load all the time.
#'               \item \code{cacheLoad} : it loads the most expensive internal memory structure from filesystem into memory. Pay attention because it override every possible previous content of such list.
#'               \item \code{cacheDrop} : it simply delete the most expensive internal memory structure withouth saving it on filesystem. It is useful if you have previously stored it and you didn't change it. So, in order to use the geoLet object you can simpli load che cache with the \code{cacheLoad} and then drop it without updating the cache on filesystem (in order to save time). Obviously it can be smarty used only if you didn't make any changes...
#'               \item \code{getAlignedStructureAndVoxelCube( double ps.x, double ps.y, double ps.z, string ROIName )} this function resample the image voxel cube and give back the interested ROI rotated on such cube ad adapted according with the new geometry. ps.x,y and z are optional: if not specified the normal pixel Spacing are used.
#'               }
#' @export
#' @import stringr XML misc3d rgl Rvcg
geoLet<-function(ROIVoxelMemoryCache=TRUE,folderCleanUp=FALSE) {
  
  dataStorage<-list();          # Attribute with ALL the data
  dataChache<-list();           # Cache (internal use)
  SOPClassUIDList<-list();      # SOPClassUIDList
  attributeList<-list()
  logObj<-list()               # log handler
  objServ<-list();
  ROIVoxelMemoryCacheArray<-list() # ROIVoxelMemoryCache Array
  mainFrameOfReferenceUID<-NA
  
  #=================================================================================
  # openDICOMFolder
  # Loads a Folder containing one or more DICOM Studies
  #=================================================================================
  # Open a folder and load the content
  openDICOMFolder<-function(pathToOpen, setValidCTRMNSeriesInstanceUID='any', setValidRTPLANSOPInstanceUID='first') {
    
    if( attributeList$verbose$lv1 == TRUE ) logObj$sendLog(pathToOpen)
    
    # get the dcm file type
    SOPClassUIDList<<-getFolderContent(pathToOpen);
    # Load RTPLan (if exists)
    if( attributeList$verbose$lv2 == TRUE ) logObj$sendLog("---------------------------------------")
    if( attributeList$verbose$lv2 == TRUE ) logObj$sendLog("Load RTPlan")
    if( attributeList$verbose$lv2 == TRUE ) logObj$sendLog("---------------------------------------")
    loadRTPlan(SOPClassUIDList,setValidRTPLANSeriesInstanceUID);

    # Load CT/RMN Scans
    if( attributeList$verbose$lv2 == TRUE ) logObj$sendLog("---------------------------------------")
    if( attributeList$verbose$lv2 == TRUE ) logObj$sendLog("Load Images")
    if( attributeList$verbose$lv2 == TRUE ) logObj$sendLog("---------------------------------------")
    loadCTRMNRDScans(SOPClassUIDList);
    
    # Load RTStruct Files
    if( attributeList$verbose$lv2 == TRUE ) logObj$sendLog("---------------------------------------")
    if( attributeList$verbose$lv2 == TRUE ) logObj$sendLog("Load RTStruct")
    if( attributeList$verbose$lv2 == TRUE ) logObj$sendLog("---------------------------------------")
    dataStorage[["structures"]]<<-loadRTStructFiles(SOPClassUIDList);   
    # Associate ROI and Images
    if(is.list(dataStorage[["structures"]])) {
      associateROIandImageSlices();
    }
    # calculate the ImageVoxelCube
    #    if( attributeList$verbose$lv2 == TRUE ) logObj$sendLog("---------------------------------------")
    #    if( attributeList$verbose$lv2 == TRUE ) logObj$sendLog("Creating image Voxel cubes")
    #    if( attributeList$verbose$lv2 == TRUE ) logObj$sendLog("---------------------------------------")
    #    createImageVoxelCube();
    # set the internal attribute indicating the path
    attributeList[["path"]]<<-pathToOpen
    changeDVHROIIDInROINames();
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
    TMP<-list()
    for(i in names(SOPClassUIDList)) {
      if(  SOPClassUIDList[[i]]$kind=="RTStructureSetStorage") {
        if( attributeList$verbose$lv2 == TRUE ) logObj$sendLog(i) 
        TMP[[i]]<-getStructuresFromXML( i );
      }
    }
    if(length(TMP)==0) {
      return( TMP );
    }
   
    # now let me use some more easy to handle variable names
    matrice2<-c(); matrice3<-c(); FORUID.m<-NA;
    for(i in names(TMP)) {
      matrice2<-cbind(matrice2,TMP[[i]]$IDROINameAssociation)
      matrice3<-rbind(matrice3,TMP[[i]]$tableROIPointList)
      # .im
      # Aggiungi le informazioni relative al FrameOfReferenceUID delle ROI caricate
      if(!is.list(dataStorage$info[["structures"]])) dataStorage$info[["structures"]]<-list();
      for( nomeROI in TMP[[i]]$IDROINameAssociation[2,] ) {
        dataStorage$info[["structures"]][[nomeROI]]<<-list();
        dataStorage$info[["structures"]][[nomeROI]]$FrameOfReferenceUID<<-TMP[[i]]$FORUID.m
        dataStorage$info[["structures"]][[nomeROI]]$SeriesInstanceUID<<-TMP[[i]]$RTStructSeriesInstanceUID
      }      
      # .fm
    }
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
        listaROI[[i]]<-list()
        
        # add properly the points to the 'listaROI' structure
        for(contatore in seq(1,quantiElementiTrovati) ) {
          
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
  # getAttrFromSOPiUID
  # Create the association between ROI and images
  #=================================================================================
  getAttrFromSOPiUID<-function( serUID, sopUID , attrName) {
    return( dataStorage$info[[serUID]][[sopUID]][[attrName]] )
  }
  getROIVoxelCube.v2<-function(ROINAME,typeOfVoxels='image',SeriesInstanceUID=NA,pixelSpacing=NA) {
    if(is.na(pixelSpacing)) pixelSpacing<-getPixelSpacing();
    if(is.na(SeriesInstanceUID) & typeOfVoxels=='image')  SeriesInstanceUID<-giveBackImageSeriesInstanceUID();
    
  }
  #=================================================================================
  # associateROIandImageSlices
  # Create the association between ROI and images
  #=================================================================================
  associateROIandImageSlices<-function( relaxCoPlanarity = FALSE) {
    numeroROIDaAssegnare<-0;
    numeroROIAssegnate<-0;
    ROINames<-getROIList();
    for(nomeROI in ROINames[2,]) {
      FrameOfReferenceUID<-dataStorage$info$structures[[nomeROI]]$FrameOfReferenceUID;
      # prendi la seriesInstanceUID con quel FrameOfReferenceUID
      # (per sapere quale serie di immagini avrà associata la ROI)
      
      seriesInstanceUID<-giveBackImageSeriesInstanceUID(FrameOfReferenceUID=FrameOfReferenceUID)
      # bene, ora per ogni ROI scorri i contorni sui vari piani assiali
      listaPuntiROI<-getROIPointList(ROINumber = nomeROI );
      # frulla su ogni assiale
      for( indiciROI in seq(1,length(listaPuntiROI) )) {
        if(length(listaPuntiROI[[indiciROI]])>0) {
          numeroROIDaAssegnare<-numeroROIDaAssegnare+1
          # prendi un punto campione
          sampleROIPoint<-listaPuntiROI[[indiciROI]][[1]][1,]
          # ora cicla sulla serie di immagini identificata come pertinente 
          # l'indice è l'instance number
          for(  imgInstanceNumber  in names(dataStorage$img[[seriesInstanceUID]] ) ) {
            # prendi l'equazione del piano di quella fetta
            planeEquation<-dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]]$planeEquation
            # calcola la distanza fra il piano dell'immagine in esame ed il punto campione della ROI
            distanza<-objServ$SV.getPointPlaneDistance(Punto=sampleROIPoint, Piano=planeEquation)
            # prendi la slice Thickness
            sliceThickness<-as.numeric(dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]]$SliceThickness) 
            # RSTruct SeriesInstanceUID
            STRUCTSeriesInstanceUID<-names(listaPuntiROI)[indiciROI]
            # verifica se la distanza è inferiore ad un dato delta, se sì
            # ciò per decidere se la ROI giace sulla slice
            if( abs(distanza)<0.1 |  ( abs(distanza)<=(sliceThickness/2) & relaxCoPlanarity==TRUE   ) ) {
              if( nomeROI %in% dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]][,1] ) {
                stop("ERRORE: La ROI sembra essere inaspettatamente associata due volte alla stessa slice (problema di posizionamento dei punti nello spazio?)");
              } else {
                # se la tabella manco c'era
                if(length(dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]])==0) {
                  dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]]<<-matrix(0,ncol=4,nrow=1)
                  dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]][1,1]<<-nomeROI
                  dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]][1,2]<<-seriesInstanceUID
                  dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]][1,3]<<-FrameOfReferenceUID
                  dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]][1,4]<<-STRUCTSeriesInstanceUID
                  numeroROIAssegnate<-numeroROIAssegnate+1
                } else {
                  # se invece c'era serve solo aggiungere una riga
                  riga<-nrow(dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]])+1
                  dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]]<<-rbind(dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]],c("","","",""))
                  dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]][riga,1]<<-nomeROI
                  dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]][riga,2]<<-seriesInstanceUID
                  dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]][riga,3]<<-FrameOfReferenceUID
                  dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]][riga,4]<<-STRUCTSeriesInstanceUID                
                  numeroROIAssegnate<-numeroROIAssegnate+1
                }
              }
            }
          }
        }
      }
    }
    if(numeroROIAssegnate!=numeroROIDaAssegnare) {
      cat("ROI non associate correttamente",numeroROIAssegnate,"/",numeroROIDaAssegnare);
    }
  }
  #=================================================================================
  # loadCTRMRDNScans
  # Loads a DICOM CT/MR Scans
  #=================================================================================  
  loadCTRMNRDScans<-function(SOPClassUIDList,setValidCTRMNSeriesInstanceUID='any') {   
    imageSerie<-list()
    objServ<-services()
    
    # loop over the list    
    for(i in names(SOPClassUIDList)) {
      if(  (
        SOPClassUIDList[[i]]$kind=="CTImageStorage" ||
        SOPClassUIDList[[i]]$kind=="MRImageStorage"  ||
        SOPClassUIDList[[i]]$kind=="PositronEmissionTomographyImageStorage")  ) {
        # get the Series number
        #              seriesInstanceUID<-getDICOMTag(i,"0020,000e")
        seriesInstanceUID<-getDICOMTagFromXML(i,"0020,000e")
        i<-i
        if(seriesInstanceUID == setValidCTRMNSeriesInstanceUID || 
           setValidCTRMNSeriesInstanceUID=='any') {
          
          # carica il FrameOfReferenceUID e cerca di capire se è uguale a quello evenutalmente
          # già caricato: se no le geometrie son troppo diverse! (e skippa la serie)
          FrameOfReferenceUID<-getDICOMTag(i,"0020,0052")
          if(  is.na(mainFrameOfReferenceUID) 
               || ( !is.na(mainFrameOfReferenceUID) & FrameOfReferenceUID==mainFrameOfReferenceUID )  ) {
            
            # se son qui significa che il FrameOfReferenceUID  è compatibile e le geometrie sono coerenti
            
            # get the Instance number (number of CT in the serie)
            instanceNumber<-getDICOMTag(i,"0020,0013")              
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]]<<-list()
            # get the Patient Position
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["PatientPosition"]]<<-getAttribute(attribute<-"PatientPosition",fileName=i)
            # get the image data
            immagine<-getDICOMTag(i,"7fe0,0010");
            # Do I have to rotate the image?
            imageToBeRotated <- -1
            if( dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["PatientPosition"]] == "HFS" ) imageToBeRotated<-1
            if( dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["PatientPosition"]] == "FFS" ) imageToBeRotated<-1
            if( dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["PatientPosition"]] == "HFP" ) imageToBeRotated<-1
            
            if( dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["PatientPosition"]] == "FFP" ) imageToBeRotated<-1
            if ( imageToBeRotated != -1 ) {
              immagine <- objServ$SV.rotateMatrix( immagine, imageToBeRotated )
            }
            
            # now update the structure in memory
            dataStorage[["img"]][[seriesInstanceUID]][[instanceNumber]]<<-immagine              
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["SOPClassUID"]]<<-SOPClassUIDList[[i]]$kind
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["SOPInstanceUID"]]<<-getDICOMTag(i,"0008,0018")
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["FrameOfReferenceUID"]]<<-FrameOfReferenceUID
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["fileName"]]<<-i
            pixelSpacing<-getAttribute(attribute<-"PixelSpacing",fileName=i)
            #pixSeriesInstanceUIDelSpacing<-getDICOMTag(i,"0028,0030")
            oM<-getAttribute(attribute<-"orientationMatrix",fileName=i)              
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["ImagePositionPatient"]]<<-getAttribute(attribute<-"ImagePositionPatient",fileName=i)
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["ImageOrientationPatient"]]<<-getAttribute(attribute<-"ImageOrientationPatient",fileName=i)
            oM[1,1]<-oM[1,1]*pixelSpacing[1]
            oM[2,1]<-oM[2,1]*pixelSpacing[1]
            oM[3,1]<-oM[3,1]*pixelSpacing[1]
            oM[1,2]<-oM[1,2]*pixelSpacing[2]
            oM[2,2]<-oM[2,2]*pixelSpacing[2]
            oM[3,2]<-oM[3,2]*pixelSpacing[2]
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["orientationMatrix"]]<<-oM
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["pixelSpacing"]]<<-pixelSpacing
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["Rows"]]<<-getDICOMTag(i,"0028,0010")
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["Columns"]]<<-getDICOMTag(i,"0028,0011")
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["SliceThickness"]]<<-getDICOMTag(i,"0018,0050")
            if(SOPClassUIDList[[i]]$kind=="MRImageStorage") {
              dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["RepetitionTime"]]<<-getDICOMTag(i,"0018,0080")
              dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["EchoTime"]]<<-getDICOMTag(i,"0018,0081")
            }
            instanceNumber<-getDICOMTag(i,"0020,0013")
            
            # three points to find out plane equation
            Pa<-c(oM[1,4],oM[2,4],oM[3,4])  
            Pb<-objServ$SV.get3DPosFromNxNy(1000,0,oM)
            Pc<-objServ$SV.get3DPosFromNxNy(0,1000,oM)
            
            abcd<-objServ$SV.getPlaneEquationBetween3Points(Pa,Pb,Pc) 
            piano<-matrix(abcd,nrow=1)
            colnames(piano)<-c("a","b","c","d")
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["planeEquation"]]<<-piano
            if( attributeList$verbose$lv2 == TRUE ) logObj$sendLog(i)
          }
        }
      }
      
      if(  SOPClassUIDList[[i]]$kind=="RTDoseStorage") {
        
        # get the mainFrameOfReferenceUID
        FrameOfReferenceUID<-getDICOMTag(i,"0020,0052")
        
        # carica il FrameOfReferenceUID e cerca di capire se è uguale a quello evenutalmente
        # già caricato: se no le geometrie son troppo diverse! (e skippa la serie)        
        if(  is.na(mainFrameOfReferenceUID) 
             || ( !is.na(mainFrameOfReferenceUID) & FrameOfReferenceUID==mainFrameOfReferenceUID )  ) {
            
          # carica l'XML ed estraine i campi
          datiDiDose<-getDoseFromXML( i )

          # verifica prima di tutto se la dose caricata fa riferimento all'RTPLAN
          # presente in memoria (se c'è!)
          if( is.null(dataStorage$info$plan$SOPInstanceUID) || ( 
            !is.null(dataStorage$info$plan$SOPInstanceUID) &
            dataStorage$info$plan$SOPInstanceUID == datiDiDose$ReferencedRTPlanSequence_ReferencedSOPInstanceUID )) {
          
            if(length(dataStorage[["info"]])==0) dataStorage[["info"]]<<-list();
            if(length(dataStorage[["info"]][["doses"]])==0) dataStorage[["info"]][["doses"]]<<-list();
            SOPInstanceUID<-datiDiDose$SOPInstanceUID;
            dataStorage[["info"]][["doses"]][[SOPInstanceUID]]<<-list();
            dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["imagePositionPatient"]]<<-splittaTAG(datiDiDose$ImagePositionPatient)
            dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["seriesInstanceUID"]]<<-datiDiDose$SeriesInstanceUID
            dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["FrameOfReferenceUID"]]<<-FrameOfReferenceUID
            dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["ImageOrientationPatient"]]<<-splittaTAG(datiDiDose$ImageOrientationPatient)
            dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["Rows"]]<<-datiDiDose$Rows
            dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["Columns"]]<<-datiDiDose$Columns
            dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["doseType"]]<<-datiDiDose$DoseType
            dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["GridFrameOffsetVector"]]<<-splittaTAG(datiDiDose$GridFrameOffsetVector)
            dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["DoseGridScaling"]]<<-datiDiDose$DoseGridScaling
            dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["ReferencedSOPInstanceUID"]]<<-datiDiDose$ReferencedRTPlanSequence_ReferencedSOPInstanceUID 
            dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["SOPClassUID"]]<<-SOPClassUIDList[[i]]$kind
            dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["pixelSpacing"]]<<-splittaTAG(datiDiDose$PixelSpacing)
            dataStorage[["info"]][["DVHs"]][[SOPInstanceUID]][["DVHFromFile"]]<<-datiDiDose$DVHList

            # estrai l'immagine
            immagine<-getDICOMTag(i,"7fe0,0010");
            if(length(dataStorage[["dose"]])==0) dataStorage[["dose"]]<<-list();
            dataStorage[["dose"]][[SOPInstanceUID]]<<-immagine * as.numeric( dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["DoseGridScaling"]] )
            
            # controlli 'formali' di congruità
            # Verifica che le DOSI abbiano tutte la stessa geometria
            if(length(dataStorage[["info"]][["doses"]])>1) {
              if(!all(dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["imagePositionPatient"]] %in% dataStorage[["info"]][["doses"]][[1]][["imagePositionPatient"]]))
                stop("ERRORE: due dosi differiscono per 'imagePositionPatient");
              if(dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["FrameOfReferenceUID"]]!=dataStorage[["info"]][["doses"]][[1]][["FrameOfReferenceUID"]])
                stop("ERRORE: due dosi differiscono per 'FrameOfReferenceUID");  
              if(!all(dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["ImageOrientationPatient"]] %in% dataStorage[["info"]][["doses"]][[1]][["ImageOrientationPatient"]]))
                stop("ERRORE: due dosi differiscono per 'ImageOrientationPatient");  
              if(dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["Rows"]]!=dataStorage[["info"]][["doses"]][[1]][["Rows"]])
                stop("ERRORE: due dosi differiscono per 'Rows");      
              if(dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["Columns"]]!=dataStorage[["info"]][["doses"]][[1]][["Columns"]])
                stop("ERRORE: due dosi differiscono per 'Columns");    
              if(dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["pixelSpacing"]][1]!=dataStorage[["info"]][["doses"]][[1]][["pixelSpacing"]][1])
                stop("ERRORE: due dosi differiscono per 'pixelSpacing' x-dim");    
              if(dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["pixelSpacing"]][2]!=dataStorage[["info"]][["doses"]][[1]][["pixelSpacing"]][2])
                stop("ERRORE: due dosi differiscono per 'pixelSpacing' y-dim");            
               if(!all(dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["GridFrameOffsetVector"]] %in% dataStorage[["info"]][["doses"]][[1]][["GridFrameOffsetVector"]]))
                 stop("ERRORE: due dosi differiscono per 'GridFrameOffsetVector");        
              if(dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["doseType"]]!="PHYSICAL")
                stop("ERRORE: la dose non è 'PHYSICAL' (attributo 'DoseType')");
            }
          }
        }
      }      
    }
  }
  #=================================================================================
  # loadRTPlan
  # load the RTPlan
  #=================================================================================    
  loadRTPlan<-function(SOPClassUIDList,setValidRTPLANSOPInstanceUID='first') {
    setValidRTPLANSOPInstanceUID<-'first'
    # loop over the list    
    for(i in names(SOPClassUIDList)) {
      if(  SOPClassUIDList[[i]]$kind=="RTPlanStorage" ) {
        
        strutture<-getPlanFromXML(fileName = i);
        
        # prendi la SOPInstanceUID e la seriesInstanceUID del piano
        SOPInstanceUID<-getDICOMTag(i,"0008,0018")
        seriesInstanceUID<-getDICOMTag(i,"0020,000e") 
        toConsider<-FALSE;

        if(SOPInstanceUID == setValidRTPLANSOPInstanceUID) toConsider<-TRUE;
        if(setValidRTPLANSOPInstanceUID == 'first' & length( dataStorage$info$plan )==0) toConsider<-TRUE;
        # se è il primo che trovi (ed era stato dichiarato di prendere il primo)
        # o se la SOPInstanceUID coincide con quanto cercato, allora caricalo
        if(toConsider==TRUE) {
          
          seriesInstanceUID<-strutture$SeriesInstanceUID
          RTPlanGeometry<-strutture$RTPlanGeometry
          ReferencedStructureSetSequence<-strutture$ReferencedStructureSetSequence_ReferencedSOPInstanceUID
          FrameOfReferenceUID<-strutture$FrameOfReferenceUID
          DoseReferenceSequence_DoseReferenceUID<-strutture$DoseReferenceSequence_DoseReferenceUID
          
          if(RTPlanGeometry!="PATIENT") stop("RTPlan has RTPlanGeometry set to something different than 'PATIENT'");
          
          if(length(dataStorage)==0) dataStorage<<-list();
          if(length(dataStorage$info)==0) dataStorage$info<<-list();
          if(length(dataStorage$info$plan)==0) dataStorage$info$plan<<-list();
    
          dataStorage$info$plan$SOPInstanceUID<<-SOPInstanceUID
          dataStorage$info$plan$seriesInstanceUID<<-seriesInstanceUID
          dataStorage$info$plan$RTPlanGeometry<<-RTPlanGeometry
          dataStorage$info$plan$FrameOfReferenceUID<<-FrameOfReferenceUID
          dataStorage$info$plan$ReferencedStructureSetSequence<<-ReferencedStructureSetSequence
          dataStorage$info$plan$DoseReferenceSequence_DoseReferenceUID<<-DoseReferenceSequence_DoseReferenceUID
          
          # setta il mainFrameOfReferenceUID  
          # ( importante: ce ne può essere solo uno)
          mainFrameOfReferenceUID<<-FrameOfReferenceUID
          
          # print the line on terminal
          if( attributeList$verbose$lv2 == TRUE ) logObj$sendLog(i)
        }
      }
    }
  }  
  #=================================================================================
  # createImageVoxelCube
  # create the imageVoxelCube for the current obj and for the image stored
  #=================================================================================   
  createImageVoxelCube<-function() {
    # ge the series Instance UID of the images
    seriesInstanceUID<-giveBackImageSeriesInstanceUID();
    # order them according with the Instance Number
    listaSeqImages<-as.character(sort(as.numeric(names( dataStorage$img[[seriesInstanceUID]]) )))
    # get dimensional data
    Rows<-dataStorage$info[[seriesInstanceUID]][[1]]$Rows
    Columns<-dataStorage$info[[seriesInstanceUID]][[1]]$Columns
    Slices<-length(listaSeqImages)
    
    cubone<-array(data = 0,dim = c(Columns,Rows,Slices))
    numSlice<-1
    # add the slices and build the cube
    for(i in listaSeqImages) {
      cubone[,,numSlice]<-dataStorage$img[[seriesInstanceUID]][[as.character(i)]]
      numSlice<-numSlice+1
    }
    return(cubone)
  }  
  #=================================================================================
  # getImageVoxelCube
  # give back the greyLevel voxel cube. If no ps.x/y/z are specified it gives back 
  # the voxelCube of the original dimensions, otherwise it gives back the interpolated
  # voxelCube according to the wished pixelSpacing along x,y or z
  #=================================================================================     
  getImageVoxelCube<-function( ps.x=NA, ps.y=NA, ps.z=NA) {
    objS<-services();
    # prendi il cubone
    voxelCube<-createImageVoxelCube()
    
    # se non  server interpolare
    if(is.na(ps.x) && is.na(ps.y) && is.na(ps.z) ) return(voxelCube)
    
    # se invece serve interpolare: prendi i pixelSpacing lungo la X, la Y e la Z (slice thickness)
    #     seriesInstanceUID<-giveBackImageSeriesInstanceUID();   
    #     oldPixelSpacing<-c();
    #     oldPixelSpacing[1]<-as.numeric(dataStorage$info[[seriesInstanceUID]][[1]]$pixelSpacing[1])
    #     oldPixelSpacing[2]<-as.numeric(dataStorage$info[[seriesInstanceUID]][[1]]$pixelSpacing[2])
    #     oldPixelSpacing[3]<-as.numeric(dataStorage$info[[seriesInstanceUID]][[1]]$SliceThickness)
    oldPixelSpacing<-getPixelSpacing();
    
    if(is.na(ps.x))  ps.x <- oldPixelSpacing[1];
    if(is.na(ps.y))  ps.y <- oldPixelSpacing[2];
    if(is.na(ps.z))  ps.z <- oldPixelSpacing[3];
    
    voxelCube<-objS$new.SV.trilinearInterpolator(
      voxelCube = voxelCube,
      pixelSpacing.new = c(ps.x,ps.y,ps.z),
      pixelSpacing.old = oldPixelSpacing )    
    
    return( voxelCube )
  }  
  #=================================================================================
  # getDoseVoxelCube
  # restituisce il voxel cube della dose
  #=================================================================================    
  getDoseVoxelCube<-function() {
    doseVC<-dataStorage$dose[[1]]
    doseInfo<-dataStorage$info$doses[[1]]
    return( list("voxelCube"=doseVC,"info"=doseInfo)  )
  }
  #=================================================================================
  # getPixelSpacing
  # una funzione specifica per il pixelSpacing, visto quanto è usato!
  #================================================================================= 
  getPixelSpacing<-function() {
    seriesInstanceUID<-giveBackImageSeriesInstanceUID();   
    ps<-c();
    ps[1]<-as.numeric(dataStorage$info[[seriesInstanceUID]][[1]]$pixelSpacing[1])
    ps[2]<-as.numeric(dataStorage$info[[seriesInstanceUID]][[1]]$pixelSpacing[2])
    ps[3]<-as.numeric(dataStorage$info[[seriesInstanceUID]][[1]]$SliceThickness)
    return(ps);
  }  
  getAlignedStructureAndVoxelCube<-function(  ps.x=NA, ps.y=NA, ps.z=NA, ROIName ) {
    voxelCube<-getImageVoxelCube( ps.x = ps.x, ps.y = ps.y, ps.z = ps.z ) 
    ROI<-rotateToAlign(ROIName = ROIName)
    old.ps<-getPixelSpacing();
    delta.x<-old.ps[1] / ps.x;
    delta.y<-old.ps[2] / ps.y;
    delta.z<-old.ps[3] / ps.z;
    for(i in names(ROI$pointList)) {
      for( ct in seq(1,length(ROI$pointList[[i]]) ) ) {
        ROI$pointList[[i]][[ct]][,3]<-ROI$pointList[[i]][[ct]][,3] * delta.z;
        ROI$pointList[[i]][[ct]][,2]<-ROI$pointList[[i]][[ct]][,2] * delta.y;
        ROI$pointList[[i]][[ct]][,1]<-ROI$pointList[[i]][[ct]][,1] * delta.x;
      }
    }
    return( list( "voxelCube"=voxelCube, "ROI"=ROI$pointList )  )
  }
  #=================================================================================
  # giveBackImageSeriesInstanceUID
  # from dataStorage it gives back the SOPInstanceUID of the series which has a 
  # SOPClassUID as 'CTImageStorage' or 'MRImageStorage' or 'PET'. 
  # Frame
  #=================================================================================
  giveBackImageSeriesInstanceUID<-function(FrameOfReferenceUID=NA, ROIName=NA) {
    # se è per FrameOfReferenceUID

    if(!is.na(FrameOfReferenceUID)) return(searchIMGSeriesForFrameOfReferenceUID(FrameOfReferenceUID=FrameOfReferenceUID));
    # se le vuoi tutte...
    list.index<-''
    for(whichIdentifier in seq(1,length(dataStorage$info) )) {
      if(names(dataStorage$info)[whichIdentifier]!="structures" &
         names(dataStorage$info)[whichIdentifier]!="plan"  &
         names(dataStorage$info)[whichIdentifier]!="DVHs"  &
         names(dataStorage$info)[whichIdentifier]!="doses" ) {
        if(dataStorage$info[[whichIdentifier]][[1]]$SOPClassUID == 'CTImageStorage' ||
           dataStorage$info[[whichIdentifier]][[1]]$SOPClassUID == 'MRImageStorage' ||
           dataStorage$info[[whichIdentifier]][[1]]$SOPClassUID == 'PositronEmissionTomographyImageStorage'
        ) {
          #list.index<-whichIdentifier
            list.index<-names(dataStorage$info)[whichIdentifier]
          }
      }
    }
    return(list.index)
  }
  #=================================================================================
  # getGeometricalInformationOfImage
  # give back pixelspacing and other little stuff about the CT/MR
  #=================================================================================  
  getGeometricalInformationOfImage<-function() {
    serieInstanceUID<-giveBackImageSeriesInstanceUID();
    ddd<-dataStorage$info[[ serieInstanceUID ]][[1]]
    return(list(
      "PatientPosition"=ddd$PatientPosition,
      "SOPClassUID"=ddd$MRImageStorage,
      "pixelSpacing"=ddd$pixelSpacing,
      "ImagePositionPatient"=ddd$ImagePositionPatient,
      "Rows"=ddd$Rows,
      "Columns"=ddd$Columns,
      "SliceThickness"=ddd$SliceThickness,
      "supposedNumberOfSlices"=length(dataStorage$info[[ serieInstanceUID ]]),
      "randomSliceImageOrientationPatient"=ddd$ImageOrientationPatient,
      "randomSlicePlaneEquation"=ddd$planeEquation
    ))
  }
  #=================================================================================
  # getPlanFromXML
  # get the interesting tag via XML instead of via normal dump (more robust)
  #=================================================================================   
  getDoseFromXML<-function(fileName) { 
    obj.S<-services();
    
    # build the XML file and get the XML structure
    doc<-obj.S$getXMLStructureFromDICOMFile(fileName = fileName, folderCleanUp = folderCleanUp)
    
    ImagePositionPatient<-xpathApply(doc,'/file-format/data-set/element[@tag="0020,0032" and @name="ImagePositionPatient"]',xmlValue)[[1]]
    ImageOrientationPatient<-xpathApply(doc,'/file-format/data-set/element[@tag="0020,0037" and @name="ImageOrientationPatient"]',xmlValue)[[1]]
    Rows<-xpathApply(doc,'/file-format/data-set/element[@tag="0028,0010" and @name="Rows"]',xmlValue)[[1]]
    Columns<-xpathApply(doc,'/file-format/data-set/element[@tag="0028,0011" and @name="Columns"]',xmlValue)[[1]]
    PixelSpacing<-xpathApply(doc,'/file-format/data-set/element[@tag="0028,0030" and @name="PixelSpacing"]',xmlValue)[[1]]
    PixelRepresentation<-xpathApply(doc,'/file-format/data-set/element[@tag="0028,0103" and @name="PixelRepresentation"]',xmlValue)[[1]]
    SOPInstanceUID<-xpathApply(doc,'/file-format/data-set/element[@tag="0008,0018" and @name="SOPInstanceUID"]',xmlValue)[[1]]
    SOPInstanceUID<-xpathApply(doc,'/file-format/data-set/element[@tag="0008,0018" and @name="SOPInstanceUID"]',xmlValue)[[1]]
    SeriesInstanceUID<-xpathApply(doc,'/file-format/data-set/element[@tag="0020,000e" and @name="SeriesInstanceUID"]',xmlValue)[[1]]
    GridFrameOffsetVector<-xpathApply(doc,'/file-format/data-set/element[@tag="3004,000c" and @name="GridFrameOffsetVector"]',xmlValue)[[1]]
    DoseUnits<-xpathApply(doc,'/file-format/data-set/element[@tag="3004,0002" and @name="DoseUnits"]',xmlValue)[[1]]
    DoseType<-xpathApply(doc,'/file-format/data-set/element[@tag="3004,0004" and @name="DoseType"]',xmlValue)[[1]]
    DoseGridScaling<-xpathApply(doc,'/file-format/data-set/element[@tag="3004,000e" and @name="DoseGridScaling"]',xmlValue)[[1]]
    ReferencedRTPlanSequence_ReferencedSOPInstanceUID<-xpathApply(doc,'/file-format/data-set/sequence[@tag="300c,0002" and @name="ReferencedRTPlanSequence"]//element[@tag="0008,1155" and @name="ReferencedSOPInstanceUID"]',xmlValue)[[1]]    
    
    # now look for some DVHs
    #a<-xpathApply(doc,'/file-format/data-set/sequence[@tag="3004,0050" and @name="DVHSequence"]',xmlValue)
    a<-getNodeSet(doc,'/file-format/data-set/sequence[@tag="3004,0050" and @name="DVHSequence"]/item')

    DVHList<-list();
    ct<-1
    if(length(a)>0) {
      for(ct in seq(1,length(a))) {

        ReferencedROINumber<-xpathApply(a[[ct]],'//item/sequence[@tag="3004,0060"]//element[@tag="3006,0084"]',xmlValue)[[ct]]
        
        ROIName<-as.character(ReferencedROINumber);

        DVHList[[ROIName]]<-list();
        
        DVHList[[ROIName]][["DVHType"]]<-xpathApply(a[[ct]],'//item/element[@tag="3004,0001"]',xmlValue)[[ct]]
        DVHList[[ROIName]][["DoseUnits"]]<-xpathApply(a[[ct]],'//item/element[@tag="3004,0002"]',xmlValue)[[ct]]
        DVHList[[ROIName]][["DoseType"]]<-xpathApply(a[[ct]],'//item/element[@tag="3004,0004"]',xmlValue)[[ct]]
        DVHList[[ROIName]][["DVHDoseScaling"]]<-xpathApply(a[[ct]],'//item/element[@tag="3004,0052"]',xmlValue)[[ct]]
        DVHList[[ROIName]][["DVHVolumeUnits"]]<-xpathApply(a[[ct]],'//item/element[@tag="3004,0054"]',xmlValue)[[ct]]
        DVHList[[ROIName]][["DVHNumberOfBins"]]<-xpathApply(a[[ct]],'//item/element[@tag="3004,0056"]',xmlValue)[[ct]]
        DVHList[[ROIName]][["DVHMeanDose"]]<-xpathApply(a[[ct]],'//item/element[@tag="3004,0074"]',xmlValue)[[ct]]
        DVHList[[ROIName]][["DVHMaximumDose"]]<-xpathApply(a[[ct]],'//item/element[@tag="3004,0072"]',xmlValue)[[ct]]
        DVHList[[ROIName]][["ReferencedROINumber"]]<-xpathApply(a[[ct]],'//item/sequence[@tag="3004,0060"]//element[@tag="3006,0084"]',xmlValue)[[ct]]
        
        dvhString<-xpathApply(a[[ct]],'//item/element[@tag="3004,0058"]',xmlValue)[[ct]];
        dvhArr<-strsplit(dvhString,"\\\\");
        DVHList[[ROIName]][["DVHData.volume"]]<-as.numeric(dvhArr[[1]][seq(2,length(dvhArr[[1]]),by=2 )])
        DVHList[[ROIName]][["DVHData.dose"]]<-cumsum(  as.numeric(dvhArr[[1]][seq(1,length(dvhArr[[1]]),by=2 )])  )
        dvh.type<-''
        final.matrix<-as.matrix(  cbind( DVHList[[ROIName]][["DVHData.volume"]],DVHList[[ROIName]][["DVHData.dose"]] ) )
        if(DVHList[[ROIName]][["DVHType"]]=="CUMULATIVE") dvh.type<-"cumulative";
        if(DVHList[[ROIName]][["DVHType"]]=="DIFFERENTIAL") dvh.type<-"differential";
        if(dvh.type=='') stop("ERROR: type of DVH not yet supported");
        if(DVHList[[ROIName]][["DoseUnits"]]!='GY') stop("ERROR: only 'GY' are supported as DoseUnit");
        if(DVHList[[ROIName]][["DoseType"]]!='PHYSICAL') stop("ERROR: only DoseType 'physical' is supported");
        if(DVHList[[ROIName]][["DVHDoseScaling"]]!='1') stop("ERROR: only DVHDoseScaling equal to 1 is supported");
        if(DVHList[[ROIName]][["DVHVolumeUnits"]]!='CM3') stop("ERROR: only DVHVolumeUnits equal to 'CM3' is supported");
        final.matrix<-as.matrix(cbind(DVHList[[ROIName]][["DVHData.dose"]],DVHList[[ROIName]][["DVHData.volume"]]))
        DVHObj<-new("dvhmatrix", dvh=final.matrix, dvh.type=dvh.type, vol.distr='absolute', volume=final.matrix[1,1])
        DVHList[[ROIName]][["DVHObj"]]<-DVHObj
      }
    }    

    return(list(
      "ImagePositionPatient"=ImagePositionPatient,"ImageOrientationPatient"=ImageOrientationPatient,
      "Rows"=Rows,"Columns"=Columns,"PixelSpacing"=PixelSpacing,"PixelRepresentation"=PixelRepresentation,
      "SOPInstanceUID"=SOPInstanceUID,"SeriesInstanceUID"=SeriesInstanceUID,"GridFrameOffsetVector"=GridFrameOffsetVector,
      "DoseUnits"=DoseUnits,"DoseType"=DoseType,"DoseGridScaling"=DoseGridScaling,
      "ReferencedRTPlanSequence_ReferencedSOPInstanceUID"=ReferencedRTPlanSequence_ReferencedSOPInstanceUID,
      "DVHList"=DVHList
    ))
  }
  #=================================================================================
  # changeDVHROIIDInROINames
  # at the end of the computation, change the ROIId in the ROINames. This cannot be 
  # done at the beginning because DicomRT object can be loaded in any order
  #=================================================================================   
  changeDVHROIIDInROINames<-function() {
    matriceNomiROI<-getROIList();
    if(!is.list(dataStorage$info$DVHs)) return;
    if(!is.list(dataStorage$info$DVHs[[1]])) return;
    for( SOPInstanceUID in names(dataStorage$info$DVHs) ) {
      listaDaSostituire<-list();
      for(indColonna in seq(1,dim(matriceNomiROI)[2] )) {
        for(numericID in names( dataStorage$info$DVHs[[SOPInstanceUID]]$DVHFromFile )) {
          nuovoNome<-matriceNomiROI[2,which(matriceNomiROI[1,]==numericID,arr.ind = T)]
          if(length(nuovoNome)>0) {
            listaDaSostituire[[nuovoNome]]<-dataStorage$info$DVHs[[SOPInstanceUID]]$DVHFromFile[[numericID]]
          }
        }
      }
      dataStorage$info$DVHs[[SOPInstanceUID]]$DVHFromFile<<-list();
      dataStorage$info$DVHs[[SOPInstanceUID]]$DVHFromFile<<-listaDaSostituire
    }
  }
  #=================================================================================
  # getPlanFromXML
  # get the interesting tag via XML instead of via normal dump (more robust)
  #=================================================================================   
  getPlanFromXML<-function(fileName) { 
    obj.S<-services();
    
    # Load the XML file
    doc<-obj.S$getXMLStructureFromDICOMFile(fileName = fileName, folderCleanUp = folderCleanUp)
    
    SOPInstanceUID<-xpathApply(doc,'//element[@tag="0008,0018" and @name="SOPInstanceUID"]',xmlValue)[[1]]
    StudyDate<-xpathApply(doc,'//element[@tag="0008,0020" and @name="StudyDate"]',xmlValue)[[1]]
    SeriesInstanceUID<-xpathApply(doc,'//element[@tag="0020,000e" and @name="SeriesInstanceUID"]',xmlValue)[[1]]
    FrameOfReferenceUID<-xpathApply(doc,'//element[@tag="0020,0052" and @name="FrameOfReferenceUID"]',xmlValue)[[1]]
    RTPlanGeometry<-xpathApply(doc,'//element[@tag="300a,000c" and @name="RTPlanGeometry"]',xmlValue)[[1]]
    ReferencedStructureSetSequence_ReferencedSOPInstanceUID<-xpathApply(doc,'//sequence[@tag="300c,0060" and @name="ReferencedStructureSetSequence"]//element[@tag="0008,1155" and @name="ReferencedSOPInstanceUID"]',xmlValue)[[1]]
    DoseReferenceSequence_DoseReferenceUID<-xpathApply(doc,'//sequence[@tag="300a,0010" and @name="DoseReferenceSequence"]//element[@tag="300a,0013" and @name="DoseReferenceUID"]',xmlValue)
    
    return( list("SOPInstanceUID"=SOPInstanceUID,"StudyDate"=StudyDate,"SeriesInstanceUID"=SeriesInstanceUID,
                 "FrameOfReferenceUID"=FrameOfReferenceUID,"RTPlanGeometry"=RTPlanGeometry,
                 "ReferencedStructureSetSequence_ReferencedSOPInstanceUID"=ReferencedStructureSetSequence_ReferencedSOPInstanceUID,
                 "DoseReferenceSequence_DoseReferenceUID" = DoseReferenceSequence_DoseReferenceUID))
  }
  getXYZFromNxNyNzOfDoseVolume<-function(Nx, Ny, Nz, startFromZero=FALSE) {
    objS<-services();
    if(length(dataStorage$info$doses)==1) {seriesInstanceUID<-names(dataStorage$info$doses)[[1]];}
    else {stop(" caso non ancora previsto ( più volumi di dose a questo livello #njj9)");} 
    pixelSpacing<-dataStorage$info$doses[[seriesInstanceUID]]$pixelSpacing
    imagePositionPatient<-dataStorage$info$doses[[seriesInstanceUID]]$imagePositionPatient
    ImageOrientationPatient<-dataStorage$info$doses[[seriesInstanceUID]]$ImageOrientationPatient
    mat<-array(dataStorage$info$doses[[seriesInstanceUID]]$ImageOrientationPatient,dim=c(3,2))
    mat<-rbind(mat,c(0,0));    mat<-cbind(mat,c(0,0,0,0));
    mat[,1]<-mat[,1]*pixelSpacing[1];   mat[,2]<-mat[,2]*pixelSpacing[2];
    mat<-cbind(mat,c(0,0,0,1)); mat[1:3,4]<-imagePositionPatient;
    mat[3,4]<-dataStorage$info$doses[[seriesInstanceUID]]$GridFrameOffsetVector[Nz]
    punto<-objS$SV.get3DPosFromNxNy(Nx,Ny,mat);
    return(punto);
  }
  getXYZFromNxNyNzOfImageVolume<-function(Nx, Ny, Nz, startFromZero=FALSE) {
    objS<-services();
    
    serieInstanceUID<-giveBackImageSeriesInstanceUID();
    istanceNumbers<-as.character(sort(as.numeric(names(dataStorage$info[[serieInstanceUID]]))))
    if(startFromZero==FALSE) ct<-1;
    if(startFromZero==TRUE) ct<-0;
    for(slice in istanceNumbers) {
      if(ct==Nz) {
        imagePosition<-dataStorage$info[[serieInstanceUID]][[slice]]$ImagePositionPatient;
        ImageOrientationPatient<-dataStorage$info[[serieInstanceUID]][[slice]]$ImageOrientationPatient;
        pixelSpacing<-dataStorage$info[[serieInstanceUID]][[slice]]$pixelSpacing;
        orientationMatrix<-dataStorage$info[[serieInstanceUID]][[slice]]$orientationMatrix;
        punto<-objS$SV.get3DPosFromNxNy(Nx,Ny,orientationMatrix);
#        print(punto);
        return(punto);
      }
      ct<-ct+1;
    }
    stop();
    print(c(Nx,Ny,Nz))
  }
  #=================================================================================
  # getStructuresFromXML
  # get the interesting tag via XML instead of via normal dump (more robust)
  #================================================================================= 
  getStructuresFromXML<-function(fileName) {    
    obj.S<-services();
    massimo<-0

    # Load the XML file
    doc<-obj.S$getXMLStructureFromDICOMFile(fileName = fileName, folderCleanUp = folderCleanUp)
    
    # prima di tutto controlla che la FrameOfReferenceUID sia la stessa OVUNQUE e che punti
    # ad una serie di immagini ESISTENTE!
    # E' un chiodo ma .... ragionevole, almeno per ora
    RTStructSeriesInstanceUID<-xpathApply(doc,'//element[@tag="0020,000e" and @name="SeriesInstanceUID"]',xmlValue)[[1]]
    FORUID.m<-xpathApply(doc,'//element[@tag="0020,0052" and @name="FrameOfReferenceUID"]',xmlValue)[[1]]
    FORUID.d<-xpathApply(doc,'//element[@tag="3006,0024" and @name="ReferencedFrameOfReferenceUID"]',xmlValue)
    for(FORUID.d_index in seq(1,length(FORUID.d))) {
      if( FORUID.d[[ FORUID.d_index ]] !=  FORUID.m ) {
        stop("ERRORE: FrameOfReferenceUID non allineati nel file RTStruct")
      }
    }
    if(is.na(giveBackImageSeriesInstanceUID(FrameOfReferenceUID = FORUID.m))) {
      stop("ERRORE: FrameOfReferenceUID del file RTStruct non associato a nessuna serie di immagini caricata")
    }

    # SEQUENCES: the one with the attribute  tag="3006,0020"  and name="StructureSetROISequence" 
    # is the one with association NAME<->ID
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
    massimo<-0
    
    matrice3<-c()
    for(i in n3XML) {
      ROINumber<-xpathApply(xmlDoc(i),'/item/element[@tag="3006,0084"]',xmlValue)[[1]]
      ROIName<-matrice2[2,which(matrice2[1,]==ROINumber)]
      listaPuntiDaRavanare<-getNodeSet(xmlDoc(i),'/item/sequence/item')
      
      for(i2 in listaPuntiDaRavanare)   {
        # ReferencedSOPInstanceUID<-xpathApply(xmlDoc(i2),'//element[@tag="0008,1155"]',xmlValue)[[1]]
        
        ROIPointList<-xpathApply(xmlDoc(i2),'/item/element[@tag="3006,0050"]',xmlValue)        
        splittedROIPointList<-as.numeric(strsplit(ROIPointList[[1]],split = "\\\\")[[1]])
        fPoint.x<-splittedROIPointList[1]
        fPoint.y<-splittedROIPointList[2]
        fPoint.z<-splittedROIPointList[3]
        
        assegnato<-FALSE
        ReferencedSOPInstanceUID<-''
        #        list.index<-''
        list.index<-giveBackImageSeriesInstanceUID()

        if(list.index=='') stop("list.index is empty!")
        for(slice.index in seq(1,length(dataStorage$info[[list.index]]))) {
          distanza<-obj.S$SV.getPointPlaneDistance(c(fPoint.x,fPoint.y,fPoint.z),dataStorage$info[[list.index]][[slice.index]]$planeEquation)
          if( abs(distanza)<0.2 ) {
            ReferencedSOPInstanceUID<-dataStorage$info[[list.index]][[slice.index]]$SOPInstanceUID
            assegnato<-TRUE
          }
        }
        if( assegnato == FALSE ) cat("WARNING: ROI ",ROIName,": the point (",fPoint.x,",",fPoint.y,",",fPoint.z,") has no image slice!\n"); 
        matrice3<-rbind(matrice3,c(ROINumber,ROIName,ROIPointList,ReferencedSOPInstanceUID))
      }
    }
    return(list("IDROINameAssociation"=matrice2,"tableROIPointList"=matrice3,"FORUID.m"=FORUID.m,"RTStructSeriesInstanceUID"=RTStructSeriesInstanceUID))
  }
  #=================================================================================
  # searchIMGSeriesForFrameOfReferenceUID
  # restituisce la seriesInstanceUID della serie di immagini che fa riferimento 
  # ad un dato FrameOfReferenceUID
  #=================================================================================
  searchIMGSeriesForFrameOfReferenceUID<-function(FrameOfReferenceUID) {
    listaSeriesInstanceUID<-names(dataStorage$img)
    for( index4Info in names(dataStorage$img)) {
      if(index4Info!="structures") {
        if(dataStorage$info[[index4Info]][[1]]$FrameOfReferenceUID == FrameOfReferenceUID) {
          return(index4Info);
        }
      }
    }
    return(NA)
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
    if(!file.exists( fileNameRAW )  | folderCleanUp==TRUE) {
      stringa1<-"dcmdump";
      #fileNameFS<-str_replace_all(string = fileName,pattern = "/",replacement = "\\\\")
      if ( Sys.info()["sysname"] == "Windows") {
        fileNameFS<-chartr("\\","/",fileName);
        stringa2<-chartr("/","\\\\",stringa2)
      }
      else fileNameFS<-fileName;
      stringa2<-paste(" +W  ",pathToStore,fileNameFS,collapse='')
      options(warn=-1)
      system2(stringa1,stringa2,stdout=NULL)
      options(warn=0)
    }
    rowsDICOM<-as.numeric(getDICOMTag(fileName,'0028,0010'))
    columnsDICOM<-as.numeric(getDICOMTag(fileName,'0028,0011'))
    bitsAllocated<-as.numeric(getDICOMTag(fileName,'0028,0100'))
    if(SOPClassUIDList[[fileName]]$kind!="RTDoseStorage"){
      if(bitsAllocated!=16) stop("16bit pixel are allowed only for non-RTDoseStorage")
      #fileNameRAWFS<-str_replace_all(string = fileNameRAW,pattern = "/",replacement = "\\\\")
      if ( Sys.info()["sysname"] == "Windows") {
        fileNameRAWFS<-chartr("\\","/",fileNameRAW);
        fileNameRAWFS<-chartr("/","\\\\",fileNameRAWFS);
        }
      else fileNameRAWFS<-fileNameRAW;
      
      rn<-readBin(con = fileNameRAWFS, what="integer", size=2, endian="little",n=rowsDICOM*columnsDICOM)    
      rn<-matrix(rn,ncol=columnsDICOM, byrow = TRUE)
      #rn<-matrix(rn,ncol=columnsDICOM)
    }
    if(SOPClassUIDList[[fileName]]$kind=="RTDoseStorage"){
      if(bitsAllocated==32) {
        if(SOPClassUIDList[[fileName]]$kind!="RTDoseStorage") stop("32bit pixel are allowed only for RTDoseStorage")
        numberOfFrames<-as.numeric(getDICOMTag(fileName,'0028,0008'))
        # fileNameRAWFS<-str_replace_all(string = fileNameRAW,pattern = "/",replacement = "\\\\")
        if ( Sys.info()["sysname"] == "Windows") {
          fileNameRAWFS<-chartr("\\","/",fileNameRAW);
          fileNameRAWFS<-chartr("/","\\\\",fileNameRAWFS);
        }
        else fileNameRAWFS<-fileNameRAW;
        
        rn<-readBin(con = fileNameRAWFS, what="integer", size=4, endian="little",n=rowsDICOM*columnsDICOM*numberOfFrames) 
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
        #fileNameRAWFS<-str_replace_all(string = fileNameRAW,pattern = "/",replacement = "\\\\")
        if ( Sys.info()["sysname"] == "Windows") {
          fileNameRAWFS<-chartr("\\","/",fileNameRAW);
          fileNameRAWFS<-chartr("/","\\\\",fileNameRAWFS)
        }
        else fileNameRAWFS<-fileNameRAW;
        rn<-readBin(con = fileNameRAWFS, what="integer", size=2, endian="little",n=rowsDICOM*columnsDICOM*numberOfFrames)
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
  #=================================================================================
  # getROIList
  # restituisce la lista delle ROI
  #=================================================================================  
  getROIList<-function() {
    mat2Ret<-matrix( c(seq(1,length(names(dataStorage$structures))),names(dataStorage$structures)),nrow=2 ,byrow=T )
    return(mat2Ret)
  }
  #=================================================================================
  # getROIPointList
  # restituisce la lista delle coordinate dei vertici di una data ROI
  #=================================================================================    
  getROIPointList<-function(ROINumber) {    
    return(dataStorage$structures[ROINumber][[names(dataStorage$structures[ROINumber])]])
  }
  #=================================================================================
  # getAttribute
  # getAttribute: what else?
  #=================================================================================
  getAttribute<-function(attribute,seriesInstanceUID="",fileName="") {
    if(fileName == "" ) {
      primoIndice<-names(dataStorage$info[[  giveBackImageSeriesInstanceUID()  ]])[1]
      fileName<-dataStorage$info[[seriesInstanceUID]][[primoIndice]]$fileName
    }
    
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
    if(attribute=="ROIVoxelMemoryCache") return (ROIVoxelMemoryCache);
    if(attribute=="DVHsFromFile") return (giveBackDVHsFromFile());
    
    if(attribute=="orientationMatrix")  return( buildOrientationMatrix(fileName)  )
    return(getDICOMTag(fileName,attribute))  
  }
  #=================================================================================
  # giveBackDVHsFromFile
  # give back the DVH loaded from File
  #=================================================================================  
  giveBackDVHsFromFile<-function() {
    listaDVH<-list();
    for(SOPInstanceUID in names(dataStorage$info$DVHs)) {
      for(ROINumber in names(dataStorage$info$DVHs[[SOPInstanceUID]]$DVHFromFile)) {
        if(is.list(listaDVH[[ROINumber]])) stop("ERRORE: lo stesso ROINUMBER compare in diversi RTDose");
        listaDVH[[ROINumber]]<-dataStorage$info$DVHs[[SOPInstanceUID]]$DVHFromFile[[ROINumber]]
      }
    }
    return(listaDVH);
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
  # NAME: getDICOMTag
  # it queries a DICOM file for a specific tag (short content, no image data)
  #=================================================================================
  getDICOMTag<-function(fileName="",tag=tag) {
    stringa1<-"dcmdump"
    
    # if he want an image, grab it by a raw dump
    if(tag == "7fe0,0010") return( getImageFromRAW(fileName) );
    
    # otherwise go on
    stringa2<-objServ$adjCommandLinePar(paste(" +L -Un +P '",tag,"'  ",fileName,collapse=''))
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
  # NAME: getDICOMTagFromXML
  # estrae una tag dall'XML
  #=================================================================================  
  getDICOMTagFromXML<-function(fileName="",tag=tag) {    
    obj.S<-services();
    massimo<-0
    fileNameXML<-paste(fileName,".xml")    
    fileNameXML<-str_replace_all(string = fileNameXML , pattern = " .xml",replacement = ".xml")
    pathToStore<-substr(fileName,1,tail(which(strsplit(fileName, '')[[1]]=='/'),1)-1)
    
    # dcmodify -i "(0008,0005)=ISO_IR 100" ./RTSTRUCT999.61977.8337.20150525122026787.dcm
    if(!file.exists( fileNameXML )  | folderCleanUp==TRUE ) {
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
    # SEQUENCES: the one with the attribute  tag="3006,0020"  and name="StructureSetROISequence" 
    # is the one with association NAME<->ID
    stringaQuery<-paste(c('/file-format/data-set/element[@tag="',tag,'"]'),collapse='');
    valore<-xpathApply(doc,stringaQuery,xmlValue);
    if(length(valore)==2) stop("ERRORE: due valori sono inattesi");
    valore<-valore[[1]]
    return(valore);
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
        if( SOPClassUIDList[[fileNameWithPath]]$tag == "1.2.840.10008.5.1.4.1.1.128" ) SOPClassUIDList[[fileNameWithPath]]$kind<-"PositronEmissionTomographyImageStorage"
      }
    } 
    return(SOPClassUIDList);
  }  
  #=================================================================================
  # NAME: setAttribute
  #=================================================================================  
  setAttribute<-function(attribute, value) {
    if(attribute=="verbose") {
      if(!is.list(value)) return;
      for(i in names(value)) {
        attributeList$verbose[[ i ]] <- value[[i]]
      }
      logObj$setOutput( attributeList$verbose )
      return;  
    }
    if(attribute=="cacheDir") {
      attributeList$virtualMemory$path<<-value;
      return;
    }
    attributeList[[ attribute ]]<<-value    
  }
  #=================================================================================
  # NAME: getDoseInformation
  # give back the information about DOSES in the space
  #=================================================================================    
  getDoseInformation<-function() {
    # da modificare (ora aggancia solo il primo!)
    informazioniDose<-ds$info$doses[[1]]
    return(informazioniDose)
  }  
  getSpatialInformationOfVoxelCube<-function(SeriesInstanceUID) {
    zpos<-c()
    objS<-services();
    for(istanceNumber in dataStorage$info[[SeriesInstanceUID]]) {
      m.DOM<-dataStorage$info[[SeriesInstanceUID]][[istanceNumber]]$orientationMatrix
      zpos<-c(zpos,m.DOM[3,4])
    }
    min.z<-min(zpos)
    max.z<-max(zpos)
    m.DOM.min<-m.DOM; m.DOM.min[3,4]<-min.z;
    m.DOM.max<-m.DOM; m.DOM.max[3,4]<-max.z;
    xyz.min<-objS$SV.get3DPosFromNxNy(0,0,m.DOM.min);
    xyz.max<-objS$SV.get3DPosFromNxNy(0,0,m.DOM.max);
    return( list(
      "xyz.min"=xyz.min,"xyz.max"=xyz.max
    )  );
  }
  #=================================================================================
  # NAME: getROIVoxels
  # restituisce i voxel interni ad una data ROI
  #=================================================================================    
  getROIVoxels<-function( Structure = Structure ) {
    objS<-services();

    # cerca nella cache di memoria, se attivata, la presenza della ROI già estratta
    if(ROIVoxelMemoryCache==TRUE  & !is.null(ROIVoxelMemoryCacheArray[[Structure]]) ) {
      return(ROIVoxelMemoryCacheArray[[Structure]]);
    }
    
    # try to find out which Series is the CT/MR serie
    SeriesInstanceUID<-giveBackImageSeriesInstanceUID()
    if(SeriesInstanceUID == '' ) stop("ERROR: missing CT/MR series")
    res<-getROIVoxelsFromCTRMN( Structure = Structure, SeriesInstanceUID = SeriesInstanceUID)
    
    croppedRes<-list()
    croppedRes$DOM<-res$DOM
    croppedRes$geometricalInformationOfImages<-res$geometricalInformationOfImages
    croppedRes$masked.images<-objS$cropCube( bigCube = res$masked.images)
    croppedRes$masked.images$location$fe<-dim(res$masked.images)[1]
    croppedRes$masked.images$location$se<-dim(res$masked.images)[2]
    croppedRes$masked.images$location$te<-dim(res$masked.images)[3]
    croppedRes$geometricalInformationOfImages$koc<-"littleCube"   
    class(croppedRes)<-"geoLetStructureVoxelList"
    
    # se la cache in memoria è attiva, salvane una copia
    if(ROIVoxelMemoryCache==TRUE) ROIVoxelMemoryCacheArray[[Structure]]<<-croppedRes;
    return( croppedRes )
  }
  old_getROIVoxels<-function( Structure = Structure ) {
    
    # cerca nella cache di memoria, se attivata, la presenza della ROI già estratta
    if(ROIVoxelMemoryCache==TRUE  & !is.null(ROIVoxelMemoryCacheArray[[Structure]]) ) {
      return(ROIVoxelMemoryCacheArray[[Structure]]);
    }
    
    # try to find out which Series is the CT/MR serie
    SeriesInstanceUID<-giveBackImageSeriesInstanceUID()
    if(SeriesInstanceUID == '' ) stop("ERROR: missing CT/MR series")
    res<-getROIVoxelsFromCTRMN( Structure = Structure, SeriesInstanceUID = SeriesInstanceUID)
    class(res)<-"geoLetStructureVoxelList"
    
    
    
    
    # se la cache in memoria è attiva, salvane una copia
    if(ROIVoxelMemoryCache==TRUE) ROIVoxelMemoryCacheArray[[Structure]]<<-res;
    return( res )
  }  
  getROIVoxelsFromCTRMN<-function( Structure = Structure, SeriesInstanceUID = SeriesInstanceUID) {
    objService<-services()  
    if ( (Structure %in% getROIList()[2,]) == FALSE )  return(NA)
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
    #    browser();
    contatoreROI<-1; indiceDOM<-1;
    # for each instance number
    for (n in index) {
      # check if there is a ROI for such slice
      for (m in which(dataStorage$info[[SeriesInstanceUID]][[n]]$ROIList[,1]==Structure)) {
        
        # find the slice and gets the key for accessing at coordinates vectors
        #key<-dataStorage$info[[SeriesInstanceUID]][[n]]$ROIList[m,2]
        key<-dataStorage$info[[SeriesInstanceUID]][[n]]$ROIList[m,4]
        
        # calculate how many ROIs are co-planar
        numeroROIComplanari<-length(dataStorage$structures[[Structure]][[key]])
        # for each one of them concat the array
        for(indiceROI in seq(1,numeroROIComplanari)) {    
#          browser();
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
    return(list(
      "DOM"=array(DOM, dim = c(3,3,length(index))), 
      "final.array"=final.array, 
      "masked.images"=final.array*image.arr,
      "geometricalInformationOfImages"=getGeometricalInformationOfImage()
    )
    )
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
  cacheLoad<-function() {
    if( !(dataStorage=='') ) return;
    fileName<-attributeList$virtualMemory$fileName
    filePath<-attributeList$virtualMemory$path;
    completeFileName<-paste( c(filePath,"/",fileName), collapse=''  );
    dataStorage<<-readRDS( completeFileName )
    attributeList$virtualMemory$isOn<<-FALSE
    return();
  }
  cacheSave<-function() {
    if( is.null(names(dataStorage)) ) return();
    fileName<-attributeList$virtualMemory$fileName
    filePath<-attributeList$virtualMemory$path;
    completeFileName<-paste( c(filePath,"/",fileName), collapse=''  );
    dataStorage<<-saveRDS( object = dataStorage, file = completeFileName )
    dataStorage<<-'';
    attributeList$virtualMemory$isOn<<-TRUE
    return();
  }
  cacheDrop<-function() {
    dataStorage<<-'';
    return();
  }  
  # ...............................................................
  #  functions to support rotation
  # ...............................................................
  
  getAssociationTable<-function( tipoTabella="SOPInstance_vs_SliceLocation", ROIName ) {
    
    SeriesInstanceUID<-giveBackImageSeriesInstanceUID();
    
    if(tipoTabella=="SOPInstance_vs_SliceLocation") {
      matrice<-c()
      involvedCT<-names(dataStorage$structures[[ROIName]]);
      
      for(index in names(dataStorage$info[[SeriesInstanceUID]]) ) {
        matrice<-rbind(matrice,cbind(dataStorage$info[[SeriesInstanceUID]][[index]]$ROIList,index) );
      }
      matrice<-matrice[which(matrice[,1]==ROIName),]
      return(matrice)
    }
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
  rotate3dMatrix<-function( point.coords , angle.x = 0, angle.y = 0, angle.z = 0 ) {
    
    rotation.x<-matrix(c( 1, 0, 0, 0, cos(angle.x), sin(angle.x), 0,-sin(angle.x), cos(angle.x)),ncol=3)
    rotation.y<-matrix(c( cos(angle.y), 0, -sin(angle.y), 0, 1, 0, sin(angle.y), 0, cos(angle.y)),ncol=3)
    rotation.z<-matrix(c( cos(angle.z), sin(angle.z), 0, -sin(angle.z), cos(angle.z), 0, 0, 0, 1),ncol=3)
    
    R = (rotation.x %*% rotation.y %*% rotation.z ) %*% point.coords
    return(R);
    
  }
  rotateToAlign<-function(ROIName) {
    SeriesInstanceUID<-giveBackImageSeriesInstanceUID();
    tabella1<-getAssociationTable( "SOPInstance_vs_SliceLocation" , ROIName )
    pointList<-getROIPointList( ROIName )
    newPointList<-pointList
    iterazione<-1
    
    for(index in names( pointList )) {
      for(internalIndex in seq(1,length(pointList[[index]]))) {
        #IMGsliceInfo<-dataStorage$info[[SeriesInstanceUID]][[ tabella1[ which(tabella1[,2]==index)  ,3  ]  ]]
        #sliceLocation<-as.numeric(tabella1[which(tabella1[,2]==index),3])
        IMGsliceInfo<-dataStorage$info[[SeriesInstanceUID]][[ tabella1[ which(tabella1[,4]==index)  ,5  ]  ]]
        sliceLocation<-as.numeric(tabella1[which(tabella1[,4]==index),5])
        DOM<-IMGsliceInfo$orientationMatrix
        
        m<-pointList[[index]][[internalIndex]];
        
        numeroRighe<-dim(m)[1]
        for(i in seq(1,numeroRighe)) {
          valori<-calcolaNX(m[i,],DOM )
          newPointList[[index]][[internalIndex]][i,1]<-valori$Nx
          newPointList[[index]][[internalIndex]][i,2]<-valori$Ny
          newPointList[[index]][[internalIndex]][i,3]<-sliceLocation
          #newPointList[[index]][[internalIndex]][i,3]<-valori$Nz
        }
      }
      iterazione<-iterazione+1
    }
    return( list("pointList"=newPointList) ) 
  }  
  importStructures<-function(strutturaDataStorage) {
    dataStorage$structures<<-strutturaDataStorage
    # Associate ROI and Images
    if(is.list(dataStorage[["structures"]])) {
      associateROIandImageSlices(relaxCoPlanarity = TRUE);
    }    
  }
  
  #=================================================================================
  # Constructor
  #=================================================================================
  constructor<-function( ROIVoxelMemoryCache = TRUE , folderCleanUp = FALSE) {
    dataStorage <<- list();   
    dataChache <<- list();
    attributeList<<-list()
    attributeList$verbose<<-list("lv1"=TRUE,"lv2"=TRUE,"lv3"=FALSE,"onScreen"=TRUE,"onFile"=FALSE)
    logObj<<-logHandler()
    logObj$setOutput( list("onScreen" = attributeList$verbose$onScreen,   "onFile" = attributeList$verbose$onFile )  )
    logObj$do("clearOutputFile")
    objServ<<-services()
    attributeList$virtualMemory<<-list();
    attributeList$virtualMemory$isOn<<-FALSE;
    attributeList$virtualMemory$path<<-'./cache';
    attributeList$virtualMemory$kindOfCache<<-'dataStorage';
    attributeList$virtualMemory$fileName<<-paste(c(format(Sys.time(), "%a%b%d_%H%M%S%Y_"),as.integer(runif(1)*10000),"_",as.integer(runif(1)*10000) ) , collapse='');
    ROIVoxelMemoryCache<<-ROIVoxelMemoryCache;
    folderCleanUp<<-folderCleanUp;
    mainFrameOfReferenceUID<<-NA
    
  }
  constructor( ROIVoxelMemoryCache , folderCleanUp)
  return(list(openDICOMFolder=openDICOMFolder,getAttribute=getAttribute,
              getDICOMTag=getDICOMTag,getROIList=getROIList,getROIPointList=getROIPointList,
              setAttribute=setAttribute,getFolderContent=getFolderContent,getROIVoxels=getROIVoxels,
              getGeometricalInformationOfImage=getGeometricalInformationOfImage,
              getImageVoxelCube=getImageVoxelCube,
              cacheLoad=cacheLoad, cacheSave=cacheSave, cacheDrop = cacheDrop, 
              getAlignedStructureAndVoxelCube = getAlignedStructureAndVoxelCube,
              getPixelSpacing = getPixelSpacing,
              importStructures=importStructures,
              rotateToAlign = rotateToAlign,
              getDoseVoxelCube = getDoseVoxelCube,
              getDoseInformation = getDoseInformation,
              giveBackImageSeriesInstanceUID = giveBackImageSeriesInstanceUID,
              getXYZFromNxNyNzOfImageVolume = getXYZFromNxNyNzOfImageVolume,
              getXYZFromNxNyNzOfDoseVolume = getXYZFromNxNyNzOfDoseVolume
  ))
}
