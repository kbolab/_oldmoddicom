#' class for handling DICOM files (load, etc)
#' 
#' @description  Instantiate an object of the class \code{geoLet}.
#' @export
logHandler<-function() {
  onScreen<-TRUE;
  onFile<-FALSE;
  returnOnEOL<-TRUE;
  
  setOutput<-function(onScreenPar = TRUE, onFilePar = FALSE) {
    onScreen<<-onScreenPar
    onFile<<-onFilePar
  }  
  sendLog<-function(message) {
    if( onFile == TRUE ) { cat("\n Not yet implemented"); stop(); }  
    if( onScreen == TRUE ) { 
      if( returnOnEOL == TRUE) cat('\n', message) 
      else cat( message )    
    }
  }  
  constructor<-function() {
    onScreen<<-TRUE;
    onFile<<-FALSE;
    returnOnEOL<<-TRUE;
  }
  constructor();  
  return(list(setOutput=setOutput,sendLog=sendLog))  
}