#' class for handling logs/warnings/errors
#' 
#' @description  It handles messages from script to a chosen output (screen, file, etc.)
#' 
#'               Many methods are available for objects built by this class:
#'               \itemize{
#'               \item \code{setOutput(attributeList)} 
#'               \cr\cr sets the attributes of the object. Many attributes can be specified in the list:
#'               \cr
#'               \itemize{
#'                 \item \code{onFilePar} : \emph{logical}, default=\code{FALSE}. Set it to \code{TRUE} if you want the msgs  written on file;
#'                 \item \code{onScreenPar} : \emph{logical}, default=\code{TRUE} Set it to \code{TRUE} if you want the msgs written on screen;
#'                 \item \code{fileName} : \emph{string}, default='./defLogHandler.txt'. Set the file name where msgs should be written;
#'                 \item \code{returnOnEOL} : \emph{logical}, default=\code{TRUE} Does it has to add a '\\n' at the end of the msg?
#'               }
#'               \cr\cr
#'               \item \code{sendLog(message)} 
#'               \cr\cr print on screen, on file (or on a void output) the given message 
#'               \cr\cr
#'               \item \code{do(what2Do, arguments="")} 
#'               \cr\cr it send a command to the object. Many commands are available:
#'               \cr
#'               \itemize{
#'                 \item \code{clearOutputFile} : it truncates the output file. No arguments are needed.
#'               }
#'               }
#' @export
logHandler<-function() {
  onScreen<-TRUE;
  onFile<-FALSE;
  returnOnEOL<-TRUE;
  fileName<-""
  
  setOutput<-function(attributeList) {
    # onScreenPar (BOOLEAN)
    if(length(attributeList$onScreenPar)>0) onScreen<<-attributeList$onScreenPar
    # onFilePar (BOOLEAN)
    if(length(attributeList$onFilePar)>0) onFile<<-attributeList$onFilePar
    # fileName (string)
    if(length(attributeList$fileName)>0) fileName<<-attributeList$fileName
    # returnOnEOL (BOOLEAN)
    if(length(attributeList$returnOnEOL)>0) returnOnEOL<<-attributeList$returnOnEOL    
  }  
  do<-function(what2Do, arguments="") {
    if( what2Do == "clearOutputFile" ) { if (file.exists(fileName)) { file.remove(fileName) }; return(); }
  }
  sendLog<-function(message) {
    if( onFile == TRUE ) { 
      write( message, file=fileName, append=TRUE)
    }  
    if( onScreen == TRUE ) { 
      if( returnOnEOL == TRUE) cat('\n', message) 
      else cat( message )    
    }
  }  
  constructor<-function() {
    onScreen<<-TRUE;
    onFile<<-FALSE;
    returnOnEOL<<-TRUE;
    fileName<<-"./defLogHandler.txt";
  }
  constructor();  
  return(list(setOutput=setOutput,sendLog=sendLog,do=do))  
}