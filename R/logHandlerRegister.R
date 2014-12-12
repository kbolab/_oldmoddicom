#' class for handling logs/warnings/errors
#' 
#' @description  It handles messages from script to a chosen output (screen, file, etc.)
#' @export
logHandlerRegister<-function() {
  
  registry<-list()
  logObj<-list()
  
  # ------------------------------------------------
  # subscription
  # request a subscription
  # ------------------------------------------------  
  subscription<-function(ticket="") {
    # if ticket is not specified it means that it is requested
    if( ticket=="" ) {
      # build a new ticket
      ticket<-as.character(as.integer(runif(1)*100000000))
      # if it is already extant, well, bye!
      if( length(registry[[ticket]]) != 0) {
        logObj$sendLog("Ticket Creation failed. Duplicated.")
        return();
      }
      # subscribe it and return
      registry[[ticket]]<<-list();
      registry[[ticket]][["objHandler"]]<<-logHandler()
      return(ticket)
    }
    # if a registry is specified and already exist
    if( length(registry[[ticket]]) > 0 )  {
      logObj$sendLog("Subscription failed: extant ticket")
      return();    
    } else {
      # subscribe it and return
      registry[[ticket]]<<-list();
      return(ticket)  
    }    
  }
  sendLog<-function(ticket,message) {
    if( length(registry[[ticket]]) == 0  ) {
      logObj$sendLog( paste(c("Missing subscription for ticket '",ticket,"'"), collapse = "") )
      return();    
    }    
    registry[[ticket]][["objHandler"]]$sendLog( message )
    return()    
  }
  do<-function(ticket,what2Do, arguments="") {
    if( length(registry[[ticket]]) == 0  ) {
      logObj$sendLog( paste(c("Missing subscription for ticket '",ticket,"'"), collapse = "") )
      return();    
    }    
    registry[[ticket]][["objHandler"]]$do( what2Do , arguments )
    return()    
  }
  setOutput<-function(ticket,attributeList) {
    if( length(registry[[ticket]]) == 0  ) {
      logObj$sendLog( paste(c("Missing subscription for ticket '",ticket,"'"), collapse = "") )
      return();    
    }    
    registry[[ticket]][["objHandler"]]$setOutput( attributeList )
    return()     
  }
  
  # ------------------------------------------------
  # constructor
  # ------------------------------------------------  
  constructor<-function() {
    registry<<-list()
    logObj<<-logHandler();    
  }
  constructor();  
  return(list(subscription=subscription , sendLog=sendLog, do=do , setOutput=setOutput))  
}