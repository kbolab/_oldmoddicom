#' class for handling a REGISTER for logs/warnings/errors
#' 
#' @description  It provide a REGISTER for handling messages from script to a chosen output (screen, file, etc.)
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
      registry[[ticket]][["objHandler"]]<<-logHandler()
      return(ticket)  
    }    
  }
  # ------------------------------------------------
  # sendLog
  # send a log to a subscripted logHandler
  # ------------------------------------------------    
  sendLog<-function(ticket,message) {
    if( length(registry[[ticket]]) == 0  ) {
      logObj$sendLog( paste(c("Missing subscription for ticket '",ticket,"'"), collapse = "") )
      return();    
    }    
    registry[[ticket]][["objHandler"]]$sendLog( message )
    return()    
  }
  # ------------------------------------------------
  # do
  # send a do to a subscripted logHandler
  # ------------------------------------------------   
  do<-function(ticket,what2Do, arguments="") {
    if( length(registry[[ticket]]) == 0  ) {
      logObj$sendLog( paste(c("Missing subscription for ticket '",ticket,"'"), collapse = "") )
      return();    
    }    
    registry[[ticket]][["objHandler"]]$do( what2Do , arguments )
    return()    
  }
  # ------------------------------------------------
  # setOutput
  # send a setOutput to a subscripted logHandler
  # ------------------------------------------------   
  setOutput<-function(ticket,attributeList) {
    if( length(registry[[ticket]]) == 0  ) {
      logObj$sendLog( paste(c("Missing subscription for ticket '",ticket,"'"), collapse = "") )
      return();    
    }    
    registry[[ticket]][["objHandler"]]$setOutput( attributeList )
    return()     
  }
  # ------------------------------------------------
  # getSubscriptions
  # get ALL the subscriptions! (for debugging issues only!)
  # ------------------------------------------------     
  getSubscriptions<-function() {
    return( registry )
  }
  
  # ------------------------------------------------
  # constructor
  # ------------------------------------------------  
  constructor<-function() {
    registry<<-list()
    logObj<<-logHandler();    
  }
  constructor();  
  return(list(
      subscription=subscription, 
      sendLog=sendLog, 
      do=do, 
      setOutput=setOutput,
      getSubscriptions=getSubscriptions))  
}

#rm(a)
#a<-logHandlerRegister()
#a$subscription("geoLet")
#a$subscription("mmButo")
#a$subscription("mmButo")
#a$subscription()
#b[["geoLet"]][["objHandler"]]$sendLog("prova")
#b[["geoLet"]][["objHandler"]]$setOutput(list("lv1"=TRUE,"lv2"=TRUE,"onFilePar"=TRUE))
#b[["geoLet"]][["objHandler"]]$sendLog("prova")
#b[["mmButo"]][["objHandler"]]$sendLog("prova")

