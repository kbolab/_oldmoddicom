
  
    # create matrices of normalized intensity values
    y.pos<-matrix(data = unlist(lapply(X = pos$BAVA$details$interpolatedD, FUN = function(x) return(x$y))), nrow = length(pos$BAVA$details$interpolatedD), byrow = TRUE)
    y.neg<-matrix(data = unlist(lapply(X = neg$BAVA$details$interpolatedD, FUN = function(x) return(x$y))), nrow = length(neg$BAVA$details$interpolatedD), byrow = TRUE)
    
      MW.test<-c()
      for ( n in 2:min(c(ncol(y.neg), ncol(y.pos)))) {
        MW.test<-c(MW.test, wilcox.test(x = y.pos[,n], y = y.neg[,n], conf.int = FALSE, alternative = 'less')$p.value)  
      }
      # matrix of p-values (y) and signal values (x)
      MW.test<-cbind(c(2:n), MW.test)

    
      KS.test.greater<-ks.test(x = colMeans(x = y.pos), y = colMeans(x = y.neg), alternative = 'greater')
    print(KS.test.greater)
      KS.test.less<-ks.test(x = colMeans(x = y.pos), y = colMeans(x = y.neg), alternative = 'less')
    print(KS.test.less)
      KS.test.dual<-ks.test(x = colMeans(x = y.pos), y = colMeans(x = y.neg))    
    print(KS.test.dual)


      plot(x = MW.test[500:1000, 1], y = MW.test[500:1000,2], type = 'l', log ='y', lwd = 2, xlab = 'Urine Signal Normalized MRI', ylab = 'Mann Withney P-Value')
      abline (h = .05, lty = 2, col ='red')


      plot(x = c(1:2221), y = colMeans(x = y.pos), col = 'red', lwd = 2, xlim = c(1,1500),
           xlab = 'Urine Signal Normalized MRI', ylab = 'Mean KDF', type = 'l', 
           ylim = c(0, 0.0025))
      lines(x = c(1:2506), y = colMeans(x = y.neg), col = 'blue', lwd = 2)
      legend(x = 1000, y = .002, legend = c('Positive', 'Negative'), col = c('red', 'blue'), lty = 1)    
