# # create vector of urine values
# message('Vectorize mmButo data...')
# Vurine <- sapply(X = Urine, FUN = function(x) return (as.vector(x$masked.images$voxelCube[x$masked.images$voxelCube!=0])))
# 
# # create normalization indeces
# norm.ind <- sapply(X = Vurine, FUN = function(x) return(mean(x) / max(sapply(X = Vurine, FUN = mean))))
# 
# # create normalized GTV values
# message('Normalize GTV data...')
# norm.GTV <- list()
# unnorm.GTV <- list()
# for (n in 1:length(norm.ind)) {
#   norm.GTV[[n]] <- GTV[[n]]$masked.images$voxelCube / norm.ind[n]
#   norm.GTV[[n]] <- as.vector(norm.GTV[[n]][norm.GTV[[n]]!=0])
#   unnorm.GTV[[n]] <- GTV[[n]]$masked.images$voxelCube
#   unnorm.GTV[[n]] <- as.vector(unnorm.GTV[[n]][unnorm.GTV[[n]]!=0])
# }
# 
# # create normalized Urine values
# message('Normalize urine data...')
# Vurine.norm <- list()
# for (n in 1:length(norm.ind)) {
#   Vurine.norm[[n]] <- Urine[[n]]$masked.images$voxelCube / norm.ind[n]
#   Vurine.norm[[n]] <- as.vector(Vurine.norm[[n]][Vurine.norm[[n]]!=0])
# }
# 
# #########################################################
# # create density of urine signals without normalization #
# #########################################################
# message('Create density functions of unnormalized urine data...')
# Vurine.density <- sapply(X = Vurine, FUN = density, from = ceiling(min(sapply(X = Vurine, FUN = min))), to = floor(max(sapply(X = Vurine, FUN = max))), simplify = FALSE)
# # draws the plot of urine densities without normalization
# plot(x = Vurine.density[[1]]$x, y = Vurine.density[[2]]$y, ylim = c(0, max(sapply(X = Vurine.density, FUN = function(x) return(max(x$y))))), 
#      type = 'l', col = 'grey', xlab = 'Unnormalized MRI signal', ylab = 'KDE')
# for (n in 2:length(norm.ind)) 
#   lines(x = Vurine.density[[1]]$x, y = Vurine.density[[n]]$y, col = 'grey')
# # create the matrix of y levels
# Vurine.unnorm.Y <- c()
# for (n in 1:length(norm.ind)) 
#   Vurine.unnorm.Y <- rbind(Vurine.unnorm.Y, Vurine.density[[n]]$y)
# # plot the mean value of the signals
# lines(x = Vurine.density[[1]]$x, y = apply(X = Vurine.unnorm.Y, MARGIN = 2, FUN = mean), lwd = 2, col = 'black')
# # create qqplot
# z.norm <- (unname(unlist(Vurine)) - mean(unlist(Vurine))) / sd(unlist(Vurine))
# qqnorm(z.norm, pch = '.')
# abline(0,1, lwd = 2, col = 'red')
# 
# #########################################################
# # create density of urine signals with normalization    #
# #########################################################
# message('Create density functions of normalized urine data...')
# Vurine.norm.density <- sapply(X = Vurine.norm, FUN = density, from = ceiling(min(sapply(X = Vurine.norm, FUN = min))), to = floor(max(sapply(X = Vurine.norm, FUN = max))), simplify = FALSE)
# # draws the plot of urine densities without normalization
# plot(x = Vurine.norm.density[[1]]$x, y = Vurine.norm.density[[2]]$y, ylim = c(0, max(sapply(X = Vurine.norm.density, FUN = function(x) return(max(x$y))))), 
#      type = 'l', col = 'grey', xlab = 'Normalized MRI signal', ylab = 'KDE')
# for (n in 2:length(norm.ind)) 
#   lines(x = Vurine.norm.density[[1]]$x, y = Vurine.norm.density[[n]]$y, col = 'grey')
# # create the matrix of y levels
# Vurine.norm.Y <- c()
# for (n in 1:length(norm.ind)) 
#     Vurine.norm.Y <- rbind(Vurine.norm.Y, Vurine.norm.density[[n]]$y)
# # plot the mean value of the signals
# lines(x = Vurine.norm.density[[1]]$x, y = apply(X = Vurine.norm.Y, MARGIN = 2, FUN = mean), lwd = 2, col = 'black')
# # creating bootstrap analysis for Confidence Intervals
# message('Bootstrapping confidence intervals...')
# indeces <- sample(x = c(1:182), size = 182*5000, replace = TRUE)
# big.matrix <- rbind(Vurine.norm.Y[indeces, ])
# Vurine.norm.CI <- apply(X = big.matrix, MARGIN = 2, FUN = quantile, probs = c(.025, .975))
# lines(x = Vurine.norm.density[[1]]$x, Vurine.norm.CI[1,])
# lines(x = Vurine.norm.density[[1]]$x, Vurine.norm.CI[2,])
# rm(big.matrix)
# # creating QQ plot
# z.norm.norm <- (unname(unlist(Vurine.norm)) - mean(unlist(Vurine.norm))) / sd(unlist(Vurine.norm))
# qqnorm(z.norm.norm, pch = '.')
# abline(0,1, lwd = 2, col = 'red')
# 
# 
# #########################################################
# # create density of GTV signals with normalization      #
# #########################################################
# message('Create density functions of normalized GTV data...')
# GTV.norm.density <- sapply(X = norm.GTV, FUN = density, from = ceiling(min(sapply(X = norm.GTV, FUN = min))), to = floor(max(sapply(X = norm.GTV, FUN = max))), simplify = FALSE)
# # draws the plot of GTV densities without normalization
# plot(x = GTV.norm.density[[1]]$x, y = GTV.norm.density[[2]]$y, ylim = c(0, max(sapply(X = GTV.norm.density, FUN = function(x) return(max(x$y))))), 
#      type = 'l', col = 'grey', xlab = 'Normalized MRI signal', ylab = 'KDE')
# for (n in 2:length(norm.ind)) 
#   lines(x = GTV.norm.density[[1]]$x, y = GTV.norm.density[[n]]$y, col = 'grey')
# # create the matrix of y levels
# GTV.norm.Y <- c()
# for (n in 1:length(norm.ind)) 
#   GTV.norm.Y <- rbind(GTV.norm.Y, GTV.norm.density[[n]]$y)
# # plot the mean value of the signals
# lines(x = GTV.norm.density[[1]]$x, y = apply(X = GTV.norm.Y, MARGIN = 2, FUN = mean), lwd = 2, col = 'black')
# # creating bootstrap analysis for Confidence Intervals
# message('Bootstrapping confidence intervals...')
# indeces <- sample(x = c(1:182), size = 182*5000, replace = TRUE)
# big.matrix <- rbind(GTV.norm.Y[indeces, ])
# GTV.norm.CI <- apply(X = big.matrix, MARGIN = 2, FUN = quantile, probs = c(.025, .975))
# lines(x = GTV.norm.density[[1]]$x, GTV.norm.CI[1,])
# lines(x = GTV.norm.density[[1]]$x, GTV.norm.CI[2,])
# rm(big.matrix)
# # creating QQ plot
# z.GTV.norm <- (unname(unlist(norm.GTV)) - mean(unlist(norm.GTV))) / sd(unlist(norm.GTV))
# qqnorm(z.GTV.norm, pch = '.')
# abline(0,1, lwd = 2, col = 'red')
# 
# #########################################################
# #    Test normality for 2 sd of normalized urine        #
# #########################################################
# Vurine.norm.v <-unlist(Vurine.norm)
# st.subset <- Vurine.norm.v[(Vurine.norm.v > (- 2 * sd(Vurine.norm.v) + mean(Vurine.norm.v))) & (Vurine.norm.v < (2 * sd(Vurine.norm.v)) + mean(Vurine.norm.v))]
# shapiro.test(x = sample(st.subset, size = 2000, replace = F))
# 
# #########################################################
# # correlation between GTV signal and Urine signal for   #
# # justifying the use of normalization by urine          #
# #########################################################
# cat('\nPearson correlation test between unnormalized urine and unnormalized GTV values')
# cor.test(x = sapply(X = unnorm.GTV, FUN = mean), y = sapply(Vurine, FUN = mean))
# cat('\nPearson correlation test between normalized urine and normalized GTV values')
# cor.test(x = sapply(X = norm.GTV, FUN = mean), y = sapply(Vurine.norm, FUN = mean))
# # crete dataframe of meanvalues
# meanvalues<-as.data.frame(cbind('unnorm.GTV' = sapply(X = unnorm.GTV, FUN = mean, simplify = TRUE), 'unnorm.Urine' = sapply(X = Vurine, FUN = mean, simplify = TRUE),
#                                 'norm.GTV' = sapply(X = norm.GTV, FUN = mean, simplify = TRUE), 'norm.Urine' = sapply(X = Vurine.norm, FUN = mean, simplify = TRUE)))
# # create linear model
# summary(lm(meanvalues$unnorm.GTV ~ meanvalues$unnorm.Urine))
# plot(lm(meanvalues$unnorm.GTV ~ meanvalues$unnorm.Urine))
