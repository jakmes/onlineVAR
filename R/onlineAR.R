onlineAR <- function(y, nu0 = 1, lags = 1, ahead = 1, trace = FALSE) {
  d <- duplicatedata(y, lags, ahead = ahead)
  names <- rownames(d)
  d <- cbind(d[,1], 1, d[,-1])
  d <- d[(lags+ahead):(nrow(d)-lags+1),]
  R <- 10*diag(ncol(d)-1)
  coef <- rep(0, ncol(d)-1)
  pred <- coef2 <- list()
  for(n in 1:(nrow(d)-ahead)) {
    ytk <- d[n,-1]
    R <- nu0*R + ytk%*%t(ytk)
    pred2 <- if(ahead == 1 & n>1) pred[[n]] else t(coef)%*%ytk
    epsilon <- d[n,1]- pred2
    coef <- coef + solve(R)%*%ytk*as.numeric(epsilon)
    if(trace) coef2[[n]] <- coef
    pred[[n+ahead]] <- t(coef)%*%d[n+ahead, -1]
  }

  pred <- unlist(pred)
  pred <- c(rep(NA, lags+ahead*2-1), pred)
  names(pred) <- names[1:length(pred)]

  coef <- as.numeric(coef)
  names(coef) <- c("intercept", paste0("lag-", rep(1:lags)))

  if(trace) {
    trace <- as.data.frame(t(do.call(cbind, coef2)))
    trace <- rbind(matrix(NA, ncol = ncol(trace), nrow = lags+ahead-1), trace)
    rownames(trace) <- names[1:nrow(trace)]
    colnames(trace) <- names(coef)
  }
  

  list(coef = coef, 
    pred = pred, 
    nu0 = nu0,
    lags = lags,
    ahead = ahead,
    trace = trace) 
}



#areg <- function(y, lags, ahead, nu0 = 1) {
#  y <- as.numeric(y)
#  N <- length(y)
#  d <- cbind(y[(ahead + lags):N], 1)
#  for(j in 1:lags) d <- cbind(d, y[(lags-j+1):(N-ahead-j+1)])
#  R <- 10*diag(ncol(d)-1)
##  for(n in 1:100) {  
##    R <- nu0*R + d[n,-1]%*%t(d[n,-1])
##  }
#  coef <- rep(0, ncol(d)-1)
#  pred <- list()
#  for(n in 1:nrow(d)) {
#    ytk <- d[n,-1]
#    R <- nu0*R + ytk%*%t(ytk)
#    pred[[n]] <- t(coef)%*%ytk
#    epsilon <- d[n,1]- pred[[n]]
#    coef <- coef + solve(R)%*%ytk*as.numeric(epsilon)
#    
#  }
##  pred[[n+1]] <- t(coef)%*%ytk
#  pred <- unlist(pred)[-1]
#  list(coef = coef, pred = pred)
#  
#}
