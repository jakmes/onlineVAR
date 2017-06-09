onlineVAR.control <- function(nlambda = NULL, lambda.min.ratio = NULL, abstol = 0.001, 
  seq = NULL, trace = FALSE, batch = 0, batch.cv = TRUE, start = NULL, 
  parallel = FALSE, ignore = 0, ...) 
{
  if(is.null(start)) {
    if(!is.null(seq)) {
      if(!is.null(nlambda) | !is.null(lambda.min.ratio)) warning("if seq is supplied nlambda and lambda.min.ratio are ignored")
      nlambda <- length(seq)
      lambda.min.ratio <- tail(seq, 1)
    } else {
      if(is.null(nlambda)) nlambda <- 10
      if(is.null(lambda.min.ratio)) lambda.min.ratio <- 0.0001
      ## lambda sequence relative to lambda max (minimum lambda such that all coef are 0)
      seq <- exp(-log(lambda.min.ratio)/(nlambda-1)*(nlambda-c(1:nlambda))) * 
        lambda.min.ratio
    }     
  } else {
    if(!is.null(nlambda) | !is.null(lambda.min.ratio) | !is.null(seq)) warning("nlambda, lambda.min.ratio, seq taken from start. Input is ignored")
    seq <- start$seq    
    nlambda <- length(seq)
    lambda.min.ratio <- tail(seq, 1)
    ignore <- 0
    batch <- 0
    batch.cv <- 0
    
  }
  rval <- list(nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, 
    abstol = abstol, seq = seq, trace = trace, batch = batch, 
    batch.cv = batch.cv, start = start, parallel = parallel, 
    ignore = ignore, ...)
  rval
}





onlineVAR <- function(y, nu0 = 1, lags = NULL, ahead = NULL,
   control = onlineVAR.control(...), ...) {

  start <- control$start
  if(is.null(start)) {
    if(is.null(lags)) lags <- 1
    if(is.null(ahead)) ahead <- 1
  } else {
    if(!is.null(lags) | !is.null(ahead)) warning("lags and ahead taken from start. Input is ignored")
    lags <- start$lags
    ahead <- start$ahead
  }

  batch <- control$batch
  
  if(batch >= nrow(y)) stop("batch must be lower than the length of y. 
    Use batchVAR for pure batch estimation")

  if(batch < 1) batch <- batch * nrow(y)

  Q <- ncol(y)
  P <- Q*lags
  
 

  ## create starting values 
  if(is.null(start)) {
    if(batch) {
      weights <- nu0^(batch-1:batch)
      start <- batchVAR(y[1:batch, ], lags = lags, ahead = ahead, control = control, 
        weights = weights, standardize = FALSE, cv = control$batch.cv)
      nlambda <- length(start$seq)
      seq <- start$seq
    } else {
      coef2 <- matrix(0, ncol = P, nrow = P/lags)
      if(lags > 1) coef2 <- extendmatrix(coef2, lags)
      zeromatrix <- matrix(0, ncol = P, nrow = P)
      coef <- werror <- R <- r <- list()
      for(i in 1:control$nlambda) {
        coef[[i]] <- coef2
        R[[i]] <- zeromatrix
        r[[i]] <- zeromatrix
      }
      fit <- list(coef = coef, R = R, r = r, mu = rep(0, P), wN = 0, 
        werror = rep(0, control$nlambda))
      start <- list(fit = fit)
    }
  } 

  ## perform onlineVAR on data not used in batch estimation

  ## create matrix y such that model can be fitted as VAR1
  names <- colnames(y)
  y <- duplicatedata(y, lags)
  if(batch) y <- y[(batch + 1):nrow(y), ]
  names2 <- c(rownames(y), rep("NA", ahead))
  y <- na.omit(y)

  x <- head(y, -ahead)
  y <- tail(y, -ahead)
  
  N <- nrow(y)



  ## call the actual workhorse lassovar.fit()
  rval <- onlineVAR.fit(x = x, y = y, nu0 = nu0, lags = lags, ahead = ahead, 
    fit = start$fit, control = control)
    
  ## fill up prediction matrix
  ## NA predictions for the first lags+ahead rows
  rval$pred <- rbind(matrix(NA, nrow = (batch==0)*(lags-1) + ahead, ncol = Q), rval$pred)
  rval$residuals <- rbind(matrix(NA, nrow = (batch==0)*(lags-1) + ahead, ncol = Q), rval$residuals)
  ## row and column names
  rownames(rval$pred) <- names2[1:nrow(rval$pred)]
#  rval$pred <- rbind(rval$pred, "NA" = pred2)
  colnames(rval$pred) <- colnames(y)[1:Q]


  rval$batch <- batch
  return(rval)
}



###############################################################################
## Batch VAR -- wrapper function for glmnet with functions from sparsevar #####
###############################################################################

batchVAR <- function(y, lags = 1, ahead = 1, nlambda = 100, 
  lambda.min.ratio = 0.0001, seq = NULL, weights = NULL, standardize = TRUE, 
  cv = TRUE, ...) {

  y <- as.matrix(y)
  N <- nrow(y)
  if(is.null(weights)) weights <- rep(1, N)
  
  Q <- ncol(y)
  P <- ncol(y)*lags

  ## center
  mu <- colSums(weights*y)
  wN <- sum(weights)
  y <- y - matrix(rep(mu/wN, N), ncol = Q, byrow = TRUE)

  ## normalize
  if(standardize) {
    sigma <- colSums(weights*y^2)
    std <- sqrt(sigma/wN)
    y <- y/matrix(rep(std, N), ncol = Q, byrow = TRUE)
  } else sigma <- NA


  ## derive r and R to update fit
  y3 <- na.omit(duplicatedata(y, lags))
  x2 <- head(y3, -ahead)
  y3 <- tail(y3, -ahead)
  


  
  r <- R <- 0
  nu0 <- weights[length(weights)-1]
  for(i in 1:nrow(y3)) {
    R <- nu0*R + x2[i,]%*%t(x2[i,])
    r <- nu0*r + y3[i,]%*%t(x2[i,])
  }
  

 if(is.null(seq)) 
    seq <- exp(-log(lambda.min.ratio)/(nlambda-1)*(nlambda-c(1:nlambda))) * 
      lambda.min.ratio

  lambdaseq <- max(abs(r))*seq/N/Q
  


  ## fit batch VAR with glmnet and transformData function from sparseVAR

  if(ahead > 1) {
    X <- transformData(head(y, -(ahead - 1)), p = lags, 
      list(method = "other", center = FALSE))
    x <- X$X
    y2 <- as.vector(tail(y, -(ahead+lags-1)))
    weights <- tail(weights, -(ahead - 1))
  } else {
    X <- transformData(y, p = lags, list(method = "other", center = FALSE))
    x <- X$X
    y2 <- X$y
  }
  lasso <- if(cv) {
      foldid <- sort(sample(10, N, replace = TRUE))
      foldid <- as.vector(matrix(rep(foldid, ncol(y)), ncol = ncol(y))[(ahead+lags):N,])

      cv.glmnet(x, y2, intercept = FALSE, standardize = FALSE, 
        weights = rep(tail(weights, -lags), Q), foldid = foldid, 
        lambda = lambdaseq, type.measure = "mse", ...)
    } else {
      glmnet(x, y2, intercept = FALSE, standardize = FALSE, 
        weights = rep(tail(weights, -lags), Q), lambda = lambdaseq, ...)
    }
    

  ## prepare rval
  coef <- list()
  ## extract coefficient matrices
  nlambda <- length(lasso$lambda)
  for(frac in 1:nlambda) {
    coef[[frac]] <- matrix(coef(lasso, s = lasso$lambda[frac])[-1], nrow = Q, ncol = P, 
      byrow = TRUE)
    if(lags > 1) coef[[frac]] <- extendmatrix(coef[[frac]], lags)
  }

  
  mu <- rep(mu, lags)
  if(standardize) {
    sigma <- rep(sigma, lags)
    std <- rep(std, lags)
  }
  
  if(cv) {
    minerror <- which(lasso$lambda == lasso$lambda.min)   
    coef2 <- if(standardize) coef[[minerror]]*std%*%t(1/std) else coef[[minerror]]
    intercept <- (mu - coef2%*%mu)/wN
    coef3 <- list(intercept = intercept[1:Q])
    for(l in 1 : lags) coef3[[paste0("lag", l)]] <- coef2[1:Q, ((l-1)*Q+1):(l*Q)]
  } else coef3 <- NULL

  werror <- if(cv) Q * wN * lasso$cvm/9*10 else rep(0, nlambda)

  fit <- list(coef = coef, R = R, r = r, mu = mu, sigma = sigma, wN = wN, 
    werror = werror)

  rval <- list(
    coef = coef3, 
    fit = fit,
    nu0 = nu0,
    lags = lags,
    ahead = ahead,
    nlambda = nlambda,
    lambda.min.ratio = lambda.min.ratio,
    weights = weights,
    standardize = standardize,
    cv = cv,
    seq = lasso$lambda/lasso$lambda[1])
  class(rval) <- c("batchVAR", "onlineVAR")
  return(rval)
}


###############################################################################
## Online VAR fitting function ################################################
###############################################################################




onlineVAR.fit <- function(x, y, nu0, lags, ahead, fit, control = control) {
  seq <- control$seq
  abstol <- control$abstol
  trace <- control$trace
  parallel <- control$parallel
  ignore <- control$ignore
  nlambda <- control$nlambda

  N <- nrow(y)
  P <- ncol(y)
  Q <- P/lags

  pred <- coefall <- intercept <- pred2 <- list()
  
  fitfun <- function(s) {
    fit2 <- .Call("onlineVARfit", x = x, y = y, nu0 = nu0, mui = fit$mu,  
      wNi=fit$wN, Ri = fit$R[[s]], ri = fit$r[[s]], seq = seq[s], coefi = fit$coef[[s]], 
      q = as.integer(Q), abstol = abstol, trace = as.integer(trace))
    fit2
  }

  fit2 <- if(parallel) mclapply(1:length(seq), fitfun) else 
    lapply(1:nlambda, fitfun)

  error <- lapply(1:length(seq), 
    function(s) c(fit$werror[s], rowSums((fit2[[s]][[6]] - y[,1:Q])^2)))

  werrorfun <- function(error) {
    werror0 <- error[1]
    error <- error[-1]
    rval <- rep(0, N)
    rval[ignore+1] <- nu0*werror0 + error[ignore+1] 
    for(n in (ignore+2):N) rval[n] <- nu0*rval[n-1] + error[n]
    rval
  }
  werror <- sapply(error, werrorfun)
#browser()
  pred <- matrix(NA, nrow = N, ncol = Q)
  for(n in (ignore+1):N) {
    pred[n,] <- fit2[[which.min(werror[n,])]][[6]][n,]
  }


  if(trace) {
    mu <- fit$mu
    wN <- fit$wN
    coefall <- list()
    for(n in 1:(N-ahead)+1) {
      mu <- nu0*mu + y[n,]
      wN <- nu0*wN + 1    
      coefn <- fit2[[which.min(werror[n-1,])]][[7]][((n-1)*P+1):(n*P),]
      interceptn <- (mu - coefn%*%mu)/wN
      coefall[[n+ahead]] <- list(intercept = interceptn[1:Q])
      for(l in 1:lags) coefall[[n+ahead]][[paste0("lag", l)]] <- coefn[1:Q, ((l-1)*Q+1):(l*Q)]
    }
  }


  coef <- R <- r <- list()
  for(i in 1:nlambda) {
    coef[[i]] <- fit2[[i]][[5]]
    R[[i]] <- fit2[[i]][[3]]
    r[[i]] <- fit2[[i]][[4]]
  }
  mu <- fit2[[1]][[1]]
  wN <- fit2[[1]][[2]]

  fit <- list(coef = coef, R = R, r = r, mu = mu, wN = wN, 
    werror = werror[N-ahead,])


  s <- which.min(werror[N,])

  
  intercept <- (mu - coef[[s]]%*%mu)/wN
  coef2 <- list(intercept = intercept[1:Q])
  for(l in 1 : lags) coef2[[paste0("lag", l)]] <- coef[[s]][1:Q, ((l-1)*Q+1):(l*Q)]



  

  if(trace) {
    trace <- coefall
  } 
  rval <- list(coef = coef2, 
    fit = fit, 
    nu0 = nu0,
    pred = pred,
    residuals = y[,1:Q] - pred,
    lags = lags,
    ahead = ahead,
    weights = weights,
    cv = TRUE,
    seq = seq,
    trace = trace)

  class(rval) <- "onlineVAR"
  return(rval)

}





###############################################################################
## Extractor functions ########################################################
###############################################################################
predict.onlineVAR <- function(object, newdata, s = NULL) {
  if(is.null(s)) {
    if(object$cv) {
      s <- which.min(object$fit$werror)
    } else {
      error("s has to be specified for cv = FALSE")
    } 
  } 
  wN <- object$fit$wN
  std <- object$fit$sigma/wN
  mu <- object$fit$mu
  coef <- object$fit$coef
  coef2 <-  coef[[s]]
  intercept <- (mu - coef2%*%mu)/wN
  
  Q <- (ncol(coef2)/object$lags)
  y <- duplicatedata(newdata, object$lags, object$ahead)
  x <- y[, (Q + 1):ncol(y)]
  if(object$lags > 1) x <- head(x, - (object$lags - 1))
  rval <- matrix(rep(intercept, nrow(x)), nrow = nrow(x), byrow = TRUE) + t(coef2 %*% t(x))
  rval <- rval[,1:Q]
  colnames(rval) <- colnames(y)[1:Q]
  return(rval)
}

coef.onlineVAR <- function(object, s = NULL) {
  if(is.null(s)) {
    if(object$cv) {
      s <- which.min(object$fit$werror)
      coef3 <- object$coef
    } else {
      error("s has to be specified for cv = FALSE")
    } 
  } else {
    Q <- ncol(object$pred)
    lags <- object$lags
    wN <- object$fit$wN
    std <- object$fit$sigma/wN
    mu <- object$fit$mu
    coef <- object$fit$coef
    coef2 <-  coef[[s]]
    intercept <- (mu - coef2%*%mu)/wN

    coef3 <- list(intercept = intercept[1:Q])
    for(l in 1 : lags) coef3[[paste0("lag", l)]] <- coef2[1:Q, ((l-1)*Q+1):(l*Q)]
  }

  rval <- coef3
  return(rval)
}




#plot.onlineVAR <- function(x, lags = 1, col = NULL, ...)  {
#  if(is.null(col)) col <- rev(sequential_hcl(16, l = c(30, 100),power = 0.5))
#  layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(4,1) )
#  image(1:10, 1:10, a[, ncol(a):1], col=col)

#  lim <- range(a, na.rm=TRUE)
#  breaks <- seq(lim[1], lim[2], length.out=(length(col)+1))

#  par(mar = c(5,2,4,2))
#  plot(1, 1, type = "n", xlim = c(0,1), ylim = range(breaks), 
#    xaxt = "n", yaxt = "n", xaxs = "i", yaxs="i", xlab = "", ylab = "")
#  axis(4)
#  for(i in 1:length(col)) {
#    polygon(c(0,0,1,1), breaks[c(i, i+1, i+1, i)], col = col[i], border = NA)
#  }
#}


#  for(i in lags){
#    a <- abs(x$coef[[i+1]])
#    plot1 <- levelplot(a[, ncol(a):1], col.regions=pal, xlab = "site", 
#  ylab = "site",main = expression(paste("Coefficient matrix ", A[1])))

#}



###############################################################################
## additional functions required for other functions ##########################
###############################################################################


## function to get response matrix
duplicatedata <- function(y, lags, ahead = 0) {
  
  y2 <- y3 <- as.matrix(y)
  names <- rownames(y2)
  if(is.null(names)) names <- 1:nrow(y2)
  if(ahead) {
    for(i in 1:ahead) {
      y2 <- rbind(y2, NA)
      y3 <- rbind(NA, y3)
      names <- c(names, "NA")
    }
    y2 <- cbind(y2, y3)
  }

  if(lags>1) {
    for(i in 1:(lags-1)) {
      y3 <- rbind(NA, y3)
      y2 <- cbind(rbind(y2, NA), y3)
      names <- c(names, "NA")
#      y2 <- cbind(tail(y2, -1), as.matrix(head(y, -i)))
    }
  } 
  rownames(y2) <- names
  
  return(y2)
}

## function to extend coefficientmatrix to transform VAR-p into VAR-1 model
extendmatrix <- function(coef, lags) {
  P <- ncol(coef)/lags
  coef <- rbind(coef, matrix(0, ncol = P*lags, nrow = P*(lags-1)))
  for(j in 1:(lags-1)) {
    coef[(j*P) + 1:P, ((j-1)*P) + 1:P] <- diag(P)
  }
  coef
}



 
