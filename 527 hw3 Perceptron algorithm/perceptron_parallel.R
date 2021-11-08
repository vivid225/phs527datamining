percep_function <- function(id){
  library(mvtnorm)
  
  
  ## GD -----
  k = 100
  n = 100000
  converge <- c()
  start_time <- Sys.time()
  for (i in 1:k){
    # n=10000
    X = rmvnorm(n, mean=c(0.5,1.5), sigma=matrix(c(1,0.5,0.5,1), ncol = 2))
    x1 = X[,1]
    x2 = X[,2]
    X = rmvnorm(n, mean=c(5,10), sigma=matrix(c(1,0.5,0.5,1), ncol = 2))
    x11 = X[,1]
    x22 = X[,2]
    
    X <- matrix(c(x1,x11,x2,x22),ncol=2)
    X <- cbind(rep(1,nrow(X)), X)
    
    Y <- matrix(rep(c(-1,1),each=n), ncol=1)
    
    rho = 0.5
    
    j = 1
    bmat <- matrix(c(1,1,1),ncol=1) # initial bmat = matrix(b0,b1,b2,ncol=1)
    m = 10000
    while (j <= m){
      
      xb <- X%*%bmat
      # ypred <- ifelse(xb <=0, -1, 1)
      # Y * xb
      if (all((Y * xb) > 0)==FALSE){
        misind <- which(Y * xb < 0)
        ### update
        Ymis <- matrix(Y[misind,], ncol=1)
        Xmis <- matrix(X[misind,], ncol=3)
        
        bmat.new <- bmat + rho * (t(Xmis) %*% Ymis)
        bmat <- bmat.new
        j = j+1
        if (j == m+1){
          converge <- c(converge, 0)
        }
      } else {
        # print(bmat)
        # print(j)
        converge <- c(converge, 1)
        break
      }
    }
  }
  end_time <- Sys.time()
  time_gd <- end_time - start_time
  convgd <- mean(converge)
  
  
  ## Stochastic gradient descent -----
  
  converge_sgd <- c()
  start_time <- Sys.time()
  for (i in 1:k){
    # n=10000
    X = rmvnorm(n, mean=c(0.5,1.5), sigma=matrix(c(1,0.5,0.5,1), ncol = 2))
    x1 = X[,1]
    x2 = X[,2]
    X = rmvnorm(n, mean=c(5,10), sigma=matrix(c(1,0.5,0.5,1), ncol = 2))
    x11 = X[,1]
    x22 = X[,2]
    
    X <- matrix(c(x1,x11,x2,x22),ncol=2)
    X <- cbind(rep(1,nrow(X)), X)
    
    Y <- matrix(rep(c(-1,1),each=n), ncol=1)
    
    rho = 0.5
    
    j = 1
    bmat <- matrix(c(1,1,1),ncol=1) # initial bmat = matrix(b0,b1,b2,ncol=1)
    m = 10000
    while (j <= m){
      
      xb <- X%*%bmat
      # ypred <- ifelse(xb <=0, -1, 1)
      # Y * xb
      if (all((Y * xb) > 0)==FALSE){
        misind <- which(Y * xb < 0)
        misind <- sample(misind,1)
        ### update
        Ymis <- matrix(Y[misind,], ncol=1)
        Xmis <- matrix(X[misind,], ncol=3)
        
        bmat.new <- bmat + rho * (t(Xmis) %*% Ymis)
        bmat <- bmat.new
        j = j+1
        if (j == m+1){
          converge_sgd <- c(converge_sgd, 0)
        }
      } else {
        # print(bmat)
        # print(j)
        converge_sgd <- c(converge_sgd, 1)
        break
      }
    }
  }
  end_time <- Sys.time()
  time_sgd <- end_time - start_time
  convsgd <- mean(converge_sgd)
  
  res <- data.frame(converge_gd = convgd,
                    time_gd = time_gd,
                    convsgd = convsgd,
                    time_sgd = time_sgd)
  
  rdsname <- paste("percep_",id,"_k1w_n10w", sep = "")
  saveRDS(res, rdsname)
}