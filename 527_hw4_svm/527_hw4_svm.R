

## Algorithm of SVM -----
library(mvtnorm)
n=10
X = rmvnorm(n, mean=c(7,7), sigma=matrix(c(1,0.5,0.5,1), ncol = 2))
x1 = X[,1]
x2 = X[,2]
X = rmvnorm(n, mean=c(5,5), sigma=matrix(c(1,0.5,0.5,1), ncol = 2))
x11 = X[,1]
x22 = X[,2]

p = plot(x1,x2,col="red",xlim=c(-5,20),ylim=c(-5,20))
    + points(x11,x22,col="blue")

X <- matrix(c(x1,x11,x2,x22),ncol=2)

Y <- matrix(rep(c(1,-1),each=n), ncol=1)

X <- scale(X)

Ymat <- Y %*% t(Y)
Xmat <- X %*% t(X)


library(quadprog)
library(Matrix)
Dmat       <- Ymat * Xmat
Dmatnear <- nearPD(Dmat)
Dp <- Dmatnear$mat

# dvec       <- t(matrix(1,2*n,1))
dvec       <- as.vector(matrix(1,2*n,1))
firstrow   <- t(Y)
diagmat <- diag(x = 2*n)
diagmat2 <- -diag(x = 2*n)
Amat_t <- rbind(firstrow, diagmat, diagmat2)
Amat <- t(Amat_t)
cval <- 10
# bvec <- t(as.matrix(c(rep(0,(2*n+1)),rep(-cval,(2*n)))))
bvec <- c(rep(0,(2*n+1)),rep(-cval,(2*n)))
res <- solve.QP(Dmat=Dmatnear$mat,dvec=dvec,Amat=Amat,bvec=bvec, meq=1)

alpha <- res$solution
betas <- t(as.matrix(alpha*Y)) %*% as.matrix(X);betas

alpha = matrix(alpha, nrow=n)
nonzero =  abs(alpha)>1e-4
beta = rowSums(sapply(which(nonzero), function(i) alpha[i]*Y[i]*X[i,]))
X_sv = X[nonzero,]
Y_sv = Y[nonzero]
beta0 <- -0.5*(min(X_sv[Y_sv==1,]%*%beta)+max(X_sv[Y_sv==-1,]%*%beta))
beta

library(e1071)
## Compared to xvm results
svm.model <- svm(as.factor(Y) ~ X, cost = 10, kernel="linear",scale=TRUE)

beta = drop(t(svm.model$coefs)%*%X[svm.model$index,])
beta0 = svm.model$rho
beta;beta0










