
## load data -----

gene <- read.delim("GEUVADIS_normalized_expression_chr20")

genotype <- read.delim("GEUVADIS_chr20_processed.traw")

## data set construction -----
### suppose we consider gene id ENSG00000000419 only
idx = 1
y_pre <- gene[idx,]
start <- y_pre$start - 500000
end <- y_pre$end + 500000
y <- t(y_pre[,5:362])

### Filter x variables with +/- 500,000 basepairs
x.idx <- which(genotype$POS >= start & genotype$POS <= end)
x <- t(genotype[x.idx,7:364])
colnames(x) <- genotype$SNP[x.idx]


## Build my own cyclic coordinate descent algorithm for lasso -----

soft_threshold <- function(u,v,lambda,n,intercept=FALSE){
  if (intercept == TRUE & j ==1 ){
    eqn1 <- (t(u)%*%v)/n
    s_lambda <- eqn1
  } else {
    eqn1 <- (t(u)%*%v)/n
    tmp <- abs(eqn1)-lambda
    eqn2 <- ifelse(tmp>0, tmp, 0)
    s_lambda <- sign(eqn1)*eqn2
  }
  return(s_lambda)
}

ccd_function <- function(x,y,lambda,threshold=1E-7,intercept=FALSE){
  ## inputs of x and y do not need to be standardized
  
  ## scale x and y
  xs <- scale(x)
  if (intercept == TRUE){
    xs <- cbind(1, xs) 
  }
  ys <- scale(y)
  
  nvar = ncol(xs) # number of variables
  nsubj = nrow(xs) # number of subjects
  
  ## set initial values of beta
  ## first, do not consider intercept 
  bmat <- matrix(1, nrow=nvar)
  err <- matrix(5,nrow=nvar)
  i =1
  # while(any(err>threshold) == TRUE){  
  while (i <=2000){ # 1000 is enough. The results at the end are calculated with 2000 repetation.
    start_time <- Sys.time()
    for (j in 1:nvar){
      
      bpre = bmat[j]
      xj = as.matrix(xs[,j])
      xk = xs[,-j]
      bk = bmat[-j]
      u = xj
      v = ys - xk %*% bk
      
      bpost <- soft_threshold(u=u,v=v,lambda=lambda,n=nsubj,
                              intercept=intercept)
      bmat[j] <- bpost
      j = j+1
    }
    end_time <- Sys.time()
    print(paste("Iteration", i, ", time:",end_time - start_time, sep = ""))
    i = i+1
  }
  
  bmat
  
}

## With intercept -----
### select a lambda = 0.1
bmat_withint <- ccd_function(x = x, y = y, lambda=0.1, intercept = TRUE)
bmat_withint <- readRDS("bmat_withint.rds")
x_withint <- cbind(1,x)
colnames(x_withint)[1] <- "(Intercept)"
int1 <- colnames(x_withint)[which(bmat_withint!=0)]
bmat_withint[which(bmat_withint!=0)]

# > colnames(x_withint)[which(bmat_withint!=0)]
# [1] "(Intercept)"               "20_49064671_G_A_b37"     "20_49078164_C_CAG_b37"  
# [4] "20_49101662_A_G_b37"     "20_49374183_G_GA_b37"    "20_49424564_C_T_b37"    
# [7] "20_49647371_C_T_b37"     "20_49693755_T_C_b37"     "20_49710044_G_A_b37"    
# [10] "20_49714264_G_A_b37"     "20_49715109_C_T_b37"     "20_49936124_G_A_b37"    
# [13] "20_49938899_G_A_b37"     "20_49974518_G_GTGGA_b37" "20_50064221_G_A_b37"    
# > bmat_withint[which(bmat_withint!=0)]
# [1]  5.944965e-16  3.125334e-02  5.602805e-03  2.781656e-03 -1.678444e-02
# [6] -6.548586e-03 -9.462111e-04  1.325202e-02 -2.532187e-02 -1.206171e-02
# [11] -9.559900e-03 -8.233062e-03  4.316459e-02 -4.768447e-02  2.444805e-02


### compared with glmnet
fit_withint <- glmnet(scale(x),scale(y),intercept = TRUE, lambda=0.1,
                      standardize = FALSE)
mycoef <- coef(fit_withint)
int2 <- rownames(mycoef)[which(mycoef != 0)]
mycoef[which(mycoef != 0)]

# > rownames(mycoef)[which(mycoef != 0)]
# [1] "(Intercept)"             "20_49064671_G_A_b37"     "20_49078164_C_CAG_b37"  
# [4] "20_49101662_A_G_b37"     "20_49374183_G_GA_b37"    "20_49424564_C_T_b37"    
# [7] "20_49647371_C_T_b37"     "20_49693755_T_C_b37"     "20_49710044_G_A_b37"    
# [10] "20_49714264_G_A_b37"     "20_49715109_C_T_b37"     "20_49936124_G_A_b37"    
# [13] "20_49938899_G_A_b37"     "20_49974518_G_GTGGA_b37" "20_50064221_G_A_b37"    
# > mycoef[which(mycoef != 0)]
# [1]  5.715633e-16  3.134430e-02  5.604514e-03  2.778403e-03 -1.680209e-02
# [6] -6.551715e-03 -9.215823e-04  1.327820e-02 -2.540688e-02 -1.205019e-02
# [11] -9.607439e-03 -8.240967e-03  4.330458e-02 -4.784740e-02  2.451292e-02

int1 == int2
# > int1 == int2
# [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

# With intercept, the results of my cyclic coordinate descent match the glmnet results.

## Without intercept -----
### select a lambda = 0.1
bmat_noint <- readRDS("bmat_noint.rds")
noint1 <- colnames(x)[which(bmat_noint!=0)]
bmat_noint[which(bmat_noint!=0)]
# > colnames(x)[which(bmat_noint!=0)]
# [1] "20_49064671_G_A_b37"     "20_49078164_C_CAG_b37"   "20_49101662_A_G_b37"    
# [4] "20_49374183_G_GA_b37"    "20_49424564_C_T_b37"     "20_49647371_C_T_b37"    
# [7] "20_49693755_T_C_b37"     "20_49710044_G_A_b37"     "20_49714264_G_A_b37"    
# [10] "20_49715109_C_T_b37"     "20_49936124_G_A_b37"     "20_49938899_G_A_b37"    
# [13] "20_49974518_G_GTGGA_b37" "20_50064221_G_A_b37"    
# > bmat_noint[which(bmat_noint!=0)]
# [1]  0.0312533426  0.0056028048  0.0027816559 -0.0167844356 -0.0065485858
# [6] -0.0009462111  0.0132520156 -0.0253218738 -0.0120617099 -0.0095599004
# [11] -0.0082330625  0.0431645887 -0.0476844724  0.0244480452

### compared with glmnet
fit_noint <- glmnet(scale(x),scale(y),intercept = FALSE, lambda=0.1,
                    standardize = FALSE)
mycoef <- coef(fit_noint)
noint2 <- rownames(mycoef)[which(mycoef != 0)]
mycoef[which(mycoef != 0)]

# > rownames(mycoef)[which(mycoef != 0)]
# [1] "20_49064671_G_A_b37"     "20_49078164_C_CAG_b37"   "20_49101662_A_G_b37"    
# [4] "20_49374183_G_GA_b37"    "20_49424564_C_T_b37"     "20_49647371_C_T_b37"    
# [7] "20_49693755_T_C_b37"     "20_49710044_G_A_b37"     "20_49714264_G_A_b37"    
# [10] "20_49715109_C_T_b37"     "20_49936124_G_A_b37"     "20_49938899_G_A_b37"    
# [13] "20_49974518_G_GTGGA_b37" "20_50064221_G_A_b37"    
# > mycoef[which(mycoef != 0)]
# [1]  0.0313443015  0.0056045137  0.0027784027 -0.0168020935 -0.0065517145
# [6] -0.0009215823  0.0132782017 -0.0254068816 -0.0120501940 -0.0096074391
# [11] -0.0082409674  0.0433045782 -0.0478473965  0.0245129233

noint1 == noint2
# > noint1 == noint2
# [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

# Without intercept, the results of my cyclic coordinate descent match the glmnet results.
