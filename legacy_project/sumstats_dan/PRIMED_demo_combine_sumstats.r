## create test data of 3 groups with PRS (nprs) and covariates (sex, age, cov1, cov2)

source("../PRSmix/PRIMED_summary_stats.R")

## simulate training data

nprs <- 4000

n1 <- 1000
n2 <- 500
n3 <- 50

dat1 <- sim_test_dat(n1, nprs, prev=.1, beta.sd=.01)
dat2 <- sim_test_dat(n2, nprs, prev=.1, beta.sd=.01)
dat3 <- sim_test_dat(n3, nprs, prev=.1, beta.sd=.01)


x_test <- dat1$x[1:5,]
y_test <- dat1$pheno[1:5]


## make some pheno's missing
dat1$pheno[3:20] <- NA
dat2$pheno[15:30] <- NA
dat3$pheno[1:16] <- NA

## make some x's missing
dat1$x[3:5, 4:5] <- NA
dat2$x[20:30, 5:8] <- NA
dat3$x[2:10, 6:9] <- NA

############ Biobanks Create Summary Statistics (3 demo biobanks) ###########

## need to discuss:
##
##  1) add cluster to make_sumstats to make xx,xy for all clusters and return
##  as list ?
##
##  2) consider adding 10-fold xx/xy subsets for xval ? (not always feasible)
##
##  3) need to exlude AoU PGS
##
##  4) need to consider how to merge (sum) xx,xy across biobanks and
##  account for missing PGS from some biobanks. One possibility is
##  to add row/col of 0's for missing PRS, and then sum across 
##  biobanks

sumstats <- list()
sumstats[[1]] <- make_sumstats_V2(dat1$x, dat1$pheno)
sumstats[[2]] <- make_sumstats_V2(dat2$x, dat2$pheno)
sumstats[[3]] <- make_sumstats_V2(dat3$x, dat3$pheno)

  
###### CC: sum over biobanks ################################################

total_stats <- combine_sumstats(sumstats)

xx <- total_stats$xx
xy <- total_stats$xy

yvar<- total_stats$yvar
ysum <- total_stats$ysum
nobs <- total_stats$nobs

## vmat used later for simulations to create simulated validation data
## used to choose lambda that gives best fit
## after standardizing x and y, vmat is a correlation matrix of y and x
vmat <- make_varxy(xx, xy, yvar)

## note that there is no intercept because x and y are centered on their means

penalty_factor <- rep(1,ncol(xx))
## if don't penalize adjusting covariates: should they be penalized??
## penalty_factor[1:4] <- 0

############### fit a grid of lambda's ##################
library(glmnet)
## alpha weights abs(beta) and (1-alpha) weights beta^2
alpha <- 0.5
lambda_max <- max(abs(xy))/alpha
lambda_frac <- seq(from=1, to=.05, by= -.05)
lambda_vec <- lambda_max * lambda_frac
nfit <- length(lambda_vec)

beta_init <- rep(0, ncol(xx))
fit_sumstats <- list()

cat("================ lambda.frac = ", lambda_frac[1], " ==================\n")

fit_sumstats[[1]] <-  glmnet_sumstats(xx, xy, beta_init, alpha=alpha, lambda=lambda_vec[1], 
                                penalty_factor, maxiter=500, tol=1e-7,f.adj=2.0, verbose=TRUE)

for(j in 2:nfit){
  ptm <- proc.time()
  cat("================ lambda.frac = ", lambda_frac[j], " ==================\n")
  ## use warm start for next fit
  beta_init <- fit_sumstats[[j-1]]$beta
  fit_sumstats[[j]] <-  glmnet_sumstats(xx, xy, beta_init, alpha=alpha, lambda=lambda_vec[j], 
                                    penalty_factor, maxiter=500, tol=1e-7,f.adj=32.0, verbose=TRUE)
  print(proc.time() - ptm)
}



##### approximate AUC and R2 for case/control

ncase <- ysum
ncont <- nobs - ncase

ncase; ncont
auc <- nbeta <- r2 <- rep(0, nfit)
for(i in 1:nfit){
  tmp <- auc_glmnet_sumstats(fit_sumstats[[i]]$beta, xx, yvar, ncase, ncont)
  auc[i] <- tmp$auc
  r2[i]  <- tmp$r2
  nbeta[i] <- sum( abs(fit_sumstats[[i]]$beta*penalty_factor) > 1e-6)
} 

plot(nbeta, r2)
plot(lambda_frac, auc)
plot(lambda_frac, nbeta)
max(auc)

loss <- rep(0, nfit)

for (i in 1:nfit) {
  beta <- as.vector(fit_sumstats[[i]]$beta)
  loss[i] <- yvar- 2*t(beta) %*% xy +  t(beta) %*% xx %*% beta
}

plot(lambda_frac, loss)
plot(nbeta, loss)
index <- (1:nfit)[loss == min(loss)]
beta <- as.vector(fit_sumstats[[index]]$beta)
hist(beta)
nbeta[index]

################ choose lambda penalty #################################
## now choose best lambda based on simulated data, using
## the covariance matrix of Y with X to simulate validation data

## not sure if this is best, or consider alternative with cross-validation
## with setup from make_sumstats to create xval data sets: need to discuss

library(MASS)

set.seed(123)

simulate_yx <- function(n, vmat) {
  ## simulate y and X matrix from covariance
  ## matrix:
  ##        | var(y),  cov(y,x) |
  ## vmat = | cov(x,y) cov(x,x) |
  ## and return matrix [Y|X]
  p <- ncol(vmat)
  svd_decomp <- svd(vmat)
  U <- svd_decomp$u
  D <- diag(sqrt(pmax(svd_decomp$d, 0)))
  Z <- matrix(rnorm(n * p), nrow = n)
  return(Z %*% U %*% D)
}

nsim <- 5
nfit <- length(fit_sumstats)

## compute correlation of predicted using beta from training data 
## with observed validation data (simulated)
## see equation (12) p 472 of Mak et al.
cormat <- lossmat <- matrix(0, nrow = nsim, ncol = nfit)
cr <- rep(0, nfit)
loss <- rep(0, nfit)

for(isim in 1:nsim) {
  yx <- simulate_yx(nobs, vmat)
  y.sim <- yx[, 1]
  x.sim <- yx[, -1]
   
  ## center y
  y.c <- y.sim - mean(y.sim)
  ## center and scale x's
  x.c <- scale(x.sim) / sqrt(nobs - 1)
  
  for (i in 1:nfit) {
    beta <- as.vector(fit_sumstats[[i]]$beta)
    xb <- (x.c %*% beta)
    num <- t(xb) %*% y.c
    den <- as.vector(sqrt((t(xb) %*% xb) * t(y.c) %*% y.c))
    if(abs(den) < 1e-6){
      cr[i] <- 0
    } else{
      cr[i] <- as.numeric(num / den)
    }
    loss[i] <- sum( (y.sim - xb)^2)
  }
  
  cormat[isim, ] <- cr
  lossmat[isim,] <- loss
}


cor_mean <- apply(cormat, 2, mean)
loss_mean <- apply(lossmat, 2, mean)
plot(lambda_vec, loss_mean)
plot(lambda_vec, cor_mean)

index_best <- (1:nfit)[loss_mean == min(loss_mean)][1]
index_best
loss_mean
lambda_vec[index_best]

## now look at results from best fit
beta <- as.vector(fit_sumstats[[index_best]]$beta)
hist(beta)
range(beta)
table(abs(beta) > 1e-6)
auc[index_best]
plot(rep(lambda_vec, nsim), as.vector(lossmat))

