
## summary stat functions for glmnet with summary stats

make_sumstats <- function(x, y){
  ## important: x and y must have the same number of subjects
  ## and ordered the same way according to subject ID, but ID 
  ## should not be a column of X
  
  ## create sum stats xx and xy 
  xcol.names <-  colnames(x)
  names.miss <- xcol.names == "" | is.null(xcol.names)
  if(any(names.miss)){
    stop("col names of x are missing")
  }
  
  nsubj <- nrow(x)
  is.miss.x <- apply(is.na(x), 1, any)
  is.miss.y <-  is.na(y)
  is.miss <- is.miss.x | is.miss.y
  z <- scale(x[!is.miss, ], scale=FALSE)
  xx <- (t(z) %*% z)
  rownames(xx) <- xcol.names
  colnames(xx) <- xcol.names
  
  y <- y[!is.miss]
  ysum <- sum(y, na.rm=TRUE)
  y <- scale(y, scale=FALSE)
  xy <- t(z) %*% y
  names(xy) <-  xcol.names
  
  attr(xx, "nsubj") <- nsubj
  attr(xx, "nmiss") <- sum(is.miss)
  
  attr(xy, "ysum") <-  ysum
  attr(xy, "yssq") <-  sum(y^2, na.rm=TRUE)
  attr(xy, "nsubj") <- nsubj
  attr(xy, "nmiss") <- sum(is.miss)
  
  return(list(xx=xx, xy=xy))
}

make_varxy <- function(xx, xy, yssq, nsubj){
  ## create var mat for later simulations
  vx <- xx
  vxy <- xy
  v12 <- c(yssq, vxy)
  vmat <- rbind(v12, cbind(vxy, vx))
  vmat <- vmat/ nsubj
  return(vmat)
}


#' Compute elastic net from summary statistics
#' @description
#' Fit an elastic net model based on summary statistics for specified 
#' penalty parameters alpha and lambda
#' 
#' @param xx matrix of x'x where x is nxp design matrix
#' @param xy vector of x'y where y is nx1 trait vector
#' @param beta_init is px1 vector of starting values
#' @param alpha (range 0-1) is fraction of penalty for L1 penalty
#' and (1-alpha) is fraction for L2 penalty
#' @param lambda is penalty parameter
#' @param penalty_factor is px1 vector for weighting penalties of design matrix columns.
#' Typically used as value of 0 for terms not to be penalized and 1 for terms to penalize
#' @param maxiter is maximum number of cyclic coordinate descent iterations.
#' One iteration is over all p parameters.
#' @param tol is tolerance to check convergence
#' @param f.adj (f.adj >=1) is an adjustment factor to guard against diverging betas, which can 
#' occur when n is small relative to p. The default value of 2 should work most of the time.
#' If beta parameters are found to diverge, larger values can be used.
#' Note that larger values of f.adj cause more iterations.
#' See Yang, Y., & Zou, H. (2012). A coordinate majorization descent algorithm for L1 penalized learning. 
#' Journal of Statistical Computation and Simulation, 84(1), 84–95. https://doi.org/10.1080/00949655.2012.695374 
#' @param verbose prints summary results at each interation if verbose=TRUE 
#' @details
#' Cyclic coordinate descent is used to fit an elastic net model based on minimizing penalized least squared. 
#' @returns list with fitted beta, number of iterations, convegence (true/false), and input penalty parameters
#' alpha and lambda
#' @author Dan Schaid (schaid@mayo.edu)
#' @export
glmnet_sumstats <- function(xx, xy, beta_init, alpha, lambda, penalty_factor,obj_check=FALSE, 
                            maxiter=500, tol=1e-5, f.adj=2, verbose=FALSE){
  eps <- 1e-3 ## used to check whether a beta is penalized by penalty_factor
  
  lambda_pen <- lambda * penalty_factor
  np <- ncol(xx)
  if(is.null(beta_init)){
    beta_init <- rep(0, np)
  }
  converge <- FALSE
  betaold <- betanew <- beta_init
  
  floss_old <-  elastic_net_loss(betaold, xx, xy, alpha,lambda_pen)
  active <- rep(FALSE, np)
  
  
  ## determine active set
  for(j in 1:np){
    diff <- xy[j] - xx[j,] %*% betaold
    if(abs(diff) > lambda_pen[j]*alpha)
    {
      active[j] <- TRUE
    }
  }
  
  ## iterate over active set
  for(iter in 1:maxiter){
    
    ## random order of beta updates
    index <- sample(1:np, np)
    for(i in 1:np){
      j <- index[i]
      if(!active[j]) next
      
      
      if(abs(lambda_pen[j]) < eps){
        betanew[j] <- (xy[j] - xx[j,] %*% betaold + xx[j,j]*betaold[j])/xx[j,j]
      } else{
        u <- (xy[j] - xx[j,] %*% betaold +  f.adj*xx[j,j]*betaold[j])
        betanew[j] <- soft(u, alpha*lambda_pen[j])/(f.adj * xx[j,j] + lambda_pen[j]*(1-alpha))
      }
    }
   
    floss_new <- elastic_net_loss(betanew, xx, xy, alpha,lambda_pen)
    fdelta <- floss_new - floss_old
    
    if(verbose){
      cat("active set: iter = ", iter,", f.adj = ", f.adj, ", f = ", floss_new, "fdelta = ", fdelta , ", beta range = ", range(betanew),"\n")
    }
    
    if(abs(fdelta) < tol){
      converge <- TRUE
      break
    }
    
    floss_old <- floss_new
    betaold <- betanew
    
  }
  
  ## now check all beta's, not just active
  if(converge){
    
    converge <- FALSE
    
    for(iter in 1:maxiter){
      
      ## random order of beta updates
      index <- sample(1:np, np)
      for(i in 1:np){
        j <- index[i]
        
        if(abs(lambda_pen[j]) < eps){
          betanew[j] <- (xy[j] - xx[j,] %*% betaold + xx[j,j]*betaold[j])/xx[j,j]
        } else{
          u <- (xy[j] - xx[j,] %*% betaold +  f.adj*xx[j,j]*betaold[j])
          betanew[j] <- soft(u, alpha*lambda_pen[j])/(f.adj * xx[j,j] + lambda_pen[j]*(1-alpha))
        }
      }
      
     
      floss_new <- elastic_net_loss(betanew, xx, xy, alpha,lambda_pen)
      
      fdelta <- floss_new - floss_old
    
      
      if(verbose){
        cat("beta check: iter = ", iter,", f.adj = ", f.adj, ", f = ", floss_new, ", fdelta = ", fdelta , ", beta range = ", range(betanew),"\n")
      }
      
      if(abs(fdelta) < tol){
        converge <- TRUE
        break
      }
      
      floss_old <- floss_new
      betaold <- betanew
    }
    
  }
  
  return(list(beta=betanew, iter=iter, converge=converge, alpha=alpha, lambda=lambda))
}

soft <- function(x, gamma){
  if(x > gamma){
    s <- x - gamma
  } else if (x < -gamma){
    s <- x + gamma
  }else{
    s <- 0
  }
  return(s)
}

elastic_net_loss <- function(beta, xx, xy, alpha,lambda_pen){
  f <-  sum( t(beta) %*% xx %*% beta - 2*t(beta) %*% xy) + alpha*sum(lambda_pen * abs(beta)) +
    (1-alpha)*sum(lambda_pen*beta^2)
  return(f)
}

gradient <- function(index, beta, xx, xy, alpha, lambda_pen){
  d <- 2*(xx %*% beta)[index] -2*xy[index]  + 
    alpha*lambda_pen[index]*sign(beta[index]) + 2*(1-alpha)*lambda_pen[index]*beta[index]
}




auc_glmnet_sumstats <- function(beta, xx, vary, nsubj,  ncase, ncont){
  ssr <- t(beta) %*% xx %*% beta
  sst <- vary * nsubj
  r2 <- ssr/sst
  a <- (ncase + ncont)^2/(ncase * ncont)
  d <- sqrt(a*r2)/sqrt(1-r2)
  auc <- pnorm(d/sqrt(2))
  return(list(auc=auc, r2=r2))
}

sim_test_dat <- function(nsubj, nprs, prev=.1, beta.sd){
 ## large beta.sd allows larger beta's
  sex <-  rbinom(nsubj,size=1,prob=.5)
  age <- round(runif(nsubj, min=40, max=70),0)
  cov1 <- rnorm(nsubj)
  cov2 <- rnorm(nsubj)
  x <- cbind(age, sex, cov1, cov2, matrix(rnorm(nsubj*nprs), nrow=nsubj))
  beta <- rnorm(nprs+4, sd=beta.sd)
  xb <- x %*% beta + log(prev/(1-prev))
  p <- exp(xb) / (1 + exp(xb))
  pheno <- 1*(runif(nsubj) <= p)
  colnames(x) <- c("age","sex","cov1","cov2", paste0("PRS00", 1:nprs))
  return(list(pheno=pheno, x=x))
}

