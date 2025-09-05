combine_sumstats <- function(sumstats){
  ngroups <- length(sumstats)
  
  xx <- sumstats[[1]]$xx
  xy <- sumstats[[1]]$xy
  
  nobs <- attr(sumstats[[1]]$xx, "nobs")
  colsum <- attr(sumstats[[1]]$xx, "colsum")  
  ysum <- attr(sumstats[[1]]$xy, "ysum")
  yssq <- attr(sumstats[[1]]$xy, "yssq")
  
  for(index in 2:ngroups){
    xx <- xx + sumstats[[index]]$xx
    xy <- xy + sumstats[[index]]$xy
    nobs <- nobs + attr(sumstats[[index]]$xx, "nobs")
    colsum <- colsum + attr(sumstats[[index]]$xx, "colsum")  
    ysum <- ysum + attr(sumstats[[index]]$xy, "ysum")
    yssq <- yssq + attr(sumstats[[index]]$xy, "yssq")
  }
  
  ## center and scale xx
  xx <- xx - colsum %o% (colsum/nobs)
  ## scale xx matrix (x'x) so that diag = 1 and off-diag serve as correlations of cols of x's
  sdx <- sqrt(diag(xx)/nobs)
  xx <- xx / (sdx %o% sdx)
  ## divide by n because loss has 1/n
  xx <- xx / nobs
  
  ## center and scale xy
  yssq <- yssq - ysum^2/nobs
  yvar <- yssq / nobs
  sdy <- sqrt(yvar)
  xy <- xy - colsum * (ysum/nobs)
  xy <- xy/ (sdx * sdy)
  ## divide by n because loss has 1/n
  xy <- xy / nobs
  
  ## after standardize y
  yvar <- yvar/yvar
  
  ## to put beta's on original scale:
  ## if   y = x * b
  ## then y/sdy = x/sdx * a
  ## and b = a * (sdy/sdx)
  
  beta_multiplier <- sdy / sdx
  

  ## need to remove attr of xx and xy
  ## because they were copied from sumstats[[1]]
  ## and not relevant after summing across all
  ## sumstats
  
  attr(xx, "nsubj") <- NULL
  attr(xx, "nmiss") <- NULL
  attr(xx, "nobs") <-  NULL
  attr(xx, "colsum") <- NULL
  
  attr(xy, "ysum") <-  NULL
  attr(xy, "yssq") <-  NULL
  attr(xy, "nsubj") <- NULL
  attr(xy, "nmiss") <- NULL
  attr(xy, "nobs") <- NULL
  
  
  return(list(xx=xx, xy=xy, yvar=yvar, ysum=ysum, nobs=nobs,
              beta_multiplier=beta_multiplier))
  
}

