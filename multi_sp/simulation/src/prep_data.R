prep.data <- function(dd, limit.to.visits, limit.to.range) {
  
  ## keep only detected species:
  sp.keep <- which(apply(dd$Z, 1, sum, na.rm=TRUE)>0)
  
  ## which sites to keep will depend on limit.to.visits
  ## keep all sites, even those without any detections
  if(limit.to.visits=='all') {
      site.keep <- 1:dd$nsite
  }
  ## keep only sites that were visited
  if(limit.to.visits=='visits') {
      site.keep <- which(apply(dd$vis.arr, 2, sum, na.rm = TRUE)>0)
  }
  ## keep only sites that yielded a detection of at least one species
  if(limit.to.visits=='detected') {
      site.keep <- which(apply(dd$X, 'site', sum, na.rm = TRUE)>0)
  }
  if(limit.to.visits=='community') {
      site.keep <- which(apply(dd$X, 'site', sum, na.rm = TRUE)>0)
  }
  
  ## subset based on the above
  dd$Z        <- dd$Z[sp.keep,site.keep,,drop=FALSE]
  dd$X        <- dd$X[sp.keep,site.keep,,,drop=FALSE]
  dd$sp.range <- dd$sp.range[sp.keep,site.keep,drop=FALSE]
  dd$vis.arr  <- dd$vis.arr[sp.keep,site.keep,,,drop=FALSE]
  dd$nsp      <- length(sp.keep)
  dd$nsite    <- length(site.keep)
  dd$vis.col <- dd$vis.col[site.keep,,,drop=FALSE]
  
  ## generate master index (to improve model efficiency (this prevents
  ## unnecessary iterating through all irrelevant sites and visits)
  get.indices <- function(sp) {
    ## visited array
    vis.arr <- dd$vis.arr[sp,,,]
    ## if modelling all visits, set visit array to 1 everywhere
    if(limit.to.visits=='all') 
      vis.arr[TRUE] <- 1
    ## if modelling sites with detections only, create new visit array
    if(limit.to.visits=='detected') {
      nsp.detected <- apply(dd$X, 2:4, sum, na.rm = TRUE)
      vis.arr[TRUE] <- 1
      vis.arr[nsp.detected==0] <- 0
    }
    ## if inferring non-detections only for community visits 
    if(limit.to.visits=='community') {
      
      # figure out sites where at least 1 species was detected
      nsp.detected <- apply(dd$X, 2:4, sum, na.rm = TRUE)
      # multiply by sites where we know it was a community collection
      sp.det.com.col <- nsp.detected*dd$vis.col
      
      vis.arr[TRUE] <- 0
      vis.arr[sp.det.com.col!=0] <- 1  ## visits occurred where at least 1 species was detected and there was a community collection
      vis.arr[dd$X[sp,,,]==1] <- 1 ## visits for sp also occurred where each species was found 
    }
    
    if(limit.to.range=='yes'){ # restrict to within species' ranges
      vis.arr[!dd$sp.range[sp,],,] <- 0
      tmp <- which(vis.arr==1, arr.ind=TRUE)
      indices <- cbind(rep(sp,nrow(tmp)),tmp)
    }
    if(limit.to.range=='no'){ # unrestricted with respect to species' ranges
      vis.arr[!dd$sp.range[sp,],,] <- 1
      tmp <- which(vis.arr==1, arr.ind=TRUE)
      indices <- cbind(rep(sp,nrow(tmp)), tmp)
    }
    
    return(indices)

  }
  master.index <- do.call(rbind, lapply(1:dd$nsp, get.indices))
  colnames(master.index) <- c('sp','site','yr','visit')
  
  ## data structures to be returned
  
  my.data <- list(X=dd$X[master.index])
  
  my.constants <- list(
    nsp=dim(dd$X)['nsp'],
    nsite=dim(dd$X)['nsite'],
    nyr=dim(dd$X)['nyr'],
    nind=nrow(master.index),
    yrv=master.index[,'yr'],
    sitev=master.index[,'site'],
    spv=master.index[,'sp'])
  
  return(list(my.constants = my.constants, my.data = my.data))
}
