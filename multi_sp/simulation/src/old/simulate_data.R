## simulate data-set

## nsp: number of species
##
## nsite: number of sites the species might actually occupy
##
## nyr: number of years (note: need to examine the special case of
## nyr=2, which is what Soroye had - this might actually be a
## worst-scase scenario).
##
## nvisit: number of visits within a year
##
## mu.psi:
## sigma.psi.sp:
## mu.psi.yr:
## sigma.psi.yr:
## mu.p:
## sigma.p.sp:
## sigma.p.yr:

## Function to simulate ranges.
##
## Note that, as we include more species, we will be including more
## 'irrelevant sites' for each species, and this is a pattern that
## will likely be true of real data.  The structure of the species
## across space may also matter here (e.g., continents, vs cts
## connected landscape) and is something we could investigate.
##
## For now, the below uses 'sample', so every range size is equally
## likely.  We may want to consider something more realistic, such as
## a log-normal distribution, or explicit space.  This likely matters.
make.ranges <- function(nsp, nsite, type.range) {
  if(type.range == 'equal') {
    nsite.by.sp <- sample.int(n=nsite, size=nsp, replace=TRUE)
  } else if(type.range == 'logn') {
    prop.sites <- rbeta(nsp, shape1=1, shape2=3)
    nsite.by.sp <- round(nsite * prop.sites)
  }
  
  get.sites.within.range <- function(ii) {
    sites <- rep(0,nsite)
    sites[sample(x=1:nsite, size=ii, replace=FALSE)] <- 1
    sites
  }
  res <- t(sapply(nsite.by.sp, get.sites.within.range))
  names(dim(res)) <- c('nsp','nsite')
  res==1
}

make.data <- function(nsp=25,
                      nsite=100,
                      nyr=5,
                      nvisit=3,
                      mu.psi=0,
                      sigma.psi.sp=0.1,
                      mu.psi.yr=0,
                      sigma.psi.yr=0.1,
                      mu.p=0,
                      sigma.p.sp=0.1,
                      p.yr=0,
                      sigma.p.site=0.1,
                      mu.v=-0.5,
                      mu.v.yr=-0.5,
                      type.range='equal',
                      missing.visits=FALSE,
                      sp.range=NULL) {

  ## ------------------------------------------------------------
  ## If species' ranges are not passed in, simulate them for the
  ## specified number of species and sites.
  if(is.null(sp.range)) {
    sp.range <- make.ranges(nsp=nsp,
                            nsite=nsite,
                            type.range=type.range)
  }
  ## ------------------------------------------------------------

  ## ------------------------------------------------------------
  ## specify species-specific occupancy and detection probabilities
  ## and also site visit probabilities
  
  ## species-specific random intercepts
  psi.sp <- rnorm(n=nsp, mean=0, sd=sigma.psi.sp)
  p.sp   <- rnorm(n=nsp, mean=0, sd=sigma.p.sp)

  ## effect of year on occupancy (species-specific random slopes)
  psi.yr <- rnorm(n=nsp, mean=mu.psi.yr, sd=sigma.psi.yr)

  ## effect of site on detection (year-specific)
  p.site <- matrix(rnorm(n=nsite*nyr, mean=0, sd=sigma.p.site),
                   nrow=nsite,
                   ncol=nyr)
  
  ## create empty occupancy, detection, and visitation probability
  ## matrices
  psi.mat <- array(NA, dim=c(nsp=nsp,
                             nsite=nsite,
                             nyr=nyr))
  p.mat <- array(NA, dim=c(nsp=nsp,
                           nsite=nsite,
                           nyr=nyr,
                           nvisit=nvisit))
  v.mat <- array(NA, dim=c(nsite=nsite,
                           nyr=nyr,
                           nvisit=nvisit))

  ## fill in these matrices
  for(site in 1:nsite) {
    for(yr in 1:nyr) {
      for(sp in 1:nsp) {
        psi.mat[sp,site,yr] <- expit(mu.psi +
                                     psi.sp[sp] +
                                     psi.yr[sp]*(yr-1))
      }
      for(visit in 1:nvisit) {
        v.mat[site,yr,visit] <- expit(mu.v + mu.v.yr*(yr-1))
        for(sp in 1:nsp) {
          p.mat[sp,site,yr,visit] <- expit(mu.p +
                                           p.sp[sp] +
                                           p.site[site,yr] +
                                           p.yr*(yr-1))
        }
      }
    }
  }
  ## ------------------------------------------------------------
  
  ## ------------------------------------------------------------
  ## Create occupancy and detectability matrices, etc

  ## To do this, we will construct a range matrix (TRUE for all sites
  ## in range and FALSE for all sites outside of range).  Array
  ## dimensions will also include yr and visit, so that we can use
  ## array multiplication to easily set non-relevant entries to zero.
  range.arr <- array(sp.range, dim=c(nsp=nsp,
                                     nsite=nsite,
                                     nyr=nyr,
                                     nvisit=nvisit))
  
  occ.arr <- array(rbinom(n=length(psi.mat),
                          size=1,
                          prob=psi.mat),
                   dim=c(nsp=nsp,
                         nsite=nsite,
                         nyr=nyr,
                         nvisit=nvisit))
  det.arr <- array(rbinom(n=length(p.mat),
                          size=1,
                          prob=p.mat),
                   dim=c(nsp=nsp,
                         nsite=nsite,
                         nyr=nyr,
                         nvisit=nvisit))
  vis.arr <- array(rbinom(n=length(v.mat),
                          size=1,
                          prob=v.mat),
                   dim=c(nsite=nsite,
                         nyr=nyr,
                         nvisit=nvisit))
  
  ## subset occ.arr down only to sites within each species' range
  occ.arr <- occ.arr*range.arr
  ## subset detection to only sites where the species was present
  det.arr <- det.arr*occ.arr
  ## subset actual detections by incorporating visits
  if(missing.visits) {
    vis.arr.with.sp <- aperm(array(vis.arr, dim=dim(det.arr)[c(2:4,1)]),
                             c(4,1:3))
    det.arr <- det.arr*vis.arr.with.sp
  }
  ## ------------------------------------------------------------
  
  ## ------------------------------------------------------------
  ## create objects to return

  ## true Z matrix is useful for later on when we evaluate performance
  ## of various models
  Z <- (apply(occ.arr, 1:3, sum)>0)*1
  
  ## observation matrix
  X <- det.arr
  
  ## add dimension names to arrays
  arr.names <- list(sp=sprintf("sp_%04d", 1:nsp),
                    site=sprintf("site_%04d", 1:nsite),
                    yr=paste('y',1:nyr,sep='_'),
                    visit=paste('v',1:nvisit,sep='_'))
  dimnames(sp.range) <- arr.names[1:2]
  dimnames(Z)        <- arr.names[1:3]
  dimnames(X)        <- arr.names[1:4]
  
  list(sp.range=sp.range,
       vis.arr=vis.arr,
       Z=Z,
       X=X,
       nsp=nsp,
       nsite=nsite,
       nyr=nyr,
       nvisit=nvisit,
       mu.psi=mu.psi,
       sigma.psi.sp=sigma.psi.sp,
       mu.psi.yr=mu.psi.yr,
       sigma.psi.yr=sigma.psi.yr,
       mu.p=mu.p,
       sigma.p.site=sigma.p.site,
       mu.v=mu.v,
       mu.v.yr=mu.v.yr,
       p.yr=p.yr,
       psi.sp=psi.sp,
       p.sp=p.sp,
       psi.yr=psi.yr,
       p.site=p.site)
  ## ------------------------------------------------------------
}


# Code from Vaughn (needs proper/clean integration). #
library(sf); library(rgeos);
library(concaveman); library(tidyverse)
library(reshape2); library(locfit)

# ----------------------------------------------------
# Function: GenPoly()
# Generates a specified number of species
# ranges with a given complexity which
# specifies the mean number of vertices
# for the set of polygons. Grid size is used to
# generate an artificial landscape on which to place
# randomly generated polygons. Grids must
# be square. Additionally, the number of
# polygons to generate are also specified.
# ----------------------------------------------------
GenPoly <- function(numpoly=100,
                    complexity=50,
                    gridsize=100){
  
  # make a placeholder polygon to create
  # a grid.
  sfc <- st_sfc(st_polygon(list(rbind(c(0,0),
                                      c(1,0),
                                      c(1,1),
                                      c(0,0)))))
  
  # create the grid using the specified grid size.
  grid <- st_make_grid(sfc, cellsize=1/gridsize, square=TRUE)
  
  # create a list to hold all of the generated polygons.
  # iterate while generating the number of specified polygons.
  poly_ls <- list()
  for(i in 1:numpoly){
    
    # draw random number of vertices from normal distribution
    n_vert <- round(rnorm(1, mean=complexity, sd=5))
    if(n_vert >= 3){
      verts <- st_sample(grid, size=n_vert)
      verts <- st_as_sf(verts)
      poly <- concaveman(verts) # generate the concave polygon for these vertices
      
      # randomly scale and translate the polygons to have different size ranges.
      scale_flag <- 1-(0.1*rpois(1, 1))
      # print(paste("Scaling by:", scale_flag))
      
      shift_flag_x <- 0.2*rnorm(1, 0, 1) # translate by x factor
      shift_flag_y <- 0.2*rnorm(1, 0, 1) # translate by y factor
      # print(paste("Translating by:", shift_flag_x, shift_flag_y))
      
      poly <- (poly*scale_flag+c(shift_flag_x, shift_flag_y))$geometry
      poly <- poly %>% st_crop(grid)
    }
    
    poly_ls[[i]] <- poly
  }
  
  return(list(poly_ls, grid))
}

# ----------------------------------------------------
# Function: PolyToMatrix()
# Takes a list of polygons and converts them to a
# binary matrix of grid intersections. 
# A list of matrices is returned. Takes only a 
# large list object from the function GenPoly.
# ----------------------------------------------------
PolyToMatrix <- function(polys){
  matrix_ls <- list()
  n_poly <- length(polys[[1]])
  grid_size <- sqrt(length(polys[[2]]))
  
  for(i in 1:n_poly){
    mtrx <- st_intersects(st_as_sf(polys[[2]]), 
                          st_as_sf(polys[[1]][[i]]))
    matrix_ls[[i]] <- matrix(as.matrix(mtrx), nrow=grid_size, ncol=grid_size)
  }
  
  return(matrix_ls)
}

## simulate data-set

## nsp: number of species
##
## nsite: number of sites the species might actually occupy
##
## nyr: number of years (note: need to examine the special case of
## nyr=2, which is what Soroye had - this might actually be a
## worst-scase scenario).
##
## nvisit: number of visits within a year
##
## mu.psi:
## sigma.psi.sp:
## mu.psi.yr:
## sigma.psi.yr:
## mu.p:
## sigma.p.sp:
## sigma.p.yr:

## Function to simulate ranges.
##
## Note that, as we include more species, we will be including more
## 'irrelevant sites' for each species, and this is a pattern that
## will likely be true of real data.  The structure of the species
## across space may also matter here (e.g., continents, vs cts
## connected landscape) and is something we could investigate.
##
## For now, the below uses 'sample', so every range size is equally
## likely.  We may want to consider something more realistic, such as
## a log-normal distribution, or explicit space.  This likely matters.
make.data <- function(nsp=25,
                      nsite=100,
                      nyr=5,
                      nvisit=3,
                      mu.psi=0,
                      sigma.psi.sp=0.1,
                      mu.psi.yr=0,
                      sigma.psi.yr=0.1,
                      mu.p=0,
                      sigma.p.sp=0.1,
                      p.yr=0,
                      sigma.p.site=0.1,
                      mu.v=-0.5,
                      mu.v.yr=-0.5,
                      type.range='equal',
                      missing.visits=FALSE,
                      sp.range=NULL) {
  
  ## ------------------------------------------------------------
  ## If species' ranges are not passed in, simulate them for the
  ## specified number of species and sites.
  if(is.null(sp.range)) {
    polys <- GenPoly(numpoly=nsp, gridsize=sqrt(nsite))
    sp.range <- array(unlist(PolyToMatrix(polys)), dim=c(nsp, nsite))
  }
  ## ------------------------------------------------------------
  
  ## ------------------------------------------------------------
  ## specify species-specific occupancy and detection probabilities
  ## and also site visit probabilities
  
  ## species-specific random intercepts
  psi.sp <- rnorm(n=nsp, mean=0, sd=sigma.psi.sp)
  p.sp   <- rnorm(n=nsp, mean=0, sd=sigma.p.sp)
  
  ## effect of year on occupancy (species-specific random slopes)
  psi.yr <- rnorm(n=nsp, mean=mu.psi.yr, sd=sigma.psi.yr)
  
  ## effect of site on detection (year-specific)
  p.site <- matrix(rnorm(n=nsite*nyr, mean=0, sd=sigma.p.site),
                   nrow=nsite,
                   ncol=nyr)
  
  ## create empty occupancy, detection, and visitation probability
  ## matrices
  psi.mat <- array(NA, dim=c(nsp=nsp,
                             nsite=nsite,
                             nyr=nyr))
  p.mat <- array(NA, dim=c(nsp=nsp,
                           nsite=nsite,
                           nyr=nyr,
                           nvisit=nvisit))
  v.mat <- array(NA, dim=c(nsite=nsite,
                           nyr=nyr,
                           nvisit=nvisit))
  
  ## fill in these matrices
  for(site in 1:nsite) {
    for(yr in 1:nyr) {
      for(sp in 1:nsp) {
        psi.mat[sp,site,yr] <- expit(mu.psi +
                                       psi.sp[sp] +
                                       psi.yr[sp]*(yr-1))
      }
      for(visit in 1:nvisit) {
        v.mat[site,yr,visit] <- expit(mu.v + mu.v.yr*(yr-1))
        for(sp in 1:nsp) {
          p.mat[sp,site,yr,visit] <- expit(mu.p +
                                             p.sp[sp] +
                                             p.site[site,yr] +
                                             p.yr*(yr-1))
        }
      }
    }
  }
  ## ------------------------------------------------------------
  
  ## ------------------------------------------------------------
  ## Create occupancy and detectability matrices, etc
  
  ## To do this, we will construct a range matrix (TRUE for all sites
  ## in range and FALSE for all sites outside of range).  Array
  ## dimensions will also include yr and visit, so that we can use
  ## array multiplication to easily set non-relevant entries to zero.
  range.arr <- array(sp.range, dim=c(nsp=nsp,
                                     nsite=nsite,
                                     nyr=nyr,
                                     nvisit=nvisit))
  
  occ.arr <- array(rbinom(n=length(psi.mat),
                          size=1,
                          prob=psi.mat),
                   dim=c(nsp=nsp,
                         nsite=nsite,
                         nyr=nyr,
                         nvisit=nvisit))
  det.arr <- array(rbinom(n=length(p.mat),
                          size=1,
                          prob=p.mat),
                   dim=c(nsp=nsp,
                         nsite=nsite,
                         nyr=nyr,
                         nvisit=nvisit))
  vis.arr <- array(rbinom(n=length(v.mat),
                          size=1,
                          prob=v.mat),
                   dim=c(nsite=nsite,
                         nyr=nyr,
                         nvisit=nvisit))
  
  ## subset occ.arr down only to sites within each species' range
  occ.arr <- occ.arr*range.arr
  ## subset detection to only sites where the species was present
  det.arr <- det.arr*occ.arr
  ## subset actual detections by incorporating visits
  if(missing.visits) {
    vis.arr.with.sp <- aperm(array(vis.arr, dim=dim(det.arr)[c(2:4,1)]),
                             c(4,1:3))
    det.arr <- det.arr*vis.arr.with.sp
  }
  ## ------------------------------------------------------------
  
  ## ------------------------------------------------------------
  ## create objects to return
  
  ## true Z matrix is useful for later on when we evaluate performance
  ## of various models
  Z <- (apply(occ.arr, 1:3, sum)>0)*1
  
  ## observation matrix
  X <- det.arr
  
  ## add dimension names to arrays
  arr.names <- list(sp=sprintf("sp_%04d", 1:nsp),
                    site=sprintf("site_%04d", 1:nsite),
                    yr=paste('y',1:nyr,sep='_'),
                    visit=paste('v',1:nvisit,sep='_'))
  dimnames(sp.range) <- arr.names[1:2]
  dimnames(Z)        <- arr.names[1:3]
  dimnames(X)        <- arr.names[1:4]
  
  list(sp.range=sp.range,
       vis.arr=vis.arr,
       Z=Z,
       X=X,
       nsp=nsp,
       nsite=nsite,
       nyr=nyr,
       nvisit=nvisit,
       mu.psi=mu.psi,
       sigma.psi.sp=sigma.psi.sp,
       mu.psi.yr=mu.psi.yr,
       sigma.psi.yr=sigma.psi.yr,
       mu.p=mu.p,
       sigma.p.site=sigma.p.site,
       mu.v=mu.v,
       mu.v.yr=mu.v.yr,
       p.yr=p.yr,
       psi.sp=psi.sp,
       p.sp=p.sp,
       psi.yr=psi.yr,
       p.site=p.site)
  ## ------------------------------------------------------------
}

