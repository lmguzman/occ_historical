## ************************************************************
## simulate data for dynamic model 
## ************************************************************
library(stringr)
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
  i <- 1
  while(i <= numpoly){
    
    # draw random number of vertices from normal distribution
    n_vert <- round(rnorm(1, mean=complexity, sd=4))
    if(n_vert >= 3){
      verts <- st_sample(grid, size=n_vert)
      verts <- st_as_sf(verts)
      poly <- concaveman(verts) # generate the concave polygon for these vertices
      
      # randomly scale and translate the polygons to have different size ranges.
      scale_flag <- 1-(0.1*rpois(1, 1))
      # print(paste("Scaling by:", scale_flag))
      
      shift_flag_x <- 0.1*runif(1, -2, 2) # translate by x factor
      shift_flag_y <- 0.1*runif(1, -2, 2) # translate by y factor
      # print(paste("Translating by:", shift_flag_x, shift_flag_y))
      
      poly_geom <- st_geometry(poly)
      poly_shift <- (poly_geom*scale_flag+matrix(data=c(shift_flag_x, shift_flag_y), ncol=2))
      
      # check to make sure that the polygon is intersecting something on the grid
      if(length(st_intersects(grid,poly_shift)) <= 1){
        poly_shift <- concaveman(verts)
      }
      
      poly <- poly_shift %>% st_crop(grid)
      poly_ls[[i]] <- poly
      i <- i+1
    }
    else{
      n_vert <- round(rnorm(1, mean=complexity, sd=4))
    }
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

## expit and logit functions
expit <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x/(1-x))

make.ranges <- function(nsp, nsite, type.range) {

  if(type.range == 'all'){
    nsite.by.sp <- rep(nsite, nsp)
  }
   if(type.range == 'equal') {
     nsite.by.sp <- sample.int(n=nsite, size=nsp, replace=TRUE)
  } else if(type.range == 'logn') {
    prop.sites <- rbeta(nsp, shape1=1, shape2=3)
    nsite.by.sp <- round(nsite * prop.sites)
  } else if(type.range == 'polys'){
    polys <- GenPoly(numpoly=nsp, gridsize=sqrt(nsite))
    res <- array(unlist(PolyToMatrix(polys)), dim=c(nsp, nsite))
    return(res)
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


visit.history <- function(nsp, nsite, nyr, nvisit, type.visit, mu.v.0, mu.v.yr, prop.visits.same){
  
  ## for later: add a spatial component that maybe visits cluster in space through time?
  
  if(type.visit == 'visit_all'){
    
    ## all of the visits for all of the species
    
    vis.arr <- array(1,
                     dim=c(nsp=nsp, 
                           nsite=nsite,
                           nyr=nyr,
                           nvisit=nvisit))
    visit_collector <- array(1,
                             dim=c(nsite=nsite,
                                   nyr=nyr,
                                   nvisit=nvisit))
    
  }else if(type.visit == 'visit_miss'){
    
    ## there are missing visits
    
    ## Here some of the visits are the same for all species or some of the visits are different for species
    ## Depends on the proportion of visits that are the same
    ## If prop.visits.same =0 then all of the species have separate visits and if prop.visits.same=1 all of the species have the same visits 
    
    ## visit_collector is like the "Collector ID" where we identify which visits were community visits = 1 and which visits were independent sampling visits = 0
    
    ## provide parameters for the change of visits through time
    
    sigma.v.0 <- 1
    sigma.v.yr <- 1
    
    mu.v.0.sp <- rnorm(nsp, mu.v.0, sd = sigma.v.0)
    mu.v.yr.sp <- rnorm(nsp, mu.v.yr, sd = sigma.v.yr)
    
    
    vis.arr <- array(NA,
                     dim=c(nsp=nsp, 
                           nsite=nsite,
                           nyr=nyr,
                           nvisit=nvisit))
    
    visit_collector <- array(NA,
                             dim=c(nsite=nsite,
                                   nyr=nyr,
                                   nvisit=nvisit))
    
      for(site in 1:nsite) {
        for(yr in 1:nyr) {
          for(visit in 1:nvisit) {
            
            # for each visit sample whether its a community or individual visit
            
            same_vis <- rbinom(1, 1, prob = prop.visits.same)
            
            visit_collector[site,yr,visit] <- same_vis
            
            if(same_vis == 1){
              
              # if it is a community visit, then all of the species get the same
              
              vis.arr[,site,yr,visit] <- rbinom(1,1,expit(mu.v.0 + mu.v.yr*(yr-1)))
              
            }else(
              for(sp in 1:nsp){
                
                # If it is an individual species visit then not all of the species get the same
                
                vis.arr[sp,site,yr,visit] <- rbinom(1,1, expit(mu.v.0.sp[sp] + mu.v.yr.sp[sp]*(yr-1)))
                  
              }
            )
          }
        }
      }
  }
  
  return(list(vis.arr = vis.arr, visit_collector = visit_collector))
}


make.data <- function(## data structure set up 
                      nsp=20,
                      nsite=18, ## number of sites
                      nyr=10, ## number of years
                      nvisit=6, ## number of samples per year
                      ## detection
                      mu.p.0 = 0,
                      p.yr = 0.5,
                      sigma.p.site = 0.3,
                      sigma.p.sp = 0.3,
                      ## persistence
                      mu.phi.0 = 0,
                      sigma.phi.sp = 0.5,
                      mu.phi.yr = 0.5,
                      sigma.phi.yr = 0.2,
                      ## colonization
                      mu.gam.0 = 0,
                      sigma.gam.sp = 0.5,
                      mu.gam.yr = 0.5, 
                      sigma.gam.yr = 0.2,
                      ## visit
                      mu.v.0 = 0, 
                      mu.v.yr = -0.5,
                      ## type sym
                      type.range = "equal",
                      type.visit = 'visit_miss',
                      prop.visits.same = 1){

  
  ## ------------------------------------------------------------
  ## Simulate species ranges them for the
  ## specified number of species and sites.
  
  sp.range <- make.ranges(nsp, nsite, type.range)
  
  ## ------------------------------------------------------------
  ## Get visit history
  
  vis.info <- visit.history(nsp=nsp, nsite=nsite, 
                           nyr=nyr, nvisit=nvisit, type.visit=type.visit, mu.v.0, mu.v.yr,
                           prop.visits.same=prop.visits.same)
  vis.arr <- vis.info$vis.arr
  
  ## --------------------------------------------------
  ## specify species-specific colonization and persistence probabilities

  ## species-specific random intercepts
  phi.sp <- rnorm(n=nsp, mean=0, sd=sigma.phi.sp)
  gam.sp <- rnorm(n=nsp, mean=0, sd=sigma.gam.sp)

  ## effect of year on colonization and persistence (species-specific random slopes)
  phi.yr <- rnorm(n=nsp, mean=mu.phi.yr, sd=sigma.phi.yr)
  gam.yr <- rnorm(n=nsp, mean=mu.gam.yr, sd=sigma.gam.yr)
  
  ### detection
  
  ## effect of site on detection (year-specific)
  p.site <- matrix(rnorm(n=nsite*nyr, mean=0, sd=sigma.p.site),
                   nrow=nsite,
                   ncol=nyr)
  p.sp   <- rnorm(n=nsp, mean=0, sd=sigma.p.sp)
  
  
  
  
  ## Create matrices for gam and phi
  
  phi.mat <- array(NA, dim =c(nsp, nsite, nyr))
  gam.mat <- array(NA, dim =c(nsp, nsite, nyr))

  
  p.mat <- array(NA, dim =c(nsp, nsite, nyr, nvisit))
  
  for(site in 1:nsite) {
    for(yr in 1:nyr) {
      for(sp in 1:nsp) {
        phi.mat[sp,site,yr] <- expit(mu.phi.0 +
                                       phi.sp[sp] +
                                       phi.yr[sp]*(yr-1))
        
        gam.mat[sp,site,yr] <- expit(mu.gam.0 +
                                       gam.sp[sp] +
                                       gam.yr[sp]*(yr-1))
        
        for(visit in 1:nvisit) {
          
            p.mat[sp,site,yr,visit] <- expit(mu.p.0 +
                                               p.sp[sp] +
                                               p.site[site,yr] +
                                               p.yr*(yr-1))
          }
        }
      }
  }
  
  ## create initial occupancy 
  
  psi1 <- runif(nsp)
  

  ## --------------------------------------------------
  ## Generate presence and absence
  
  ## generate initial presence absence
  
  Z <- array(NA, dim=c(nsp=nsp,
                             nsite=nsite,
                             nyr=nyr))
  
  
  ## generate year 1 occupancy only in sites within the range
  
  for(sp in 1:nsp){
    Z[sp,which(sp.range[sp,]),1] <- rbinom(n=sum(sp.range[sp,]), size=1, prob=psi1[sp])
  }
  
  ### generate true presence absence for subsequent years
  
  for(yr in 2:nyr){
    for(sp in 1:nsp){
     
      exp.Z <- Z[sp,,yr-1] * phi.mat[sp,,yr-1] + (1 - Z[sp,,yr-1]) * gam.mat[sp,,yr-1]
      Z[sp, which(sp.range[sp,]) ,yr] <- rbinom(n = sum(sp.range[sp,]), size = 1, prob = exp.Z[which(sp.range[sp,])])

    }
  }

  
  ## --------------------------------------------------
  ## Generate detection non detection data 
  
  V <- array(NA, dim=c(nsp=nsp,
                       nsite=nsite,
                       nyr=nyr,
                       nvisit = nvisit))
  
  for(sp in 1:nsp){
    for(yr in 1:nyr){
      for(v in 1:nvisit){
        V[sp, which(sp.range[sp,]),yr,v] <- rbinom(n = sum(sp.range[sp,]), size = 1, prob = Z[sp, which(sp.range[sp,]),yr] * p.mat[sp,which(sp.range[sp,]),yr, v])
      }
    }
  }
  
  ## --------------------------------------------------
  ## Generate visit data 
  
  X <-V*vis.arr
  

  ## add dimension names to arrays
  
  arr.names <- list(sp=paste('sp',str_pad(1:nsp,4,pad='0'),sep='_'),
                    site=paste('site',str_pad(1:nsite,4,pad='0'),sep='_'),
                    yr=paste('y',1:nyr,sep='_'),
                    visit=paste('v',1:nvisit,sep='_'))
  dimnames(sp.range) <- arr.names[1:2]
  dimnames(Z)        <- arr.names[1:3]
  dimnames(X)        <- arr.names[1:4]
  
  
  return(list(sp.range=sp.range,
       vis.arr=vis.arr,
       vis.col=vis.info$visit_collector,
       Z=Z,
       X=X,
       nsp=nsp,
       nsite=nsite, ## number of sites
       nyr=nyr, ## number of years
       nvisit=nvisit, ## number of samples per year
       ## detection
       mu.p.0 = mu.p.0,
       p.yr = p.yr,
       sigma.p.site = sigma.p.site,
       sigma.p.sp = sigma.p.sp,
       ## persistence
       mu.phi.0 = mu.phi.0,
       sigma.phi.sp = sigma.phi.sp,
       mu.phi.yr = mu.phi.yr,
       sigma.phi.yr = sigma.phi.yr,
       ## colonization
       mu.gam.0 = mu.gam.0,
       sigma.gam.sp = sigma.gam.sp,
       mu.gam.yr = mu.gam.yr, 
       sigma.gam.yr = sigma.gam.yr,
       ## visit
       mu.v.0 = mu.v.0, 
       mu.v.yr = mu.v.yr,
       ## type sym
       type.range = type.range,
       type.visit = type.visit,
       prop.visits.same = prop.visits.same))

}
