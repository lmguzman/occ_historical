
prep.data <- function(dd, limit.to.visits, limit.to.range, time.interval.yr, time.interval.visit) {
  
  ## keep only detected species:
  sp.keep <- which(apply(dd$X, 1, sum, na.rm=TRUE)>0)
  
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
  
  ## subset based on the above
  dd$nsite.old  <- dim(dd$X)[2]
  dd$Z        <- dd$Z[sp.keep,site.keep,,drop=FALSE]
  dd$X        <- dd$X[sp.keep,site.keep,,,drop=FALSE]
  dd$sp.range <- dd$sp.range[sp.keep,site.keep,drop=FALSE]
  dd$vis.arr  <- dd$vis.arr[sp.keep,site.keep,,,drop=FALSE]
  dd$nsp      <- length(sp.keep)
  dd$nsite    <- length(site.keep)
  dd$vis.col <- dd$vis.col[site.keep,,,drop=FALSE]
  
  
  ####### time interval ###
  
  ## convert data into long format
  
  presence.only <- which(dd$X==1, arr.ind=TRUE) %>% data.frame()
  
  visits.only <- which(dd$vis.arr==1, arr.ind=TRUE) %>% data.frame()
  
  names(presence.only) <- c("spn", "siten", "yr", "visit")
  
  names(visits.only) <- c("spn", "siten", "yr", "visit")
  
  new_year <- expand.grid(yr = 1:dd$nyr, visit = 1:dd$nvisit) %>% 
    arrange(yr) %>% 
    mutate(syr = 1:n())
  
  presence.new.year <- presence.only %>% 
    left_join(new_year)
  
  visits.new.year <- visits.only %>% 
    left_join(new_year)
  
  #### use given time interval to bin ###
  
  # time.interval.yr <- 2
  # time.interval.visit <- 3
  
  unique.syr <- unique(sort(presence.new.year$syr))
  
  time.cut.yr <- as.numeric(cut(unique.syr, time.interval.yr))
  
  time.cut.visit <- as.numeric(cut(unique.syr, time.interval.visit*time.interval.yr))
  
  time.cut.visit.2 <-  rep(rep(1:time.interval.visit, each = max(unique.syr)/(time.interval.yr*time.interval.visit)), 
                           time.interval.yr)
  
  new.time.interval <-data.frame(syr = unique.syr, time.cut.yr, time.cut.visit, time.cut.visit.2)
  
  sp_names <- data.frame(spn = 1:length(dimnames(dd$X)$sp), sp = dimnames(dd$X)$sp)
  site_names <-data.frame(siten = 1:length(dimnames(dd$X)$site), site = dimnames(dd$X)$site)
  
  presence.new.time.interval <- presence.new.year %>% 
    left_join(new.time.interval) %>% 
    group_by(spn, siten, time.cut.yr, time.cut.visit.2) %>% 
    summarise(presence = 1) %>% 
    left_join(sp_names) %>% left_join(site_names) %>% 
    mutate(yr = paste('yr',time.cut.yr,sep='_'),
           visit = paste('v',time.cut.visit.2,sep='_'))
  
  visits.new.time.interval <- visits.new.year %>% 
    left_join(new.time.interval) %>% 
    group_by(spn, siten, time.cut.yr, time.cut.visit.2) %>% 
    summarise(presence = 1) %>% 
    left_join(sp_names) %>% left_join(site_names) %>% 
    mutate(yr = paste('yr',time.cut.yr,sep='_'),
           visit = paste('v',time.cut.visit.2,sep='_'))
  
  
  # presence.new.time.interval$time.cut.yr %>% table()
  # presence.new.time.interval$time.cut.visit %>% table()
  
  occ.arr <- array(0, dim = c(dd$nsp, dd$nsite, time.interval.yr, time.interval.visit), 
                   dimnames = list(sp=dimnames(dd$X)$sp,
                                   site=dimnames(dd$X)$site,
                                   year= paste0("yr_", 1:time.interval.yr),
                                   visit=paste0("v_", 1:time.interval.visit)))
  
  occ.arr[cbind(match(presence.new.time.interval$sp, dimnames(dd$X)$sp), match(presence.new.time.interval$site, dimnames(dd$X)$site), 
                match(presence.new.time.interval$yr, paste0("yr_", 1:time.interval.yr)), 
                match(presence.new.time.interval$visit,paste0("v_", 1:time.interval.visit)))] <- 1 
  
  
  vis.arr <- array(0, dim = c(dd$nsp, dd$nsite, time.interval.yr, time.interval.visit), 
                   dimnames = list(sp=dimnames(dd$X)$sp,
                                   site=dimnames(dd$X)$site,
                                   year= paste0("yr_", 1:time.interval.yr),
                                   visit=paste0("v_", 1:time.interval.visit)))
  
  names(dim(occ.arr)) <- c("nsp", "nsite", 'nyr', 'nvisit')
  
  vis.arr[cbind(match(visits.new.time.interval$sp, dimnames(dd$X)$sp), match(visits.new.time.interval$site, dimnames(dd$X)$site), 
                match(visits.new.time.interval$yr, paste0("yr_", 1:time.interval.yr)), 
                match(visits.new.time.interval$visit,paste0("v_", 1:time.interval.visit)))] <- 1 
  
  dd$X2 <- occ.arr
  dd$vis.arr2 <- vis.arr
  
  
  ## generate master index (to improve model efficiency (this prevents
  ## unnecessary iterating through all irrelevant sites and visits)
  get.indices <- function(sp) {
    ## visited array
    vis.arr <- dd$vis.arr2[sp,,,]
    ## if modelling all visits, set visit array to 1 everywhere
    if(limit.to.visits=='all') 
      vis.arr[TRUE] <- 1
    ## if modelling sites with detections only, create new visit array
    if(limit.to.visits=='detected') {
      nsp.detected <- apply(dd$X2, 2:4, sum, na.rm = TRUE)
      vis.arr[TRUE] <- 1
      vis.arr[nsp.detected==0] <- 0
    }
    
    # either restrict or not restrict inference to inferred species ranges
    if(limit.to.range=='yes'){ # restrict to within species' ranges
      # infer species ranges from affirmative detections
      # make a placeholder polygon to create
      # a grid.
      sfc <- st_sfc(st_polygon(list(rbind(c(0,0),
                                          c(1,0),
                                          c(1,1),
                                          c(0,0)))))
      
      # create the grid using the specified number of simulated sites
      grid <- st_make_grid(sfc, cellsize=1/sqrt(dd$nsite.old), square=TRUE) %>%
        st_as_sf() %>%
        dplyr::mutate(GID=row_number())
      
      # iteratively add points for each species based on the centroids of
      # each grid cell in the matrix, then construct a convex polygon for each
      # species and add to a list
      
      
      sp_GID <- rowSums(dd$X2[sp,,,], c(2,3,4))
      sp_centroids <- c()
      for(j in 1:dd$nsite){
        if(sp_GID[j] >= 1){
          sp_centroids <- sp_centroids %>% base::append(st_centroid(grid[j,])$x)
          
        }
      }
      
      sp_poly <- st_convex_hull(st_union(sp_centroids))
      grid_size <- sqrt(dd$nsite.old)
      mtrx <- st_intersects(st_as_sf(grid),
                            st_as_sf(sp_poly))
      matrix_ls <- matrix(as.matrix(mtrx), nrow=grid_size, ncol=grid_size)
      
      vector_range <- as.vector(matrix_ls)
      vis.arr[!vector_range[site.keep],,] <- 0
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
  
  my.data <- list(X=dd$X2[master.index])
  
  my.constants <- list(
    nsp=dim(dd$X2)['nsp'],
    nsite=dim(dd$X2)['nsite'],
    nyr=dim(dd$X2)['nyr'],
    nind=nrow(master.index),
    yrv=master.index[,'yr'],
    sitev=master.index[,'site'],
    spv=master.index[,'sp'])
  
  my.info <- list(
    time.interval.yr=time.interval.yr,
    time.interval.visit=time.interval.visit,
    limit.to.visits=limit.to.visits,
    limit.to.range=limit.to.range,
    sp.keep = names(sp.keep), 
    site.keep = names(site.keep)
  )
  
  return(list(my.constants = my.constants, my.data = my.data, my.info = my.info))
}