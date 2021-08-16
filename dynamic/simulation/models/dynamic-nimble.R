library(nimble)
multi_season_nimble <- nimbleCode({
  
  ### priors
  
  #### detection #####
  
  mu.p.0   ~ dnorm(mean = 0, sd= 10)
  
  ## effect of yr on detection
  p.yr ~ dnorm(mean = 0, sd = 10)
  
  ## random effect of species on detection
  sigma.p.sp   ~ dunif(min = 0, max = 10)
  
  for(sp in 1:nsp) {
    p.sp[sp]   ~ dnorm(mean = 0, sd = sigma.p.sp)
  }
  # 
  ## random effect of site on detection
  sigma.p.site ~ dunif(min = 0, max = 10)
  
  for(site in 1:nsite) {
    for(yr in 1:nyr) {
      p.site[site,yr]   ~ dnorm(0, sd = sigma.p.site)
    }
  }
  
  #### colonization and extinction. 
  
  ## intercepts
  
  mu.phi.0 ~ dnorm(mean = 0, sd = 10)
  mu.gam.0 ~ dnorm(mean = 0, sd = 10)
  
  # species specific intercepts
  
  sigma.phi.sp ~ dunif(0,10)
  sigma.gam.sp ~ dunif(0,10)

  for(sp in 1:nsp){
    phi.sp[sp] ~ dnorm(0,sd = sigma.phi.sp)
    gam.sp[sp] ~ dnorm(0,sd = sigma.gam.sp)
  }
  
  ## species specific slopes 
  
  mu.phi.yr ~ dnorm(0, sd = 10)
  mu.gam.yr ~ dnorm(0, sd = 10)
  
  sigma.phi.yr ~ dunif(0,10)
  sigma.gam.yr ~ dunif(0,10)
  
  for(sp in 1:nsp){
    phi.yr[sp] ~ dnorm(mu.phi.yr, sd = sigma.phi.yr)
    gam.yr[sp] ~ dnorm(mu.gam.yr, sd = sigma.gam.yr)
  }
  
  
  ### True occupancy 
  
  ## occupancy in year 1
  
  for(sp in 1:nsp){
    
    mu.psi.1[sp] ~ dunif(0,1)
    
    for(site in 1:nsite){
      
      psi[site,1,sp] <- mu.psi.1[sp]
      
    }
  }
  
  ## occupancy in subsequent years
  
  for(sp in 1:nsp) {
    for(site in 1:nsite) {
      for(yr in 1:(nyr-1)) {
        logit(phi[site,yr,sp]) <-
          mu.phi.0 +
          phi.sp[sp] +
          phi.yr[sp]*yr 
        
        logit(gam[site,yr,sp]) <-
          mu.gam.0 +
          gam.sp[sp] +
          gam.yr[sp]*yr 
        
        psi[site,yr+1,sp] <-
          (Z[site,yr,sp] * phi[site,yr,sp] +
             (1-Z[site,yr,sp]) * gam[site,yr,sp])
      }
    }
  }
  
  ### detection 
  
  for(sp in 1:nsp){
    for(site in 1:nsite){
      for(yr in 1:nyr){
        logit(p[site,yr,sp]) <-
          mu.p.0 +
          p.yr*(yr-1) +
          p.sp[sp] +
          p.site[site,yr]
        
      }
    }  
  }
  
  ## latent state and likelihood    
  
  for(sp in 1:nsp){
    for(site in 1:nsite){
      for(yr in 1:nyr) {
        
        #latent state
        Z[site,yr,sp] ~ dbern(psi[site,yr,sp])
        
      }
      
    }
    
  }
  
  for(ind in 1:nind) {
    mu.p[ind] <-
      Z[sitev[ind],yrv[ind],spv[ind]] *
      p[sitev[ind],yrv[ind],spv[ind]]
    X[ind] ~ dbern(mu.p[ind])
  }
  
  ### derived quantities ##
  
  ## species occupancy per year
  for(sp in 1:nsp){
    for(yr in 1:nyr){
      psi.sp.yr[yr,sp] <- mean(psi[1:nsite,yr,sp])
    }
  }
  
})


params <- c('mu.p.0',
            'p.yr',
            'sigma.p.sp',
            'sigma.p.site',
            'mu.phi.0',
            'mu.gam.0',
            'sigma.phi.sp',
            'sigma.gam.sp',
            'mu.phi.yr',
            'mu.gam.yr',
            'sigma.phi.yr',
            'sigma.phi.yr')














