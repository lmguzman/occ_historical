library(nimble)
ms_nimble <- nimbleCode({
  
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
  
  #### occupancy. 
  
  ## intercepts
  
  mu.psi.0 ~ dnorm(mean = 0, sd = 10)

  # species specific intercepts
  
  sigma.psi.sp ~ dunif(0,10)

  for(sp in 1:nsp){
    psi.sp[sp] ~ dnorm(0,sd = sigma.psi.sp)
  }
  
  ## species specific slopes 
  
  mu.psi.yr ~ dnorm(0, sd = 10)

  sigma.psi.yr ~ dunif(0,10)

  for(sp in 1:nsp){
    psi.yr[sp] ~ dnorm(mu.psi.yr, sd = sigma.psi.yr)
  }
  
  
  
  for(sp in 1:nsp) {
    for(yr in 1:nyr) {
      ## occupancy
      logit(psi[yr,sp]) <-
        mu.psi.0 +
        psi.yr[sp]*(yr-1) +
        psi.sp[sp]
      ## detection
      for(site in 1:nsite) {
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
        Z[site,yr,sp] ~ dbern(psi[yr,sp])
        
      }
      
    }
    
  }
  
  for(ind in 1:nind) {
    mu.p[ind] <-
      Z[sitev[ind],yrv[ind],spv[ind]] *
      p[sitev[ind],yrv[ind],spv[ind]]
    X[ind] ~ dbern(mu.p[ind])
  }
  
})


params <- c('mu.p.0',
            'p.yr',
            'sigma.p.sp',
            'sigma.p.site',
            'mu.psi.0',
            'sigma.psi.sp',
            'mu.psi.yr',
            'sigma.psi.yr')














