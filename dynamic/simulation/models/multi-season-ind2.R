model.jags <- function() {
  
  ### priors
  
  # species specific detection
  for(sp in 1:nsp){
    p[sp] ~ dunif(0,1)
  }
  
  # random effect of year on initial colonization and extinction
  
  sigma.phi.sp ~ dunif(0,10)
  tau.phi.sp <- 1/(sigma.phi.sp*sigma.phi.sp)
  sigma.gam.sp ~ dunif(0,10)
  tau.gam.sp <- 1/(sigma.gam.sp*sigma.gam.sp)
  
  mu.phi.0 ~ dnorm(0,0.01)
  mu.gam.0 ~ dnorm(0,0.01)
  
  for(sp in 1:nsp){
    phi.sp[sp] ~ dnorm(0,tau.phi.sp)
    gam.sp[sp] ~ dnorm(0,tau.gam.sp)
  }
  
  ## effect of temperature on colonization and extinction
  
  phi.tmax ~ dnorm(0,0.01)
  gam.tmax ~ dnorm(0,0.01)
  
  ## occupancy in year 1
  
  for(sp in 1:nsp){
    for(site in 1:nsite){
      
      mu.psi.1[site,1,sp] ~ dunif(0,1)
      psi[site,1,sp] <- mu.psi.1[site,1,sp]
    
    }
  }
  
  ## occupancy and detectability in subsequent years
  
  for(sp in 1:nsp) {
    for(site in 1:nsite) {
      for(yr in 1:(nyr-1)) {
        logit(phi[site,yr,sp]) <-
          mu.phi.0 +
          phi.sp[sp] +
          phi.tmax*tmax[site,yr]
        
        logit(gam[site,yr,sp]) <-
          mu.gam.0 +
          gam.sp[sp] +
          gam.tmax*tmax[site,yr]
        
        psi[site,yr+1,sp] <-
          (Z[site,yr,sp] * phi[site,yr,sp] +
             (1-Z[site,yr,sp]) * gam[site,yr,sp])
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
      Z[site[ind],yr[ind],sp[ind]] *
      p[sp[ind]]
    X[ind] ~ dbern(mu.p[ind])
  }
  
}
get.params <- function()
  c('mu.phi.0',
    'mu.gam.0',
    'phi.sp',
    'gam.sp',
    'phi.tmax',
    'gam.tmax',
    'p')















