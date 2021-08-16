model.jags <- function() {

  ## mean occupancy and detection
  mu.psi ~ dnorm(0,0.01)
  mu.p   ~ dnorm(0,0.01)
  
  ## random effect of species on occupancy and detection
  sigma.psi.sp ~ dunif(0,10)
  sigma.p.sp   ~ dunif(0,10)
  tau.psi.sp  <- 1/(sigma.psi.sp*sigma.psi.sp)
  tau.p.sp    <- 1/(sigma.p.sp*sigma.p.sp)
  for(sp in 1:nsp) {
    psi.sp[sp] ~ dnorm(0, tau.psi.sp)
    p.sp[sp]   ~ dnorm(0, tau.p.sp)
  }

  ## random effect of site on detection
  sigma.p.site   ~ dunif(0,10)
  tau.p.site    <- 1/(sigma.p.site*sigma.p.site)
  for(site in 1:nsite) {
    for(yr in 1:nyr) {
      p.site[site,yr]   ~ dnorm(0, tau.p.site)
    }
  }

  ## effect of year on occupancy and detection (random slopes)
  mu.psi.yr    ~ dnorm(0,0.01)
  sigma.psi.yr ~ dunif(0,10)
  tau.psi.yr  <- 1/(sigma.psi.yr*sigma.psi.yr)
  for(sp in 1:nsp) {
    psi.yr[sp] ~ dnorm(mu.psi.yr, tau.psi.yr)
  }

  ## effect of yr on detection
  p.yr ~ dnorm(0,0.01)

  for(sp in 1:nsp) {
    for(site in 1:nsite) {
      for(yr in 1:nyr) {
        ## occupancy
        logit(psi[sp,site,yr]) <-
          mu.psi +
          psi.yr[sp]*(yr-1) +
          psi.sp[sp]
        ## create modified psi that is set to zero for sites outside
        ## each species' range
        psi.range[sp,site,yr] <- psi[sp,site,yr] * in.range[sp,site]
        ## detection
        logit(p[sp,site,yr]) <-
          mu.p +
          p.yr*(yr-1) +
          p.sp[sp] +
          p.site[site,yr]
      }
    }
  }

  ## latent Z state
  for(sp in 1:nsp) {
    for(site in 1:nsite) {
      for(yr in 1:nyr) {
        Z[sp,site,yr] ~ dbern(psi.range[sp,site,yr])
      }
    }
  }
  
  ## likelihood
  for(sp in 1:nsp) {
    for(site in 1:nsite) {
      for(yr in 1:nyr) {
        p.eff[sp,site,yr] <- Z[sp,site,yr] * p[sp,site,yr]
        for(visit in 1:nvisit) {
          X[sp,site,yr,visit] ~ dbern(p.eff[sp,site,yr])
        }
      }
    }
  }
}


## specify the parameters to be monitored
get.params <- function()
  c('mu.psi',
    'mu.p',
    'psi.sp',
    'sigma.psi.sp',
    'p.sp',
    'sigma.p.sp',
    'p.site',
    'sigma.p.site',
    'psi.yr',
    'mu.psi.yr',
    'sigma.psi.yr',
    'p.yr')
