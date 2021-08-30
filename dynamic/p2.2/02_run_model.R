#!/usr/bin/env Rscript
args <-
  commandArgs(trailingOnly=TRUE)


library(parallel)
library(stringr)
source('~/scratch/occ_historical/dynamic/simulation/src/prep_data.R')

file_sim <- list.files("~/scratch/occ_historical/dynamic/p2.2/outputs/sim.data/")

run_prep_model <- function(file, case, model){
  
  load(paste0("~/scratch/occ_historical/dynamic/p2.2/outputs/sim.data/", file), verbose = TRUE)
  
  ## source model
  source(sprintf('~/scratch/occ_historical/dynamic/simulation/models/%s.R', model))

  data.prepped <- prep.data(limit.to.visits = case, sim.data)
  
  ## set model parameters 
  
  n.iter=1e4
  n.burnin=1e2
  n.adapt=1e2
  n.thin=3e1
  n.chains=3
  
  ## Initial values 
  #Zst <- array(1,dim=c(sim.data$nsite,sim.data$nyr,sim.data$nsp))

  Zint <- aperm(sim.data$X[,,,1], c(2,3,1))
  Zint[is.na(Zint)] <- 0
  
  inits <- list(Z = Zint,
                mu.p.0 = 0,
                p.yr = 0,
                sigma.p.sp = 0.1,
                p.sp = rep(0, sim.data$nsp),
                p.site = matrix(0, nrow = sim.data$nsite, ncol = sim.data$nyr),  
                sigma.p.site = 0.1,
                mu.phi.0 = 0,
                mu.gam.0 = 0,
                sigma.phi.sp = 0.1,
                sigma.gam.sp = 0.1,
                phi.sp = rep(0, sim.data$nsp),
                gam.sp = rep(0, sim.data$nsp),
                mu.phi.yr = 0,
                mu.gam.yr = 0,
                sigma.phi.yr = 0.1,
                sigma.phi.yr = 0.1,
                phi.yr = rep(0, sim.data$nsp),
                gam.yr = rep(0, sim.data$nsp),
                mu.psi.1 = rep(0, sim.data$nsp))
  
  ### nimble model 
  
  R.model <- nimbleModel(multi_season_nimble,
                         constants= data.prepped$my.constants,
                         data= data.prepped$my.data,
                         inits=inits,
                         check=FALSE)
  
  mcmc.spec <- configureMCMC(R.model,
                             print=TRUE,
                             monitors=params,
                             thin=n.thin)
  
  mcmc <- buildMCMC(mcmc.spec)
  
  ## compile model in C++
  
  C.model <- compileNimble(R.model)
  
  C.mcmc <- compileNimble(mcmc, project=R.model)
  
  ## run model
  
  message('running model')
  
  res <- runMCMC(C.mcmc,
                 niter=n.iter,
                 nburnin=n.burnin,
                 nchains=n.chains)
  
  
  
  #### fix file name saving 
  save(res, data.prepped, sim.data,
       file=paste0("~/scratch/occ_historical/dynamic/p2.2/outputs/model.res/", case, file))
}

run_id <-
  as.numeric(args)
print(run_id)


  run_prep_model(file_sim[run_id], case = 'all', model = 'dynamic-nimble')
  
  run_prep_model(file_sim[run_id], case = 'detected', model = 'dynamic-nimble')
  
  run_prep_model(file_sim[run_id], case = 'visits', model = 'dynamic-nimble')
  
  run_prep_model(file_sim[run_id], case = 'community', model = 'dynamic-nimble')


