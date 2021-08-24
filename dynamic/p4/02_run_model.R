#!/usr/bin/env Rscript
args <-
  commandArgs(trailingOnly=TRUE)

library(parallel)
library(stringr)
source('dynamic/simulation/src/prep_data.R')
print("Loaded prep_data from source...")

file_sim <- list.files("p3/outputs/sim.data/")
print("Obtained simulated data list...")

run_prep_model <- function(file, case, model){
  
  load(paste0("p4/outputs/sim.data/", file), verbose = TRUE)
  
  ## source model
  source(sprintf('dynamic/simulation/models/%s.R', model))

  data.prepped <- prep.data(limit.to.visits = case, sim.data)
  
  ## set model parameters 
  
  n.iter=1e5
  n.burnin=1e3
  n.adapt=1e3
  n.thin=3e2
  n.chains=3
  
  ## Initial values 
  Zst <- array(1,dim=c(sim.data$nsite,sim.data$nyr,sim.data$nsp))
  
  inits <- list(Z = Zst,
                mu.p.0 = 0,
                p.yr = 0,
                sigma.p.sp = 0.1,
                sigma.p.site = 0.1,
                mu.phi.0 = 0,
                mu.gam.0 = 0,
                sigma.phi.sp = 0.1,
                sigma.gam.sp = 0.1,
                mu.phi.yr = 0,
                mu.gam.yr = 0,
                sigma.phi.yr = 0.1,
                sigma.phi.yr = 0.1)
  
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
       file=paste0("p4/outputs/model.res/", case, file))
}

run_id <-
  as.numeric(args)
print(run_id)

run_prep_model(file_sim[run_id], case = 'all', model = 'dynamic-nimble')

run_prep_model(file_sim[run_id], case = 'detected', model = 'dynamic-nimble')

run_prep_model(file_sim[run_id], case = 'visits', model = 'dynamic-nimble')

#mclapply(file_sim, run_prep_model, case = 'all', model = 'dynamic-nimble', mc.cores = 40)

#mclapply(file_sim, run_prep_model, case = 'detected', model = 'dynamic-nimble', mc.cores = 40)

#mclapply(file_sim, run_prep_model, case = 'visits', model = 'dynamic-nimble', mc.cores = 40)