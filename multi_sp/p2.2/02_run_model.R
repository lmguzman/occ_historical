#!/usr/bin/env Rscript
args <-
  commandArgs(trailingOnly=TRUE)


library(parallel)
library(stringr)
source('multi_sp/simulation/src/prep_data.R')

file_sim <- list.files("multi_sp/p2.2/outputs/sim.data/")

run_prep_model <- function(file, case, model){
  
  load(paste0("multi_sp/p2.2/outputs/sim.data/", file), verbose = TRUE)
  
  ## source model
  source(sprintf('multi_sp/simulation/models/%s.R', model))

  data.prepped <- prep.data(limit.to.visits = case, sim.data)
  
  ## set model parameters 
  
  n.iter=1e4
  n.burnin=1e2
  n.adapt=1e2
  n.thin=3e1
  n.chains=3
  
  ## Initial values 
  Zst <- array(1,dim=c(sim.data$nsite,sim.data$nyr,sim.data$nsp))
  
  inits <- list(Z = Zst,
                mu.p.0 = 0,
                p.yr = 0,
                sigma.p.sp = 0.1,
                p.sp = rep(0,sim.data$nsp),
                sigma.p.site = 0.1,
                p.site = matrix(0, nrow = sim.data$nsite, ncol = sim.data$nyr), 
                mu.psi.0 = 0,
                sigma.psi.sp = 0.1,
                psi.sp = rep(0,sim.data$nsp),
                mu.psi.yr = 0,
                sigma.psi.yr = 0.1,
                psi.yr = rep(0,sim.data$nsp))
  
  ### nimble model 
  
  R.model <- nimbleModel(ms_nimble,
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
       file=paste0("multi_sp/p2.2/outputs/model.res/", case, file))
}

run_id <-
  as.numeric(args)
print(run_id)


  run_prep_model(file_sim[run_id], case = 'all', model = 'ms-nimble')
  
  run_prep_model(file_sim[run_id], case = 'detected', model = 'ms-nimble')
  
  run_prep_model(file_sim[run_id], case = 'visits', model = 'ms-nimble')
  
  run_prep_model(file_sim[run_id], case = 'community', model = 'ms-nimble')


