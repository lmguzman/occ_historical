#p2
#!/usr/bin/env Rscript
args <-
  commandArgs(trailingOnly=TRUE)


library(parallel)
library(stringr)
library(dplyr)
source('~/scratch/occ_historical/multi_sp/simulation/src/prep_data3.R')

file_sim <- list.files("~/scratch/occ_historical/multi_sp/p2.2/outputs/sim.data/")

run_prep_model <- function(file, case, range, model, time.interval.yr, time.interval.visit=3){
  
  load(paste0("~/scratch/occ_historical/multi_sp/p2.2/outputs/sim.data/", file), verbose = TRUE)

    ## source model
  source(sprintf('~/scratch/occ_historical/multi_sp/simulation/models/%s.R', model))


  data.prepped <- prep.data(dd = sim.data, limit.to.visits = case, limit.to.range=range, time.interval.yr=time.interval.yr, time.interval.visit=time.interval.visit)
  
  ## set model parameters 
  
  n.iter=1e5
  n.burnin=1e3
  n.adapt=1e2
  n.thin=1e2
  n.chains=3

  ## Initial values 
  Zst <- array(1,dim=c(data.prepped$my.constants$nsite,data.prepped$my.constants$nyr,data.prepped$my.constants$nsp))
  
  inits <- list(Z = Zst,
                mu.p.0 = 0,
                p.yr = 0,
                sigma.p.sp = 0.1,
                p.sp = rep(0,data.prepped$my.constants$nsp),
                sigma.p.site = 0.1,
                p.site = matrix(0, nrow = data.prepped$my.constants$nsite, ncol = data.prepped$my.constants$nyr), 
                mu.psi.0 = 0,
                sigma.psi.sp = 0.1,
                psi.sp = rep(0,data.prepped$my.constants$nsp),
                mu.psi.yr = 0,
                sigma.psi.yr = 0.1,
                psi.yr = rep(0,data.prepped$my.constants$nsp),
                mu.p = rep(0, data.prepped$my.constants$nind))
  
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
       file=paste0("~/scratch/occ_historical/multi_sp/p2.2/outputs/model.res/", case, range, "_eras", time.interval.yr, file))
}

run_id <-
  as.numeric(args)
print(run_id)

  run_prep_model(file_sim[run_id], case = 'all', range="no", model = 'ms-nimble', time.interval.yr = 2)
  
  run_prep_model(file_sim[run_id], case = 'detected', range="no", model = 'ms-nimble', time.interval.yr = 2)
  
  run_prep_model(file_sim[run_id], case = 'visits', range="no", model = 'ms-nimble', time.interval.yr = 2)
  


