library(nimble)
load("multi_sp/odes_caseStudy/output/ODE_env.RData")

n.iter=1e5
n.burnin=1e3
n.adapt=1e3
n.thin=7.5e1
n.chains=3

# DETECTED ALL MODEL
Zst <- array(1,dim=c(dd_det_all$nsite,dd_det_all$nyr,dd_det_all$nsp))
inits <- list(Z = Zst,
              mu.p.0 = 0,
              p.yr = 0,
              sigma.p.sp = 0.1,
              p.sp = rep(0,dd_det_all$nsp),
              sigma.p.site = 0.1,
              p.site = matrix(0, nrow = dd_det_all$nsite, ncol = dd_det_all$nyr), 
              mu.psi.0 = 0,
              sigma.psi.sp = 0.1,
              psi.sp = rep(0,dd_det_all$nsp),
              mu.psi.yr = 0,
              sigma.psi.yr = 0.1,
              psi.yr = rep(0,dd_det_all$nsp))

# nimble model 
R.model <- nimbleModel(ms_nimble,
                       constants= dd_det_all_prep$my.constants,
                       data= dd_det_all_prep$my.data,
                       inits=inits,
                       check=FALSE)
mcmc.spec <- configureMCMC(R.model,
                           print=TRUE,
                           monitors=params,
                           thin=n.thin)

mcmc <- buildMCMC(mcmc.spec)

# compile model in C++
C.model <- compileNimble(R.model)
C.mcmc <- compileNimble(mcmc, project=R.model)

# run model
message('running model')

res_det_all <- runMCMC(C.mcmc,
                       niter=n.iter,
                       nburnin=n.burnin,
                       nchains=n.chains)

saveRDS(res_det_all, "ODE_res_det_all.rds")