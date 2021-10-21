# run the detected_range model
Zst <- array(1,dim=c(dd_det_range$nsite,dd_det_range$nyr,dd_det_range$nsp))
inits <- list(Z = Zst,
              mu.p.0 = 0,
              p.yr = 0,
              sigma.p.sp = 0.1,
              p.sp = rep(0,dd_det_range$nsp),
              sigma.p.site = 0.1,
              p.site = matrix(0, nrow = dd_det_range$nsite, ncol = dd_det_range$nyr), 
              mu.psi.0 = 0,
              sigma.psi.sp = 0.1,
              psi.sp = rep(0,dd_det_range$nsp),
              mu.psi.yr = 0,
              sigma.psi.yr = 0.1,
              psi.yr = rep(0,dd_det_range$nsp))

# nimble model 
R.model <- nimbleModel(ms_nimble,
                       constants= dd_det_range_prep$my.constants,
                       data= dd_det_range_prep$my.data,
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

res_det_range <- runMCMC(C.mcmc,
                       niter=n.iter,
                       nburnin=n.burnin,
                       nchains=n.chains)

saveRDS(res_det_range, "ODE_res_det_range.rds")
