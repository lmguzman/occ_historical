setwd('~/Dropbox/occ_historical/multi_sp')
rm(list=ls())

source('src/initialize.R')

set.seed(1)

## VISITS: select only one of the below
## case.visits <- 'visits'
## case.visits <- 'detected'
case.visits <- 'all'

## SITES: select only one of the below
case.sites <- 'range'
## case.sites <- 'all'

## MODEL: select only one of the below
method <- 'ind' ## index method (only index over relevant values)
## method <- 'delta' ## delta method (multiple by 0/1 for sites in/out of range)

## ------------------------------------------------------------
## create simulated data (if using simulated data)
sim.data <- make.data(nsp=25,
                      nsite=250,
                      nyr=5,
                      nvisit=3,
                      mu.psi=0,
                      sigma.psi.sp=0.8,
                      mu.psi.yr=0,
                      sigma.psi.yr=0.7,
                      mu.p=-0.5,
                      sigma.p.sp=0.6,
                      p.yr=0,
                      sigma.p.site=1.5,
                      mu.v=-0.5,
                      mu.v.yr=-0.5,
                      missing.visits=TRUE)
## ------------------------------------------------------------

## ------------------------------------------------------------
## model

my.data <- prep.data(dd=sim.data,
                     limit.to.visits=case.visits,
                     limit.to.sites=case.sites,
                     method=method)

## source JAGS model
model <- sprintf('ms_%s', method)
source(sprintf('src/models/%s.R', model))
write.model(model.jags, con=sprintf('src/models/%s.txt', model))

scale <- 1e2
source('src/run_model.R')
model.out <- run.model(dd=my.data,
                       model=model,
                       n.burnin=1e3*scale,
                       n.adapt=1e3*scale,
                       n.iter=1e3*scale,
                       n.thin=1e1*scale)

fn <- sprintf('%s_%s', model, case.visits)
save(model.out,
     my.data,
     sim.data,
     file=sprintf('saved/%s.RData', fn))
## ------------------------------------------------------------


## --------------------------------------------------
## load data and create summary (skip if summary exists)
load(sprintf('saved/%s-%s.RData', model, case.visits), verbose=TRUE)

## ## create my own summary without gelman statistics, etc
## summ <- make.summary(model.out$jags.out)

## create full summary
model.summary <-  add.summary(model.out$jags.out)
get.summ <- function(pars) {
  summ <- round(cbind(
    model.summary$summary$statistics[pars,'Mean',drop=FALSE],
    model.summary$summary$quantiles[pars,c('2.5%', '97.5%'),drop=FALSE],
    Rhat=model.summary$psrf$psrf[pars,1]
  ), digits=3)
  colnames(summ)[1] <- 'mean'
  summ
}
vars <- rownames(model.summary$psrf$psrf)
summ <- get.summ(vars)

save(model.out,
     my.data,
     sim.data,
     model.summary,
     file=sprintf('saved/%s_with_summ.RData', fn))
## --------------------------------------------------

## ## --------------------------------------------------
## load(sprintf('%s/%s_with_summ.RData', dir, fn), verbose=TRUE)

## sims.arr <- aperm(sapply(res$mcmc, I, simplify='array'), c(1,3,2))
## vars <- dimnames(sims.arr)[[3]]
## ## --------------------------------------------------

## --------------------------------------------------
## run-time info
model.out$jags.out$burnin
model.out$jags.out$sample
model.out$jags.out$thin
sprintf('%2.1f minutes', as.numeric(model.out$jags.out$timetaken/60))
sprintf('%2.1f hour', as.numeric(model.out$jags.out$timetaken/60/60))
## --------------------------------------------------

## --------------------------------------------------
## examine summary
summ.paper <- summ[c('mu.psi',
                     'sigma.psi.sp',
                     'mu.psi.yr',
                     'sigma.psi.yr',
                     'mu.p',
                     'sigma.p.sp',
                     'p.yr',
                     'sigma.p.site'),]
summ.paper
