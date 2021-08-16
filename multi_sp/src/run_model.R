run.model <- function(dd,
                      model,
                      n.iter=1e5,
                      n.burnin=1e3,
                      n.adapt=1e3,
                      n.thin=1e2) {

  model.txt <- sprintf('src/models/%s.txt', model)
  
  if(model=='ms_ind') {
    jags.data <- list(X=dd$X,
                      yr=dd$master.index[,'yr'],
                      site=dd$master.index[,'site'],
                      sp=dd$master.index[,'sp'],
                      ## visit=dd$master.index[,'visit'],
                      nsp=dd$nsp,
                      nsite=dd$nsite,
                      nyr=dd$nyr,
                      ## nvisit=dd$nvisit,
                      nind=dd$nind)
    ## Initial values 
    Zst <- array(1,dim=c(dd$nsp,dd$nsite,dd$nyr))
  }
  if(model=='ms_delta') {
    jags.data <- list(X=dd$X,
                      nsp=dd$nsp,
                      nsite=dd$nsite,
                      nyr=dd$nyr,
                      nvisit=dd$nvisit,
                      in.range=dd$in.range)
    ## Initial values 
    Zst <- array(1,dim=c(dd$nsp,dd$nsite,dd$nyr))
    in.range.expanded <- array(dd$in.range, dim=dim(Zst))
    Zst <- Zst*in.range.expanded
  }

  make.inits <- function() {
    RNG <- parallel.seeds("base::BaseRNG", 1)
    c(list(Z=Zst), RNG[[1]])
  }
  inits1 <- make.inits()
  inits2 <- make.inits()
  inits3 <- make.inits()

  jags.out <- run.jags(model=model.txt,
                       monitor=get.params(),
                       data=jags.data,
                       inits=list(inits1,inits2,inits3),
                       n.chains=3,
                       burnin=n.burnin,
                       sample=floor(n.iter/n.thin),
                       thin=n.thin,
                       adapt=n.adapt,
                       method='parallel')
  
  list(jags.data=jags.data, jags.out=jags.out)
}
