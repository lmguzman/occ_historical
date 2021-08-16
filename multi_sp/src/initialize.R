library(abind)
library(parallel)
library(rjags)
library(R2jags)
library(runjags)
library(R2WinBUGS)

source('src/simulate_data.R')
source('src/prep_data.R')
source('src/run_model.R')

## expit and logit functions
expit <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x/(1-x))

## standard error function
se <- function(x) sqrt(var(x)/(length(x)-1))

## nice little pdf function
pdf.f <- function(f, file, ...) {
  cat(sprintf('Writing %s\n', file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}

## load and return loaded object
load.local <- function(file) {
 v <- load(file)
 stopifnot(length(v) == 1)
 get(v)
}

id <- function(x) unique(sort(x))

## make a colour transparent
make.transparent <- function(..., alpha=0.5) {
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  alpha = floor(255*alpha)  
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  apply(newColor, 2, .makeTransparent, alpha=alpha)
}

## make a summary for a run.jags model
make.summary <- function(jags.out) {
  ## create sims.matrix
  sims.mat <- do.call(rbind, jags.out$mcmc)
  vars <- colnames(sims.mat)
  mean <- colMeans(sims.mat)
  quantiles <- apply(sims.mat, 2, quantile, p=c(0.025,0.975))
  cbind(mean, t(quantiles))
}
