# Odonate figures

# Load libraries
library(tidyverse); library(cowplot); library(ggpubr);
library(stringr); library(mcmcr); library(coda);
library(purrr); library(data.table); library(tidybayes);
library(modelr)

# Load in .RDS files from Compute Canada runs
ode_det_all <- readRDS("../output/ODE_res_det_all.rds")
ode_det_range <- readRDS("../output/ODE_res_det_range.rds")

# Detected-All Model
samplech1 <- as.mcmc(ode_det_all$chain1)
samplech2 <- as.mcmc(ode_det_all$chain2)
samplech3 <- as.mcmc(ode_det_all$chain3)

coda.samples <- as.mcmc.list(samplech1, samplech2, samplech3)

summary_mcmc <- coef(coda.samples, simplify = TRUE) %>%
  dplyr::mutate(param=term)
rhat_mcmc <- map_df(rhat(coda.samples, by = "term"), ~data.frame(rhat = .x), .id = 'term')
which(rhat_mcmc[,2:ncol(rhat_mcmc)]>=1.1)

# Detected-Range Model
samplech1 <- as.mcmc(ode_det_range$chain1)
samplech2 <- as.mcmc(ode_det_range$chain2)
samplech3 <- as.mcmc(ode_det_range$chain3)

coda.samples <- as.mcmc.list(samplech1, samplech2, samplech3)
sims.mat <- as.matrix(coda.samples)

summary_mcmc <- coef(coda.samples, simplify = TRUE) %>%
  dplyr::mutate(param=term)
rhat_mcmc <- map_df(rhat(coda.samples, by = "term"), ~data.frame(rhat = .x), .id = 'term')
which(rhat_mcmc[,2:ncol(rhat_mcmc)]>=1.1)

plot(NA,
     xlim=c(1,5),
     ylim=c(0,1),
     xlab='Era',
     ylab='Occupancy Probability',
     xaxt='n',
     yaxt='n',
     bty="n",
     las=1)

ss <- 1
get.y.val <- function(dd){
  chains <- expit(sims.mat[,"mu.psi.0"]+
                    sims.mat[,sprintf("psi.yr[%d]", ss)]*dd+
                    sims.mat[,sprintf("psi.sp[%d]", ss)])
  
  c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
}

x.vals <- seq(from=1, to=5, length=5)

y.vals <- lapply(x.vals, get.y.val)
y.vals2 <- do.call(rbind, y.vals)
lines(y.vals2[,"mean"], type="l", ylim=c(0,1), col="grey") # plot a single species line

# plot a community of lines
for(ss in 1:195){
  get.y.val <- function(dd){
    chains <- expit(sims.mat[,"psi.0"]+
                      sims.mat[,sprintf("psi.yr[%d]", ss)]*dd+
                      sims.mat[,sprintf("psi.sp[%d]", ss)])
  }
  x.vals <- seq(from=1, to=5, length=5)
  
  y.vals <- lapply(x.vals, get.y.val)
  y.vals2 <- do.call(rbind, y.vals)
  
  lines(y.vals2[,"mean"], type="l", ylim=c(0,1), col="grey")
}

# plot the mean community effect
get.y.val <- function(xx){
  chains <- expit(sims.mat[,"mu.psi.0"]+
                    sims.mat[,sprintf("psi.yr[%d]", ss)]*dd+
                    sims.mat[,sprintf("psi.sp[%d]", ss)])
  
  c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
}
x.vals <- seq(from=1, to=5, length=5)
y.vals <- lapply(x.vals, get.y.val)
y.vals2 <- do.call(rbind, y.vals)

lines(y.vals2[,"mean"], type="l", ylim=c(0,1), col="black", lwd=3)
lines(y.vals2[,"2.5%"], lty=2, lwd=3, col='black')
lines(y.vals2[,"97.5%"], lty=2, lwd=3, col='black')
