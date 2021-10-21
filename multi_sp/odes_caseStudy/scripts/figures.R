# Odonate figures

# load libraries
library(tidyverse); library(cowplot); library(ggpubr);

# load in .RDS files from Compute Canada runs
ode_det_all <- readRDS("../output/ODE_res_det_all.rds")
ode_det_range <- readRDS("../output/ODE_res_det_range.rds")

samplech1 <- as.mcmc(ode_det_all$chain1)
samplech2 <- as.mcmc(ode_det_all$chain2)
samplech3 <- as.mcmc(ode_det_all$chain3)

coda.samples <- as.mcmc.list(samplech1, samplech2, samplech3)

summary_mcmc <- coef(coda.samples, simplify = TRUE)
rhat_mcmc <- map_df(rhat(coda.samples, by = "term"), ~data.frame(rhat = .x), .id = 'term')