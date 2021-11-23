# Odonate figures

# Load libraries
library(tidyverse); library(cowplot); library(ggpubr);
library(stringr); library(mcmcr); library(coda);
library(purrr); library(data.table); library(tidybayes);
library(modelr); library(colorspace)

# Load in .RDS files from Compute Canada runs
ode_det_all <- readRDS("../output/ODE_res_det_all.rds")
ode_det_range <- readRDS("../output/ODE_res_det_range.rds")

# Detected-All Model (MODEL 5)
samplech1 <- as.mcmc(ode_det_all$chain1)
samplech2 <- as.mcmc(ode_det_all$chain2)
samplech3 <- as.mcmc(ode_det_all$chain3)

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
     bty="n",
     las=1)

# plot a community of lines
spp_lines <- list()
for(ss in 1:195){
  get.y.val <- function(dd){
    chains <- plogis(sims.mat[,"mu.psi.0"]+
                       sims.mat[,sprintf("psi.yr[%d]", ss)]*dd+
                       sims.mat[,sprintf("psi.sp[%d]", ss)])
    
    c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
  }
  x.vals <- seq(from=1, to=5, length=5)
  
  y.vals <- lapply(x.vals, get.y.val)
  y.vals2_spp <- do.call(rbind, y.vals)
  
  lines(y.vals2_spp[,"mean"], type="l", ylim=c(0,1), col="grey")
  
  y.vals2_spp <- y.vals2_spp %>% as.data.frame() %>%
    dplyr::mutate(spp=ss, x=row_number())
  spp_lines[[ss]] <- y.vals2_spp
}
y.vals2_spp <- do.call(rbind, spp_lines)

# plot the mean community effect
get.y.val <- function(dd){
  chains <- plogis(sims.mat[,"mu.psi.0"]+
                     sims.mat[,sprintf("psi.yr[%d]", ss)]*dd+
                     sims.mat[,sprintf("psi.sp[%d]", ss)])
  
  c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
}
x.vals <- seq(from=1, to=5, length=5)
y.vals <- lapply(x.vals, get.y.val)
y.vals2_com <- do.call(rbind, y.vals) %>% as.data.frame()

lines(y.vals2_com[,"mean"], type="l", ylim=c(0,1), col="black", lwd=3)
lines(y.vals2_com[,"2.5%"], lty=2, lwd=3, col='black')
lines(y.vals2_com[,"97.5%"], lty=2, lwd=3, col='black')

y.vals2_spp <- y.vals2_spp %>%
  group_by(spp) %>%
  dplyr::mutate(min=dplyr::first(mean), max=dplyr::last(mean)) %>%
  dplyr::mutate(diff=max-min) %>%
  dplyr::mutate(sgn=sign(diff)) %>%
  ungroup()

## ggplot
M5_plot <- ggplot()+
  geom_line(y.vals2_spp, mapping=aes(x=x, y=mean, group=spp, color=diff), 
            alpha=0.8)+
  scale_color_continuous_divergingx(palette="RdYlBu", mid=0.0, name="Slope")+
  geom_line(y.vals2_com, mapping=aes(x=c(1:5), y=mean), size=1, linetype=1,
            color="grey10", alpha=0.8)+
  geom_line(y.vals2_com, mapping=aes(x=c(1:5), y=`2.5%`), size=1, linetype=2,
            color="grey10", alpha=0.8)+
  geom_line(y.vals2_com, mapping=aes(x=c(1:5), y=`97.5%`), size=1, linetype=2,
            color="grey10", alpha=0.8)+
  ylab("Occupancy Probability")+
  xlab("Era")+
  theme_cowplot()+
  background_grid()
M5_plot

# All-Range Model (MODEL 7)
samplech1 <- as.mcmc(ode_all_range$chain1)
samplech2 <- as.mcmc(ode_all_range$chain2)
samplech3 <- as.mcmc(ode_all_range$chain3)

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
     bty="n",
     las=1)

# plot a community of lines
spp_lines <- list()
for(ss in 1:195){
  get.y.val <- function(dd){
    chains <- plogis(sims.mat[,"mu.psi.0"]+
                       sims.mat[,sprintf("psi.yr[%d]", ss)]*dd+
                       sims.mat[,sprintf("psi.sp[%d]", ss)])
    
    c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
  }
  x.vals <- seq(from=1, to=5, length=5)
  
  y.vals <- lapply(x.vals, get.y.val)
  y.vals2_spp <- do.call(rbind, y.vals)
  
  lines(y.vals2_spp[,"mean"], type="l", ylim=c(0,1), col="grey")
  
  y.vals2_spp <- y.vals2_spp %>% as.data.frame() %>%
    dplyr::mutate(spp=ss, x=row_number())
  spp_lines[[ss]] <- y.vals2_spp
}
y.vals2_spp <- do.call(rbind, spp_lines)

# plot the mean community effect
get.y.val <- function(dd){
  chains <- plogis(sims.mat[,"mu.psi.0"]+
                     sims.mat[,sprintf("psi.yr[%d]", ss)]*dd+
                     sims.mat[,sprintf("psi.sp[%d]", ss)])
  
  c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
}
x.vals <- seq(from=1, to=5, length=5)
y.vals <- lapply(x.vals, get.y.val)
y.vals2_com <- do.call(rbind, y.vals) %>% as.data.frame()

lines(y.vals2_com[,"mean"], type="l", ylim=c(0,1), col="black", lwd=3)
lines(y.vals2_com[,"2.5%"], lty=2, lwd=3, col='black')
lines(y.vals2_com[,"97.5%"], lty=2, lwd=3, col='black')

y.vals2_spp <- y.vals2_spp %>%
  group_by(spp) %>%
  dplyr::mutate(min=dplyr::first(mean), max=dplyr::last(mean)) %>%
  dplyr::mutate(diff=max-min) %>%
  dplyr::mutate(sgn=sign(diff)) %>%
  ungroup()

## ggplot
M7_plot <- ggplot()+
  geom_line(y.vals2_spp, mapping=aes(x=x, y=mean, group=spp, color=diff), 
            alpha=0.8)+
  scale_color_continuous_divergingx(palette="RdYlBu", mid=0.0, name="Slope")+
  geom_line(y.vals2_com, mapping=aes(x=c(1:5), y=mean), size=1, linetype=1,
            color="grey10", alpha=0.8)+
  geom_line(y.vals2_com, mapping=aes(x=c(1:5), y=`2.5%`), size=1, linetype=2,
            color="grey10", alpha=0.8)+
  geom_line(y.vals2_com, mapping=aes(x=c(1:5), y=`97.5%`), size=1, linetype=2,
            color="grey10", alpha=0.8)+
  ylab("Occupancy Probability")+
  xlab("Era")+
  theme_cowplot()+
  background_grid()
M7_plot

## summary table
y.vals2_sum <- y.vals2_spp %>%
  dplyr::select()

# Detected-Range Model (MODEL 8)
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
     bty="n",
     las=1)

# plot a community of lines
spp_lines <- list()
for(ss in 1:195){
  get.y.val <- function(dd){
    chains <- plogis(sims.mat[,"mu.psi.0"]+
                      sims.mat[,sprintf("psi.yr[%d]", ss)]*dd+
                      sims.mat[,sprintf("psi.sp[%d]", ss)])
    
    c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
  }
  x.vals <- seq(from=1, to=5, length=5)
  
  y.vals <- lapply(x.vals, get.y.val)
  y.vals2_spp <- do.call(rbind, y.vals)
  
  lines(y.vals2_spp[,"mean"], type="l", ylim=c(0,1), col="grey")
  
  y.vals2_spp <- y.vals2_spp %>% as.data.frame() %>%
    dplyr::mutate(spp=ss, x=row_number())
  spp_lines[[ss]] <- y.vals2_spp
}
y.vals2_spp <- do.call(rbind, spp_lines)

# plot the mean community effect
get.y.val <- function(dd){
  chains <- plogis(sims.mat[,"mu.psi.0"]+
                    sims.mat[,sprintf("psi.yr[%d]", ss)]*dd+
                    sims.mat[,sprintf("psi.sp[%d]", ss)])
  
  c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
}
x.vals <- seq(from=1, to=5, length=5)
y.vals <- lapply(x.vals, get.y.val)
y.vals2_com <- do.call(rbind, y.vals) %>% as.data.frame()

lines(y.vals2_com[,"mean"], type="l", ylim=c(0,1), col="black", lwd=3)
lines(y.vals2_com[,"2.5%"], lty=2, lwd=3, col='black')
lines(y.vals2_com[,"97.5%"], lty=2, lwd=3, col='black')

y.vals2_spp <- y.vals2_spp %>%
  group_by(spp) %>%
  dplyr::mutate(min=dplyr::first(mean), max=dplyr::last(mean)) %>%
  dplyr::mutate(diff=max-min) %>%
  dplyr::mutate(sgn=sign(diff)) %>%
  ungroup()

## ggplot
M8_plot <- ggplot()+
  geom_line(y.vals2_spp, mapping=aes(x=x, y=mean, group=spp, color=diff), 
            alpha=0.8)+
  scale_color_continuous_divergingx(palette="RdYlBu", mid=0.0, name="Slope")+
  geom_line(y.vals2_com, mapping=aes(x=c(1:5), y=mean), size=1, linetype=1,
            color="grey10", alpha=0.8)+
  geom_line(y.vals2_com, mapping=aes(x=c(1:5), y=`2.5%`), size=1, linetype=2,
            color="grey10", alpha=0.8)+
  geom_line(y.vals2_com, mapping=aes(x=c(1:5), y=`97.5%`), size=1, linetype=2,
            color="grey10", alpha=0.8)+
  ylab("Occupancy Probability")+
  xlab("Era")+
  theme_cowplot()+
  background_grid()
M8_plot

## summary table
y.vals2_sum <- y.vals2_spp %>%
  dplyr::select(spp, diff, sgn) %>% unique()

sp_list <- sp_list %>%
  dplyr::mutate(spp=row_number())

model8.sum <- inner_join(sp_list, y.vals2_sum)

write.csv(model8.sum, "../output/speciesTrends.csv")

## combined figure
ggarrange(M5_plot, M7_plot, M8_plot, ncol=3, labels=c("(a)", "(b)"), common.legend=TRUE,
          legend="bottom")