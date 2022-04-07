# Load libraries
library(tidyverse); library(cowplot); library(ggpubr);
library(stringr); library(mcmcr); library(coda);
library(purrr); library(data.table); library(tidybayes);
library(modelr); library(colorspace); library(gridExtra)

# Load in .RDS files from Compute Canada runs
ode_all_range <- readRDS("../output/ODE_res_all_range.rds")
ode_det_all <- readRDS("../output/ODE_res_det_all.rds")
ode_det_range <- readRDS("../output/ODE_res_det_range.rds")
ode_det_range2 <- readRDS("../output/ODE_res_det_range_com.rds")

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
  x.vals <- seq(from=1, to=10, length=10)
  
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
                     sims.mat[,"mu.psi.yr"]*dd)
  
  c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
}
x.vals <- seq(from=1, to=10, length=10)
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
  scale_color_continuous_divergingx(palette="RdYlBu", mid=0.0, name="Total Occupancy Change",
                                    guide=guide_colorbar(
                                      direction = "horizontal",
                                      title.position = "top", barwidth=8, 
                                      label.theme=element_text(size=8, angle=0)
                                    ), limits=c(-0.4, 0.4))+
  geom_line(y.vals2_com, mapping=aes(x=c(1:10), y=mean), size=1, linetype=1,
            color="grey10", alpha=0.8)+
  geom_line(y.vals2_com, mapping=aes(x=c(1:10), y=`2.5%`), size=1, linetype=2,
            color="grey10", alpha=0.8)+
  geom_line(y.vals2_com, mapping=aes(x=c(1:10), y=`97.5%`), size=1, linetype=2,
            color="grey10", alpha=0.8)+
  ylab("Occupancy Probability")+
  xlab("")+
  scale_y_continuous(limits=c(0,1))+
  scale_x_continuous(breaks=c(1:10))+
  theme_cowplot()+
  theme(legend.position="bottom", legend.box = "horizontal")+
  background_grid()
M5_plot

diff.spp <- y.vals2_spp %>%
  dplyr::select(spp, diff) %>%
  unique() %>%
  dplyr::filter(diff >= 0.25 | diff <= -0.25) %>%
  nrow()
diff.spp

y.vals2_sum <- y.vals2_spp %>%
  dplyr::select(spp, diff, sgn) %>% unique()

M5_plot_inner <- ggplot()+
  geom_vline(xintercept=0, linetype=2)+
  geom_histogram(y.vals2_sum, mapping=aes(x=diff), binwidth=0.025)+
  xlab("Shift in Occupancy Probability")+
  ylab("Frequency")+
  theme_half_open()+
  theme(text = element_text(size = 6),
        axis.text = element_text(size = 6))

M5_plot <- M5_plot +
  annotation_custom(ggplotGrob(M5_plot_inner),
                    xmin=6, xmax=10,
                    ymin=0, ymax=0.5)

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
  x.vals <- seq(from=1, to=10, length=10)
  
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
                     sims.mat[,"mu.psi.yr"]*dd)
  
  c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
}
x.vals <- seq(from=1, to=10, length=10)
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
  scale_color_continuous_divergingx(palette="RdYlBu", mid=0.0, name="Total Occupancy Change",
                                    guide=guide_colorbar(
                                      direction = "horizontal",
                                      title.position = "top", barwidth=8, 
                                      label.theme=element_text(size=8, angle=0)
                                    ), limits=c(-0.4, 0.4))+
  geom_line(y.vals2_com, mapping=aes(x=c(1:10), y=mean), size=1, linetype=1,
            color="grey10", alpha=0.8)+
  geom_line(y.vals2_com, mapping=aes(x=c(1:10), y=`2.5%`), size=1, linetype=2,
            color="grey10", alpha=0.8)+
  geom_line(y.vals2_com, mapping=aes(x=c(1:10), y=`97.5%`), size=1, linetype=2,
            color="grey10", alpha=0.8)+
  scale_y_continuous(limits=c(0,1))+
  scale_x_continuous(breaks=c(1:10))+
  ylab("")+
  xlab("")+
  theme_cowplot()+
  background_grid()
M7_plot

## summary table
y.vals2_sum <- y.vals2_spp %>%
  dplyr::select()

y.vals2_sum <- y.vals2_spp %>%
  dplyr::select(spp, diff, sgn) %>% unique()

M7_plot_inner <- ggplot()+
  geom_vline(xintercept=0, linetype=2)+
  geom_histogram(y.vals2_sum, mapping=aes(x=diff), binwidth=0.025)+
  xlab("Shift in Occupancy Probability")+
  ylab("Frequency")+
  theme_half_open()+
  theme(text = element_text(size = 6),
        axis.text = element_text(size = 6))

M7_plot <- M7_plot +
  annotation_custom(ggplotGrob(M7_plot_inner),
                    xmin=6, xmax=10,
                    ymin=0, ymax=0.5)

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
  x.vals <- seq(from=1, to=10, length=10)
  
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
                     sims.mat[,"mu.psi.yr"]*dd)
  
  c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
}
x.vals <- seq(from=1, to=10, length=10)
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
  scale_color_continuous_divergingx(palette="RdYlBu", mid=0.0, name="Total Occupancy Change",
                                    guide=guide_colorbar(
                                      direction = "horizontal",
                                      title.position = "top", barwidth=8, 
                                      label.theme=element_text(size=8, angle=0)
                                    ), limits=c(-0.4, 0.4))+
  geom_line(y.vals2_com, mapping=aes(x=c(1:10), y=mean), size=1, linetype=1,
            color="grey10", alpha=0.8)+
  geom_line(y.vals2_com, mapping=aes(x=c(1:10), y=`2.5%`), size=1, linetype=2,
            color="grey10", alpha=0.8)+
  geom_line(y.vals2_com, mapping=aes(x=c(1:10), y=`97.5%`), size=1, linetype=2,
            color="grey10", alpha=0.8)+
  ylab("Occupancy Probability")+
  xlab("Era")+
  scale_y_continuous(limits=c(0,1))+
  scale_x_continuous(breaks=c(1:10))+
  theme_cowplot()+
  background_grid()
M8_plot

diff.spp <- y.vals2_spp %>%
  dplyr::select(spp, diff) %>%
  unique() %>%
  dplyr::filter(diff >= 0.25 | diff <= -0.25) %>%
  nrow()
diff.spp

y.vals2_sum <- y.vals2_spp %>%
  dplyr::select(spp, diff, sgn) %>% unique()

M8_plot_inner <- ggplot()+
  geom_vline(xintercept=0, linetype=2)+
  geom_histogram(y.vals2_sum, mapping=aes(x=diff), binwidth=0.025)+
  xlab("Shift in Occupancy Probability")+
  ylab("Frequency")+
  theme_half_open()+
  theme(text = element_text(size = 6),
        axis.text = element_text(size = 6))

M8_plot <- M8_plot +
  annotation_custom(ggplotGrob(M8_plot_inner),
                    xmin=6, xmax=10,
                    ymin=0, ymax=0.5)

# Detected-Range Model (MODEL 8)
samplech1 <- as.mcmc(ode_det_range2$chain1)
samplech2 <- as.mcmc(ode_det_range2$chain2)
samplech3 <- as.mcmc(ode_det_range2$chain3)

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
  x.vals <- seq(from=1, to=10, length=10)
  
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
                     sims.mat[,"mu.psi.yr"]*dd)
  
  c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
}
x.vals <- seq(from=1, to=10, length=10)
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
M8_plot2 <- ggplot()+
  geom_line(y.vals2_spp, mapping=aes(x=x, y=mean, group=spp, color=diff), 
            alpha=0.8)+
  scale_color_continuous_divergingx(palette="RdYlBu", mid=0.0, name="Total Occupancy Change",
                                    guide=guide_colorbar(
                                      direction = "horizontal",
                                      title.position = "top", barwidth=8, 
                                      label.theme=element_text(size=8, angle=0)
                                    ), limits=c(-0.4, 0.4))+
  geom_line(y.vals2_com, mapping=aes(x=c(1:10), y=mean), size=1, linetype=1,
            color="grey10", alpha=0.8)+
  geom_line(y.vals2_com, mapping=aes(x=c(1:10), y=`2.5%`), size=1, linetype=2,
            color="grey10", alpha=0.8)+
  geom_line(y.vals2_com, mapping=aes(x=c(1:10), y=`97.5%`), size=1, linetype=2,
            color="grey10", alpha=0.8)+
  ylab("")+
  xlab("Era")+
  scale_y_continuous(limits=c(0,1))+
  scale_x_continuous(breaks=c(1:10))+
  theme_cowplot()+
  background_grid()
M8_plot2

diff.spp <- y.vals2_spp %>%
  dplyr::select(spp, diff) %>%
  unique() %>%
  dplyr::filter(diff >= 0.25 | diff <= -0.25) %>%
  nrow()
diff.spp

## summary table
y.vals2_sum <- y.vals2_spp %>%
  dplyr::select(spp, diff, sgn) %>% unique()

M8_plot2_inner <- ggplot()+
  geom_vline(xintercept=0, linetype=2)+
  geom_histogram(y.vals2_sum, mapping=aes(x=diff), binwidth=0.025)+
  xlab("Shift in Occupancy Probability")+
  ylab("Frequency")+
  theme_half_open()+
  theme(text = element_text(size = 6),
        axis.text = element_text(size = 6))

M8_plot2 <- M8_plot2 +
  annotation_custom(ggplotGrob(M8_plot2_inner),
                    xmin=6, xmax=10,
                    ymin=0, ymax=0.5)

sp_list <- sp_list %>%
  dplyr::mutate(spp=row_number())

model8.sum <- inner_join(sp_list, y.vals2_sum)

write.csv(model8.sum, "../output/speciesTrends.csv")

## combined figure
blankPlot <- ggplot()+geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

odesPlot <- plot_grid(M5_plot+theme(legend.position="none"), 
                      M7_plot+theme(legend.position="none"), 
                      M8_plot+theme(legend.position="none"), 
                      M8_plot2+theme(legend.position="none"),
                      labels=c("a.", "b.", "c.", "d."), label_size=12)
bottomPlot <- plot_grid(blankPlot, get_legend(M8_plot2), blankPlot, blankPlot, ncol=4,
                        rel_widths=c(0.1, 1, 1, 1))

odesPlotFinal <- plot_grid(odesPlot, bottomPlot, nrow=2,
                           rel_heights=c(0.9,0.12))
odesPlotFinal

ggsave("../output/odesFigure.pdf", odesPlotFinal, dpi=350, height=8, width=8)