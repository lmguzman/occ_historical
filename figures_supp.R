library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(tidyr)
library(Metrics)
library(dplyr)

case <- 'p5'

output_p <- readRDS(paste0("all_outputs/p5_out.rds"))

fact <- nrow(dplyr::filter(output_p, term=="mu.psi.yr") %>% dplyr::select(term, estimate))

output_p[rhat > 1.1]

plot_mupsiyr <- function(mupsiyr){
  
  summary_1 <- output_p[mu.psi.yr == mupsiyr & term == 'mu.psi.yr', .(mean_es = mean(estimate), sd_es = sd(estimate)), by = .(visit_mod, p.yr, mu.v.yr, prop.visits.same)]
  
  summary_1 <- summary_1 %>% 
    mutate(p.yr = factor(p.yr, levels = c(-0.1, -0.05, 0, 0.05, 0.1)),
           mu.v.yr = factor(mu.v.yr, levels = c(-0.1, 0, 0.1)))
  
  zero_mu.psi.yr1 <- summary_1 %>% 
    filter(visit_mod == 'detectedno') %>% 
    ggplot() +
    geom_point(aes(x = prop.visits.same, y = mean_es)) +
    geom_errorbar(aes(x = prop.visits.same, ymin = mean_es -sd_es, ymax = mean_es +sd_es)) +
    geom_hline(yintercept = mupsiyr, colour = 'grey', linetype = 'dashed') +
    facet_grid(mu.v.yr ~ p.yr) +
    theme_cowplot() +
    theme(strip.background = element_blank()) +
    xlab(expression(rho[com])) + 
    ylab(expression("Estimated change in occupancy(" ~mu[psi~",era"]~")")) +
    theme(plot.margin = unit(c(1, 0.3, 0.1, 0.1), "cm"), 
          axis.title = element_text(size = 18),
          axis.text.y = element_text(size = 15))
  
  zero_mu.psi.yr1.1 <- ggdraw(zero_mu.psi.yr1) +
    draw_label(expression(p["era"]), 
               x = 0.5, y = 0.97, size = 18) +
    #draw_label(expression(mu[nu~",era"]), 
    #           x = 0.96, y = 0.5, angle = 270) +
    draw_label(expression("WF"["range,detected"]), 
               x = 0.1, y = 0.98, fontface = "bold", size = 18)
  
  zero_mu.psi.yr2 <- summary_1 %>% 
    filter(visit_mod == 'visitsno') %>% 
    ggplot() +
    geom_point(aes(x = prop.visits.same, y = mean_es)) +
    geom_errorbar(aes(x = prop.visits.same, ymin = mean_es -sd_es, ymax = mean_es +sd_es)) +
    geom_hline(yintercept = mupsiyr, colour = 'grey', linetype = 'dashed') +
    facet_grid(mu.v.yr ~ p.yr) +
    theme_cowplot() +
    theme(strip.background = element_blank()) +
    xlab(expression(rho[com])) + 
    #ylab(expression("Estimated change in occupancy(" ~mu[psi~",era"]~")")) +
    ylab("") +
    theme(plot.margin = unit(c(1, 1, 0.1, 0.1), "cm"),
          axis.title = element_text(size = 18),
          axis.text.y = element_text(size = 15)) 
  
  zero_mu.psi.yr2.1 <- ggdraw(zero_mu.psi.yr2) +
    draw_label(expression(p["era"]), 
               x = 0.5, y = 0.97, size = 18) +
    draw_label(expression(mu[nu~",era"]), 
               x = 0.97, y = 0.5, angle = 270, size = 18) +
    draw_label(expression("WF"["range,visits"]), 
               x = 0.1, y = 0.98, fontface = "bold", size = 18)
  
  zero_mu.psi.yr <- plot_grid(zero_mu.psi.yr1.1, zero_mu.psi.yr2.1)
  
  ggsave(zero_mu.psi.yr, filename = paste0("figures_clean/p5_mu.psi.yr_",mupsiyr,".pdf"), height = 10, width = 18)
  
  
}

plot_mupsiyr(mupsiyr = 0)

plot_mupsiyr(mupsiyr = -0.1)

plot_mupsiyr(mupsiyr = 0.1)

