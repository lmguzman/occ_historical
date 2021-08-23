library(data.table)
library(ggplot2)
library(cowplot)
library(tidyr)
library(Metrics)
library(dplyr)


### main plot code ##

main_plot <- function(case, load_tf, intervals,...){
  
  if(load_tf == TRUE){
    output_p <- readRDS(paste0(case,"/outputs/model.summary/all_outputs.rds"))
    
    output_p$nyr <- factor(output_p$nyr, levels = c('2', '5', '10'))
    
    output_p$sim <- rep(1:(nrow(output_p)/11), each = 11)
    
    scenarios <- output_p[term %in% c('p.yr', 'mu.phi.yr', 'mu.gam.yr'), .(term, true_val, sim)] %>% 
      pivot_wider(names_from = 'term', values_from = 'true_val') 
    
    output_p_s <- output_p %>% 
      left_join(scenarios)
  }
  
  ## varying detection ##
  
  p_p.yr <- output_p_s[mu.phi.yr == 0 & mu.gam.yr == 0 & term == 'p.yr' & nyr != 2] %>% 
    ggplot() +
    geom_point(aes(x = p.yr, y = estimate, colour = nyr, group = r), position = position_dodge(width = 0.2)) +
    geom_linerange(aes(x = p.yr, ymin = lower, ymax = upper, colour = nyr, group = r), position = position_dodge(width = 0.2)) +
    theme_cowplot() +
    facet_wrap(~visit_mod) +
    geom_abline(intercept = 0, slope = 1, colour = 'grey', linetype = 'dashed') +
    theme(strip.background = element_blank()) +
    ylab('Estimated p.yr') + xlab('p.yr')
  
  p_mu.phi.yr <- output_p_s[mu.phi.yr == 0 & mu.gam.yr == 0 & term == 'mu.phi.yr' & nyr != 2] %>% 
    ggplot() +
    geom_point(aes(x = p.yr, y = estimate, colour = nyr, group = r), position = position_dodge(width = 0.2)) +
    geom_linerange(aes(x = p.yr, ymin = lower, ymax = upper, colour = nyr, group = r), position = position_dodge(width = 0.2)) +
    theme_cowplot() +
    facet_wrap(~visit_mod) +
    geom_hline(yintercept = 0,  colour = 'grey', linetype = 'dashed') +
    theme(strip.background = element_blank()) +
    ylab('Estimated mu.phi.yr') + xlab('p.yr')
  
  p_mu.gam.yr <- output_p_s[mu.phi.yr == 0 & mu.gam.yr == 0 & term == 'mu.gam.yr' & nyr != 2] %>% 
    ggplot() +
    geom_point(aes(x = p.yr, y = estimate, colour = nyr, group = r), position = position_dodge(width = 0.2)) +
    geom_linerange(aes(x = p.yr, ymin = lower, ymax = upper, colour = nyr, group = r), position = position_dodge(width = 0.2)) +
    theme_cowplot() +
    facet_wrap(~visit_mod) +
    geom_hline(yintercept = 0,  colour = 'grey', linetype = 'dashed') +
    theme(strip.background = element_blank()) +
    ylab('Estimated mu.gam.yr') + xlab('p.yr')
  
  p_v_p <- plot_grid(p_p.yr, p_mu.phi.yr, p_mu.gam.yr, ncol = 1)
  
  ggsave(p_v_p, filename = paste0("figures/",case,"_v_p.jpeg"), height = 12)
  
  
  
  
  ## varying persistence ##
  
  p_p.yr_2 <- output_p_s[p.yr == 0 & mu.gam.yr == 0 & term == 'p.yr' & nyr != 2] %>% 
    ggplot() +
    geom_point(aes(x = mu.phi.yr, y = estimate, colour = nyr, group = r), position = position_dodge(width = 0.2)) +
    geom_linerange(aes(x = mu.phi.yr, ymin = lower, ymax = upper, colour = nyr, group = r), position = position_dodge(width = 0.2)) +
    theme_cowplot() +
    facet_wrap(~visit_mod) +
    geom_hline(yintercept = 0,  colour = 'grey', linetype = 'dashed') +
    theme(strip.background = element_blank()) +
    ylab('Estimated p.yr') + xlab('mu.phi.yr')
  
  p_mu.phi.yr_2 <- output_p_s[p.yr == 0 & mu.gam.yr == 0 & term == 'mu.phi.yr' & nyr != 2] %>% 
    ggplot() +
    geom_point(aes(x = mu.phi.yr, y = estimate, colour = nyr, group = r), position = position_dodge(width = 0.2)) +
    geom_linerange(aes(x = mu.phi.yr, ymin = lower, ymax = upper, colour = nyr, group = r), position = position_dodge(width = 0.2)) +
    theme_cowplot() +
    facet_wrap(~visit_mod) +
    geom_abline(intercept = 0, slope = 1, colour = 'grey', linetype = 'dashed') +
    theme(strip.background = element_blank()) +
    ylab('Estimated mu.phi.yr') + xlab('mu.phi.yr')
  
  p_mu.gam.yr_2 <- output_p_s[p.yr == 0 & mu.gam.yr == 0 & term == 'mu.gam.yr' & nyr != 2] %>% 
    ggplot() +
    geom_point(aes(x = mu.phi.yr, y = estimate, colour = nyr, group = r), position = position_dodge(width = 0.2)) +
    geom_linerange(aes(x = mu.phi.yr, ymin = lower, ymax = upper, colour = nyr, group = r), position = position_dodge(width = 0.2)) +
    theme_cowplot() +
    facet_wrap(~visit_mod) +
    geom_hline(yintercept = 0,  colour = 'grey', linetype = 'dashed') +
    theme(strip.background = element_blank()) +
    ylab('Estimated mu.gam.yr') + xlab('mu.phi.yr')
  
  p_v_mu.phi.yr <- plot_grid(p_p.yr_2, p_mu.phi.yr_2, p_mu.gam.yr_2, ncol = 1)
  
  ggsave(p_v_mu.phi.yr, filename = paste0("figures/",case,"_v_mu.phi.yr.jpeg"), height = 12)
  
  
  ## varying colonization ##
  
  p_p.yr_3 <- output_p_s[p.yr == 0 & mu.phi.yr == 0 & term == 'p.yr' & nyr != 2] %>% 
    ggplot() +
    geom_point(aes(x = mu.gam.yr, y = estimate, colour = nyr, group = r), position = position_dodge(width = 0.2)) +
    geom_linerange(aes(x = mu.gam.yr, ymin = lower, ymax = upper, colour = nyr, group = r), position = position_dodge(width = 0.2)) +
    theme_cowplot() +
    facet_wrap(~visit_mod) +
    geom_hline(yintercept = 0,  colour = 'grey', linetype = 'dashed') +
    theme(strip.background = element_blank()) +
    ylab('Estimated p.yr') + xlab('mu.gam.yr')
  
  p_mu.phi.yr_3 <- output_p_s[p.yr == 0 & mu.phi.yr == 0 & term == 'mu.phi.yr' & nyr != 2] %>% 
    ggplot() +
    geom_point(aes(x = mu.gam.yr, y = estimate, colour = nyr, group = r), position = position_dodge(width = 0.2)) +
    geom_linerange(aes(x = mu.gam.yr, ymin = lower, ymax = upper, colour = nyr, group = r), position = position_dodge(width = 0.2)) +
    theme_cowplot() +
    facet_wrap(~visit_mod) +
    geom_hline(yintercept = 0,  colour = 'grey', linetype = 'dashed') +
    theme(strip.background = element_blank()) +
    ylab('Estimated mu.phi.yr') + xlab('mu.gam.yr')
  
  p_mu.gam.yr_3 <- output_p_s[p.yr == 0 & mu.phi.yr == 0 & term == 'mu.gam.yr' & nyr != 2] %>% 
    ggplot() +
    geom_point(aes(x = mu.gam.yr, y = estimate, colour = nyr, group = r), position = position_dodge(width = 0.2)) +
    geom_linerange(aes(x = mu.gam.yr, ymin = lower, ymax = upper, colour = nyr, group = r), position = position_dodge(width = 0.2)) +
    theme_cowplot() +
    facet_wrap(~visit_mod) +
    geom_abline(intercept = 0, slope = 1, colour = 'grey', linetype = 'dashed') +
    theme(strip.background = element_blank()) +
    ylab('Estimated mu.gam.yr') + xlab('mu.gam.yr')
  
  p_v_mu.gam.yr <- plot_grid(p_p.yr_3, p_mu.phi.yr_3, p_mu.gam.yr_3, ncol = 1)
  
  ggsave(p_v_mu.gam.yr, filename = paste0("figures/",case,"_v_mu.gam.yr.jpeg"), height = 12)
  
}



######## p1 ########
# All species are everywhere, no ranges
# visits the same for all species
# MS_all_all, MS_all_visits vs MS_all_detected 

main_plot(case = 'p1', load_tf = TRUE)


######## P2 ########
# All species are everywhere, no ranges
# visits are different for all species
# MS_all_all, MS_all_visits vs MS_all_detected 

main_plot(case = 'p2', load_tf = TRUE)


######## P2.1 ########
# All species are everywhere, no ranges
# visits are different for all species
# MS_all_all, MS_all_visits vs MS_all_detected 
# Same as P2 but mu.v.yr = 0 <- the number of visits through time does not decrease

main_plot(case = 'p2.1', load_tf = TRUE)


######## P2.2 ########
# All species are everywhere, no ranges
# visits vary dependding on the proportion of visits that are the same 
# MS_all_all, MS_all_visits vs MS_all_detected vs MS_all_community
# MS_all_community is where only community visits are used to infer non-detections


case <- 'p2.2'

output_p <- readRDS(paste0(case,"/outputs/model.summary/all_outputs.rds"))

output_p$nyr <- factor(output_p$nyr, levels = c('5', '10'))

output_p$s <- rep(1:(nrow(output_p)/11), each = 11)

output_p %>% head()

scenarios <- output_p[term %in% c('p.yr', 'mu.phi.yr', 'mu.gam.yr'), .(term, true_val, s)] %>% 
  pivot_wider(names_from = 'term', values_from = 'true_val') 

output_p_s <- output_p %>% 
  left_join(scenarios)

summarised_rmse <- output_p_s[,.(rmse_vals = rmse(true_val, estimate), N = .N), by = .(visit_mod, term, nyr, prop.visits.same)]

summarised_rmse %>% 
  filter(term %in% c('mu.gam.yr', 'mu.phi.yr', 'p.yr')) %>% 
  ggplot(aes(x = prop.visits.same, y = rmse_vals, colour =  visit_mod, group = visit_mod)) + 
  geom_point() +
  geom_line() +
  facet_grid(nyr~term) +
  theme_cowplot()

##### by prop.visits.same ###

output_p.v.same.0 <- output_p_s[prop.visits.same == 0]

main_plot(case = '2.2_p.v.same.0', load_tf = FALSE, output_p_s = output_p.v.same.0)

output_p.v.same.1 <- output_p_s[prop.visits.same == 1]

main_plot(case = '2.2_p.v.same.1', load_tf = FALSE, output_p_s = output_p.v.same.1)

output_p.v.same.0.5 <- output_p_s[prop.visits.same == 0.5]

main_plot(case = '2.2_p.v.same.0.5', load_tf = FALSE, output_p_s = output_p.v.same.0.5)

output_p_s[estimate > 10]
