library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(tidyr)
library(Metrics)
library(dplyr)
library(stringr)
### main plot code ##
main_plot <- function(case, censor, load_tf, intervals,...){
  
   case <- "p2.2"
  if(load_tf == TRUE){
    output_p <- readRDS(paste0("multi_sp/",case,"/outputs/model.summary/all_outputs.rds"))
    print("Main file loading successful...")
    
    output_p1 <- dplyr::filter(output_p, mu.v.yr==-0.1)  
    output_p1$eras <- factor(output_p1$eras, levels = c('eras2', 'eras5', 'eras10'))
    output_p1$sim <- rep(1:(nrow(output_p1)/108), each = 108)
    scenarios1 <- output_p1[term %in% c('p.yr', 'mu.psi.yr'), .(term, true_val, sim)] %>% 
      pivot_wider(names_from = 'term', values_from = 'true_val') 
    output_p_s1 <- output_p1 %>% 
      left_join(scenarios1)
    print("First file parsed...")
    
    output_p2 <- dplyr::filter(output_p, mu.v.yr==0)  
    output_p2$eras <- factor(output_p2$eras, levels = c('eras2', 'eras5', 'eras10'))
    output_p2$sim <- rep(1:(nrow(output_p2)/108), each =108)
    scenarios2 <- output_p2[term %in% c('p.yr', 'mu.psi.yr'), .(term, true_val, sim)] %>% 
      pivot_wider(names_from = 'term', values_from = 'true_val') 
    output_p_s2 <- output_p2 %>% 
      left_join(scenarios2)
    print("First file parsed...")
    
    output_p3 <- dplyr::filter(output_p, mu.v.yr==0.1)  
    output_p2$eras <- factor(output_p2$eras, levels = c('eras2', 'eras5', 'eras10'))
    output_p3$sim <- rep(1:(nrow(output_p3)/108), each = 108)
    scenarios3 <- output_p3[term %in% c('p.yr', 'mu.psi.yr'), .(term, true_val, sim)] %>% 
      pivot_wider(names_from = 'term', values_from = 'true_val') 
    output_p_s3 <- output_p3 %>% 
      left_join(scenarios3)
    print("First file parsed...")
    
  }
  
  ## RMSE Values ##
  summarised_rmse <- output_p_s1[,.(rmse_vals = rmse(true_val, estimate), N = .N), 
                                 by = .(visit_mod, term, eras, prop.visits.same)] %>%
    dplyr::mutate(eras=paste(str_extract(eras, "\\d+"), 'eras'),
                  visit_mod=case_when(visit_mod=="allno" ~ "All",
                                      visit_mod=="detectedno" ~ "Detected",
                                      visit_mod=="visitsno" ~ "True Visits"))
  summarised_rmse$eras <- factor(summarised_rmse$eras, levels=c("2 eras",
                                                              "5 eras",
                                                              "10 eras"))
  summarised_rmse$visit_mod <- factor(summarised_rmse$visit_mod, 
                                      levels=c("All",
                                               "Detected",
                                               "True Visits"))
  
  saveRDS(summarised_rmse, "multi_sp/p2.2/outputs/model.summary/RMSE_summary_1.rds")
  
  p1 <- summarised_rmse %>% 
    filter(term %in% c('mu.psi.yr', 'p.yr'), !is.na(visit_mod)) %>% 
    ggplot(aes(x = prop.visits.same, y = rmse_vals, colour =  visit_mod, 
               group = visit_mod)) + 
    geom_point(alpha=0.9) +
    geom_line(alpha=0.9) +
    scale_color_viridis_d(end=0.8, name="Type of Modelling Approach")+
    geom_hline(yintercept=0, linetype=2)+
    xlab("Proportion of Community Visits")+
    ylab("RMSE")+
    facet_grid(eras~term, scales = 'free_y') +
    theme_cowplot()
  
  summarised_rmse2 <- output_p_s2[,.(rmse_vals = rmse(true_val, estimate), N = .N), 
                                  by = .(visit_mod, term, eras, prop.visits.same)] %>%
    dplyr::mutate(eras=paste(str_extract(eras, "\\d+"), 'eras'),
                  visit_mod=case_when(visit_mod=="allno" ~ "All",
                                      visit_mod=="detectedno" ~ "Detected",
                                      visit_mod=="visitsno" ~ "True Visits"))
  summarised_rmse2$eras <- factor(summarised_rmse2$eras, levels=c("2 eras",
                                                                "5 eras",
                                                                "10 eras"))
  summarised_rmse2$visit_mod <- factor(summarised_rmse2$visit_mod, 
                                       levels=c("All",
                                                "Detected",
                                                "True Visits"))
  p2 <- summarised_rmse2 %>% 
    filter(term %in% c('mu.psi.yr', 'p.yr'), !is.na(visit_mod)) %>% 
    ggplot(aes(x = prop.visits.same, y = rmse_vals, colour =  visit_mod, 
               group = visit_mod)) + 
    geom_point(alpha=0.9) +
    geom_line(alpha=0.9) +
    scale_color_viridis_d(end=0.8, name="Type of Modelling Approach")+
    geom_hline(yintercept=0, linetype=2)+
    xlab("Proportion of Community Visits")+
    ylab("RMSE")+
    facet_grid(eras~term, scales = 'free_y') +
    theme_cowplot()
  
  summarised_rmse3 <- output_p_s3[,.(rmse_vals = rmse(true_val, estimate), N = .N), 
                                  by = .(visit_mod, term, eras, prop.visits.same)] %>%
    dplyr::mutate(eras=paste(str_extract(eras, "\\d+"), 'eras'),
                  visit_mod=case_when(visit_mod=="allno" ~ "All",
                                      visit_mod=="detectedno" ~ "Detected",
                                      visit_mod=="visitsno" ~ "True Visits"))
  summarised_rmse3$eras <- factor(summarised_rmse3$eras, levels=c("2 eras",
                                                                "5 eras",
                                                                "10 eras"))
  summarised_rmse3$visit_mod <- factor(summarised_rmse3$visit_mod, 
                                       levels=c("All",
                                                "Detected",
                                                "True Visits"))
  p3 <- summarised_rmse3 %>% 
    filter(term %in% c('mu.psi.yr', 'p.yr'), !is.na(visit_mod)) %>% 
    ggplot(aes(x = prop.visits.same, y = rmse_vals, colour =  visit_mod, 
               group = visit_mod)) + 
    geom_point(alpha=0.9) +
    geom_line(alpha=0.9) +
    scale_color_viridis_d(end=0.8, name="Type of Modelling Approach")+
    geom_hline(yintercept=0, linetype=2)+
    xlab("Proportion of Community Visits")+
    ylab("RMSE")+
    facet_grid(eras~term, scales = 'free_y') +
    theme_cowplot()
  
  p_all <- ggarrange(p1, p2, p3, ncol=3, common.legend=TRUE, labels=c("(a)", "(b)", "(c)"),
            legend="bottom")
  
  ggsave(p_all, filename = "multi_sp/Figure_1_RMSE.pdf")
}

######## P4 ########
# Species have ranges but we don't restrict our
# indices to them.
main_plot(case = 'P4', censor="all_uncensored_outputs", load_tf = TRUE)



summarised_rmse %>% 
  filter(visit_mod  !='All') %>% 
  filter(term %in% c('mu.psi.yr', 'p.yr'), !is.na(visit_mod)) %>% 
  ggplot(aes(x = prop.visits.same, y = rmse_vals, colour =  visit_mod, 
             group = visit_mod)) + 
  geom_point(alpha=0.9) +
  geom_line(alpha=0.9) +
  scale_color_viridis_d(end=0.8, name="Type of Modelling Approach")+
  geom_hline(yintercept=0, linetype=2)+
  xlab("Proportion of Community Visits")+
  ylab("RMSE")+
  facet_grid(eras~term) +
  theme_cowplot()

