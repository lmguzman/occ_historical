### main plot code ##
main_plot <- function(case, censor, load_tf, intervals,...){
  
  if(load_tf == TRUE){
    output_p <- readRDS(paste0("multi_sp/",case,"/outputs/model.summary/all_outputs.rds"))
    print("Main file loading successful...")
    
    output_p1 <- dplyr::filter(output_p, mu.v.yr==-0.1)  
    output_p1$nyr <- factor(output_p1$nyr, levels = c('2', '5', '10'))
    output_p1$sim <- rep(1:(nrow(output_p1)/8), each = 8)
    scenarios1 <- output_p1[term %in% c('p.yr', 'mu.psi.yr'), .(term, true_val, sim)] %>% 
      pivot_wider(names_from = 'term', values_from = 'true_val') 
    output_p_s1 <- output_p1 %>% 
      left_join(scenarios1)
    print("First file parsed...")
    
    output_p2 <- dplyr::filter(output_p, mu.v.yr==0)  
    output_p2$nyr <- factor(output_p2$nyr, levels = c('2', '5', '10'))
    output_p2$sim <- rep(1:(nrow(output_p2)/8), each = 8)
    scenarios2 <- output_p2[term %in% c('p.yr', 'mu.psi.yr'), .(term, true_val, sim)] %>% 
      pivot_wider(names_from = 'term', values_from = 'true_val') 
    output_p_s2 <- output_p2 %>% 
      left_join(scenarios2)
    print("First file parsed...")
    
    output_p3 <- dplyr::filter(output_p, mu.v.yr==0.1)  
    output_p3$nyr <- factor(output_p3$nyr, levels = c('2', '5', '10'))
    output_p3$sim <- rep(1:(nrow(output_p3)/8), each = 8)
    scenarios3 <- output_p3[term %in% c('p.yr', 'mu.psi.yr'), .(term, true_val, sim)] %>% 
      pivot_wider(names_from = 'term', values_from = 'true_val') 
    output_p_s3 <- output_p3 %>% 
      left_join(scenarios3)
    print("First file parsed...")
    
  }
  
  ## RMSE Values ##
  summarised_rmse <- output_p_s1[,.(rmse_vals = rmse(true_val, estimate), N = .N), 
                                 by = .(visit_mod, term, nyr, prop.visits.same)] %>%
    dplyr::mutate(nyr=paste(nyr, "eras"),
                  visit_mod=case_when(visit_mod=="allno" ~ "All",
                                      visit_mod=="detectedno" ~ "Detected",
                                      visit_mod=="visitsno" ~ "True Visits"))
  summarised_rmse$nyr <- factor(summarised_rmse$nyr, levels=c("2 eras",
                                                              "5 eras",
                                                              "10 eras"))
  summarised_rmse$visit_mod <- factor(summarised_rmse$visit_mod, 
                                      levels=c("All",
                                               "Detected",
                                               "True Visits"))
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
    facet_grid(nyr~term) +
    theme_cowplot()
  
  summarised_rmse2 <- output_p_s2[,.(rmse_vals = rmse(true_val, estimate), N = .N), 
                                  by = .(visit_mod, term, nyr, prop.visits.same)] %>%
    dplyr::mutate(nyr=paste(nyr, "eras"),
                  visit_mod=case_when(visit_mod=="allno" ~ "All",
                                      visit_mod=="detectedno" ~ "Detected",
                                      visit_mod=="visitsno" ~ "True Visits"))
  summarised_rmse2$nyr <- factor(summarised_rmse2$nyr, levels=c("2 eras",
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
    facet_grid(nyr~term) +
    theme_cowplot()
  
  summarised_rmse3 <- output_p_s3[,.(rmse_vals = rmse(true_val, estimate), N = .N), 
                                  by = .(visit_mod, term, nyr, prop.visits.same)] %>%
    dplyr::mutate(nyr=paste(nyr, "eras"),
                  visit_mod=case_when(visit_mod=="allno" ~ "All",
                                      visit_mod=="detectedno" ~ "Detected",
                                      visit_mod=="visitsno" ~ "True Visits"))
  summarised_rmse3$nyr <- factor(summarised_rmse3$nyr, levels=c("2 eras",
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
    facet_grid(nyr~term) +
    theme_cowplot()
  
  p_all <- ggarrange(p1, p2, p3, ncol=3, common.legend=TRUE, labels=c("(a)", "(b)", "(c)"),
            legend="bottom")
  
  ggsave(p_all, filename = "multi_sp/Figure_1_RMSE.pdf")
}

######## P4 ########
# Species have ranges but we don't restrict our
# indices to them.
main_plot(case = 'P4', censor="all_uncensored_outputs", load_tf = TRUE)
