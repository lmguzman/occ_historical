library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(tidyr)
library(Metrics)
library(dplyr)

######## P2 ########
p2_output <- readRDS("fin_sim/p2_out.rds")

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
print("Second file parsed...")

output_p3 <- dplyr::filter(output_p, mu.v.yr==0.1)  
output_p3$nyr <- factor(output_p3$nyr, levels = c('2', '5', '10'))
output_p3$sim <- rep(1:(nrow(output_p3)/8), each = 8)
scenarios3 <- output_p3[term %in% c('p.yr', 'mu.psi.yr'), .(term, true_val, sim)] %>% 
  pivot_wider(names_from = 'term', values_from = 'true_val') 
output_p_s3 <- output_p3 %>% 
  left_join(scenarios3)
print("Third file parsed...")

summarised_rmse <- output_p_s1[,.(rmse_vals = rmse(true_val, estimate), N = .N), 
                               by = .(visit_mod, term, nyr, prop.visits.same)] %>%
  dplyr::mutate(nyr=paste(nyr, "eras"),
                visit_mod=case_when(visit_mod=="allno" ~ "MODEL 1",
                                    visit_mod=="detectedno" ~ "MODEL 2",
                                    visit_mod=="visitsno" ~ "MODEL 3"))
summarised_rmse$nyr <- factor(summarised_rmse$nyr, levels=c("2 eras",
                                                            "5 eras",
                                                            "10 eras"))
summarised_rmse$visit_mod <- factor(summarised_rmse$visit_mod, 
                                    levels=c("MODEL 1",
                                             "MODEL 2",
                                             "MODEL 3"))
p1 <- summarised_rmse %>% 
  filter(term %in% c('mu.psi.yr'), !is.na(visit_mod)) %>% 
  ggplot(aes(x=prop.visits.same, y=rmse_vals, colour=nyr, 
             group = nyr)) + 
  geom_point(alpha=0.9) +
  geom_line(alpha=0.9) +
  scale_color_viridis_d(end=0.8, name="Number of Eras in Model")+
  geom_hline(yintercept=0, linetype=2)+
  xlab("Probability of Community Visits")+
  ylab("RMSE")+
  facet_grid(visit_mod~term) +
  theme_cowplot()

summarised_rmse <- output_p_s2[,.(rmse_vals = rmse(true_val, estimate), N = .N), 
                               by = .(visit_mod, term, nyr, prop.visits.same)] %>%
  dplyr::mutate(nyr=paste(nyr, "eras"),
                visit_mod=case_when(visit_mod=="allno" ~ "MODEL 1",
                                    visit_mod=="detectedno" ~ "MODEL 2",
                                    visit_mod=="visitsno" ~ "MODEL 3"))
summarised_rmse$nyr <- factor(summarised_rmse$nyr, levels=c("2 eras",
                                                            "5 eras",
                                                            "10 eras"))
summarised_rmse$visit_mod <- factor(summarised_rmse$visit_mod, 
                                    levels=c("MODEL 1",
                                             "MODEL 2",
                                             "MODEL 3"))
p2 <- summarised_rmse %>% 
  filter(term %in% c('mu.psi.yr'), !is.na(visit_mod)) %>% 
  ggplot(aes(x=prop.visits.same, y=rmse_vals, colour=nyr, 
             group = nyr)) + 
  geom_point(alpha=0.9) +
  geom_line(alpha=0.9) +
  scale_color_viridis_d(end=0.8, name="Number of Eras in Model")+
  geom_hline(yintercept=0, linetype=2)+
  xlab("Probability of Community Visits")+
  ylab("RMSE")+
  facet_grid(visit_mod~term) +
  theme_cowplot()

summarised_rmse <- output_p_s3[,.(rmse_vals = rmse(true_val, estimate), N = .N), 
                               by = .(visit_mod, term, nyr, prop.visits.same)] %>%
  dplyr::mutate(nyr=paste(nyr, "eras"),
                visit_mod=case_when(visit_mod=="allno" ~ "MODEL 1",
                                    visit_mod=="detectedno" ~ "MODEL 2",
                                    visit_mod=="visitsno" ~ "MODEL 3"))
summarised_rmse$nyr <- factor(summarised_rmse$nyr, levels=c("2 eras",
                                                            "5 eras",
                                                            "10 eras"))
summarised_rmse$visit_mod <- factor(summarised_rmse$visit_mod, 
                                    levels=c("MODEL 1",
                                             "MODEL 2",
                                             "MODEL 3"))
p3 <- summarised_rmse %>% 
  filter(term %in% c('mu.psi.yr'), !is.na(visit_mod)) %>% 
  ggplot(aes(x=prop.visits.same, y=rmse_vals, colour=nyr, 
             group = nyr)) + 
  geom_point(alpha=0.9) +
  geom_line(alpha=0.9) +
  scale_color_viridis_d(end=0.8, name="Number of Eras in Model")+
  geom_hline(yintercept=0, linetype=2)+
  xlab("Probability of Community Visits")+
  ylab("RMSE")+
  facet_grid(visit_mod~term) +
  theme_cowplot()

p2_plot <- ggarrange(p1, p2, p3, ncol=3, common.legend=TRUE, legend="bottom",
          labels=c("(a)", "(b)", "(c)"))
p2_plot
ggsave("../figures/p2_plot.pdf", p2_plot, dpi=350)

######## P4 ########
output_p <- readRDS("fin_sim/p4_out.rds")

output_p1 <- dplyr::filter(output_p, mu.v.yr==-0.1)
output_p1$nyr <- factor(output_p1$nyr, levels = c('2', '5', '10'))
output_p1$sim <- rep(1:(nrow(output_p1)/8), each = 8) %>% as.character()
scenarios1 <- output_p1[term %in% c('p.yr', 'mu.psi.yr'), .(term, true_val, sim)] %>% 
  pivot_wider(names_from = 'term', values_from = 'true_val')  %>%
  dplyr::mutate(across(everything(), as.character))
output_p_s1 <- output_p1 %>% 
  left_join(scenarios1)
print("First file parsed...")

output_p2 <- dplyr::filter(output_p, mu.v.yr==0)  %>%
  dplyr::mutate(across(everything(), as.character))  
output_p2$nyr <- factor(output_p2$nyr, levels = c('2', '5', '10'))
output_p2$sim <- rep(1:(nrow(output_p2)/8), each = 8)
scenarios2 <- output_p2[term %in% c('p.yr', 'mu.psi.yr'), .(term, true_val, sim)] %>% 
  pivot_wider(names_from = 'term', values_from = 'true_val') 
output_p_s2 <- output_p2 %>% 
  left_join(scenarios2)
print("Second file parsed...")

output_p3 <- dplyr::filter(output_p, mu.v.yr==0.1) %>%
  dplyr::mutate(across(everything(), as.character))
output_p3$nyr <- factor(output_p3$nyr, levels = c('2', '5', '10'))
output_p3$sim <- rep(1:(nrow(output_p3)/8), each = 8)
scenarios3 <- output_p3[term %in% c('p.yr', 'mu.psi.yr'), .(term, true_val, sim)] %>% 
  pivot_wider(names_from = 'term', values_from = 'true_val') 
output_p_s3 <- output_p3 %>% 
  left_join(scenarios3)
print("Third file parsed...")

summarised_rmse <- output_p_s1[,.(rmse_vals = rmse(true_val, estimate), N = .N), 
                               by = .(visit_mod, term, nyr, prop.visits.same)] %>%
  dplyr::mutate(nyr=paste(nyr, "eras"),
                visit_mod=case_when(visit_mod=="allno" ~ "MODEL 1",
                                    visit_mod=="detectedno" ~ "MODEL 2",
                                    visit_mod=="visitsno" ~ "MODEL 3"))
summarised_rmse$nyr <- factor(summarised_rmse$nyr, levels=c("2 eras",
                                                            "5 eras",
                                                            "10 eras"))
summarised_rmse$visit_mod <- factor(summarised_rmse$visit_mod, 
                                    levels=c("MODEL 1",
                                             "MODEL 2",
                                             "MODEL 3"))
p1 <- summarised_rmse %>% 
  filter(term %in% c('mu.psi.yr'), !is.na(visit_mod)) %>% 
  ggplot(aes(x=prop.visits.same, y=rmse_vals, colour=nyr, 
             group = nyr)) + 
  geom_point(alpha=0.9) +
  geom_line(alpha=0.9) +
  scale_color_viridis_d(end=0.8, name="Number of Eras in Model")+
  geom_hline(yintercept=0, linetype=2)+
  xlab("Probability of Community Visits")+
  ylab("RMSE")+
  facet_grid(visit_mod~term) +
  theme_cowplot()

summarised_rmse <- output_p_s2[,.(rmse_vals = rmse(true_val, estimate), N = .N), 
                               by = .(visit_mod, term, nyr, prop.visits.same)] %>%
  dplyr::mutate(nyr=paste(nyr, "eras"),
                visit_mod=case_when(visit_mod=="allno" ~ "MODEL 1",
                                    visit_mod=="detectedno" ~ "MODEL 2",
                                    visit_mod=="visitsno" ~ "MODEL 3"))
summarised_rmse$nyr <- factor(summarised_rmse$nyr, levels=c("2 eras",
                                                            "5 eras",
                                                            "10 eras"))
summarised_rmse$visit_mod <- factor(summarised_rmse$visit_mod, 
                                    levels=c("MODEL 1",
                                             "MODEL 2",
                                             "MODEL 3"))
p2 <- summarised_rmse %>% 
  filter(term %in% c('mu.psi.yr'), !is.na(visit_mod)) %>% 
  ggplot(aes(x=prop.visits.same, y=rmse_vals, colour=nyr, 
             group = nyr)) + 
  geom_point(alpha=0.9) +
  geom_line(alpha=0.9) +
  scale_color_viridis_d(end=0.8, name="Number of Eras in Model")+
  geom_hline(yintercept=0, linetype=2)+
  xlab("Probability of Community Visits")+
  ylab("RMSE")+
  facet_grid(visit_mod~term) +
  theme_cowplot()

summarised_rmse <- output_p_s3[,.(rmse_vals = rmse(true_val, estimate), N = .N), 
                               by = .(visit_mod, term, nyr, prop.visits.same)] %>%
  dplyr::mutate(nyr=paste(nyr, "eras"),
                visit_mod=case_when(visit_mod=="allno" ~ "MODEL 1",
                                    visit_mod=="detectedno" ~ "MODEL 2",
                                    visit_mod=="visitsno" ~ "MODEL 3"))
summarised_rmse$nyr <- factor(summarised_rmse$nyr, levels=c("2 eras",
                                                            "5 eras",
                                                            "10 eras"))
summarised_rmse$visit_mod <- factor(summarised_rmse$visit_mod, 
                                    levels=c("MODEL 1",
                                             "MODEL 2",
                                             "MODEL 3"))
p3 <- summarised_rmse %>% 
  filter(term %in% c('mu.psi.yr'), !is.na(visit_mod)) %>% 
  ggplot(aes(x=prop.visits.same, y=rmse_vals, colour=nyr, 
             group = nyr)) + 
  geom_point(alpha=0.9) +
  geom_line(alpha=0.9) +
  scale_color_viridis_d(end=0.8, name="Number of Eras in Model")+
  geom_hline(yintercept=0, linetype=2)+
  xlab("Probability of Community Visits")+
  ylab("RMSE")+
  facet_grid(visit_mod~term) +
  theme_cowplot()

p4.1_plot <- ggarrange(p1, p2, p3, ncol=3, common.legend=TRUE, legend="bottom",
                     labels=c("(a)", "(b)", "(c)"))
p4.1_plot
ggsave("../figures/p4.1_plot.pdf", p2_plot, dpi=350)

p2_output <- readRDS("fin_sim/p4.2_out.rds")

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
print("Second file parsed...")

output_p3 <- dplyr::filter(output_p, mu.v.yr==0.1)  
output_p3$nyr <- factor(output_p3$nyr, levels = c('2', '5', '10'))
output_p3$sim <- rep(1:(nrow(output_p3)/8), each = 8)
scenarios3 <- output_p3[term %in% c('p.yr', 'mu.psi.yr'), .(term, true_val, sim)] %>% 
  pivot_wider(names_from = 'term', values_from = 'true_val') 
output_p_s3 <- output_p3 %>% 
  left_join(scenarios3)
print("Third file parsed...")

summarised_rmse <- output_p_s1[,.(rmse_vals = rmse(true_val, estimate), N = .N), 
                               by = .(visit_mod, term, nyr, prop.visits.same)] %>%
  dplyr::mutate(nyr=paste(nyr, "eras"),
                visit_mod=case_when(visit_mod=="allno" ~ "MODEL 1",
                                    visit_mod=="detectedno" ~ "MODEL 2",
                                    visit_mod=="visitsno" ~ "MODEL 3"))
summarised_rmse$nyr <- factor(summarised_rmse$nyr, levels=c("2 eras",
                                                            "5 eras",
                                                            "10 eras"))
summarised_rmse$visit_mod <- factor(summarised_rmse$visit_mod, 
                                    levels=c("MODEL 1",
                                             "MODEL 2",
                                             "MODEL 3"))
p1 <- summarised_rmse %>% 
  filter(term %in% c('mu.psi.yr'), !is.na(visit_mod)) %>% 
  ggplot(aes(x=prop.visits.same, y=rmse_vals, colour=nyr, 
             group = nyr)) + 
  geom_point(alpha=0.9) +
  geom_line(alpha=0.9) +
  scale_color_viridis_d(end=0.8, name="Number of Eras in Model")+
  geom_hline(yintercept=0, linetype=2)+
  xlab("Probability of Community Visits")+
  ylab("RMSE")+
  facet_grid(visit_mod~term) +
  theme_cowplot()

summarised_rmse <- output_p_s2[,.(rmse_vals = rmse(true_val, estimate), N = .N), 
                               by = .(visit_mod, term, nyr, prop.visits.same)] %>%
  dplyr::mutate(nyr=paste(nyr, "eras"),
                visit_mod=case_when(visit_mod=="allno" ~ "MODEL 1",
                                    visit_mod=="detectedno" ~ "MODEL 2",
                                    visit_mod=="visitsno" ~ "MODEL 3"))
summarised_rmse$nyr <- factor(summarised_rmse$nyr, levels=c("2 eras",
                                                            "5 eras",
                                                            "10 eras"))
summarised_rmse$visit_mod <- factor(summarised_rmse$visit_mod, 
                                    levels=c("MODEL 1",
                                             "MODEL 2",
                                             "MODEL 3"))
p2 <- summarised_rmse %>% 
  filter(term %in% c('mu.psi.yr'), !is.na(visit_mod)) %>% 
  ggplot(aes(x=prop.visits.same, y=rmse_vals, colour=nyr, 
             group = nyr)) + 
  geom_point(alpha=0.9) +
  geom_line(alpha=0.9) +
  scale_color_viridis_d(end=0.8, name="Number of Eras in Model")+
  geom_hline(yintercept=0, linetype=2)+
  xlab("Probability of Community Visits")+
  ylab("RMSE")+
  facet_grid(visit_mod~term) +
  theme_cowplot()

summarised_rmse <- output_p_s3[,.(rmse_vals = rmse(true_val, estimate), N = .N), 
                               by = .(visit_mod, term, nyr, prop.visits.same)] %>%
  dplyr::mutate(nyr=paste(nyr, "eras"),
                visit_mod=case_when(visit_mod=="allno" ~ "MODEL 1",
                                    visit_mod=="detectedno" ~ "MODEL 2",
                                    visit_mod=="visitsno" ~ "MODEL 3"))
summarised_rmse$nyr <- factor(summarised_rmse$nyr, levels=c("2 eras",
                                                            "5 eras",
                                                            "10 eras"))
summarised_rmse$visit_mod <- factor(summarised_rmse$visit_mod, 
                                    levels=c("MODEL 1",
                                             "MODEL 2",
                                             "MODEL 3"))
p3 <- summarised_rmse %>% 
  filter(term %in% c('mu.psi.yr'), !is.na(visit_mod)) %>% 
  ggplot(aes(x=prop.visits.same, y=rmse_vals, colour=nyr, 
             group = nyr)) + 
  geom_point(alpha=0.9) +
  geom_line(alpha=0.9) +
  scale_color_viridis_d(end=0.8, name="Number of Eras in Model")+
  geom_hline(yintercept=0, linetype=2)+
  xlab("Probability of Community Visits")+
  ylab("RMSE")+
  facet_grid(visit_mod~term) +
  theme_cowplot()

p4.2_plot <- ggarrange(p1, p2, p3, ncol=3, common.legend=TRUE, legend="bottom",
                       labels=c("(a)", "(b)", "(c)"))
p4.2_plot
ggsave("../figures/p4.2_plot.pdf", p2_plot, dpi=350)
  
ggarrange(p2_plot, p4.1_plot, p4.2_plot, nrow=3, common.legend=TRUE)