library(dplyr)
library(cowplot)
library(ggplot2)
library(viridis)
library(data.table)
library(Metrics)
library(stringr)
library(ggthemes)



###### Main effect data processing workflow x interaction with Main effect of prob of community visits

p2<- readRDS("all_outputs/p2_out.rds")

summarised_rmse <- p2[mu.v.yr == 0 & term %in% c("mu.psi.yr", "p.yr") & nyr == 10][,.(rmse_vals = rmse(true_val, estimate), N = .N), 
                                                    by = .(visit_mod, term, nyr, prop.visits.same)] %>%
  dplyr::mutate(nyr=paste(nyr, "eras"))

summarised_rmse$term <- factor(summarised_rmse$term, 
                               levels=c("mu.psi.yr", "p.yr"),
                               labels=c(expression(mu[psi~",era"]),expression(p["era"])))

plot_wf_prob_com <- summarised_rmse %>% 
  ggplot(aes(x = prop.visits.same, y = rmse_vals, colour = visit_mod, group = visit_mod)) +
  geom_point() +
  geom_line() +
  facet_wrap(~term, labeller=label_parsed) +
  theme_cowplot() +
  xlab(expression(rho["com"])) +
  ylab("RMSE") +
  theme(strip.background = element_blank(),
        legend.position = "bottom") +
  geom_point(alpha=0.9) +
  geom_line(alpha=0.9) +
  scale_color_viridis_d(end=0.8, name="", 
                        labels = c(expression("WF"["all,all"]), 
                                   expression("WF"["all,detected"]),
                                   expression("WF"["all,visits"])))
  


  
##### Effect of nyr ####

p2<- readRDS("all_outputs/p2_out.rds")


summarised_rmse <- p2[mu.v.yr == 0 & term %in% c("mu.psi.yr", "p.yr")][,.(rmse_vals = rmse(true_val, estimate), N = .N), 
                                                                                   by = .(visit_mod, term, nyr, prop.visits.same)]  %>% 
  mutate(prop.visits.same = as.numeric(prop.visits.same))


summarised_rmse$nyr <- factor(summarised_rmse$nyr, 
                              levels = c("2", "5", "10"), 
                              labels = c(expression("2"~" Eras"),expression("5"~" Eras"),expression("10"~" Eras")))

summarised_rmse$term <- factor(summarised_rmse$term, 
                        levels=c("mu.psi.yr", "p.yr"),
                        labels=c(expression(mu[psi~",era"]),expression(p["era"])))


data_text <- data.frame(label = c("e.", "a.", "c.", "f.", "b.", 'd.'),  # Create data for text
                        x = c(-0.04),
                        y = c(0.38)) %>% 
  bind_cols(expand.grid(nyr = unique(summarised_rmse$nyr), term = unique(summarised_rmse$term)))


plot_wf_prob_com_nyr <-  ggplot() +
  geom_point(data = summarised_rmse, aes(x = prop.visits.same, y = rmse_vals, colour = visit_mod, group = visit_mod)) +
  geom_line(data = summarised_rmse, aes(x = prop.visits.same, y = rmse_vals, colour = visit_mod, group = visit_mod)) +
  facet_grid(nyr~term, labeller=label_parsed, scales = "free") +
  xlab(expression("Probability that visits are community samples vs. targeted samples" ~ rho["com"])) +
  ylab("RMSE") +
  geom_point(alpha=0.9) +
  geom_line(alpha=0.9) +
  scale_color_viridis_d(end=0.8, name="", 
                        labels = c(expression("WF"["all,all"]), 
                                   expression("WF"["all,detected"]),
                                   expression("WF"["all,visits"]))) +
  theme_few()+
  theme(strip.background = element_blank(),
        legend.position = "bottom", 
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 17),
        strip.text = element_text(size = 20), 
        legend.text = element_text(size = 20),
        axis.line=element_line()) + 
  scale_x_continuous(limits=c(-0.05,1)) + 
  scale_y_continuous(limits=c(0,0.4)) +
  geom_text(data = data_text, aes(x = x, y = y,label = label), size = 6)

plot_wf_prob_com_nyr

ggsave(plot_wf_prob_com_nyr, filename = "figures_clean/nyr_pcom.pdf", height = 8)


####### effect of range #####

p2<- readRDS("all_outputs/p2_out.rds")

summarised_rmse_p2 <- p2[mu.v.yr == 0 & term %in% c("mu.psi.yr", "p.yr") & nyr == 5][,.(rmse_vals = rmse(true_val, estimate), N = .N), 
                                                                                   by = .(visit_mod, term, nyr, prop.visits.same)] %>% 
  mutate(scenario = "p2")


p4_ignored<- readRDS("all_outputs/p4_all_uncensored.rds")

summarised_rmse_p4_ig <- p4_ignored[visit_sim == 0 & term %in% c("mu.psi.yr", "p.yr") & nyr ==5][,.(rmse_vals = rmse(true_val, estimate), N = .N), 
                                                                                      by = .(visit_mod, term, nyr, prop.visits.same)] %>% 
  mutate(scenario = "p4_ig")


p4_cen<- readRDS("all_outputs/p4_all_censored.rds")

summarised_rmse_p4_cen <- p4_cen[visit_sim == 0 & term %in% c("mu.psi.yr", "p.yr") & nyr == 5][,.(rmse_vals = rmse(true_val, estimate), N = .N), 
                                                                                                   by = .(visit_mod, term, nyr, prop.visits.same)] %>% 
  mutate(scenario = "p4_cen")


range_together <- bind_rows(summarised_rmse_p2, summarised_rmse_p4_ig, summarised_rmse_p4_cen) %>% 
  mutate(visit_mod = str_remove(visit_mod, "no|yes")) %>% 
  mutate(prop.visits.same = as.numeric(prop.visits.same))

range_together$term <- factor(range_together$term, 
                               levels=c("mu.psi.yr", "p.yr"),
                               labels=c(expression(mu[psi~",era"]),expression(p["era"])))

range_together$scenario <- factor(range_together$scenario, 
                              levels=c("p2", "p4_ig", "p4_cen"),
                              labels=c(expression("Sim"["all"]~"WF"["all,X"]), expression("Sim"["range"]~"WF"["all,X"]), expression("WF"["range,X"])))

data_text <- data.frame(label = c("a.", "c.", "e.", "b.", "d.", 'f.'),  # Create data for text
                        x = c(-0.04),
                        y = c(0.38)) %>% 
  bind_cols(expand.grid(scenario = unique(range_together$scenario), term = unique(range_together$term)))


plot_wf_prob_com_range <- ggplot() +
  geom_point(data = range_together, aes(x = prop.visits.same, y = rmse_vals, colour = visit_mod, group = visit_mod)) +
  geom_line(data = range_together, aes(x = prop.visits.same, y = rmse_vals, colour = visit_mod, group = visit_mod)) +
  facet_grid(scenario~term, labeller=label_parsed) +
  xlab(expression("Probability that visits are community samples vs. targeted samples" ~ rho["com"])) +
  ylab("RMSE") +
  geom_point(alpha=0.9) +
  geom_line(alpha=0.9) +
  scale_color_viridis_d(end=0.8, name="", 
                        labels = c(expression("WF"["X,all"]), 
                                   expression("WF"["X,detected"]),
                                   expression("WF"["X,visits"]))) +
  theme_few()+
  theme(strip.background = element_blank(),
        legend.position = "bottom", 
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 17),
        strip.text = element_text(size = 20), 
        legend.text = element_text(size = 20),
        axis.line=element_line()) +
  scale_x_continuous(limits=c(-0.05,1)) + 
  scale_y_continuous(limits=c(0,0.4)) +
  geom_text(data = data_text, aes(x = x, y = y,label = label), size = 6)

plot_wf_prob_com_range

ggsave(plot_wf_prob_com_range, filename = "figures_clean/range_pcom.pdf", height = 8)


###### Effect of mu.v.yr


p4_cen<- readRDS("all_outputs/p4_all_censored.rds")

summarised_rmse_p4_cen <- p4_cen[term %in% c("mu.psi.yr", "p.yr") & nyr == 10][,.(rmse_vals = rmse(true_val, estimate), N = .N), 
                                                                                               by = .(visit_mod, term, visit_sim, prop.visits.same)] %>% 
  mutate(prop.visits.same = as.numeric(prop.visits.same))


summarised_rmse_p4_cen$term <- factor(summarised_rmse_p4_cen$term, 
                              levels=c("mu.psi.yr", "p.yr"),
                              labels=c(expression(mu[psi~",era"]),expression(p["era"])))

summarised_rmse_p4_cen$visit_sim <- factor(summarised_rmse_p4_cen$visit_sim,
                                  levels=c("-0.1", "0", "0.1"),
                                  labels=c(expression(mu[nu~",era"]~" = -0.1"),expression(mu[nu~",era"]~" = 0"),expression(mu[nu~",era"]~" = 0.1")))


data_text <- data.frame(label = c("a.", "e.", "c.", "b.", "f.", 'd.'),  # Create data for text
                        x = c(-0.04),
                        y = c(0.42)) %>% 
  bind_cols(expand.grid(visit_sim = unique(summarised_rmse_p4_cen$visit_sim), term = unique(summarised_rmse_p4_cen$term)))


plot_wf_prob_com_visit <-  ggplot() +
  geom_point(data = summarised_rmse_p4_cen, aes(x = prop.visits.same, y = rmse_vals, colour = visit_mod, group = visit_mod)) +
  geom_line(data = summarised_rmse_p4_cen, aes(x = prop.visits.same, y = rmse_vals, colour = visit_mod, group = visit_mod)) +
  facet_grid(visit_sim~term, labeller=label_parsed) +
  xlab(expression("Probability that visits are community samples vs. targeted samples" ~ rho["com"])) +
  ylab("RMSE") +
  geom_point(alpha=0.9) +
  geom_line(alpha=0.9) +
  scale_color_viridis_d(end=0.8, name="", 
                        labels = c(expression("WF"["range,all"]), 
                                   expression("WF"["range,detected"]),
                                   expression("WF"["range,visits"]))) +
  theme_few()+
  theme(strip.background = element_blank(),
        legend.position = "bottom", 
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 17),
        strip.text = element_text(size = 20), 
        legend.text = element_text(size = 20),
        axis.line=element_line()) +
  scale_x_continuous(limits=c(-0.05,1)) + 
  scale_y_continuous(limits=c(0,0.45)) +
  geom_text(data = data_text, aes(x = x, y = y,label = label), size = 6)

plot_wf_prob_com_visit

ggsave(plot_wf_prob_com_visit, filename = "figures_clean/visit_pcom.pdf", height = 8)

