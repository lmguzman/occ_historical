library(stringr)
library(mcmcr)
library(coda)
library(purrr)
library(dplyr)
library(data.table)

compiled_res <- list()

for(i in 1:length(file_sim)){
  
  load(paste0("sim_all/outputs/model.res/model.res/",file_sim[i]), verbose = TRUE)
  
  file_params <- unlist(str_split(file_sim[i], "_"))
  
  samplech1 <- as.mcmc(res$chain1)
  samplech2 <- as.mcmc(res$chain2)
  samplech3 <- as.mcmc(res$chain3)
  
  coda.samples <- as.mcmc.list(samplech1, samplech2, samplech3)
  
  summary_mcmc <- coef(coda.samples, simplify = TRUE)
  
  ## rhat
  
  rhat_int <- rhat(coda.samples, by = "term")
  
  rhat_mcmc_1 <- map_df(rhat_int[summary_mcmc$term[!str_detect(summary_mcmc$term, "psi.yr\\[|psi.sp\\[")]], ~data.frame(rhat = .x), .id = 'term')
  
  rhat_psi.yr <- data.frame(term = paste0("psi.yr[",1:length(data.prepped$my.info$sp.keep) ,"]" ), rhat = rhat_int$psi.yr)
  
  rhat_psi.sp <- data.frame(term = paste0("psi.sp[", 1:length(data.prepped$my.info$sp.keep),"]" ), rhat = rhat_int$psi.sp)
  
  rhat_mcmc <- bind_rows(rhat_mcmc_1, rhat_psi.yr, rhat_psi.sp)
  
  ## true values
  
  sp_kept <- as.numeric(str_remove(data.prepped$my.info$sp.keep, "^sp_0+"))
  
  true_value_1 <- map_df(sim.data[summary_mcmc$term[!str_detect(summary_mcmc$term, "psi.yr\\[|psi.sp\\[")]], ~data.frame(true_val = .x), .id = 'term')
  
  true_value_psi.yr <- data.frame(term = paste0("psi.yr[", 1:length(data.prepped$my.info$sp.keep),"]" ), true_val = sim.data$psi.yr[sp_kept])
  
  true_value_psi.sp <- data.frame(term = paste0("psi.sp[", 1:length(data.prepped$my.info$sp.keep),"]" ), true_val = sim.data$psi.sp[sp_kept])
  
  true_value <- bind_rows(true_value_1, true_value_psi.yr, true_value_psi.sp)
  
  compiled_res[[i]] <- full_join(full_join(summary_mcmc, rhat_mcmc), true_value) %>% 
    mutate(visit_sim = file_params[16], visit_mod = file_params[1], eras = file_params[2],
           r = file_params[4], nyr = file_params[10], prop.visits.same = file_params[12], mu.v.yr = file_params[19],
           s = i) %>% 
    data.table()
  
}

all_outputs <- rbindlist(compiled_res)

saveRDS(all_outputs, 'sim_range/outputs/model.summary/sim_all_outputs.rds')