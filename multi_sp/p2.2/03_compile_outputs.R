library(stringr)
library(mcmcr)
library(coda)
library(purrr)
library(dplyr)
library(data.table)

file_sim <- list.files("~/scratch/occ_historical/multi_sp/p2.2/outputs/model.res/")

compiled_res <- list()

for(i in 1:length(file_sim)){
  
  load(paste0("~/scratch/occ_historical/multi_sp/p2.2/outputs/model.res/",file_sim[i]), verbose = TRUE)
  
  file_params <- unlist(str_split(file_sim[i], "_"))
  
  samplech1 <- as.mcmc(res$chain1)
  samplech2 <- as.mcmc(res$chain2)
  samplech3 <- as.mcmc(res$chain3)
  
  coda.samples <- as.mcmc.list(samplech1, samplech2, samplech3)
  
  summary_mcmc <- coef(coda.samples, simplify = TRUE)
  rhat_mcmc <- map_df(rhat(coda.samples, by = "term"), ~data.frame(rhat = .x), .id = 'term')
  
  true_value <- map_df(sim.data[summary_mcmc$term], ~data.frame(true_val = .x), .id = 'term')
  
  compiled_res[[i]] <- full_join(full_join(summary_mcmc, rhat_mcmc), true_value) %>% 
    mutate(range_sim = file_params[15], range_mod = NA, visit_sim = file_params[18], visit_mod = file_params[1], 
           r = file_params[3], nyr = file_params[11], prop.visits.same = file_params[13], s = i) %>% 
    data.table()
  
}


all_outputs <- rbindlist(compiled_res)

saveRDS(all_outputs, 'p2.2/outputs/model.summary/all_outputs.rds')
