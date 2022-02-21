library(stringr)
library(mcmcr)
library(coda)
library(purrr)
library(dplyr)
library(data.table)

file_sim <- list.files("~/scratch/occ_historical/analysis/sim_all/outputs/model.res/")

compiled_res <- list()

for(i in 1:length(file_sim)){
  
  load(paste0("~/scratch/occ_historical/analysis/sim_all/outputs/model.res/",file_sim[i]), verbose = TRUE)
  
  file_params <- unlist(str_split(file_sim[i], "_"))
  
  samplech1 <- as.mcmc(res$chain1)
  samplech2 <- as.mcmc(res$chain2)
  samplech3 <- as.mcmc(res$chain3)
  
  coda.samples <- as.mcmc.list(samplech1, samplech2, samplech3)
  
  summary_mcmc <- coef(coda.samples, simplify = TRUE)
  
  ## rhat
  
  rhat_int <- rhat(coda.samples, by = "term")
  
  rhat_mcmc_1 <- map_df(rhat_int[summary_mcmc$term[!str_detect(summary_mcmc$term, "psi.yr\\[|psi.sp\\[")]], ~data.frame(rhat = .x), .id = 'term')
  
  rhat_psi.yr <- data.frame(term = paste0("psi.yr[", 1:length(sim.data$psi.yr),"]" ), rhat = rhat_int$psi.yr)
  
  rhat_psi.sp <- data.frame(term = paste0("psi.sp[", 1:length(sim.data$psi.yr),"]" ), rhat = rhat_int$psi.sp)
  
  rhat_mcmc <- bind_rows(rhat_mcmc_1, rhat_psi.yr, rhat_psi.sp)
  
  ## true values
  
  true_value_1 <- map_df(sim.data[summary_mcmc$term[!str_detect(summary_mcmc$term, "psi.yr\\[|psi.sp\\[")]], ~data.frame(true_val = .x), .id = 'term')
  
  true_value_psi.yr <- data.frame(term = paste0("psi.yr[", 1:length(sim.data$psi.yr),"]" ), true_val = sim.data$psi.yr)
  
  true_value_psi.sp <- data.frame(term = paste0("psi.sp[", 1:length(sim.data$psi.yr),"]" ), true_val = sim.data$psi.sp)
  
  true_value <- bind_rows(true_value_1, true_value_psi.yr, true_value_psi.sp)
  
  compiled_res[[i]] <- full_join(full_join(summary_mcmc, rhat_mcmc), true_value) %>% 
    mutate(visit_sim = file_params[18], visit_mod = file_params[1], eras = file_params[2],
           r = file_params[4], nyr = file_params[10], prop.visits.same = file_params[12], mu.v.yr = file_params[14],
           s = i) %>% 
    data.table()
  
}


all_outputs <- rbindlist(compiled_res)

saveRDS(all_outputs, '~/scratch/occ_historical/analysis/sim_all/outputs/model.summary/all_outputs.rds')







### compile species rarity

library(stringr)
library(mcmcr)
library(coda)
library(purrr)
library(dplyr)
library(data.table)

file_sim <- list.files("~/scratch/occ_historical/analysis/sim_all/outputs/model.res/")

compiled_res <- list()

for(i in 1:length(file_sim)){
  
  load(paste0("~/scratch/occ_historical/analysis/sim_all/outputs/model.res/",file_sim[i]), verbose = TRUE)

  file_params <- unlist(str_split(file_sim[i], "_"))
  
  n_site <- apply(sim.data$Z, c(1,3), sum)
  
  compiled_res[[i]] <- data.frame(n_site) %>% 
    mutate(visit_sim = file_params[18], visit_mod = file_params[1], eras = file_params[2],
           r = file_params[4], nyr = file_params[10], prop.visits.same = file_params[12], mu.v.yr = file_params[14],
           s = i) %>% 
    data.table()
  
}


all_outputs <- rbindlist(compiled_res)

saveRDS(all_outputs, '~/scratch/occ_historical/analysis/sim_all/outputs/model.summary/species_orig_occupancy.rds')




### find missing runs ##

library(stringr)

file_sim <- list.files("~/scratch/occ_historical/analysis/sim_all/outputs/model.res/")

x <- file_sim[1]

model_res_name <- function(x){
  all_name <- unlist(str_split(x, "_"))
  
  data.frame(model = all_name[1], eras = all_name[2], r = as.numeric(all_name[4]), p.yr = as.numeric(all_name[6]),
             mu.psi.yr = as.numeric(all_name[8]), nyr = as.numeric(all_name[10]),
             prop.visits.same = as.numeric(all_name[12]), mu.v.yr = as.numeric(all_name[14]),
             s = as.numeric(str_remove(all_name[21], ".RData")))
}

model_res <- rbindlist(lapply(file_sim, model_res_name))

s_need <- c(1:150)[!(1:150 %in% model_res[model == 'detectedno' & eras == 'eras10']$s)]

file_sim <- list.files("p2.2/outputs/sim.data")

model_res_mis <- model_res[!s %in% s_need, .(r, nyr, p.yr, mu.psi.yr, prop.visits.same, mu.v.yr)] %>% unique()

news <- paste0("_r_",model_res_mis$r ,"_p.yr_",model_res_mis$p.yr ,"_mu.psi.yr_",model_res_mis$mu.psi.yr,"_nyr_",model_res_mis$nyr,"_prop.visits.same_",model_res_mis$prop.visits.same,"_mu.v.yr_",model_res_mis$mu.v.yr,"_range_all_visit_visit_miss_")

#str_detect(file_sim, news)

s_to_run <- sapply(news, FUN = function(x) which(str_detect(file_sim, x)))
names(s_to_run) <- NULL

c(1:150)[!(1:150 %in% s_to_run)]

paste0(c(1:150)[!(1:150 %in% s_to_run)], collapse = ',')

## 5 eras

# visitsno # detectedno

4,5,12,13,75,81,99,100,109,125,130

## 10 eras

## detectedno

1,3,5,6,7,8,9,10,11,12,20,21,46,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,106,107,108,109,110,121,122,123,124,125,126,127,128
