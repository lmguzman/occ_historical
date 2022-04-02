library(stringr)

n_runs <- 150



file_sim <- list.files("analysis/sim_range/outputs/model.res/")

model_res_name <- function(x){
  all_name <- unlist(str_split(x, "_"))
  
  data.frame(model = all_name[1], eras = all_name[2], r = as.numeric(all_name[4]), p.yr = as.numeric(all_name[6]),
             mu.psi.yr = as.numeric(all_name[8]), nyr = as.numeric(all_name[10]),
             prop.visits.same = as.numeric(all_name[12]), mu.v.yr = as.numeric(all_name[19]),
             s = as.numeric(str_remove(all_name[21], ".RData")))
}

model_res <- rbindlist(lapply(file_sim, model_res_name))

s_need <- c(1:n_runs)[!(1:n_runs %in% model_res[model == "visitsno" & eras == "eras5"]$s)]

file_sim2 <- list.files("analysis/sim_range/outputs/sim.data")

model_res_mis <- model_res[!s %in% s_need, .(r, nyr, p.yr, mu.psi.yr, prop.visits.same, mu.v.yr)] %>% unique()

news <- paste0("_r_",model_res_mis$r ,"_p.yr_",model_res_mis$p.yr ,"_mu.psi.yr_",model_res_mis$mu.psi.yr,"_nyr_",model_res_mis$nyr,"_prop.visits.same_",model_res_mis$prop.visits.same,"_range_polys_visit_visit_miss_","mu.v.yr_",model_res_mis$mu.v.yr)

#str_detect(file_sim, news)

s_to_run <- sapply(news, FUN = function(x) which(str_detect(file_sim2, x)))
names(s_to_run) <- NULL

c(1:n_runs)[!(1:n_runs %in% s_to_run)]

paste0(c(1:n_runs)[!(1:n_runs %in% s_to_run)], collapse = ',')