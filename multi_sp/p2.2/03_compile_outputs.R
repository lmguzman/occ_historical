library(stringr)
library(mcmcr)
library(coda)
library(purrr)
library(dplyr)
library(data.table)

file_sim <- list.files("~/scratch/occ_historical/multi_sp/p2.2/outputs/model.res/")
#file_sim <- list.files("multi_sp/p2.2/outputs/model.res/")

compiled_res <- list()

for(i in 1:length(file_sim)){
  
  load(paste0("~/scratch/occ_historical/multi_sp/p2.2/outputs/model.res/",file_sim[i]), verbose = TRUE)
  #load(paste0("multi_sp/p2.2/outputs/model.res/",file_sim[i]), verbose = TRUE)
  
  file_params <- unlist(str_split(file_sim[i], "_"))
  
  samplech1 <- as.mcmc(res$chain1)
  samplech2 <- as.mcmc(res$chain2)
  samplech3 <- as.mcmc(res$chain3)
  
  coda.samples <- as.mcmc.list(samplech1, samplech2, samplech3)
  
  summary_mcmc <- coef(coda.samples, simplify = TRUE)
  rhat_mcmc <- map_df(rhat(coda.samples, by = "term"), ~data.frame(rhat = .x), .id = 'term')
  
  true_value <- map_df(sim.data[summary_mcmc$term], ~data.frame(true_val = .x), .id = 'term')
  
  compiled_res[[i]] <- full_join(full_join(summary_mcmc, rhat_mcmc), true_value) %>% 
    mutate(visit_sim = file_params[17], visit_mod = file_params[1], 
           r = file_params[3], nyr = file_params[9], prop.visits.same = file_params[11], mu.v.yr = file_params[13],
           s = i) %>% 
    data.table()
  
}


all_outputs <- rbindlist(compiled_res)

saveRDS(all_outputs, '~/scratch/occ_historical/multi_sp/p2.2/outputs/model.summary/all_outputs.rds')




x <- file_sim[1]

model_res_name <- function(x){
  all_name <- unlist(str_split(x, "_"))
  
  data.frame(model = all_name[1], r = as.numeric(all_name[3]), p.yr = as.numeric(all_name[5]),
             mu.psi.yr = as.numeric(all_name[7]), nyr = as.numeric(all_name[9]),
             prop.visits.same = as.numeric(all_name[11]), mu.v.yr = as.numeric(all_name[13]),
             s = as.numeric(str_remove(all_name[20], ".RData")))
}

model_res <- rbindlist(lapply(file_sim, model_res_name))


r <- data.frame(r = 1:5)
nyr <- data.frame(nyr = c(2, 5, 10))
prop.visits.same <- data.frame(prop.visits.same = c(0,0.25,0.5, 0.75, 1))
mu.v.yr <- data.frame(mu.v.yr = c(-0.1,0,0.1))


model_res$model %>% unique()

paste0(c(1:225)[!(1:225 %in% model_res[model == 'detectedno']$s)], collapse = ',')

7,11,12,13,14,15,22,26,27,28,29,30,42,44,45,57,59,60,82,86,87,90,97,101,102,103,104,105,117,119,157,161,162,163,164,165,176,177,178,179,191,192,194,195


#### real 
"_r_2_p.yr_0_mu.psi.yr_0_nyr_10_prop.visits.same_0_mu.v.yr_-0.1_range_all_visit_visit_miss_"

model_res %>% 
  group_by(model) %>% 
  summarise(n = n())

s_need <- model_res %>% 
  group_by(s) %>% 
  summarise(n = n()) %>% 
  filter(n<2) %>% 
  dplyr::select(s) %>% unlist()

model_res_mis <- model_res[s %in% s_need, .(r, nyr, prop.visits.same, mu.v.yr)] %>% unique()

news <- paste0("_r_",model_res_mis$r,"_p.yr_0_mu.psi.yr_0_nyr_",model_res_mis$nyr,"_prop.visits.same_",model_res_mis$prop.visits.same,"_mu.v.yr_",model_res_mis$mu.v.yr,"_range_all_visit_visit_miss_")
which(str_detect(file_sim, news))
s_to_run <- sapply(news, FUN = function(x) which(str_detect(file_sim, x)))
names(s_to_run) <- NULL

paste0(s_to_run, collapse = ",")
[1] "12,1,2,3,4,5,6,8,9,10,39,46,47,48,49,50,51,52,53,54,55,56,58,63,64,65,66,67,91,92,93,94,95,96,98,99,121,122,123,125,126,127,128,129,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,168,169,170,181,182,183,184,185,186,187,188,189,190,193,211,213,214,215,218"


## for ones that never finished 


model_res %>% 
  group_by(model) %>% 
  summarise(n = n())

s_need <- c(1:1875)[!(1:1875 %in% model_res[model == 'visitsno']$s)]

file_sim <- list.files("../sim.data")

model_res_mis <- model_res[!s %in% s_need, .(r, nyr, p.yr, mu.psi.yr, prop.visits.same, mu.v.yr)] %>% unique()

news <- paste0("_r_",model_res_mis$r,"_p.yr_",model_res_mis$p.yr ,"_mu.psi.yr_",model_res_mis$mu.psi.yr,"_nyr_",model_res_mis$nyr,"_prop.visits.same_",model_res_mis$prop.visits.same,"_mu.v.yr_",model_res_mis$mu.v.yr,"_range_polys_visit_visit_miss_")

#str_detect(file_sim, news)

s_to_run <- sapply(news, FUN = function(x) which(str_detect(file_sim, x)))
names(s_to_run) <- NULL

c(1:1875)[!(1:1875 %in% s_to_run)]

paste0(c(1:1875)[!(1:1875 %in% s_to_run)], collapse = ',')

## visitsno

45,55,86,87,90,92,93,94,96,97,125,126,127,133,134,135,143,144,145,146,173,174,175,176,190,191,192,206,245,275,284,303,312,313,314,315,318,372,373,374,396,411,412,413,414,415,418,419,441,466,467,511,512,513,514,717,718,788,797,825,846,866,867,915,920,1115,1116,1187,1250,1251,1252,1253,1338,1378,1392,1529,1530,1531,1538,1546,1553,1556,1557,1558,1572,1575,1609,1610,1611,1612,1613,1654,1662,1664,1685,1688,1691,1736,1737,1740,1768,1772,1790,1793,1796,1797,1801,1818,1851,1858,1870,187

## detectedno

852,945,1019,1020,1021,1022,1030,1031,1043,1064,1114,1115,1124,1137,1235,1241,1242,1243,1287,1288,1289,1331,1332,1333,1360,1375,1384,1389,1416,1418,1421,1425,1454,1464,1481,1482,1486,1487,1488,1490,1491,1527,1536,1543,1544,1545,1548,1549,1574,1575,1579,1604,1607,1609,1610,1643,1653,1656,1666,1667,1668,1669,1692,1711,1714,1716,1745,1746,1747,1748,1749,1750,1753,1754,1755,1756,1757,1758,1763,1787,1788,1789,1790,1791,1799,1800,1801,1812,1820

