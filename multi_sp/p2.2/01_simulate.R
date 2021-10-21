###### Part 1 #####

## 1. 
# Simulation:
# - All species are everywhere, no ranges
# - Missing visits but same visit for all species 

# Modeling: 
# MS_all_all, MS_all_visits vs MS_all_detected for more than 2 eras 
# this would be useful for a scenario that has small spatial scale 

library(parallel)
library(dplyr)
source('multi_sp/simulation/src/simulate_ms.R')

expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))


## Make a data frame of the parameters 

## Change either p.yr, or mu.psi.yr

base_scenario <- data.frame(nsp          = 50,
                             nsite        = 100,
                             nvisit       = 3,
                             ## detection
                             mu.p.0 = -0.5,
                             sigma.p.site = 1.5,
                             sigma.p.sp = 0.5,
                             ## occupancy
                             mu.psi.0 = 0,
                             sigma.psi.sp = 0.5,
                             sigma.psi.yr = 0.2,
                             ## visit
                             mu.v.0 = 0, 
                             mu.v.yr = -0.1,
                             ## type sym
                             type.range = 'all',
                             type.visit = 'visit_miss')

r <- data.frame(r = c(1,2))
p.yr <- data.frame(p.yr = c(-0.1, -0.05, 0, 0.05, 0.1))
mu.psi.yr <- data.frame(mu.psi.yr = c(-0.1, -0.05, 0, 0.05, 0.1))
nyr <- data.frame(nyr = c(2,5,10))
prop.visits.same <- data.frame(prop.visits.same = c(0,0.25,0.5,0.75,1))

all_scenarios <- bind_rows(expand.grid.df(base_scenario, p.yr, data.frame(mu.psi.yr = 0), r, nyr, prop.visits.same),
                           expand.grid.df(base_scenario, mu.psi.yr, data.frame(p.yr = 0), r, nyr, prop.visits.same)) %>% 
  unique()


run_all_simulation_1 <- function(s, all_scenarios){

  
  sim.data <- make.data(## data structure set up 
                        nsp = all_scenarios$nsp[s],
                        nsite = all_scenarios$nsite[s], 
                        nyr = all_scenarios$nyr[s], 
                        nvisit = all_scenarios$nvisit[s], 
                        ## detection
                        mu.p.0 = all_scenarios$mu.p.0[s],
                        p.yr = all_scenarios$p.yr[s],
                        sigma.p.site = all_scenarios$sigma.p.site[s],
                        sigma.p.sp = all_scenarios$sigma.p.sp[s],
                        ## occupancy
                        mu.psi.0 = all_scenarios$mu.psi.0[s],
                        sigma.psi.sp = all_scenarios$sigma.psi.sp[s],
                        mu.psi.yr = all_scenarios$mu.psi.yr[s], 
                        sigma.psi.yr = all_scenarios$sigma.psi.yr[s],
                        ## visit
                        mu.v.0 = all_scenarios$mu.v.0[s], 
                        mu.v.yr = all_scenarios$mu.v.yr[s],
                        ## type sym
                        type.range = all_scenarios$type.range[s],
                        type.visit = all_scenarios$type.visit[s], 
                        prop.visits.same = all_scenarios$prop.visits.same[s])
  
  save(sim.data,
       file=paste0("multi_sp/p2.2/outputs/", "sim.data/",
                   
                   "_r_",all_scenarios[s,'r'],"_p.yr_",all_scenarios[s,'p.yr'],
                   "_mu.psi.yr_",all_scenarios[s,'mu.psi.yr'], "_nyr_",all_scenarios[s,'nyr'],
                   "_prop.visits.same_",all_scenarios[s,'prop.visits.same'],
                   "_range_", all_scenarios[s,'type.range'], "_visit_", all_scenarios[s,'type.visit'],
                   "_s_", s,
                   ".RData"))

    
}


lapply(1:nrow(all_scenarios), run_all_simulation_1, all_scenarios = 
         all_scenarios)


