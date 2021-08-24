###### Part 3 #####

## 1. 
# Simulation:
# - Species have spatially explicit ranges
# - Missing visits but same visit for all species 

# Modeling: 
# MS_all_all, MS_all_visits vs MS_all_detected for more than 2 eras 
# this would be useful for a scenario that has small spatial scale 

library(parallel)
library(dplyr)
source('dynamic/simulation/src/simulate_dynamic.R')
expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))


## Make a data frame of the parameters 

## Change either p.yr, mu.phi.yr or mu.gam.yr

base_scenario <- data.frame(nsp          = 50,
                             nsite        = 100,
                             nvisit       = 3,
                             ## detection
                             mu.p.0 = -0.5,
                             sigma.p.site = 1.5,
                             sigma.p.sp = 0.5,
                             ## persistence
                             mu.phi.0 = 0,
                             sigma.phi.sp = 0.5,
                             sigma.phi.yr = 0.2,
                             ## colonization
                             mu.gam.0 = 0,
                             sigma.gam.sp = 0.5,
                             sigma.gam.yr = 0.2,
                             ## visit
                             mu.v.0 = 0, 
                             mu.v.yr = -0.5,
                             ## type sym
                             type.range = "polys",
                             type.visit = 'visit_same_sp',
                             prop.visits.same = 1)

r <- data.frame(r = 1:5)
p.yr <- data.frame(p.yr = c(-1, -0.5, 0, 0.5, 1))
mu.phi.yr <- data.frame(mu.phi.yr = c(-1, -0.5, 0, 0.5, 1))
mu.gam.yr <- data.frame(mu.gam.yr = c(-1, -0.5, 0, 0.5, 1))
nyr <- data.frame(nyr = c(2,5,10))

all_scenarios <- bind_rows(expand.grid.df(base_scenario, p.yr, data.frame(mu.phi.yr = 0, mu.gam.yr = 0), r, nyr),
                           expand.grid.df(base_scenario, mu.phi.yr, data.frame(p.yr= 0, mu.gam.yr = 0), r, nyr),
                           expand.grid.df(base_scenario, mu.gam.yr, data.frame(mu.phi.yr = 0,  p.yr = 0), r, nyr))

run_all_simulation_1 <-function(s, all_scenarios){

  
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
                        ## persistence
                        mu.phi.0 = all_scenarios$mu.phi.0[s],
                        sigma.phi.sp = all_scenarios$sigma.phi.sp[s],
                        mu.phi.yr = all_scenarios$mu.phi.yr[s],
                        sigma.phi.yr = all_scenarios$sigma.phi.yr[s],
                        ## colonization
                        mu.gam.0 = all_scenarios$mu.gam.0[s],
                        sigma.gam.sp = all_scenarios$sigma.gam.sp[s],
                        mu.gam.yr = all_scenarios$mu.gam.yr[s], 
                        sigma.gam.yr = all_scenarios$sigma.gam.yr[s],
                        ## visit
                        mu.v.0 = all_scenarios$mu.v.0[s], 
                        mu.v.yr = all_scenarios$mu.v.yr[s],
                        ## type sym
                        type.range = all_scenarios$type.range[s],
                        prop.visits.same = all_scenarios$prop.visits.same[s])
  
  save(sim.data,
       file=paste0("p3/outputs/", "sim.data/",
                   
                   "_r_",all_scenarios[s,'r'],"_p.yr_",all_scenarios[s,'p.yr'],"_mu.phi.yr_",all_scenarios[s,'mu.phi.yr'],
                   "_mu.gam.yr_",all_scenarios[s,'mu.gam.yr'], "_nyr_",all_scenarios[s,'nyr'],"_range_", all_scenarios[s,'type.range'], "_visit_", all_scenarios[s,'prop.visits.same'],
                   ".RData"))

    
}


lapply(1:nrow(all_scenarios), run_all_simulation_1, all_scenarios = 
         all_scenarios)
