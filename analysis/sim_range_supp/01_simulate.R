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
source('analysis/simulation/src/simulate_ms.R')

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
                            ## type sym
                            type.range = 'polys',
                            type.visit = 'visit_miss')

r <- data.frame(r = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
p.yr <- data.frame(p.yr = c(-0.1, -0.05, 0, 0.05, 0.1))
mu.psi.yr <- data.frame(mu.psi.yr = c(-0.1, -0.05, 0, 0.05, 0.1))
mu.v.yr <- data.frame(mu.v.yr = c(-0.1, 0, 0.1))
nyr <- data.frame(nyr = c(10))
prop.visits.same <- data.frame(prop.visits.same = c(0, 0.25, 0.5, 0.75, 1))

all_scenarios <- bind_rows(expand.grid.df(base_scenario, 
                                          mu.v.yr, 
                                          mu.psi.yr,
                                          p.yr, 
                                          r, 
                                          nyr, 
                                          prop.visits.same)) %>% 
  unique()