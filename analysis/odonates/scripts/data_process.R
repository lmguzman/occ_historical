# This script ingests GBIF and iDigBio data to explore real-world
# scenarios of museum data collection, community collection, visit
# probabilities and multispecies occupancy modeling for Odonates
# in eastern North America.

# load libraries
library(tidyverse); library(taxotools); library(data.table)
library(sf); library(cowplot); library(maps); library(mapdata)
library(ggpubr); library(nimble); library(raster); library(lubridate)
library(reclin); library(sp); library(rgdal); library(geosphere)

# set the project crs to North American Equal Area Conic (Albers)
crs_1 <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"

# source prewritten functions from the simulation studies, modified slightly to deal with 
# real world data
source('prep_data_real.R')

# load a basemap
basemap <- st_read("ne_10m_land.shp")
basemap <- basemap %>% 
  st_crop(xmin=-90, xmax=-50,
          ymin=10, ymax=50) %>%
  st_transform(crs_1)

# read in GBIF and iDigBio data
gbif_dat <- fread("../data/gbif.csv", 
                  sep="\t", 
                  stringsAsFactors=FALSE, 
                  quote="") %>%
  dplyr::filter(coordinateUncertaintyInMeters < 100*1000 | is.na(coordinateUncertaintyInMeters))

occ_dat <- gbif_dat

# filter the occurrence data to a sensible time frame (1970s-2010s)
sp_list <- as.data.frame(table(occ_dat$species)) %>%
  dplyr::filter(Freq >= 100)

occ_dat <- occ_dat %>%
  filter(between(year, 1970, 2019), !is.na(species), species!="", species %in% sp_list$Var1) %>%
  dplyr::select(species, year, decimalLatitude, decimalLongitude)
head(occ_dat)
plot(occ_dat$decimalLongitude, occ_dat$decimalLatitude)


# approximate the number of "community visits" using event dates, collectors, and locations
comm_coll <- gbif_dat %>%
  dplyr::select(species, eventDate, recordedBy, decimalLatitude, decimalLongitude )%>% 
  distinct() %>% 
  filter(!is.na(eventDate), !is.na(decimalLatitude), !is.na(decimalLongitude)) %>% 
  filter(!is.na(recordedBy)) %>%
  mutate(date_clean = ymd(eventDate)) %>% 
  filter(date_clean > ymd("1970-01-01")) %>% 
  data.table()

comm_coll$recordedBy %>% table() %>% sort(decreasing = TRUE) %>% head()

# Cluster observations by same day within 1km square areas
comm_coll <- gbif_dat %>%
  dplyr::select(species, eventDate, recordedBy, decimalLatitude, decimalLongitude )%>% 
  distinct() %>% 
  filter(!is.na(eventDate), !is.na(decimalLatitude), !is.na(decimalLongitude)) %>% 
  filter(!is.na(recordedBy)) %>%
  mutate(date_clean = ymd(eventDate)) %>% 
  filter(date_clean > ymd("1970-01-01")) %>% 
  data.table()

n_obs_day <- table(comm_coll$date_clean)
single_obs_dates <- names(n_obs_day[n_obs_day == 1])
unique_dates <- names(n_obs_day[n_obs_day > 1])

cluster_lists <- list()
sin_obs_data <- comm_coll[date_clean %in% ymd(single_obs_dates)]
sin_obs_data[, "cluster" := paste0(date_clean, "-", 1)]
cluster_lists[[1]] <- sin_obs_data
counter <- 2

for(date_use in unique_dates){
  
  cur_date <- comm_coll[date_clean == date_use]
  
  lat_lon <- cur_date[,.(decimalLongitude, decimalLatitude)]
  
  xy <- SpatialPointsDataFrame(
    lat_lon, data.frame(ID=seq(1:nrow(lat_lon))),
    proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
  
  # use the distm function to generate a geodesic distance matrix in meters
  mdist <- distm(xy)
  
  # cluster all points using a hierarchical clustering approach
  hc <- hclust(as.dist(mdist), method="complete")
  
  # define the distance threshold
  d=1000
  
  cur_date[, "cluster" := paste0(date_clean, "-", cutree(hc, h=d))]
  
  cluster_lists[[counter]] <- cur_date
  
  counter <- counter+1
}

all_clusters <- rbindlist(cluster_lists)

size_clusters <- table(all_clusters$cluster) %>% table()

# percentage of community samplings
size_clusters[1]/sum(size_clusters[1:length(size_clusters)])

###############################################################################
# convert the occurrence data to a spatial object
occ_spat <- st_as_sf(occ_dat, 
                     coords=c("decimalLongitude", "decimalLatitude"),
                     crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
occ_spat <- occ_spat %>% st_transform(crs_1)
head(occ_spat)

grid_1 <- st_make_grid(extent(basemap), cellsize=100*1000) %>%
  st_as_sf(crs=crs_1) %>%
  dplyr::mutate(GID=row_number())

# format occurrence data for data.prep function
# there are two cases here, all and detection. All includes
# all species at all sites, detection sites where one other species
# was detected. Range maps are inferred from convex polygons that
# are derived from occurrence data

# make range maps from convex polygons
range <- list()
sp_list <- unique(occ_spat$species)
for(i in 1:length(sp_list)){
  range[[i]] <- st_convex_hull(st_union(filter(occ_spat,
                                               species==sp_list[i])))
}

grid1_plot <- ggplot()+
  geom_polygon(data=basemap, mapping=aes(x=long, y=lat))+
  geom_sf(grid_1, mapping=aes(), fill=NA, color="grey")+
  geom_sf(range[[1]], mapping=aes(), fill="cyan", color="cyan", alpha=0.35)+
  geom_sf(filter(occ_spat, species==sp_list[1]), mapping=aes(), color="cyan")+
  geom_sf(range[[2]], mapping=aes(), fill="red", color="red", alpha=0.35)+
  geom_sf(filter(occ_spat, species==sp_list[2]), mapping=aes(), color="red")+
  geom_sf(range[[10]], mapping=aes(), fill="gold", color="gold", alpha=0.35)+
  geom_sf(filter(occ_spat, species==sp_list[10]), mapping=aes(), color="gold")+
  scale_x_continuous(limits=c(-93, -60))+
  theme_map()

# intersect range maps with the grid
grid_int <- list()
for(i in 1:length(range)){
  grid_int[[i]] <- matrix(as.matrix(st_intersects(grid_1, range[[i]])), nrow=24, ncol=26) 
}
res <- array(unlist(grid_int), dim=c(length(unique(occ_dat$species)), nrow(grid_1)))

# assign occurrences to the grid, create era codes for each
# occurrence, then take only unique era/visit/grid/species occurrences
sp_list <- sp_list %>% as.data.frame()
colnames(sp_list) <- c("species")

sp_list <- sp_list %>%
  dplyr::mutate(SPID=row_number())

occ_grid <- occ_spat %>%
  st_intersection(grid_1) %>%
  dplyr::mutate(era=(year-year%%10)) %>%
  dplyr::mutate(year=(year-era)+1) %>%
  dplyr::mutate(era=(era-1960)/10) %>%
  left_join(sp_list) %>%
  dplyr::select(SPID, era, year, GID) %>%
  unique()
head(occ_grid)

# merge clusters to occ_grid for community sampling investigation
com_clusters <- all_clusters %>%
  dplyr::mutate(year=year(date_clean)) %>%
  dplyr::mutate(era=(year-year%%10)) %>%
  dplyr::mutate(year=(year-era)+1) %>%
  dplyr::mutate(era=(era-1960)/10) %>%
  left_join(sp_list) %>%
  filter(!is.na(SPID))

com_clusters <- st_as_sf(com_clusters, coords=c('decimalLongitude', 'decimalLatitude'),
                         crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
com_clusters <- com_clusters %>% st_transform(crs=crs_1)

com_clusters <- grid_1 %>%
  st_intersection(com_clusters) %>%
  dplyr::select(SPID, era, year, GID, cluster) %>%
  unique()

com_clusters_count <- com_clusters %>%
  dplyr::select(GID, era, year, cluster) %>%
  group_by(GID, era, year) %>%
  dplyr::mutate(n=n()) %>%
  ungroup() %>%
  arrange(n)

com_clusters_count <- com_clusters_count %>%
  inner_join(data.frame(table(all_clusters$cluster)), by=c("cluster"="Var1"))

com_filter <- com_clusters_count %>%
  dplyr::group_by(GID, era, year) %>%
  dplyr::mutate(com_vis=length(Freq[Freq>=2])) %>% ungroup() %>%
  dplyr::mutate(prop=com_vis/n) %>%
  dplyr::filter(prop >= 0.5) %>%
  dplyr::select(era, year, GID, prop) %>%
  dplyr::mutate(key=paste0(era, year, GID))

com_cluster_test <- com_clusters_count %>%
  dplyr::select(GID, era) %>%
  unique()
era_grid <- table(com_cluster_test$GID, com_cluster_test$era) %>% rowSums()
sum(era_grid==1)
era_grid <- era_grid %>%
  as.data.frame() %>%
  dplyr::filter(. > 1)

com_cluster_visit <- com_clusters_count %>%
  dplyr::select(GID, era, year) %>%
  unique()

pre_filter <- com_cluster_visit %>%
  st_drop_geometry() %>%
  group_by(era, GID) %>%
  dplyr::mutate(n=n()) %>%
  ungroup() %>%
  dplyr::select(era, n) %>%
  unique() %>%
  group_by(era) %>%
  dplyr::mutate(ave=mean(n)) %>%
  dplyr::select(era, ave) %>%
  unique()
mean(pre_filter$ave)

com_filter2 <- dplyr::select(com_filter, era, year, GID, prop) %>% unique()

ggplot()+
  geom_histogram(com_filter2, mapping=aes(x=prop), binwidth=0.05)+
  scale_x_continuous(breaks=seq(0.5,1,by=0.05))+
  labs(x="Proportion of Community Visits",
       y="Frequency")+
  theme_half_open()+
  background_grid()

occ_grid <- occ_grid %>%
  dplyr::mutate(key=paste0(era, year, GID))
occ_grid_count <- nrow(occ_grid)

occ_grid <- occ_grid %>%
  dplyr::filter(key %in% com_filter$key)

nrow(occ_grid)/occ_grid_count

size_clusters <- table(all_clusters$cluster) %>% table()

# percentage of single samplings
size_clusters[1]/nrow(all_clusters)

# create the master matrix of occurrences
X <- array(0, dim=c(nsp=nrow(sp_list),
                    nsite=nrow(grid_1),
                    nyr=5,
                    nvisit=10))

for(i in 1:nrow(occ_grid)){
  X[occ_grid$SPID[i],
    occ_grid$GID[i], 
    occ_grid$era[i], 
    occ_grid$year[i]] <- 1
}
Z <- apply(X, c(1,2,3), sum)
dd <- list(Z=Z, X=X, nsite=nrow(grid_1), nsp=nrow(sp_list), nyr=5, nvisit=10,
           vis.arr=array(1, dim=c(nsp=nrow(sp_list),
                                  nsite=nrow(grid_1),
                                  nyr=5,
                                  nvisit=10)),
           sp.range=res)

# all_range model
dd_all_range <- dd
dd_all_range_prep <- prep.data(dd_all_range, limit.to.visits="all", 
                               limit.to.range="yes")

# detected_all modeling case
dd_det_all <- dd
dd_det_all_prep <- prep.data(dd_det_all, limit.to.visits="detected", 
                             limit.to.range="no")

# detected_range modeling case
dd_det_range <- dd
dd_det_range_prep <- prep.data(dd_det_range, limit.to.visits="detected", 
                               limit.to.range="yes")


# model code
ms_nimble <- nimbleCode({
  
  ### priors
  
  #### detection #####
  
  mu.p.0   ~ dnorm(mean = 0, sd= 10)
  
  ## effect of yr on detection
  p.yr ~ dnorm(mean = 0, sd = 10)
  
  ## random effect of species on detection
  sigma.p.sp   ~ dunif(min = 0, max = 10)
  
  for(sp in 1:nsp) {
    p.sp[sp]   ~ dnorm(mean = 0, sd = sigma.p.sp)
  }
  # 
  ## random effect of site on detection
  sigma.p.site ~ dunif(min = 0, max = 10)
  
  for(site in 1:nsite) {
    for(yr in 1:nyr) {
      p.site[site,yr]   ~ dnorm(0, sd = sigma.p.site)
    }
  }
  
  #### occupancy. 
  
  ## intercepts
  
  mu.psi.0 ~ dnorm(mean = 0, sd = 10)
  
  # species specific intercepts
  
  sigma.psi.sp ~ dunif(0,10)
  
  for(sp in 1:nsp){
    psi.sp[sp] ~ dnorm(0,sd = sigma.psi.sp)
  }
  
  ## species specific slopes 
  
  mu.psi.yr ~ dnorm(0, sd = 10)
  
  sigma.psi.yr ~ dunif(0,10)
  
  for(sp in 1:nsp){
    psi.yr[sp] ~ dnorm(mu.psi.yr, sd = sigma.psi.yr)
  }
  
  
  
  for(sp in 1:nsp) {
    for(yr in 1:nyr) {
      ## occupancy
      logit(psi[yr,sp]) <-
        mu.psi.0 +
        psi.yr[sp]*(yr) +
        psi.sp[sp]
      ## detection
      for(site in 1:nsite) {
        logit(p[site,yr,sp]) <-
          mu.p.0 +
          p.yr*(yr) +
          p.sp[sp] +
          p.site[site,yr]
      }
    }
  }
  
  ## latent state and likelihood
  for(sp in 1:nsp){
    for(site in 1:nsite){
      for(yr in 1:nyr) {
        
        #latent state
        Z[site,yr,sp] ~ dbern(psi[yr,sp])
        
      }
      
    }
    
  }
  
  for(ind in 1:nind) {
    mu.p[ind] <-
      Z[sitev[ind],yrv[ind],spv[ind]] *
      p[sitev[ind],yrv[ind],spv[ind]]
    X[ind] ~ dbern(mu.p[ind])
  }
  
})
params <- c('mu.p.0',
            'p.yr',
            "p.sp",
            "p.site",
            "psi",
            "psi.yr",
            "psi.sp",
            'sigma.p.sp',
            'sigma.p.site',
            'mu.psi.0',
            'sigma.psi.sp',
            'mu.psi.yr',
            'sigma.psi.yr')

# save workspace to port to ComputeCanada
save.image("ODE_env_100x100km_5occInt_com.RData")

