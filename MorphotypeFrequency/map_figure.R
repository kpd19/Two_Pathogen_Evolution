library(tidyverse)
library(raster)
library(scatterpie)
library(ggnewscale)
library(geosphere)
library(sf)
library(geodata)

# British Columbia and US State data:

latlong <- read_csv("data/morphotye_distribution_data.csv")
latlong$total2 <- latlong$SNPV + latlong$MNPV

canada_data = gadm(country="CAN", level = 1, path = tempdir())
bc_data = canada_data[canada_data$NAME_1 %in% c("British Columbia",'Alberta','Saskatchewan'), ]
bc_data =st_as_sf(bc_data)

us_data = gadm(country="USA", level = 1, path = tempdir())
state_data = us_data[us_data$NAME_1 %in% c("Washington", "Oregon", "Idaho", "California", "Nevada", "Arizona", "New Mexico", "Utah", "Colorado","Montana","Wyoming"), ]
state_data = st_as_sf(state_data)

st_crs(bc_data) <- NA
st_crs(state_data) <- NA

# Tree Species data:

douglas_fir = read_sf("tree_polygons/pseumenz/", "pseumenz")
douglas_fir <- douglas_fir %>% rename(sp1 = 'PSEUMENZ_',sp2 = 'PSEUMENZ_I' ) %>% mutate(species = "Douglas-fir",
                                                                                        genus = 'Douglas-fir')

bigcone_doug = read_sf("tree_polygons/pseumacr/", "pseumacr")
bigcone_doug <- bigcone_doug %>% rename(sp1 = 'PSEUMACR_',sp2 = 'PSEUMACR_I' ) %>% mutate(species = "Bigcone Douglas-fir",
                                                                                          genus = 'Douglas-fir') 

grand_fir = read_sf("tree_polygons/abiegran/", "abiegran")
grand_fir <- grand_fir %>% rename(sp1 = 'ABIEGRAN_',sp2 = 'ABIEGRAN_I' ) %>% mutate(species = "Grand fir",
                                                                                    genus = 'Abies spp.')

white_fir = read_sf("tree_polygons/abieconc/", "abieconc")
white_fir <- white_fir %>% rename(sp1 = 'ABIECONC_',sp2 = 'ABIECONC_I' ) %>% mutate(species = "White fir",
                                                                                    genus = 'Abies spp.') 

subalpine_fir = read_sf("tree_polygons/abielasi/", "abielasi")
subalpine_fir <- subalpine_fir %>% rename(sp1 = 'ABIELASI_',sp2 = 'ABIELASI_I' ) %>% mutate(species = "Subalpine fir",
                                                                                            genus = 'Abies spp.') 

pacific_silver = read_sf("tree_polygons/abieamab/", "abieamab")
pacific_silver <- pacific_silver %>% rename(sp1 = 'ABIEAMAB_',sp2 = 'ABIEAMAB_I' ) %>% mutate(species = "Pacific silver fir",
                                                                                              genus = 'Abies spp.') 

red_fir = read_sf("tree_polygons/abiemagn/", "abiemagn")
red_fir <- red_fir %>% rename(sp1 = 'ABIEMAGN_',sp2 = 'ABIEMAGN_I' ) %>% mutate(species = "Red fir",
                                                                                genus = 'Abies spp.') 

firs_int <- rbind(douglas_fir,grand_fir, white_fir,subalpine_fir,pacific_silver,
                  red_fir)
firs_int <- firs_int %>% filter(CODE == 1)
#firs_int <- firs_int %>% sf::st_set_crs(4326)

firs_sf <- dplyr::bind_rows(list(grand_fir, white_fir,subalpine_fir,pacific_silver,
                                 red_fir))

all_firs <- st_cast(st_union(firs_sf))

clust <- unlist(st_intersects(firs_sf, all_firs))

diss <- cbind(firs_sf, clust) %>% filter(CODE == 1) %>% 
  group_by(clust) %>%
  summarize() %>% mutate(genus = "Abies spp.")
#dissolve_firs <- st_combine(firs_sf)

douglas_fir2 <- douglas_fir %>% filter(CODE == 1) %>% dplyr::select(geometry,genus) %>% mutate(clust = 0)

tree_data = rbind(douglas_fir2, diss)

#tree_data <- tree_data %>% sf::st_set_crs(4326)
tree_data$genus <- factor(tree_data$genus,levels = c('Abies spp.','Douglas-fir'))

pdf("figures/tree_distributions.pdf",height = 10, width = 15)
#mapcols = c("#134F5C", "#B1B27A")
ggplot() + geom_sf(data = firs_int, aes(fill = species, color = species), alpha = 0.75) +
  # scale_fill_manual(name="", values = mapcols) +
  # scale_color_manual(name="", values = mapcols) +
  # scale_alpha_manual(name = "",values = c("Abies spp." = 0.8,'Douglas-fir' = 1)) +
  geom_sf(data = state_data, aes(geometry = geometry), color = "black", fill = NA) + 
  geom_sf(data = bc_data, aes(geometry = geometry), color = "black", fill = NA) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = 'none') + 
  facet_wrap(~species) +
  coord_sf(ylim = c(32.25,55),xlim = c(-125.5,-103.0)) 
dev.off()

llround <- latlong %>% mutate(lat_round = round(latitude/1)*1,
                              lon_round = round(longitude/1)*1)

xy <- SpatialPointsDataFrame(
  cbind(llround[3],llround[2]), data.frame(ID=seq(1:dim(llround)[1])),
  proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))

# use the distm function to generate a geodesic distance matrix in meters
mdist <- distm(xy)

# cluster all points using a hierarchical clustering approach
hc <- hclust(as.dist(mdist), method="complete")

d = 100000
xy$clust <- cutree(hc, h=d)

cluster_centroids <- st_as_sf(xy) %>% group_by(clust) %>% summarize(geometry = st_union(geometry)) %>% st_centroid()

cc <- data.frame(st_coordinates(cluster_centroids), clust = cluster_centroids$clust)
colnames(cc) <- c('lon_cent','lat_cent','clust')

llround$clust <- xy$clust

llround <- merge(llround,cc)

llround2 <- llround %>% group_by(clust,lat_cent,lon_cent) %>%
  summarize(MNPV = sum(MNPV),
            SNPV = sum(SNPV),
            total = sum(total)) %>% ungroup() %>%  mutate(max_total = max(total))

llround2 <- llround2 %>% mutate(radius = total/(max_total*3) + 0.2)

llround2 <- llround2 %>% rename(`Single-capsid morphotype` = "SNPV", `Multi-capsid morphotype` = "MNPV")

llround2 <- llround2 %>% mutate(nudge = case_when(total <= 2 ~ 0.25,
                                                  total >2 & total <= 10 ~ 0.325,
                                                  total >10 & total <= 50 ~ 0.5,
                                                  total >50 & total <= 160 ~ 0.6,
                                                  total >60 & total <= 400 ~ 0.75,
                                                  total >400 ~ 0.85))

get_x <- function(y){
  log(y)/20+ 0.10
}

llround2 <- llround2 %>% mutate(radius = get_x(total))

llround2 <- llround2 %>% mutate(nudge = case_when(total <= 2 ~ 0.275,
                                                  total >2 & total <= 10 ~ 0.35,
                                                  total >10 & total <= 50 ~ 0.5,
                                                  total >50 & total <= 160 ~ 0.6,
                                                  total >60 & total <= 400 ~ 0.75,
                                                  total >400 ~ 0.85))

mapcols = c("#275D3C", "#CAB462")

plt1 <- ggplot() + geom_sf(data = tree_data, aes(fill = genus, color = genus,alpha = genus)) +
  scale_fill_manual(name="", values = mapcols,labels = c(expression(italic("Abies spp.")), 'Douglas-fir')) +
  scale_color_manual(name="", values = mapcols,labels = c(expression(italic("Abies spp.")), 'Douglas-fir')) +
  scale_alpha_manual(name = "",values = c("Abies spp." = 0.7,'Douglas-fir' = 0.8),
                     labels = c(expression(italic("Abies spp.")), 'Douglas-fir')) +
  geom_sf(data = state_data, aes(geometry = geometry), color = "black", fill = NA) + 
  geom_sf(data = bc_data, aes(geometry = geometry), color = "black", fill = NA) +
  theme_classic(base_size = 20) +
  theme(#axis.title = element_blank(),
  #      axis.text = element_blank(),
        #axis.ticks = element_blank(),
        #axis.line = element_blank(),
        legend.key = element_rect(color = 'transparent', fill = 'transparent'))  +
  labs(y = "Latitude", x = "Longitude")  +
  coord_sf(ylim = c(34,51),xlim = c(-124,-104.9)) 


pdf("figures/morph_dist_map.pdf",height = 12,width = 16)
plt1 + new_scale_fill() + 
  geom_scatterpie(data = llround2, aes(x = lon_cent, y = lat_cent, r = radius),
                  cols = c('Single-capsid morphotype','Multi-capsid morphotype'), color = 'white', alpha = 1) + 
  scale_fill_manual("", values = c("Single-capsid morphotype" = "#ee8800", "Multi-capsid morphotype" = '#5D65C5')) + 
  geom_scatterpie_legend(breaks = c(get_x(1), get_x(10), get_x(100),get_x(550)), 
                         x = -123, y = 34, labeller = function(x) exp((x - 0.1)*20), label_position = 'right') + 
  geom_text(data = llround2, aes(x = lon_cent + nudge, y = lat_cent, label = total), color = 'grey5', size = 5)
dev.off()

