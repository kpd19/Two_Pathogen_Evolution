library(raster)
library(tidyverse)
library(sp)

forest <- raster("data/conus_foresttype.img")
latlong <- read_csv("data/morphotye_distribution_data.csv")
tree_labels <- read_csv("data/tree_labels.csv")

sites <- latlong %>% select(latitude,longitude) %>% rename(Latitude = latitude, Longitude = longitude)
coordinates(sites)  <-  c("Longitude",  "Latitude")
proj4string(sites)  <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") # this is what google is
sites_transformed<-spTransform(sites, CRS("+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD27 +units=m +no_defs"))

sp_cover = raster::extract(forest,sites_transformed, buffer = 5000, factors = TRUE)

sp_cover_df <- c()
for(i in 1:length(latlong$site)){
  temp <- data.frame(site = latlong$site[i], number = sp_cover[[i]])
  sp_cover_df <- rbind(sp_cover_df,temp)
}

sp_cover_df <- merge(sp_cover_df,tree_labels,all_x = TRUE)

sp_cover_df <- sp_cover_df %>% drop_na(number)

forest_counts <- sp_cover_df %>% group_by(site) %>% count(name) %>% group_by(site) %>% mutate(sum_n = sum(n)) 

latlong2 <- merge(latlong,forest_counts,all = TRUE)

lat_long_wide <- latlong2 %>% drop_na(name) %>% mutate(prop = n/sum_n) %>% select(-n) %>% select(-sum_n) %>% pivot_wider(names_from = name, values_from = prop)

lat_long_wide <- lat_long_wide %>% mutate_at(colnames(lat_long_wide)[17:38], ~replace_na(.,0))

#write_csv(lat_long_wide,"forest_comp_usa.csv")

##################

`%ni%` <- Negate(`%in%`)

successful_sites <- unique(lat_long_wide$site_no)
all_sites <- unique(latlong$site_no)
unsuccessful_sites <- all_sites[all_sites %ni% successful_sites]

missing_sites <- latlong %>% filter(site_no %in% unsuccessful_sites)

usa_sites <- missing_sites %>% filter(state != "BC") 

sites2 <- usa_sites %>% select(latitude,longitude) %>% rename(Latitude = latitude, Longitude = longitude)
coordinates(sites2)  <-  c("Longitude",  "Latitude")
proj4string(sites2)  <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
sites_transformed2<-spTransform(sites2, CRS("+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD27 +units=m +no_defs"))

sp_cover2 = raster::extract(forest,sites_transformed2, buffer = 10000, factors = TRUE)

sp_cover_df2 <- c()
for(i in 1:length(usa_sites$site)){
  temp <- data.frame(site = usa_sites$site[i], number = sp_cover2[[i]])
  sp_cover_df2 <- rbind(sp_cover_df2,temp)
}

sp_cover_df2 <- merge(sp_cover_df2,tree_labels,all_x = TRUE)

sp_cover_df2 <- sp_cover_df2 %>% drop_na(number)

forest_counts2 <- sp_cover_df2 %>% group_by(site) %>% count(name) %>% group_by(site) %>% mutate(sum_n = sum(n)) 

f1 <- unique(forest_counts$site)
f2 <- unique(forest_counts2$site)

forest_counts4 <- rbind(forest_counts, forest_counts2)

latlong2 <- merge(latlong,forest_counts4,all = TRUE)

all_good <- latlong2[!is.na(latlong2$name),]
cr3 <- all_good %>% filter(site == "CRA3")
all_good <- all_good %>% filter(site != "CRA3")
problems <- latlong2[is.na(latlong2$name),]

# canadian sites not included in FIA dataset, they are predominantly Douglas-fir
problems[problems$state == "BC",]$name <- "Douglas-fir"
problems[problems$state == "BC",]$n <- 1
problems[problems$state == "BC",]$sum_n <- 1

# sites are Craters of the Moon National Monument in Idaho (CTR and CRA3) are Douglas-fir
problems[problems$site == "CTR",]$name <- "Douglas-fir"
problems[problems$site == "CTR",]$n <- 1
problems[problems$site == "CTR",]$sum_n <- 1

cr3 <- cr3 %>% mutate(name = "Douglas-fir", n = 1, sum_n = 1)

all_sites2 <- rbind(problems, all_good,cr3)

lat_long_wide <- all_sites2 %>% mutate(prop = n/sum_n) %>% select(-n) %>% select(-sum_n) %>% drop_na(name)
lat_long_wide <- lat_long_wide %>% pivot_wider(names_from = name, values_from = prop)

lat_long_wide <- lat_long_wide %>% mutate_at(colnames(lat_long_wide)[11:38], ~replace_na(.,0))

write_csv(lat_long_wide,"data/forest_composition_wide.csv")

df_vals <- unique(lat_long_wide$`Douglas-fir`)
df_options <- seq(0,1,1/37)

df_round <- c()
for(i in 1:length(df_vals)){
  df_round[i] <- df_options[which.min(abs(df_options - df_vals[i]))]
}

df_round_df <- data.frame(Douglas_fir = df_vals,df_round)

ll_small <- lat_long_wide %>% select(source,site,MNPV,total,`Douglas-fir`) %>%
  rename(Douglas_fir = `Douglas-fir`) %>%
  mutate(id = rownames(lat_long_wide)) 

ll_small <- merge(ll_small,df_round_df)

ll_small$n_trees <- ll_small$df_round*37

ll_small %>% arrange(n_trees,site)

ll_no_ext <- ll_small %>% mutate(n_trees = ifelse(n_trees == 0, 1,n_trees)) %>% 
  mutate(n_trees = ifelse(n_trees == 37, 36, n_trees))

write_csv(ll_no_ext,"data/data_for_ll.csv")

