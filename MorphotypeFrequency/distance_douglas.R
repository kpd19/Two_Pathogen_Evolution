library(tidyverse)
library(gridExtra)
library(geosphere)

morph <- read_csv("MorphotypeFrequency/data/2021_morphotype_dist.csv")
site_data <- read_csv("MorphotypeFrequency/data/forest_composition_wide.csv")
site_data <- site_data %>% rename(doug_fir = `Douglas-fir`)

dist_mat <- distm(cbind(site_data$longitude,site_data$latitude), fun = distGeo)

dist_df <- c()
for(i in 1:dim(dist_mat)[1]){
  temp <- data.frame(site_a = site_data$site_no[i],site_b = site_data$site_no,
                     dist = dist_mat[i,], pd_a = site_data$doug_fir[i], pd_b = site_data$doug_fir)
  dist_df <- rbind(dist_df,temp)
}

dist_df <- dist_df %>% filter(site_a != site_b)

pdf("MorphotypeFrequency/figures/dist_doug.pdf",height = 4, width = 6)
dist_df %>% ggplot() + aes(x = dist/1000,y = abs(pd_a - pd_b)) + geom_point() + theme_classic() +
  coord_cartesian(xlim = c(0,50)) +
  #scale_x_log10() +
  xlab("Distance (km)") + 
  ylab("Difference between % Douglas-fir")
dev.off()

dist_df %>% arrange(dist) %>% mutate(dist = dist/1000)

site_data %>% filter(site_no %in% c(1,100))
