library(tidyverse)
library(ggforce)

pd12 <- read_csv("HostMobility/data/prop_disp_12.csv", col_names = FALSE)
pd13 <- read_csv("HostMobility/data/prop_disp_13.csv", col_names = FALSE)
pd14 <- read_csv("HostMobility/data/prop_disp_14.csv", col_names = FALSE)
pd15 <- read_csv("HostMobility/data/prop_disp_15.csv", col_names = FALSE)
normal <- read_csv("HostMobility/data/prop_disp_normal.csv", col_names = FALSE)

coord_dist <- read_csv("ModelFitting/data/coord_distances_R3.csv")
coord_dist2 <- read_csv("HostMobility/data/coord_distances_meta12.csv")

pd12$pop2 <- paste0("X",1:74)
pd12 <- pd12 %>% pivot_longer(cols = paste0("X",1:74), names_to = 'pop1')

pd12 %>% filter(pop1 != pop2) %>% mutate(group1 = ifelse(pop1 %in% c(paste0("X",1:37)),1,2),
                                         group2 = ifelse(pop2 %in% c(paste0("X",1:37)),1,2)) %>% 
  filter(group1 == 1) %>% group_by(group2) %>% summarize(mean_value = sum(value)/37*100)

pd12 <- pd12 %>% mutate(pop1 = factor(pop1, levels = c(paste0("X",1:74))),
                        pop2 = factor(pop2, levels = c(paste0("X",1:74))))

pd12 %>% filter(pop1 == "X5") %>% ggplot() +
  aes(x = pop2, y = value) + geom_point() + theme_classic(base_size = 15)

pd13$pop2 <- paste0("X",1:74)
pd13 <- pd13 %>% pivot_longer(cols = paste0("X",1:74), names_to = 'pop1')

pd13 %>% filter(pop1 != pop2) %>% mutate(group1 = ifelse(pop1 %in% c(paste0("X",1:37)),1,3),
                                         group2 = ifelse(pop2 %in% c(paste0("X",1:37)),1,3)) %>% 
  filter(group1 == 1) %>% group_by(group2) %>% summarize(mean_value = sum(value)/37*100)

pd13 %>% filter(pop1 != pop2) %>% mutate(group1 = ifelse(pop1 %in% c(paste0("X",1:37)),1,3),
                                         group2 = ifelse(pop2 %in% c(paste0("X",1:37)),1,3)) %>% 
  filter(group1 == 1, group2 == 3) %>% summarize(max_value = max(value)*100)

pd13 <- pd13 %>% mutate(pop1 = factor(pop1, levels = c(paste0("X",1:74))),
                        pop2 = factor(pop2, levels = c(paste0("X",1:74))))

pd13 %>% filter(pop1 == "X5") %>% ggplot() +
  aes(x = pop2, y = value) + geom_point() + theme_classic(base_size = 15)


pd14$pop2 <- paste0("X",1:74)
pd14 <- pd14 %>% pivot_longer(cols = paste0("X",1:74), names_to = 'pop1')

pd14 %>% filter(pop1 != pop2) %>% mutate(group1 = ifelse(pop1 %in% c(paste0("X",1:37)),1,4),
                                         group2 = ifelse(pop2 %in% c(paste0("X",1:37)),1,4)) %>% 
  filter(group1 == 1) %>% group_by(group2) %>% summarize(mean_value = sum(value)/37*100)

pd14 <- pd14 %>% mutate(pop1 = factor(pop1, levels = c(paste0("X",1:74))),
                        pop2 = factor(pop2, levels = c(paste0("X",1:74))))

pd14 %>% filter(pop1 == "X5") %>% ggplot() +
  aes(x = pop2, y = value) + geom_point() + theme_classic(base_size = 15)


pd15$pop2 <- paste0("X",1:74)
pd15 <- pd15 %>% pivot_longer(cols = paste0("X",1:74), names_to = 'pop1')

pd15 %>% filter(pop1 != pop2) %>% mutate(group1 = ifelse(pop1 %in% c(paste0("X",1:37)),1,4),
                                         group2 = ifelse(pop2 %in% c(paste0("X",1:37)),1,4)) %>% 
  filter(group1 == 1) %>% group_by(group2) %>% summarize(mean_value = sum(value)/37*100)

pd15 <- pd15 %>% mutate(pop1 = factor(pop1, levels = c(paste0("X",1:74))),
                        pop2 = factor(pop2, levels = c(paste0("X",1:74))))

pd15 %>% filter(pop1 == "X5") %>% ggplot() +
  aes(x = pop2, y = value) + geom_point() + theme_classic(base_size = 15)


dist12 <- pd12 %>% filter(pop1 != pop2) %>% mutate(group1 = ifelse(pop1 %in% c(paste0("X",1:37)),1,2),
                                                   group2 = ifelse(pop2 %in% c(paste0("X",1:37)),1,2)) %>% 
  filter(group1 == 1) %>% group_by(group2) %>% summarize(mean_value = sum(value)/37*100) %>% 
  mutate(grid_dist = 7)
dist13 <- pd13 %>% filter(pop1 != pop2) %>% mutate(group1 = ifelse(pop1 %in% c(paste0("X",1:37)),1,2),
                                                   group2 = ifelse(pop2 %in% c(paste0("X",1:37)),1,2)) %>% 
  filter(group1 == 1) %>% group_by(group2) %>% summarize(mean_value = sum(value)/37*100) %>% 
  mutate(grid_dist = 14)
dist14 <- pd14 %>% filter(pop1 != pop2) %>% mutate(group1 = ifelse(pop1 %in% c(paste0("X",1:37)),1,2),
                                                   group2 = ifelse(pop2 %in% c(paste0("X",1:37)),1,2)) %>% 
  filter(group1 == 1) %>% group_by(group2) %>% summarize(mean_value = sum(value)/37*100) %>% 
  mutate(grid_dist = 21)
dist15 <- pd15 %>% filter(pop1 != pop2) %>% mutate(group1 = ifelse(pop1 %in% c(paste0("X",1:37)),1,2),
                                                   group2 = ifelse(pop2 %in% c(paste0("X",1:37)),1,2)) %>% 
  filter(group1 == 1) %>% group_by(group2) %>% summarize(mean_value = sum(value)/37*100) %>% 
  mutate(grid_dist = 28)

distance_df <- rbind(dist12,dist13,dist14,dist15)
distance_df <- distance_df %>% filter(group2 == 2)

write_csv(distance_df,"HostMobility/data/prop_dispersing.csv")

normal$pop2 <- paste0("X",1:37)
normal <- normal %>% pivot_longer(cols = paste0("X",1:37), names_to = 'pop1')
normal %>% group_by(pop2) %>% mutate(total = sum(value))

coord_dist <- coord_dist %>% mutate(pop2 = paste0("X",num), pop1 = paste0("X",ref))

coord_dist <- merge(coord_dist,normal)

coord_dist %>% filter(distance != 0) %>%
  group_by(pop2) %>% summarize(sum_n = sum(value))

coord_dist <- coord_dist %>% mutate(distance2 = distance/(sqrt(3/4)*2))

avg_dist <- coord_dist %>% filter(distance != 0) %>% group_by(pop2) %>% 
  summarize(mean_d = sum(distance2*value)) 

distm <- mean(avg_dist$mean_d)

pdf("spatial_model/figures/pdisp2.pdf",height = 5, width = 6)
coord_dist %>% filter(distance != 0) %>% ggplot() + aes(x = distance2, y = value) + geom_point() + theme_classic() +
  xlab("Distance (Grid Units)") + ylab("Proportion dispersing from i to j") +
  annotate(geom = 'text', x = 2, y = 0.020, label = paste0("Avg. distance = ", round(mean(avg_dist$mean_d),2)), 
           color = 'blue') + 
  geom_point(data = coord_dist[coord_dist$num == 1 & coord_dist$distance != 0,], aes(x = distance2, y = value), color = 'red') 
dev.off()

50/distm

ref <- coord_dist %>% filter(ref == 1) %>% filter(num == 1)

pdf('HostMobility/figures/distances.pdf')
coord_dist %>% filter(ref == 1) %>% mutate(distance = distance/(sqrt(3/4)*2)) %>% 
  ggplot() + aes(x0 = xvals, y0 = yvals, sides = 6, angle = 0, r = 1) +
  geom_regon(color = 'grey55', fill = 'white') + theme_classic() + 
  geom_regon(data = ref, aes(x0 = xvals, y0 = yvals, sides = 6, angle = 0, r = 1), color = 'red', fill = 'white', size = 1.5) + 
  coord_fixed(xlim = c(-7,7), ylim = c(-7,7)) + 
  theme(legend.position = 'none') +
  geom_text(aes(x = xvals, y = yvals, label = round(distance,2)), color = 'black')
dev.off()
