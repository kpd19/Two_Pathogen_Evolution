library(tidyverse)
library(ggforce)

`%ni%` <- Negate(`%in%`)

grid.r <-17

rows <- grid.r*2 + 1
cols <- rows

r = 1
c = sqrt(3/4)

n.cells = 3*(grid.r + 1)^2 - 3*(grid.r + 1) + 1

row_id <- seq(1,rows)
coord_0 <- seq(-grid.r,grid.r)
nums <- c(seq(grid.r + 1,rows),seq(rows -1, grid.r + 1))

coord_df <- c()
for(i in row_id){
  print(i)
  xvals <- rep(1.5*coord_0[i],nums[i])
  yvals <- seq(-(nums[i]-1)*c,(nums[i]-1)*c,2*c)
  
  temp <- data.frame(xvals,yvals)
  coord_df <- rbind(coord_df,temp)
}

coord_df$num <- seq(1,n.cells)

trees <- sample(x = c("DO","GR"), size = dim(coord_df)[1], prob = c(0.5,0.5), replace = TRUE)

coord_df$trees <- trees

pdf("/Users/katherinedixon/Documents/StuffINeed/_Research/_DFTM_2021/Field_2021/_plots/_eeid/hexagons_linkedin2.pdf",height = 5, width = 5)
coord_df %>% mutate(trees = ifelse(trees == "DO",'Douglas fir','Grand fir')) %>%
  ggplot() + aes(x0 = xvals, y0 = yvals, sides = 6, angle = 0, r = 1, fill = trees) +
  geom_regon(color = 'palegreen4') + theme_classic() + 
  coord_fixed() + 
  # geom_text(aes(x = xvals,y = yvals, label = num), color = 'white') + 
  scale_fill_manual("", values = c("Douglas fir" = 'palegreen4', 'Grand fir' = 'darkgreen')) + 
  theme_void() +
  #theme(legend.position = c(0.9,0.93))
  theme(legend.position = 'none')
dev.off()

num_vals <- unique(coord_df$num)

coord_dist <- c()
for(i in 1:length(num_vals)){
  point.x <- coord_df %>% filter(num == num_vals[i]) %>% pull(xvals)
  point.y <- coord_df %>% filter(num == num_vals[i]) %>% pull(yvals)
  
  temp <- coord_df %>% mutate(distance = sqrt(abs(xvals - point.x)^2 + abs(yvals - point.y)^2))
  temp$ref <- num_vals[i]
  
  coord_dist <- rbind(coord_dist,temp)
  
}

get_ids <- function(l1){
  ids <- c(l1[1]:(l1[1]+3),l1[2]:(l1[2]+4),l1[3]:(l1[3]+5),l1[4]:(l1[4]+6),l1[5]:(l1[5]+5),l1[6]:(l1[6]+4),l1[7]:(l1[7]+3))
  return(ids)
}

id1 <- get_ids(c(1,19,38,58,80,103,127))
id2 <- get_ids(c(155,180,206,233,262,292,323))
id3 <- get_ids(c(358,390,423,457,492,526,559))
id4 <- get_ids(c(594,624,653,681,709,736,762))
id5 <- get_ids(c(790,813,835,856,877,897,916))


coord_dist <- coord_dist %>% mutate(group = case_when(num %in% id1 ~ 1,
                                                      num %in% id2 ~ 2,
                                                      num %in% id3 ~ 3,
                                                      num %in% id4 ~ 4,
                                                      num %in% id5 ~ 5,
                                                      num %ni% c(id1,id2,id3,id4,id5) ~ 0))

coord_dist %>% filter(ref == 1, group != 0) %>%
  ggplot() + aes(x0 = xvals, y0 = yvals, sides = 6, angle = 0, r = 1, fill = group) +
  geom_regon(color = 'grey55') + theme_classic() + 
  coord_fixed(xlim = c(-30,30), ylim = c(-30,30)) + 
  scale_fill_viridis_c(option = 'viridis') +
  theme(legend.position = 'none')+ 
  geom_text(aes(x = xvals,y = yvals, label = num), color = 'white')

pdf("AnalyzeSimulations/figures/grid_hm2.pdf",height = 10, width = 10)
coord_dist %>% filter(ref == 1, group %in% c(1,3)) %>%
  ggplot() + aes(x0 = xvals, y0 = yvals, sides = 6, angle = 0, r = 1, fill = as.factor(group)) +
  geom_regon(color = 'grey55') + theme_classic() + 
  coord_fixed(xlim = c(-30,30), ylim = c(-20,20)) + 
  scale_fill_viridis_c(option = 'viridis') +
  theme(legend.position = 'none') +
  scale_fill_manual(values = c(`1` = 'grey95', `2` = '#ffffcc', `3` = '#a1dab4', `4` = '#41b6c4', `5` ='#225ea8'))
  #geom_text(aes(x = xvals,y = yvals, label = group), color = 'white')
dev.off()

coord_dist %>% filter(num %in% c(id1,id3),ref %in% c(id1,id4)) %>% 
  filter(ref == 118, group !=0) %>%
  ggplot() + aes(x0 = xvals, y0 = yvals, sides = 6, angle = 0, r = 1, fill = distance) +
  geom_regon() + theme_classic() + 
  coord_fixed() + 
  scale_fill_viridis_c(option = 'viridis') + 
  geom_text(aes(x = xvals,y = yvals, label = num), color = 'white')

spat_12 <- coord_dist %>% filter(num %in% c(id1,id2),ref %in% c(id1,id2))
spat_13 <- coord_dist %>% filter(num %in% c(id1,id3),ref %in% c(id1,id3))
spat_14 <- coord_dist %>% filter(num %in% c(id1,id4),ref %in% c(id1,id4))
spat_15 <- coord_dist %>% filter(num %in% c(id1,id5),ref %in% c(id1,id5))

recode_12 <- data.frame(num = c(id1,id2),num2 = 1:74)
recode_13 <- data.frame(num = c(id1,id3),num2 = 1:74)
recode_14 <- data.frame(num = c(id1,id4),num2 = 1:74)
recode_15 <- data.frame(num = c(id1,id5),num2 = 1:74)

spat_12 <- merge(spat_12,recode_12)
spat_12 <- spat_12 %>% select(-num) %>% rename(num = num2)

spat_12 <- merge(spat_12, recode_12, by.x = c('ref'), by.y = c('num'))
spat_12 <- spat_12 %>% select(-ref) %>% rename(ref = num2)

spat_13 <- merge(spat_13,recode_13)
spat_13 <- spat_13 %>% select(-num) %>% rename(num = num2)

spat_13 <- merge(spat_13, recode_13, by.x = c('ref'), by.y = c('num'))
spat_13 <- spat_13 %>% select(-ref) %>% rename(ref = num2)

spat_14 <- merge(spat_14,recode_14)
spat_14 <- spat_14 %>% select(-num) %>% rename(num = num2)

spat_14 <- merge(spat_14, recode_14, by.x = c('ref'), by.y = c('num'))
spat_14 <- spat_14 %>% select(-ref) %>% rename(ref = num2)

spat_15 <- merge(spat_15,recode_15)
spat_15 <- spat_15 %>% select(-num) %>% rename(num = num2)

spat_15 <- merge(spat_15, recode_15, by.x = c('ref'), by.y = c('num'))
spat_15 <- spat_15 %>% select(-ref) %>% rename(ref = num2)

spat_15 %>% 
  filter(ref == 3, group !=0) %>%
  ggplot() + aes(x0 = xvals, y0 = yvals, sides = 6, angle = 0, r = 1, fill = distance) +
  geom_regon() + theme_classic() + 
  coord_fixed() + 
  scale_fill_viridis_c(option = 'viridis') + 
  geom_text(aes(x = xvals,y = yvals, label = num), color = 'white')

spat_12 <- spat_12 %>% select(xvals,yvals,num,distance,ref) %>% arrange(ref,num)
spat_13 <- spat_13 %>% select(xvals,yvals,num,distance,ref) %>% arrange(ref,num)
spat_14 <- spat_14 %>% select(xvals,yvals,num,distance,ref) %>% arrange(ref,num)
spat_15 <- spat_15 %>% select(xvals,yvals,num,distance,ref) %>% arrange(ref,num)

write.csv(spat_12, paste0("HostMobility/data/coord_distances_meta12.csv"))
write.csv(spat_13, paste0("HostMobility/data/coord_distances_meta13.csv"))
write.csv(spat_14, paste0("HostMobility/data/coord_distances_meta14.csv"))
write.csv(spat_15, paste0("HostMobility/data/coord_distances_meta15.csv"))

