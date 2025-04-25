library(tidyverse)
library(ggforce)

grid.r <-6

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

coord_dist %>% filter(ref == 1) %>%
  ggplot() + aes(x0 = xvals, y0 = yvals, sides = 6, angle = 0, r = 1, fill = distance) +
  geom_regon() + theme_classic() + 
  coord_fixed() + 
  scale_fill_viridis_c(option = 'turbo') + 
  geom_text(aes(x = xvals,y = yvals, label = num), color = 'white')


write.csv(coord_dist, paste0("_data/coord_distances_R",grid.r,".csv"))

