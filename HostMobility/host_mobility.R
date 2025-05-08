library(tidyverse)
library(gridExtra)

df_means <- read_csv("HostMobility/data/host_mobility_means.csv")

df_means <- df_means %>% mutate(patch1 = ifelse(pd1 == 2, "5% Douglas-fir", "95% Douglas-fir"),
                             patch2 = ifelse(pd2 == 2, "5% Douglas-fir", "95% Douglas-fir")) %>% 
  mutate(type = paste0("Patch 1: ", patch1, ", Patch 2: ", patch2)) %>% 
  filter(type != "Patch 1: 95% Douglas-fir, Patch 2: 5% Douglas-fir") %>% 
  mutate(type2 = case_when(pd1 == 2 & pd2 == 2 ~ '5% Douglas-fir',
                           pd1 == 2 & pd2 == 35 ~ "Patch 1: 5% Douglas-fir, Patch 2: 95% Douglas-fir",
                           pd1 == 35 & pd2 == 35 ~ '95% Douglas-fir'))


#hm12 <- read_csv("data/host_mobility_patch12.csv")
#hm13 <- read_csv("data/host_mobility_patch13.csv")
#hm14 <- read_csv("data/host_mobility_patch14.csv")

df_means <- df_means %>% mutate(patch3 = ifelse(patch == 1, patch1,patch2))

df_means <- merge(df_means,distance_df)

df_means <- df_means %>% mutate(pdisp = factor(pdisp,levels = distance_df$pdisp))

plt1 <- df_means %>% filter(type2 == "5% Douglas-fir", grid_dist == 7) %>% 
  ggplot() + aes(x = patch3, y = nu1, group = interaction(patch3,pdisp)) +
  geom_boxplot(outliers = FALSE, fill = 'grey95') + theme_classic(base_size = 15) + 
  facet_wrap(~type2, drop = TRUE, scales = 'free_x') +
  scale_y_log10() + 
  xlab("") + ylab(expression(bar(nu)[SNPV])) +
  theme(legend.position = 'none') +
  coord_cartesian(ylim = c(5e-5,500))

plt2 <- df_means %>% filter(type2 ==  "Patch 1: 5% Douglas-fir, Patch 2: 95% Douglas-fir") %>% 
  ggplot() + aes(x = patch3, y = nu1, fill = as.factor(pdisp), group = interaction(patch3,pdisp)) +
  geom_boxplot(outliers = FALSE) + theme_classic(base_size = 15) + 
  facet_wrap(~type2, drop = TRUE, scales = 'free_x') +
  scale_y_log10() + 
  scale_fill_brewer("Proportion dispersing \nbetweeen patches", palette = "YlGnBu") + 
  xlab("Patch") + ylab('') +
  theme(legend.position = 'none', axis.text.y = element_blank())+
  coord_cartesian(ylim = c(5e-5,500)) 

plt3 <- df_means %>% filter(type2 == "95% Douglas-fir", grid_dist == 7) %>% 
  ggplot() + aes(x = patch3, y = nu1, group = interaction(patch3,pdisp)) +
  geom_boxplot(outliers = FALSE, fill = 'grey95') + theme_classic(base_size = 15) + 
  facet_wrap(~type2, drop = TRUE, scales = 'free_x') +
  scale_y_log10() + 
  xlab("") + ylab('') +
  theme(legend.position = 'none', axis.text.y = element_blank()) +
  coord_cartesian(ylim = c(5e-5,500))

plt4 <- df_means %>% filter(type2 == "5% Douglas-fir", grid_dist == 7) %>% 
  ggplot() + aes(x = patch3, y = nu2, group = interaction(patch3,pdisp)) +
  geom_boxplot(outliers = FALSE, fill = 'grey95') + theme_classic(base_size = 15) + 
  facet_wrap(~type2, drop = TRUE, scales = 'free_x') +
  scale_y_log10() + 
  xlab("") + ylab(expression(bar(nu)[MNPV])) +
  theme(legend.position = 'none') +
  coord_cartesian(ylim = c(5e-5,500))

plt5 <- df_means %>% filter(type2 ==  "Patch 1: 5% Douglas-fir, Patch 2: 95% Douglas-fir") %>% 
  ggplot() + aes(x = patch3, y = nu2, fill = as.factor(pdisp), group = interaction(patch3,pdisp)) +
  geom_boxplot(outliers = FALSE) + theme_classic(base_size = 15) + 
  facet_wrap(~type2, drop = TRUE, scales = 'free_x') +
  scale_y_log10() + 
  scale_fill_brewer("Proportion dispersing \nbetweeen patches", palette = "YlGnBu") + 
  xlab("Patch") + ylab('') +
  theme(legend.position = 'none', axis.text.y = element_blank())+
  coord_cartesian(ylim = c(5e-5,500))

plt6 <- df_means %>% filter(type2 == "95% Douglas-fir", grid_dist == 7) %>% 
  ggplot() + aes(x = patch3, y = nu2, group = interaction(patch3,pdisp)) +
  geom_boxplot(outliers = FALSE, fill = 'grey95') + theme_classic(base_size = 15) + 
  facet_wrap(~type2, drop = TRUE, scales = 'free_x') +
  scale_y_log10() + 
  xlab("") + ylab('') +
  theme(legend.position = 'none', axis.text.y = element_blank()) +
  coord_cartesian(ylim = c(5e-5,500))

pdf("AnalyzeSimulations/figures/host_mobility_transmission.pdf",height = 10, width = 12)
grid.arrange(plt1,plt2,plt3,plt4,plt5,plt6,nrow = 2, widths = c(0.6,2,0.5))
dev.off()

