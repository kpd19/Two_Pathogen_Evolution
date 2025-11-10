library(tidyverse)

ll_data <- read_csv("CompareModels/data/morphotype_dist_data.csv")

ll_data <- ll_data %>% mutate(n_trees_jit = jitter(n_trees,0.1),
                              pmnpv = MNPV/total)

ll_data <- ll_data %>% mutate(n = 1) %>%
  group_by(n_trees,pmnpv) %>% mutate(sum_n = sum(n))
ll_data <- ll_data %>% mutate(n_trees_jit = ifelse(sum_n == 1,n_trees,n_trees_jit))
ll_data <- ll_data %>% group_by(n_trees,pmnpv) %>%
  mutate(jit2 = ifelse(sum_n >1, seq(-0.1,0.1, length.out = sum_n), 0))

pdf("MorphotypeFrequency/figures/conceptual_mnpv2.pdf",height = 2.5, width = 15)
ll_data %>% filter(n_trees %in% c(1,2,3,4,5,36)) %>%
  ggplot() + aes(x =n_trees_jit, y = MNPV/total*100) +
  geom_point(size = 3, alpha = 0.75) +
  theme_bw(base_size = 22) + 
  facet_wrap(~n_trees, scales = 'free_x',nrow = 1) +
  scale_color_viridis_c(option = 'turbo') +
  theme(strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ylab("% Multi-capsid") +
  scale_x_continuous(breaks = c(1,2,3,4,5,36))
dev.off()

