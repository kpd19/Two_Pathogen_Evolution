library('tidyverse')
library('gridExtra')
library('zoo')

dose_2020_wide <- read.csv("Dose/data/dose_response_2020_wide.csv")
dose_2021_wide <- read_csv('Dose/data/dose_response_2021_wide.csv')

starting_values <- read.csv("Dose/data/starting_values.csv")

dose_2020_wide <- merge(dose_2020_wide,starting_values)

dose_2020_wide %>% filter(morphotype == "MNPV") %>%
  pivot_longer(cols = c('second','virus'), names_to = 'category', values_to = 'number') %>%
  ggplot() +
  aes(x= day, y = number, color = factor(category)) + geom_line(size= 2) + geom_point(size = 3) + 
  facet_wrap(~dose, scales = 'free_y') + 
  theme_classic(base_size = 15) + 
  scale_color_brewer(palette = "Set1", direction = -1)

dose_response %>% filter(morphotype == "SNPV") %>%
  pivot_longer(cols = c('second','virus'), names_to = 'category', values_to = 'number') %>%
  ggplot() +
  aes(x= day, y = number, color = factor(category)) + geom_line(size= 2) + geom_point(size = 3) + 
  facet_wrap(~dose, scales = 'free_y') + 
  theme_classic(base_size = 15) + 
  scale_color_brewer(palette = "Set1", direction = -1)

pdf("Dose/figures/dose_response_2020.pdf", height = 7, width = 15)
dose_response %>% drop_na(alive) %>% arrange(day) %>% group_by(morphotype,strain,dose) %>% mutate(csum_early = cumsum(dead_early),
                                                                       csum_second = cumsum(second),
                                                                       csum_virus = cumsum(virus),
                                                                       csum_dead = cumsum(dead)) %>% 
  select(morphotype,strain,dose,day,alive,csum_early,csum_second,csum_virus,csum_dead,start) %>% 
  pivot_longer(cols = c('alive','csum_early','csum_second','csum_virus','csum_dead'),
               names_to = 'category', values_to = 'number') %>%
  mutate(cat2 = case_when(category == 'csum_dead' ~ 'non-virus dead',
                          category == 'csum_early' ~ 'died too early',
                          category == 'csum_second' ~ 'second instar',
                          category == 'csum_virus' ~ 'virus killed',
                          category == 'alive' ~ 'first instar')) %>% 
  mutate(cat2 = factor(cat2, levels = c('first instar', 'second instar', 'virus killed','non-virus dead', 'died too early'))) %>% 
  ggplot() +
  aes(x= day, y = number/start, fill = factor(cat2)) + geom_bar(stat = 'identity') +  
  facet_grid(morphotype~dose, scales = 'free_y') + 
  theme_classic(base_size = 15) + 
  scale_fill_manual("", values = c('non-virus dead' = '#6a3d9a','died too early' = 'grey65',
                                   'second instar' = '#B2DF8A', 'virus killed' = '#ff7f00','first instar' = '#22801C')) +
  xlab("Day Post Infection") + ylab("Proportion") + 
  scale_x_continuous(breaks = c(7,9,11,13,15))
dev.off()


pdf("Dose/figures/dose_response_2021.pdf", height = 6, width = 18)
dose_2021_wide %>% 
  group_by(dose,day,morphotype) %>% summarize(alive = sum(alive),
                                              second = sum(second),
                                              virus = sum(virus),
                                              dead = sum(dead),
                                              dead_early = sum(dead_early),
                                              start = sum(start)) %>%
  arrange(day) %>% group_by(morphotype,dose) %>% mutate(csum_early = cumsum(dead_early),
                                                                            csum_second = cumsum(second),
                                                                            csum_virus = cumsum(virus),
                                                                            csum_dead = cumsum(dead)) %>% 
  select(morphotype,dose,day,alive,csum_early,csum_second,csum_virus,csum_dead,start) %>% 
  pivot_longer(cols = c('alive','csum_early','csum_second','csum_virus','csum_dead'),
               names_to = 'category', values_to = 'number') %>%
  mutate(cat2 = case_when(category == 'csum_dead' ~ 'non-virus dead',
                          category == 'csum_early' ~ 'died too early',
                          category == 'csum_second' ~ 'second instar',
                          category == 'csum_virus' ~ 'virus killed',
                          category == 'alive' ~ 'first instar')) %>% 
  mutate(cat2 = factor(cat2, levels = c('first instar', 'second instar', 'virus killed','non-virus dead', 'died too early'))) %>% 
  ggplot() +
  aes(x= day, y = number/start, fill = factor(cat2)) + geom_bar(stat = 'identity') +  
  facet_grid(morphotype~dose, scales = 'free_y') + 
  theme_classic(base_size = 15) + 
  scale_fill_manual("", values = c('non-virus dead' = '#6a3d9a','died too early' = 'grey65',
                                   'second instar' = '#B2DF8A', 'virus killed' = '#ff7f00','first instar' = '#22801C')) +
  xlab("Day Post Infection") + ylab("Proportion") + 
  scale_x_continuous(breaks = c(5,7,9,11,13,15))
dev.off()
