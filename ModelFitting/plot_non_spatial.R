library(tidyverse)

simulation <- read_csv("/Users/katherinedixon/Documents/StuffINeed/_Research/Two_Pathogen_Evolution/spatial_model/realization_op/non_spatial_evolution/non_spatial_sig0.5.csv")

snpv_col = "#ee8800"
mnpv_col = '#5D65C5'

simulation <- simulation %>% mutate(tree_sp = ifelse(tree_sp == "DO",'Douglas-fir','Grand fir'))

pdf("figures/population_time_series.pdf",height = 6, width = 10)
simulation %>% filter(time >= 25, time <= 75, rep <= 5) %>% ggplot() +
  geom_line(aes(x = time, y = Host, color = "Host")) + 
  geom_line(aes(x = time, y = SNPV, color = "SNPV")) + 
  geom_line(aes(x = time, y = MNPV, color = "MNPV")) + 
  theme_classic(base_size = 15) + 
  scale_y_log10() + 
  scale_color_manual("", values = c('Host' = 'darkgreen','SNPV' = snpv_col,'MNPV' = mnpv_col)) +
  facet_grid(tree_sp~rep) + 
  xlab("Time (years)") + ylab(expression(log[10]~"Population size"))
dev.off()

pdf("figures/mean_infectiousness_time_series.pdf",height = 6, width = 10)
simulation %>% filter(time >= 25, time <= 75, rep <= 5) %>%
  ggplot() +
  geom_line(aes(x = time, y = SNPV, color = "SNPV")) + 
  geom_line(aes(x = time, y = MNPV, color = "MNPV")) + 
  theme_classic(base_size = 15) + 
  scale_color_manual(expression(bar(nu)), values = c('SNPV' = snpv_col,'MNPV' = mnpv_col)) +
  scale_y_log10() + 
  facet_grid(tree_sp~rep) + 
  xlab("Time (years)") + ylab(expression(log[10]~"Mean Infectiousness ("~bar(nu)~")"))
dev.off()

pdf("figures/fraction_infected_time_series.pdf",height = 6, width = 10)
simulation %>% filter(time >= 25, time <= 75, rep <= 5) %>% 
  pivot_longer(cols = c('frac_SNPV','frac_MNPV')) %>% 
  mutate(name = ifelse(name == 'frac_SNPV','SNPV','MNPV')) %>% 
  ggplot() + 
  geom_bar(aes(x = time, y = value, fill = name), stat = 'identity', position = 'stack') + 
  scale_fill_manual("", values = c('SNPV' = snpv_col, 'MNPV' = mnpv_col)) +
  theme_classic(base_size = 15) + 
  facet_grid(tree_sp~rep) + 
  xlab("Time (years)") + ylab("Fraction Infected")
dev.off()

  