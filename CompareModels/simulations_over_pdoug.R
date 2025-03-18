library(tidyverse)
library(gridExtra)
library(gg3D)

ll_data <- read_csv("data/data_for_ll.csv")
id1 <- ll_data %>% filter(n_trees %in% 1:36) %>%
  arrange(n_trees) %>% group_by(n_trees) %>%
  mutate(n = seq(1,length(total))) %>% filter(n == 1) %>% pull(id)

all1 <- read_csv("CompareModels/data/t6s_4k.csv")
deterministic <- read_csv("CompareModels/data/t6s_2k_deterministic.csv")

pd_pick <- c(2,7,18,30,35)

all1 <- all1 %>% pivot_longer(cols = c('S','Z1',"Z2",'nu1','nu2','FracI1','FracI2'),names_to = 'Species',values_to = "Column1")

deterministic <- deterministic %>% pivot_longer(cols = c('S','Z1',"Z2",'nu1','nu2','FracI1','FracI2'),names_to = 'Species',values_to = "Column1")
 
all_pop <- all1  %>% filter(pdoug %in% pd_pick, rep == 5, Year <= 4000)
all_pop <- all_pop %>% mutate(pd = recode(pdoug, '2' = '5% Douglas-fir', '7' = "20% Douglas-fir", '18' = '50% Douglas-fir', '30' = '80% Douglas-fir', '35' = '95% Douglas-fir'),
                              pop = recode(Species, "S" = "Host", "Z1" = "SNPV", "Z2" = "MNPV"))

all_pop$pd <- factor(all_pop$pd, levels = c('5% Douglas-fir',"20% Douglas-fir", '50% Douglas-fir', '80% Douglas-fir', '95% Douglas-fir'))

all_pop_wide <- all_pop %>% select(-Species) %>% 
  pivot_wider(names_from = pop, values_from = Column1)
all_pop_wide$pd <- factor(all_pop_wide$pd, levels = c('5% Douglas-fir',"20% Douglas-fir", '50% Douglas-fir', '80% Douglas-fir', '95% Douglas-fir'))

det_pop <- deterministic %>% filter(pdoug %in% pd_pick, rep == 5, Year <= 4000)
det_pop <- det_pop %>% mutate(pd = recode(pdoug, '2' = '5% Douglas-fir', '7' = "20% Douglas-fir", '18' = '50% Douglas-fir', '30' = '80% Douglas-fir', '35' = '95% Douglas-fir'),
                              pop = recode(Species, "S" = "Host", "Z1" = "SNPV", "Z2" = "MNPV"))

det_pop$pd <- factor(det_pop$pd, levels = c('5% Douglas-fir',"20% Douglas-fir", '50% Douglas-fir', '80% Douglas-fir', '95% Douglas-fir'))
det_pop_wide <- det_pop %>% select(-Species) %>% 
  pivot_wider(names_from = pop, values_from = Column1)
det_pop_wide$pd <- factor(det_pop_wide$pd, levels = c('5% Douglas-fir',"20% Douglas-fir", '50% Douglas-fir', '80% Douglas-fir', '95% Douglas-fir'))

year1 <- 150
year2 <- 190

snpv_col <- "#ee8800"
mnpv_col <-'#5D65C5'

plt_lines <- all_pop %>% filter(Species %in% c("S","Z1","Z2"), Year >= year1, Year <= year2) %>% 
  ggplot() + aes(x = Year, y= Column1, color = pop) + geom_line() + theme_classic(base_size = 15) + 
  facet_wrap(~pd,nrow = 1) + 
  scale_y_log10() +
  theme(axis.title.x = element_blank()) + 
  scale_color_manual("", values = c("Host" = "forestgreen", 'SNPV' = snpv_col, "MNPV" = mnpv_col)) + 
  ylab(expression(log[10] ~ "Population Size"))

legend_labels <- c(expression(bar(nu)["SNPV"]),expression(bar(nu)["MNPV"]))

nuSDO = 10.5e-3
nuSGR = 4.84e-3
nuMDO = 8.11e-3
nuMGR = 8.11e-3

nu_df <- data.frame(name =c('nuSDO','nuSGR','nuMDO','nuMGR'), values = c(nuSDO,nuSGR,nuMDO,nuMGR),
                    x1 = 140, x2 = 160,Species =c('nu1','nu1','nu2','nu2'), tree = c("DO","GR","DO","GR"))

plt_transm <- all_pop %>% filter(Species %in% c("nu1","nu2"), Year >= year1, Year <= year2) %>% 
  ggplot() + aes(x = Year, y= Column1, color = Species) + geom_line() + theme_classic(base_size = 15) + 
  facet_wrap(~pd,nrow = 1) + 
  scale_y_log10() +
  theme(axis.title.x = element_blank(),strip.text = element_blank(), strip.background = element_blank()) +
  scale_color_manual("", values = c(snpv_col, mnpv_col), labels = parse(text = legend_labels)) + 
  ylab(expression(log[10] ~ "Transmission Risk")) +
  geom_segment(data = nu_df, aes(x = x1, xend = x2, y = values, yend = values, linetype = tree), size = 1) +
  coord_cartesian(xlim = c(150,190), ylim = c(0.0001, 24)) + 
  scale_linetype_manual(values = c('DO' = 'solid','GR' = 'dashed'))

all_pop_wide$type <- "stochastic"
det_pop_wide$type <- "deterministic"

both_wide <- rbind(all_pop_wide,det_pop_wide)

plt_dots3 <- ggplot(both_wide[both_wide$Year <= 4000 & both_wide$Year <= 1000,], aes(x=log10(SNPV), y=log10(MNPV), z=log10(Host), color = type, size = type, alpha = type)) + 
  stat_3D(theta = 80, phi = 0) + 
  axes_3D(theta = 80, phi = 0) +
  #labs_3D(labs=c("SNPV", "MNPV", "Host")) +
  theme_void() + 
  facet_wrap(~pd,nrow = 1) + 
  theme(strip.text = element_blank()) + 
  scale_color_manual("", values = c("deterministic" = 'red', 'stochastic'= 'black')) + 
  scale_size_manual("", values = c("deterministic" = 0.5, 'stochastic'= 0.05)) + 
  scale_alpha_manual("", values = c("deterministic" = 1, 'stochastic'= 0.5))


plt_fracs <- all_pop_wide %>% mutate(MNPV = (FracI2)/(FracI1 + FracI2),
                                     SNPV = (FracI1)/(FracI1 + FracI2)) %>% 
  select(Year,pdoug,pd,MNPV,SNPV) %>% pivot_longer(cols = c("MNPV","SNPV"), names_to = "Species", values_to = 'val') %>% 
  filter(Year >= year1, Year <= year2) %>% 
  ggplot() + aes(x = Year, y= val*100, fill = Species, color = Species) + geom_bar(stat = 'identity', position = 'stack') + theme_classic(base_size = 15) + 
  facet_wrap(~pd,nrow = 1) + 
  theme(axis.title.x = element_blank(),strip.text = element_blank(), strip.background = element_blank()) + 
  scale_fill_manual("", values = c('SNPV' = snpv_col, "MNPV" = mnpv_col)) +
  scale_color_manual("", values = c('SNPV' = snpv_col, "MNPV" = mnpv_col)) + 
  ylab("% of Fraction Infected")

plt_tfrac <- all_pop_wide %>% mutate(MNPV = FracI2,
                        SNPV = FracI1) %>% 
  select(Year,pdoug,pd,MNPV,SNPV) %>%
  pivot_longer(cols = c("MNPV","SNPV"), names_to = "Species", values_to = 'val') %>% 
  filter(Year >= year1, Year <= year2) %>% 
  ggplot() + aes(x = Year, y= val, fill = Species, color = Species) + geom_bar(stat = 'identity', position = 'stack') + theme_classic(base_size = 15) + 
  facet_wrap(~pd,nrow = 1) + 
  theme(axis.title.x = element_blank(),strip.text = element_blank(), strip.background = element_blank()) + 
  scale_fill_manual("", values = c('SNPV' = snpv_col, "MNPV" = mnpv_col)) +
  scale_color_manual("", values = c('SNPV' = snpv_col, "MNPV" = mnpv_col)) + 
  ylab("Fraction Infected") + 
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 0, ymax = 0.3, fill = 'grey70', alpha = 0.025, color = NA)

high_frac <- read_csv('CompareModels/data/high_frac.csv')
high_frac$y1 <- 0
high_frac$y2 <- 0.5

plt_tfrac2 <- all_pop_wide %>% mutate(MNPV = FracI2,
                        SNPV = FracI1) %>% 
  select(Year,pdoug,pd,MNPV,SNPV) %>%
  pivot_longer(cols = c("MNPV","SNPV"), names_to = "Species", values_to = 'val') %>% 
  filter(Year >= year1, Year <= year2) %>% 
  ggplot() + aes(x = Year, y= val, fill = Species, color = Species) + geom_bar(stat = 'identity', position = 'stack') +
  theme_classic(base_size = 15) + 
  facet_wrap(~pd,nrow = 1) + 
  theme(axis.title.x = element_blank(),strip.text = element_blank(), strip.background = element_blank()) + 
  scale_fill_manual("", values = c('SNPV' = snpv_col, "MNPV" = mnpv_col)) +
  scale_color_manual("", values = c('SNPV' = snpv_col, "MNPV" = mnpv_col)) + 
  ylab("Fraction Infected") + 
  geom_rect(data = high_frac,
            aes(xmin = start2, xmax = end2, ymin = y1, ymax =y2, y = NULL), fill = 'white', alpha = 0.75, color = NA) +
  coord_cartesian(xlim = c(150,190)) + 
  geom_hline(yintercept = 0.3, linetype = 'dashed', color = 'grey65')

plt_box <- all_pop_wide %>% mutate(total = FracI1 + FracI2,
                                   MNPV = (FracI2)/(FracI1 + FracI2),
                                   SNPV = (FracI1)/(FracI1 + FracI2)) %>% filter(total >= 0.3) %>% 
  select(Year,pdoug,pd,MNPV,SNPV) %>% pivot_longer(cols = c("MNPV","SNPV"), names_to = "Species", values_to = 'val') %>% 
  filter(Year >= 50, Year <= 200) %>% 
  ggplot() + aes(x = Species, y= val*100, fill = Species) + geom_boxplot() + theme_classic(base_size = 15) + 
  facet_wrap(~pd,nrow = 1) + 
  theme(axis.title.x = element_blank(),strip.text = element_blank(), strip.background = element_blank()) + 
  scale_fill_manual("", values = c('SNPV' = snpv_col, "MNPV" = mnpv_col)) +
  ylab("% of Fraction Infected")


pdf("figures/simulations_over_pd_noext_fix2.pdf",height = 14, width = 14)
grid.arrange(plt_lines,plt_transm,plt_fracs,plt_tfrac2,plt_box,plt_dots3,nrow=6,heights = c(1,0.9,0.9,0.9,0.9,0.9))
dev.off()

small1 <- all1 %>% filter(pdoug %in% pd_pick, Year <= 200, Year >= 150, Species %in% c("S","Z1","Z2"))
small1 <- small1 %>% mutate(sp2 = recode(Species, "Z1" = "SNPV", "Z2" = "MNPV"))
small1$sp2 <- factor(small1$sp2, levels = c("S", "SNPV","MNPV"))

nu1 <- all1 %>% filter(pdoug %in% pd_pick, Year <= 200, Year >= 150, Species %in% c("nu1","FracI1"))

nu1 %>% ggplot() + aes(x = Year, y = Column1, color = Species, group = Species) +
  geom_line() + theme_classic(base_size = 15) +
  scale_color_manual("", values = c('nu1' = 'tomato2', "nu2" = 'dodgerblue3'),
                     labels = c(bquote(bar(nu)['SNPV']), bquote(bar(nu)['MNPV']))) +
  scale_y_log10() +
  facet_grid(rep~pdoug)

all1 %>% filter(Year >=150, Year <= 200, rep == 5) %>% pivot_wider(names_from = Species,values_from = Column1) %>% 
  ggplot() + aes(x = nu1,y = FracI1, color = Year) + geom_path() + theme_classic(base_size = 15) +
  scale_y_log10() + scale_x_log10() +
  facet_wrap(~pdoug,nrow = 1) + scale_color_viridis_c(option = 'magma')

all1 %>% filter(Year >=150, Year <= 200, rep == 5) %>% pivot_wider(names_from = Species,values_from = Column1) %>% 
  ggplot() + aes(x = nu2,y = FracI2, color = Year) + geom_path() + theme_classic(base_size = 15) +
  scale_y_log10() + scale_x_log10() +
  facet_wrap(~pdoug,nrow = 1) + scale_color_viridis_c(option = 'magma')

