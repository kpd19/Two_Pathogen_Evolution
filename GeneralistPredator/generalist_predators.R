library(tidyverse)

gen_pred <- read_csv("GeneralistPredator/data/generalist_predator_simulations.csv")
gen_pred_pmnpv <- read_csv("GeneralistPredator/data/generalist_predator_pMNPV.csv")
ll_data <- read_csv("GeneralistPredator/data/morphotype_dist_data.csv")
t6s_pmnpv <- read_csv("CompareModels/data/pmnpv_t6s_fix2noext.csv")
t6s_pmnpv <- t6s_pmnpv %>% filter(sigma == 0.5,pset == 7)

field_data_round <- ll_data %>% mutate(round_pd = round(Douglas_fir*9)/9) %>% group_by(round_pd) %>% 
  summarize(sum_MNPV = sum(MNPV), sum_tot = sum(total),
            pMNPV = mean(MNPV/total),
            sd_pMNPV = sd(MNPV/total),
            len = length(MNPV),
            pdoug = mean(Douglas_fir),
            sd_pd = sd(Douglas_fir),
            trees= mean(n_trees),
            sd_tree = sd(n_trees)) %>% 
  mutate(se_pMNPV = sd_pMNPV/sqrt(len),
         se_pd = sd_pd/sqrt(len))

pd_pick <- c(2,7,18,30,35)

snpv_col <- "#ee8800"
mnpv_col <-'#5D65C5'

gen_pred_long <- gen_pred %>% mutate(FI = FracI1 + FracI2) %>% 
  pivot_longer(cols = c('S','Z1',"Z2",'nu1','nu2','FracI1','FracI2'),names_to = 'Species',values_to = "Column1") %>% 
  mutate(pd = recode(pdoug, '2' = '5% Douglas-fir', '7' = "20% Douglas-fir", '18' = '50% Douglas-fir', '30' = '80% Douglas-fir', '35' = '95% Douglas-fir'))

gen_pred_long %>% filter(Species %in% c("S","Z1","Z2"), rep == 1, omega == 0.14,a %in% c(0.1,0.9)) %>%
  mutate(pop = recode(Species, "S" = "Host", "Z1" = "SNPV", "Z2" = "MNPV")) %>% 
  ggplot() + aes(x = Year, y = Column1, group = pop, color = pop) + geom_line() +
  theme_classic(base_size = 15) + 
  facet_grid(pd~a) + 
  scale_y_log10() +
  scale_color_manual("", values = c("Host" = "forestgreen", 'SNPV' = snpv_col, "MNPV" = mnpv_col)) 

legend_labels <- c(expression(bar(nu)["SNPV"]),expression(bar(nu)["MNPV"]))

gen_pred_long %>% filter(Species %in% c("nu1","nu2"), rep == 1, omega == 0.14,a %in% c(0.1,0.9)) %>%
  ggplot() + aes(x = Year, y = Column1, group = Species, color = Species) + geom_line() +
  scale_color_manual("", values = c(snpv_col, mnpv_col), labels = parse(text = legend_labels)) + 
  theme_classic(base_size = 15) + 
  facet_grid(pd~a) + 
  scale_y_log10()

gen_pred_long %>% filter(Species %in% c("FracI1",'FracI2'), rep == 1, omega == 0.14, a %in% c(0.1,0.9)) %>% 
  mutate(Species = recode(Species, "FracI1" = "SNPV", "FracI2" = 'MNPV')) %>% 
  ggplot() + aes(x = Year, y = Column1/FI, group = Species, fill = Species, color = Species) + geom_bar(stat = 'identity') +
  theme_classic(base_size = 15) + 
  facet_grid(a~pd) +
  scale_fill_manual("", values = c('SNPV' = snpv_col, "MNPV" = mnpv_col))  +
  scale_color_manual("", values = c('SNPV' = snpv_col, "MNPV" = mnpv_col)) 

gen_pred_avgs <- gen_pred_pmnpv %>% group_by(pdoug,omega,a) %>% summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE))

pdf("GeneralistPredator/figures/generalist_pred_pmnpv.pdf",height = 6, width = 10)
gen_pred_avgs %>% 
  ggplot() +
  theme_classic(base_size = 15) + 
  geom_point(data = field_data_round, aes(x = trees/37*100, y = pMNPV*100, size = sum_tot), show.legend = FALSE) +
  geom_errorbar(data = field_data_round, aes(x = trees/37*100, ymin = pMNPV*100 - se_pMNPV*100, ymax = pMNPV*100 + se_pMNPV*100), size = 1) +
  geom_errorbar(data = field_data_round, aes(y = pMNPV*100, xmin = trees/37*100 - se_pd*100, xmax = trees/37*100 + se_pd*100), size = 1) +
  geom_line(aes(x = pdoug/37*100, y = mean_MNPV*100, group = omega, color = as.factor(omega)),size = 2) +
  xlab("% Douglas-fir") + ylab("% Multi-capsid morphotype") +
  ylim(0,100) +
  scale_color_viridis_d(expression(omega), option = 'viridis') + 
  facet_wrap(~a,nrow = 2, labeller = label_bquote(rows = a == .(a))) +
  geom_line(data = t6s_pmnpv, aes(x = pdoug/37*100, y = mean_MNPV*100), linetype = 'dashed', color = 'red', size = 1) 
dev.off()
