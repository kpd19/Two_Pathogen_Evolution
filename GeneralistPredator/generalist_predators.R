library(tidyverse)

pred <- read_csv("data/predator_test1.csv")
pred2 <- read_csv("spatial_model/data/predator_all2.csv")
lls <- read_csv("spatial_model/data/predator_ll1.csv")


ll_data <- read_csv("CompareModels/data/morphotype_dist_data.csv")

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


small <- pred %>% filter(pdoug %in% c(18))

small %>% filter(Species %in% c("S","Z1","Z2"), rep == 1) %>%
  ggplot() + aes(x = Year, y = Column1, group = Species, color = Species) + geom_line() +
  theme_classic(base_size = 15) + 
  facet_grid(omega~a) + 
  scale_y_log10()

stats <- pred %>% filter(Species %in% c('FracI1','FracI2'))  %>% filter(Year>=50) %>% 
  pivot_wider(names_from = Species, values_from = Column1) %>% mutate(FI = FracI1 + FracI2, n = 1) %>% 
  filter(FI >= 0.3) %>% group_by(pdoug,omega,a,rep) %>%
  summarize(mean_MNPV = mean(FracI2/FI),
            sum_n = sum(n),
            mean_FI = mean(FI))

stats_all <- stats %>% group_by(pdoug,omega,a) %>% summarize(mean_reps = mean(mean_MNPV),
                                                             mean_FI = mean(mean_FI))

pdf("spatial_model/figures/gen_pred_test.pdf",height = 6, width = 8)
stats_all %>% ggplot() + aes(x = pdoug, y = mean_reps, group = omega, color = omega) +
  geom_line() + theme_classic(base_size = 15) + 
  facet_wrap(~a,nrow = 2, labeller = label_bquote(rows = a == .(a))) + 
  scale_color_viridis_c(expression(omega), option = 'viridis') + 
  ylab("% Multi-capsid morphotype") + xlab("% Douglas-fir")

stats_all %>% ggplot() + aes(x = pdoug, y = mean_FI, group = omega, color = omega) +
  geom_line() + theme_classic(base_size = 15) + 
  facet_wrap(~a,nrow = 2, labeller = label_bquote(rows = a == .(a))) + 
  scale_color_viridis_c(expression(omega), option = 'viridis') + 
  ylab("% Multi-capsid morphotype") + xlab("% Douglas-fir")
dev.off()

lls %>% ggplot() + aes(x = pdoug, y = mean_pMNPV) + geom_point() + 
  theme_classic(base_size = 15)

lls %>% group_by(pdoug,omega,a) %>% summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE)) %>% 
  ggplot() + aes(x = pdoug, y = mean_MNPV, group = omega, color = omega) +
  geom_line() + theme_classic(base_size = 15) + 
  facet_wrap(~a,nrow = 2, labeller = label_bquote(rows = a == .(a))) + 
  scale_color_viridis_c(expression(omega), option = 'viridis') + 
  ylab("% Multi-capsid morphotype") + xlab("% Douglas-fir")
ll_avgs <- lls %>% group_by(pdoug,omega,a) %>% summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE))

pdf("spatial_model/figures/gen_pred_test2.pdf",height = 6, width = 10)
ll_avgs %>% 
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
