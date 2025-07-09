library(tidyverse)

phi <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/cycle_test/vary_phi_all.csv")
phi2 <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/cycle_test/vary_phi2_all.csv")

phi <- rbind(phi,phi2)

phi_lls <- read_csv("/Volumes/My Book/Spatial_DFTM_modeling/Julia_sims/cycle_test/vary_phi_many_ll.csv")

ll_data <- read_csv("data/data_for_ll.csv")
field_data_round <- ll_data %>% mutate(round_pd = round(Douglas_fir*9)/9) %>% group_by(round_pd) %>% 
  summarize(sum_MNPV = sum(MNPV), sum_tot = sum(total),
            mean_MNPV = mean(MNPV/total),
            sd_pMNPV = sd(MNPV/total),
            len = length(MNPV),
            pdoug = mean(Douglas_fir),
            sd_pd = sd(Douglas_fir),
            pdoug= mean(n_trees),
            sd_tree = sd(n_trees)) %>% 
  mutate(se_pMNPV = sd_pMNPV/sqrt(len),
         se_pd = sd_pd/sqrt(len))

id1 <- ll_data %>% filter(n_trees %in% 1:36) %>%
  arrange(n_trees) %>% group_by(n_trees) %>%
  mutate(n = seq(1,length(total))) %>% filter(n == 1) %>% pull(id)

size_pick <- 6

q1 <- 0.025
q2 <- 0.975

phi_lls <- phi_lls %>% filter(id %in% id1) %>% select(quality,pdoug,mean_pMNPV,phi)

write_csv(phi_lls, "AnalyzeSimulations/data/phi_pmnpv.csv")

phi_pmnpv_2k <- read_csv("AnalyzeSimulations/data/phi_pmnpv.csv")

phi_avg <- phi_pmnpv_2k %>% filter(quality == 1) %>% group_by(pdoug,phi) %>%
  summarize(mean_MNPV = mean(mean_pMNPV,na.rm=TRUE),
            sd = sd(mean_pMNPV,na.rm=TRUE),
            q1 = quantile(mean_pMNPV, probs = c(q1), na.rm=TRUE),
            q9 = quantile(mean_pMNPV, probs = c(q2), na.rm=TRUE))

pdf("AnalyzeSimulations/figures/phi.pdf",height = 6, width = 12)
phi_avg %>% 
  filter(phi >= 5, phi <=90) %>% 
  ggplot()+
  geom_ribbon(aes(x = pdoug/37*100, ymax = q9*100, ymin = q1*100),
              fill = 'blue', color = 'black', alpha = 0.25,linetype = 'dashed') +
  geom_point(data = field_data_round, aes(x = pdoug/37*100, y = mean_MNPV*100, size = sum_tot), color = 'black') +
  geom_errorbar(data = field_data_round, aes(x = pdoug/37*100, ymin = mean_MNPV*100 - se_pMNPV*100, ymax = mean_MNPV*100 + se_pMNPV*100), size = 1, color = 'black') +
  geom_errorbar(data = field_data_round, aes(y = mean_MNPV*100, xmin = pdoug/37*100 - se_pd*100, xmax = pdoug/37*100 + se_pd*100), size = 1, color = 'black') +
  aes(x = pdoug/37*100, y = mean_MNPV*100) + geom_line(size = 1, color = 'black') +
  theme_classic(base_size = 15) + 
  facet_wrap(~phi, labeller = label_bquote(rows = phi == .(phi)),nrow = 2) +
    scale_color_brewer(expression(phi), palette = "Spectral") +
  xlim(0,100) + ylim(0,100) + 
  xlab("% Douglas-fir") + ylab("% Multi-capsid morphotype") +
  guides(size = 'none')
dev.off()

snpv_col <- "#ee8800"
mnpv_col <-'#5D65C5'

phi <- phi %>% filter(pdoug %in% c(2,7,18,30,35), rep == 9)

write_csv(phi,"AnalyzeSimulations/data/phi_dynamics.csv")

phi_dynamics <- read_csv("AnalyzeSimulations/data/phi_dynamics.csv")

pdf("AnalyzeSimulations/figures/cycles_phi_pd.pdf",height = 15, width = 10)
phi_dynamics %>% filter(pdoug %in% c(2,7,18,30,35),rep == 9, Year >= 130, Year <= 170, phi <= 90) %>% select(Year,Z1, Z2, S,pdoug,phi,rep) %>% 
  mutate(pd = recode(pdoug, '2' = '5% Douglas-fir', '7' = "20% Douglas-fir",'13' = '35% Douglas-fir',
                     '18' = '50% Douglas-fir', '24' = '65% Douglas-fir', 
                     '30' = '80% Douglas-fir', '35' = '95% Douglas-fir')) %>% 
  mutate(pd = factor(pd, levels = c('5% Douglas-fir',"20% Douglas-fir",'35% Douglas-fir', '50% Douglas-fir','65% Douglas-fir', '80% Douglas-fir', '95% Douglas-fir'))) %>% 
  pivot_longer(cols = c('S',"Z1","Z2"), names_to = 'Species', values_to = 'Column1') %>% 
  mutate(Species = recode(Species, "S" = "Host", "Z1" = "Single-capsid morphotype", "Z2" = "Multi-capsid morphotype")) %>% 
  ggplot() + aes(x = Year, y = Column1, color = Species, group = Species) + geom_line() + 
  theme_classic(base_size = 15) + 
  facet_grid(phi~pd, scales = 'free', labeller = label_bquote(rows = phi == .(phi))) + 
  scale_y_log10() + 
  scale_color_manual(values = c("Host" = 'forestgreen', 'Single-capsid morphotype' = snpv_col, "Multi-capsid morphotype" = mnpv_col)) +
  ylab("Population density") + 
  theme(legend.position = 'top')
dev.off()


