library(tidyverse)
library(gridExtra)
library(gg3D)

snpv_col <- "#ee8800"
mnpv_col <-'#5D65C5'

biocontrol_grid <- read_csv("Biocontrol/data/biocontrol_data.csv")
biocontrol_t0.5p150 <- read_csv("Biocontrol/data/biocontrol_data_t0.5_p150.csv")

bicontrol_long <- biocontrol_grid %>% pivot_longer(cols = c('S','nu1','nu2'),
                                   names_to = 'Species',values_to = "Column1") 

pdf("Biocontrol/figures/transmission_risk_MNPV_bio.pdf",height = 8, width = 12)
bicontrol_long %>% filter(type == "100% multi-capsid morphotype", rep == 1, pdoug == 18, Species %in% c("nu2",'nu1')) %>% 
  group_by(rep,thresh,pdoug) %>%
  ggplot() + aes(x = Year,y = Column1, group = Species, color = Species) + geom_path() + theme_classic(base_size = 15) +
  scale_y_log10() + 
  facet_grid(biocontrol~thresh, labeller = label_bquote(rows = 'TMB' == .(biocontrol),
                                                        cols = 'Threshold' == .(thresh))) + 
  geom_vline(aes(xintercept = 150),linetype = 'dashed', color = 'grey55') + 
  geom_vline(aes(xintercept = 250),linetype = 'dashed', color = 'grey55') + 
  scale_color_manual("", values = c("nu1" = snpv_col, 'nu2' = mnpv_col)) +
  ylab("Transmission risk")
dev.off()

pdf("Biocontrol/figures/transmission_risk_SNPV_bio.pdf",height = 8, width = 12)
bicontrol_long %>% filter(type == "100% single-capsid morphotype", rep == 1, pdoug == 18, Species %in% c("nu2",'nu1')) %>% 
  group_by(rep,thresh,pdoug) %>%
  ggplot() + aes(x = Year,y = Column1, group = Species, color = Species) + geom_path() + theme_classic(base_size = 15) +
  scale_y_log10() + 
  facet_grid(biocontrol~thresh, labeller = label_bquote(rows = 'TMB' == .(biocontrol),
                                                        cols = 'Threshold' == .(thresh))) + 
  geom_vline(aes(xintercept = 150),linetype = 'dashed', color = 'grey55') + 
  geom_vline(aes(xintercept = 250),linetype = 'dashed', color = 'grey55') + 
  scale_color_manual("", values = c("nu1" = snpv_col, 'nu2' = mnpv_col)) 
dev.off()

thresh_pick <- 0.5
biocontrol_pick <- 150

legend_labels <- c(expression(bar(nu)["SNPV"]),expression(bar(nu)["MNPV"]))

pd_multi <- bicontrol_long %>% filter(type == "100% multi-capsid morphotype", rep == 1, thresh == thresh_pick, biocontrol == biocontrol_pick,
                  Species %in% c("nu2",'nu1'), Year >= 100, Year <= 300) %>% 
  mutate(pd = recode(pdoug, '2' = '5% Douglas-fir', '7' = "20% Douglas-fir",'13' = '35% Douglas-fir',
                     '18' = '50% Douglas-fir', '24' = '65% Douglas-fir', 
                     '30' = '80% Douglas-fir', '35' = '95% Douglas-fir')) %>% 
  mutate(pd_num = recode(pdoug, '2' = 1, '7' = 2,'13' = 3,
                         '18' = 4, '24' = 5, 
                         '30' = 6, '35' = 7)) %>% 
  mutate(pd = factor(pd, levels = c('5% Douglas-fir',"20% Douglas-fir",'35% Douglas-fir', '50% Douglas-fir','65% Douglas-fir', '80% Douglas-fir', '95% Douglas-fir'))) %>% 
  ggplot() + aes(x = Year,y = Column1, group = Species, color = Species) + geom_path() + theme_classic(base_size = 10) +
  scale_y_log10() + 
  facet_grid(~pd) + 
  geom_vline(aes(xintercept = 150),linetype = 'dashed', color = 'grey55') +
  geom_vline(aes(xintercept = 250),linetype = 'dashed', color = 'grey55') + 
  scale_color_manual("", values = c(snpv_col, mnpv_col), labels = parse(text = legend_labels)) + 
  ggtitle("100% multi-capsid morphotype")  +
  ylab("Transmission risk")
pd_single <- bicontrol_long %>% filter(type == "100% single-capsid morphotype", rep == 1, thresh == thresh_pick, biocontrol == biocontrol_pick,
                    Species %in% c("nu2",'nu1'), Year >= 100, Year <= 300) %>% 
  mutate(pd = recode(pdoug, '2' = '5% Douglas-fir', '7' = "20% Douglas-fir",'13' = '35% Douglas-fir',
                     '18' = '50% Douglas-fir', '24' = '65% Douglas-fir', 
                     '30' = '80% Douglas-fir', '35' = '95% Douglas-fir')) %>% 
  mutate(pd_num = recode(pdoug, '2' = 1, '7' = 2,'13' = 3,
                         '18' = 4, '24' = 5, 
                         '30' = 6, '35' = 7)) %>% 
  mutate(pd = factor(pd, levels = c('5% Douglas-fir',"20% Douglas-fir",'35% Douglas-fir', '50% Douglas-fir','65% Douglas-fir', '80% Douglas-fir', '95% Douglas-fir'))) %>% 
  ggplot() + aes(x = Year,y = Column1, group = Species, color = Species) + geom_path() + theme_classic(base_size = 10) +
  scale_y_log10() + 
  facet_grid(~pd) + 
  geom_vline(aes(xintercept = 150),linetype = 'dashed', color = 'grey55') + 
  geom_vline(aes(xintercept = 250),linetype = 'dashed', color = 'grey55') + 
  scale_color_manual("", values = c(snpv_col, mnpv_col), labels = parse(text = legend_labels)) + 
  ggtitle("100% single-capsid morphotype")  +
  ylab("Transmission risk")

pd_mix <- bicontrol_long %>% filter(type == "50% multi-capsid, 50% single-capsid", rep == 1, thresh == thresh_pick, biocontrol == biocontrol_pick,
                                 Species %in% c("nu2",'nu1'), Year >= 100, Year <= 300) %>% 
  mutate(pd = recode(pdoug, '2' = '5% Douglas-fir', '7' = "20% Douglas-fir",'13' = '35% Douglas-fir',
                     '18' = '50% Douglas-fir', '24' = '65% Douglas-fir', 
                     '30' = '80% Douglas-fir', '35' = '95% Douglas-fir')) %>% 
  mutate(pd_num = recode(pdoug, '2' = 1, '7' = 2,'13' = 3,
                         '18' = 4, '24' = 5, 
                         '30' = 6, '35' = 7)) %>% 
  mutate(pd = factor(pd, levels = c('5% Douglas-fir',"20% Douglas-fir",'35% Douglas-fir', '50% Douglas-fir','65% Douglas-fir', '80% Douglas-fir', '95% Douglas-fir'))) %>% 
  ggplot() + aes(x = Year,y = Column1, group = Species, color = Species) + geom_path() + theme_classic(base_size = 10) +
  scale_y_log10() + 
  facet_grid(~pd) + 
  geom_vline(aes(xintercept = 150),linetype = 'dashed', color = 'grey55') + 
  geom_vline(aes(xintercept = 250),linetype = 'dashed', color = 'grey55') + 
  scale_color_manual("", values = c(snpv_col, mnpv_col), labels = parse(text = legend_labels)) + 
  ggtitle("50% multi-capsid, 50% single-capsid")  +
  ylab("Transmission risk")

pdf("Biocontrol/figures/transmission_risk_pdoug_t0.5_p150.pdf",height = 9, width = 12)
grid.arrange(pd_multi,pd_mix,pd_single,nrow = 3)
dev.off()

plt_snpv <- biocontrol_t0.5p150 %>% mutate(type = factor(type,
                                 levels = c('No biocontrol applied','100% multi-capsid morphotype',
                                            '75% multi-capsid, 25% single-capsid', '50% multi-capsid, 50% single-capsid',
                                            '25% multi-capsid, 75% single-capsid',
                                            "100% single-capsid morphotype"))) %>% 
  ggplot() + aes(x = pdoug/37*100, y = nu1, group = interaction(pdoug,type), fill = type) +
  geom_boxplot(outliers = FALSE, width = 10) + theme_classic() + 
  scale_y_log10() +
  ylab(expression(bar(nu)[SNPV])) + 
  scale_fill_manual("", values = c('No biocontrol applied' = 'grey65',
                                   '100% multi-capsid morphotype' = '#08519c',
                                   '75% multi-capsid, 25% single-capsid' = '#3182bd',
                                   '50% multi-capsid, 50% single-capsid' ='#fef0d9',
                                   '25% multi-capsid, 75% single-capsid' = '#fc8d59',
                                   "100% single-capsid morphotype" = '#e34a33'))+ 
  xlab("% Douglas-fir") + ggtitle("Single-capsid morphotype transmission risk") + 
  theme(plot.title = element_text(hjust = 0.5))

plt_mnpv <- biocontrol_t0.5p150 %>% mutate(type = factor(type,
                               levels = c('No biocontrol applied','100% multi-capsid morphotype',
                                          '75% multi-capsid, 25% single-capsid', '50% multi-capsid, 50% single-capsid',
                                          '25% multi-capsid, 75% single-capsid',
                                          "100% single-capsid morphotype"))) %>% 
  ggplot() + aes(x = pdoug/37*100, y = nu2, group = interaction(pdoug,type), fill = type) +
  geom_boxplot(outliers = FALSE, width = 10) + theme_classic() + 
  scale_y_log10() +
  ylab(expression(bar(nu)[MNPV])) + 
  scale_fill_manual("", values = c('No biocontrol applied' = 'grey65',
                                   '100% multi-capsid morphotype' = '#08519c',
                                   '75% multi-capsid, 25% single-capsid' = '#3182bd',
                                   '50% multi-capsid, 50% single-capsid' ='#fef0d9',
                                   '25% multi-capsid, 75% single-capsid' = '#fc8d59',
                                   "100% single-capsid morphotype" = '#e34a33'))+ 
  xlab("% Douglas-fir") + ggtitle("Multi-capsid morphotype transmission risk") + 
  theme(plot.title = element_text(hjust = 0.5))

plt_host <- biocontrol_t0.5p150 %>% mutate(type = factor(type,
                               levels = c('No biocontrol applied','100% multi-capsid morphotype',
                                          '75% multi-capsid, 25% single-capsid', '50% multi-capsid, 50% single-capsid',
                                          '25% multi-capsid, 75% single-capsid',
                                          "100% single-capsid morphotype"))) %>% 
  ggplot() + aes(x = pdoug/37*100, y = S, group = interaction(pdoug,type), fill = type) +
  geom_boxplot(outliers = FALSE, width = 10) + theme_classic() + 
  scale_y_log10() +
  ylab(expression(log[10]~"Host Population Size")) + 
  scale_fill_manual("", values = c('No biocontrol applied' = 'grey65',
                                   '100% multi-capsid morphotype' = '#08519c',
                                   '75% multi-capsid, 25% single-capsid' = '#3182bd',
                                   '50% multi-capsid, 50% single-capsid' ='#fef0d9',
                                   '25% multi-capsid, 75% single-capsid' = '#fc8d59',
                                   "100% single-capsid morphotype" = '#e34a33'))+ 
  xlab("% Douglas-fir") + ggtitle("Host Population") + 
  theme(plot.title = element_text(hjust = 0.5))

pdf("Biocontrol/figures/host_transmission.pdf",height = 10, width = 12)
grid.arrange(plt_snpv, plt_mnpv,plt_host,nrow = 3)
dev.off()


