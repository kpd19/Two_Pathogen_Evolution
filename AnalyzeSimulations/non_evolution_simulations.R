library(tidyverse)
library(gridExtra)

snpv_col <- "#ee8800"
mnpv_col <-'#5D65C5'

rlzn_wide <- read_csv('AnalyzeSimulations/data/simulation_no_evo.csv')

rlzn_long <- rlzn_wide %>% pivot_longer(cols = c("S","Z1","Z2","FracI1","FracI2"), names_to = "Species", values_to = 'Column1') 

rlzn_long <- rlzn_long %>% mutate(pop = recode(Species, "S" = "Host", "Z1" = "Single-capsid morphotype", "Z2" = "Multi-capsid morphotype",
                                  'FracI1' = 'FracS',"FracI2" = "FracM"))

rlzn_long$pd <- factor(rlzn_long$pd, levels = c('5% Douglas-fir',"20% Douglas-fir", '50% Douglas-fir', '80% Douglas-fir', '95% Douglas-fir'))

snpv_col <- "#ee8800"
mnpv_col <-'#5D65C5'

pdf("AnalyzeSimulations/figures/best_param_noevo_sig.pdf",height = 7, width = 10)
rlzn_long %>% filter(Species %in% c('Z1','Z2','S'), rep == 6, pdoug %in% c(2,7,18,30,35), sigma < 2) %>%
  ggplot() + aes(x = Year, y = Column1, color = pop) + 
  geom_line() + theme_classic(base_size = 15) + 
  scale_y_log10() + 
  scale_color_manual("", values = c("Host" = "forestgreen", 'Single-capsid morphotype' = snpv_col, "Multi-capsid morphotype" = mnpv_col)) + 
  facet_grid(sigma~pd, scales = 'free_y', labeller = label_bquote(rows = sigma == .(sigma))) +
  ylab("Population Density") + 
  theme(legend.position = 'top') + 
  xlim(0,100)
dev.off()
