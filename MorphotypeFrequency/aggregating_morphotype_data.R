library(tidyverse)
library(gridExtra)

site_inf <- read_csv("data/site_infection_data_2019-2020.csv")
ll_data <- read_csv("data/literature_dist_data.csv")

# calculating the infection totals for all the natural sites that weren't included in biopesticide spray projects
site_inf_totals <- site_inf %>% filter(sprayed == FALSE) %>% group_by(state,latitude,longitude,site,year) %>%
  summarize(sum_MNPV = sum(MNPV,na.rm=TRUE), sum_SNPV = sum(SNPV,na.rm=TRUE),
            coinf = sum(coinfections,na.rm=TRUE),sum_total = sum(infected),
            sum_insect = sum(total),
            sum_infected = sum(infected)) %>% rename(MNPV = sum_MNPV,SNPV = sum_SNPV, Site = site) %>%
  mutate(total_iso = MNPV + SNPV,
         total_known = MNPV + SNPV - coinf,
         pSNPV = SNPV/total_iso,
         pMNPV = MNPV/total_iso,
         p_coinf = coinf/total_known) %>% 
  mutate(frac_inf = sum_infected/sum_insect,
         dom = max(pMNPV,pSNPV), ndom = min(pMNPV,pSNPV)) %>% mutate(rat = ndom/dom)

# analysing coinfections

mod1 <- glm(p_coinf ~ pMNPV, data = site_inf_totals, family = 'quasibinomial')
mod2 <- glm(p_coinf ~ pSNPV, data = site_inf_totals, family = 'quasibinomial')
mod3 <- glm(p_coinf ~ rat, data = site_inf_totals, family = 'quasibinomial')
mod4 <- glm(p_coinf ~ 1, data = site_inf_totals, family = 'quasibinomial')
summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)

anova(mod1,mod4, test = "Chisq")
anova(mod2,mod4, test = "Chisq")
anova(mod3,mod4, test = "Chisq")

pdf("figures/prop_coinf_ratio.pdf",height = 10, width = 10)
plt1 <- site_inf_totals %>% ungroup() %>% mutate(pred = pSNPV*pMNPV) %>%
  ggplot() + aes(x = rat,y = p_coinf) + geom_point() + theme_classic(base_size = 15) + 
  xlab("Non-dominant morph: Dominant morph") + ylab("Prop. coinf") +
  geom_smooth(method = 'glm')
plt2 <- site_inf_totals %>% ungroup() %>% mutate(pred = pSNPV*pMNPV) %>%
  ggplot() + aes(x = frac_inf,y = p_coinf) + geom_point() + theme_classic(base_size = 15) + 
  xlab("Fraction infected") + ylab("Prop. coinf") +
  geom_smooth(method = 'glm')
plt3 <- site_inf_totals %>% ungroup() %>% mutate(pred = pSNPV*pMNPV) %>%
  ggplot() + aes(x = pSNPV,y = p_coinf) + geom_point() + theme_classic(base_size = 15) + 
  xlab("pSNPV") + ylab("Prop. coinf") +
  geom_smooth(method = 'gam')
plt4 <- site_inf_totals %>% ungroup() %>% mutate(pred = pSNPV*pMNPV) %>%
  ggplot() + aes(x = pMNPV,y = p_coinf) + geom_point() + theme_classic(base_size = 15) + 
  xlab("pMNPV") + ylab("Prop. coinf") +
  geom_smooth(method = 'gam')

grid.arrange(plt1,plt2,plt3,plt4,nrow=2)
dev.off()

site_summary <- site_inf_totals %>% filter(total_iso >0) %>% arrange(Site) %>%
  select(Site,MNPV,SNPV,total_iso,state,latitude,longitude) %>% 
  rename(site = Site, total = total_iso) %>% mutate(elevation = NA, source = 'Dwyer Lab PCR Collections')
site_summary$site_no <- 111:128

ll_data <- rbind(ll_data,site_summary)

write_csv(ll_data, "data/morphotye_distribution_data.csv")
