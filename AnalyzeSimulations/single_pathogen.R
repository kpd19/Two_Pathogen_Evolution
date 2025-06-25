library(tidyverse)
library(gridExtra)
library(WaveletComp)

find_peaks2 <- function (x, m = 1, thresh=1){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  
  if (!missing(thresh)) {
    pks <- pks[which(x[pks]>=thresh)]
  }
  else pks
  
  return(pks)
}

snpv_col <- "#ee8800"
mnpv_col <-'#5D65C5'

trio <- read_csv('AnalyzeSimulations/data/single_pathogen.csv')

trio_long <- trio %>% select(Year,S,nu1,nu2,fecundity,type,rep,pdoug) %>% pivot_longer(cols = c('S','nu1','nu2','fecundity'))

trio_long <- trio_long %>% mutate(pd = recode(pdoug, '2' = '5% Douglas-fir', '7' = "20% Douglas-fir",'13' = '35% Douglas-fir',
                                                '18' = '50% Douglas-fir', '24' = '65% Douglas-fir', 
                                                '30' = '80% Douglas-fir', '35' = '95% Douglas-fir')) %>% 
  mutate(pd = factor(pd, levels = c('5% Douglas-fir', "20% Douglas-fir", '35% Douglas-fir',
                                    '50% Douglas-fir', '65% Douglas-fir', 
                                    '80% Douglas-fir', '95% Douglas-fir')))

single_path <- trio_long

rp_pick <- 2

plt1 <- single_path %>% filter(rep == rp_pick, pdoug %in% c(2,7,18,30,35), name == 'nu1', Year >= 160, Year <=200) %>%
  ggplot() + aes(x = Year, y = value, color = type) + 
  geom_line() + theme_classic(base_size = 15) + 
  scale_y_log10() + 
  facet_wrap(~pd,nrow=1) + 
  ggtitle("Single-capsid morphotype") + 
  ylab("Average Infectiousness")+
  scale_color_manual("", values = c("Both pathogens" = "forestgreen", 'Single-capsid morphotype only' = snpv_col)) +
  theme(legend.position = 'none')

plt2 <- single_path %>% filter(rep == rp_pick, pdoug %in% c(2,7,18,30,35), name == 'nu1') %>%
  ggplot() + aes(x = Year, y = value, color = type) + 
  geom_boxplot(alpha = 0.75) + theme_classic(base_size = 15) + 
  scale_y_log10() + 
  facet_wrap(~pd,nrow=1) + 
  ylab("Average Infectiousness")+
  scale_fill_manual("", values = c("Both pathogens" = "forestgreen", 'Single-capsid morphotype only' = snpv_col)) + 
  scale_x_continuous(breaks = c(55,145), labels = c('Both',"SNPV only")) + 
  theme(axis.title.x = element_blank())

plt3 <- single_path %>% filter(rep == rp_pick, pdoug %in% c(2,7,18,30,35), name == 'nu2', Year >= 160, Year <=200) %>%
  ggplot() + aes(x = Year, y = value, color = type) + 
  geom_line() + theme_classic(base_size = 15) + 
  scale_y_log10() + 
  facet_wrap(~pd,nrow=1) + 
  ggtitle("Multi-capsid morphotype") + 
  ylab("Average Infectiousness")+
  scale_color_manual("", values = c("Both pathogens" = "forestgreen", "Multi-capsid morphotype only" = mnpv_col)) +
  theme(legend.position = 'none')

plt4 <- single_path %>% filter(rep == rp_pick, pdoug %in% c(2,7,18,30,35), name == 'nu2') %>%
  ggplot() + aes(x = Year, y = value, fill = type) + 
  geom_boxplot(alpha = 0.75) + theme_classic(base_size = 15) + 
  scale_y_log10() + 
  facet_wrap(~pd,nrow=1) + 
  ylab("Average Infectiousness")+
  scale_color_manual("", values = c("Both pathogens" = "forestgreen", "Multi-capsid morphotype only" = mnpv_col)) +
  scale_x_continuous(breaks = c(55,145), labels = c('Both',"MNPV only")) + 
  theme(axis.title.x = element_blank())

plt5 <- single_path %>% filter(rep == rp_pick, pdoug %in% c(2,7,18,30,35), name %in% c('nu1','nu2')) %>%
  ggplot() + aes(x = name, y = value, fill = type, group = interaction(type,name)) + 
  geom_boxplot(alpha = 0.75, outliers = FALSE) + theme_classic(base_size = 15) + 
  scale_y_log10() + 
  facet_wrap(~pd,nrow=1) + 
  ylab("Average Infectiousness")+
  scale_fill_manual("", values = c("Both pathogens" = "forestgreen", "Multi-capsid morphotype only" = mnpv_col,
                                   'Single-capsid morphotype only' = snpv_col)) +
  scale_x_discrete(labels = c("SNPV",'MNPV')) + 
  theme(axis.title.x = element_blank(),
        legend.position = 'bottom') 

plt6 <- single_path %>% filter(rep == rp_pick, pdoug %in% c(2,7,18,30,35), name == 'fecundity', Year >= 50, Year <=200) %>%
  ggplot() + aes(x = type, y = value, fill = type, group = type) + 
  geom_boxplot(alpha = 0.75, outliers = FALSE) + theme_classic(base_size = 15) + 
  scale_y_log10() + 
  facet_wrap(~pd,nrow=1) + 
  ylab("Fecundity") +
  scale_fill_manual("", values = c("Both pathogens" = "forestgreen", "Multi-capsid morphotype only" = mnpv_col,
                                   'Single-capsid morphotype only' = snpv_col)) +
  scale_x_discrete(labels = c('Both',"MNPV\nonly",'SNPV\nonly')) + 
  theme(axis.title.x = element_blank())

pdf("AnalyzeSimulations/figures/single_path_nu_bar.pdf",height = 10, width = 10)
grid.arrange(plt1,plt3,plt5,nrow=3,heights = c(1,1,1.1))
dev.off()

a.time <- Sys.time()

host_pop <- single_path %>% filter(name == "S")
sim_types <- unique(host_pop$type)

all_amps <- c()
all_periods <- c()
all_peaks <- c()
for(s in 1:length(sim_types)){
  mini_S <- host_pop %>% filter(Year >= 50, Year <= 200, rep <= 20, type == sim_types[s]) %>%
    mutate(logP = log10(value)) %>% select(Year,logP,rep,pdoug)
  
  periods_S <- c()
  amps_S <- c()
  peaks_S <- c()
  
  periods <- c()
  amps <- list()
  
  pdougs <- sort(unique(mini_S$pdoug))
  
  for(j in 1:length(pdougs)){
    mini_S2 <- mini_S %>% filter(pdoug == pdougs[j])
    for(i in 1:20){
      test <- mini_S2 %>% filter(rep ==i)
      invisible(capture.output(my_wave <- analyze.wavelet(test,'logP', loess.span=0,dt =1,upperPeriod = 50, n.sim = 50, verbose = FALSE)))
      periods[i] <- my_wave$Period[which.max(my_wave$Power.avg)]
      
      tpeak <- find_peaks2(test$logP,m=2)
      peaks <- test$logP[tpeak]
      
      tmin <- find_peaks2(-test$logP,m=2)
      mins <- test$logP[tmin]
      
      if(is.null(tpeak) == TRUE || is.null(tmin) == TRUE){
        temp_amps <- data.frame(rep = i, amp = NA,pdoug = pdougs[j])
      }else{
        mm_df <- data.frame(t = c(tpeak,tmin), lpop = c(peaks,mins),type = c(rep(1,length(tpeak)),rep(0,length(tmin))))
        amps <- mm_df %>% arrange(t) %>% mutate(pp_diff = c(1,diff(type))) %>% filter(pp_diff !=0) %>%
          mutate(pop = 10^lpop) %>% mutate(diff = c(0,diff(pop))) %>% filter(diff >0) %>% pull(diff)
        
        temp_amps <- data.frame(rep = i, amp = amps,pdoug = pdougs[j])
      }
      
      amps_S <- rbind(amps_S,temp_amps)
      
      if(is.null(tpeak) == TRUE || is.null(tmin) == TRUE){
        temp_peaks <- data.frame(tpeak = NA, peak = NA, pdoug = pdougs[j])
      }else{
        temp_peaks <- data.frame(tpeak = tpeak, peak = peaks,pdoug = pdougs[j])
      }
      peaks_S <- rbind(peaks_S,temp_peaks)    
    }
    
    periods_S <- rbind(periods_S, data.frame(periods,rep=1:20,pdoug = pdougs[j]))
    print(paste0("Done pdoug: ", pdougs[j]))

  }
  
  amps_S$Species <-"S"
  periods_S$Species <-"S"
  peaks_S$Species <-"S"
  
  amps_S$type <- sim_types[s]
  periods_S$type <-sim_types[s]
  peaks_S$type <-sim_types[s]
  
  all_amps <- rbind(all_amps,amps_S)
  all_periods <- rbind(all_periods, periods_S)
  all_peaks <- rbind(all_peaks, peaks_S)
}

b.time <- Sys.time()
print(b.time-a.time)

all_peaks <- all_peaks %>% filter(pdoug %in% c(2,7,13,18,24,30,35)) %>%
  mutate(pd = recode(pdoug, '2' = '5% Douglas-fir', '7' = "20% Douglas-fir",'13' = '35% Douglas-fir',
                     '18' = '50% Douglas-fir', '24' = '65% Douglas-fir', 
                     '30' = '80% Douglas-fir', '35' = '95% Douglas-fir')) %>% 
  mutate(pd_num = recode(pdoug, '2' = 1, '7' = 2,'13' = 3,
                     '18' = 4, '24' = 5, 
                     '30' = 6, '35' = 7))

all_peaks$pd <- factor(all_peaks$pd, levels = c('5% Douglas-fir',"20% Douglas-fir",'35% Douglas-fir', '50% Douglas-fir','65% Douglas-fir', '80% Douglas-fir', '95% Douglas-fir'))

pdf("AnalyzeSimulations/figures/host_peaks.pdf",height = 6, width = 10)
all_peaks %>% 
  ggplot() + aes(x = pdoug/37*100, y = peak, group=interaction(pdoug,type),fill = type) +
  geom_boxplot(alpha = 0.75, outliers = FALSE, width = 8) + theme_classic(base_size = 15) + 
  ylab(expression(log[10] ~ "Peak Host Density")) +
  scale_fill_manual("", values = c("Both pathogens" = "forestgreen", 'Single-capsid morphotype only' = snpv_col,
                                   "Multi-capsid morphotype only" = mnpv_col)) + 
  xlab("% Douglas-fir") + 
  theme(legend.position = 'top')
dev.off()




find_peaks2 <- function (x, m = 1, thresh=1){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  
  if (!missing(thresh)) {
    pks <- pks[which(x[pks]>=thresh)]
  }
  else pks
  
  return(pks)
}

x <- seq(0,50,0.1) # x is your time, but you don't really need it except to plot
y <- sin(x)*3 # y is your time series, so probably comparejulia[comparejulia$tree == "Trembling Aspen",]$host

plot(x,y)

time_of_peak <- find_peaks2(y,m=2) # finds the locations of the peaks in the array
amplitudes <- y[time_of_peak] # finds the actual peak values
periods <- diff(time_of_peak) # this works because every time step of yours is a year

avg_period = mean(periods)
avg_amp <- mean(amplitudes)
avg_pop = mean(y)





