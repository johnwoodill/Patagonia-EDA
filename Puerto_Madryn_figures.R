library(lubridate)
library(tidyverse)
library(feather)
library(ggthemes)


dat = as.data.frame(read_feather('~/Data/GFW_point/Patagonia_Shelf/complete/Puerto_Madryn_2016-03-10_2016-03-20.feather'))
# dat = as.data.frame(read_feather('~/test.feather'))
head(dat)

dat <- filter(dat, rank <= 5)

dat <- dat %>% 
  mutate(day = day(timestamp),
         hour = hour(timestamp),
         ln_distance = log(1 + distance)) 

bdat <- data.frame()
# Bootstrap n times
for (i in 1:1000){
  ships <- sample(unique(dat$vessel_A), 45, )
  outdat <- dat %>% 
    filter(vessel_A %in% ships) %>% 
    group_by(day, hour) %>% 
    summarise(mean_dist = mean(ln_distance, na.rm=TRUE),
              med_dist = median(ln_distance, na.rm=TRUE),
              sd_dist = sd(ln_distance, na.rm=TRUE)) %>% 
  ungroup
  outdat$bstrap <- i
  bdat <- rbind(bdat, outdat)
}

ggplot(bdat, aes(mean_dist)) +
  geom_density() +
  facet_wrap(~day)

mdat <- bdat %>% 
  group_by(day) %>% 
  summarise(mean = mean(mean_dist),
            median = median(mean_dist),
            kurtosis = kurtosis(mean_dist),
            skewness = skewness(mean_dist)) %>% 
  gather(moment, value, -day)
  
mdat
ggplot(mdat, aes(x=day, y=value)) +
  theme_tufte(11) + 
  geom_line()+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "grey") +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "grey") +
  scale_x_continuous(breaks = seq(9, 20, 1)) + 
  geom_vline(xintercept = 15, color = "red") +
  facet_wrap(~moment)



pdat <- bdat %>% 
  group_by(day, hour) %>% 
  summarise(bmean_dist = mean(mean_dist),
            bmedian_dist = median(med_dist),
            bmean_sd = sd(mean_dist),
            bmed_sd = sd(med_dist))


ggplot(bdat) +
  theme_tufte(11) +
  geom_tufteboxplot(aes(x=factor(day), y=mean_dist)) +
  #geom_point(aes(x=hour, y=bmean_dist)) +
  #geom_line(aes(x=hour, y=bmean_dist)) +
  #geom_line(aes(x=hour, y=(bmean_dist + 1.96*bmean_sd)), linetype = "dashed") +
  #geom_line(aes(x=hour, y=(bmean_dist - 1.96*bmean_sd)), linetype = "dashed") +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "grey") +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "grey") +
  # scale_x_continuous(breaks = seq(9, 20, 1)) +
  ylab("Mean Distance") +
  xlab("March")
  
  
# ggplot(pdat) +
#   theme_tufte(11) +
#   geom_point(aes(x=day, y=bmedian_dist)) +
#   geom_line(aes(x=day, y=bmedian_dist)) +
#   geom_line(aes(x=day, y=(bmedian_dist + 1.96*bmed_sd)), linetype = "dashed") +
#   geom_line(aes(x=day, y=(bmedian_dist - 1.96*bmed_sd)), linetype = "dashed") +
#   annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "grey") +
#   annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "grey") +
#   scale_x_continuous(breaks = seq(9, 20, 1)) +
#   ylab("Mean Distance") +
#   xlab("March")
#                           
#                           
# 
# 
# 
# 
# ggplot(pdat, aes(log(distance))) +
#   geom_density() +
#   annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "grey") +
#   annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "grey") +
#   facet_wrap(~day)
# 
# ships <- dat %>%
#   group_by(day) %>%
#   summarise(nships = length(unique(vessel_A)))
# 
# ggplot(ships, aes(x=day, y=nships)) + 
#   geom_line()
# 
# # Sample n ships per day
# 
