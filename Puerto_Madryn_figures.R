library(lubridate)
library(tidyverse)
library(feather)
library(ggthemes)
library(cowplot)

dat = as.data.frame(read_feather('~/Data/GFW_point/Patagonia_Shelf/complete/Puerto_Madryn_2016-03-10_2016-03-20.feather'))


# dat = as.data.frame(read_feather('~/test.feather'))
head(dat)

#dat <- filter(dat, rank <= 5)
dat$date <- paste0(year(dat$timestamp), "-", month(dat$timestamp), "-", day(dat$timestamp), "-", hour(dat$timestamp))
#dat$date <- paste0(year(dat$timestamp), "-", month(dat$timestamp), "-", day(dat$timestamp), "-")
dat <- dat %>% 
  mutate(day = day(timestamp),
         hour = hour(timestamp),
         ln_distance = log(distance)) 

bdat <- data.frame()
# Bootstrap n times
for (i in 1:1000){
  ships <- sample(unique(dat$vessel_A), 45)
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
  ylab("Log(Distance)") +
  xlab("March")
  

dat = as.data.frame(read_feather('~/Data/GFW_point/Patagonia_Shelf/complete/Puerto_Madryn_2016-03-10_2016-03-20.feather'))
dat <- filter(dat, rank <= 5)

#dat <- filter(dat, timestamp == '2016-03-09 16:00:00')
dat$date <- paste0(year(dat$timestamp), "-", month(dat$timestamp), "-", day(dat$timestamp), "-", hour(dat$timestamp))
# dat$date <- paste0(year(dat$timestamp), "-", month(dat$timestamp), "-", day(dat$timestamp))

dat <- dat %>% 
  mutate(day = day(timestamp),
         hour = hour(timestamp)) %>% 
  group_by(timestamp, date, vessel_A) %>% 
  summarise(distance = mean(distance)) %>% 
  mutate(ln_distance = log(1 + distance))
head(dat)  


outdat <- data.frame()
for (i in unique(dat$date)){
  indat <- filter(dat, date == i)
  efunc <- ecdf(indat$ln_distance)
  cdat <- cut(indat$ln_distance, breaks = 30)
  indat$breaks <- cdat 
  indat$n_breaks <- as.numeric(indat$breaks)
  
  retdat <- indat %>% 
    group_by(n_breaks) %>% 
    summarise(prob = efunc(max(ln_distance)) - efunc(min(ln_distance))) %>% 
    mutate(prob = prob/sum(prob)) %>% 
    ungroup
  
  retdat$timestamp <- indat$timestamp[1]
  retdat$date <- indat$date[1]
  outdat <- rbind(outdat, retdat)
}

outdat %>% 
  group_by(timestamp) %>% 
  summarise(sum = sum(prob))

ggplot(outdat, aes(timestamp, factor(n_breaks))) + 
  geom_tile(aes(fill = prob)) + 
  theme_tufte(11) +
  #scale_x_datetime(date_breaks = "1 day", rotate=45) +     ### changed
  #geom_raster(aes(fill = prob)) +
  scale_fill_gradient(low = "white", high = "red") +
  #scale_x_continuous(labels = seq(10, 21, 1)) +
  xlab("March") +
  ylab(paste0("Distance Breaks"))

ggsave("~/Projects/Patagonia-EDA/figures/heatmap_NN5.pdf", width = 6, height = 4)
{
# Medoids data outcomes
meddat1 <- as.data.frame(read_feather('~/Projects/Patagonia-EDA/data/k_medoids_day_NN1_3-10-2016-3-20-2016.feather'))
meddat2 <- as.data.frame(read_feather('~/Projects/Patagonia-EDA/data/k_medoids_day_NN5_3-10-2016-3-20-2016.feather'))
meddat3 <- as.data.frame(read_feather('~/Projects/Patagonia-EDA/data/k_medoids_day_hour_NN1_3-10-2016-3-20-2016.feather'))
meddat4 <- as.data.frame(read_feather('~/Projects/Patagonia-EDA/data/k_medoids_day_hour_NN5_3-10-2016-3-20-2016.feather'))
meddat1$group <- meddat1$group + 1
meddat2$group <- meddat2$group + 1
meddat3$group <- meddat3$group + 1
meddat4$group <- meddat4$group + 1

meddat1$value <- meddat1$value + 10
meddat2$value <- meddat2$value + 10


p1 <- ggplot(meddat1, aes(x=value, y=group)) +
  geom_point(size=.5) + 
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "grey") +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "grey") +
  geom_vline(xintercept = 15, color='red') +
  theme_tufte(11) +
  xlab(NULL) +
  ylab("Cluster") +
  annotate("text", x=10.5, y=3, label = "Daily", size=2.5) +
  annotate("text", x=10.5, y=2.8, label = "NN=1", size=2.5) +
  scale_y_continuous(breaks = seq(1, length(unique(meddat1$group)))) +
  scale_x_continuous(breaks = meddat1$value)
p1
p2 <- ggplot(meddat2, aes(x=value, y=group)) +
  geom_point(size=.5) + 
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "grey") +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "grey") +
  geom_vline(xintercept = 15, color='red') +
  theme_tufte(11) +
  xlab(NULL) +
  ylab("Cluster") +
  annotate("text", x=10.5, y=3, label = "Daily", size=2.5) +
  annotate("text", x=10.5, y=2.8, label = "NN=5", size=2.5) +
  scale_y_continuous(breaks = seq(1, length(unique(meddat1$group)))) +
  scale_x_continuous(breaks = meddat1$value)

plot_grid(p1, p2, ncol=2)

p3 <- ggplot(meddat3, aes(x=value, y=group)) +
  geom_point(size=.5) + 
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "grey") +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "grey") +
  geom_vline(xintercept = 120, color='red') +
  theme_tufte(11) +
  xlab("March") +
  ylab("Cluster") +
  annotate("text", x=12, y=3, label = "Hourly", size=2.5) +
  annotate("text", x=12, y=2.8, label = "NN=1", size=2.5) +
  scale_x_continuous(breaks = c(0,  24,  48,  72,  96, 120, 144, 168, 192, 216, 240, 264),
                     labels = c('10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21')) +
  scale_y_continuous(breaks = seq(1, length(unique(meddat2$group)), 1))
p3
p4 <- ggplot(meddat4, aes(x=value, y=group)) +
  geom_point(size=.5) + 
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "grey") +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "grey") +
  geom_vline(xintercept = 120, color='red') +
  theme_tufte(11) +
  xlab("March") +
  ylab("Cluster") +
  annotate("text", x=12, y=3, label = "Hourly", size=2.5) +
  annotate("text", x=12, y=2.8, label = "NN=5", size=2.5) +
  scale_x_continuous(breaks = c(0,  24,  48,  72,  96, 120, 144, 168, 192, 216, 240, 264),
                     labels = c('10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21')) +
  scale_y_continuous(breaks = seq(1, length(unique(meddat2$group)), 1))
p4
plot_grid(p3, p4)

plot_grid(p1, p2, p3, p4, ncol=2)
ggsave("~/Projects/Patagonia-EDA/figures/clustering_results_March10-20.pdf", width=6, height=4)
plot_grid(p1, p2, p3, p4, ncol=2)
}

# Clustering March5-25
{
  # Medoids data outcomes
  meddat1 <- as.data.frame(read_feather('~/Projects/Patagonia-EDA/data/k_medoids_day_NN1_3-05-2016-3-25-2016.feather'))
  meddat2 <- as.data.frame(read_feather('~/Projects/Patagonia-EDA/data/k_medoids_day_NN5_3-05-2016-3-25-2016.feather'))
  meddat3 <- as.data.frame(read_feather('~/Projects/Patagonia-EDA/data/k_medoids_day_hour_NN1_3-05-2016-3-25-2016.feather'))
  meddat4 <- as.data.frame(read_feather('~/Projects/Patagonia-EDA/data/k_medoids_day_hour_NN5_3-05-2016-3-25-2016.feather'))
  meddat1$group <- meddat1$group + 1
  meddat2$group <- meddat2$group + 1
  meddat3$group <- meddat3$group + 1
  meddat4$group <- meddat4$group + 1
  
  meddat1$value <- meddat1$value + 5
  meddat2$value <- meddat2$value + 5
  
  
  p1 <- ggplot(meddat1, aes(x=value, y=group)) +
    geom_point(size=.5) + 
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "grey") +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "grey") +
    geom_vline(xintercept = 15, color='red') +
    theme_tufte(11) +
    xlab(NULL) +
    ylab("Cluster") +
    annotate("text", x=5.5, y=3, label = "Daily", size=2.5) +
    annotate("text", x=5.5, y=2.8, label = "NN=1", size=2.5) +
    scale_y_continuous(breaks = seq(1, length(unique(meddat1$group)))) +
    scale_x_continuous(breaks = (meddat1$value))
  p1
  p2 <- ggplot(meddat2, aes(x=value, y=group)) +
    geom_point(size=.5) + 
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "grey") +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "grey") +
    geom_vline(xintercept = 15, color='red') +
    theme_tufte(11) +
    xlab(NULL) +
    ylab("Cluster") +
    annotate("text", x=5.5, y=3, label = "Daily", size=2.5) +
    annotate("text", x=5.5, y=2.8, label = "NN=5", size=2.5) +
    scale_y_continuous(breaks = seq(1, length(unique(meddat1$group)))) +
    scale_x_continuous(breaks = (meddat1$value))
  p2
  plot_grid(p1, p2, ncol=2)
  
  p3 <- ggplot(meddat3, aes(x=value, y=group)) +
    geom_point(size=.5) + 
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "grey") +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "grey") +
    geom_vline(xintercept = 240, color='red') +
    theme_tufte(11) +
    xlab("March") +
    ylab("Cluster") +
    annotate("text", x=12, y=3, label = "Hourly", size=2.5) +
    annotate("text", x=12, y=2.8, label = "NN=1", size=2.5) +
    scale_x_continuous(breaks = seq(0, 21*24, 24),
                       labels = seq(5, 26, 1)) +
    scale_y_continuous(breaks = seq(1, length(unique(meddat2$group)), 1))
  p3
  p4 <- ggplot(meddat4, aes(x=value, y=group)) +
    geom_point(size=.5) + 
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "grey") +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "grey") +
    geom_vline(xintercept = 240, color='red') +
    theme_tufte(11) +
    xlab("March") +
    ylab("Cluster") +
    annotate("text", x=12, y=3, label = "Hourly", size=2.5) +
    annotate("text", x=12, y=2.8, label = "NN=5", size=2.5) +
    scale_x_continuous(breaks = seq(0, 21*24, 24),
                       labels = seq(5, 26, 1)) +
    scale_y_continuous(breaks = seq(1, length(unique(meddat2$group)), 1))
  p4
  plot_grid(p3, p4)
  
  plot_grid(p1, p2, p3, p4, ncol=2)
  sggsave("~/Projects/Patagonia-EDA/figures/clustering_results_March05-25.pdf", width=6, height=4)
}


s# ggplot(pdat) +
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
