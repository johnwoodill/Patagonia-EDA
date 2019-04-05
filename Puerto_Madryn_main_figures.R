library(tidyverse)
library(lubridate)
library(ggmap)
library(cowplot)
library(gridExtra)
library(ggthemes)
library(feather)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

GAPI_Key <- file("~/Projects/Patagonia-EDA/Google_api_key.txt", "r")
GAPI_Key <- readLines(GAPI_Key)
register_google(key=GAPI_Key)

dat <- read_feather('~/Data/GFW_point/Patagonia_Shelf/complete/Puerto_Madryn_region3_2016-3-1_2016-3-31.feather')
dat$month <- month(dat$timestamp)
dat$day <- day(dat$timestamp)
dat$hour <- hour(dat$timestamp)
dat$ln_distance <- log(1 + dat$distance)


# Figure 1a
fig1_dat <- dat %>% 
  filter(day == 15) %>% 
  group_by(vessel_A, month, day, hour) %>% 
  summarise(vessel_A_lat = mean(vessel_A_lat),
            vessel_A_lon = mean(vessel_A_lon))
fig1_dat  

lat <- mean(fig1_dat$vessel_A_lat)
lon <- mean(fig1_dat$vessel_A_lon)

# Patagonia shelf
lon1 = -63
lon2 = -55
lat1 = -45
lat2 = -40

# Puerto Madryn
#-42.7694° S, -65.0317° W

myLocation <- c(left=lon1 - 5, bottom=lat1 - 5, right=lon2 + 5, top=lat2 + 5)

p1 <- ggmap(get_stamenmap(myLocation, zoom=5)) +
  theme_tufte(11) +
  theme(legend.position = 'none') +
  annotate("text", x=-53.5, y = -49.5, label="March 15, 2016", size = 4) +
  annotate("text", x=-65.0317, y = -42.1, label="Puerto Madryn", size = 3) +
  geom_point(aes(x=-65.0317, y=-42.7694), size = .5) +
  # coord_fixed(ratio = 2) +
  xlab("Longitude") +
  ylab("Latitude") +
  geom_path(aes(x=vessel_A_lon, y=vessel_A_lat, color=factor(vessel_A)), size = .5, alpha=.5, data = fig1_dat) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "grey") +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=-Inf, y=Inf, yend=Inf, color = "grey")

p1
myMap <- get_stamenmap(location=myLocation,
          source="stamen", color="bw")

myMap
ggmap(myMap, source="stamen", color="bw", zoom=5) +
  theme(legend.position = 'none') +
  geom_path(aes(x=vessel_A_lon, y=vessel_A_lat, color=factor(vessel_A)), size = .5, alpha=.5, data = fig1_dat) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "grey") +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=-Inf, y=Inf, yend=Inf, color = "grey")
  
p1
stamen: terrain


# Figure 1b

fig2_dat <- dat %>% 
  group_by(rank) %>% 
  summarise(ln_distance = mean(ln_distance)) %>% 
  filter(rank > 0) %>% 
  ungroup()

p2 <- ggplot(fig2_dat, aes(x=rank, y=ln_distance)) +
  # geom_bar(stat='identity', width=0.1, space=1) +
  geom_point(space=1) +
  ylab("Log(NN Distance)") + xlab("Nearest Neighbor (NN)") +
  ylim(0, 3) + 
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "grey") +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=-Inf, y=Inf, yend=Inf, color = "grey") +
  geom_linerange(aes(ymin=0, ymax = ln_distance)) +
  theme_tufte(11) +
  scale_x_continuous(breaks = seq(1, 10, 1))
p2
#ggsave("~/Projects/Patagonia-EDA/figures/Main Figures/figure2.pdf")

# Figure 1c

fig3_dat <- dat %>%
  filter(rank <= 5) %>% 
  group_by(day, hour) %>% 
  summarise(mean_distance = mean(ln_distance))

p3 <- ggplot(fig3_dat, aes(mean_distance)) + 
  geom_density(aes(y=..density../nrow(fig3_dat))) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "grey") +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=-Inf, y=Inf, yend=Inf, color = "grey") +
  theme_tufte(11) +
  ylab("P(Nearest Neighbor Distances)") +
  xlab("Mean Distance (NN=5)")
p3


# Panel plot
figure1 <- ggdraw() + draw_plot(p1, .0, .00, height = 1, width = .55)
figure1 <- figure1 + draw_plot(p2, .55, .5, height = .5, width = .45)
figure1 + draw_plot(p3, .55, .0, height = .5, width = .45)

ggsave(filename = "~/Projects/Patagonia-EDA/figures/Main Figures/figure1.pdf", width = 7, height = 5)

# Figure 2

# (A)
{
dat <- read_feather('~/Data/GFW_point/Patagonia_Shelf/complete/Puerto_Madryn_region1_2016-03-10_2016-03-20.feather')
dat <- filter(dat, distance != 0)
dat <- filter(dat, rank <= 5)

dat$date <- paste0(year(dat$timestamp), "-", month(dat$timestamp), "-", day(dat$timestamp), "-", hour(dat$timestamp))

dat <- dat %>% 
  mutate(day = day(timestamp),
         hour = hour(timestamp)) %>% 
  group_by(timestamp, date, vessel_A) %>% 
  summarise(distance = mean(distance)) %>% 
  ungroup() %>% 
  mutate(ln_distance = log(1 + distance))
head(dat)  

breaks <- 50
outdat <- data.frame()
for (i in unique(dat$timestamp)){
  indat <- filter(dat, timestamp == i)
  efunc <- ecdf(indat$ln_distance)
  cdat <- cut(indat$ln_distance, breaks = breaks)
  
  indat$breaks <- cdat 
  indat$bin <- as.numeric(indat$breaks)
  
  sumdat <- indat %>% 
    group_by(bin) %>% 
    summarise(mean_ln_distance = mean(ln_distance))
  
  char <- unique(cdat)
  retdat <- data.frame()
  for (j in char){
    ldat <- data.frame(a = as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", j), ",")[[1]][1]),
                       b = as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", j), ",")[[1]][2]))
    retdat <- rbind(retdat, ldat)
  }
  
  retdat <- arrange(retdat, a, b)
  retdat$bin <- seq(1, length(unique(cdat)))
  retdat$prob <- efunc(retdat$b) - efunc(retdat$a)
  retdat$a <- NULL; retdat$b <- NULL
  retdat$timestamp <- indat$timestamp[1]
  retdat <- left_join(retdat, sumdat, by = 'bin')
  #retdat$date <- indat$date[1]
  outdat <- rbind(outdat, retdat) 
  print(i)
}
  
outdat %>% 
  group_by(timestamp) %>% 
  summarise(sum = sum(prob))

outdat <- filter(outdat, bin >= 5)
outdat$day <- day(outdat$timestamp)


outdat$timestamp
as.Date.POSIXct(outdat$timestamp, c("%Y-%m-%d %h:00:00 PST"))
outdat$timestamp <- as.POSIXct(outdat$timestamp)

#as.Date(outdat$timestamp, c("%Y-%m-%d %h:00:00 PST"))

#outdat$timestamp <- as.Date(outdat$timestamp)

fig2a <- ggplot(outdat, aes(timestamp, factor(bin))) + 
  geom_raster(aes(fill = prob)) +
  theme_tufte(11) +
  annotate("text", x=as.POSIXct("2016-03-20 6:00:00"), y = 44, label='(a)', size = 5) +
  labs(x="March", y="Binned Distance", fill="P(d)") +
  scale_x_datetime(date_breaks = "1 day", 
                   date_labels = "%d", 
                   labels = c(seq(10, 21, 1)), 
                   breaks = seq(10, 21, 1), 
                   expand = expand_scale(mult = c(0, 0))) +
  scale_y_discrete(breaks = floor(seq(1, 50, length.out=11)), 
                   labels = floor(seq(1, 50, length.out=11)), 
                   expand = c(0, 0)) +
  scale_fill_gradient(low = "white", high = "red", 
                      breaks = c(0, 0.025, 0.050, 0.075, 0.1), 
                      labels = c('0', '0.025', '0.050', '0.075', '0.1'),
                      na.value="transparent",
                      lim = c(0, 0.1)) +
  theme(legend.margin=margin(l = 0, unit='cm'),
        panel.border = element_rect(colour = "grey", fill=NA, size=1),
        panel.grid = element_blank()
        ) +
  guides(fill = guide_colorbar(label.hjust = unit(0, 'cm'),
                               frame.colour = "black",
                               barwidth = .5,
                               barheight = 16))

fig2a

ggsave("~/Projects/Patagonia-EDA/figures/Main Figures/figure2a.png", width = 6, height = 4)
ggsave("~/Projects/Patagonia-EDA/figures/Main Figures/figure2a.pdf", width = 6, height = 4)

}


# (B)

# Figure 2 (B) 

dat <- as.data.frame(read_feather("/home/server/pi/homes/woodilla/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region1_NN1_day-hour_2016-03-01_2016-03-31.feather"))
dat$index <- NULL
names(dat)[1] <- "day"


pdat <- gather(dat, key = day, value = value, -day)
pdat$day <- as.numeric(pdat$day)
pdat$day2 <- seq(1, length(unique(pdat$day)))

fig2b <- ggplot(pdat, aes(x=day, y=day2)) + 
  theme_tufte(11) + 
  geom_tile(aes(fill=value)) +
  labs(y="March", x="March", fill="J-S Metric") +
  #geom_raster(aes(fill=value), interpolate = TRUE) +
  scale_fill_gradient(low='white', high='red') +
  scale_y_continuous(trans = "reverse", expand=expand_scale(mult = c(0, 0))) +
  scale_x_continuous(expand=expand_scale(mult = c(0, 0))) +
  # geom_vline(xintercept = 6) +
  # geom_vline(xintercept = 20) +
  # geom_hline(yintercept = 6) +
  # geom_hline(yintercept = 20) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "black") +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "black") +
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, color = "black") +
  annotate("segment", x=Inf, xend=-Inf, y=Inf, yend=Inf, color = "black") +
  
  # Cluster 1
  annotate("segment", x=6, xend=6, y=0, yend=6, color = "black") +
  annotate("segment", x=0, xend=6, y=6, yend=6, color = "black") +
  
  # Cluster 2
  annotate("segment", x=6, xend=20, y=6, yend=6, color = "black") +
  annotate("segment", x=6, xend=6, y=6, yend=20, color = "black") +
  annotate("segment", x=20, xend=20, y=6, yend=20, color = "black") +
  annotate("segment", x=6, xend=20, y=20, yend=20, color = "black") +
  
  # Cluster 3
  annotate("segment", x=20, xend=20, y=20, yend=32, color = "black") +
  annotate("segment", x=20, xend=32, y=20, yend=20, color = "black") +
  
  annotate("text", x=4, y = 1, label='Cluster 1', size = 3) +
  annotate("text", x=18, y = 7, label='Cluster 2', size = 3) +
  annotate("text", x=29, y = 21, label='Cluster 3', size = 3) +
  annotate("text", x=30.5, y = 1.5, label='(b)', size = 5) +
  guides(fill = guide_colorbar(label.hjust = unit(0, 'cm'),
                               frame.colour = "black",
                               barwidth = .5,
                               barheight = 16)) +
  NULL
fig2b
ggsave("~/Projects/Patagonia-EDA/figures/Main Figures/figure2b.png", width=6, height = 4)
ggsave("~/Projects/Patagonia-EDA/figures/Main Figures/figure2b.pdf", width=6, height = 4)

plot_grid(fig2a, fig2b, ncol=2)
ggsave("~/Projects/Patagonia-EDA/figures/Main Figures/figure2.pdf", width = 18, height = 6)


# Figure 3
{
#dat <- as.data.frame(read_feather("/home/server/pi/homes/woodilla/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region1_NN1_day-hour_2016-03-01_2016-03-31.feather"))
dat <- read_feather('~/Data/GFW_point/Patagonia_Shelf/complete/Puerto_Madryn_region1_2016-3-1_2016-3-31.feather')
dat <- filter(dat, distance != 0)
dat <- filter(dat, rank <= 5)
dat$ln_distance <- log(1 + dat$distance)
dat$month <- month(dat$timestamp)
dat$day <- day(dat$timestamp)
dat$hour <- hour(dat$timestamp)

dat <- filter(dat, month == 3)

fig3_dat <- dat %>% 
  group_by(day, hour, vessel_A) %>% 
  summarise(medoid = as.numeric(pam(ln_distance, 1)[1]),
            mean = mean(ln_distance),
            med = median(ln_distance),
            mode = Mode(ln_distance))
fig3_dat <- ungroup(fig3_dat)
pdat <- gather(fig3_dat, key = moment, value=value, -day, -hour, -vessel_A)

ggplot(pdat, aes(factor(day), value, color=moment)) + geom_boxplot() + facet_wrap(~moment, scale='free')

head(fig3_dat)

dat$day <- day(dat$timestamp)

dat <- dat %>% 
  group_by(timestamp, vessel_A) %>% 
  summarise(ln_distance = mean(ln_distance),
            med_ln_distance = median(ln_distance),
            ) %>% 
  ungroup()
dat
# dat$index <- NULL
# rownames(dat) <- marg_dat$`0`
# dat
# dat$x <- seq(1, 31, 1)
# 
# 
# dat <- gather(dat, key = x, value = value)
# dat <- filter(dat, x > 0)
# dat$y <- seq(1, 31,1 )
# dat$x <- as.numeric(dat$x)

# fig3_dat <- dat %>% 
#   group_by(x) %>% 
#   summarise(mean = mean(value),
#             q25 = quantile(value)[2],
#             q75 = quantile(value)[4])
dat$day <- day(dat$timestamp)


fig3_dat <- fig3_dat %>% 
  group_by(day) %>% 
  summarise(mean = mean(med),
            med = median(mean))
,
            q25 = quantile(ln_distance)[2],
            q75 = quantile(ln_distance)[4])

fig3_dat$cluster <- 3
fig3_dat$cluster <- ifelse(fig3_dat$day <= 6, 1, fig3_dat$cluster)
fig3_dat$cluster <- ifelse(fig3_dat$day > 6 & fig3_dat$day < 20, 2, fig3_dat$cluster)
fig3_dat

fig3 <- ggplot(fig3_dat, aes(day, med)) + 
  theme_tufte(11) +
  # ylim(0, 3) +
  labs(y = "Log(Average Distance) \n (NN=5)", x="March", shape="Cluster") +
  # geom_linerange(aes(ymin=q25, ymax = q75), color = 'grey') +
  # geom_point(aes(day, q25), size = .25) +
  # geom_point(aes(day, q75), size = .25) +
  geom_point(aes(shape = factor(cluster))) +
  # geom_line() +
  scale_x_continuous(breaks = seq(1, 31), labels = seq(1, 31)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "grey") +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=-Inf, y=Inf, yend=Inf, color = "grey") +
  theme(legend.position = c(.5, .90),
        legend.direction = 'horizontal',
        legend.justification = 'center',
        legend.title=element_text(size=8),
        legend.text = element_text(size=8)) +
  guides(shape = guide_legend(title.position="top", title.hjust = 0.5))

fig3

# k-medoids cluster analysis
#[0:123], [124:491], [492, 673]

fig3_dat$cluster <- 3
fig3_dat$cluster <- ifelse(fig3_dat$x <= 6, 1, fig3_dat$cluster)
fig3_dat$cluster <- ifelse(fig3_dat$x > 6 & fig3_dat$x < 20, 2, fig3_dat$cluster)

dens1 <- filter(dat, day <= 6)
dens2 <- filter(dat, day >= 6 & day <= 20)
dens3 <- filter(dat, day > 20)

dp1 <- ggplot(dens1, aes(ln_distance)) +  
  ylim(0, 0.55) +
  geom_density(color = 'darkgrey') +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "grey") +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=-Inf, y=Inf, yend=Inf, color = "grey") +
  annotate("text", x=.5, y=0.5, label='Cluster 1', size=2.5) +
  theme_tufte(11) +
  xlab(NULL) +
  ylab(NULL)
dp1

dp2 <- ggplot(dens2, aes(ln_distance)) +  
  ylim(0, 0.55) +
  geom_density(color = 'darkgrey') +
  theme_tufte(11) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "grey") +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=-Inf, y=Inf, yend=Inf, color = "grey") +
  annotate("text", x=0.5, y=0.5, label='Cluster 2', size=2.5) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  xlab("Log(Distance)") +
  ylab(NULL)
dp2

dp3 <- ggplot(dens3, aes(ln_distance)) +  
  ylim(0, 0.55) +
  geom_density(color = 'darkgrey') +
  theme_tufte(11) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "grey") +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=-Inf, y=Inf, yend=Inf, color = "grey") +
  annotate("text", x=0.3, y=0.5, label='Cluster 3', size=2.5) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  xlab(NULL) +
  ylab(NULL)

dp3


# Panel plot
{
  ggdraw() + draw_plot(fig3, 0 ,.35, height = .65, width = 1) +
    draw_plot(dp1, .05, .043, height=.305, width = .3) + 
    draw_plot(dp2, .35, 0, height=.35, width = .3) +
    draw_plot(dp3, .65, .043, height=.305, width = .3) 
    
ggsave("~/Projects/Patagonia-EDA/figures/Main Figures/figure3.png", width = 8, height = 4)
ggsave("~/Projects/Patagonia-EDA/figures/Main Figures/figure3.pdf", width = 8, height = 4)
}

# Figure 4 - time-series of the trailing rate of change in JS divergence as an index of behavioral change

dat <- as.data.frame(read_feather("/home/server/pi/homes/woodilla/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region1_NN1_day-hour_2016-03-01_2016-03-31.feather"))
dat$index <- NULL
names(dat)[1] <- "day"


pdat <- gather(dat, key = day, value = value, -day)
pdat$day <- as.numeric(pdat$day)
pdat$day2 <- seq(1, length(unique(pdat$day)))
pdat <- filter(pdat, day <= 31 & day2 <= 31)
pdat <- filter(pdat, day == 31)
pdat$value <- 1/pdat$value
plot(x=pdat$day2, y=pdat$value)

pdat <- arrange(pdat, -day, -day2)
pdat$effect <- (lead(pdat$value) - pdat$value)/(lead(pdat$day2) - pdat$day2)

plot(pdat$day2, pdat$effect)

p1 <- ggplot(pdat, aes(x=day2, y=value)) + 
  geom_line(color='grey') +
  geom_point(size=.5) +
  theme_tufte(11) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "grey") +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=-Inf, y=Inf, yend=Inf, color = "grey") +
  labs(y="1 / JS-Metric", x="March") +
  scale_x_continuous(breaks = seq(1, 31, 1)) +
  NULL
  
  
p2 <- ggplot(pdat, aes(x=day2, y=effect)) + 
  geom_line(color='grey') +
  geom_point(size=.5) +
  # ylim(0, 3) +
  theme_tufte(11) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "grey") +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=-Inf, y=Inf, yend=Inf, color = "grey") +
  labs(y="JS-Metric Rate of Change", x="March") +
  scale_x_continuous(breaks = seq(2, 31, 1), labels = seq(2, 31, 1)) +
NULL

plot_grid(p1, p2, ncol=1)
ggsave("~/Projects/Patagonia-EDA/figures/Main Figures/figure4.pdf", width=6, height=4)
