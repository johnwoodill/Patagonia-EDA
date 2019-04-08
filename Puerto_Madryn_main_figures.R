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
  scale_x_continuous(breaks = seq(1, 31, 1), expand=expand_scale(mult = c(0, 0))) +
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
JSD <- function(p, q){
  pfun <- approxfun(density(p))
  p <- pfun(p)/sum(pfun(p))
  qfun <- approxfun(density(q))
  q <- qfun(q)/sum(qfun(q))
  m <- 0.5 * (p + q)
  JS <- 0.5 * (sum(p * log(p / m)) + sum(q * log(q / m)))
  return(JS)
}
library(MASS)
dat <- as.data.frame(read_feather("~/Data/GFW_point/Patagonia_Shelf/Puerto_Madryn/dist_matrices/dmat_Puerto_Madryn_region1_NN1_day-hour_2016-03-01_2016-03-31.feather"))
dat <- dat[-1, -1]
d <- as.matrix(dat)
fit <- isoMDS(d, k=3)

isoMDS_dat <- data.frame(x = fit$points[,1], y = fit$points[,2])

# Hourly cluster
#Clusters: 

cluster1 <- c(43, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 
              28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 
              55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80,
              81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 
              106, 107, 108, 109, 110, 116, 118, 122, 123)
cluster2 <- c(209, 111, 112, 113, 114, 115, 117, 119, 120, 121, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 
  138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 
  162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 
  186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 210, 211,
  212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 
  236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 
  260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 
  284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 
  308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 
  332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 
  356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 
  380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 
  404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 
  428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 
  452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 
  476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491)
cluster3 <- c(673, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 
              513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 
              535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 
              557, 558, 559, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569, 570, 571, 572, 573, 574, 575, 576, 577, 578, 
              579, 580, 581, 582, 583, 584, 585, 586, 587, 588, 589, 590, 591, 592, 593, 594, 595, 596, 597, 598, 599, 600,
              601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622, 
              623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 
              645, 646, 647, 648, 649, 650, 651, 652, 653, 654, 655, 656, 657, 658, 659, 660, 661, 662, 663, 664, 665, 666, 
              667, 668, 669, 670, 671, 672, 674, 675, 676, 677, 678, 679, 680, 681, 682, 683, 684, 685, 686, 687, 688, 689, 
              690, 691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 
              712, 713, 714, 715, 716, 717, 718, 719, 720, 721, 722, 723, 724, 725, 726, 727, 728, 729, 730, 731, 732, 733, 
              734, 735, 736, 737, 738, 739, 740, 741, 742, 743)
isoMDS_dat$row <- seq(1, nrow(isoMDS_dat))

clustdat <- data.frame(cluster = c(rep(1, length(cluster1)),
                                   rep(2, length(cluster2)),
                                   rep(3, length(cluster3))),
                       row = c(cluster1, cluster2, cluster3))
isoMDS_dat <- left_join(isoMDS_dat, clustdat, by = 'row')

isoMDS_dat$dist <- sqrt((isoMDS_dat$x - isoMDS_dat$y)^2 + (lead(isoMDS_dat$x) - lead(isoMDS_dat$y))^2)
library(RColorBrewer)
ggplot(isoMDS_dat, aes(x, y, color = dist)) + geom_point() +
  theme_tufte(11) +
  # ylim(-0.12, .12) +
  # xlim(-0.12, 0.12) +
  scale_color_gradientn(colours=brewer.pal(7,"YlGnBu")) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "grey") +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=-Inf, y=Inf, yend=Inf, color = "grey") +
  labs(x = "Dimension 1", y= "Dimension 2", color = "Cluster", title = "Non-metric Multidimensional Scaling of JS-Distance \n March 1 - 31 2016 (Clustered by hour using k-medoids)") +
  theme(legend.position = c(.5, .90),
        legend.direction = 'horizontal',
        legend.justification = 'center',
        legend.title=element_text(size=8),
        legend.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -0.08, y = -0.07, label = "stress=6.91", size=2.5)
  NULL
ggsave("~/Projects/Patagonia-EDA/figures/Main Figures/NMDS_k=2_by_cluster.pdf", width = 6, height = 4)

# Color by speed (distance from day to lead(day))
ggplot(isoMDS_dat, aes(x, y, color = dist)) + geom_point() +
  theme_tufte(11) +
  # ylim(-0.12, .12) +
  # xlim(-0.12, 0.12) +
  scale_color_gradientn(colours=brewer.pal(7,"YlGnBu")) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "grey") +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=-Inf, y=Inf, yend=Inf, color = "grey") +
  labs(x = "Dimension 1", y= "Dimension 2", color = "Speed", title = "Non-metric Multidimensional Scaling of JS-Distance \n March 1 - 31 2016 (Clustered by hour using k-medoids)") +
  theme(legend.position = c(.5, .90),
        legend.direction = 'horizontal',
        legend.justification = 'center',
        legend.title=element_text(size=8),
        legend.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = -0.08, y = -0.06, label = "k=2", size=2.5) + 
  annotate("text", x = -0.08, y = -0.07, label = "stress=6.91", size=2.5)
NULL
ggsave("~/Projects/Patagonia-EDA/figures/Main Figures/NMDS_k=2_by_speed.pdf", width = 6, height = 4)


# 3-D plot
library(plotly)

isoMDS_dat <- data.frame(x = fit$points[,1], y = fit$points[,2], z = fit$points[, 3])
isoMDS_dat$row <- seq(1, nrow(isoMDS_dat))
isoMDS_dat <- left_join(isoMDS_dat, clustdat, by = 'row')


day15 <- isoMDS_dat[336:360, ]
day15
isoMDS_dat$event <- ifelse(isoMDS_dat$row >= 336, ifelse(isoMDS_dat$row <= 360, 1, 0), 0)
ggplot(isoMDS_dat, aes(x=x, y=y, z=z, color=factor(cluster))) + 
  # scale_shape_manual(values = c(1, 16)) +
  labs(x = "Dimension 1", y= "Dimension 2", color = "Cluster", title = "Non-metric Multidimensional Scaling of JS-Distance \n March 1 - 31 2016 (Clustered by hour using k-medoids)") +
  theme_void() +
  theme_tufte(11) +
  axes_3D() +
  stat_3D() +
  theme(legend.position = c(.75, .90),
        legend.direction = 'horizontal',
        legend.justification = 'center',
        legend.title=element_text(size=8),
        legend.text = element_text(size=8),
        plot.title = element_text(hjust = 0.5)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "grey") +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=-Inf, y=Inf, yend=Inf, color = "grey") +
  annotate("text", x = -.25, y = -0.28, label = "stress=4.86", size=2.5) +
  annotate("text", x = -.25, y = -0.25, label = "k=3", size=2.5) +

NULL
ggsave("~/Projects/Patagonia-EDA/figures/Main Figures/NMDS_k=3_by_cluster.pdf", width = 6, height = 4)


# Speed plot from above
isoMDS_dat$day <- (isoMDS_dat$row/24)
isoMDS_dat$day <- lead(isoMDS_dat$day) + 1
ggplot(isoMDS_dat, aes(x=day, y=dist)) + 
  labs(x="March", y="Hourly JS-Divergence Distance") +
  theme_tufte(11) +
  geom_line() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, color = "grey") +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=Inf, y=-Inf, yend=Inf, color = "grey") +
  annotate("segment", x=Inf, xend=-Inf, y=Inf, yend=Inf, color = "grey") +
  scale_x_continuous(breaks = seq(1, 31, 1), labels = seq(1, 31, 1)) +
  # geom_vline(xintercept = 1)
  NULL
ggsave("~/Projects/Patagonia-EDA/figures/Main Figures/Divergence_distance_vplot.pdf", width=6, height=4)


########### END OF MAIN FIGURES
dat2 <- read_feather('~/Data/GFW_point/Patagonia_Shelf/complete/Puerto_Madryn_region1_2016-3-1_2016-3-31.feather')
dat$index <- NULL
names(dat)[1] <- "day"

dat2 <- filter(dat2, rank == 1)
dat2 <- filter(dat2, distance != 0)
dat2$day <- day(dat2$timestamp)
dat2$hour <- hour(dat2$timestamp) + 1
dat2$month <- month(dat2$timestamp)
dat2$day_hour <- dat2$day * dat2$hour
dat2 <- filter(dat2, month == 3)

jsd_dat <- data.frame()
for (i in unique(dat2$timestamp)){
  p <- filter(dat2, timestamp == i )
  q <- filter(dat2, timestamp == i + 60*60)
  date = p$timestamp[1]
  # p <- filter(dat2, day_hour == i)
  # q <- filter(dat2, day_hour == i + 1)
  p <- p$distance
  q <- q$distance
  qmin <- length(q)
  pmin <- length(p)
  minobs <- min(qmin, pmin)
  lst <- data.frame()
  for (j in 1:10){
    pp <- sample(p, minobs, replace = TRUE)
    qq <- sample(q, minobs, replace = TRUE)
    js <- JSD(qq, pp)
    indat <- data.frame(js = js)
    lst <- rbind(lst, indat)
  }
  outdat <- data.frame(t = date, jsd_mean = mean(lst$js), jsd_sd = sd(lst$js))
  jsd_dat <- rbind(jsd_dat, outdat)
}

jsd_dat$day <- day(jsd_dat$t)
jsd_dat2 <- jsd_dat %>% 
  group_by(day) %>% 
  summarise(jsd_mean = mean(jsd_mean))
  

ggplot(jsd_dat2, aes(x=day, y=jsd_mean)) + 
  geom_point() +
  #geom_point(aes(x=t, y=jsd_mean + 1.96*jsd_sd)) +
  ylim(0, 0.15) +
  geom_vline(xintercept = 6) +
  geom_vline(xintercept = 20)

geom_point(aes(x=t, y=jsd_mean - 1.96*jsd_sd))

jsd_dat <- data.frame()
for (i in 10:20){
  
  q <- dat[, i+1]
  p <- dat[, i]
  r <- dat[, i-1]
  js1 <- sqrt(JSD(q, p))
  js2 <- sqrt(JSD(r, p))
  outdat <- data.frame(t = i, jsd1 = js1, jsd2 = js2)
  jsd_dat <- rbind(jsd_dat, outdat)
}


qplot(x=10:20, y=(jsd_dat$jsd2 - jsd_dat$jsd1))

plot(x=10:20, y=(jsd_dat$jsd - lag(jsd_dat$jsd))/(lag(jsd_dat$jsd)))
plot(jsd_dat$t[2:31], jsd_dat$jsd[2:31])

p <- dat$`2`
q <- dat$`1`

JSD(p, q)

pdat <- gather(dat, key = day, value = value, -day)
pdat$day <- as.numeric(pdat$day)
pdat$day2 <- seq(1, length(unique(pdat$day)))

pdat2 <- pdat %>% group_by(day) %>% 
  summarise(mvalue = mean(value)) %>% 
  filter(day > 0) %>% 
  ungroup() %>% 
  mutate(change = (mvalue - lag(mvalue))/lag(mvalue))

plot(x=pdat2$day, y=pdat2$mvalue)
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
