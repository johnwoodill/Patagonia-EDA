library(tidyverse)
library(lubridate)
library(ggmap)
library(cowplot)
library(gridExtra)
library(ggthemes)

GAPI_Key <- file("~/Projects/Patagonia-EDA/Google_api_key.txt", "r")
GAPI_Key <- readLines(GAPI_Key)
register_google(key=GAPI_Key)

dat <- read_feather('~/Data/GFW_point/Patagonia_Shelf/complete/Puerto_Madryn_region3_2016-3-1_2016-3-31.feather')
dat$month <- month(dat$timestamp)
dat$day <- day(dat$timestamp)
dat$hour <- hour(dat$timestamp)
dat$ln_distance <- log(1 + dat$distance)


# Figure 1
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
#42.7694° S, 65.0317° W

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


# Figure 2

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
#ggsave("~/Projects/Patagonia-EDA/figures/Main Figures/figure2.pdf")

# Figure 3

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
