setwd("~/chinook_male_age/map")

library(ggplot2); library(tidyverse); library(sf)
library(bcmaps); library(ggspatial); library(sp) 
library(ggrepel); library(cowplot); library(geodata)

bch <- st_transform(bcmaps::bc_bound_hres(), crs = 4326)
sites <- read.csv("../data/hatchery_locations.csv")
USA <- sf::st_as_sf(geodata::gadm(country = "USA", level = 0, path = "."))
fraser <- bcmaps::watercourses_5M() %>% 
  filter(!name_en %in% c("Lillooet River"))
rivers <- st_transform(st_read(dsn = "BC_WATER_LINES_500M"), crs = 4326) %>% 
  filter(FCODE %in% c("GB15300000", "GA24850000")) %>% 
  mutate(lon = as.numeric(gsub(",.*", "", gsub("c", "", gsub("*\\(", "", geometry )))),
         lat = as.numeric(gsub("[^0-9.-]","", str_sub(string = geometry, -18 , -2))))

punt <- rivers %>% 
  filter(FEAT_LEN > 1e3) %>% 
  filter(between(lon, -125.2, -124.6) & between(lat, 49.65, 49.77)) 

qual <- rivers %>% 
  filter(between(lon, -124.73, -124.5) & between(lat, 49.36, 49.42))

chil <- rivers %>% 
  filter(between(lon, -122.5, -121.3) & between(lat, 49.05, 49.15))

lakes <- st_transform(st_read(dsn = "FWA_LAKES_POLY"), crs = 4326) %>% 
  filter(GNSNM1 %in% c("Cultus Lake", "Chilliwack Lake", "Comox Lake"))

(b <- ggplot() +
  geom_sf(data = USA) + theme_bw() +
  geom_sf(data = bch, fill = "gray90") +
  theme(panel.background = element_rect(fill = alpha("skyblue", 1/10)),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid = element_blank()) +
  labs(x = NULL, y = NULL) +
  geom_sf(data = fraser, colour = "skyblue", linewidth = 1/2) +
  geom_sf(data =  punt, colour = "skyblue", linewidth = 1/4) +
  geom_sf(data =  qual, colour = "skyblue", linewidth = 1/4) +
  geom_sf(data =  chil, colour = "skyblue", linewidth = 1/4) +
  geom_sf(data = lakes, color = "skyblue", fill = "navyblue") +
  geom_point(data = sites, aes(x = Lon, y = Lat), size = 1) +
  geom_label_repel(data = sites, aes(x = Lon, y = Lat, label = Site),
                   size = 2, min.segment.length = 0,
                   box.padding = 0, nudge_y = 0.1,
                   point.padding = 1/10, segment.size = 0.2) +
  coord_sf(xlim = c(-121.5, -125.4), ylim = c(49, 50)))

ggsave("plots/map.tiff", dpi = 300, height = 6, width = 6)

# Inset ------------------------------------------------------------------------

# Download low-res country outlines.
ca <- map_data("world", "Canada")
us <- map_data("world", "USA") 

# Make the inset plot by itself. 
(ins <- ggplot() +
    geom_polygon(data = us, aes(x = long, y = lat, group = group),
                 fill = "grey85", colour = "black", linewidth = 1/8) +
    geom_polygon(data = ca, aes(x = long, y = lat, group = group),
                 fill = "grey95", colour = "black", linewidth = 1/8) +
    theme_void() +
    theme(panel.border = element_rect(colour = "black", 
                                      fill = NA, linewidth = 1/4),
          panel.background = element_rect(fill = alpha("skyblue", 1/4)))  +
    annotate("rect", fill = NA, colour = "black",
             linewidth = 1/2,
             xmin = -125.4, xmax = -121.5,
             ymin = 49, ymax = 50) +
    annotate(geom = "text", label = "British 
Columbia", y = 55, x = -123.8, size = 3/2) +
    annotate(geom = "text", label = "Pacific 
Ocean", y = 50, x = -140, size = 3/2) +
    # Important to maintain accurate proportions/orientations. 
    # Plot is cartesian otherwise and appears distorted.
    coord_map(ylim = c(60, 45),
              xlim = c(-120, -150)))


# Add inset --------------------------------------------------------------------

ggdraw(plot = b) +
  draw_plot({
    ins
  },
  x = 0.77,
  y = 0.35,
  width = 0.2,
  height = 0.5)

ggsave("plots/map_winset.tiff", dpi = 300, 
       width = 6, height = 6)
