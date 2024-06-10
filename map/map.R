setwd("~/chinook_male_age/map")

library(ggplot2); library(tidyverse); library(sf)
library(bcmaps); library(ggspatial); library(sp) 
library(ggrepel); library(cowplot); library(geodata)
library(ggsci)

################################################################################
# To do:
# Change x-axis from Qualicum River to Big Qualicum River
# Choose a more accessible colour palette 
################################################################################


# High-resolution outline of BC. Convert to common coordinates.
bch <- st_transform(bcmaps::bc_bound_hres(), crs = 4326)

# Read in site information. 
sites <- read.csv("../data/hatchery_locations.csv")

# Shape file of USA. 
USA <- sf::st_as_sf(geodata::gadm(country = "USA", level = 0, path = "."))

# Relatively low-resolution data of BC watercourses. 
# At this scale, basically just Fraser and a few others.
# Remove Lillooet so it doesn't bleed through inset plot later on. 
fraser <- bcmaps::watercourses_5M() %>% 
  filter(!name_en %in% c("Lillooet River"))

# Downloaded high-res shapefiles of water lines via iMap BC.
rivers <- st_transform(st_read(dsn = "BC_WATER_LINES_500M"), crs = 4326) %>% 
  # Filter out coastlines other small features. 
  filter(FCODE %in% c("GB15300000", "GA24850000")) %>% 
  # Awkwardly create lat/lon features for some filtering.
  # Just prevents a super busy and messy figure later on.
  mutate(lon = as.numeric(gsub(",.*", "", gsub("c", "", gsub("*\\(", "", geometry )))),
         lat = as.numeric(gsub("[^0-9.-]","", str_sub(string = geometry, -18 , -2))))

# Isolate waterways near Puntledge River site.
punt <- rivers %>% 
  filter(FEAT_LEN > 1e3) %>% 
  filter(between(lon, -125.2, -124.6) & between(lat, 49.65, 49.77)) 

# Isolate waterways near Qualicum River site.
qual <- rivers %>% 
  filter(between(lon, -124.73, -124.5) & between(lat, 49.36, 49.42))

# Isolate waterways near Chilliwack site.
chil <- rivers %>% 
  filter(between(lon, -122.5, -121.3) & between(lat, 49.05, 49.15))

# Polygon shape files of various lakes in BC courtesy of iMap BC download.
# Only include those directly relevant to above sample sites. 
lakes <- st_transform(st_read(dsn = "FWA_LAKES_POLY"), crs = 4326) %>% 
  filter(GNSNM1 %in% c("Cultus Lake", "Chilliwack Lake", "Comox Lake"))

# Plot sample sites on map. 
# A lot of trial and error in placing objects on plot below. 
(b <- ggplot() +
  geom_sf(data = USA) + theme_bw() +
  geom_sf(data = bch, fill = "gray90") +
  theme(panel.background = element_rect(fill = alpha("skyblue", 1/10)),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid = element_blank()) +
  labs(x = NULL, y = NULL) +
  geom_sf(data = fraser, colour = "skyblue", linewidth = 1/4) +
  geom_sf(data =  punt, colour  = "skyblue", linewidth = 1/4) +
  geom_sf(data =  qual, colour  = "skyblue", linewidth = 1/4) +
  geom_sf(data =  chil, colour  = "skyblue", linewidth = 1/4) +
  geom_sf(data = lakes, colour  = "skyblue", fill = "skyblue2") +
  geom_point(data = sites, aes(x = Lon, y = Lat), size = 1) +
  geom_label_repel(data = sites, aes(x = Lon, y = Lat, label = Site),
                   size = 2, min.segment.length = 0,
                   box.padding = 0, nudge_y = 0.1,
                   point.padding = 1/10, segment.size = 0.2) +
  ggspatial::annotation_scale(location = "tr",
                              width_hint = 1/10,
                              pad_x = unit(0.32, "cm"),
                              pad_y = unit(3.20, "cm")) +
  # Add custom N arrow because the ggspatial options are all very nautical.
  geom_segment(aes(x = -122.4, xend = -122.4,
                     y = 49.60, yend = 50),
                 arrow = arrow(length = unit(1/5, "cm"))) +
  annotate("text", label = "N", x = -122.4, y = 49.55) +
  annotate("text", label = "Straight of Georgia", 
             x = -123.7, y = 49.3, size = 3) +
  coord_sf(xlim = c(-121.5, -125.4), ylim = c(49, 50)))


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
    # Awkward parsing below to ensure text is centered over two lines.
    annotate(geom = "text", label = "British 
Columbia", y = 55, x = -123.8, size = 3/2) +
    annotate(geom = "text", label = "Pacific 
Ocean", y = 50, x = -140, size = 3/2) +
    geom_segment(aes(x = -135, xend = -126,
                     y = 46, yend = 48.8),
                 arrow = arrow(length = unit(1/5, "cm"))) +
    # Important to maintain accurate proportions/orientations. 
    # Plot is Cartesian otherwise and appears distorted.
    coord_map(ylim = c(60, 45),
              xlim = c(-120, -150)))


# Add inset --------------------------------------------------------------------


# Use cowplot to insert inset plot of NA. 
(mapwinset <- ggdraw(plot = b) +
  draw_plot({
    ins
  },
  x = 0.77,
  y = 0.50,
  width = 0.2,
  height = 0.5))


# Sample info ------------------------------------------------------------------


# Read in sample information.
dat <- read.csv("../data/chRADseq_samples.csv") %>% 
  select(c("fishID", "pop", "year", "age")) %>% 
  mutate(pop = as.factor(str_to_title(gsub("_", " ", pop)))) 

# Count samples by grouping factors and count.
tab <- dat %>% 
  group_by(pop, age, year) %>% 
  tally() %>% 
  mutate(year = factor(year,
         levels = c("2016", "2017", "2018", "2019", "2020", "2021")))

# Plot barplots of sample composition by site. 
(samples <- ggplot(data = tab %>% 
         mutate(year = as.factor(year))) +
  geom_bar(aes(x = age, y = n, fill = year), 
           color = "black", stat = "identity",
           linewidth = 1/5) +
  scale_fill_npg(palette = c("nrc")) +
  facet_grid(~ pop, scales = "free_x", space = "free_x", switch = "x") +
  theme_bw() + labs(x = NULL, y = "Samples") +
  theme(strip.placement = "outside", 
        strip.text = element_text(size = 12, vjust = 2.5),
        strip.background = element_rect(fill = NA, color = NA),
        legend.position = "right", legend.title = element_blank(),
        legend.margin = unit(-2, "cm")) +
  guides(fill = guide_legend(ncol = 1, byrow = T)))

# Join map and sample info figures together.
ggdraw(cowplot::plot_grid(mapwinset, samples, ncol = 1, 
                align = "VH", rel_heights = c(1, 0.75),
                labels = c('a)', 'b)'),
                label_size = 10, label_y = 1.02)) +
  draw_label("Sea age:", x = 0.025, y = 0.095, size = 8, hjust = 0)

ggsave("../plots/map_wsamples.tiff", dpi = 300, 
       width = 7, height = 5)



