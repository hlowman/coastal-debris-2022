# Figure 1 Map Creation
# Heili Lowman
# 11/15/20
# Edited: 3/19/21

# The following script will create the map to be used in the lignin manuscript.
# It will use the sf package rather than creating the map by hand in powerpoint.

# Load packages
library(tidyverse)
library(ggplot2)
library(ggmap)
library(ggrepel)
library(sf)
library(USAboundaries)
library(USAboundariesData)
library(ggspatial)

# Load data
map_data <- read_csv("rapid_sites.csv")

map_df <- map_data[1:19,1:5] %>% # Removing odd rows/columns of NAs.
  mutate(lon = Lon) %>%
  mutate(lat = Lat) %>%
  mutate(env_f = factor(Environment, levels = c("Stream", "Estuary", "Marine near Stream", "Marine far from Stream"))) %>%
  mutate(site_f = factor(LTER_code, levels = c("GV01", "HO00", "AHND", "AQUE", "RG01", "REFU", "BC02", "IVEE", "TE03", "SP02", "GOSL", "GOLB", "ATMY", "ABURE", "AB00", "ABUR", "MOHK", "MC00", "MICR"))) %>%
  mutate(n_x = c(0.00,-0.01,0.015,-0.01,0.0,-0.01,0.015,0.0,0.0,-0.03,-0.035,-0.02,0.01,0.00,0.0,0.01,0.01,-0.01,0.0)) %>%
  mutate(n_y = c(0.02,0.02,0.02,0.02,0.02,0.02,0.015,0.02,0.02,0.0,0.0,-0.02,-0.02,-0.02,-0.02,-0.02,-0.015,-0.02,0.0)) %>%
  mutate(fill_color = ifelse(env_f == "Marine near Stream", 2, 1))

# Create data sf object
map_sf <- st_as_sf(map_df,
                   coords = c("lon", "lat"),
                   remove = F,
                   crs = 4326) # WGS 84 projection

# Base plot to see how things are looking...
plot(map_sf$geometry)

# ggplot + sf map

# create CA county dataset for use in map below using USAboundaries
CA_counties <- us_counties(states = "California")

SB_county <- CA_counties %>%
  filter(name == "Santa Barbara")

# create base terrain map tile

# create bounding box
lats <- c(34.368063, 34.569875)
lons <- c(-120.287304, -119.634991)
bb <- make_bbox(lon = lons, lat = lats, f = 0.05)

sb_basemap <- get_stamenmap(bb, 
                      zoom = 12,
                      maptype = 'terrain-background')

ggmap(sb_basemap)
  
fullmap <- ggmap(sb_basemap) + # base google maps tile
  geom_point(data = map_sf, aes(x = lon, y = lat, shape = env_f, fill = env_f), 
          size = 4,
          shape = 21,
          inherit.aes = FALSE) + # adds points
  scale_shape_manual(values = c(21, 24, 22, 22)) +
  scale_fill_manual(values = c("white", "lightgrey", "grey40", "black")) +
  geom_text_repel(data = map_sf, 
                     aes(x = lon, 
                         y = lat, 
                         label = site_f),
                  nudge_x = map_sf$n_x,
                  nudge_y = map_sf$n_y,
                  segment.size = 0.2,
                  size = 5) +
  ggspatial::annotation_north_arrow(location = "tr") + # adds compass due north
  ggspatial::annotation_scale() + # adds scale
  geom_text(x = -120, y = 34.40504, label = "Santa Barbara Channel", color = "gray40", size = 4, fontface = "italic") +
  geom_text(x = -119.95, y = 34.5, label = "Santa Ynez Mountains", color = "gray10", size = 4, fontface = "italic") +
  labs(x = "Longitude (WGS84)",
       y = "Latitude",
       shape = "Environment",
       fill = "Environment") +
  theme_bw() +
  theme(legend.position = c(0.1,0.25),
        legend.background = element_rect(fill = "white", size = 0.5, linetype = "solid")) +
  coord_sf(crs = st_crs(4326))

fullmap

# Export map to desktop.

# ggsave(("Figure_1.tiff"),
#        path = "/Users/heililowman/Desktop/R_figures/Lignin",
#        width = 30,
#        height = 15,
#        units = "cm"
#        )

# End of R script.
