# Figure 1 Map Creation
# Heili Lowman
# 11/15/20
# Edited: 3/19/21

# The following script will create the map to be used in the manuscript.
# It will use the sf package rather than creating the map by hand in powerpoint.

# Load packages
library(tidyverse)
library(ggplot2)
library(ggmap)
library(ggrepel)
library(sf)
library(USAboundaries)
#library(USAboundariesData)
library(ggspatial)
library(calecopal)
library(patchwork)
library(jpeg)
library(ggpubr)

# Load data
map_data <- read_csv("data_raw/site_locations.csv") %>%
  mutate(print_names = c("Disposal Site", "Goleta Beach", "Goleta Slough", "10m", "20m", "5m", "10m", "20m")) # add in names to actually be displayed

map_df <- map_data %>%
  mutate(lon = Long) %>%
  mutate(lat = Lat) %>%
  mutate(site_f = factor(Site, levels = c("Deposition Site", "Goleta Beach", "Goleta Slough", "Goleta Bay West 5m", "Goleta Bay West 10m", "Goleta Bay West 20m", "Goleta Bay East 10m", "Goleta Bay East 20m"))) %>%
  # adding nudging distances for each site's label
  mutate(n_x = c(-0.006,0.002,0.004,0.006,0.005,-0.004,-0.004,-0.004)) %>%
  mutate(n_y = c(0.002,0.004,0.002,0.00,0.00,0.00,0.00,0.00))

# Create data sf object
map_sf <- st_as_sf(map_df,
                   coords = c("lon", "lat"),
                   remove = F,
                   crs = 4326) # WGS 84 projection

# Base plot to see how things are looking...
plot(map_sf$geometry)

#### INSET ####

# create base terrain map tile
# create bounding box
lats <- c(34.395, 34.425)
lons <- c(-119.855, -119.80)
bb <- make_bbox(lon = lons, lat = lats, f = 0.05)

sb_basemap <- get_stamenmap(bb, 
                      zoom = 13,
                      maptype = 'terrain-background')

ggmap(sb_basemap)
  
(insetmap <- ggmap(sb_basemap) + # base google maps tile
  geom_point(data = map_sf, aes(x = lon, y = lat), 
          size = 3,
          shape = 16,
          inherit.aes = FALSE) + # adds points
  geom_text_repel(data = map_sf, 
                     aes(x = lon, 
                         y = lat, 
                         label = print_names),
                  nudge_x = map_sf$n_x,
                  nudge_y = map_sf$n_y,
                  segment.size = 0.2,
                  size = 4) +
  geom_text(x = -119.835, y = 34.4095, label = "W Goleta Bay", color = "gray40", size = 4, fontface = "italic") + # labels the west side
    geom_text(x = -119.8238, y = 34.414, label = "E Goleta Bay", color = "gray40", size = 4, fontface = "italic") + # labels the east side
  ggspatial::annotation_scale() + # adds scale
  labs(x = "Longitude (WGS84)",
       y = "Latitude") +
  theme_bw() +
  theme(legend.position = "none") +
  coord_sf(crs = st_crs(4326)))

#### SBC FULL MAP ####

# Create a new dataframe for plotting the cities on a larger map
cities <- c("Santa Barbara", "Montecito", "Goleta")
lat <- c(34.419118, 34.420693, 34.430335)
lon <- c(-119.698701, -119.643306, -119.869271)
city_map_data <- data.frame(cities, lat, lon)

city_map_data <- city_map_data %>%
  # adding nudging distances for each site's label
  mutate(n_x = c(0.00, 0.00, 0.00)) %>%
  mutate(n_y = c(0.01, 0.01, 0.01))

# Create data sf object
city_sf <- st_as_sf(city_map_data,
                   coords = c("lon", "lat"),
                   remove = F,
                   crs = 4326) # WGS 84 projection

# Base plot to see how things are looking...
#plot(city_sf$geometry)

# create bounding box
city_lats <- c(34.35, 34.50)
city_lons <- c(-119.9, -119.4)
city_bb <- make_bbox(lon = city_lons, lat = city_lats, f = 0.05)

city_basemap <- get_stamenmap(city_bb, 
                            zoom = 12,
                            maptype = 'terrain-background')

ggmap(city_basemap)

# Downloaded Thomas Fire perimeter from NIFC database.
# https://data-nifc.opendata.arcgis.com/datasets/historic-perimeters-combined-2000-2018/explore?filters=eyJpbmNpZGVudG5hbWUiOlsiVGhvbWFzIl19&location=34.462055%2C-119.400630%2C10.00&showTable=true

tfp <- read_sf("data_raw/Historic_Perimeters_Combined_2000-2018_Thomas_Fire/US_HIST_FIRE_PERIMTRS_2000_2018_DD83.shp")

# select for only the perimeter determined on Jan 6, 2018
final_perimeter <- tfp[1,]
# transform to EPCG 3857
fp_3857 <- st_transform(final_perimeter, 3857)
# more info on why this is necessary here: https://stackoverflow.com/questions/47749078/how-to-put-a-geom-sf-produced-map-on-top-of-a-ggmap-produced-raster

(fullmap <- ggmap(city_basemap) + # base google maps tile
    coord_sf(crs = st_crs(3857)) + # force the ggplot2 map to be in 3857
    geom_sf(data = fp_3857, color = "firebrick1", fill = "firebrick3",
            alpha = 0.65, inherit.aes = FALSE) +
  geom_point(data = city_sf, aes(x = lon, y = lat), 
             size = 3, shape = 15, inherit.aes = FALSE) + # adds points
  geom_text_repel(data = city_sf, 
                  aes(x = lon, y = lat, label = cities),
                  nudge_x = city_sf$n_x,
                  nudge_y = city_sf$n_y,
                  segment.color = 'transparent',
                  size = 4) +
  ggspatial::annotation_north_arrow(location = "tr") + # adds compass due north
  ggspatial::annotation_scale() + # adds scale
  geom_text(x = -119.6, y = 34.37, label = "Santa Barbara Channel", color = "gray40", size = 6, fontface = "italic") + # labels the channel
  geom_text(x = -119.8, y = 34.485, label = "Santa Ynez Mountains", color = "gray10", size = 6, fontface = "italic") + # labels the mountains
  geom_rect(xmin = -119.8578, xmax = -119.7972, ymin = 34.3935, ymax = 34.4265,
            color = "gray10", fill = NA) +
  labs(x = "Longitude (WGS84)",
       y = "Latitude") + # labels axes
  theme_bw() +
  theme(legend.position = "none") + # removes the legend
  coord_sf(crs = st_crs(4326)))

#### Disposal image ####
bulldozer <- readJPEG('figures/thumbnail_IMG_1087.jpeg')

bulldozer_plot <- ggplot() +
  background_image(bulldozer)
  # coord_fixed()

after_beach <- readJPEG('figures/thumbnail_IMG_6677.jpeg')

beach_plot <- ggplot() +
  background_image(after_beach)

g_beach <- readJPEG('figures/thumbnail_IMG_6680.jpeg')

gbeach_plot <- ggplot() +
  background_image(g_beach)

g_slough <- readJPEG('figures/thumbnail_IMG_6686.jpeg')

gslough_plot <- ggplot() +
  background_image(g_slough)

# Not playing nice when in it's own in a line, so
# have to put it on a line with the inset map.

# Combine the maps together.
(fig1 <- (fullmap + insetmap) +
    plot_annotation(tag_levels = 'A'))

# Export map.

# ggsave(fig1,
#        filename = "figures/Fig1_map_0622.png",
#        width = 40,
#        height = 20,
#        units = "cm"
#        )

# End of R script.
