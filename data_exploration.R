# Coastal Fund Data Exploration
# October 29, 2021
# Heili Lowman

# The following script will do a bit of data exploration with the newly
# curated dataset combining both lignin and pyC measures.

# Load packages
library(tidyverse)
library(here)
library(lubridate)
library(calecopal)
library(patchwork)

# Load data
dat <- read_csv("data_raw/lignin_pyc_results.csv") 

# Tidy data
dat_ed <- dat %>%
  mutate(Date = mdy(Collection_Date)) %>% # properly format date
  mutate(Location_f = factor(Location_id, levels = c("MD", "DD", "CD", "GBEA", 
                                                     "GOSL", "GOBW", "GOBE"))) %>% # level locations
  mutate(Environment_f = factor(case_when(Location %in% c("Montecito Debris",
                                                 "Goleta Beach",
                                                 "Goleta Slough",
                                                 "Debris Deposition Site",
                                                 "Clay Debris") ~ "Beach",
                                 Location %in% c("Goleta Bay West",
                                                 "Goleta Bay East") ~ "Ocean"))) # create environment column
# Create summary tables
# First, parsing out core section as a level
summary1 <- dat_ed %>%
  group_by(Environment_f, Location_f, Date, Water_Depth, Core_Section) %>%
  summarize(D13C = mean(d13C, na.rm = TRUE), sdD13C = sd(d13C, na.rm = TRUE),
                      OC = mean(perc_OC, na.rm = TRUE), sdOC = sd(perc_OC, na.rm = TRUE),
                      PYC = mean(pyC, na.rm = TRUE), sdPYC = sd(pyC, na.rm = TRUE),
                      LAM = mean(Lambda, na.rm = TRUE), sdLAM = sd(Lambda, na.rm = TRUE),
                      SIG = mean(Sigma8, na.rm = TRUE), sdSIG = sd(Sigma8, na.rm = TRUE),
                      S = mean(S, na.rm = TRUE), sds = sd(S, na.rm = TRUE),
                      V = mean(V, na.rm = TRUE), sdv = sd(V, na.rm = TRUE),
                      C = mean(C, na.rm = TRUE), sdc = sd(C, na.rm = TRUE),
                      P = mean(P, na.rm = TRUE), sdp = sd(P, na.rm = TRUE),
                      SV = mean(SV, na.rm = TRUE), sdsv = sd(SV, na.rm = TRUE),
                      CV = mean(CV, na.rm = TRUE), sdcv = sd(CV, na.rm = TRUE),
                      PVS = mean(PVS, na.rm = TRUE), sdpvs = sd(PVS, na.rm = TRUE),
                      BDV = mean(BdV, na.rm = TRUE), sdbdv = sd(BdV, na.rm = TRUE),
                      AAS = mean(AcAds, na.rm = TRUE), sdbdv = sd(AcAds, na.rm = TRUE),
                      AAV = mean(AcAdv, na.rm = TRUE), sdbdv = sd(AcAdv, na.rm = TRUE)) %>%
  ungroup()

# Then, looking at cores as a whole
summary2 <- dat_ed %>%
  group_by(Environment_f, Location_f, Date, Water_Depth) %>%
  summarize(D13C = mean(d13C, na.rm = TRUE), sdD13C = sd(d13C, na.rm = TRUE),
            OC = mean(perc_OC, na.rm = TRUE), sdOC = sd(perc_OC, na.rm = TRUE),
            PYC = mean(pyC, na.rm = TRUE), sdPYC = sd(pyC, na.rm = TRUE),
            LAM = mean(Lambda, na.rm = TRUE), sdLAM = sd(Lambda, na.rm = TRUE),
            SIG = mean(Sigma8, na.rm = TRUE), sdSIG = sd(Sigma8, na.rm = TRUE),
            S = mean(S, na.rm = TRUE), sds = sd(S, na.rm = TRUE),
            V = mean(V, na.rm = TRUE), sdv = sd(V, na.rm = TRUE),
            C = mean(C, na.rm = TRUE), sdc = sd(C, na.rm = TRUE),
            P = mean(P, na.rm = TRUE), sdp = sd(P, na.rm = TRUE),
            SV = mean(SV, na.rm = TRUE), sdsv = sd(SV, na.rm = TRUE),
            CV = mean(CV, na.rm = TRUE), sdcv = sd(CV, na.rm = TRUE),
            PVS = mean(PVS, na.rm = TRUE), sdpvs = sd(PVS, na.rm = TRUE),
            BDV = mean(BdV, na.rm = TRUE), sdbdv = sd(BdV, na.rm = TRUE),
            AAS = mean(AcAds, na.rm = TRUE), sdbdv = sd(AcAds, na.rm = TRUE),
            AAV = mean(AcAdv, na.rm = TRUE), sdbdv = sd(AcAdv, na.rm = TRUE)) %>%
  ungroup()

# Export this for inclusion in the manuscript as table 1
# write_csv(summary2, path = "data_working/table1.csv")

# Generate plots of select data
dat_plots <- summary2 %>%
  # replacing the Montecito debris category for ease of plotting
  mutate(loc_plot = factor(case_when(Location_f == "MD" | Location_f == "DD" ~ "DDS",
                                     TRUE ~ as.character(Location_f)),
                           levels = c("DDS", "CD", "GBEA", "GOSL", "GOBW", "GOBE"))) %>%
  # and log-transforming pyC
  mutate(log_pyc = log10(PYC))

# Pyrogenic C results through time at all sites
(fig1 <- dat_plots %>%
  ggplot(aes(x = Date, y = log_pyc, color = loc_plot)) +
  geom_point() +
  scale_color_manual(name = "Sampling Location", 
                     values = cal_palette("fire", n = 6, type = "continuous")) + # custom colors
  facet_grid(rows = vars(loc_plot)) +
  theme_bw() +
  theme(legend.position = "none"))

# Pyrogenic C results vs. d13C values at all sites
(fig2 <- dat_plots %>%
  ggplot(aes(x = D13C, y = log_pyc, color = loc_plot)) +
  geom_point() +
  scale_color_manual(name = "Sampling Location", 
                     values = cal_palette("fire", n = 6, type = "continuous")) +
  theme_bw() +
  theme(legend.position = "none"))

# Lambda results through time at all sites
(fig3 <- dat_plots %>%
    ggplot(aes(x = Date, y = LAM, color = loc_plot)) +
    geom_point() +
    scale_color_manual(name = "Sampling Location", 
                       values = cal_palette("fire", n = 6, type = "continuous")) + # custom colors
    facet_grid(rows = vars(loc_plot)) +
    theme_bw() +
    theme(legend.position = "none"))

# Lambda results vs. d13C values at all sites
(fig4 <- dat_plots %>%
    ggplot(aes(x = D13C, y = LAM, color = loc_plot)) +
    geom_point() +
    scale_color_manual(name = "Sampling Location", 
                       values = cal_palette("fire", n = 6, type = "continuous")) +
    theme_bw())

# Combine above plots into a single figure and export for easier viewing.
(fig1234 <- (fig1 + fig3) /
  (fig2 + fig4))

# ggsave(("figures/exploration/pyc_lam_fig.png"),
#        width = 16,
#        height = 20,
#        units = "cm"
# )

# SV vs. CV values at all sites
(fig5 <- dat_plots %>%
    ggplot(aes(x = log10(CV), y = log10(SV), color = loc_plot)) +
    geom_point() +
    scale_color_manual(name = "Sampling Location", 
                       values = cal_palette("fire", n = 6, type = "continuous")) +
    theme_bw())

# PVS vs. 35BDV values at all sites
(fig6 <- dat_plots %>%
    ggplot(aes(x = log10(PVS), y = log10(BDV), color = loc_plot)) +
    geom_point() +
    scale_color_manual(name = "Sampling Location", 
                       values = cal_palette("fire", n = 6, type = "continuous")) +
    theme_bw())

# Need to examine same measures in core depths at all marine sites sampled.

#### Results calculations ####

# % OC measures

min(dat_ed$perc_OC) # 0.05
max(dat_ed$perc_OC) # 3.59

a <- dat_ed %>%
  group_by(Environment_f) %>%
  summarize(meanOC = mean(perc_OC, na.rm = TRUE),
            sdOC = sd(perc_OC, na.rm = TRUE))

a1 <- dat_ed %>%
  filter(Location_id %in% c("GBEA", "GOSL")) %>%
  group_by(Date) %>%
  summarize(meanOC = mean(perc_OC, na.rm = TRUE),
            sdOC = sd(perc_OC, na.rm = TRUE))

a2 <- dat_ed %>%
  filter(Environment_f == "Ocean") %>%
  group_by(Location_id) %>%
  summarize(meanOC = mean(perc_OC, na.rm = TRUE),
            sdOC = sd(perc_OC, na.rm = TRUE))

a3 <- dat_ed %>%
  filter(Environment_f == "Ocean") %>%
  group_by(Water_Depth) %>%
  summarize(meanOC = mean(perc_OC, na.rm = TRUE),
            sdOC = sd(perc_OC, na.rm = TRUE))

a4 <- dat_ed %>%
  filter(Environment_f == "Ocean") %>%
  group_by(Core_Section) %>%
  summarize(meanOC = mean(perc_OC, na.rm = TRUE),
            sdOC = sd(perc_OC, na.rm = TRUE))

# % OC figure
(fig7 <- ggplot(dat_plots) +
  geom_point(aes(x = Date, y = OC, color = loc_plot)) +
  scale_color_manual(name = "Sampling Location", 
                     values = cal_palette("fire", n = 6, type = "continuous")) +
  theme_bw())

# δ13C measures

min(dat_ed$d13C) # -27.05
max(dat_ed$d13C) # -21.52

b <- dat_ed %>%
  group_by(Environment_f) %>%
  summarize(mean13C = mean(d13C, na.rm = TRUE),
            sd13C = sd(d13C, na.rm = TRUE))

b1 <- dat_ed %>%
  filter(Location_id %in% c("GBEA", "GOSL")) %>%
  group_by(Date) %>%
  summarize(meand13C = mean(d13C, na.rm = TRUE),
            sdd13C = sd(d13C, na.rm = TRUE))

b2 <- dat_ed %>%
  filter(Environment_f == "Ocean") %>%
  group_by(Location_id) %>%
  summarize(mean13C = mean(d13C, na.rm = TRUE),
            sd13C = sd(d13C, na.rm = TRUE))

b3 <- dat_ed %>%
  filter(Environment_f == "Ocean") %>%
  group_by(Water_Depth) %>%
  summarize(mean13C = mean(d13C, na.rm = TRUE),
            sd13C = sd(d13C, na.rm = TRUE))

b4 <- dat_ed %>%
  filter(Environment_f == "Ocean") %>%
  group_by(Core_Section) %>%
  summarize(mean13C = mean(d13C, na.rm = TRUE),
            sd13C = sd(d13C, na.rm = TRUE))

# d13C figure
(fig8 <- ggplot(dat_plots) +
    geom_point(aes(x = Date, y = D13C, color = loc_plot)) +
    scale_color_manual(name = "Sampling Location", 
                       values = cal_palette("fire", n = 6, type = "continuous")) +
    theme_bw())

# End of script.
