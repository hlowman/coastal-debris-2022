# Coastal Fund Data Analysis & Figure Creation
# October 29, 2021
# Heili Lowman

# The following script will do a bit of data exploration and analysis with the newly
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
                      s = mean(S, na.rm = TRUE), sds = sd(S, na.rm = TRUE),
                      v = mean(V, na.rm = TRUE), sdv = sd(V, na.rm = TRUE),
                      c = mean(C, na.rm = TRUE), sdc = sd(C, na.rm = TRUE),
                      p = mean(P, na.rm = TRUE), sdp = sd(P, na.rm = TRUE),
                      sv = mean(SV, na.rm = TRUE), sdsv = sd(SV, na.rm = TRUE),
                      cv = mean(CV, na.rm = TRUE), sdcv = sd(CV, na.rm = TRUE),
                      pvs = mean(PVS, na.rm = TRUE), sdpvs = sd(PVS, na.rm = TRUE),
                      bdv = mean(BdV, na.rm = TRUE), sdbdv = sd(BdV, na.rm = TRUE),
                      aas = mean(AcAds, na.rm = TRUE), sdbdv = sd(AcAds, na.rm = TRUE),
                      aav = mean(AcAdv, na.rm = TRUE), sdbdv = sd(AcAdv, na.rm = TRUE)) %>%
  ungroup()

# Then, looking at cores as a whole
summary2 <- dat_ed %>%
  group_by(Environment_f, Location_f, Date, Water_Depth) %>%
  summarize(D13C = mean(d13C, na.rm = TRUE), sdD13C = sd(d13C, na.rm = TRUE),
            OC = mean(perc_OC, na.rm = TRUE), sdOC = sd(perc_OC, na.rm = TRUE),
            PYC = mean(pyC, na.rm = TRUE), sdPYC = sd(pyC, na.rm = TRUE),
            LAM = mean(Lambda, na.rm = TRUE), sdLAM = sd(Lambda, na.rm = TRUE),
            SIG = mean(Sigma8, na.rm = TRUE), sdSIG = sd(Sigma8, na.rm = TRUE),
            s = mean(S, na.rm = TRUE), sds = sd(S, na.rm = TRUE),
            v = mean(V, na.rm = TRUE), sdv = sd(V, na.rm = TRUE),
            c = mean(C, na.rm = TRUE), sdc = sd(C, na.rm = TRUE),
            p = mean(P, na.rm = TRUE), sdp = sd(P, na.rm = TRUE),
            sv = mean(SV, na.rm = TRUE), sdsv = sd(SV, na.rm = TRUE),
            cv = mean(CV, na.rm = TRUE), sdcv = sd(CV, na.rm = TRUE),
            pvs = mean(PVS, na.rm = TRUE), sdpvs = sd(PVS, na.rm = TRUE),
            bdv = mean(BdV, na.rm = TRUE), sdbdv = sd(BdV, na.rm = TRUE),
            aas = mean(AcAds, na.rm = TRUE), sdbdv = sd(AcAds, na.rm = TRUE),
            aav = mean(AcAdv, na.rm = TRUE), sdbdv = sd(AcAdv, na.rm = TRUE)) %>%
  ungroup()

# Export this for inclusion in the manuscript as table 1
#write_csv(summary2, path = "data_working/table1.csv")

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

#### Manuscript Figure 2 ####
dat_F1A <- dat_ed %>%
  group_by(Location_f, Date) %>%
  summarize(D13C = mean(d13C, na.rm = TRUE), 
            sdD13C = sd(d13C, na.rm = TRUE),
            PYC = mean(pyC, na.rm = TRUE), 
            sdPYC = sd(pyC, na.rm = TRUE),
            LAM = mean(Lambda, na.rm = TRUE), 
            sdLAM = sd(Lambda, na.rm = TRUE)) %>%
  ungroup() %>%
  # replacing the Montecito debris category for ease of plotting
  mutate(loc_plot = factor(case_when(Location_f == "MD" | Location_f == "DD" ~ "DDS",
                                     TRUE ~ as.character(Location_f)),
                           levels = c("DDS", "GBEA", "GOSL", "GOBW", "GOBE", "CD"))) %>%
  mutate(date_plot = factor(case_when(Date == "2018-02-02" ~ "A",
                               Date == "2018-02-23" ~ "B",
                               Location_f == "CD" ~ NA_character_,
                               TRUE ~ "C")))

# Pyrogenic C results vs. d13C values at all sites
(F1A <- dat_F1A %>%
    ggplot(aes(x = D13C, y = PYC, fill = date_plot, shape = loc_plot)) +
    geom_point(size = 3) +
    scale_fill_manual(name = "Sampling Date", 
                      values = cal_palette("sierra2"),
                      labels = c("Feb 2", "Feb 23", "Apr 23-24"),
                      guide = 'none' # remove legend
                      ) +
    scale_shape_manual(name = "Sampling Site",
                       values = c(21, 22, 23, 24, 25, 19),
                       labels = c("Deposition Site", "Goleta Beach", "Goleta Slough", "Goleta Bay West", "Goleta Bay East", "Clay Debris"),
                       guide = 'none' # remove legend
                       ) +
    labs(y = "% Pyrogenic Carbon",
         x = expression("δ"^{13}*"C (‰)")) +
    theme_bw())

# For clarity, I'm going to plot these through time, rather than against d13C
date_labels <- c("February 2", "February 23", "April 24")
site_labels <- c(`DDS` = "Disposal Site",
                 `GBEA` = "Goleta Beach",
                 `GOSL` = "Goleta Slough")
(fig_beachpyc <- dat_plots %>%
    filter(Environment_f == "Beach") %>% # only beach samples
    filter(Location_f != "CD") %>% # remove the clay debris layer
    mutate(Location_b = factor(loc_plot)) %>%
    mutate(Date_f = factor(case_when(Date == "2018-02-02" ~ "A",
                                     Date == "2018-02-23" ~ "B",
                                     Date == "2018-04-24" ~ "C"))) %>% # add new column for date as factor
    ggplot(aes(x = Date_f, y = PYC)) +
    geom_point() +
    geom_pointrange(aes(ymin=PYC-sdPYC, ymax=PYC+sdPYC)) +
    labs(x = "Date", y = "% Pyrogenic Carbon") +
    scale_x_discrete(labels = date_labels) +
    facet_wrap(.~Location_b, scales = "free",
               nrow = 3,
               labeller = as_labeller(site_labels)) +
    theme_bw() +
    theme(strip.background = element_blank(),
    legend.position = "none"))

# Lambda results vs. d13C values at all sites
(F1B <- dat_F1A %>%
    ggplot(aes(x = D13C, y = LAM, fill = date_plot, shape = loc_plot)) +
    geom_point(size = 3) +
    scale_fill_manual(name = "Sampling Date", 
                      values = cal_palette("sierra2"),
                      labels = c("Feb 2", "Feb 23", "Apr 23-24"),
                      na.translate = F) +
    scale_shape_manual(name = "Sampling Site",
                       values = c(21, 22, 23, 24, 25, 19),
                       labels = c("Deposition Site", "Goleta Beach", "Goleta Slough", "Goleta Bay West", "Goleta Bay East", "Low Tide Line")) +
    labs(y = "Λ (mg/100 mg OC)",
         x = expression("δ"^{13}*"C (‰)")) +
    # have to override the aes for the legend
    guides(fill = guide_legend(override.aes=list(shape=21))) +
    theme_bw())

# create similar paneled lambda figure as above
(fig_beachlam <- dat_plots %>%
    filter(Environment_f == "Beach") %>% # only beach samples
    filter(Location_f != "CD") %>% # remove the clay debris layer
    mutate(Location_b = factor(loc_plot)) %>%
    mutate(Date_f = factor(case_when(Date == "2018-02-02" ~ "A",
                                     Date == "2018-02-23" ~ "B",
                                     Date == "2018-04-24" ~ "C"))) %>% # add new column for date as factor
    ggplot(aes(x = Date_f, y = LAM)) +
    geom_point() +
    geom_pointrange(aes(ymin=LAM-sdLAM, ymax=LAM+sdLAM)) +
    labs(x = "Date", y = "Λ (mg/100 mg OC)") +
    scale_x_discrete(labels = date_labels) +
    facet_wrap(.~Location_b, scales = "free",
               nrow = 3,
               labeller = as_labeller(site_labels)) +
    theme_bw() +
    theme(strip.background = element_blank(),
          legend.position = "none"))

# Combine above plots into a single figure and export for easier viewing.
(fig1_manuscript <- (F1A + F1B) +
    plot_annotation(tag_levels = 'A'))

# ggsave(("figures/Fig2_pyc_lam.png"),
#        width = 20,
#        height = 8,
#        units = "cm"
# )

(fig_panels_manuscript <- (fig_beachpyc + fig_beachlam) +
    plot_annotation(tag_levels = 'A'))

# ggsave(("figures/Fig2_paneled_freescales_pyc_lam.png"),
#        width = 15,
#        height = 15,
#        units = "cm"
# )

# SV vs. CV values at all sites
(fig5 <- dat_plots %>%
    ggplot(aes(x = log10(cv), y = log10(sv), color = loc_plot)) +
    geom_point() +
    scale_color_manual(name = "Sampling Location", 
                       values = cal_palette("fire", n = 6, type = "continuous")) +
    theme_bw())

# PVS vs. 35BDV values at all sites
(fig6 <- dat_plots %>%
    ggplot(aes(x = log10(pvs), y = log10(bdv), color = loc_plot)) +
    geom_point() +
    scale_color_manual(name = "Sampling Location", 
                       values = cal_palette("fire", n = 6, type = "continuous")) +
    theme_bw())

# Since many of the degradation measures were low in beach sediments,
# going to focus those figures on marine sediment.

# Pyrogenic C results vs. d13C values at marine sites
(F2A <- dat_plots %>%
    filter(Environment_f == "Ocean") %>% # removing beach data
    mutate(Depth_f = factor(Water_Depth)) %>%
    ggplot(aes(x = D13C, y = PYC, fill = Depth_f, shape = Location_f)) +
    geom_point(size = 3) +
    scale_fill_manual(name = "Sampling Depth", 
                      values = c("#A1CAF6", "#4C6FA1", "#1E2F46"),
                      labels = c("5 m", "10 m", "20 m"),
                      guide = 'none' # remove legend
    ) +
    scale_shape_manual(name = "Sampling Site",
                       values = c(24, 25),
                       labels = c("Goleta Bay West", "Goleta Bay East"),
                       guide = 'none' # remove legend
    ) +
    labs(y = "% Pyrogenic Carbon",
         x = expression("δ"^{13}*"C (‰)")) +
    theme_bw())

# making clearer plots after Andy's suggestions
(fig_marpyc <- dat_plots %>%
    filter(Environment_f == "Ocean") %>% # only beach samples
    mutate(Location_m = factor(loc_plot)) %>%
    mutate(Depth_f = factor(Water_Depth)) %>%
    ggplot(aes(x = Depth_f, y = PYC, color = Location_m)) +
    geom_point(position=position_dodge(width=0.5)) +
    geom_pointrange(aes(ymin=PYC-sdPYC, ymax=PYC+sdPYC),
                    position=position_dodge(width=0.5)) +
    scale_color_manual(values = c("#0FB2D3", "#026779")) +
    labs(x = "Water Depth (m)", y = "% Pyrogenic Carbon", color = "Site") +
    theme_bw() +
    theme(legend.position = "none"))

# Lambda results vs. d13C values at marine sites
(F2B <- dat_plots %>%
    filter(Environment_f == "Ocean") %>% # removing beach data
    mutate(Depth_f = factor(Water_Depth)) %>%
    ggplot(aes(x = D13C, y = LAM, fill = Depth_f, shape = Location_f)) +
    geom_point(size = 3) +
    scale_fill_manual(name = "Sampling Depth", 
                      values = c("#A1CAF6", "#4C6FA1", "#1E2F46"),
                      labels = c("5 m", "10 m", "20 m")) +
    scale_shape_manual(name = "Sampling Site",
                       values = c(24, 25),
                       labels = c("Goleta Bay West", "Goleta Bay East")) +
    labs(y = "Λ (mg/100 mg OC)",
         x = expression("δ"^{13}*"C (‰)")) +
    # have to override the aes for the legend
    guides(fill = guide_legend(override.aes=list(shape=21))) +
    theme_bw())

# making clearer plots after Andy's suggestions
(fig_marlam <- dat_plots %>%
    filter(Environment_f == "Ocean") %>% # only beach samples
    mutate(Location_m = factor(loc_plot)) %>%
    mutate(Depth_f = factor(Water_Depth)) %>%
    ggplot(aes(x = Depth_f, y = LAM, color = Location_m)) +
    geom_point(position=position_dodge(width=0.5)) +
    geom_pointrange(aes(ymin=LAM-sdLAM, ymax=LAM+sdLAM),
                    position=position_dodge(width=0.5)) +
    scale_color_manual(values = c("#0FB2D3", "#026779"),
                       labels = c("West Goleta Bay", "East Goleta Bay")) +
    labs(x = "Water Depth (m)", y = "Λ (mg/100 mg OC)", color = "Site") +
    theme_bw())

(fig_panels_marine_manuscript <- (fig_marpyc + fig_marlam) +
    plot_annotation(tag_levels = 'A'))

# ggsave(("figures/Fig3_paneled_marine_pyc_lam.png"),
#        width = 18,
#        height = 8,
#        units = "cm"
# )

# S/V vs. C/V values at marine sites
(F2C <- dat_plots %>%
    filter(Environment_f == "Ocean") %>% # removing beach data
    mutate(Depth_f = factor(Water_Depth)) %>%
    ggplot(aes(x = cv, y = sv, fill = Depth_f, shape = Location_f)) +
    geom_point(size = 3) +
    scale_fill_manual(name = "Sampling Depth", 
                      values = c("#A1CAF6", "#4C6FA1", "#1E2F46"),
                      labels = c("5 m", "10 m", "20 m"),
                      guide = 'none' # remove legend
    ) +
    scale_shape_manual(name = "Sampling Site",
                       values = c(24, 25),
                       labels = c("Goleta Bay West", "Goleta Bay East"),
                       guide = 'none' # remove legend
    ) +
    labs(y = "Syringyl / Vanillyl",
         x = "Cinnamyl / Vanillyl") +
    theme_bw())

# making clearer plots after Andy's suggestions
(fig_mar_svcv <- dat_plots %>%
    filter(Environment_f == "Ocean") %>% # only beach samples
    mutate(Location_m = factor(loc_plot)) %>%
    mutate(Depth_f = factor(Water_Depth)) %>%
    ggplot(aes(x = cv, y = sv, color = Location_m)) +
    geom_point(aes(shape = Depth_f), size = 3) +
    geom_linerange(aes(ymin=sv-sdsv, ymax=sv+sdsv)) +
    geom_linerange(aes(xmin=cv-sdcv, xmax=cv+sdcv)) +
    scale_color_manual(values = c("#0FB2D3", "#026779"),
                       labels = c("West Goleta Bay", "East Goleta Bay")) +
    labs(x = "Cinnamyl / Vanillyl", 
         y = "Syringyl / Vanillyl", 
         color = "Site",
         shape = "Water Depth (m)") +
    theme_bw() +
    theme(legend.position = "none"))

# P/V+S vs. 3,5Bd/V values at marine sites
(F2D <- dat_plots %>%
    filter(Environment_f == "Ocean") %>% # removing beach data
    mutate(Depth_f = factor(Water_Depth)) %>%
    ggplot(aes(x = bdv, y = pvs, fill = Depth_f, shape = Location_f)) +
    geom_point(size = 3) +
    scale_fill_manual(name = "Sampling Depth", 
                      values = c("#A1CAF6", "#4C6FA1", "#1E2F46"),
                      labels = c("5 m", "10 m", "20 m"),
                      guide = 'none' # remove legend
    ) +
    scale_shape_manual(name = "Sampling Site",
                       values = c(24, 25),
                       labels = c("Goleta Bay West", "Goleta Bay East"),
                       guide = 'none' # remove legend
    ) +
    labs(y = "P-hydroxyl / (Vanillyl + Syringyl)",
         x = "3,5-Bd / Vanillyl") +
    theme_bw())

# making clearer plots after Andy's suggestions
(fig_mar_pvsbdv <- dat_plots %>%
    filter(Environment_f == "Ocean") %>% # only beach samples
    mutate(Location_m = factor(loc_plot)) %>%
    mutate(Depth_f = factor(Water_Depth)) %>%
    ggplot(aes(x = bdv, y = pvs, color = Location_m)) +
    geom_point(aes(shape = Depth_f), size = 3) +
    geom_linerange(aes(ymin=pvs-sdpvs, ymax=pvs+sdpvs)) +
    geom_linerange(aes(xmin=bdv-sdbdv, xmax=bdv+sdbdv)) +
    scale_color_manual(values = c("#0FB2D3", "#026779"),
                       labels = c("West Goleta Bay", "East Goleta Bay")) +
    labs(x = "3,5-Bd / Vanillyl", 
         y = "P-hydroxyl / (Vanillyl + Syringyl)", 
         color = "Site",
         shape = "Water Depth (m)") +
    theme_bw())

# Combine above plots into a single figure and export for easier viewing.
(fig2_manuscript <- (F2A + F2B) / (F2C + F2D) +
    plot_annotation(tag_levels = 'A'))

# ggsave(("figures/Fig3_pyc_lig.png"),
#        width = 20,
#        height = 16,
#        units = "cm"
# )

(fig_panels_marine2_manuscript <- (fig_mar_svcv + fig_mar_pvsbdv) +
    plot_annotation(tag_levels = 'A'))

# ggsave(("figures/Fig4_paneled_marine_deg.png"),
#        width = 20,
#        height = 8,
#        units = "cm"
# )

# For sake of figure clarity, going to focus on site and depth, but leave discussion
# of within core results to the results section of the text.

#### Results calculations ####
# Summary tables
# by environment (beach/ocean)
m <- dat_ed %>%
  group_by(Environment_f) %>%
  summarize(meanOC = mean(perc_OC, na.rm = TRUE),
            sdOC = sd(perc_OC, na.rm = TRUE),
            mean13C = mean(d13C, na.rm = TRUE),
            sd13C = sd(d13C, na.rm = TRUE),
            meanpyC = mean(pyC, na.rm = TRUE),
            sdpyC = sd(pyC, na.rm = TRUE),
            meanLAM = mean(Lambda, na.rm = TRUE),
            sdLAM = sd(Lambda, na.rm = TRUE),
            meanSV = mean(SV, na.rm = TRUE),
            sdSV = sd(SV, na.rm = TRUE),
            meanCV = mean(CV, na.rm = TRUE),
            sdCV = sd(CV, na.rm = TRUE),
            meanPVS = mean(PVS, na.rm = TRUE),
            sdPVS = sd(PVS, na.rm = TRUE),
            meanBDV = mean(BdV, na.rm = TRUE),
            sdBDV = sd(BdV, na.rm = TRUE))

# by beaches only (not including deposition site or low-tide line)
dlist <- c("DD", "MD")

m1 <- dat_ed %>%
  filter(Location_id %in% c("GBEA", "GOSL")) %>%
  group_by(Date) %>%
  summarize(meanOC = mean(perc_OC, na.rm = TRUE),
            sdOC = sd(perc_OC, na.rm = TRUE),
            meand13C = mean(d13C, na.rm = TRUE),
            sdd13C = sd(d13C, na.rm = TRUE),
            meanpyC = mean(pyC, na.rm = TRUE),
            sdpyC = sd(pyC, na.rm = TRUE),
            meanLAM = mean(Lambda, na.rm = TRUE),
            sdLAM = sd(Lambda, na.rm = TRUE),
            meanSV = mean(SV, na.rm = TRUE),
            sdSV = sd(SV, na.rm = TRUE),
            meanCV = mean(CV, na.rm = TRUE),
            sdCV = sd(CV, na.rm = TRUE),
            meanPVS = mean(PVS, na.rm = TRUE),
            sdPVS = sd(PVS, na.rm = TRUE),
            meanBDV = mean(BdV, na.rm = TRUE),
            sdBDV = sd(BdV, na.rm = TRUE))

# by beaches only (not including low-tide line)
m1.2 <- dat_ed %>%
  filter(Environment_f == "Beach") %>%
  mutate(Location_f = factor(case_when(Location_id %in% dlist ~ "DD",
                                       TRUE ~ Location_id))) %>%
  filter(Location_f != "CD") %>%
  group_by(Date) %>%
  summarize(meanOC = mean(perc_OC, na.rm = TRUE),
            sdOC = sd(perc_OC, na.rm = TRUE),
            meand13C = mean(d13C, na.rm = TRUE),
            sdd13C = sd(d13C, na.rm = TRUE),
            meanpyC = mean(pyC, na.rm = TRUE),
            sdpyC = sd(pyC, na.rm = TRUE),
            meanLAM = mean(Lambda, na.rm = TRUE),
            sdLAM = sd(Lambda, na.rm = TRUE),
            meanSV = mean(SV, na.rm = TRUE),
            sdSV = sd(SV, na.rm = TRUE),
            meanCV = mean(CV, na.rm = TRUE),
            sdCV = sd(CV, na.rm = TRUE),
            meanPVS = mean(PVS, na.rm = TRUE),
            sdPVS = sd(PVS, na.rm = TRUE),
            meanBDV = mean(BdV, na.rm = TRUE),
            sdBDV = sd(BdV, na.rm = TRUE))

# by marine sites only
m2 <- dat_ed %>%
  filter(Environment_f == "Ocean") %>%
  group_by(Location_id) %>%
  summarize(meanOC = mean(perc_OC, na.rm = TRUE),
            sdOC = sd(perc_OC, na.rm = TRUE),
            mean13C = mean(d13C, na.rm = TRUE),
            sd13C = sd(d13C, na.rm = TRUE),
            meanpyC = mean(pyC, na.rm = TRUE),
            sdpyC = sd(pyC, na.rm = TRUE),
            meanLAM = mean(Lambda, na.rm = TRUE),
            sdLAM = sd(Lambda, na.rm = TRUE),
            meanSV = mean(SV, na.rm = TRUE),
            sdSV = sd(SV, na.rm = TRUE),
            meanCV = mean(CV, na.rm = TRUE),
            sdCV = sd(CV, na.rm = TRUE),
            meanPVS = mean(PVS, na.rm = TRUE),
            sdPVS = sd(PVS, na.rm = TRUE),
            meanBDV = mean(BdV, na.rm = TRUE),
            sdBDV = sd(BdV, na.rm = TRUE))

# by water sampling depth at marine sites
m3 <- dat_ed %>%
  filter(Environment_f == "Ocean") %>%
  group_by(Water_Depth) %>%
  summarize(meanOC = mean(perc_OC, na.rm = TRUE),
            sdOC = sd(perc_OC, na.rm = TRUE),
            mean13C = mean(d13C, na.rm = TRUE),
            sd13C = sd(d13C, na.rm = TRUE),
            meanpyC = mean(pyC, na.rm = TRUE),
            sdpyC = sd(pyC, na.rm = TRUE),
            meanLAM = mean(Lambda, na.rm = TRUE),
            sdLAM = sd(Lambda, na.rm = TRUE),
            meanSV = mean(SV, na.rm = TRUE),
            sdSV = sd(SV, na.rm = TRUE),
            meanCV = mean(CV, na.rm = TRUE),
            sdCV = sd(CV, na.rm = TRUE),
            meanPVS = mean(PVS, na.rm = TRUE),
            sdPVS = sd(PVS, na.rm = TRUE),
            meanBDV = mean(BdV, na.rm = TRUE),
            sdBDV = sd(BdV, na.rm = TRUE))

# by section of core at marine sites
m4 <- dat_ed %>%
  filter(Environment_f == "Ocean") %>%
  group_by(Core_Section) %>%
  summarize(meanOC = mean(perc_OC, na.rm = TRUE),
            sdOC = sd(perc_OC, na.rm = TRUE),
            mean13C = mean(d13C, na.rm = TRUE),
            sd13C = sd(d13C, na.rm = TRUE),
            meanpyC = mean(pyC, na.rm = TRUE),
            sdpyC = sd(pyC, na.rm = TRUE),
            meanLAM = mean(Lambda, na.rm = TRUE),
            sdLAM = sd(Lambda, na.rm = TRUE),
            meanSV = mean(SV, na.rm = TRUE),
            sdSV = sd(SV, na.rm = TRUE),
            meanCV = mean(CV, na.rm = TRUE),
            sdCV = sd(CV, na.rm = TRUE),
            meanPVS = mean(PVS, na.rm = TRUE),
            sdPVS = sd(PVS, na.rm = TRUE),
            meanBDV = mean(BdV, na.rm = TRUE),
            sdBDV = sd(BdV, na.rm = TRUE))

# % OC range
min(dat_ed$perc_OC) # 0.05
max(dat_ed$perc_OC) # 3.59

# % OC figure
(fig7 <- ggplot(dat_plots) +
  geom_point(aes(x = Date, y = OC, color = loc_plot)) +
  scale_color_manual(name = "Sampling Location", 
                     values = cal_palette("fire", n = 6, type = "continuous")) +
  theme_bw())

# δ13C range
min(dat_ed$d13C) # -27.05
max(dat_ed$d13C) # -21.52

# d13C figure
(fig8 <- ggplot(dat_plots) +
    geom_point(aes(x = Date, y = D13C, color = loc_plot)) +
    scale_color_manual(name = "Sampling Location", 
                       values = cal_palette("fire", n = 6, type = "continuous")) +
    theme_bw())

# %PyC range
min(dat_ed$pyC) # 0.26
max(dat_ed$pyC) # 4.75

# pyC figure
(fig9 <- ggplot(dat_plots) +
    geom_point(aes(x = Date, y = PYC, color = loc_plot)) +
    scale_color_manual(name = "Sampling Location", 
                       values = cal_palette("fire", n = 6, type = "continuous")) +
    theme_bw())

# lambda range
min(dat_ed$Lambda, na.rm = TRUE) # 0.03
max(dat_ed$Lambda, na.rm = TRUE) # 4.6

# lambda figure
(fig10 <- ggplot(dat_plots) +
    geom_point(aes(x = Date, y = LAM, color = loc_plot)) +
    scale_color_manual(name = "Sampling Location", 
                       values = cal_palette("fire", n = 6, type = "continuous")) +
    theme_bw())

# sv & cv ranges
marine <- dat_ed %>% filter(Environment_f == "Ocean")
min(marine$SV, na.rm = TRUE) # 0.93
max(marine$SV, na.rm = TRUE) # 3.88

min(marine$CV, na.rm = TRUE) # 0.05
max(marine$CV, na.rm = TRUE) # 0.48

# pvs & bdv ranges
min(marine$PVS, na.rm = TRUE) # 0.04
max(marine$PVS, na.rm = TRUE) # 0.44

min(marine$BdV, na.rm = TRUE) # 0.11
max(marine$BdV, na.rm = TRUE) # 1.35

# End of script.
