#################################################
## Paper Figure 2 ###############################
#################################################

## Here analysis and Plotting

## SISAL Karst Plot

#################################################
PLOTTING_VARIABLES <- list()
PLOTTING_VARIABLES$WIDTH = 8.3
PLOTTING_VARIABLES$HEIGHT = 5.5

source("Functions/Plotting/karst_map_plot.R")
library(tidyverse)

#aus 1_5_create_dataset_SISAL.R --> site_past1000_min1

#ALL_SITES <- read.csv("/stacywork/ginnyweasley/02_SISAL/SISAL_v2_CARLA/site_countries.csv")
ALL_SITES <- read.csv("/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/02_SISAL/SISAL_v2/site_countries.csv")

sites_spec <- DATA_past1000$CAVES$entity_info$site_id[mask_spec]
sites_var  <- DATA_past1000$CAVES$entity_info$site_id[mask_var]
sites_mean <- DATA_past1000$CAVES$entity_info$site_id[mask_mean]

USED_SITES_spec <- ALL_SITES %>% filter(site_id %in% DATA_past1000$CAVES$site_info$site_id) %>% filter(site_id %in% sites_spec) %>% distinct(site_id, longitude, latitude)
USED_SITES_spec <- data.frame(
  lon = USED_SITES_spec$longitude,
  lat = USED_SITES_spec$latitude,
  value = USED_SITES_spec$site_id
)
USED_SITES_var <- ALL_SITES %>% filter(site_id %in% DATA_past1000$CAVES$site_info$site_id) %>% filter(!site_id %in% sites_spec) %>% 
  filter(site_id %in% sites_var) %>% distinct(site_id, longitude, latitude)
USED_SITES_var <- data.frame(
  lon = USED_SITES_var$longitude,
  lat = USED_SITES_var$latitude,
  value = USED_SITES_var$site_id
)
USED_SITES_mean <- ALL_SITES %>% filter(site_id %in% DATA_past1000$CAVES$site_info$site_id) %>% filter(!site_id %in% sites_spec) %>% 
  filter(!site_id %in% sites_var) %>% 
  filter(site_id %in% sites_mean) %>% distinct(site_id, longitude, latitude)
USED_SITES_mean <- data.frame(
  lon = USED_SITES_mean$longitude,
  lat = USED_SITES_mean$latitude,
  value = USED_SITES_mean$site_id
)
NOT_SITES <- ALL_SITES %>% filter(!site_id %in% sites_mean) %>% distinct(site_id, longitude, latitude)
NOT_SITES <- data.frame(
  lon = NOT_SITES$longitude,
  lat = NOT_SITES$latitude,
  value = NOT_SITES$site_id
)

#################################################
## PLOTTING #####################################

GLOBAL_STACY_OPTIONS$GLOBAL_FONT_SIZE = 10
plot <- karst_map_plot(USED_SITES_spec = USED_SITES_spec,
                       USED_SITES_var = USED_SITES_var,
                       USED_SITES_mean = USED_SITES_mean,
                       NOT_SITES = NOT_SITES, pt_size = 2.5) + 
  #theme(legend.position = "bottom", legend.direction = "horizontal")#c(0.0, 0.02))
  theme(legend.position = c(-0.01, 0), legend.justification = c(0, 0), legend.box = 'vertical',
        axis.text = element_blank())
#plot
plot <- plot + theme(panel.border = element_blank()) #,legend.text = element_text(size = 8)

#plot
plot %>% ggsave(filename = paste('ECHAM_Plot_2_SISAL_database', 'pdf', sep = '.'), plot = ., path = 'Plots', 
                width = 2*PLOTTING_VARIABLES$WIDTH, height = 2*PLOTTING_VARIABLES$HEIGHT-2, units = 'cm', dpi = 'print', device = "pdf")

remove(plot, NOT_SITES, USED_SITES_mean, USED_SITES_spec, USED_SITES_var, ALL_SITES, sites_mean, sites_spec, sites_var)
remove(karst_map)
