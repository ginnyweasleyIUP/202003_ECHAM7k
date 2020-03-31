#################################################
## Paper Figure 1 ###############################
#################################################

## Here analysis and Plotting

## Einführungsplot (GENERAL --> von Kira kopieren)

library(dplyr)
library(latex2exp)
source("Functions/Plotting/STACYmap_5.R")
source("Functions/Plotting/STACYmap_5_1_NAgrid.R")
source("Functions/projection_ptlyr.R")

# this is Regression with downsampled temp and prec

#################################################
## PLOTS ########################################
#################################################

r.squared = 0.5

Plot_lyr_temp <- ANALYSIS$REG$FIELD$REG_TEMP_ISOT
Plot_lyr_temp_r <- ANALYSIS$REG$FIELD$REG_TEMP_ISOT_R
Plot_lyr_temp[Plot_lyr_temp_r > r.squared] <- NA
Plot_lyr_temp[abs(Plot_lyr_temp)>2.5] = NA
Plot_lyr_prec <- ANALYSIS$REG$FIELD$REG_PREC_ISOT
Plot_lyr_prec_r <- ANALYSIS$REG$FIELD$REG_PREC_ISOT_R
Plot_lyr_prec[Plot_lyr_prec_r > r.squared] <- NA
Plot_lyr_prec[abs(Plot_lyr_prec)>2] = NA

Plot_lyr_temp <- rbind(Plot_lyr_temp[49:96,1:48],
                       Plot_lyr_temp[1:48,1:48])
Plot_lyr_prec <- rbind(Plot_lyr_prec[49:96,1:48],
                       Plot_lyr_prec[1:48,1:48])

##### Point Layer

Point_Lyr_temp <- list(lon = list(), lat = list(), value = list())
Point_Lyr_prec <- list(lon = list(), lat = list(), value = list())

length_cave = length(DATA_past7000$CAVES$entity_info$site_id)

for(ii in 1:length_cave){
  site <- DATA_past7000$CAVES$entity_info$site_id[ii]
  print(ii)
  if(!mask_mean[ii]){next}
  # 1) sortiert aus, was nicht signifikant ist
  if(!is.na(ANALYSIS$REG$POINTS$p_TEMP[ii]) & ANALYSIS$REG$POINTS$p_TEMP[ii] > r.squared){
  #if(!is.na(ANALYSIS$REG$POINTS$p_TEMP[ii])){
    Point_Lyr_temp$lon = c(Point_Lyr_temp$lon, DATA_past7000$CAVES$site_info$longitude[DATA_past7000$CAVES$site_info$site_id == site])
    Point_Lyr_temp$lat = c(Point_Lyr_temp$lat, DATA_past7000$CAVES$site_info$latitude[DATA_past7000$CAVES$site_info$site_id == site])
    Point_Lyr_temp$value = c(Point_Lyr_temp$value, ANALYSIS$REG$POINTS$REG_TEMP[ii])
    # 2) betrachte signifikante Korrelationen:
  }
  if(!is.na(ANALYSIS$REG$POINTS$p_PREC[ii]) & ANALYSIS$REG$POINTS$p_PREC[ii] > r.squared){
    Point_Lyr_prec$lon = c(Point_Lyr_prec$lon, DATA_past7000$CAVES$site_info$longitude[DATA_past7000$CAVES$site_info$site_id == site])
    Point_Lyr_prec$lat = c(Point_Lyr_prec$lat, DATA_past7000$CAVES$site_info$latitude[DATA_past7000$CAVES$site_info$site_id == site])
    Point_Lyr_prec$value = c(Point_Lyr_prec$value, ANALYSIS$REG$POINTS$REG_PREC[ii])
    # 2) betrachte signifikante Korrelationen:
  }
}



Point_Lyr_temp$lon = as.numeric(Point_Lyr_temp$lon)
Point_Lyr_temp$lat = as.numeric(Point_Lyr_temp$lat)
Point_Lyr_temp$value = as.numeric(Point_Lyr_temp$value)

Point_Lyr_prec$lon = as.numeric(Point_Lyr_prec$lon)
Point_Lyr_prec$lat = as.numeric(Point_Lyr_prec$lat)
Point_Lyr_prec$value = as.numeric(Point_Lyr_prec$value)

mask_china <- logical(length(Point_Lyr_temp$lon))

for(ii in 1:length(Point_Lyr_temp$lon)){
  if(is.na(Point_Lyr_temp$lon[ii])){next}
  if(Point_Lyr_temp$lon[ii] > 100 & Point_Lyr_temp$lon[ii] < 120){
    if(Point_Lyr_temp$lat[ii] < 35 & Point_Lyr_temp$lat[ii] > 22){
      mask_china[ii] = T}
  }
}

ptlyr_china_temp <- data.frame(lon = Point_Lyr_temp$lon[mask_china], lat = Point_Lyr_temp$lat[mask_china], value = Point_Lyr_temp$value[mask_china])
ptlyr_china_prec <- data.frame(lon = Point_Lyr_prec$lon[mask_china], lat = Point_Lyr_prec$lat[mask_china], value = Point_Lyr_prec$value[mask_china])

ptlyr_rest_temp <- data.frame(lon = Point_Lyr_temp$lon[!mask_china], lat = Point_Lyr_temp$lat[!mask_china], value = Point_Lyr_temp$value[!mask_china])
ptlyr_rest_prec <- data.frame(lon = Point_Lyr_prec$lon[!mask_china], lat = Point_Lyr_prec$lat[!mask_china], value = Point_Lyr_prec$value[!mask_china])

ptlyr_china_temp_p <- projection_ptlyr(ptlyr_china_temp, projection = as.character('+proj=robin +datum=WGS84'))
ptlyr_rest_temp_p <- projection_ptlyr(ptlyr_rest_temp, projection = as.character('+proj=robin +datum=WGS84'))
ptlyr_china_prec_p <- projection_ptlyr(ptlyr_china_prec, projection = as.character('+proj=robin +datum=WGS84'))
ptlyr_rest_prec_p <- projection_ptlyr(ptlyr_rest_prec, projection = as.character('+proj=robin +datum=WGS84'))

GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE <- 3

NA_plot_lyr <- Plot_lyr_temp
NA_plot_lyr[!is.na(NA_plot_lyr)] <- 0
NA_plot_lyr[is.na(NA_plot_lyr)] <- 1


plot_temp <- STACYmap_NA(gridlyr = Plot_lyr_temp, centercolor = 0, graticules = T,
                         NA_gridlyr = NA_plot_lyr, NA_color = "grey",
                         legend_names = list(grid = TeX("$Reg(T, \\delta^{18}O), r^2>0.5$"))) +
  geom_point(data = ptlyr_china_temp_p, aes(x = long, y = lat, fill = layer), shape = 21, alpha = 0.7, color = "black",
             size = (GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE), show.legend = c(color =TRUE), position = position_jitter(width = 1000000, height = 500000))+
  geom_point(data = ptlyr_rest_temp_p, aes(x = long, y = lat, fill = layer), shape = 21, alpha = 0.7, color = "black",
             size = (GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE), show.legend = c(color =TRUE)) +
  theme(panel.border = element_blank(),
        legend.background = element_blank(),
        axis.text = element_blank(),
        text = element_text(size = 12),
        legend.title = element_text(size = 12))

plot_temp

NA_plot_lyr <- Plot_lyr_prec
NA_plot_lyr[!is.na(NA_plot_lyr)] <- 0
NA_plot_lyr[is.na(NA_plot_lyr)] <- 1

plot_prec <- STACYmap_NA(gridlyr = Plot_lyr_prec, centercolor = 0, graticules = T,
                         NA_gridlyr = NA_plot_lyr, NA_color = "grey",
                         legend_names = list(grid = TeX("$Reg(P, \\delta^{18}O), r^2>0.5$"))) +
  geom_point(data = ptlyr_china_prec_p, aes(x = long, y = lat, fill = layer), shape = 21, alpha = 0.7, color = "black",
             size = (GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE), show.legend = c(color =TRUE), position = position_jitter(width = 1000000, height = 500000))+
  geom_point(data = ptlyr_rest_prec_p, aes(x = long, y = lat, fill = layer), shape = 21, alpha = 0.7, color = "black",
             size = (GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE), show.legend = c(color =TRUE)) +
  theme(panel.border = element_blank(),
        legend.background = element_blank(),
        axis.text = element_blank(),
        text = element_text(size = 12),
        legend.title = element_text(size = 12))

#plot_prec

library(ggpubr)
plot <- ggarrange(plot_temp, plot_prec,
                  labels = c("(a)", "(b)"),
                  ncol = 2, nrow = 1)

plot  %>% ggsave(filename = paste0('ECHAM7k_Plot_5_Regression.pdf'), plot = ., path = 'Plots', 
                 width = 2*12, height = 12/8.3*5, units = 'cm', dpi = 'print', device = "pdf")
plot  %>% ggsave(filename = paste0('ECHAM7k_Plot_5_Regression.png'), plot = ., path = 'Plots', 
                 width = 2*12, height = 12/8.3*5, units = 'cm', dpi = 'print', device = "png")

remove(COR, double_time, plot, Plot_lyr_prec, Plot_lyr_prec_p, Plot_lyr_temp, Plot_lyr_temp_p, plot_prec, plot_temp, Point_Lyr_prec, Point_Lyr_temp, s, test)
remove(allmax, allmax_real, diff_dt, entity, ii, length_cave, record, run, sim, site, var)
