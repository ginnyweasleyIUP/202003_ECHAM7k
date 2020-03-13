#################################################
## Paper Figure 3 ###############################
#################################################

## Here analysis and Plotting

## Mean State

## DATA Acquisition #############################

# mask_mean = logical(length = length(DATA_past1000$CAVES$entity_info$entity_id))
# mask_var  = logical(length = length(DATA_past1000$CAVES$entity_info$entity_id))
# mask_spec = logical(length = length(DATA_past1000$CAVES$entity_info$entity_id))
# 
# for(entity in 1:length(DATA_past1000$CAVES$entity_info$entity_id)){
#   if(DATA_past1000$CAVES$entity_info$n[entity] > 10 & DATA_past1000$CAVES$entity_info$period[entity] > 600){mask_mean[entity] = T}
#   if(DATA_past1000$CAVES$entity_info$n[entity] > 20 & DATA_past1000$CAVES$entity_info$period[entity] > 600){mask_var[entity] = T}
#   if(DATA_past1000$CAVES$entity_info$n[entity] > 30 & DATA_past1000$CAVES$entity_info$period[entity] > 600){mask_spec[entity] = T}
# }
# 
PLOTTING_DATA <- list()
PLOTTING_DATA$FIG3 <- list()

#################################################

library(plyr)
library(dplyr)
library(rgdal)

## PLOTTING #####################################

source("Functions/Plotting/STACYmap_5.R")
#source("Functions/Plotting/STACYmap_5_1_NAgrid.R")
source("Functions/Plotting/STACYmap_5_2_logscale.R")
source("Functions/aw_mean.R")

GLOBAL_STACY_OPTIONS$GLOBAL_FONT_SIZE = 10
color_slpr = c("#edf7fb", "#ddedf4", "#cfe2ef", "#c1d6e8", "#b3cde2", "#a8bfdb", "#9eb1d4", "#95a3cd",
               "#8c95c6", "#8a85bd", "#8975b5", "#8966ae", "#8756a7", "#85449c", "#811e84", "#800e7c")

PLOTTING_VARIABLES$COLORS$SLPR <- color_slpr
remove(color_slpr)

for(run in c("a")){
  for(var in c("ISOT", "ITPC")){
    # TEMP LAYER
    temp_lyr <- apply(DATA_past1000[[paste0("SIM_yearly_",run)]]$TEMP, c(1,2), mean, na.rm = T)
    temp_lyr_hadcm3 <- apply(DATA_past1000_HadCM3[[paste0("SIM_yearly_",run)]]$TEMP, c(1,2), mean, na.rm = T)
    plot_temp <- STACYmap(gridlyr = rbind(temp_lyr[49:96,1:48],temp_lyr[1:48,1:48]),
                          #zoom = c(-180, -60, 180, 48),
                          legend_names = list(grid = "Temperature (°C)"),
                          graticules = TRUE,
                          colorscheme = "temp", 
                          centercolor = 0, 
                          allmax_forced = c(-max(abs(temp_lyr_hadcm3)), max(abs(temp_lyr_hadcm3)))) +
      theme(panel.border = element_blank(),
            legend.background = element_blank(),
            axis.text = element_blank(),
            legend.text = element_text(size = 8)) 
    
    #plot_temp
    # PREC LAYER
    prec_lyr <- apply(DATA_past1000[[paste0("SIM_yearly_",run)]]$PREC, c(1,2), mean, na.rm = T)
    prec_lyr_hadcm3 <- 8.6148e4*apply(DATA_past1000_HadCM3[[paste0("SIM_yearly_",run)]]$PREC, c(1,2), mean, na.rm = T)
    plot_prec <- STACYmap(gridlyr = rbind(prec_lyr[49:96,1:48],prec_lyr[1:48,1:48]),
                          #zoom = c(-180, -60, 180, 48),
                          legend_names = list(grid = "Precipitation (mm/day)"),
                          graticules = TRUE,
                          colorscheme = "prcp_grd",
                          allmax_forced = c(0, max(abs(prec_lyr_hadcm3)))) +
      theme(panel.border = element_blank(),
            legend.background = element_blank(),
            axis.text = element_blank(),
            legend.text = element_text(size = 8)) 
    
    plot_prec
    
    # ISOT POINTS
    PLOTTING_DATA$FIG3$CAVElyr_isot <- data.frame(
      lon = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
      lat = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)), 
      value = numeric(length(DATA_past1000$CAVES$entity_info$entity_id))
    )
    
    
    for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id)){
      site = DATA_past1000$CAVES$entity_info$site_id[ii]
      entity = DATA_past1000$CAVES$entity_info$entity_id[ii]
      PLOTTING_DATA$FIG3$CAVElyr_isot$lon[ii]   = DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site]
      PLOTTING_DATA$FIG3$CAVElyr_isot$lat[ii]   = DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site]
      PLOTTING_DATA$FIG3$CAVElyr_isot$value[ii] = mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]][[paste0("d18O_dw_eq_",run)]], na.rm = T)
    }
    
    PLOTTING_DATA$FIG3$CAVElyr_isot_used <- data.frame(
      lon = PLOTTING_DATA$FIG3$CAVElyr_isot$lon[mask_mean],
      lat = PLOTTING_DATA$FIG3$CAVElyr_isot$lat[mask_mean],
      value = PLOTTING_DATA$FIG3$CAVElyr_isot$value[mask_mean]
    )
    
    remove(entity, site, ii)
    
    #ISOT LAYER
    isot_lyr <- apply(DATA_past1000[[paste0("SIM_yearly_",run)]][[var]], c(1,2), mean, na.rm = T)
    Plot_lyr1 <-rbind(isot_lyr[49:96,1:48],isot_lyr[1:48,1:48])
    
    Plot_lyr3 <- Plot_lyr1
    Plot_lyr3[is.na(Plot_lyr3)] = 1000
    Plot_lyr3[Plot_lyr3>0] <- Plot_lyr3[Plot_lyr3>0]+1 
    Plot_lyr3[Plot_lyr3<0] <- Plot_lyr3[Plot_lyr3<0]-1
    Plot_lyr3[Plot_lyr3>0] <- log10(Plot_lyr3[Plot_lyr3>0])
    Plot_lyr3[Plot_lyr3<0] <- - log10(abs(Plot_lyr3[Plot_lyr3<0]))
    Plot_lyr3[abs(Plot_lyr3)>5] <- NA
    Plot_lyr3[,1] <- NA
    Plot_lyr3[,48] <- NA
    #Plot_lyr3[,1:6] <- NA
    #Plot_lyr3[,60:48] <- NA
    #Plot_lyr1[,1:6] <- NA
    #Plot_lyr1[,60:48] <- NA
    
    
    Point_Lyr <- data.frame(
      lon = PLOTTING_DATA$FIG3$CAVElyr_isot_used$lon,
      lat = PLOTTING_DATA$FIG3$CAVElyr_isot_used$lat,
      value = - log10(abs(PLOTTING_DATA$FIG3$CAVElyr_isot_used$value -1))
    )
    
    Point_Lyr$value[[57]] <- log10(PLOTTING_DATA$FIG3$CAVElyr_isot_used$value[[57]]+1)
    
    # Plot_lyr1[Plot_lyr1 < -17] = -17 -0.1*Plot_lyr1
    # Plot_lyr2 <-rbind(DATA_past1000$SIM_mean$ISOT[49:96,1:48],DATA_past1000$SIM_mean$ISOT[1:48,1:48])
    # Plot_lyr2[Plot_lyr2 > -17] = NA
    # Plot_lyr2[Plot_lyr2 < -17] = -17
    
    allmax = - min(Plot_lyr3, na.rm = T)
    allmax_real = - min(Plot_lyr1, na.rm = T)
    
    GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE = 4
    
    leg_name = expression(paste(delta^{plain(18)}, plain(O), " (%)"))
    if(var == "ITPC"){leg_name = expression(paste(delta^{plain(18)}, plain(O)[plain(pw)], " (%)"))}
    
    plot_isot <- STACYmap_isot(gridlyr = Plot_lyr3,
                               ptlyr = Point_Lyr,
                               legend_names = list(grid =leg_name),
                               graticules = TRUE,
                               colorscheme = rev(RColorBrewer::brewer.pal(11, 'BrBG'))[c(1:8,10)],
                               #colorscheme = PLOTTING_VARIABLES$COLORS$SLPR,
                               breaks_isot = c(-log10(71),-log10(11),-log10(6), -log10(2), 0, log10(2), log10(6)),
                               labels_isot = c(-70, -10,-5, -1, 0, 1, 5),
                               min_grid = -log10(71), max_grid = log10(12)) +
      theme(panel.border = element_blank(),
            legend.background = element_blank(),
            axis.text = element_blank(),
            text = element_text(size = 8)) 
    
    plot_isot
    
    # PLOT MEAN BIAS
    
    PLOTTING_DATA$FIG3$CAVElyr_diff <- data.frame(
      lon = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
      lat = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)), 
      value = numeric(length(DATA_past1000$CAVES$entity_info$entity_id))
    )
    
    
    for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id)){
      site = DATA_past1000$CAVES$entity_info$site_id[ii]
      entity = DATA_past1000$CAVES$entity_info$entity_id[ii]
      PLOTTING_DATA$FIG3$CAVElyr_diff$lon[ii]   = DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site]
      PLOTTING_DATA$FIG3$CAVElyr_diff$lat[ii]   = DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site]
      PLOTTING_DATA$FIG3$CAVElyr_diff$value[ii] = mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]][[paste0(var,"_",run)]], na.rm = T) - 
        mean(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]][[paste0("d18O_dw_eq_",run)]], na.rm = T)
    }
    
    PLOTTING_DATA$FIG3$CAVElyr_diff_used <- data.frame(
      lon = PLOTTING_DATA$FIG3$CAVElyr_diff$lon[mask_mean],
      lat = PLOTTING_DATA$FIG3$CAVElyr_diff$lat[mask_mean],
      value = PLOTTING_DATA$FIG3$CAVElyr_diff$value[mask_mean]
    )
    PLOTTING_DATA$FIG3$CAVElyr_diff_used$value[1] = NA
    PLOTTING_DATA$FIG3$CAVElyr_diff_used$value[105] = NA
    PLOTTING_DATA$FIG3$CAVElyr_diff_used$value[103] = NA
    PLOTTING_DATA$FIG3$CAVElyr_diff_used$value[104] = NA
    
    ptlyr <- projection_ptlyr(PLOTTING_DATA$FIG3$CAVElyr_diff_used, as.character('+proj=robin +datum=WGS84'))
    
    leg_name = expression(paste(delta^{plain(18)}, plain(O)[plain(sim)]," - ", delta^{plain(18)}, plain(O)[plain(rec)], " (%)"))
    if(var == "ITPC"){leg_name = expression(paste(delta^{plain(18)}, plain(O)["sim,pw"]," - ", delta^{plain(18)}, plain(O)[plain(rec)], " (%)"))}
    
    plot_diff <- STACYmap(coastline = TRUE) +
      geom_point(data = ptlyr, aes(x = long, y = lat, fill = layer), shape = 21, alpha = 0.7, color = "black",
                 size = 4, show.legend = c(color =TRUE)) +
      scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, 'RdBu')), 
                           limits = c(-max(abs(ptlyr$layer)), max(abs(ptlyr$layer))),
                           name = leg_name, 
                           guide = guide_colorbar(barwidth = 10, barheight = 0.3)) +
      theme(legend.direction = "horizontal", 
            panel.border = element_blank(),
            legend.background = element_blank(),
            axis.text = element_blank(),
            text = element_text(size = 10),
            legend.title = element_text(size = 10))
    
    plot_diff

    
    library(ggpubr)
    plot <- ggarrange(plot_temp, plot_prec, plot_isot, plot_diff,
                      labels = c("(a)", "(b)", "(c)", "(d)"),
                      ncol = 2, nrow = 2)
    
    plot  %>% ggsave(filename = paste0('ECHAM_Plot_3_Mean_',var,'_',run, '.pdf'), plot = ., path = 'Plots', 
                     width = 2*12, height = 2*12/8.3*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "pdf")
    plot  %>% ggsave(filename = paste0('ECHAM_Plot_3_Mean_',var,'_',run, '.png'), plot = ., path = 'Plots', 
                     width = 2*12, height = 2*12/8.3*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "png")
  }
}

# temp <- simpleawmean(DATA_past1000$SIM_mean$TEMP, seq(from = -90, to = 90, length.out = 48))
# temp_sd <- simpleawsd(DATA_past1000$SIM_mean$TEMP, seq(from = -90, to = 90, length.out = 48))
# prec <- simpleawmean(DATA_past1000$SIM_mean$PREC, seq(from = -90, to = 90, length.out = 48))*8.6148e4
# prec_sd <- simpleawsd(DATA_past1000$SIM_mean$PREC, seq(from = -90, to = 90, length.out = 48))*8.6148e4
# 
# isot <- simpleawmean(DATA_past1000$SIM_mean$ISOT, seq(from = -90, to = 90, length.out = 48))
# isot_sd <- simpleawsd(DATA_past1000$SIM_mean$ISOT, seq(from = -90, to = 90, length.out = 48))
# slpr <- simpleawmean(DATA_past1000$SIM_mean$SLPR, seq(from = -90, to = 90, length.out = 48))
# slpr_sd <- simpleawsd(DATA_past1000$SIM_mean$SLPR, seq(from = -90, to = 90, length.out = 48))




remove(plot_temp, plot_prec, plot_isot, plot_diff, leg_name, Plot_lyr1, Plot_lyr3, plot, Point_Lyr)
remove(temp_lyr, prec_lyr)