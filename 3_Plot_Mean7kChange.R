#################################################
## TREND COMPARE ################################
#################################################

ANALYSIS$MEAN <- list()
source("Functions/aw_mean.R")

for(var in c("TEMP","PREC","ISOT","ITPC")){
  value <- list("echam" = numeric(abs(diff(DATA_past7000$time))))
  
  for(ii in 1:length(value$echam)){
    value$echam[ii] <- simpleawmean(DATA_past7000$SIM_yearly[[var]][,,ii], lats = seq(from = -90, to = 90, length.out = 48))
  }
  
  #mean(value$a)*500000/(4419*mean(value$a)-500000) Temp Änderung pro Delta Änderung
  
  
  ANALYSIS$MEAN$TS[[var]] <- zoo(x = value$echam, order.by = seq(from = DATA_past7000$time[1], to = DATA_past7000$time[2], by = -1))
}

plot(ANALYSIS$MEAN$TS$ITPC)

ANALYSIS$MEAN$FIRST30 <- list(TEMP = array(dim = c(96,48)),
                              PREC = array(dim = c(96,48)),
                              ISOT = array(dim = c(96,48)),
                              ITPC = array(dim = c(96,48)))
ANALYSIS$MEAN$LAST30 <- list(TEMP = array(dim = c(96,48)),
                             PREC = array(dim = c(96,48)),
                             ISOT = array(dim = c(96,48)),
                             ITPC = array(dim = c(96,48)))

for(lon in 1:96){
  print(lon)
  for(lat in 1:48){
    for(var in c("TEMP", "PREC", "ISOT", "ITPC")){
      ANALYSIS$MEAN$FIRST30[[var]][lon,lat] = mean(DATA_past7000$SIM_yearly[[var]][lon,lat,1:100], na.rm = T)
      ANALYSIS$MEAN$LAST30[[var]][lon,lat] = mean(DATA_past7000$SIM_yearly[[var]][lon,lat,6897:6997], na.rm = T)
    }
  }
}

source("Functions/Plotting/STACYmap_5.R")

ANALYSIS$MEAN$POINTS <- data.frame(
  lon = numeric(length(DATA_past7000$CAVES$entity_info$entity_id[mask_var])),
  lat = numeric(length(DATA_past7000$CAVES$entity_info$entity_id[mask_var])),
  value = numeric(length(DATA_past7000$CAVES$entity_info$entity_id[mask_var]))
)

for(ii in 1:length(DATA_past7000$CAVES$entity_info$entity_id[mask_var])){
  entity = DATA_past7000$CAVES$entity_info$entity_id[mask_var][ii]
  site = DATA_past7000$CAVES$entity_info$site_id[DATA_past7000$CAVES$entity_info$entity_id == entity]
  ANALYSIS$MEAN$POINTS$lon[ii] <- DATA_past7000$CAVES$site_info$longitude[DATA_past7000$CAVES$site_info$site_id == site]
  ANALYSIS$MEAN$POINTS$lat[ii] <- DATA_past7000$CAVES$site_info$latitude[DATA_past7000$CAVES$site_info$site_id == site]
  
  s <- DATA_past7000$CAVES$record_data[[paste0("ENTITY", entity)]]
  year_start = PaleoSpec::FirstElement(s$interp_age)
  year_stop = PaleoSpec::LastElement(s$interp_age)
  s_first100 <- s %>% filter(interp_age > year_stop-100)
  s_last100 <- s %>% filter(interp_age < year_start+100)
  
  ANALYSIS$MEAN$POINTS$value[ii] <- mean(s_last100$d18O_dw_eq, na.rm = T)-mean(s_first100$d18O_dw_eq, na.rm = T)
}


point_lyr <- ANALYSIS$MEAN$POINTS
plot_lyr <- ANALYSIS$MEAN$LAST30$ITPC-ANALYSIS$MEAN$FIRST30$ITPC
plot_lyr[plot_lyr>2.0] = max(plot_lyr[plot_lyr<2.0], na.rm = T)

plot <- STACYmap(gridlyr = rbind(plot_lyr[49:96,1:48],plot_lyr[1:48,1:48]),
                 ptlyr = point_lyr,
                 centercolor = 0,
                 legend_names = list(grid = "d18O change 0BP-7000BP (%)"))
plot %>% ggsave(filename = paste0('ECHAM7k_MeanChange_100y_.pdf'), plot = ., path = 'Plots', 
                width = 2*12, height = 2*12/8.3*5, units = 'cm', dpi = 'print', device = "pdf")


#################################################
## TRY same with a trend analysis ###############
#################################################

ANALYSIS$MEAN$TREND <- list(TEMP = array(dim = c(96,48)),
                            PREC = array(dim = c(96,48)),
                            ISOT = array(dim = c(96,48)),
                            ITPC = array(dim = c(96,48)))

for(lon in 1:96){
  for(lat in 1:48){
    for(var in c("TEMP", "PREC", "ISOT","ITPC")){
      ANALYSIS$MEAN$TREND[[var]][lon,lat] <- lm(DATA_past7000$SIM_yearly[[var]][lon,lat,] ~ c(6997:1))$coefficients[[2]]
    }
  }
}

ANALYSIS$MEAN$TREND$CAVES <- data.frame(
  lon = numeric(length(DATA_past7000$CAVES$entity_info$entity_id[mask_var])),
  lat = numeric(length(DATA_past7000$CAVES$entity_info$entity_id[mask_var])),
  value = numeric(length(DATA_past7000$CAVES$entity_info$entity_id[mask_var])),
  value_sim = numeric(length(DATA_past7000$CAVES$entity_info$entity_id[mask_var]))
)

for(ii in 1:length(DATA_past7000$CAVES$entity_info$entity_id[mask_var])){
  entity = DATA_past7000$CAVES$entity_info$entity_id[mask_var][ii]
  site = DATA_past7000$CAVES$entity_info$site_id[DATA_past7000$CAVES$entity_info$entity_id == entity]
  ANALYSIS$MEAN$TREND$CAVES$lon[ii] <- DATA_past7000$CAVES$site_info$longitude[DATA_past7000$CAVES$site_info$site_id == site]
  ANALYSIS$MEAN$TREND$CAVES$lat[ii] <- DATA_past7000$CAVES$site_info$latitude[DATA_past7000$CAVES$site_info$site_id == site]
  ANALYSIS$MEAN$TREND$CAVES$value[ii] <- -1000*lm(DATA_past7000$CAVES$record_data[[paste0("ENTITY",entity)]]$d18O_measurement ~ DATA_past7000$CAVES$record_data[[paste0("ENTITY",entity)]]$interp_age)$coefficients[[2]]
  ANALYSIS$MEAN$TREND$CAVES$value_sim[ii] <- -1000*lm(DATA_past7000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$ITPC ~ DATA_past7000$CAVES$sim_data_downsampled[[paste0("ENTITY",entity)]]$interp_age)$coefficients[[2]]
}



plot_lyr <- -1000*ANALYSIS$MEAN$TREND$ITPC
plot_lyr[plot_lyr>0.5] = max(plot_lyr[plot_lyr<0.5], na.rm = T)
point_lyr <- ANALYSIS$MEAN$TREND$CAVES[1:3]


plot <- STACYmap(gridlyr = rbind(plot_lyr[49:96,1:48],plot_lyr[1:48,1:48]),
                 ptlyr = point_lyr,
                 centercolor = 0,
                 legend_names = list(grid = "d18O/ky since 7000BP"))+
  theme(panel.border = element_blank(),
        legend.background = element_blank(),
        axis.text = element_blank(),
        legend.text = element_text(size = 8)) 

plot %>% ggsave(filename = paste0('ECHAM7k_Trend_ITPC.pdf'), plot = ., path = 'Plots', 
                width = 12, height = 12/8.3*5, units = 'cm', dpi = 'print', device = "pdf")



#################################################

colZ <- RColorBrewer::brewer.pal(12, "Set3")
range_d18O <- range(c(as.numeric(ANALYSIS$MEAN$TS$ITPC)), na.rm = T)
range_d18O <- c(-10,0)
# Start plot
plot(xlimz,range_d18O,type="n",xlab="",ylab="")

lines(ANALYSIS$MEAN$TS$ITPC,type="l",col=adjustcolor(col_hadcm3_temp,0.3))
lines(gaussdetr(ANALYSIS$MEAN$TS$ITPC,tsc.in=100)$Xsmooth,type="l",col=col_hadcm3_temp)

lines(DATA_past7000$CAVES$record_data$ENTITY20$interp_age, DATA_past7000$CAVES$record_data$ENTITY20$d18O_dw_eq,type="l",col=colZ[1])
lines(DATA_past7000$CAVES$record_data$ENTITY27$interp_age, DATA_past7000$CAVES$record_data$ENTITY27$d18O_dw_eq,type="l",col=colZ[2])
lines(DATA_past7000$CAVES$record_data$ENTITY32$interp_age, DATA_past7000$CAVES$record_data$ENTITY32$d18O_dw_eq,type="l",col=colZ[3])
lines(DATA_past7000$CAVES$record_data$ENTITY48$interp_age, DATA_past7000$CAVES$record_data$ENTITY48$d18O_dw_eq,type="l",col=colZ[4])
lines(DATA_past7000$CAVES$record_data$ENTITY51$interp_age, DATA_past7000$CAVES$record_data$ENTITY51$d18O_dw_eq,type="l",col=colZ[5])
lines(DATA_past7000$CAVES$record_data$ENTITY52$interp_age, DATA_past7000$CAVES$record_data$ENTITY52$d18O_dw_eq,type="l",col=colZ[6])
lines(DATA_past7000$CAVES$record_data$ENTITY54$interp_age, DATA_past7000$CAVES$record_data$ENTITY54$d18O_dw_eq,type="l",col=colZ[7])
lines(DATA_past7000$CAVES$record_data$ENTITY58$interp_age, DATA_past7000$CAVES$record_data$ENTITY58$d18O_dw_eq,type="l",col=colZ[8])
lines(DATA_past7000$CAVES$record_data$ENTITY90$interp_age, DATA_past7000$CAVES$record_data$ENTITY90$d18O_dw_eq,type="l",col=colZ[9])
lines(DATA_past7000$CAVES$record_data$ENTITY95$interp_age, DATA_past7000$CAVES$record_data$ENTITY95$d18O_dw_eq,type="l",col=colZ[10])
lines(DATA_past7000$CAVES$record_data$ENTITY109$interp_age, DATA_past7000$CAVES$record_data$ENTITY109$d18O_dw_eq,type="l",col=colZ[11])
lines(DATA_past7000$CAVES$record_data$ENTITY115$interp_age, DATA_past7000$CAVES$record_data$ENTITY115$d18O_dw_eq,type="l",col=colZ[12])



#################################################

pdf(file = "Plots/ECHAM7k_Trend-Ratio_Scatter_lat.pdf", width = 8, height = 5)
plot(ANALYSIS$MEAN$TREND$CAVES$lat, ANALYSIS$MEAN$TREND$CAVES$value/ANALYSIS$MEAN$TREND$CAVES$value_sim, 
     xlab = "", ylab = "", panel.first = grid(),
     pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2)
abline(h=1)
lines(lowess(ANALYSIS$MEAN$TREND$CAVES$lat[order(ANALYSIS$MEAN$TREND$CAVES$lat)], ANALYSIS$MEAN$TREND$CAVES$value[order(ANALYSIS$MEAN$TREND$CAVES$lat)], f = 2/3, delta = 0.01*180), lwd = 4, col = "#B2182B")
mtext(text = "latitude",side = 1,line = 2)
mtext(text = "Trend Ratio Caves/Sim d18O/ky since 7000BP",side = 2,line = 2)
dev.off()


#################################################
## Caves in China ###############################

pdf(file = "Plots/ECHAM7k_Trend_ChinaCaves.pdf", width = 8, height = 6)
plot(DATA_past7000$CAVES$record_data$ENTITY115$interp_age, DATA_past7000$CAVES$record_data$ENTITY115$d18O_measurement,  
     col = "black", type = "l", xlab = "time (y BP)", ylab = expression(paste(delta^"18",O)), ylim = c(-12,0),
     main = "Caves in China past 7000y")
text(7000, 0, "e_115", col = "black", adj = 1)
lines(DATA_past7000$CAVES$record_data$ENTITY253$interp_age, DATA_past7000$CAVES$record_data$ENTITY253$d18O_measurement,
      col = "red")
text(7000, -0.7, "e_253", col = "red", adj = 1)
lines(DATA_past7000$CAVES$record_data$ENTITY296$interp_age, DATA_past7000$CAVES$record_data$ENTITY296$d18O_measurement,
      col = "green")
text(7000, -1.4, "e_296", col = "green", adj = 1)
lines(DATA_past7000$CAVES$record_data$ENTITY297$interp_age, DATA_past7000$CAVES$record_data$ENTITY297$d18O_measurement,
      col = "blue")
text(7000, -2.1, "e_297", col = "blue", adj = 1)
lines(DATA_past7000$CAVES$record_data$ENTITY329$interp_age, DATA_past7000$CAVES$record_data$ENTITY329$d18O_measurement,
      col = "cyan")
text(7000, -2.8, "e_329", col = "cyan", adj = 1)
lines(DATA_past7000$CAVES$record_data$ENTITY420$interp_age, DATA_past7000$CAVES$record_data$ENTITY420$d18O_measurement,
      col = "darkgreen")
text(7000, -3.4, "e_420", col = "darkgreen", adj = 1)
lines(DATA_past7000$CAVES$record_data$ENTITY496$interp_age, DATA_past7000$CAVES$record_data$ENTITY496$d18O_measurement,
      col = "darkblue")
text(7000, -4.2, "e_496", col = "darkblue", adj = 1)
dev.off()