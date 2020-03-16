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
point_lyr$value[74] =NA
point_lyr$value[75] =NA
plot_lyr <- ANALYSIS$MEAN$LAST30$ITPC-ANALYSIS$MEAN$FIRST30$ITPC
plot_lyr[plot_lyr>2.0] = NA

plot <- STACYmap(gridlyr = rbind(plot_lyr[49:96,1:48],plot_lyr[1:48,1:48]),
                 ptlyr = point_lyr,
                 centercolor = 0)
plot
#git comment