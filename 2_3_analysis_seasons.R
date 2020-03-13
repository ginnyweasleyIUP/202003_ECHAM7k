#################################################
## ANALYSIS SEASONS #############################
#################################################

library(plyr)
library(dplyr)
library(tidyverse)

#################################################

ANALYSIS$SEASONS <- list()

source("Functions/SubsampleTimeseriesBlock_highresNA.R")

## We need Correlation of rec with winter, sommer, spring, autumn temp, prec, isot, itpc?

entity_list <- DATA_past1000$CAVES$entity_info$entity_id[mask_var]

for(run in c("a", "b")){
  
  
  for(var in c("TEMP", "PREC", "ISOT")){
    ANALYSIS$SEASONS[[paste0("Plot_Lyr_",var,"_",run)]] <- data.frame(
        lon = numeric(length(entity_list)),
        lat = numeric(length(entity_list)),
        value = numeric(length(entity_list)),
        season = numeric(length(entity_list))
    )
      
    counter = 1
    for(entity in entity_list){
      print(entity)
        
      TEST <- list(
        corr_TEMP = numeric(5),
        corr_TEMP_p = numeric(5),
        corr_PREC = numeric(5),
        corr_PREC_p = numeric(5),
        corr_ISOT = numeric(5),
        corr_ISOT_p = numeric(5)
      )
      site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
      ANALYSIS$SEASONS[[paste0("Plot_Lyr_", var,"_",run)]]$lon[counter] = DATA_past1000$CAVES$site_info$longitude[DATA_past1000$CAVES$site_info$site_id == site]
      ANALYSIS$SEASONS[[paste0("Plot_Lyr_", var,"_",run)]]$lat[counter] = DATA_past1000$CAVES$site_info$latitude[DATA_past1000$CAVES$site_info$site_id == site]
      
      TS <- list()
      s <- DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]
      # zoo cannot handle objects where order.by has two elements which is why they are sorted out here (no better option found)
      double_time <- s %>% group_by(interp_age) %>% count() %>% filter(n>1)
      s <- s %>% filter(!interp_age %in% double_time$interp_age) %>% filter(!is.na(d18O_measurement))
      
      TS[["Record"]] <- zoo(x = s[[paste0("d18O_dw_eq_",run)]], order.by = s$interp_age)
      season_num = 1
      for(season in c("WINTER", "SPRING", "SUMMER", "AUTUMN")){
        TS[["SIM"]] <- zoo( x= SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_seasonal[[run]][[paste0("CAVE", site)]][[season]][[paste0(tolower(var), "_mean")]]), 
                                                                    start = 1950-DATA_past1000$time[2],
                                                                    end = 1950-DATA_past1000$time[1]),
                                                                  s$interp_age), order.by = s$interp_age)
        
        corr <- cor.test(TS$Record, TS$SIM, conflevel = 0.1)
        TEST[[paste0("corr_", var)]][season_num] = corr$estimate[[1]]
        TEST[[paste0("corr_", var, "_p")]][season_num] = corr$p.value
        season_num = season_num + 1
      }

      TS[["SIM"]] <- zoo( x= SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_yearly[[paste0("CAVE", site)]][[paste0(var,"_",run)]]), 
                                                                  start = 1950-DATA_past1000$time[2],
                                                                  end = 1950-DATA_past1000$time[1]),
                                                                s$interp_age), order.by = s$interp_age)
      
      corr <- cor.test(TS$Record, TS$SIM, conflevel = 0.1)
      TEST[[paste0("corr_", var)]][season_num] = corr$estimate[[1]]
      TEST[[paste0("corr_", var, "_p")]][season_num] = corr$p.value
      
      ANALYSIS$SEASONS[[paste0("Plot_Lyr_",var,"_",run)]]$value[counter] = max(abs(TEST[[paste0("corr_", var)]]))
      if(all(is.na(TEST[[paste0("corr_", var)]]))){
        ANALYSIS$SEASONS[[paste0("Plot_Lyr_",var,"_",run)]]$season[counter] = NA
      } else {
        ANALYSIS$SEASONS[[paste0("Plot_Lyr_",var,"_",run)]]$season[counter] = which.max(abs(TEST[[paste0("corr_", var)]]))
      }
      counter = counter + 1
    }#höhle
  }#var
}#run
