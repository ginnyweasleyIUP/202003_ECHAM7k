#################################################
## ANALYSIS SEASONS #############################
#################################################

library(plyr)
library(dplyr)
library(tidyverse)
library(PaleoSpec)

ANALYSIS$REG <- list()

#################################################
## CALCULATION ##################################
#################################################

# 1) FIELD (TEMP-ISOT and PREC-ISOT)

ANALYSIS$REG$FIELD$REG_TEMP_ISOT = array(dim = c(96,48))
ANALYSIS$REG$FIELD$REG_TEMP_ISOT_R = array(dim = c(96,48))
ANALYSIS$REG$FIELD$REG_PREC_ISOT = array(dim = c(96,48))
ANALYSIS$REG$FIELD$REG_PREC_ISOT_R = array(dim = c(96,48))
  
for (lon in 1:96){
  for (lat in 1:48){
    #TEMP ISOT
    if(!any(is.na(DATA_past7000$SIM_yearlyISOT[lon,lat,]))){
      REG_TI = lm(DATA_past7000$SIM_yearly$TEMP[lon,lat,] ~ DATA_past7000$SIM_yearly$ISOT[lon,lat,])
      ANALYSIS$REG$FIELD$REG_TEMP_ISOT[lon,lat] = REG_TI$coefficients[[2]]
      ANALYSIS$REG$FIELD$REG_TEMP_ISOT_R[lon,lat] = summary(REG_TI)$r.squared
    }else{
      ANALYSIS$REG$FIELD$REG_TEMP_ISOT[lon,lat] = NA
      ANALYSIS$REG$FIELD$REG_TEMP_ISOT_R[lon,lat] = NA
    }
    
    # PREC ISOT
    if(!any(is.na(DATA_past7000$SIM_yearly$ISOT[lon,lat,]))){
      REG_PI = lm(DATA_past7000$SIM_yearly$PREC[lon,lat,] ~ DATA_past7000$SIM_yearly$ISOT[lon,lat,])
      ANALYSIS$REG$FIELD$REG_PREC_ISOT[lon,lat] = REG_PI$coefficients[[2]]
      ANALYSIS$REG$FIELD$REG_PREC_ISOT_R[lon,lat] = summary(REG_PI)$r.squared
    }else{
      ANALYSIS$REG$FIELD$REG_PREC_ISOT[lon,lat] = NA
      ANALYSIS$REG$FIELD$REG_PREC_ISOT_R[lon,lat] = NA
    }
    
  }
}

remove(lon,lat, REG_TI, REG_PI)


# 2) POINT (TEMP-d18O_dw_eq and PREC-d18O_dw_eq)

length_cave = length(DATA_past7000$CAVES$entity_info$entity_id)

ANALYSIS$REG$POINTS <- list(
  entity_id = numeric(length_cave)
)

ANALYSIS$REG$POINTS <- list()
for(var in c("TEMP","PREC","ISOT","ITPC")){
  print(var)
  ANALYSIS$REG$POINTS[[paste0("REG_",var)]] <- numeric(length_cave)
  ANALYSIS$REG$POINTS[[paste0("p_", var)]] <- numeric(length_cave)
  for(ii in 1:length_cave){
    entity <- DATA_past7000$CAVES$entity_info$entity_id[ii]
    site <- DATA_past7000$CAVES$entity_info$site_id[ii]
    ANALYSIS$REG$POINTS$entity_id[ii] <- entity
    s <- DATA_past7000$CAVES$record_data[[paste0("ENTITY", entity)]]
    # zoo cannot handle objects where order.by has two elements which is why they are sorted out here (no better option found)
    double_time <- s %>% group_by(interp_age) %>% count() %>% filter(n>1)
    s <- s %>% filter(!interp_age %in% double_time$interp_age)
    if(length(DATA_past7000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age)>4){
      record <- zoo(x = s$d18O_measurement,
                    order.by = s$interp_age)
      sim <- DATA_past7000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]
      double_time <- sim %>% group_by(interp_age) %>% count() %>% filter(n>1)
      sim <- sim %>% filter(!interp_age %in% double_time$interp_age)
      sim <- zoo(x = sim[[var]],
                 order.by = sim$interp_age)
      REG <- lm(as.numeric(sim) ~ as.numeric(record))
        ANALYSIS$REG$POINTS[[paste0("REG_",var)]][ii] = REG$coefficients[[2]]
      ANALYSIS$REG$POINTS[[paste0("p_",var)]][ii] = summary(REG)$r.squared
    }else{
      ANALYSIS$REG$POINTS[[paste0("REG_",var)]][ii] = NA
      ANALYSIS$REG$POINTS[[paste0("p_",var)]][ii] = NA
    }
  }
}


## EQUIDISTANT OPTION
# for(run in c("a","b","c")){
#   print(run)
#   ANALYSIS$REG$POINTS[[run]] <- list()
#   for(var in c("TEMP","PREC","ISOT","ITPC")){
#     print(var)
#     ANALYSIS$REG$POINTS[[run]][[paste0("REG_",var)]] <- numeric(length_cave)
#     ANALYSIS$REG$POINTS[[run]][[paste0("p_", var)]] <- numeric(length_cave)
#     for(ii in 1:length_cave){
#       entity <- DATA_past7000$CAVES$entity_info$entity_id[ii]
#       site <- DATA_past7000$CAVES$entity_info$site_id[ii]
#       ANALYSIS$REG$POINTS$entity_id[ii] <- entity
#       diff_dt = mean(diff(DATA_past7000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age), na.rm = T)
#       if(length(DATA_past7000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age)>4 & ii != 95 & ii != 53 & ii != 109 & ii != 158){
#         record <- PaleoSpec::MakeEquidistant(DATA_past7000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age,DATA_past7000$CAVES$record_data[[paste0("ENTITY", entity)]][[paste0("d18O_dw_eq_",run)]],
#                                              time.target = seq(from = head(DATA_past7000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age, n = 1),
#                                                                to = tail(DATA_past7000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age, n = 1),
#                                                                by = diff_dt))
#         sim <- PaleoSpec::MakeEquidistant(DATA_past7000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, DATA_past7000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]][[paste0("ISOT_",run)]],
#                                           time.target = seq(from = FirstElement(DATA_past7000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age),
#                                                             to = LastElement(DATA_past7000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age),
#                                                             by = diff_dt))
#         COR <- cor.test(record, sim)
#         
#         ANALYSIS$REG$POINTS[[run]][[paste0("REG_",var)]][ii] = COR$estimate[[1]]
#         ANALYSIS$REG$POINTS[[run]][[paste0("p_",var)]][ii] = COR$p.value
#       }else{
#         ANALYSIS$REG$POINTS[[run]][[paste0("REG_",var)]][ii] = NA
#         ANALYSIS$REG$POINTS[[run]][[paste0("p_",var)]][ii] = NA
#       }
#     }
#   }
# }
#git comment