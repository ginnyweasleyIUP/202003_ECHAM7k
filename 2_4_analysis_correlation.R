#################################################
## ANALYSIS SEASONS #############################
#################################################

library(plyr)
library(dplyr)
library(tidyverse)
library(PaleoSpec)

ANALYSIS$CORR <- list()

#################################################
## CALCULATION ##################################
#################################################

# 1) FIELD (TEMP-ISOT and PREC-ISOT)

ANALYSIS$CORR$FIELD <- list("a" = list(), "b" = list())
for(run in c("a","b")){
  ANALYSIS$CORR$FIELD[[run]]$CORR_TEMP_ISOT = array(dim = c(96,48))
  ANALYSIS$CORR$FIELD[[run]]$CORR_TEMP_ISOT_P = array(dim = c(96,48))
  ANALYSIS$CORR$FIELD[[run]]$CORR_PREC_ISOT = array(dim = c(96,48))
  ANALYSIS$CORR$FIELD[[run]]$CORR_PREC_ISOT_P = array(dim = c(96,48))
  
  for (lon in 1:96){
    for (lat in 1:48){
      #TEMP ISOT
      if(!any(is.na(DATA_past1000[[paste0("SIM_yearly_",run)]]$ISOT[lon,lat,]))){
        COR_TI = cor.test(DATA_past1000[[paste0("SIM_yearly_",run)]]$TEMP[lon,lat,], DATA_past1000[[paste0("SIM_yearly_",run)]]$ISOT[lon,lat,], na.rm = TRUE)
        ANALYSIS$CORR$FIELD[[run]]$CORR_TEMP_ISOT[lon,lat] = COR_TI$estimate[[1]]
        ANALYSIS$CORR$FIELD[[run]]$CORR_TEMP_ISOT_P[lon,lat] = COR_TI$p.value
      }else{
        ANALYSIS$CORR$FIELD[[run]]$CORR_TEMP_ISOT[lon,lat] = NA
        ANALYSIS$CORR$FIELD[[run]]$CORR_TEMP_ISOT_P[lon,lat] = NA
      }
      
      # PREC ISOT
      if(!any(is.na(DATA_past1000[[paste0("SIM_yearly_",run)]]$ISOT[lon,lat,]))){
        COR_PI = cor.test(DATA_past1000[[paste0("SIM_yearly_",run)]]$PREC[lon,lat,], DATA_past1000[[paste0("SIM_yearly_",run)]]$ISOT[lon,lat,], na.rm = TRUE)
        ANALYSIS$CORR$FIELD[[run]]$CORR_PREC_ISOT[lon,lat] = COR_PI$estimate[[1]]
        ANALYSIS$CORR$FIELD[[run]]$CORR_PREC_ISOT_P[lon,lat] = COR_PI$p.value
      }else{
        ANALYSIS$CORR$FIELD[[run]]$CORR_PREC_ISOT[lon,lat] = NA
        ANALYSIS$CORR$FIELD[[run]]$CORR_PREC_ISOT_P[lon,lat] = NA
      }
      
    }
  }
}




remove(lon,lat, COR_TI, COR_PI)


# 2) POINT (TEMP-d18O_dw_eq and PREC-d18O_dw_eq)

length_cave = length(DATA_past1000$CAVES$entity_info$entity_id)

ANALYSIS$CORR$POINTS <- list(
  entity_id = numeric(length_cave)
)

for(run in c("a","b")){
  print(run)
  ANALYSIS$CORR$POINTS[[run]] <- list()
  for(var in c("TEMP","PREC","ISOT","ITPC")){
    print(var)
    ANALYSIS$CORR$POINTS[[run]][[paste0("CORR_",var)]] <- numeric(length_cave)
    ANALYSIS$CORR$POINTS[[run]][[paste0("p_", var)]] <- numeric(length_cave)
    for(ii in 1:length_cave){
      entity <- DATA_past1000$CAVES$entity_info$entity_id[ii]
      site <- DATA_past1000$CAVES$entity_info$site_id[ii]
      ANALYSIS$CORR$POINTS$entity_id[ii] <- entity
      s <- DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]
      # zoo cannot handle objects where order.by has two elements which is why they are sorted out here (no better option found)
      double_time <- s %>% group_by(interp_age) %>% count() %>% filter(n>1)
      s <- s %>% filter(!interp_age %in% double_time$interp_age)
      if(length(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age)>4 & ii != 95 & ii != 53 & ii != 109){
        record <- zoo(x = s$d18O_measurement,
                      order.by = s$interp_age)
        sim <- DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]
        double_time <- sim %>% group_by(interp_age) %>% count() %>% filter(n>1)
        sim <- sim %>% filter(!interp_age %in% double_time$interp_age)
        sim <- zoo(x = sim[[paste0(var,"_", run)]],
                   order.by = sim$interp_age)
        COR <- nexcf_ci(record, sim)

        ANALYSIS$CORR$POINTS[[run]][[paste0("CORR_",var)]][ii] = COR$rxy
        ANALYSIS$CORR$POINTS[[run]][[paste0("p_",var)]][ii] = COR$pval
      }else{
        ANALYSIS$CORR$POINTS[[run]][[paste0("CORR_",var)]][ii] = NA
        ANALYSIS$CORR$POINTS[[run]][[paste0("p_",var)]][ii] = NA
      }
    }
  }
}


## EQUIDISTANT OPTION
# for(run in c("a","b","c")){
#   print(run)
#   ANALYSIS$CORR$POINTS[[run]] <- list()
#   for(var in c("TEMP","PREC","ISOT","ITPC")){
#     print(var)
#     ANALYSIS$CORR$POINTS[[run]][[paste0("CORR_",var)]] <- numeric(length_cave)
#     ANALYSIS$CORR$POINTS[[run]][[paste0("p_", var)]] <- numeric(length_cave)
#     for(ii in 1:length_cave){
#       entity <- DATA_past1000$CAVES$entity_info$entity_id[ii]
#       site <- DATA_past1000$CAVES$entity_info$site_id[ii]
#       ANALYSIS$CORR$POINTS$entity_id[ii] <- entity
#       diff_dt = mean(diff(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age), na.rm = T)
#       if(length(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age)>4 & ii != 95 & ii != 53 & ii != 109 & ii != 158){
#         record <- PaleoSpec::MakeEquidistant(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age,DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]][[paste0("d18O_dw_eq_",run)]],
#                                              time.target = seq(from = head(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age, n = 1),
#                                                                to = tail(DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age, n = 1),
#                                                                by = diff_dt))
#         sim <- PaleoSpec::MakeEquidistant(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]][[paste0("ISOT_",run)]],
#                                           time.target = seq(from = FirstElement(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age),
#                                                             to = LastElement(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age),
#                                                             by = diff_dt))
#         COR <- cor.test(record, sim)
#         
#         ANALYSIS$CORR$POINTS[[run]][[paste0("CORR_",var)]][ii] = COR$estimate[[1]]
#         ANALYSIS$CORR$POINTS[[run]][[paste0("p_",var)]][ii] = COR$p.value
#       }else{
#         ANALYSIS$CORR$POINTS[[run]][[paste0("CORR_",var)]][ii] = NA
#         ANALYSIS$CORR$POINTS[[run]][[paste0("p_",var)]][ii] = NA
#       }
#     }
#   }
# }
