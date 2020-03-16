#################################################
## Paper Figure SPECTRUM ########################
#################################################

library(plyr)
library(dplyr)
library(tidyverse)
library(zoo)
library(PaleoSpec)
library(nest)
library(latex2exp)

#################################################

## CALC

# We also use d18O measurement as it is not influenced by simulation

# We need 
# [ ] full yearly spectrum    TEMP, PREC, TEMP (prec weighted), ISOT, ITPC
# [ ] downsamples spectrum    TEMP, PREC, TEMP (prec weighted), ISOT, ITPC
# [ ] record spectrum
# [ ] gridbox weighing

ANALYSIS <- list()
ANALYSIS$SPECTRA <- list()
ANALYSIS$SPECTRA$RECORDS <- list()
ANALYSIS$SPECTRA$SIM_ds <- list()
ANALYSIS$SPECTRA$SIM_full <- list()
ANALYSIS$SPECTRA$MEAN_SPEC <- list()
ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH <- list()

## RECORDS:

length_cave <-length(DATA_past7000$CAVES$entity_info$entity_id)
entities_spec <- list()

# PaleoSpec::SpecMTM needs equally spaced data...

print("Spectrum of Records")

for(ii in 1:length_cave){
  if(ii%%10 == 0){
    print(ii)
  }
  entity = DATA_past7000$CAVES$entity_info$entity_id[ii]
  name = paste0("ENTITY", entity)
  if(length(DATA_past7000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age) > 8 & mask_spec[ii] & entity != 548 & entity != 546 & entity != 547){
    #85 -> eID351 der Quatsch macht allgemein
    entities_spec = c(entities_spec, entity)
    start_ts = ceiling(head(DATA_past7000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age, n = 1))
    length = length(DATA_past7000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age)
    stop_ts = floor(tail(DATA_past7000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age, n = 1))
    if(length > (stop_ts-start_ts)){length = (stop_ts-start_ts)}
    stop_ts = floor((stop_ts-start_ts)/length)*length+start_ts
    
    record <- PaleoSpec::MakeEquidistant(DATA_past7000$CAVES$record_data[[paste0("ENTITY", entity)]]$interp_age,DATA_past7000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_measurement,
                                         time.target = seq(from = start_ts, to = stop_ts, by = floor((stop_ts-start_ts)/length)))
    record = na.omit(record)
    data <- ts(data = record, start = start_ts, end = stop_ts, deltat = floor((stop_ts-start_ts)/length))
    
    ANALYSIS$SPECTRA$RECORDS[[name]] = SpecMTM(data)
    
  }
  
}

ANALYSIS$SPECTRA$entities_spec_rec <- as.numeric(entities_spec)

ANALYSIS$SPECTRA$MEAN_SPEC[["Record"]] <- MeanSpectrum(ANALYSIS$SPECTRA$RECORDS)


# ## WEIGHING
# print("weighing calculations")
# 
# entity_gridbox = array(dim = c(length(ANALYSIS$SPECTRA$entities_spec_rec), 4))
# counter = 1
# 
# for(entity in ANALYSIS$SPECTRA$entities_spec_rec){
#   entity_gridbox[counter, 1] = entity 
#   entity_gridbox[counter, 2] = ANALYSIS$NETWORK$entity_meta$gridbox_id[ANALYSIS$NETWORK$entity_meta$entity_id == entity]
#   entity_gridbox[counter, 3] = ANALYSIS$NETWORK$entity_meta$gridbox_lat[ANALYSIS$NETWORK$entity_meta$entity_id == entity]
#   counter = counter+1
# }
# 
# entity_gridbox <- as.data.frame(entity_gridbox)
# colnames(entity_gridbox) <- c("entity_id", "gridbox_id", "gridbox_lat", "weighing")
# total_gb <- entity_gridbox %>% select(gridbox_id, gridbox_lat) %>% group_by(gridbox_id) %>% count()
# total_lat <- entity_gridbox %>% select(gridbox_lat) %>% group_by(gridbox_lat) %>% count()
# 
# 
# for(entity in entity_gridbox$entity_id){
#   gridbox = entity_gridbox$gridbox_id[entity_gridbox$entity_id == entity]
#   lat = entity_gridbox$gridbox_lat[entity_gridbox$entity_id == entity]
#   lat_weight = cos((90-lat*3.75)*pi/180)/sum(cos(seq(from = -90, to = 90, length.out = 48)*pi/180))/total_lat$n[total_lat$gridbox_lat == lat]
#   entity_gridbox$weighing[entity_gridbox$entity_id == entity] = 1/dim(total_gb)[1]/total_gb$n[total_gb$gridbox_id == gridbox]*lat_weight
# }
# 
# entity_gridbox$weighing <- entity_gridbox$weighing/sum(entity_gridbox$weighing)



# ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[["Record"]] <- MeanSpectrum(ANALYSIS$SPECTRA$RECORDS, weights = entity_gridbox$weighing)


## SIMULATION DOWNSAMPLED

ANALYSIS$SPECTRA$SIM_ds <- list(
  TEMP = list(),
  PREC = list(),
  ISOT = list(), 
  ITPC = list()
)
print("Spectra of downsampled simulation")
for(ii in 1:length(ANALYSIS$SPECTRA$entities_spec_rec)){
  if(ii%%10 == 0){
    print(ii)
  }
  entity = ANALYSIS$SPECTRA$entities_spec_rec[ii]
  name = paste0("ENTITY", entity)

  start_ts = ceiling(head(DATA_past7000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, n = 1))
  length = length(DATA_past7000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age)
  stop_ts = floor(tail(DATA_past7000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, n = 1))
  if(length > (stop_ts-start_ts)){length = (stop_ts-start_ts)}
  stop_ts = floor((stop_ts-start_ts)/length)*length+start_ts
  #ISOT
  record <- PaleoSpec::MakeEquidistant(DATA_past7000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age,
                                       DATA_past7000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ISOT,
                                       time.target = seq(from = start_ts, to = stop_ts, by = floor((stop_ts-start_ts)/length)))
  record = na.omit(record)
  data  = ts(data = record, start = start_ts, end   = stop_ts, deltat = floor((stop_ts-start_ts)/length))
  
  ANALYSIS$SPECTRA$SIM_ds$ISOT[[name]] = SpecMTM(data)
  
  #ITPC
  record <- PaleoSpec::MakeEquidistant(DATA_past7000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age,
                                       DATA_past7000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ITPC,
                                       time.target = seq(from = start_ts, to = stop_ts, by = floor((stop_ts-start_ts)/length)))
  record = na.omit(record)
  data  = ts(data = record, start = start_ts, end   = stop_ts, deltat = floor((stop_ts-start_ts)/length))
  
  ANALYSIS$SPECTRA$SIM_ds$ITPC[[name]] = SpecMTM(data)
  
  #TEMP
  record <- PaleoSpec::MakeEquidistant(DATA_past7000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age,
                                       DATA_past7000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$TEMP,
                                       time.target = seq(from = start_ts, to = stop_ts, by = floor((stop_ts-start_ts)/length)))
  record = na.omit(record)
  data  = ts(data = record, start = start_ts, end   = stop_ts, deltat = floor((stop_ts-start_ts)/length))
  
  ANALYSIS$SPECTRA$SIM_ds$TEMP[[name]] = SpecMTM(data)
  
  #PREC
  record <- PaleoSpec::MakeEquidistant(DATA_past7000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age,
                                       DATA_past7000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$PREC,
                                       time.target = seq(from = start_ts, to = stop_ts, by = floor((stop_ts-start_ts)/length)))
  record = na.omit(record)
  data  = ts(data = record, start = start_ts, end   = stop_ts, deltat = floor((stop_ts-start_ts)/length))
  
  ANALYSIS$SPECTRA$SIM_ds$PREC[[name]] = SpecMTM(data)
  
}

ANALYSIS$SPECTRA$MEAN_SPEC$TEMP$ds <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_ds$TEMP)
ANALYSIS$SPECTRA$MEAN_SPEC$PREC$ds <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_ds$PREC)
ANALYSIS$SPECTRA$MEAN_SPEC$ISOT$ds <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_ds$ISOT)
ANALYSIS$SPECTRA$MEAN_SPEC$ITPC$ds <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_ds$ITPC)

## SIMULATION FULL 
print("Spectra of full simulation")
ANALYSIS$SPECTRA$SIM_full <- list(
  TEMP = list(),
  PREC = list(),
  ISOT = list(), 
  ITPC = list()
)

for(ii in 1:length(ANALYSIS$SPECTRA$entities_spec_rec)){
  if(ii%%10 == 0){
    print(ii)
  }
  entity = ANALYSIS$SPECTRA$entities_spec_rec[ii]
  site = DATA_past7000$CAVES$entity_info$site_id[DATA_past7000$CAVES$entity_info$entity_id == entity]
  name = paste0("ENTITY", entity)
  
  record = na.omit(DATA_past7000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ISOT)
  data  = ts(data = rev(record), 
             start = DATA_past7000$time[2], 
             end   = DATA_past7000$time[1], 
             deltat = 1)
  
  ANALYSIS$SPECTRA$SIM_full$ISOT[[name]] = SpecMTM(data)
  
  record = na.omit(DATA_past7000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$ITPC)
  data  = ts(data = rev(record), 
             start = DATA_past7000$time[2], 
             end   = DATA_past7000$time[1], 
             deltat = 1)
  
  ANALYSIS$SPECTRA$SIM_full$ITPC[[name]] = SpecMTM(data)
  
  record = na.omit(DATA_past7000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$TEMP)
  data  = ts(data = rev(record), 
             start = DATA_past7000$time[2], 
             end   = DATA_past7000$time[1],
             deltat = 1)
  
  ANALYSIS$SPECTRA$SIM_full$TEMP[[name]] = SpecMTM(data)
  
  record = na.omit(DATA_past7000$CAVES$sim_data_yearly[[paste0("CAVE", site)]]$PREC)
  data  = ts(data = rev(record), 
             start = DATA_past7000$time[2], 
             end   = DATA_past7000$time[1], 
             deltat = 1)
  
  ANALYSIS$SPECTRA$SIM_full$PREC[[name]] = SpecMTM(data)
    
}

ANALYSIS$SPECTRA$MEAN_SPEC$TEMP$full <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_full$TEMP)
ANALYSIS$SPECTRA$MEAN_SPEC$PREC$full <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_full$PREC)
ANALYSIS$SPECTRA$MEAN_SPEC$ISOT$full <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_full$ISOT)
ANALYSIS$SPECTRA$MEAN_SPEC$ITPC$full <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_full$ITPC)


# #################################################
# ## Sim full FILTER Spectra ######################
# # filtere full Sim signal auf down sim signal
# # filtere full sim signal auf rec signal
# # filtere down sim signal auf rec signal
# 
# print("Filter process start")
# 
# source('Functions/Filter/EASY_Sensor_WM4.R')
# source('Functions/Filter/filter_function3.R')
# 
# 
# ANALYSIS$SPECTRA$SIM_filter_full_down <- list(TEMP = list(), PREC = list(), ISOT = list(), ITPC = list())
# ANALYSIS$SPECTRA$SIM_filter_full_rec <- list(TEMP = list(), PREC = list(), ISOT = list(), ITPC = list())
# ANALYSIS$SPECTRA$SIM_filter_down_rec <- list(TEMP = list(), PREC = list(), ISOT = list(), ITPC = list())
# 
# filter <- list(
#   ISOT = c(2.0, 20.0, 10.0),
#   ITPC = c(3.0, 12.0, 4.0)
# )
# 
# for(ii in 1:length(ANALYSIS$SPECTRA$entities_spec_rec)){
#   if(ii%%10 == 0){
#     print(ii)
#   }
#   entity = ANALYSIS$SPECTRA$entities_spec_rec[ii]
#   site = DATA_past7000$CAVES$entity_info$site_id[DATA_past7000$CAVES$entity_info$entity_id == entity]
#   name = paste0("ENTITY", entity)
#   
#   #full -> down
#   if(ii%%10 == 0){
#     print("full -> down")
#   }
#   
#   for(run in c("a", "b")){
#     for(var in c("ISOT","ITPC")){
#       Results <- easy_sensor_wm4(1.0, na.omit(DATA_past7000$CAVES$sim_data_yearly[[paste0("CAVE", site)]][[paste0(var,"_",run)]]), filter[[var]][1])
#       data  = ts(data = rev(Results), start = 1950-DATA_past7000$time[2]-29, end   = 1950-DATA_past7000$time[1], deltat = 1)
#       
#       ANALYSIS$SPECTRA$SIM_filter_full_down[[var]][[run]][[name]] = SpecMTM(data)
#     }
#   }
#   
#   #full-> rec
#   if(ii%%10 == 0){
#     print("full -> rec")
#   }
#   for(run in c("a", "b")){
#     for(var in c("ISOT","ITPC")){
#       Results <- easy_sensor_wm4(1.0, na.omit(DATA_past7000$CAVES$sim_data_yearly[[paste0("CAVE", site)]][[paste0(var,"_", run)]]), filter[[var]][2])
#       data  = ts(data = rev(Results), start = 1950-DATA_past7000$time[2]-144, end   = 1950-DATA_past7000$time[1], deltat = 1)
#       
#       ANALYSIS$SPECTRA$SIM_filter_full_rec[[var]][[run]][[name]] = SpecMTM(data)
#     }
#   }
#   
#   
#   #down->rec
#   if(ii%%10 == 0){
#     print("down -> rec")
#   }
#   start_ts = ceiling(head(DATA_past7000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, n = 1))
#   length = length(DATA_past7000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age)
#   stop_ts = floor(tail(DATA_past7000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age, n = 1))
#   diff = floor((stop_ts-start_ts)/length)
#   if(diff<1){diff = 1}
#   
#   for(run in c("a", "b")){
#     for(var in c("ISOT","ITPC")){
#       record <- PaleoSpec::MakeEquidistant(DATA_past7000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age,
#                                            DATA_past7000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]][[paste0(var,"_", run)]],
#                                            time.target = seq(from = start_ts, to = stop_ts, by = diff))
#       Results <- easy_sensor_wm4(diff, na.omit(rev(record)), filter[[var]][3])
#       data = ts(data = Results, start = start_ts-length(Results)+length(record), deltat = diff)
#       ANALYSIS$SPECTRA$SIM_filter_down_rec[[var]][[run]][[name]] = SpecMTM(data)
#     }
#   }
# }
# 
# print("Mean Spectra for filter")
# for(var in c("ISOT", "ITPC")){
#   for(run in c("a", "b")){
#     ANALYSIS$SPECTRA$MEAN_SPEC[[var]][[paste0("full_down_", run)]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_full_down[[var]][[run]])
#     ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[[var]][[paste0("full_down_", run)]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_full_down[[var]][[run]], weights = entity_gridbox$weighing)
#     ANALYSIS$SPECTRA$MEAN_SPEC[[var]][[paste0("full_rec_", run)]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_full_rec[[var]][[run]])
#     ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[[var]][[paste0("full_rec_", run)]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_full_rec[[var]][[run]], weights = entity_gridbox$weighing)
#     ANALYSIS$SPECTRA$MEAN_SPEC[[var]][[paste0("down_rec_", run)]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_down_rec[[var]][[run]])
#     ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[[var]][[paste0("down_rec_", run)]] <- MeanSpectrum(ANALYSIS$SPECTRA$SIM_filter_down_rec[[var]][[run]], weights = entity_gridbox$weighing)
#   }
#   
#   ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_down_",var)]]<- MeanSpectrum(list(ANALYSIS$SPECTRA$MEAN_SPEC[[var]]$full_down_a$spec, 
#                                                                                  ANALYSIS$SPECTRA$MEAN_SPEC[[var]]$full_down_b$spec))
#   ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[[paste0("SIM_full_down_",var)]]<- MeanSpectrum(list(ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[[var]]$full_down_a$spec,
#                                                                                        ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[[var]]$full_down_b$spec))
#   ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_full_rec_",var)]]<- MeanSpectrum(list(ANALYSIS$SPECTRA$MEAN_SPEC[[var]]$full_rec_a$spec, 
#                                                                                  ANALYSIS$SPECTRA$MEAN_SPEC[[var]]$full_rec_b$spec))
#   ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[[paste0("SIM_full_rec_",var)]]<- MeanSpectrum(list(ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[[var]]$full_rec_a$spec,
#                                                                                        ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[[var]]$full_rec_b$spec))
#   ANALYSIS$SPECTRA$MEAN_SPEC[[paste0("SIM_down_rec_",var)]]<- MeanSpectrum(list(ANALYSIS$SPECTRA$MEAN_SPEC[[var]]$down_rec_a$spec, 
#                                                                                 ANALYSIS$SPECTRA$MEAN_SPEC[[var]]$down_rec_b$spec))
#   ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[[paste0("SIM_down_rec_",var)]]<- MeanSpectrum(list(ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[[var]]$down_rec_a$spec,
#                                                                                       ANALYSIS$SPECTRA$MEAN_SPEC_WEIGH[[var]]$down_rec_b$spec))
# }
# 
# remove(dist_matrix, entities_spec, entity_gridbox, total_gb, total_lat, counter, data, diff, entity, filter, gridbox,
#        ii, lat, lat_weight, length, length_cave, name,record, Results, run, site, start_ts, stop_ts, var, easy_sensor_wm4, filter_function3,
#        simpleawmean, simpleawsd)

#git comment