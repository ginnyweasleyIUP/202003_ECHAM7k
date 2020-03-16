#################################################
## CREATE DATASET ###############################
#################################################

## 0.1) Set Data-Structure? #####################
## 1.0) read in simulation data from xnapa ######
## 2.0) SISAL TIME SERIES #######################
## 3.0) Extract data from Caves #################
## 4.0) SIMULATION YEARLY #######################
## 5.0) Seasonal Data ###########################
## 6.0) CLIMATE INFO ############################
## 7.0) DOWNSAMPELING ###########################
## 8.0) SIMULATION SEASONS ######################
## 9.0)SIMULATION MEAN ##########################
## 10 )ENTITY in LAT band #######################
## 11 )ARAGONITE, CALCITE, SMOW #################

#################################################
#################################################

library(plyr)
library(dplyr)
library(stacy.hadcm.tools)
library(PaleoSpec)
library(nest)
library(tidyverse)



#################################################
##0.1) Set Data-Structure?#######################
#################################################

DATA_past7000 <- list()
DATA_past7000$SIM_yearly <- list(
  TEMP = list(),
  PREC = list(), 
  ISOT = list(),
  SLPR = list()
)


DATA_past7000$CAVES <- list()
DATA_past7000$CAVES$site_info <- read.csv("/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/02_SISAL/SISAL_v2/site_countries.csv")
DATA_past7000$CAVES$entity_info <- list()
DATA_past7000$CAVES$record_data <- list()
DATA_past7000$CAVES$sim_data_yearly <- list()
DATA_past7000$CAVES$sim_data_downsampled <- list()
DATA_past7000$time <- c(6997, 0)

#################################################
##2) SISAL TIME SERIES ##########################
#################################################

# needs to be imported first, as then only the relevant cave sites will be extracted and calculated further

source("Functions/fun_with_SISALv2_Janica_ECHAM7k.R")

data <- load_sisal_data_janica_echam(year_start = 6997, year_stop = 0)
DATA_past7000$CAVES$entity_info <- data[[1]]
DATA_past7000$CAVES$entity_dating <- data[[3]]
DATA_past7000$CAVES$site_to_entity <- data[[4]]

#Schmeißt alle Höhlen raus, die nicht gebraucht werden
DATA_past7000$CAVES$site_info <- DATA_past7000$CAVES$site_info %>% filter(site_id %in% DATA_past7000$CAVES$entity_info$site_id)
for (ii in DATA_past7000$CAVES$entity_info$entity_id){
  name = paste0("ENTITY", ii)
  #if(ii%%10 == 0){
  print(name)
  #}
  site <- as.numeric(DATA_past7000$CAVES$entity_info %>% filter(entity_id == ii) %>% distinct(site_id))
  DATA_past7000$CAVES$record_data[[name]] <- data[[2]] %>% filter(entity_id == ii) %>% distinct(entity_id, interp_age, d18O_measurement) %>%
    mutate(site_id = site)
  #if(ii == 144 | ii == 226 | ii == 20){next}
  #DATA_past7000$CAVES$record_data[[name]]$chron <- as.tibble(data[[5]] %>% filter(entity_id == ii))
}

remove(data, site, ii, name, load_sisal_data_janica_echam)





#################################################
## 1) Read in DATA ##############################
#################################################

## 1.1) TEMP ####################################
source("Functions/clear_data_matrix.R")

DATA_past7000$SIM_yearly <- list(TEMP = array(dim = c(96,48,abs(diff(DATA_past7000$time)))),
                                 PREC = array(dim = c(96,48,abs(diff(DATA_past7000$time)))),
                                 ISOT = array(dim = c(96,48,abs(diff(DATA_past7000$time)))),
                                 ITPC = array(dim = c(96,48,abs(diff(DATA_past7000$time)))))

for(var in c("TEMP", "PREC", "ISOT")){
#for(var in c("PREC", "ISOT")){
  DATA_past7000_SIM_RAW <- list()
  
  print(var)
  if(var == "TEMP"){
    ncf <- (ncdf4::nc_open("/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/10_ECHAM5/results_Hol-T7k/Hol-T_echam5_wiso_tsurf_0004_to_7000.nc"))
    test <- clear_data_matrix(ncdf4::ncvar_get(ncf),1)
  }
  if(var == "PREC"){
    ncf<-ncdf4::nc_open("/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/10_ECHAM5/results_Hol-T7k/Hol-T_echam5_wiso_aprt_0004_to_7000.nc")
    test <- clear_data_matrix(ncdf4::ncvar_get(ncf),2)/30
  }
  if(var == "ISOT"){
    ncf<-ncdf4::nc_open("/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/10_ECHAM5/results_Hol-T7k/Hol-T_echam5_wiso_wisoaprt_d_0004_to_7000_2.nc")
    test<- clear_data_matrix(ncdf4::ncvar_get(ncf, 'wisoaprt_d'),3)
  }
  
  print("Data in ...")
  
  test_NA <- array(dim = c(96,48,7))
  DATA_past7000_SIM_RAW$TEMP <- abind::abind(test[,,1:66557], test_NA, test[,,66558:83957], along = 3)
  DATA_past7000_SIM_RAW$lon <- ncf$dim$lon$vals
  DATA_past7000_SIM_RAW$lat <- ncf$dim$lat$vals
  ncdf4::nc_close(ncf)
  
  print("Extract grid boxes ...")
  source("Functions/extract_gridboxes.R")
  
  for (ii in 1:(dim(DATA_past7000$CAVES$site_info)[1])){
    lon_cave = DATA_past7000$CAVES$site_info$longitude[ii]
    
    if(lon_cave<0){lon_cave = 360+lon_cave}
    lat_cave = DATA_past7000$CAVES$site_info$latitude[ii]
    site_id = DATA_past7000$CAVES$site_info$site_id[ii]
    
    ratios <- extract_gridboxes_echam(lon_cave, lat_cave)
    name <- paste0("CAVE",site_id)
    DATA_past7000$CAVES$sim_data_raw[[name]][[var]] <- rowSums(cbind(ratios$E1*DATA_past7000_SIM_RAW$TEMP[ratios$E1_lon_pos, ratios$E1_lat_pos,],
                                                                     ratios$E2*DATA_past7000_SIM_RAW$TEMP[ratios$E2_lon_pos, ratios$E2_lat_pos,],
                                                                     ratios$E3*DATA_past7000_SIM_RAW$TEMP[ratios$E3_lon_pos, ratios$E3_lat_pos,],
                                                                     ratios$E4*DATA_past7000_SIM_RAW$TEMP[ratios$E4_lon_pos, ratios$E4_lat_pos,]), na.rm = T)
  }
  
}

rm(DATA_past7000_SIM_RAW)

print("Yearly Mean ...")

#TEMP
ncf <- (ncdf4::nc_open("/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/10_ECHAM5/results_Hol-T7k/Hol-T_echam5_wiso_tsurf_0004_to_7000_yearmean.nc"))
DATA_past7000$SIM_yearly$TEMP <- clear_data_matrix(ncdf4::ncvar_get(ncf),1)
ncdf4::nc_close(ncf)

ncf<-ncdf4::nc_open("/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/10_ECHAM5/results_Hol-T7k/Hol-T_echam5_wiso_aprt_0004_to_7000_yearmean.nc")
DATA_past7000$SIM_yearly$PREC <- clear_data_matrix(ncdf4::ncvar_get(ncf),2)/30
ncdf4::nc_close(ncf)

ncf<-ncdf4::nc_open("/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/10_ECHAM5/results_Hol-T7k/Hol-T_echam5_wiso_wisoaprt_d_0004_to_7000_2.nc")
DATA_past7000$SIM_yearly$ISOT <- clear_data_matrix(ncdf4::ncvar_get(ncf, 'wisoaprt_d'),3)
ncdf4::nc_close(ncf)




save.image()

ncf<-ncdf4::nc_open("/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/10_ECHAM5/results_Hol-T7k/Hol-T_echam5_wiso_aprt_0004_to_7000.nc")
test <- clear_data_matrix(ncdf4::ncvar_get(ncf),2)/30
test_NA <- array(dim = c(96,48,7))
PREC <- abind::abind(test[,,1:66557], test_NA, test[,,66558:83957], along = 3)
ncdf4::nc_close(ncf)
rm(test)
ncf<-ncdf4::nc_open("/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/10_ECHAM5/results_Hol-T7k/Hol-T_echam5_wiso_wisoaprt_d_0004_to_7000_2.nc")
test <- clear_data_matrix(ncdf4::ncvar_get(ncf, 'wisoaprt_d'),3)
ISOT <- abind::abind(test[,,1:66557], test_NA, test[,,66558:83957], along = 3)
ncdf4::nc_close(ncf)
rm(test)

save.image()

for (lon in (1:96)){
  for (lat in 1:48){
    print(paste(lon,lat))
    for(year in 1:diff(DATA_past7000$time)){
      pos_start = 12*(year-1)+1
      pos_stop  = 12*(year-1)+12
      DATA_past7000$SIM_yearly$ITPC[lon,lat,year] = sum(PREC[lon,lat, pos_start:pos_stop]*ISOT[lon,lat, pos_start:pos_stop], na.rm = T)/sum(PREC[lon,lat, pos_start:pos_stop], na.rm = T)
    }
  }
}

load("~/Dokumente/01_Promotion/07_R_Code/202003_ECHAM7k/ITPC.RData")

DATA_past7000$SIM_yearly$ITPC <- ITPC

rm(ncf, ratios, ii, ITPC, lat_cave, lon_cave, name, site_id, test, test_NA, var)
rm(clear_data_matrix, extract_gridboxes_echam)

#################################################
## CAVES YEARLY #################################
#################################################

for (ii in 1:(dim(DATA_past7000$CAVES$site_info)[1])){
  print(ii)
  site_id = DATA_past7000$CAVES$site_info$site_id[ii]
  name <- paste0("CAVE", site_id)
  DATA_past7000$CAVES$sim_data_yearly[[name]]$TEMP <- numeric(abs(diff(DATA_past7000$time)))
  DATA_past7000$CAVES$sim_data_yearly[[name]]$PREC <- numeric(abs(diff(DATA_past7000$time)))
  DATA_past7000$CAVES$sim_data_yearly[[name]]$ISOT <- numeric(abs(diff(DATA_past7000$time)))
  DATA_past7000$CAVES$sim_data_yearly[[name]]$ITPC <- numeric(abs(diff(DATA_past7000$time)))
  for(year in 1:abs(diff(DATA_past7000$time))){
    pos_start = 12*(year-1)+1
    pos_stop  = 12*(year-1)+12
    DATA_past7000$CAVES$sim_data_yearly[[name]]$TEMP[year] <- mean(DATA_past7000$CAVES$sim_data_raw[[name]]$TEMP[pos_start:pos_stop], na.rm = T)
    DATA_past7000$CAVES$sim_data_yearly[[name]]$PREC[year] <- mean(DATA_past7000$CAVES$sim_data_raw[[name]]$PREC[pos_start:pos_stop], na.rm = T)
    DATA_past7000$CAVES$sim_data_yearly[[name]]$ISOT[year] <- mean(DATA_past7000$CAVES$sim_data_raw[[name]]$ISOT[pos_start:pos_stop], na.rm = T)
    DATA_past7000$CAVES$sim_data_yearly[[name]]$ITPC[year] <- sum(DATA_past7000$CAVES$sim_data_raw[[name]]$ISOT[pos_start:pos_stop]*
                                                                    DATA_past7000$CAVES$sim_data_raw[[name]]$PREC[pos_start:pos_stop], na.rm = T)/
      sum(DATA_past7000$CAVES$sim_data_raw[[name]]$PREC[pos_start:pos_stop], na.rm = T)
  }
}

rm(ii, name, pos_start, pos_stop, site_id, year)
  

#################################################
## 7.0) DOWNSAMPELING ###########################
#################################################

source("Functions/SubsampleTimeseriesBlock_highresNA.R")


for(ii in 1:length(DATA_past7000$CAVES$entity_info$entity_id)){
  print(ii)
  nameE = paste0("ENTITY", DATA_past7000$CAVES$entity_info$entity_id[ii])
  nameC = paste0("CAVE", DATA_past7000$CAVES$entity_info$site_id[ii])
  data_temp <- SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past7000$CAVES$sim_data_yearly[[nameC]]$TEMP),
                                                     start = DATA_past7000$time[2],
                                                     end = DATA_past7000$time[1]),
                                                  DATA_past7000$CAVES$record_data[[nameE]]$interp_age)

  data_prec <- SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past7000$CAVES$sim_data_yearly[[nameC]]$PREC),
                                                     start = DATA_past7000$time[2],
                                                     end = DATA_past7000$time[1]),
                                                   DATA_past7000$CAVES$record_data[[nameE]]$interp_age)

  data_isot <- SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past7000$CAVES$sim_data_yearly[[nameC]]$ISOT),
                                                     start = DATA_past7000$time[2],
                                                     end = DATA_past7000$time[1]),
                                                  DATA_past7000$CAVES$record_data[[nameE]]$interp_age)

  data_itpc <- SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past7000$CAVES$sim_data_yearly[[nameC]]$ITPC),
                                                     start = DATA_past7000$time[2],
                                                     end = DATA_past7000$time[1]),
                                                  DATA_past7000$CAVES$record_data[[nameE]]$interp_age)
  
  data <- matrix(c(DATA_past7000$CAVES$record_data[[nameE]]$interp_age,
                   data_temp, data_prec, data_isot, data_itpc), ncol= 5)
  colnames(data) = c("interp_age", "TEMP", "PREC","ISOT", "ITPC")

  DATA_past7000$CAVES$sim_data_downsampled[[nameE]] <- as.tibble(data)

}

remove(nameE, nameC, ii, data, data_isot, data_itpc, data_prec, data_temp)
remove(SubsampleTimeseriesBlock_highresNA)


#################################################
## ARAGONITE, CALCITE, SMOW #####################
#################################################
 
## Aragonite and Calcite d18O have to be converted to drip water equivalents to be comparable
## Also in ORder to be comparable to SMOW which is the standard of the simulation we need to convert the dripwater d18O from VPDB to SMOW
 
for(ii in 1:length(DATA_past7000$CAVES$entity_info$entity_id)){
 entity = DATA_past7000$CAVES$entity_info$entity_id[ii]
 site = DATA_past7000$CAVES$entity_info$site_id[ii]
 print(entity)
 data_rec <- DATA_past7000$CAVES$record_data[[paste0("ENTITY", entity)]]
 data_sim <- DATA_past7000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]
 if(DATA_past7000$CAVES$entity_info$mineralogy[ii] == "calcite"){
   dw_eq <- 1.03092 * (data_rec$d18O_measurement - ((16.1*1000)/(data_sim$TEMP+273.15)-24.6)) + 30.92
 }else if(DATA_past7000$CAVES$entity_info$mineralogy[ii] == "aragonite"){
   dw_eq <- 1.03092 * (data_rec$d18O_measurement - ((18.34*1000)/(data_sim$TEMP+273.15)-31.954)) + 30.92
 }else{
   dw_eq <- numeric(length(data))+NA
 }

 DATA_past7000$CAVES$record_data[[paste0("ENTITY", entity)]] <- data.frame(
   site_id = data_rec$site_id,
   entity_id = data_rec$entity_id,
   interp_age = data_rec$interp_age,
   d18O_measurement = data_rec$d18O_measurement,
   d18O_dw_eq = dw_eq
   )
}

remove(entity, dw_eq, ii, site, data_rec, data_sim)

#################################################
## UMWANDLUNGEN #################################
#################################################

DATA_past7000$CAVES$site_info <- DATA_past7000$CAVES$site_info %>% mutate(elevation = as.numeric(as.character(elevation)))
#remove(DATA_past7000_SIM_RAW)

#################################################
## MASKS ########################################
#################################################

mask_mean = logical(length = length(DATA_past7000$CAVES$entity_info$entity_id))
mask_var  = logical(length = length(DATA_past7000$CAVES$entity_info$entity_id))
mask_spec = logical(length = length(DATA_past7000$CAVES$entity_info$entity_id))

for(entity in 1:length(DATA_past7000$CAVES$entity_info$entity_id)){
  if(DATA_past7000$CAVES$entity_info$n[entity] > 10 & DATA_past7000$CAVES$entity_info$period[entity] > 4000){mask_mean[entity] = T}
  if(DATA_past7000$CAVES$entity_info$n[entity] > 20 & DATA_past7000$CAVES$entity_info$period[entity] > 5000){mask_var[entity] = T}
  if(DATA_past7000$CAVES$entity_info$n[entity] > 30 & DATA_past7000$CAVES$entity_info$period[entity] > 6000){mask_spec[entity] = T}
}

#git comment