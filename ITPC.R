library(plyr)
library(dplyr)
library(stacy.hadcm.tools)
library(PaleoSpec)
library(nest)
library(tidyverse)

source("Functions/clear_data_matrix.R")

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

ITPC = array(dim = c(96,48,6997))

for (lon in (1:96)){
  for (lat in 1:48){
    print(paste(lon,lat))
    for(year in 1:6997){
      pos_start = 12*(year-1)+1
      pos_stop  = 12*(year-1)+12
      ITPC[lon,lat,year] = sum(PREC[lon,lat, pos_start:pos_stop]*ISOT[lon,lat, pos_start:pos_stop], na.rm = T)/sum(PREC[lon,lat, pos_start:pos_stop], na.rm = T)
    }
  }
}