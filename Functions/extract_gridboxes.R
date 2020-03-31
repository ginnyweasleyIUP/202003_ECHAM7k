##############################################################
## extract Data for Cave Site from surrounding grid boxes ####
##############################################################

extract_gridboxes_echam <- function(lon_cave, lat_cave){
  result_list <- list()
  
  if(lon_cave<0){lon_cave = 360+lon_cave}
  
  edges <- list(
    E1_lon = lon_cave+3.75/2, E1_lat = lat_cave+3.75/2,
    E2_lon = lon_cave-3.75/2, E2_lat = lat_cave+3.75/2,
    E3_lon = lon_cave+3.75/2, E3_lat = lat_cave-3.75/2,
    E4_lon = lon_cave-3.75/2, E4_lat = lat_cave-3.75/2
  )
  
  if(edges$E1_lon<0){edges$E1_lon = 360 + edges$E1_lon}
  if(edges$E2_lon<0){edges$E2_lon = 360 + edges$E2_lon}
  if(edges$E3_lon<0){edges$E3_lon = 360 + edges$E3_lon}
  if(edges$E4_lon<0){edges$E4_lon = 360 + edges$E4_lon}
  
  #we desire that lon_real>long_grid and lat_real<lat_grid
  edges$E1_lon_pos <- which.min(abs(DATA_past7000$SIM_mean$lon-edges$E1_lon))
  if(DATA_past7000$SIM_mean$lon[edges$E1_lon_pos]>edges$E1_lon){edges$E1_lon_pos = edges$E1_lon_pos -1}
  edges$E1_lat_pos <- which.min(abs(DATA_past7000$SIM_mean$lat-edges$E1_lat))
  if(DATA_past7000$SIM_mean$lat[edges$E1_lat_pos]>edges$E1_lat){edges$E1_lat_pos = edges$E1_lat_pos +1}
  
  edges$E2_lon_pos <- which.min(abs(DATA_past7000$SIM_mean$lon-edges$E2_lon))
  if(DATA_past7000$SIM_mean$lon[edges$E2_lon_pos]>edges$E2_lon){edges$E2_lon_pos = edges$E2_lon_pos -1}
  edges$E2_lat_pos <- which.min(abs(DATA_past7000$SIM_mean$lat-edges$E2_lat))
  if(DATA_past7000$SIM_mean$lat[edges$E2_lat_pos]>edges$E2_lat){edges$E2_lat_pos = edges$E2_lat_pos +1}
  
  edges$E3_lon_pos <- which.min(abs(DATA_past7000$SIM_mean$lon-edges$E3_lon))
  if(DATA_past7000$SIM_mean$lon[edges$E3_lon_pos]>edges$E3_lon){edges$E3_lon_pos = edges$E3_lon_pos -1}
  edges$E3_lat_pos <- which.min(abs(DATA_past7000$SIM_mean$lat-edges$E3_lat))
  if(DATA_past7000$SIM_mean$lat[edges$E3_lat_pos]>edges$E3_lat){edges$E3_lat_pos = edges$E3_lat_pos +1}
  
  edges$E4_lon_pos <- which.min(abs(DATA_past7000$SIM_mean$lon-edges$E4_lon))
  if(DATA_past7000$SIM_mean$lon[edges$E4_lon_pos]>edges$E4_lon){edges$E4_lon_pos = edges$E4_lon_pos -1}
  edges$E4_lat_pos <- which.min(abs(DATA_past7000$SIM_mean$lat-edges$E4_lat))
  if(DATA_past7000$SIM_mean$lat[edges$E4_lat_pos]>edges$E4_lat){edges$E4_lat_pos = edges$E4_lat_pos +1}
  
  ratio <- list()
  
  ratio$E1 <- (edges$E1_lon - DATA_past7000$SIM_mean$lon[edges$E1_lon_pos])*(edges$E1_lat - DATA_past7000$SIM_mean$lat[edges$E1_lat_pos])/(3.75*3.75)
  ratio$E2 <- (DATA_past7000$SIM_mean$lon[edges$E2_lon_pos]+3.75-edges$E2_lon)*(edges$E2_lat - DATA_past7000$SIM_mean$lat[edges$E2_lat_pos])/(3.75*3.75)
  ratio$E3 <- (edges$E3_lon - DATA_past7000$SIM_mean$lon[edges$E3_lon_pos])*(DATA_past7000$SIM_mean$lat[edges$E3_lat_pos-1]-edges$E3_lat)/(3.75*3.75)
  ratio$E4 <- (DATA_past7000$SIM_mean$lon[edges$E4_lon_pos]+3.75-edges$E4_lon)*(DATA_past7000$SIM_mean$lat[edges$E4_lat_pos-1]-edges$E4_lat)/(3.75*3.75)
  #ratio$E1+ratio$E2+ratio$E3+ratio$E4
  
  result_list <- c(edges, ratio)
  
  return(result_list)
}