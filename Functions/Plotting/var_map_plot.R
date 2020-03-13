library(plyr)
library(dplyr)
library(rgdal)
library(latex2exp)

source("Functions/Plotting/STACYmap_5.R")
source("Functions/projection_ptlyr.R")


var_map_plot <- function(Point_Lyr, 
                         projection = as.character('+proj=robin +datum=WGS84'),
                         pt_size,
                         txt_size){
  
  
  mask_china <- logical(length(Point_Lyr$lon))
  
  for(ii in 1:length(Point_Lyr$lon)){
    if(is.na(Point_Lyr$lon[ii])){next}
    if(Point_Lyr$lon[ii] > 100 & Point_Lyr$lon[ii] < 120){
       if(Point_Lyr$lat[ii] < 35 & Point_Lyr$lat[ii] > 22){
        mask_china[ii] = T}
    }
  }
  
  ptlyr_china <- data.frame(
    lon = Point_Lyr$lon[mask_china],
    lat = Point_Lyr$lat[mask_china],
    value = Point_Lyr$value[mask_china]
  )
  
  ptlyr_rest <- data.frame(
    lon = Point_Lyr$lon[!mask_china],
    lat = Point_Lyr$lat[!mask_china],
    value = Point_Lyr$value[!mask_china]
  )
  
  ptlyr_china_p <- projection_ptlyr(ptlyr_china, projection)
  ptlyr_rest_p <- projection_ptlyr(ptlyr_rest, projection)
  
  
  allmax = max(abs(Point_Lyr$value), na.rm = T)
  
  plot <- STACYmap(coastline = TRUE) +
    geom_point(data = ptlyr_china_p, aes(x = long, y = lat, fill = layer), shape = 21, alpha = 0.7, color = "black",
               size = (pt_size), show.legend = c(color =TRUE), position = position_jitter(width = 1000000, height = 500000)) +
    geom_point(data = ptlyr_rest_p, aes(x = long, y = lat, fill = layer), shape = 21, alpha = 0.7, color = "black",
               size = (pt_size), show.legend = c(color =TRUE)) +
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, 'RdBu')), 
                         limits = c(-allmax, allmax),
                         #limits = c(-allmax, log(0.025),log(40), allmax),
                         breaks = c(log10(0.01), log10(0.05), log10(0.1), log10(0.5), 0,log10(5), log10(10), log10(50), log10(100)),
                         labels = c(0.01, "", 0.1, "", 1, "", 10, "", 100),
                         name = TeX("$Var_{Rec}/Var_{Sim}$"), guide = guide_colorbar(barwidth = 10, barheight = 0.3)) +
    theme(legend.direction = "horizontal", 
          panel.border = element_blank(),
          legend.background = element_blank(),
          axis.text = element_blank(),
          text = element_text(size = txt_size),
          legend.title = element_text(size = txt_size))
  
  return(plot)
}

