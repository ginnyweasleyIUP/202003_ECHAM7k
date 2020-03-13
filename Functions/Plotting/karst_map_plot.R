library(plyr)
library(dplyr)
library(rgdal)
library(latex2exp)

source("Functions/Plotting/STACYmap_5.R")
source("Functions/projection_ptlyr.R")

# coastline_map <- rgdal::readOGR(dsn ="~/07_R_Code/naturalearth_10m_physical/ne_10m_coastline.shp", verbose = FALSE)
# karst_map <- rgdal::readOGR(dsn = "~/07_R_Code/202001_PaperDraft/Karst/karst_wgs.shp", verbose = FALSE)
coastline_map <- rgdal::readOGR(dsn ="/home/ginnyweasley/Dokumente/01_Promotion/07_R_Code/naturalearth_10m_physical/ne_10m_coastline.shp", verbose = FALSE)
karst_map <- rgdal::readOGR(dsn = "/home/ginnyweasley/Dokumente/01_Promotion/07_R_Code/202001_PaperDraft/Karst/karst_wgs.shp", verbose = FALSE)



karst_map_plot <- function(USED_SITES_spec,
                           USED_SITES_var,
                           USED_SITES_mean, 
                           NOT_SITES, 
                      projection = as.character('+proj=robin +datum=WGS84'),
                      pt_size){
  
  karst_map2 <- karst_map %>% spTransform(., CRSobj = CRS(projection)) %>% fortify(.)
  coastline_map2 <- coastline_map %>% spTransform(., CRSobj = CRS(projection)) %>% fortify(.)
  
  USED_SITES_spec_p <- projection_ptlyr(USED_SITES_spec, projection)
  USED_SITES_var_p <- projection_ptlyr(USED_SITES_var, projection)
  USED_SITES_mean_p <- projection_ptlyr(USED_SITES_mean, projection)
  NOT_SITES_p <- projection_ptlyr(NOT_SITES, projection)
  
  plot <- STACYmap(coastline = TRUE, filledbg = TRUE)+
    geom_polygon(data = karst_map2, aes(x=long, y = lat, group = group), fill = '#A57F7C', color = NA) +
    new_scale_color() +
    geom_point(data = NOT_SITES_p, aes(x = long, y = lat, shape = "4", color = '4'),# shape = 20,
               size = (pt_size), alpha = 0.7, show.legend = c(shape =TRUE)) +
    geom_point(data = USED_SITES_spec_p, aes(x = long, y = lat, shape = "1", color = '1'),# shape = 17,
               size = (pt_size), alpha = 0.7, show.legend = c(shape =TRUE)) +
    geom_point(data = USED_SITES_var_p, aes(x = long, y = lat, shape = "2", color = '2'),# shape = 17,
               size = (pt_size), alpha = 0.7, show.legend = c(shape =TRUE)) +
    geom_point(data = USED_SITES_mean_p, aes(x = long, y = lat, shape = "3", color = '3'),# shape = 17,
               size = (pt_size), alpha = 0.7, show.legend = c(shape =TRUE)) +
    scale_color_manual(name = "SISAL v2 sites", labels = c("d18O res > 30n", "d18O res > 20n", "d18O res > 10n", "other"), 
                        values = c('#DD2525', '#F98D11', '#8D0DC4', '#000000')) +
    scale_shape_manual(name = "SISAL v2 sites", labels = c("d18O res > 30n", "d18O res > 20n", "d18O res > 10n", "other"), 
                       values = c(17, 17, 17, 20)) +
    #geom_polygon(data = coastline_map2, aes(x=long, y=lat, group = group), color = 'black',  size = 0.2, fill = NA, alpha = 0.8) +
    theme(legend.position = c(0.02, 0.02),
          legend.justification = c(0, 0),
          legend.box = 'vertical',
          legend.box.background = element_blank(),
          legend.background = element_rect(colour = 'black'),
          panel.border = element_blank())
  
  return(plot)
}

