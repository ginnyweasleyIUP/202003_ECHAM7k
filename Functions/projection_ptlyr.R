library(dplyr)
library(tidyverse)
library(rgdal)

projection_ptlyr <- function(ptlyr, projection){

  if (!is.null(ptlyr)) {
    if (!any(class(ptlyr) %in% c('matrix', 'data.frame'))) {
      stop('class(ptlyr) has to be one of \'matrix\', \'data.frame\'')
    } else {
      prep <- list(matrix = function(x) as_tibble(as.data.frame(x)), 
                   data.frame = function(x) as_tibble(x)
      )
      ptlyr <- ptlyr %>% 
        prep[[class(ptlyr)[length(class(ptlyr))]]](.) %>% 
        rename(long = 1, lat = 2, layer = 3)
      
      ptlyr <- ptlyr %>% 
        rownames_to_column('ptid')
      ptlyr_trf <- project(cbind(ptlyr$long, ptlyr$lat),
                           proj = as.character(projection)) %>% 
        as_tibble() %>% 
        rename(long = V1, lat = V2) %>% 
        bind_cols(select(ptlyr, layer))
      ptlyr <- ptlyr_trf
      rm(ptlyr_trf)
    }
  }
  
  return(ptlyr)
}