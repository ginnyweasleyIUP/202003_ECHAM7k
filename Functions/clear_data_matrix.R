clear_data_matrix <- function(data, type){
  
  #temp, prec, isot
  max_value=array(c(60, 1000, 100, 40))
  min_value=array(c(-100, 0, -400, -40))
  
  n = 0
  
  for(lon in 1:dim(data)[1]){
    for(lat in 1:dim(data)[2]){
      for(ii in 1:dim(data)[3]){
        if(is.na(data[lon,lat,ii])){
          next
        } else if(data[lon,lat,ii] > max_value[type]){
          data[lon,lat,ii] <- NA
          n = n +1
        } else if(data[lon,lat,ii] < min_value[type]){
          data[lon,lat,ii] <- NA
          n = n + 1
        }
      }
    }
  }
  print(n)
  return(data)
}

