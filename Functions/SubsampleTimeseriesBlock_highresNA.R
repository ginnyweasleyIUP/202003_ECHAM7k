SubsampleTimeseriesBlock_highresNA <- function (ts, timepoints){
  result <- list()
  dt <- mean(diff(timepoints))
  timepoints.bound <- c(head(timepoints, 1) - dt/2, timepoints[-length(timepoints)] + 
                          diff(timepoints)/2, tail(timepoints + dt/2, 1))
  breaks <- .bincode(c(time(ts)), breaks = timepoints.bound, 
                     TRUE, TRUE)
  temp <- tapply(c(ts), breaks, mean, na.rm = TRUE)
  
  hist <- hist(breaks, length(timepoints.bound)-1, plot = FALSE)
  res_missing <- hist$breaks[hist$counts == 0]
  
  for (ii in res_missing){
    temp_new <- append(temp, NA, after = ii)
    temp <- temp_new
  }
  
  temp_new <- ts(as.vector(temp), start = 1, end = length(timepoints.bound)-1)
  
  return(temp_new)
}