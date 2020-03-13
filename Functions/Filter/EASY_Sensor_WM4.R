#' easy filter
#'
#' easier filter than the convolution filter from Silvia Dee. The function is the filter_fuction.R
#'
#' @param dt time step length in years 
#' @param d18O d18O input that comes back filtered
#' @param tau0 tau0 is the mean transit time through the karst in years
#' @param temp temperatur that fits to the data --> later important for butterworth lowpass filter
#'
#' @return list with filtered-d18O and filter function h
#' @export
#'
#' @examples
easy_sensor_wm4 <- function(dt, d18O, tau0){
  
  #source('~/201903_FirstFilter/functions_used.R')
  tau0 <- tau0/dt
  
  h_s=1E-4						# ensure kernel approaches 0
  tau_s=-tau0*log(tau0*h_s)
  tau=seq(from = 1.0, to = ceiling(tau_s), by =1.0)
  
  timeseries_len <- length(d18O)
  if (length(tau)>= (0.5*timeseries_len)){
    writeLines("Warning: mean residence time (tau0) is too close to half of the length of the timeseries.")
  }
  
  h <- (1/tau0)*exp(-tau/tau0) #kernel length as long as time serie
  hint <- DescTools::AUC(seq(0:(length(h)-1)), h)   # obtain normalization constant using Simpson's rule
  h <- h/hint# normalize transit time distribution
  
  d18OK <- filter_function3(d18O,h, tau0,dt)
  
  return (d18OK)
}
