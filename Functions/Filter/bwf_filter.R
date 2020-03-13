#' butterworth-filter
#'
#' Silvia Dees butter_lowpass filter translated to R
#'
#' @param cutoff frequency where to cut off (in which units???) --> for Dee this is usually 1/10
#' @param order order for butter filter. Order usuall = 3 (WHY???)
#' @param x time series to be filtered
#'
#' @return filtered x by butter filter
#' @export
#'
#' @examples
bwf_filter <- function(x, cutoff, order) {
  if(missing(order)){
    order <-3
  }
  
  butter_list <- signal::butter(order, cutoff)
  m <- mean(x)
  y <-signal::filtfilt (butter_list,x-m) +m
  return(y)
}