# seasonal forcing function
# seasonal_forcing <- function(times, high, low){
#   yrs <- round(length(times)/365)
#   frequency <- 2*yrs*pi/length(times)
#   days <- seq(1, length(times), by = 1)
#   ts <- (high + low)/2 + (high - low)/2 * cos(days * frequency)
#   return(ts)
# }
# 
seasonal_forcing2 <- function(times, x){
  yrs <- round(length(times)/365)
  frequency <- 2*yrs*pi/length(times)
  days <- seq(1, length(times), by = 1)
  ts <- (x + 0)/2 + (x - 0)/2 * sin(days * frequency - 60)
  return(ts)
}


seasonal_forcing <- function(times, high, low, phase = 0) {
  amplitude <- (high - low) / 2
  mean_val  <- (high + low) / 2
  frequency <- 2 * pi / 365
  # Repeat seasonal pattern every 365 days
  ts <- mean_val + amplitude * cos(frequency * (times %% 365) + phase)
  return(ts)
}
