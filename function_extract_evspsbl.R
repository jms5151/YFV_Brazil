# function to extract evspsbl data from a netcdf file
extract_climvar <- function(nc, var_name){
  # Extract latitude, longitude, and time variables
  lat <- ncvar_get(nc, 'lat')
  lon <- ncvar_get(nc, 'lon')
  time <- ncvar_get(nc, 'time') # Units may need converting
  
  # Get time units and origin
  time_units <- ncatt_get(nc, 'time', 'units')$value
  time_origin <- sub('.*since ', '', time_units)
  time_origin <- as.Date(time_origin, '%Y-%m-%d')
  
  # Convert time to standard dates
  time_converted <- time_origin + as.difftime(time, units = 'days')
  
  # Find the nearest indices for your specific lat/lon
  target_lat <- 19.9191  
  target_lon <- 43.9387 
  lat_idx <- which.min(abs(lat - target_lat))
  lon_idx <- which.min(abs(lon - target_lon))
  
  # Extract the variable for all times at the specific lat/lon
  data_subset <- ncvar_get(nc, var_name, 
                           start = c(lon_idx, lat_idx, 1), 
                           count = c(1, 1, -1))
  # Combine time, data, and model name into a data.frame
  result <- data.frame(
    time = time_converted,
    value = data_subset
  )
  
  return(result)
  
}
