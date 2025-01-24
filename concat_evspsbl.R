library(ncdf4)
library(tidyverse)
library(fitdistrplus)

# source functions to extract evspsbl data from a netcdf file, calculate return period
source('function_extract_evspsbl.R')
source('function_calculate_return_period.R')

# format SPEI data from Brazil
spei_df <- read.csv('../SPEI.csv')
spei_df <- subset(x, !is.na(SPEI.8))
spei_df$SPEI.8 <- spei_df$SPEI.8 * -1

# calculate return period for SPEI data
x_100 <- return_period(values = spei_df$SPEI.8, return_period = 100)

# format climate model data

# list netcdf file
nc_files <- list.files('../climate_data/', full.names = TRUE, pattern = '.nc$')

# create empty data frame to store results
evspsbl_df <- data.frame()

# loop through all netcdf files
for(i in 1:length(nc_files)){
  # Open the NetCDF file
  nc_temp <- nc_open(nc_files[i])
  
  # Inspect the file structure
  # print(nc)
  
  # get data
  results <- extract_evspsbl(nc = nc_temp)
  
  # Close the NetCDF file
  nc_close(nc_temp)
  
  # add model attributes
  model_name <- gsub('../climate_data/|ssp245_|historical_|_1850_2014|_2015_2100|.nc', '', nc_files[i])
  
  results$model <- model_name
  
  # combine data
  evspsbl_df <- rbind(evspsbl_df, results) 
}

# summarize climate data
evspsbl_summary <- evspsbl_df %>%
  group_by(time) %>%
  summarize(sum_evspsbl = mean(value, na.rm = TRUE))

# calculate return period for current and future periods
historical <- evspsbl_summary$sum_evspsbl[evspsbl_summary$time > '1970-10-01' & evspsbl_summary$time <= '2000-01-01']
rp_100_historical <- return_period(values = historical, return_period = 100)

future <- evspsbl_summary$sum_evspsbl[evspsbl_summary$time >= '2070-12-01']
rp_100_future <- return_period(values = future, return_period = 100)

# Plot the histograms
par(mfrow = c(2, 1))

# plot SPEI data
pdf('../figures/SPEI_updated.pdf', width = 8, height = 6)
hist((spei_df$SPEI.8), breaks = 50, col = '#ab0e4a', main = 'Standardized Precipitation-Evapotranspiration Index (SPEI)\nin Minas Gerais, Brazil', xlab = 'SPEI', xlim = c(-3, 3), ylim = c(0, 10))
abline(v = x_100, col = '#ab0e4a', lwd = 2)
peak_drought_b4_outbreak <- (spei_df$SPEI.8[x$Date == '2016-09-01'])
abline(v = peak_drought_b4_outbreak, lwd = 2)
legend('topright', c('Once in a century drought', 'Peak drought prior to outbreak'), col = c('#ab0e4a',  'black'), lty = c(1, 1), lwd = c(2, 2), bty = 'n')
dev.off()

# plot modeled climate data
hist(evspsbl_summary$sum_evspsbl[evspsbl_summary$time > '1970-10-01' & evspsbl_summary$time <= '2000-01-01'], breaks = 40, col = 'orange4', main = 'Current and future freshwater availability', xlab = '(Precipitation - Evapotranspiration)', xlim = range(evspsbl_summary$sum_evspsbl), ylim = c(0, 65))
hist(evspsbl_summary$sum_evspsbl[evspsbl_summary$time >= '2070-12-01'], breaks = 40, col = '#5a940a', add = T)
abline(v = rp_100_historical, col = 'orange4', lwd = 2)
abline(v = rp_100_future, col = '#5a940a', lwd = 2)
legend('topright', c('Historical (1970-2000)', 'Future (2070-2100)'), fill = c('orange4',  col = '#5a940a'), bty = 'n')



