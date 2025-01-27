library(ncdf4)
library(tidyverse)
library(fitdistrplus)

# source functions to extract evspsbl data from a netcdf file, calculate return period
source('function_extract_evspsbl.R')
source('function_calculate_return_period.R')

# format SPEI data from Brazil
spei_df <- read.csv('../SPEI.csv')
spei_df <- subset(spei_df, !is.na(SPEI.8))

# calculate return period for SPEI data
x_100 <- return_period(values = spei_df$SPEI.8, return_period = 100)

# format climate model data

# list netcdf file
nc_files <- list.files('../climate_data/', full.names = TRUE, pattern = '.nc$')

# create empty data frame to store results
df <- data.frame()

# loop through all netcdf files
for(i in 1:length(nc_files)){
  
  if(grepl('ET', nc_files[i])){
    variable_name <- 'evspsbl'
  } else {
    variable_name <- 'pr'
  }
  
  # Open the NetCDF file
  nc_temp <- nc_open(nc_files[i])
  
  # Inspect the file structure
  # print(nc_temp)
  
  # get data
  results <- extract_climvar(nc = nc_temp, var_name = variable_name)
  
  # Close the NetCDF file
  nc_close(nc_temp)
  
  # add model attributes
  results$model <- gsub('../climate_data/|ET_|precip_|ssp245_|historical_|_1850_2014|_2015_2100|.nc', '', nc_files[i])
  results$variable <- variable_name  

  # combine data
  df <- rbind(df, results) 
}

df_summary <- df %>%
  spread(variable, value) %>%
  mutate(diff = pr - evspsbl) %>%
  group_by(time) %>%
  summarize(diff = mean(diff, na.rm = TRUE))

# calculate return period for current and future periods
historical <- df_summary$diff[df_summary$time > '1970-10-01' & df_summary$time <= '2000-01-01']
rp_100_historical <- return_period(values = historical, return_period = 100)

future <- df_summary$diff[df_summary$time >= '2070-12-01']
rp_100_future <- return_period(values = future, return_period = 100)

# Plot the histograms

# plot SPEI data
alphaPink = adjustcolor('#ab0e4a', alpha.f = 0.6)
alphaGreen = adjustcolor('#5a940a', alpha.f = 0.8)
alphaOrange = adjustcolor('orange4', alpha.f = 0.8)

pdf('../figures/drought_histograms.pdf', width = 6, height = 8)
layout(matrix(c(1, 2, 3), nrow = 3, byrow = TRUE), heights = c(4, 4, 1))  # 3 rows, 1 column

hist((spei_df$SPEI.8), breaks = 50, col = alphaPink, main = 'Extreme drought in Minas Gerais, Brazil', xlab = 'Standardized Precipitation-Evapotranspiration Index (SPEI)', xlim = c(-3, 3), ylim = c(0, 13))
abline(v = x_100, col = '#ab0e4a', lwd = 4)
peak_drought_b4_outbreak <- (spei_df$SPEI.8[spei_df$Date == '2016-09-01'])
abline(v = peak_drought_b4_outbreak, lwd = 4)
legend(
  'topright', 
  legend = c('SPEI (2007-2020)', 'Once in a century drought', 'Peak drought prior to outbreak'),
  pch = c(22, NA, NA),  # Square symbols for first two, no symbol for lines
  pt.bg = c(alphaPink, NA, NA),  # Background fill colors for squares
  pt.cex = 2,  # Increase box size
  lty = c(NA, 1, 1),  # First two have no lines, last two have lines
  lwd = c(NA, 4, 4),  # Width of the lines
  col = c('black', '#ab0e4a', 'black'),  # Colors for the lines
  bty = 'n'
)

# plot modeled climate data
hist(df_summary$diff[df_summary$time >= '2070-12-01'], breaks = 40, col = alphaGreen, main = 'Global change increases dry conditions', xlab = 'Net precipitation (precipitation - evapotranspiration)')
hist(df_summary$diff[df_summary$time > '1970-10-01' & df_summary$time <= '2000-01-01'], breaks = 25, col = alphaOrange, add = T)
abline(v = rp_100_historical, col = 'orange4', lwd = 4)
abline(v = rp_100_future, col = '#5a940a', lwd = 4)
legend(
  'topright', 
  legend = c('Historical net precipitation (1970-2000)', 'Future net precipitation (2070-2100)', 
             'Historical 100-yr return period', 'Future 100-yr return period'),
  pch = c(22, 22, NA, NA),  # Square symbols for first two, no symbol for lines
  pt.bg = c(alphaOrange, alphaGreen, NA, NA),  # Background fill colors for squares
  pt.cex = 2,  # Increase box size
  lty = c(NA, NA, 1, 1),  # First two have no lines, last two have lines
  lwd = c(NA, NA, 4, 4),  # Width of the lines
  col = c('black', 'black', 'orange4', '#5a940a'),  # Colors for the lines
  bty = 'n'
)

# Add directional arrow
par(mar = c(0, 4, 0, 2))  # Remove margins for the arrow plot
plot(1, type = 'n', axes = FALSE, xlab = '', ylab = '', xlim = c(-3, 3), ylim = c(0, 1))
arrows(x0 = -2.8, y0 = 0.5, x1 = 2.8, y1 = 0.5, length = 0.1, col = 'black', code = 3)
text(-2.8, 0.3, 'Drier', col = "black", pos = 4)
text(2.8, 0.3, 'Wetter', col = "black", pos = 2)

dev.off()
