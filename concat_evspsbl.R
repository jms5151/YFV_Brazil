library(ncdf4)
library(tidyverse)

# source function to extract evspsbl data from a netcdf file
source('function_extract_evspsbl.R')

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
  nc_close(nc)
  
  # add model attributes
  model_name <- gsub('../climate_data/|ssp245_|historical_|_1850_2014|_2015_2100|.nc', '', nc_files[i])
  
  results$model <- model_name
  
  # combine data
  evspsbl_df <- rbind(evspsbl_df, results) 
}

# summarize
evspsbl_df$Year <- NULL

evspsbl_summary <- evspsbl_df %>%
  mutate(Year = as.numeric(format(time, '%Y'))) %>%
  group_by(Year, model) %>%
  summarize(sum_evspsbl = sum(value, na.rm = TRUE)) %>%
  summarize(sum_evspsbl = mean(sum_evspsbl))

# Install and load necessary packages
if (!requireNamespace("fitdistrplus", quietly = TRUE)) install.packages("fitdistrplus")
library(fitdistrplus)

# Fit a distribution (e.g., Gumbel, Normal)
return_period <- function(values, exceedence_prob){
  fit <- fitdist(values, 'norm')
  
  # Extract parameters
  params <- fit$estimate
  
  # Calculate the quantile for a 100-year return period
  target_probability <- 1 / exceedence_prob
  return_value <- qnorm(1 - target_probability, mean = params["mean"], sd = params["sd"])  # Adjust for your distribution
  return(return_value)  
}

pre2017 <- evspsbl_summary$sum_evspsbl[evspsbl_summary$Year < 2017]
rp_100_historical <- return_period(values = pre2017, exceedence_prob = 100)

post2017 <- evspsbl_summary$sum_evspsbl[evspsbl_summary$Year >= 2017]
rp_100_future <- return_period(values = xx, exceedence_prob = 100)

# Plot the histograms
hist(evspsbl_summary$sum_evspsbl[evspsbl_summary$Year < 2017], breaks = 20, col = 'lightblue', main = 'Histogram of evspsbl', xlab = 'evspsbl', xlim = (range(evspsbl_summary$sum_evspsbl)), ylim = c(0, 18))
hist(evspsbl_summary$sum_evspsbl[evspsbl_summary$Year >= 2017], breaks = 40, col = rgb(1,0,0,0.5), add = TRUE)
abline(v = evspsbl_summary$sum_evspsbl[evspsbl_summary$Year == 2017])
abline(v = rp_100_historical, col = 'lightblue')
abline(v = rp_100_future, col=rgb(1,0,0,0.3))
legend('topright', c('Before 2017', 'After 2017'), fill = c('lightblue',  col=rgb(1,0,0,0.3)))


x <- read.csv('../SPEI.csv')
x <- subset(x, !is.na(SPEI.8))
x$SPEI.8 <- x$SPEI.8 * -1
x_100 <- return_period(values = x$SPEI.8, exceedence_prob = 100)
hist(x$SPEI.8, breaks = 20, col = 'lightblue', main = 'Histogram of (SPEI * -1)', xlab = 'SPEI', xlim = (range(x$SPEI.8)), ylim = c(0, 18))
abline(v = x_100, col = 'lightblue')
mean_spei_outbreak <- mean(x$SPEI.8[x$Date >= '2016-01-01' & x$Date < '2017-01-01'])
abline(v = mean_spei_outbreak)
legend('topright', c('once in a century drought', '2016-2017'), col = c('lightblue',  'black'), lty = c(1, 1), bty = 'n')
