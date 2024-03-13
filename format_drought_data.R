# Expand drought data from monthly to daily values
# (repeating values for every day in month)

# load libraries
library(tidyverse)
library(lubridate)

# load data
spei <- read.csv('../SPEI.csv')

# format date
spei$Date <- as.Date(spei$Date, '%Y-%m-%d')

# Function to expand the data frame
expand_df <- function(df) {
  df %>%
    # mutate(year_month = format(Date, "%Y-%m-")) %>%
    complete(Date = seq.Date(from = min(Date), to = max(Date), by = "day")) %>%
    fill(Drought)
}

# Expand the data frame
spei <- expand_df(spei)
# plot(spei$Date, spei$Drought, type = 'l', ylab = 'Inverse SPEI', xlab = 'Date')
