# Expand drought data from monthly to daily values
# (repeating values for every day in month)

# load libraries
library(tidyverse)
library(lubridate)

# load data
spei <- read.csv('../SPEI.csv')

# format date
spei$Date <- as.Date(spei$Date, '%Y-%m-%d')

# format SPEI indexes: remove positive values and make negative values positive
spei[c('SPEI.1', 'SPEI.3', 'SPEI.8')] <- lapply(spei[c('SPEI.1', 'SPEI.3', 'SPEI.8')], function(x) ifelse(x > 0, 0, x))
spei[c('SPEI.1', 'SPEI.3', 'SPEI.8')] <- lapply(spei[c('SPEI.1', 'SPEI.3', 'SPEI.8')], function(x) abs(x))

# Function to expand the data frame
expand_df <- function(df) {
  df %>%
    # mutate(year_month = format(Date, "%Y-%m-")) %>%
    complete(Date = seq.Date(from = min(Date), to = max(Date), by = "day")) %>%
    fill(SPEI.1, SPEI.3, SPEI.8)
}

# Expand the data frame
spei <- expand_df(spei)
