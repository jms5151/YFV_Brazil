# Load required libraries
library(ncdf4)
library(tidyverse)
library(fitdistrplus)
library(ggpubr)

# Source custom functions
source('function_extract_evspsbl.R')
source('function_calculate_return_period.R')

# Load and clean SPEI data
spei_df <- read.csv('../SPEI.csv') %>%
  filter(!is.na(SPEI.8))

# Calculate return period for SPEI data
x_150 <- return_period(values = spei_df$SPEI.8, return_period = 150)

# Function to process NetCDF files
process_nc_file <- function(file) {
  variable_name <- ifelse(grepl('ET', file), 'evspsbl', 'pr')
  nc_data <- nc_open(file)
  results <- extract_climvar(nc = nc_data, var_name = variable_name)
  nc_close(nc_data)
  
  results %>%
    mutate(
      model = gsub('../climate_data/|ET_|precip_|ssp245_|historical_|_1850_2014|_2015_2100|.nc', '', file),
      variable = variable_name
    )
}

# Process climate model data
nc_files <- list.files('../climate_data/', full.names = TRUE, pattern = '.nc$')
df <- map_dfr(nc_files, process_nc_file)

# Summarize climate model data
df_summary <- df %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate(diff = evspsbl - pr) %>%
  group_by(time) %>%
  summarise(diff50 = mean(diff) * 86400, diffsd = sd(diff) * 86400, .groups = 'drop') %>%
  mutate(
    time_period = case_when(
      time >= '2070-12-01' ~ 'Future (2070-2100)',
      time > '1970-10-01' & time <= '2000-01-01' ~ 'Historical (1970-2000)',
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(time_period))

# Compute 1 - CDF
cdf_data <- df_summary %>%
  arrange(diff50) %>%
  group_by(time_period) %>%
  mutate(cdf = ecdf(diff50)(diff50), survival = 1 - cdf) %>%
  ungroup()

# Define custom theme
custom_theme <- theme_bw() +
  theme(
    legend.position = c(0.18, 0.6),
    legend.background = element_rect(fill = 'white', color = 'black'),
    legend.title = element_text(size = 10, face = 'bold'),
    legend.text = element_text(size = 9)
  )

# Define colors
alphaPink <- adjustcolor('#ab0e4a', alpha.f = 0.6)
alphaGreen <- adjustcolor('#5a940a', alpha.f = 0.8)
alphaOrange <- adjustcolor('orange4', alpha.f = 0.8)

# CDF Plot
cum_prob_threshold <- 0.01
cdf_plot <- ggplot(cdf_data, aes(x = diff50, y = pmax(survival, 1e-3), color = time_period)) +
  geom_step(size = 1.5) +
  geom_hline(yintercept = 0.01, linetype = 'dashed', color = 'black', lwd = 1) +
  scale_y_log10(breaks = c(1, 0.1, 0.01, 0.001), labels = scales::comma_format(accuracy = 0.001), limits = c(1e-3, 1)) +
  scale_color_manual(values = c('Historical (1970-2000)' = alphaOrange, 'Future (2070-2100)' = alphaGreen)) +
  labs(title = 'Global Change Increases Dry Conditions', x = 'Evapotranspiration - Precipitation (mm/day)', y = '1 - Cumulative Probability (log scale)', color = 'Time Period') +
  custom_theme +
  theme(legend.background = element_rect(color = NA)) +
  xlim(-0.3, 0.27)

# SPEI Histogram
peak_drought_b4_outbreak <- spei_df$SPEI.8[spei_df$Date == '2016-09-01']
spei_plot <- ggplot(spei_df, aes(x = SPEI.8 * -1)) +
  geom_histogram(bins = 50, fill = alphaPink, color = 'black') +
  geom_vline(aes(xintercept = x_150 * -1), color = '#ab0e4a', linewidth = 1.5) +
  geom_vline(aes(xintercept = peak_drought_b4_outbreak * -1), color = 'black', linewidth = 1.5) +
  labs(title = 'Extreme Drought in Minas Gerais, Brazil', x = 'Inverse Standardized Precipitation - Evapotranspiration Index (SPEI)', y = 'Frequency') +
  scale_x_continuous(limits = c(-3, 3)) +
  scale_y_continuous(limits = c(0, 13)) +
  # Custom Legend
  guides(fill = guide_legend(override.aes = list(fill = alphaPink, shape = 22, size = 2)),
         color = guide_legend(override.aes = list(linetype = c(1, 1), size = c(1.5, 1.5)))) +
  # Manual Legend using Dummy Data
  annotate('segment', x = -3, xend = -2.5, y = 12, yend = 12, color = '#ab0e4a', size = 1.5) +
  annotate('text', x = -2.4, y = 12, label = 'Once in a Century Drought', hjust = 0) +
  annotate('segment', x = -3, xend = -2.5, y = 11, yend = 11, color = 'black', size = 1.5) +
  annotate('text', x = -2.4, y = 11, label = 'Peak Drought Prior to Outbreak', hjust = 0) +
  annotate('rect', xmin = -3, xmax = -2.5, ymin = 9.5, ymax = 10, fill = alphaPink, color = 'black') +
  annotate('text', x = -2.4, y = 9.75, label = 'SPEI (2007-2020)', hjust = 0) +
  custom_theme

# Combine and save plots
drought_plots <- ggarrange(spei_plot, cdf_plot, ncol = 1, nrow = 2)
ggsave('../Figures/drought_plots.pdf', plot = drought_plots, width = 6, height = 8)

# Calculate differences in extremes
spei_ratio <- peak_drought_b4_outbreak / x_150

# Climate model extremes
closest_diff50_by_period <- cdf_data %>%
  group_by(time_period) %>%
  summarise(closest_diff50 = diff50[which.min(abs(survival - cum_prob_threshold))], .groups = 'drop')

climate_ratio <- closest_diff50_by_period$closest_diff50[1] / closest_diff50_by_period$closest_diff50[2]
