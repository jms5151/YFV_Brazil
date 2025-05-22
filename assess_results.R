# Create figures

# Constants
LEAD_TIME <- 30
QUANTILES <- c(0.25, 0.5, 0.75)

# --- Helper Functions ---

extract_summary_quantiles <- function(df_list) {
  third_elements <- lapply(df_list, function(x) x[[3]])
  combined_df <- do.call(rbind, third_elements)
  combined_df %>%
    drop_na(Ih, Ip) %>%
    group_by(time) %>%
    summarise(across(c(Ih, Ip), list(
      median = ~ quantile(.x, QUANTILES[2], na.rm = TRUE),
      p25 = ~ quantile(.x, QUANTILES[1], na.rm = TRUE),
      p75 = ~ quantile(.x, QUANTILES[3], na.rm = TRUE)
    )))
}

align_simulation_dates <- function(df_summary, model_name, start_date, epidemic_dates = NULL) {
  df_summary$model <- model_name
  if (grepl('reduce|limit|shift|combined', model_name)) {
    df_summary$Date <- seq.Date(from = start_date, by = 1, length.out = nrow(df_summary))
  } else {
    df_summary$Date <- epidemic_dates[1:nrow(df_summary)]
  }
  return(df_summary)
}

lead_human_cases <- function(df_summary, lead_time = LEAD_TIME) {
  df_summary[, c('Ih_median', 'Ih_p25', 'Ih_p75')] <-
    lapply(df_summary[, c('Ih_median', 'Ih_p25', 'Ih_p75')], lead, n = lead_time)
  return(df_summary)
}

get_medians <- function(dflist, model_name, start_date, epidemic_dates) {
  extract_summary_quantiles(dflist) %>%
    align_simulation_dates(model_name = model_name, start_date = start_date, epidemic_dates = epidemic_dates) %>%
    lead_human_cases()
}

format_for_plotting <- function(df, human_label = 'Infected people', monkey_label = 'Infected monkeys') {
  df %>%
    pivot_longer(cols = c(Ih_median, Ip_median), names_to = "variable", values_to = "value") %>%
    mutate(
      I_p25 = case_when(
        variable == "Ih_median" ~ Ih_p25,
        variable == "Ip_median" ~ Ip_p25
      ),
      I_p75 = case_when(
        variable == "Ih_median" ~ Ih_p75,
        variable == "Ip_median" ~ Ip_p75
      ),
      variable = ifelse(variable == "Ih_median", human_label, monkey_label)
    ) %>%
    select(Date, model, variable, value, I_p25, I_p75)
}

prepare_comparison_df <- function(model_list) {
  do.call(rbind, model_list) %>% distinct()
}

scale_by_reporting_rate <- function(df, rho_humans, rho_monkeys, human_label = 'Infected people', monkey_label = 'Infected monkeys') {
  df[df$model != 'Observed' & df$variable == human_label, c('value', 'I_p25', 'I_p75')] <-
    lapply(df[df$model != 'Observed' & df$variable == human_label, c('value', 'I_p25', 'I_p75')], function(x) x * rho_humans)
  df[df$model != 'Observed' & df$variable == monkey_label, c('value', 'I_p25', 'I_p75')] <-
    lapply(df[df$model != 'Observed' & df$variable == monkey_label, c('value', 'I_p25', 'I_p75')], function(x) x * rho_monkeys)
  return(df)
}

calc_correlation <- function(df, val_data, rho, var_model, var_obs) {
  merged <- left_join(df, val_data, by = "Date")
  merged[[var_model]] <- merged[[var_model]] * rho
  cor_df <- merged[, c(var_model, var_obs)] %>% na.omit()
  test <- cor.test(cor_df[[var_model]], cor_df[[var_obs]])
  return(data.frame(corr = round(test$estimate, 2), pvalue = round(test$p.value, 3)))
}

expand_observed_for_facet <- function(df_observed, df_model) {
  facet_combinations <- expand.grid(
    model = unique(df_model$model),
    variable = unique(df_model$variable)
  )
  merge(df_observed, facet_combinations, by = 'variable')
}

# --- Main Workflow: Load and process model data ---

library(tidyverse)
library(ggpubr)

resultsNew <- readRDS('../model_results.RData')
source('parameters.R')

IC <- read.csv('init_conditions.csv')
start_date <- as.Date(IC$nonnumeric_value[IC$variable == 'start_date'], '%m/%d/%Y')
end_date <- as.Date(IC$nonnumeric_value[IC$variable == 'end_date'], '%m/%d/%Y')

rho_humans <- 0.45
rho_monkeys <- p$value[p$variable == 'p']

val_data <- read.csv('../validation_data.csv')
val_data$Date <- as.Date(val_data$Month, '%m/%d/%Y')
val_data_long <- val_data %>%
  filter(Date >= start_date) %>%
  pivot_longer(cols = c(MG_human, MG_primate), names_to = "variable", values_to = "value") %>%
  mutate(variable = ifelse(variable == 'MG_human', 'Ih_median', 'Ip_median'),
         model = 'Observed', I_p25 = NA, I_p75 = NA)

# Correlation calculation for each model
corDF <- data.frame(matrix(nrow = 0, ncol = 5))
colnames(corDF) <- c('model', 'correlation_humans', 'pvalue_humans', 'correlation_primates', 'pvalue_primates')

for (i in 1:length(resultsNew)) {
  model_name <- names(yfv_params_list)[i]
  df_summary <- get_medians(dflist = resultsNew[[i]], model_name = model_name, start_date = start_date, epidemic_dates = yfv_epidemic)
  cor_human <- calc_correlation(df_summary, val_data, rho_humans, 'Ih_median', 'MG_human')
  cor_primates <- calc_correlation(df_summary, val_data, rho_monkeys, 'Ip_median', 'MG_primate')
  corDF[i,] <- c(model_name, cor_human$corr, cor_human$pvalue, cor_primates$corr, cor_primates$pvalue)
  
  # Store formatted output for plotting
  formatted_df <- format_for_plotting(df_summary)
  observed_df <- if (grepl('reduce|limit|shift|combined', model_name)) val_data_long else val_data_long %>% filter(Date <= end_date)
  combined_df <- rbind(formatted_df, observed_df[, colnames(formatted_df)])
  assign(model_name, combined_df)
}

corDF[,c('correlation_humans', 'correlation_primates')] <- lapply(corDF[,c('correlation_humans', 'correlation_primates')], as.numeric)
corDF$average_correlation <- rowMeans(corDF[,c('correlation_humans', 'correlation_primates')])
write.csv(corDF, 'model_validation.csv', row.names = FALSE)

# --- Plotting Functions ---

plot_theme <- function() {
  theme_bw() +
    theme(
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      strip.text = element_text(size = 12),
      legend.position = 'bottom'
    )
}

panel_plot <- function(df_model, df_observed, custom_colors, custom_labels, title = '') {
  ggplot(df_model, aes(x = Date, y = value, color = model, fill = model)) +
    geom_ribbon(aes(ymin = I_p25, ymax = I_p75), alpha = 0.3, linetype = 0) +
    geom_line(linewidth = 1.1) +
    geom_line(data = df_observed, aes(x = Date, y = value), linewidth = 0.9, inherit.aes = FALSE) +
    geom_point(data = df_observed, aes(x = Date, y = value), size = 2, inherit.aes = FALSE) +
    facet_grid(model ~ variable) +
    xlab('Date') + ylab('Infected') +
    scale_color_manual(values = custom_colors, labels = custom_labels, name = 'Model', breaks = names(custom_labels)) +
    scale_fill_manual(values = custom_colors, labels = custom_labels, name = 'Model', breaks = names(custom_labels)) +
    plot_theme() +
    guides(color = guide_legend(ncol = 1, title = title, override.aes = list(linetype = rep("solid", length(custom_labels)))),
           fill = guide_legend(ncol = 1, title = title)) +
    ggtitle(title)
}

overlaid_plot <- function(df_model, df_observed, custom_colors, custom_labels, title = '') {
  ggplot(df_model, aes(x = Date, y = value, color = model, fill = model)) +
    geom_ribbon(aes(ymin = I_p25, ymax = I_p75), alpha = 0.3, linetype = 0) +
    geom_line(linewidth = 1.1) +
    geom_line(data = df_observed, aes(x = Date, y = value), linewidth = 0.9) +
    geom_point(data = df_observed, aes(x = Date, y = value), size = 2) +
    facet_wrap(~variable, scales = 'free') +
    xlab('Date') + ylab('Infected') +
    scale_color_manual(values = custom_colors, labels = custom_labels, name = 'Model', breaks = names(custom_labels)) +
    scale_fill_manual(values = custom_colors, labels = custom_labels, name = 'Model', breaks = names(custom_labels)) +
    plot_theme() +
    guides(color = guide_legend(ncol = 1, override.aes = list(linetype = rep("solid", length(custom_labels)))),
           fill = guide_legend(ncol = 1)) +
    ggtitle(title)
}

create_comparison_plot <- function(df, custom_colors, custom_labels, title = '', plot_type = 1) {
  df <- scale_by_reporting_rate(df, rho_humans, rho_monkeys)
  df_observed <- df[df$model == 'Observed', ]
  df_model <- df[df$model != 'Observed', ]
  df_observed_expanded <- expand_observed_for_facet(df_observed, df_model)
  
  if (plot_type == 1) {
    panel_plot(df_model, df_observed_expanded, custom_colors, custom_labels, title)
  } else {
    overlaid_plot(df_model, df_observed, custom_colors, custom_labels, title)
  }
}

# --- Figure Creation Blocks ---

save_plot <- function(plot_obj, filename, width, height) {
  ggsave(filename = filename, plot = plot_obj, width = width, height = height)
}

build_and_save_comparison_plot <- function(model_list, colors, labels, filename, title = '', plot_type = 2, 
                                           start_date, end_date, facet_strip_blank = TRUE, legend_position = c(0.2, 0.7), 
                                           y_label = '', subtitle = '', y_lim = c(0, 3000)) {
  df <- prepare_comparison_df(model_list)
  df <- df[df$variable == 'Infected people' & !is.na(df$value), ]
  df$model <- gsub('_high_R0', '', df$model)
  df$model <- factor(df$model, levels = names(labels))
  
  p <- create_comparison_plot(df, custom_colors = colors, custom_labels = labels, title = title, plot_type = plot_type) +
    theme(legend.position = legend_position, legend.background = element_blank()) +
    xlim(start_date, end_date - 180) +
    ylim(y_lim[1], y_lim[2]) +
    labs(title = title, subtitle = subtitle, y = y_label, x = '')
  
  if (facet_strip_blank) {
    p <- p + theme(strip.text = element_blank())
  }
  
  save_plot(p, filename, width = 12, height = 3.5)
  return(p)
}

build_sensitivity_plot <- function(start_date, end_date) {
  mu_plot <- build_and_save_comparison_plot(
    model_list = list(full_model, low_mu_v1, high_mu_v1),
    colors = c('low_mu_v1' = '#1f324a', 'full_model' = '#1561b0', 'high_mu_v1' = '#5fa6dc', 'Observed' = 'black'),
    labels = c('low_mu_v1' = 'Low mortality rate (20%)', 'full_model' = 'Moderate mortality rate (50%)*',
               'high_mu_v1' = 'High mortality rate (80%)', 'Observed' = 'Observed'),
    filename = '../Figures/mu_comparison_plot.pdf',
    title = 'YFV-induced howler monkey mortality rate', start_date = start_date, end_date = end_date,
    y_label = 'Human YFV cases'
  )
  
  pmi_plot <- build_and_save_comparison_plot(
    model_list = list(full_model, high_pMI),
    colors = c('high_pMI' = '#1561b0', 'full_model' = '#5fa6dc', 'Observed' = 'black'),
    labels = c('high_pMI' = '2x literature values', 'full_model' = 'Literature values*', 'Observed' = 'Observed'),
    filename = '../Figures/pMI_comparison_plot.pdf',
    title = 'Probability of mosquito infection with YFV', start_date = start_date, end_date = end_date
  )
  
  move_plot <- build_and_save_comparison_plot(
    model_list = list(full_model, low_movement, high_movement),
    colors = c('low_movement' = '#1f324a', 'full_model' = '#1561b0', 'high_movement' = '#5fa6dc', 'Observed' = 'black'),
    labels = c('low_movement' = 'Low', 'full_model' = 'Moderate*', 'high_movement' = 'High', 'Observed' = 'Observed'),
    filename = '../Figures/movement_comparison_plot.pdf',
    title = 'Immigration rate of non-human primates\nand Hg mosquitos from forest to city', start_date = start_date, end_date = end_date
  )
  
  p_combined <- ggarrange(mu_plot, pmi_plot, move_plot, ncol = 1)
  save_plot(p_combined, '../Figures/Sensitivity_plot.pdf', width = 12, height = 9)
}

build_intervention_plots <- function(start_date, end_date) {
  int_models <- list(full_model, reduce_mosquitoes, limit_monkey_movement, shift_vax, combined_interventions)
  colors <- c('reduce_mosquitoes' = '#ff595e', 'limit_monkey_movement' = '#ffca3a', 'shift_vax' = '#80d819',
              'combined_interventions' = 'purple', 'full_model' = '#1abc9c', 'Observed' = 'black')
  labels <- c('reduce_mosquitoes' = 'Vector control', 'limit_monkey_movement' = 'Conservation (limit NHP movement into city)',
              'shift_vax' = 'Start vaccination earlier', 'combined_interventions' = 'All three interventions combined',
              'full_model' = 'Full model', 'Observed' = 'Observed')
  
  df_all <- prepare_comparison_df(int_models)
  df_all <- df_all[df_all$variable == 'Infected people', ]
  df_all$model <- factor(df_all$model, levels = names(labels))
  
  df_short <- df_all[df_all$Date <= as.Date('2018-07-01'), ]
  
  p_short <- create_comparison_plot(df_short, colors, labels, plot_type = 2) +
    facet_wrap(~variable, scales = 'free') +
    theme(legend.position = 'inside', legend.position.inside = c(0.5, 0.65), legend.background = element_blank(),
          strip.text = element_blank()) +
    labs(title = 'Intervention comparison', subtitle = 'Short-term interventions', y = 'Human YFV cases', x = '') +
    ylim(0, 3000)
  
  p_long <- create_comparison_plot(df_all, colors, labels, plot_type = 2) +
    theme(legend.position = 'none', strip.text = element_blank()) +
    labs(subtitle = 'Long-term interventions', y = '', x = '') +
    ylim(0, 3000)
  
  int_high_R0 <- list(full_model, reduce_mosquitoes_high_R0, limit_monkey_movement_high_R0, shift_vax_high_R0, combined_interventions_high_R0)
  df_high <- prepare_comparison_df(int_high_R0)
  df_high <- df_high[df_high$variable == 'Infected people' & !is.na(df_high$value), ]
  df_high$model <- gsub('_high_R0', '', df_high$model)
  df_high$model <- factor(df_high$model, levels = names(labels))
  
  p_high <- create_comparison_plot(df_high, colors, labels, plot_type = 2) +
    theme(legend.position = 'none', strip.text = element_blank()) +
    labs(subtitle = 'Long-term interventions, higher transmission rate', y = '', x = '') +
    ylim(0, 3000)
  
  p_all <- ggarrange(p_short, p_long, p_high, ncol = 3)
  save_plot(p_all, '../Figures/Intervention_comparison_plot.pdf', width = 12.5, height = 3.5)
}