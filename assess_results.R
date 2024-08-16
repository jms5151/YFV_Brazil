# load libraries
library(tidyverse)
library(ggpubr)

# load model results
resultsNew <- readRDS('../model_results.RData')
results2 <- readRDS('../model_results_interventions_long.RData')
resultsNew[[14]] <- results2[[1]]

# source parameters code to bring some variables into environment
source('parameters.R')

# reporting rates
rho_humans = 0.12 # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4632853/#:~:text=For%20YF%20virus%20infections%2C%20the,%25%20CI%200.31%E2%80%930.62).
rho_monkeys = p$value[p$variable == 'p']
  
# format validation data
val_data <- read.csv('../validation_data.csv')
val_data$Date <- as.Date(val_data$Month, '%m/%d/%Y')
val_data_long <- val_data %>%
  filter(Date >= start_date & Date <= end_date) %>%
  pivot_longer(
    cols = c(MG_human, MG_primate),
    names_to = "variable",
    values_to = "value"
  ) %>%
  as.data.frame() 

val_data_long$variable <- ifelse(val_data_long$variable == 'MG_human', 'I_h_median', 'I_p_median')  
val_data_long$model <- 'Observed'
val_data_long$I_25 <- NA
val_data_long$I_75 <- NA

# format simulated data
get_medians <- function(dflist){
  df <- dflist
  
  # Extract the third element from each list within grouped_results
  third_elements <- lapply(df, function(x) x[[3]])
  
  # Combine all third elements into a single data frame using do.call and rbind
  combined_df <- do.call(rbind, third_elements)
  
  # summarise
  df2 <- combined_df %>%
    na.omit() %>%
    group_by(time) %>%
    summarise(
      I_h_median = quantile(I_h, 0.5)
      , I_h_25 = quantile(I_h, 0.25)
      , I_h_75 = quantile(I_h, 0.75)
      , I_p_median = quantile(I_p, 0.5)
      , I_p_25 = quantile(I_p, 0.25)
      , I_p_75 = quantile(I_p, 0.75)
    ) %>%
    mutate('model' = names(yfv_params_list)[i])
  
  if(unique(df2$model) == 'reduce_mosquitoes' | unique(df2$model) == 'reduce_NHP_movement' | unique(df2$model) == 'shift_vax' | unique(df2$model) == 'combined_interventions'){
    end_date_new <- start_date + (length(int_times) - 2)
    df2$Date <- seq.Date(from = start_date, to = end_date_new, by = 1)
  } else {
    df2$Date <- yfv_epidemic[1:nrow(df2)]
  }
  
  return(df2)
}

# calculate correlations between simulations and observations
calc_correlation <- function(df, rho, var1, var2){
  # add observed data
  x2 <- df %>% left_join(val_data) %>% as.data.frame()
  x2[, var1] <- x2[, var1] * rho
  x3 <- x2[,c(var1, var2)]
  x3 <- subset(x3, !is.na(x3[,var2]))
  # Calculate the cross-correlation function
  ccf_result <- ccf(x3[,var1], x3[,var2], lag.max = 0, plot = FALSE)
  dist_out <- round(ccf_result$acf[which(ccf_result$lag == 0)], 2)
  return(dist_out)
}

corDF <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(corDF) <- c('model', 'correlation_humans', 'correlation_primates')

for(i in 1:10){
  x <- get_medians(dflist = resultsNew[[i]])
  corHuman <- calc_correlation(df = x, rho = rho_humans, var1 = 'I_h_median', var2 = 'MG_human')
  corPrimates <- calc_correlation(df = x, rho = rho_monkeys, var1 = 'I_p_median', var2 = 'MG_primate')
  corDF[i,] <- c(unique(x$model), corHuman, corPrimates)
}

corDF[,c("correlation_humans", "correlation_primates")] <- lapply(corDF[,c("correlation_humans", "correlation_primates")], function(x) as.numeric(x))
corDF$average_correlation <- rowMeans(corDF[,c("correlation_humans", "correlation_primates")])

write.csv(corDF, 'model_validation.csv', row.names = F)

# format simulated data for plotting
for(i in 1:length(yfv_params_list)){

  df2 <- get_medians(dflist = resultsNew[[i]])
  
  # transform to good format for plotting
  df_long_median <- df2 %>%
    pivot_longer(
      cols = c(I_h_median, I_p_median),
      names_to = "variable",
      values_to = "value"
    )
  
  # Create a key to join the percentiles back to the long format
  df_long_median <- df_long_median %>%
    mutate(
      I_25 = case_when(
        variable == "I_h_median" ~ I_h_25,
        variable == "I_p_median" ~ I_p_25
      ),
      I_75 = case_when(
        variable == "I_h_median" ~ I_h_75,
        variable == "I_p_median" ~ I_p_75
      )
    ) %>%
    select(Date, model, variable, value, I_25, I_75) %>%
    as.data.frame()
  
  # add observed data
  df_long_median <- rbind(df_long_median, val_data_long[,colnames(df_long_median)])
  
  # rename variable labels
  df_long_median$variable <- ifelse(df_long_median$variable == 'I_h_median', 'Infected people', 'Infected monkeys')
  
  # rename dataset
  assign(names(yfv_params_list)[i], df_long_median)
  
}


# Plotting function
create_comparison_plot <- function(df, custom_colors, custom_labels, titleName = '') {
  df[df$model != 'Observed' & df$variable == 'Infected people', c('value', 'I_25', 'I_75')] <- lapply(df[df$model != 'Observed' & df$variable == 'Infected people', c('value', 'I_25', 'I_75')], function(x) x * rho_humans)
  df[df$model != 'Observed' & df$variable == 'Infected monkeys', c('value', 'I_25', 'I_75')] <- lapply(df[df$model != 'Observed' & df$variable == 'Infected monkeys', c('value', 'I_25', 'I_75')], function(x) x * rho_monkeys)
  
  # Separate the observed data
  df_observed <- df[df$model == 'Observed', ]
  df_model <- df[df$model != 'Observed', ]
  
  # plot
  p <- ggplot(df_model, aes(x = Date, y = value, color = model, fill = model)) +
    geom_ribbon(aes(ymin = I_25, ymax = I_75), alpha = 0.3, linetype = 0) +
    geom_line(lwd = 1.1) +
    geom_line(data = df_observed, aes(x = Date, y = value), lwd = 1) +
    geom_point(data = df_observed, aes(x = Date, y = value), size = 3) +
    facet_wrap(~variable, scales = 'free') +
    theme_bw() +
    xlab('Date') +
    ylab('Infected') +
    scale_color_manual(values = custom_colors, labels = custom_labels, name = 'Model', breaks = names(custom_labels)) +
    scale_fill_manual(values = custom_colors, labels = custom_labels, name = 'Model', breaks = names(custom_labels)) +
    theme(legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          strip.text = element_text(size = 12)) +
    theme(legend.position = 'bottom') +
    guides(color = guide_legend(ncol = 1), fill = guide_legend(ncol = 1)) +
    ggtitle(titleName)
  
  return(p)
}

# Model comparison figures --------------
# combine datasets
mod_compare <- do.call(rbind, list(base_model, fixed_bite_rate, fixed_movement, fixed))
mod_compare <- mod_compare[!duplicated(mod_compare), ]

# Define custom colors
mod_comp_colors <- c('base_model' = '#149889'
                   , 'fixed_bite_rate' = '#82118e'
                   , 'fixed_movement' = '#2e3d54'
                   , 'fixed' = '#c5a030'
                   , 'Observed' = 'black')

# Custom labels for the legend
mod_comp_labels <- c('base_model' = 'Seasonal primate movement & seasonal biting rate'
                   , 'fixed_bite_rate' = 'No primate movement, seasonally-driven biting rate'
                   , 'fixed_movement' = 'Seasonally-driven primate movement, fixed biting rate'
                   , 'fixed' = 'No primate movement & fixed biting rate'
                   , 'Observed' = 'Observed')


model_comparison_plot <- create_comparison_plot(df = mod_compare, custom_colors = mod_comp_colors, custom_labels = mod_comp_labels)
ggsave(filename = '../Figures/Model_comparison_plot.pdf', plot = model_comparison_plot, width = 10, height = 5.5)

# Compare YFV primate mortality rates
mu_compare <- do.call(rbind, list(base_model, low_mu_v1, high_mu_v1))
mu_compare <- mu_compare[!duplicated(mu_compare), ]
mu_compare$model <- factor(mu_compare$model, levels = c('high_mu_v1', 'base_model', 'low_mu_v1', 'Observed')) 
  
mu_colors <- c('low_mu_v1' = '#1f324a', 'base_model' = '#1561b0', 'high_mu_v1' = '#5fa6dc', 'Observed' = 'black')
mu_labels <- c('low_mu_v1' = '20%', 'base_model' = '50%', 'high_mu_v1' = '80%', 'Observed' = 'Observed')

mu_comparison_plot <- create_comparison_plot(df = mu_compare, custom_colors = mu_colors, custom_labels = mu_labels, titleName = 'YFV monkey mortality rate') + theme(legend.position = c(0.9, 0.6))

# Compare Hm probability of biting NHP vs humans
p_compare <- do.call(rbind, list(base_model, low_p, mod_p))
p_compare <- p_compare[!duplicated(p_compare), ]
p_compare$model <- factor(p_compare$model, levels = c('base_model', 'mod_p', 'low_p', 'Observed')) 

p_colors <- c('low_p' = '#1f324a', 'mod_p' = '#1561b0', 'base_model' = '#5fa6dc', 'Observed' = 'black')
p_labels <- c('low_p' = '30%', 'mod_p' = '50%', 'base_model' = '70%', 'Observed' = 'Observed')

p_comparison_plot <- create_comparison_plot(df = p_compare, custom_colors = p_colors, custom_labels = p_labels, titleName = 'Proportion Hm mosquitoes biting NHP vs humans') + theme(legend.position = c(0.9, 0.6))

# Compare seasonal monkey movement rates
move_compare <- do.call(rbind, list(base_model, mod_move, high_move))
move_compare <- move_compare[!duplicated(move_compare), ]
move_compare$model <- factor(move_compare$model, levels = c('high_move', 'mod_move', 'base_model', 'Observed')) 

move_colors <- c('base_model' = '#1f324a', 'mod_move' = '#1561b0', 'high_move' = '#5fa6dc', 'Observed' = 'black')
move_labels <- c('base_model' = '0.5%', 'mod_move' = '2%', 'high_move' = '10%', 'Observed' = 'Observed')

move_comparison_plot <- create_comparison_plot(df = move_compare, custom_colors = move_colors, custom_labels = move_labels, titleName = 'Movement rate of monkeys from forest to city (%/month)') + theme(legend.position = c(0.9, 0.6))

# combine plots
sensitivity_plot <- ggarrange(mu_comparison_plot, p_comparison_plot, move_comparison_plot, ncol = 1)
ggsave(filename = '../Figures/Sensitivity_plot.pdf', sensitivity_plot, width = 12, height = 9)

# Interventions 
# add base_model
int_compare <- do.call(rbind, list(reduce_mosquitoes, reduce_NHP_movement, shift_vax, combined_interventions))
int_compare <- int_compare[!duplicated(int_compare), ]

# Define custom colors
int_comp_colors <- c('reduce_mosquitoes' = '#ff595e'
                     , 'reduce_NHP_movement' = '#ffca3a'
                     , 'shift_vax' = '#80d819'
                     , 'combined_interventions' = '#1f324a'
                     , 'Observed' = 'black')

# Custom labels for the legend
int_comp_labels <- c('reduce_mosquitoes' = 'Vector control'
                     , 'reduce_NHP_movement' = 'Limit NHP movement into city'
                     , 'shift_vax' = 'Start vaccination earlier'
                     , 'combined_interventions' = 'Combined mosquito and vaccination control'
                     , 'Observed' = 'Observed')

intervention_comparison_plot <- create_comparison_plot(df = int_compare, custom_colors = int_comp_colors, custom_labels = int_comp_labels)
ggsave(filename = '../Figures/Intervention_comparison_plot_long.pdf', plot = intervention_comparison_plot, width = 10, height = 5.5)
