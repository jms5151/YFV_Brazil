# load librariesMore actions
library(tidyverse)
library(ggbreak)
library(patchwork)
library(ggpubr)

# load model results
resultsNew <- readRDS('../model_results.RData')

# source parameters code to bring some variables into environment
source('parameters.R')

# load start_date
IC <- read.csv('init_conditions.csv')
start_date <- as.Date(IC$nonnumeric_value[IC$variable == 'start_date'], '%m/%d/%Y')
end_date <- as.Date(IC$nonnumeric_value[IC$variable == 'end_date'], '%m/%d/%Y')

# reporting rates
rho_humans = 0.45 # symptomatic rate, from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4632853/#:~:text=For%20YF%20virus%20infections%2C%20the,%25%20CI%200.31%E2%80%930.62
rho_monkeys = p$value[p$variable == 'p']

# format validation data
# allow for observations to extend for the interventions
val_data <- read.csv('../validation_data.csv')
val_data$Date <- as.Date(val_data$Month, '%m/%d/%Y')
val_data_long <- val_data %>%
  filter(Date >= start_date) %>% 
  pivot_longer(
    cols = c(MG_human, MG_primate),
    names_to = "variable",
    values_to = "value"
  ) %>%
  as.data.frame() 

val_data_long$variable <- ifelse(val_data_long$variable == 'MG_human', 'Ih_median', 'Ip_median')  
val_data_long$model <- 'Observed'
val_data_long$I_p25 <- NA
val_data_long$I_p75 <- NA

# format simulated data
get_medians <- function(dflist){
  df <- dflist
  
  # Extract the third element from each list within grouped_results
  third_elements <- lapply(df, function(x) x[[3]])
  
  # Combine all third elements into a single data frame using do.call and rbind
  combined_df <- do.call(rbind, third_elements)
  
  # summarise
  df2 <- combined_df %>%
    drop_na(Ih, Ip) %>% #, Ic
    # mutate('Ip' = Ip + Ic) %>%
    group_by(time) %>%
    summarise(
      across(c(Ih, Ip), list(
        median = ~ quantile(.x, 0.5, na.rm = TRUE),
        p25 = ~ quantile(.x, 0.25, na.rm = TRUE),
        p75 = ~ quantile(.x, 0.75, na.rm = TRUE)
      ))
    ) %>%
    mutate('model' = names(yfv_params_list)[i])
  
  if(grepl('reduce|limit|shift|combined', unique(df2$model))==T){
    end_date_new <- start_date + (nrow(df2) - 1)
    df2$Date <- seq.Date(from = start_date, to = end_date_new, by = 1)
  } else {
    df2$Date <- yfv_epidemic[1:nrow(df2)]
  }
  
  # lead modeled human cases by 1 month
  df2[,c('Ih_median', 'Ih_p25', 'Ih_p75')] <- lapply(df2[,c('Ih_median', 'Ih_p25', 'Ih_p75')], function(x) lead(x, n = 30))
  
  return(df2)
}

# calculate correlations between simulations and observations
calc_correlation <- function(df, rho, var1, var2){
  # add observed data
  x2 <- df %>% left_join(val_data) %>% as.data.frame()
  x2[, var1] <- x2[, var1] * rho
  x3 <- x2[,c(var1, var2)]
  x3 <- subset(x3, !is.na(x3[,var2]))
  # calculate correlation and p-value
  corrsum <- cor.test(x3[,var1], x3[,var2])
  corVal <- round(corrsum$estimate, 2)
  corPval <- round(corrsum$p.value, 3)
  # create output
  dist_out <- data.frame('corr' = corVal, 'pvalue' = corPval)
  return(dist_out)
}

corDF <- data.frame(matrix(nrow = 0, ncol = 5))
colnames(corDF) <- c('model', 'correlation_humans', 'pvalue_humans', 'correlation_primates', 'pvalue_primates')

for(i in 1:10){
  x <- get_medians(dflist = resultsNew[[i]])
  corHuman <- calc_correlation(df = x, rho = rho_humans, var1 = 'Ih_median', var2 = 'MG_human')
  corPrimates <- calc_correlation(df = x, rho = rho_monkeys, var1 = 'Ip_median', var2 = 'MG_primate')
  corDF[i,] <- c(unique(x$model), corHuman, corPrimates)
}

corDF[,c('correlation_humans', 'correlation_primates')] <- lapply(corDF[,c('correlation_humans', 'correlation_primates')], function(x) as.numeric(x))
corDF$average_correlation <- rowMeans(corDF[,c('correlation_humans', 'correlation_primates')])

write.csv(corDF, 'model_validation.csv', row.names = F)

# format simulated data for plotting
for(i in 1:length(yfv_params_list)){ #
  
  df2 <- get_medians(dflist = resultsNew[[i]])
  
  # transform to good format for plotting
  df_long_median <- df2 %>%
    pivot_longer(
      cols = c(Ih_median, Ip_median),
      names_to = "variable",
      values_to = "value"
    )
  
  # Create a key to join the percentiles back to the long format
  df_long_median <- df_long_median %>%
    mutate(
      I_p25 = case_when(
        variable == "Ih_median" ~ Ih_p25,
        variable == "Ip_median" ~ Ip_p25
      ),
      I_p75 = case_when(
        variable == "Ih_median" ~ Ih_p75,
        variable == "Ip_median" ~ Ip_p75
      )
    ) %>%
    select(Date, model, variable, value, I_p25, I_p75) %>%
    as.data.frame()
  
  # add observed data
  if(grepl('reduce|limit|shift|combined', unique(df2$model))==T){
    df_long_median <- rbind(df_long_median, val_data_long[,colnames(df_long_median)])
  } else {
    x <- val_data_long %>%
      filter(Date <= end_date) 
    df_long_median <- rbind(df_long_median, x[,colnames(df_long_median)])
  }
  
  # rename variable labels
  df_long_median$variable <- ifelse(df_long_median$variable == 'Ih_median', 'Infected people', 'Infected monkeys')
  
  # rename dataset
  assign(names(yfv_params_list)[i], df_long_median)
  
}

# Plotting function
panelPlot <- function(df_model, df_observed_expanded, custom_colors, custom_labels, titleName = ''){
  p <- ggplot(df_model, aes(x = Date, y = value, color = model, fill = model)) +
    geom_ribbon(aes(ymin = I_p25, ymax = I_p75), alpha = 0.3, linetype = 0) +
    geom_line(lwd = 1.1) +
    geom_line(data = df_observed_expanded, aes(x = Date, y = value), lwd = 0.9, inherit.aes = FALSE) +
    geom_point(data = df_observed_expanded, aes(x = Date, y = value), size = 2, inherit.aes = FALSE) +
    facet_grid(model ~ variable) + # ~ variable
    theme_bw() +
    xlab('Date') +
    ylab('Infected') +
    scale_color_manual(values = custom_colors, labels = custom_labels, name = 'Model', breaks = names(custom_labels)) +
    scale_fill_manual(values = custom_colors, labels = custom_labels, name = 'Model', breaks = names(custom_labels)) +
    theme(
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      strip.text = element_text(size = 12),
      strip.text.y = element_blank()
    ) +
    theme(legend.position = 'bottom') +
    guides(
      color = guide_legend(
        ncol = 1,
        title = titleName,
        override.aes = list(
          linetype = rep("solid", length(custom_labels)),
          shape = c(rep(NA, length(custom_labels)))  # Adjust based on legend requirements
        )
      ),
      fill = guide_legend(ncol = 1, title = titleName)
    ) +
    ggtitle(titleName)
  
  return(p)
}

overlaidPlot <- function(df_model, df_observed, titleName = '', custom_colors, custom_labels, plotType = 0){
  p <- ggplot(df_model, aes(x = Date, y = value, color = model, fill = model)) +
    geom_ribbon(aes(ymin = I_p25, ymax = I_p75), alpha = 0.3, linetype = 0) +
    geom_line(lwd = 1.1) +
    geom_line(data = df_observed, aes(x = Date, y = value), lwd = 0.9) +
    geom_point(data = df_observed, aes(x = Date, y = value), size = 2) +
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
    guides(color = guide_legend(ncol = 1)
           , override.aes = list(linetype = c(rep('solid', length(custom_labels))))
           , fill = guide_legend(ncol = 1)
    ) +
    ggtitle(titleName)
  
  return(p)
  
}

create_comparison_plot <- function(df, custom_colors, custom_labels, titleName = '', plotType = 1) {
  df[df$model != 'Observed' & df$variable == 'Infected people', c('value', 'I_p25', 'I_p75')] <- lapply(df[df$model != 'Observed' & df$variable == 'Infected people', c('value', 'I_p25', 'I_p75')], function(x) x * rho_humans)
  df[df$model != 'Observed' & df$variable == 'Infected monkeys', c('value', 'I_p25', 'I_p75')] <- lapply(df[df$model != 'Observed' & df$variable == 'Infected monkeys', c('value', 'I_p25', 'I_p75')], function(x) x * rho_monkeys)
  
  # Separate the observed data
  df_observed <- df[df$model == 'Observed', ]
  df_model <- df[df$model != 'Observed', ]
  
  # Create a grid of all combinations of `model` and `variable`
  facet_combinations <- expand.grid(
    model = unique(df_model$model),
    variable = unique(df_model$variable)
  )
  
  # Replicate df_observed for each combination
  df_observed_expanded <- merge(df_observed, facet_combinations, by = 'variable')
  
  # plot
  if(plotType == 1){
    p <- panelPlot(df_model = df_model, df_observed_expanded = df_observed_expanded, custom_colors = custom_colors, custom_labels = custom_labels)
  } else {
    p <- overlaidPlot(df_model = df_model, df_observed = df_observed, custom_colors = custom_colors, custom_labels = custom_labels)
  }
  
  return(p)
}

# Model comparison figures --------------
# combine datasets
mod_compare <- do.call(rbind, list(full_model, constant_bite_rate, no_movement, constant_bite_rate_and_no_movement))
mod_compare <- mod_compare[!duplicated(mod_compare), ]
mod_compare$model <- factor(mod_compare$model, levels = c('full_model', 'constant_bite_rate', 'no_movement', 'constant_bite_rate_and_no_movement', 'Observed'))

# Define custom colors
mod_comp_colors <- c('full_model' = '#149889'
                     , 'constant_bite_rate' = '#82118e'
                     , 'constant_bite_rate_and_no_movement' = '#c5a030'
                     , 'no_movement' = '#2e3d54'
)

# Custom labels for the legend
mod_comp_labels <- c('full_model' = 'Seasonal forest species movement & seasonal biting rate'
                     , 'constant_bite_rate_and_no_movement' = 'No seasonally-driven forest species movement or biting rate'
                     , 'constant_bite_rate' = 'Seasonally-driven  forest species movement'
                     , 'no_movement' = 'Seasonally-driven biting rate'
)

left_panel <- mod_compare  %>% filter(variable == 'Infected monkeys') 
right_panel <- mod_compare  %>% filter(variable == 'Infected people') 

left_plot <-  create_comparison_plot(df = left_panel, custom_colors = mod_comp_colors, custom_labels = mod_comp_labels, plotType = 1) + xlim(start_date, end_date - 180) + ylim(0,300)
right_plot <-  create_comparison_plot(df = right_panel, custom_colors = mod_comp_colors, custom_labels = mod_comp_labels, plotType = 1) + xlim(start_date, end_date - 180) + ylab('')

model_comparison_plot <- (left_plot + right_plot) + 
  plot_layout(guides = 'collect', axis_titles = 'collect') & 
  theme(legend.position = 'bottom')

ggsave(filename = '../Figures/Model_comparison_plot.pdf', plot = model_comparison_plot, width = 8, height = 8.5)

# Sensitivity plots ---------------------------------------------------
# Compare YFV primate mortality rates
mu_compare <- do.call(rbind, list(full_model, low_mu_v1, high_mu_v1))
mu_compare <- mu_compare[!duplicated(mu_compare), ]
mu_compare$model <- factor(mu_compare$model, levels = c('high_mu_v1', 'full_model', 'low_mu_v1', 'Observed')) 

mu_colors <- c('low_mu_v1' = '#1f324a', 'full_model' = '#1561b0', 'high_mu_v1' = '#5fa6dc', 'Observed' = 'black')
mu_labels <- c('low_mu_v1' = 'Low mortality rate (20%)', 'full_model' = 'Moderate mortality rate (50%)*', 'high_mu_v1' = 'High mortality rate (80%)', 'Observed' = 'Observed')

mu_comparison_plot <- create_comparison_plot(
  df = mu_compare
  , custom_colors = mu_colors
  , custom_labels = mu_labels
  , titleName = 'YFV-induced howler monkey mortality rate'
  , plotType = 2
  ) + 
  guides(fill = guide_legend(title = 'YFV-induced howler monkey mortality rate')
         , color = guide_legend(title = 'YFV-induced howler monkey mortality rate')) +
  theme(legend.position = c(0.16, 0.6)
        , legend.background=element_blank()) +  
  xlim(start_date, end_date - 180) +
  xlab('')

# Compare higher probability of mosquito infection rates
p_compare <- do.call(rbind, list(full_model, high_pMI))
p_compare <- p_compare[!duplicated(p_compare), ]
p_compare$model <- factor(p_compare$model, levels = c('full_model', 'high_pMI', 'Observed')) 

p_colors <- c('high_pMI' = '#1561b0', 'full_model' = '#5fa6dc', 'Observed' = 'black')
p_labels <- c('high_pMI' = '2x literature values', 'full_model' = 'Literature values*', 'Observed' = 'Observed')

p_comparison_plot <- create_comparison_plot(
  df = p_compare
  , custom_colors = p_colors
  , custom_labels = p_labels
  , titleName = 'Probability of mosquito infection with YFV'
  , plotType = 2
) + 
  guides(fill = guide_legend(title = 'Probability of mosquito infection with YFV')
         , color = guide_legend(title = 'Probability of mosquito infection with YFV')) +
  theme(legend.position = c(0.16, 0.6)
        , legend.background=element_blank()
        , strip.text = element_blank()) +
  xlim(start_date, end_date - 180) +
  xlab('')

# Compare seasonal monkey movement rates
move_compare <- do.call(rbind, list(full_model, low_movement, high_movement))
move_compare <- move_compare[!duplicated(move_compare), ]
move_compare$model <- factor(move_compare$model, levels = c('low_movement', 'full_model', 'high_movement', 'Observed')) 

move_colors <- c('low_movement' = '#1f324a', 'full_model' = '#1561b0', 'high_movement' = '#5fa6dc', 'Observed' = 'black')
move_labels <- c('low_movement' = 'Low', 'full_model' = 'Moderate*', 'high_movement' = 'High', 'Observed' = 'Observed')

move_comparison_plot <- create_comparison_plot(
  df = move_compare
  , custom_colors = move_colors
  , custom_labels = move_labels
  , titleName = 'Immigration rate of non-human primates\nand Hg mosquitos from forest to city'
  , plotType = 2
) + 
  guides(fill = guide_legend(title = 'Immigration rate of non-human primates\nand Hg mosquitos from forest to city')
         , color = guide_legend(title = 'Immigration rate of non-human primates\nand Hg mosquitos from forest to city')) +
  
  theme(legend.position = c(0.16, 0.6)
        , legend.background=element_blank()
        , strip.text = element_blank()) + 
  xlim(start_date, end_date - 180)

# combine plots
sensitivity_plot <- ggarrange(mu_comparison_plot, p_comparison_plot, move_comparison_plot, ncol = 1)
ggsave(filename = '../Figures/Sensitivity_plot.pdf', sensitivity_plot, width = 12, height = 9)

# Interventions ----------------------
int_compare <- do.call(rbind, list(full_model, reduce_mosquitoes, limit_monkey_movement, shift_vax, combined_interventions))
int_compare <- int_compare[!duplicated(int_compare), ]
int_compare <- subset(int_compare, variable == 'Infected people')
int_compare$model <- factor(int_compare$model, levels = c('Observed'
                                                          , 'full_model'
                                                          , 'limit_monkey_movement'
                                                          , 'reduce_mosquitoes'
                                                          , 'shift_vax'
                                                          , 'combined_interventions'
                                                          ))

# Define custom colors
int_comp_colors <- c('reduce_mosquitoes' = '#ff595e'
                     , 'limit_monkey_movement' = '#ffca3a'
                     , 'shift_vax' = '#80d819'
                     , 'combined_interventions' = 'purple'
                     ,  'full_model' = '#1abc9c'
                     , 'Observed' = 'black')

# Custom labels for the legend
int_comp_labels <- c('reduce_mosquitoes' = 'Vector control'
                     , 'limit_monkey_movement' = 'Conservation (limit NHP movement into city)'
                     , 'shift_vax' = 'Start vaccination earlier'
                     , 'combined_interventions' = 'All three interventions combined'
                     , 'full_model' = 'Full model'
                     , 'Observed' = 'Observed')

int_compare_short <- subset(int_compare, Date <= '2018-07-01')

intervention_comparison_plot_short <- create_comparison_plot(
  df = int_compare_short
  , custom_colors = int_comp_colors
  , custom_labels = int_comp_labels
  , plotType = 2
  ) +   
  facet_wrap(~variable, scales = 'free') +
  theme(legend.position = 'right', strip.text = element_blank()) +
  labs(title = 'Intervention comparison'
       , subtitle = 'Short-term interventions'
       , y = 'Human YFV cases'
       , x = '')

intervention_comparison_plot_long <- create_comparison_plot(
  df = int_compare
  , custom_colors = int_comp_colors
  , custom_labels = int_comp_labels
  , plotType = 2
  ) + 
  scale_y_break(c(200, 400), scales = 0.5, space = 0.05) +  # add the break
  theme(legend.position = 'none'
        , strip.text = element_blank()
  ) +
  labs(title = ''
       , subtitle = 'Long-term interventions'
       , y = 'Human YFV cases'
       , x = ''
  ) +
  ylim(0, 3000) + 
  theme(
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    axis.line.y.right = element_blank()
  )

# Intervention plot with higher R0 here
int_compare_high_R0 <- do.call(rbind, list(full_model, reduce_mosquitoes_high_R0, limit_monkey_movement_high_R0, shift_vax_high_R0, combined_interventions_high_R0))
int_compare_high_R0 <- int_compare_high_R0[!duplicated(int_compare_high_R0), ]
int_compare_high_R0 <- subset(int_compare_high_R0, variable == 'Infected people' & !is.na(value))
int_compare_high_R0$model <- gsub('_high_R0', '', int_compare_high_R0$model)
int_compare_high_R0$model <- factor(int_compare_high_R0$model, levels = c('Observed'
                                                                          , 'full_model'
                                                                          , 'limit_monkey_movement'
                                                                          , 'reduce_mosquitoes'
                                                                          , 'shift_vax'
                                                                          , 'combined_interventions'
                                                                          ))
df = int_compare_high_R0
df[df$model != 'Observed' & df$variable == 'Infected people', c('value', 'I_p25', 'I_p75')] <- lapply(df[df$model != 'Observed' & df$variable == 'Infected people', c('value', 'I_p25', 'I_p75')], function(x) x * rho_humans)


intervention_comparison_plot_high_R0 <- create_comparison_plot(
  df = int_compare_high_R0
  , custom_colors = int_comp_colors
  , custom_labels = int_comp_labels
  , plotType = 2
  ) + 
  scale_y_break(c(200,400), scales = 0.5, space = 0.05) +  # add the break
  theme(legend.position = 'none'
        , strip.text = element_blank()
  ) +
  labs(title = ''
       , subtitle = 'Long-term interventions, higher transmission rate'
       , y = ''
       , x = ''
  ) + 
  ylim(0, 3000) + 
  theme(
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    axis.line.y.right = element_blank()
  )

long_int_plots <- (intervention_comparison_plot_long + intervention_comparison_plot_high_R0)
intervention_plots <- ggarrange(intervention_comparison_plot_short, long_int_plots, ncol = 1) 

ggsave(filename = '../Figures/Intervention_comparison_plot.pdf', plot = intervention_plots, width = 9, height = 6)
