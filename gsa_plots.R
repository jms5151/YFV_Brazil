library(ggplot2)
library(dplyr)

# Your sensitivity results data
sensitivity_results <- read.csv('../sensitivity_analysis_output.csv')

baseline_results <- read.csv('model_validation.csv')

baseline_NRMSE_humans <- baseline_results$nrmse_humans[1]
baseline_Correlation_humans <- baseline_results$correlation_humans[1]  
baseline_NRMSE_primates <- baseline_results$nrmse_primates[1]
baseline_Correlation_primates <- baseline_results$correlation_primates[1]

# ============================================================================
# Prepare data for tornado plot
# ============================================================================

# Choose which metric to plot (change this line to switch metrics)
metric_column <- "NRMSE_humans"  # Options: "NRMSE_humans", "Correlation_humans", 
# "NRMSE_primates", "Corrleation_primates"

# Get baseline value
baseline_value <- switch(metric_column,
                         "NRMSE_humans" = baseline_NRMSE_humans,
                         "Correlation_humans" = baseline_Correlation_humans,
                         "NRMSE_primates" = baseline_NRMSE_primates,
                         "Corrleation_primates" = baseline_Correlation_primates
)

# Reshape data to get low and high values for each parameter
tornado_data <- sensitivity_results %>%
  group_by(parameter) %>%
  summarise(
    value_low = min(.data[[metric_column]]),
    value_high = max(.data[[metric_column]]),
    range = abs(value_high - value_low)
  ) %>%
  arrange(desc(range)) %>%
  mutate(parameter = factor(parameter, levels = parameter))

# ============================================================================
# Create tornado plot
# ============================================================================

# Set axis label based on metric
x_label <- switch(metric_column,
                  "NRMSE_humans" = "NRMSE - Humans (lower = better fit)",
                  "Correlation_humans" = "Correlation - Humans (higher = better fit)",
                  "NRMSE_primates" = "NRMSE - Primates (lower = better fit)",
                  "Corrleation_primates" = "Correlation - Primates (higher = better fit)"
)

# Create the plot
tornado_plot <- ggplot(tornado_data, aes(y = parameter)) +
  geom_segment(
    aes(x = value_low, xend = value_high, 
        y = parameter, yend = parameter),
    linewidth = 10, 
    color = "steelblue", 
    alpha = 0.7
  ) +
  geom_vline(
    xintercept = baseline_value, 
    linetype = "dashed", 
    color = "red", 
    linewidth = 1
  ) +
  annotate(
    "text", 
    x = baseline_value, 
    y = nrow(tornado_data) + 0.5,
    label = "Baseline", 
    color = "red", 
    hjust = -0.1, 
    size = 4
  ) +
  labs(
    title = "Tornado Plot: Parameter Sensitivity Analysis",
    subtitle = paste("Metric:", gsub("_", " ", metric_column), "| Variation: ±20%"),
    x = x_label,
    y = "Parameter"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

print(tornado_plot)

# ============================================================================
# Print summary table
# ============================================================================

cat("\n=== Sensitivity Analysis Summary ===\n")
cat("Metric:", metric_column, "\n")
cat("Baseline value:", round(baseline_value, 3), "\n\n")

summary_table <- tornado_data %>%
  mutate(
    value_low = round(value_low, 3),
    value_high = round(value_high, 3),
    range = round(range, 3)
  )

print(summary_table)

# ============================================================================
# OPTIONAL: Create all four tornado plots at once
# ============================================================================

create_all_tornado_plots <- function(sensitivity_results, 
                                     baseline_NRMSE_humans,
                                     baseline_Correlation_humans,
                                     baseline_NRMSE_primates,
                                     baseline_Correlation_primates) {
  
  library(gridExtra)
  
  metrics <- list(
    list(col = "NRMSE_humans", baseline = baseline_NRMSE_humans, 
         label = "NRMSE - Humans\n(lower = better)"),
    list(col = "Correlation_humans", baseline = baseline_Correlation_humans, 
         label = "Correlation - Humans\n(higher = better)"),
    list(col = "NRMSE_primates", baseline = baseline_NRMSE_primates, 
         label = "NRMSE - Primates\n(lower = better)"),
    list(col = "Corrleation_primates", baseline = baseline_Correlation_primates, 
         label = "Correlation - Primates\n(higher = better)")
  )
  
  plot_list <- lapply(metrics, function(m) {
    tornado_data <- sensitivity_results %>%
      group_by(parameter) %>%
      summarise(
        value_low = min(.data[[m$col]]),
        value_high = max(.data[[m$col]]),
        range = abs(value_high - value_low)
      ) %>%
      arrange(desc(range)) %>%
      mutate(parameter = factor(parameter, levels = parameter))
    
    ggplot(tornado_data, aes(y = parameter)) +
      geom_segment(
        aes(x = value_low, xend = value_high, y = parameter, yend = parameter),
        linewidth = 8, color = "steelblue", alpha = 0.7
      ) +
      geom_vline(xintercept = m$baseline, linetype = "dashed", 
                 color = "red", linewidth = 0.8) +
      labs(x = m$label, y = NULL) +
      theme_minimal(base_size = 10) +
      theme(
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 9)
      )
  })
  
  grid.arrange(grobs = plot_list, ncol = 2, 
               top = "Tornado Plots: Parameter Sensitivity Analysis (±20% variation)")
}

# Uncomment to create all four plots:
create_all_tornado_plots(sensitivity_results,
                         baseline_NRMSE_humans,
                         baseline_Correlation_humans,
                         baseline_NRMSE_primates,
                         baseline_Correlation_primates)
