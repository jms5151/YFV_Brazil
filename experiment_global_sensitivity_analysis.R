library(lhs)
library(deSolve)
library(dplyr)
library(ggplot2)
library(tidyr)

# format validation data
val_data <- read.csv('../validation_data.csv')
val_data$Date <- as.Date(val_data$Month, '%m/%d/%Y')

# ============================================================================
# STEP 1: Define parameter ranges for Latin Hypercube Sampling
# ============================================================================

source('parameters.R')
source('model.R')
source('state_variables.R')

# simulation dates
start_date <- as.Date(IC$nonnumeric_value[IC$variable == 'start_date'], '%m/%d/%Y')

# Step 1: Start with your full parameter list
params_baseline <- yfv_params_list$full_model

# Step 2: Identify seasonal vs non-seasonal parameters
seasonal_params <- c('a', 'K', 'birth_hg', 'birth_aa', 'monkey_seed', 'marmoset_seed', 'hg_seed')

p_seasonal <- p[p$variable %in% seasonal_params, ]
p_nonseasonal <- p[!(p$variable %in% seasonal_params), ]

# Step 3: Create param_ranges list
param_ranges <- list()

# Add non-seasonal parameters (from your spreadsheet rows you want to vary)
# Adjust idrange_nonseasonal to match the rows you want
idrange_nonseasonal <- which(p$variable %in% p_nonseasonal$variable & 
                               !(p$variable %in% c('sigma', 'initial_v', 'final_v'))) # exclude any you don't want

for (i in idrange_nonseasonal) {
  param_ranges[[p$variable[i]]] <- c(
    p$low_value[i],
    p$value[i] * 0.9, 
    p$value[i], 
    p$value[i] * 1.1 
    # p$high_value[i]
  )
}

# Add seasonal parameters as amplitudes
for (i in 1:nrow(p_seasonal)) {
  param_name <- p_seasonal$variable[i]
  
  baseline_mean <- (p_seasonal$high_value[i] + p_seasonal$low_value[i]) / 2
  baseline_amplitude <- (p_seasonal$high_value[i] - p_seasonal$low_value[i]) / 2
  
  # Vary amplitude (seasonality strength)
  param_ranges[[paste0(param_name, "_amplitude")]] <- c(
    baseline_amplitude * 0.9,   # 50% weaker seasonality
    baseline_amplitude,          # baseline seasonality
    baseline_amplitude * 1.1    # 50% stronger seasonality
  )
}

# ============================================================================
# STEP 2: Generate Latin Hypercube Sample
# ============================================================================

n_samples <- 1000

# Generate LHS (values between 0 and 1)
lhs_design <- randomLHS(n_samples, length(param_ranges))

# Scale to actual parameter ranges
param_samples <- as.data.frame(lhs_design)
colnames(param_samples) <- names(param_ranges)

for (i in 1:length(param_ranges)) {
  param_name <- names(param_ranges)[i]
  min_val <- param_ranges[[i]][1]
  max_val <- param_ranges[[i]][3]
  
  # Scale from [0,1] to [min, max]
  param_samples[[param_name]] <- min_val + param_samples[[param_name]] * (max_val - min_val)
}

# ============================================================================
# STEP 3: Run model for each parameter set
# ============================================================================

# Your model function wrapper
run_model_for_prcc <- function(params_row, initial_state, times, p_seasonal, p_nonseasonal, 
                               sharpened_season, trend) {
  # Start with baseline parameter list
  params_list <- params_baseline
  
  # Update non-seasonal parameters directly
  for (i in 1:nrow(p_nonseasonal)) {
    param_name <- p_nonseasonal$variable[i]
    if (param_name %in% names(params_row)) {
      params_list[[param_name]] <- params_row[[param_name]]
    }
  }
  
  # Reconstruct seasonal forcing functions from amplitude
  for (i in 1:nrow(p_seasonal)) {
    sparam <- p_seasonal$variable[i]
    amplitude_name <- paste0(sparam, "_amplitude")
    
    if (amplitude_name %in% names(params_row)) {
      # Get the varied amplitude from LHS sample
      amplitude <- params_row[[amplitude_name]]
      
      # Use fixed baseline mean (not varied)
      baseline_mean <- (p_seasonal$high_value[i] + p_seasonal$low_value[i]) / 2
      
      # Calculate high and low from mean and amplitude
      high <- baseline_mean + amplitude
      low <- baseline_mean - amplitude
      
      # Recreate the seasonal forcing based on parameter type
      if (sparam %in% c('monkey_seed', 'marmoset_seed', 'hg_seed')) {
        # These use sharpened_season * trend
        base_value <- (high + low) / 2  # Use mean value
        seasonal_vals <- base_value * sharpened_season * trend
        params_list[[sparam]] <- approxfun(times, seasonal_vals)
        
      } else if (sparam == 'a') {
        # This is alpha/biting rate - simple seasonal forcing
        a_vals <- seasonal_forcing(times, high = high, low = low)
        params_list[[sparam]] <- approxfun(times, a_vals)
        
      } else if (sparam == 'K') {
        # Carrying capacity
        k_vals <- seasonal_forcing(times, high = high, low = low)
        params_list$K <- approxfun(times, k_vals)
        
      } else if (sparam == 'birth_hg') {
        # Birth rate for hg mosquitoes
        br_vals <- seasonal_forcing(times, high = high, low = low)
        params_list$birth_hg <- approxfun(times, br_vals)
        
      } else if (sparam == 'birth_aa') {
        # Birth rate for aa mosquitoes
        br_vals <- seasonal_forcing(times, high = high, low = low)
        params_list$birth_aa <- approxfun(times, br_vals)
      }
    }
  }
  
  # Remove rainy_windows
  params_list$rainy_windows <- NULL
  
  # Set up event
  event_setting <- list(
    func = event_function_reduce_mosquitoes,
    time = 15
  )
  
  # Run model
  result <- as.data.frame(
    ode(
      y = unlist(initial_state),
      times = times,
      func = yfv_model,
      parms = params_list,
      events = event_setting
    )
  )
  
  return(result)
}

# Function to calculate fit metric
calculate_fit_metric <- function(model_output, pred_var, obs_var, testmetric) {
  predicted <- model_output[, pred_var]  
  observed_data <- model_output[, obs_var]
  
  if(testmetric == 'NRMSE'){
    rmse <- sqrt(mean((predicted - observed_data)^2))
    nrmse <- rmse / mean(observed_data)
    output <- nrmse
  } else if(testmetric == 'correlation'){
    correlation <- cor.test(predicted, observed_data)
    output <- unname(round(correlation$estimate, 2))
  }
  
  return(output)
}

# Run all simulations
cat("Running", n_samples, "model simulations...\n")

# Initialize storage for results
output_metrics <- data.frame(
  sample_id = 1:n_samples,
  NRMSE_humans = NA,
  Correlation_humans = NA,
  NRMSE_primates = NA,
  Correlation_primates = NA
)

# Initial state
initial_state <- state_start_list[[1]]

# Run simulations
for (i in 1:n_samples) {
  if (i %% 100 == 0) cat("Completed", i, "of", n_samples, "simulations\n")
  
  # Get parameter set
  params_i <- param_samples[i, ]
  
  # Run model - pass sharpened_season and trend
  result <- run_model_for_prcc(params_i, initial_state, times, p_seasonal, p_nonseasonal,
                               sharpened_season, trend)
  
  # Format output
  result$Date <- start_date + result$time
  result$Ih <- lead(result$Ih, n = 30)
  
  # Join with validation data
  x <- result %>% inner_join(val_data, by = "Date") %>% drop_na() %>% as.data.frame()
  
  # Calculate metrics
  output_metrics$NRMSE_humans[i] <- calculate_fit_metric(x, 'Ih', 'MG_human', 'NRMSE')
  output_metrics$Correlation_humans[i] <- calculate_fit_metric(x, 'Ih', 'MG_human', 'correlation')
  output_metrics$NRMSE_primates[i] <- calculate_fit_metric(x, 'Ip', 'MG_primate', 'NRMSE')
  output_metrics$Correlation_primates[i] <- calculate_fit_metric(x, 'Ip', 'MG_primate', 'correlation')
}

write.csv(output_metrics, '../sensitivity_analysis_output.csv', row.names = FALSE)
write.csv(param_samples, '../sensitivity_analysis_params.csv', row.names = FALSE)

# ============================================================================
# STEP 4: Calculate PRCC
# ============================================================================

calculate_prcc <- function(params_df, output_vector) {
  # Combine parameters and output
  data_combined <- cbind(params_df, output = output_vector)
  
  # Calculate ranks
  data_ranked <- as.data.frame(apply(data_combined, 2, rank))
  
  # Initialize results
  prcc_results <- data.frame(
    parameter = names(params_df),
    PRCC = NA,
    p_value = NA
  )
  
  # Calculate partial correlation for each parameter
  for (i in 1:ncol(params_df)) {
    param_name <- names(params_df)[i]
    
    # Other parameters (all except current one)
    other_params <- setdiff(names(params_df), param_name)
    
    # Residuals from linear model: parameter ~ other parameters
    lm_param <- lm(data_ranked[[param_name]] ~ ., 
                   data = data_ranked[, other_params, drop = FALSE])
    resid_param <- residuals(lm_param)
    
    # Residuals from linear model: output ~ other parameters
    lm_output <- lm(data_ranked$output ~ ., 
                    data = data_ranked[, other_params, drop = FALSE])
    resid_output <- residuals(lm_output)
    
    # PRCC is correlation between residuals
    prcc_test <- cor.test(resid_param, resid_output)
    
    prcc_results$PRCC[i] <- prcc_test$estimate
    prcc_results$p_value[i] <- prcc_test$p.value
  }
  
  # Add significance indicator
  prcc_results$significant <- prcc_results$p_value < 0.05
  
  return(prcc_results)
}

# Calculate PRCC for each output metric
prcc_humans_nrmse <- calculate_prcc(param_samples, output_metrics$NRMSE_humans)
prcc_humans_corr <- calculate_prcc(param_samples, output_metrics$Correlation_humans)
prcc_primates_nrmse <- calculate_prcc(param_samples, output_metrics$NRMSE_primates)
prcc_primates_corr <- calculate_prcc(param_samples, output_metrics$Correlation_primates)

# Save PRCC results
write.csv(prcc_humans_nrmse, '../prcc_humans_nrmse.csv', row.names = FALSE)
write.csv(prcc_humans_corr, '../prcc_humans_corr.csv', row.names = FALSE)
write.csv(prcc_primates_nrmse, '../prcc_primates_nrmse.csv', row.names = FALSE)
write.csv(prcc_primates_corr, '../prcc_primates_corr.csv', row.names = FALSE)

# ============================================================================
# STEP 5: Create visualizations
# ============================================================================

# Function to format parameter names (reuse from tornado plot)
format_parameter_names <- function(param_names) {
  sapply(param_names, function(p) {
    if (grepl("^mu_", p)) {
      subscript <- sub("^mu_", "", p)
      return(bquote(mu[.(subscript)]))
    }
    else if (grepl("^gamma_", p)) {
      subscript <- sub("^gamma_", "", p)
      return(bquote(gamma[.(subscript)]))
    }
    else if (grepl("^delta_", p)) {
      subscript <- sub("^delta_", "", p)
      return(bquote(delta[.(subscript)]))
    }
    else if (grepl("^pMI", p)) {
      subscript <- sub("^pMI", "", p)
      return(bquote(pMI[.(subscript)]))
    }
    else if (grepl("^PDR_", p)) {
      subscript <- sub("^PDR_", "", p)
      return(bquote(PDR[.(subscript)]))
    }
    else if (p == "b") {
      return(bquote(b))
    }
    else {
      return(p)
    }
  }, USE.NAMES = FALSE)
}

# ============================================================================
# VISUALIZATION 1: Bar plot with significance
# ============================================================================

plot_prcc_barplot <- function(prcc_df, title) {
  # Sort by absolute PRCC value
  prcc_df <- prcc_df %>%
    arrange(abs(PRCC)) %>%
    mutate(parameter = factor(parameter, levels = parameter))
  
  # Format parameter names
  formatted_labels <- format_parameter_names(levels(prcc_df$parameter))
  
  # Create plot
  p <- ggplot(prcc_df, aes(x = parameter, y = PRCC, fill = significant)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
    geom_hline(yintercept = c(-0.3, 0.3), linetype = "dashed", color = "gray50", linewidth = 0.3) +
    scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "gray70"),
                      labels = c("Not significant", "p < 0.05"),
                      name = "") +
    scale_x_discrete(labels = formatted_labels) +
    coord_flip() +
    labs(
      title = title,
      subtitle = paste("Based on", n_samples, "Latin Hypercube Samples"),
      x = "Parameter",
      y = "Partial Rank Correlation Coefficient (PRCC)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )
  
  return(p)
}

# Create individual plots
p1 <- plot_prcc_barplot(prcc_humans_nrmse, "PRCC: NRMSE (Humans)")
p2 <- plot_prcc_barplot(prcc_humans_corr, "PRCC: Correlation (Humans)")
p3 <- plot_prcc_barplot(prcc_primates_nrmse, "PRCC: NRMSE (Primates)")
p4 <- plot_prcc_barplot(prcc_primates_corr, "PRCC: Correlation (Primates)")

library(gridExtra)
prcc_plots <- grid.arrange(p1, p2, p3, p4, ncol = 2)

ggsave(filename = '../Figures/PRCC.pdf', plot = prcc_plots, width = 10, height = 8)
# ============================================================================
# VISUALIZATION 2: Heatmap showing all metrics at once
# ============================================================================

create_prcc_heatmap <- function() {
  # Combine all PRCC results
  prcc_combined <- data.frame(
    parameter = prcc_humans_nrmse$parameter,
    NRMSE_Humans = prcc_humans_nrmse$PRCC,
    Corr_Humans = prcc_humans_corr$PRCC,
    NRMSE_Primates = prcc_primates_nrmse$PRCC,
    Corr_Primates = prcc_primates_corr$PRCC
  )
  
  # Reshape to long format
  prcc_long <- prcc_combined %>%
    pivot_longer(cols = -parameter, names_to = "metric", values_to = "PRCC")
  
  # Order parameters by average absolute PRCC
  param_order <- prcc_combined %>%
    mutate(avg_abs_prcc = rowMeans(abs(.[, -1]))) %>%
    arrange(desc(avg_abs_prcc)) %>%
    pull(parameter)
  
  prcc_long$parameter <- factor(prcc_long$parameter, levels = param_order)
  
  # Format parameter names
  formatted_labels <- format_parameter_names(levels(prcc_long$parameter))
  
  # Create heatmap
  p <- ggplot(prcc_long, aes(x = metric, y = parameter, fill = PRCC)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", PRCC)), size = 3) +
    scale_fill_gradient2(
      low = "darkred", 
      mid = "white", 
      high = "darkblue",
      midpoint = 0,
      limits = c(-1, 1),
      name = "PRCC"
    ) +
    scale_y_discrete(labels = formatted_labels) +
    scale_x_discrete(labels = c("NRMSE\nHumans", "Correlation\nHumans", 
                                "NRMSE\nPrimates", "Correlation\nPrimates")) +
    labs(
      title = "Global Sensitivity Analysis: PRCC Heatmap",
      subtitle = paste("Based on", n_samples, "Latin Hypercube Samples"),
      x = "Output Metric",
      y = "Parameter"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      panel.grid = element_blank()
    )
  
  return(p)
}

heatmap_plot <- create_prcc_heatmap()
print(heatmap_plot)
ggsave(filename = '../Figures/PRCC_heatmap.pdf', plot = heatmap_plot, width = 6, height = 8)

# ============================================================================
# VISUALIZATION 3: Summary table
# ============================================================================

create_prcc_summary_table <- function() {
  # Combine all results
  summary_table <- data.frame(
    Parameter = prcc_humans_nrmse$parameter,
    NRMSE_Humans = sprintf("%.3f%s", prcc_humans_nrmse$PRCC, 
                           ifelse(prcc_humans_nrmse$significant, "*", "")),
    Corr_Humans = sprintf("%.3f%s", prcc_humans_corr$PRCC,
                          ifelse(prcc_humans_corr$significant, "*", "")),
    NRMSE_Primates = sprintf("%.3f%s", prcc_primates_nrmse$PRCC,
                             ifelse(prcc_primates_nrmse$significant, "*", "")),
    Corr_Primates = sprintf("%.3f%s", prcc_primates_corr$PRCC,
                            ifelse(prcc_primates_corr$significant, "*", ""))
  )
  
  # Order by average absolute PRCC
  prcc_combined <- data.frame(
    parameter = prcc_humans_nrmse$parameter,
    avg_abs = rowMeans(abs(cbind(
      prcc_humans_nrmse$PRCC,
      prcc_humans_corr$PRCC,
      prcc_primates_nrmse$PRCC,
      prcc_primates_corr$PRCC
    )))
  )
  
  summary_table <- summary_table[order(-prcc_combined$avg_abs), ]
  
  return(summary_table)
}

prcc_table <- create_prcc_summary_table()
cat("\n=== PRCC Summary Table (* = p < 0.05) ===\n")
print(prcc_table, row.names = FALSE)

# ============================================================================
# Export results
# ============================================================================

# Save plots
# ggsave("prcc_barplot_humans_nrmse.png", p1, width = 8, height = 6, dpi = 300)
# ggsave("prcc_heatmap.png", heatmap_plot, width = 8, height = 8, dpi = 300)

# Save table
write.csv(prcc_table, "prcc_summary_table.csv", row.names = FALSE)


