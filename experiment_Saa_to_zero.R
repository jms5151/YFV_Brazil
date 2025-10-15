# load the required packages
library(foreach)
library(doParallel)
library(deSolve)
library(tidyverse)

# model
source('model.R')

# lists of parameters
source('parameters.R')

# lists of initial conditions/state variables
source('state_variables.R')

# make Saa zero for all starting conditions
state_start_list <- lapply(state_start_list, function(x) {
  if (!is.null(x$Saa)) x$Saa <- 0
  x
})

# Set up parallel backend to use multiple processors
numCores <- detectCores() - 1
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Run the ODE model in parallel for each combination of state_start_list and yfv_params_list
Saa0 <- foreach(yfv_params_idx = 1, .packages = 'deSolve') %:%
  foreach(state_start_idx = 1:length(state_start_list), .packages = 'deSolve') %dopar% {
    state_start <- state_start_list[[state_start_idx]]
    yfv_params <- yfv_params_list[[yfv_params_idx]]
    times <- times_list[[yfv_params_idx]]
    
    if (!is.null(yfv_params$rainy_windows)) {
      rainy_days <- unlist(lapply(yfv_params$rainy_windows, function(w) w$start_day:w$end_day))
      event_times_final <- sort(unique(c(15, rainy_days)))
    } else {
      event_times_final <- 15
    }
    
    event_setting <- list(
      func = event_function_reduce_mosquitoes,
      time = event_times_final
    )
    
    result <- as.data.frame(
      ode(
        y = unlist(state_start),  # Ensure state_start is a numeric vector
        times = times,
        func = yfv_model,
        parms = yfv_params,
        events = event_setting
      )
    )
    
    list(
      yfv_params_idx = yfv_params_idx,
      state_start_idx = state_start_idx,
      result = result
    )
  }

# Stop the parallel backend
stopCluster(cl)

# save
saveRDS(object = Saa0, file = '../Saa0_results.RData')

# assess results
# reporting rates
rho_humans = 0.45 # symptomatic rate, from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4632853/#:~:text=For%20YF%20virus%20infections%2C%20the,%25%20CI%200.31%E2%80%930.62
rho_monkeys = p$value[p$variable == 'p']

# format validation data
val_data <- read.csv('../validation_data.csv')
val_data$Date <- as.Date(val_data$Month, '%m/%d/%Y')

# format simulated data
get_medians <- function(dflist){
  df <- dflist
  
  # Extract the third element from each list within grouped_results
  third_elements <- lapply(df, function(x) x[[3]])
  
  # Combine all third elements into a single data frame using do.call and rbind
  combined_df <- do.call(rbind, third_elements)
  
  # summarise
  df2 <- combined_df[[3]] %>%
    drop_na(Ih, Ip) %>% #, Ic
    group_by(time) %>%
    summarise(
      across(c(Ih, Ip), list(
        median = ~ quantile(.x, 0.5, na.rm = TRUE),
        p25 = ~ quantile(.x, 0.25, na.rm = TRUE),
        p75 = ~ quantile(.x, 0.75, na.rm = TRUE)
      ))
    ) %>%
    mutate('model' = names(yfv_params_list)[i])
  
  df2$Date <- yfv_epidemic[1:nrow(df2)]
  
  # lead modeled human cases by 1 month
  df2[,c('Ih_median', 'Ih_p25', 'Ih_p75')] <- lapply(df2[,c('Ih_median', 'Ih_p25', 'Ih_p75')], function(x) lead(x, n = 30))
  
  return(df2)
}


# calculate correlations and NRMSE between simulations and observations
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
  meanObs <- mean(x3[,var2], na.rm = TRUE)
  nrmse <- round(sqrt(mean((x3[,var1] - x3[,var2])^2, na.rm = TRUE))/meanObs)
  # create output
  dist_out <- data.frame('corr' = corVal, 'pvalue' = corPval, 'nrmse' = nrmse)
  return(dist_out)
}

x <- get_medians(dflist = Saa0)
corHuman <- calc_correlation(df = x, rho = rho_humans, var1 = 'Ih_median', var2 = 'MG_human')
corPrimates <- calc_correlation(df = x, rho = rho_monkeys, var1 = 'Ip_median', var2 = 'MG_primate')
corSaa0 <- data.frame(corHuman, corPrimates)
colnames(corSaa0) <- c('correlation_humans', 'pvalue_humans', 'nrmse_humans', 'correlation_primates', 'pvalue_primates', 'nrmse_primates')
write.csv(corSaa0, 'Saa0_validation.csv', row.names = F)
