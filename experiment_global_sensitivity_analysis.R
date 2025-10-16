# load the required packages
library(foreach)
library(doParallel)
library(deSolve)
library(tidyverse)

# model
source('model.R')

# lists of parameters
source('parameters.R')

# use full model 
yfv_params_list  <- yfv_params_list[1]

# lists of initial conditions/state variables
source('state_variables.R')

# take one starting condition since it doesn't affect results
state_start_list <- state_start_list[1]

# reporting rates
rho_humans = 0.45 # symptomatic rate, from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4632853/#:~:text=For%20YF%20virus%20infections%2C%20the,%25%20CI%200.31%E2%80%930.62
rho_monkeys = p$value[p$variable == 'p']

# format validation data
val_data <- read.csv('../validation_data.csv')
val_data$Date <- as.Date(val_data$Month, '%m/%d/%Y')

# simulation dates
start_date <- as.Date(IC$nonnumeric_value[IC$variable == 'start_date'], '%m/%d/%Y')

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
  nrmse <- sqrt(mean((x3[,var1] - x3[,var2])^2, na.rm = TRUE))/meanObs
  # create output
  dist_out <- data.frame('corr' = corVal, 'pvalue' = corPval, 'nrmse' = nrmse)
  return(dist_out)
}

# adjust parameters by +/- 20% and re-run, calculate correlation and RMSE
output_df <- data.frame(matrix(data = NA, nrow = 0, ncol = 6))
colnames(output_df) <- c('parameter', 'percentage_change', 'NRMSE_humans', 'Correlation_humans', 'NRMSE_primates', 'Corrleation_primates')

param_idx <- seq(2,21,1) # use 2 for baseline

# create a loop that goes through each parameter in param_idx and changes it by -20 or +20%, runs the model, and calculates the NRMSE and correlation for humans and primates 
for(i in param_idx){
  if(i == 2){
    percentage_change <- 1
  } else {
    percentage_change <- c(0.8, 1.2)
  }
  
  for(j in percentage_change){
    yfv_params_list_tmp <- yfv_params_list$full_model
    yfv_params_list_tmp[[i]] <- yfv_params_list_tmp[[i]] * j
    
    # Event adds 10 infected monkeys (Ic) at t=15 (same as main models)
    event_setting <- list(
      func = event_function_reduce_mosquitoes,
      time = 15
    )
    
    # Run the ODE
    result <- as.data.frame(
      ode(
        y = unlist(state_start_list),
        times = times,
        func = yfv_model,
        parms = unlist(yfv_params_list_tmp),
        events = event_setting
      )
    )
    
    # format output
    result$Date <- start_date + result$time
    # lead human infections by 30 days to account for delay between infection and detection
    result$Ih <- lead(result$Ih, n = 30)
    # join result and validation
    x <- result %>% left_join(val_data) %>% as.data.frame()
    
    # calculate correlation and NRMSE for humans and primates
    corHuman <- calc_correlation(df = x, rho = rho_humans, var1 = 'Ih', var2 = 'MG_human')
    corPrimates <- calc_correlation(df = x, rho = rho_monkeys, var1 = 'Ip', var2 = 'MG_primate')
    
    # add to output dataframe
    if(i == 2){
      names(yfv_params_list_tmp)[i] <- 'baseline'
    }
    output_df <- rbind(output_df, data.frame('parameter' = names(yfv_params_list_tmp)[i],
                                             'percentage_change' = j,
                                             'NRMSE_humans' = corHuman$nrmse,
                                             'Correlation_humans' = corHuman$corr,
                                             'NRMSE_primates' = corPrimates$nrmse,
                                             'Corrleation_primates' = corPrimates$corr))  
  }
}

# write output to csv
write.csv(output_df, '../sensitivity_analysis_output.csv', row.names = FALSE)


