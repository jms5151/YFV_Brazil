# load the required packages
library(foreach)
library(doParallel)
library(deSolve)

# model
source('model.R')

# lists of parameters
source('parameters.R')

# lists of initial conditions/state variables
source('state_variables.R')

# Set up parallel backend to use multiple processors
numCores <- detectCores() - 1
cl <- makeCluster(numCores)
registerDoParallel(cl)

#length(yfv_params_list)
# Run the ODE model in parallel for each combination of state_start_list and yfv_params_list
resultsNew <- foreach(yfv_params_idx = 1:5, .packages = 'deSolve') %:%
  foreach(state_start_idx = 1:length(state_start_list), .packages = 'deSolve') %dopar% {
    state_start <- state_start_list[[state_start_idx]]
    yfv_params <- yfv_params_list[[yfv_params_idx]]
    times <- times_list[[yfv_params_idx]]
    
    # event_setting <- list(func = event_function, time = c(400))
      
        # # Determine if the event function should be applied
    event_setting <- if (yfv_params_idx %in% specific_idx) {
      list(func = event_function, time = event_times)
    } else {
      NULL
    }

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
# saveRDS(object = resultsNew, file = '../model_results.RData')


out <- ode(
  y = state_start,
  times = times,
  func = yfv_model,
  parms = yfv_params
)

out_df <- as.data.frame(out)

out_df$Inhp <- rowSums(out_df[, c("Ip", "Ic")])
x <- out_df %>%
  mutate(Date = yfv_epidemic) %>%
  left_join(val_data, by = "Date")

plot(x$Date, x$Inhp * rho_monkeys, type = "l", lty = 1, ylab = "Infected", xlab = "Time (days)")
points(x$Date, x$MG_primate, col = "black", pch = 16)
legend("topright", legend = "NHP", lty = 1)

plot(x$Date, x$Ip * rho_monkeys, type = "l", lty = 1, ylab = "Infected", xlab = "Time (days)")
points(x$Date, x$MG_primate, col = "black", pch = 16)
legend("topright", legend = "Howler monkeys", lty = 1)

plot(x$Date, x$Ih * rho_humans, type = "l", lty = 1, ylab = "Infected", xlab = "Time (days)", ylim = c(0,400))
points(x$Date, x$MG_human, col = "black", pch = 16)
legend("topright", legend = "Humans", lty = 1)


plot(yfv_epidemic, resultsNew[[1]][[1]]$result$Ih, main = 'humans', type = 'l') # 'base_model'
lines(yfv_epidemic, resultsNew[[2]][[1]]$result$Ih, col = 'blue') # fixed_bite_rate
lines(yfv_epidemic, resultsNew[[3]][[1]]$result$Ih, col = 'red') # fixed_biting_weights
lines(yfv_epidemic, resultsNew[[4]][[1]]$result$Ih, col = 'purple') # fixed

plot(yfv_epidemic, resultsNew[[1]][[1]]$result$Ip, main = 'primates', type = 'l') # 'base_model'
lines(yfv_epidemic, resultsNew[[2]][[1]]$result$Ip, col = 'blue')
lines(yfv_epidemic, resultsNew[[3]][[1]]$result$Ip, col = 'red')
lines(yfv_epidemic, resultsNew[[4]][[1]]$result$Ip, col = 'purple')


plot.ts(resultsNew[[1]][[1]]$result$Sp)
plot.ts(resultsNew[[1]][[1]]$result$Sc)

plot.ts(resultsNew[[1]][[1]]$result$Ip)
plot.ts(resultsNew[[1]][[1]]$result$Ic)
plot.ts(resultsNew[[1]][[1]]$result$Ih)
plot.ts(resultsNew[[1]][[1]]$result$Ihg)
plot.ts(resultsNew[[1]][[1]]$result$Iaa)

plot.ts(resultsNew[[1]][[1]]$result$Saa)
plot.ts(resultsNew[[1]][[1]]$result$Nhg)
plot.ts(resultsNew[[1]][[1]]$result$Naa)
