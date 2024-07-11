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

# Run the ODE model in parallel for each combination of state_start_list and yfv_params_list
resultsNew <- foreach(yfv_params_idx = 1:length(yfv_params_list), .packages = 'deSolve') %:%
  foreach(state_start_idx = 1:length(state_start_list), .packages = 'deSolve') %dopar% {
    state_start <- state_start_list[[state_start_idx]]
    yfv_params <- yfv_params_list[[yfv_params_idx]]
    result <- as.data.frame(
      ode(
        y = unlist(state_start),  # Ensure state_start is a numeric vector
        times = times,
        func = yfv_model,
        parms = yfv_params
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
saveRDS(object = resultsNew, file = '../model_results.RData')
