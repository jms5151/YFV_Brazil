return_period <- function(values, return_period){
  # Fit a Normal distribution (change to Gamma or Weibull if needed)
  fit <- fitdist(values, 'norm')  # Adjust distribution if needed
  
  # Extract parameters
  params <- fit$estimate
  
  # Calculate the probability for the return period
  target_probability <- 1 / return_period  # Exceedance probability for rare drought
  
  # Use target_probability directly for droughts (low values)
  return_value <- qnorm(target_probability, mean = params["mean"], sd = params["sd"])  
  
  return(return_value)  
}
