# source parameter values
p <- read.csv('parameter_values.csv')
IC <- read.csv('init_conditions.csv')

# source seasonal forcing functions
source('function_seasonal_forcing.R')

# set up time sequence for simulation
start_date <- as.Date(IC$nonnumeric_value[IC$variable == 'start_date'], '%m/%d/%Y')
end_date <- as.Date(IC$nonnumeric_value[IC$variable == 'end_date'], '%m/%d/%Y')
yfv_epidemic <- seq.Date(start_date, end_date, by = 'days')
times <- seq(from = 1, to = length(yfv_epidemic), by = 1)

# vaccination rate
vax_campaign_start <- as.Date(IC$nonnumeric_value[IC$variable == 'start_date_vaccination'], '%m/%d/%Y') + 30
prevaxtimes <- as.numeric(difftime(vax_campaign_start, start_date))
postvac_times <- length(times) - prevaxtimes

start_pop_vaccinated <- p$value[p$variable == 'initial_v']
end_pop_vaccinated <- p$value[p$variable == 'final_v']

vaccination_rate <- (end_pop_vaccinated - start_pop_vaccinated)/postvac_times

v_ts <- c(rep(0, prevaxtimes), rep(vaccination_rate, postvac_times))

# carrying capacity as function of season
k <- seasonal_forcing(times = times, high = p$value[p$variable == 'K_wet'], low = p$value[p$variable == 'K_dry'])
# Plot the wave
# plot(yfv_epidemic, k, type = 'l', col = 'blue', xlab = 'Day', ylab = 'Value')

# birth rates as function of season
br_aa <- seasonal_forcing(times = times, high = p$value[p$variable == 'mu_aa_high'], low = p$value[p$variable == 'mu_aa_low'])
# plot(yfv_epidemic, br_aa, type = 'l', col = 'blue', xlab = 'Day', ylab = 'Value')

br_hm <- seasonal_forcing(times = times, high = p$value[p$variable == 'mu_hm_high'], low = p$value[p$variable == 'mu_hm_low'])
# plot(yfv_epidemic, br_hm, type = 'l', col = 'blue', xlab = 'Day', ylab = 'Value')

# biting rates as a function of season
a1vals <- seasonal_forcing(times = times, high = p$value[p$variable == 'a1_high'], low = p$value[p$variable == 'a1_low'])
# plot.ts(a1vals)
# divide wave in half for Hm on people

# Monkeys move from "R" compartment towards city at rate x
# look at range of m = from <1-15%/month? 
movement_fn <- function(x1){
  moveRate <- x1/(12*30)
  movement <- seasonal_forcing2(times = times, x = moveRate) 
  return(movement)  
}

movement <- movement_fn(x1 = 0.5)

# plot(yfv_epidemic, movement, type = 'l')
# lines(yfv_epidemic, movement, col = 'blue')
# movement <- ifelse(x$SPEI.1 > 1, movement, 0)

# create time varying functions
a1approx <- approxfun(times, a1vals)
a2approx <- approxfun(times, a1vals/2)
Vapprox <- approxfun(times, v_ts)
Kapprox <- approxfun(times, k)
br1approx <- approxfun(times, br_hm)
br2approx <- approxfun(times, br_aa)
mapprox <- approxfun(times, movement)

# list parameters
yfv_params <- list(
  a1 = a1approx
  , a2 = a2approx
  , b = p$value[p$variable == 'b']
  , pMI1 = p$value[p$variable == 'pMI1']
  , pMI2 = p$value[p$variable == 'pMI2']
  , pMI3 = p$value[p$variable == 'pMI3']
  , pMI4 = p$value[p$variable == 'pMI4']
  , PDR_hm = p$value[p$variable == 'PDR_hm']
  , PDR_aa = p$value[p$variable == 'PDR_aa']
  , mu_hm = p$value[p$variable == 'mu_hm']
  , mu_aa = p$value[p$variable == 'mu_aa']
  , mu_p = p$value[p$variable == 'mu_p']
  , mu_h = p$value[p$variable == 'mu_h']
  , mu_c = p$value[p$variable == 'mu_c']
  , mu_v1 = p$value[p$variable == 'mu_v1'] 
  , mu_v2 = p$value[p$variable == 'mu_v2']
  , mu_v3 = p$value[p$variable == 'mu_v3']
  , gamma_p = p$value[p$variable == 'gamma_p']
  , gamma_h = p$value[p$variable == 'gamma_h']
  , gamma_c = p$value[p$variable == 'gamma_c']
  , delta_h = p$value[p$variable == 'delta_h']
  , p = p$value[p$variable == 'p']
  , V = Vapprox
  , K = Kapprox
  , br1 = br1approx # Hm
  , br2 = br2approx # Aa
  , m = mapprox
)

adjust_params <- function(varName, varValue){
  yfv_params[varName] <- varValue
  return(yfv_params)
}

# remove seasonally varying biting rates
yfv_params_bite_fixed <- yfv_params
a1approx_fixed <- approxfun(times, rep(p$value[p$variable=='a1'], length(times)))
a2approx_fixed <- approxfun(times, rep(p$value[p$variable=='a2'], length(times)))
yfv_params_bite_fixed$a1 <- a1approx_fixed
yfv_params_bite_fixed$a2 <- a2approx_fixed

# remove seasonally varying movement
yfv_params_move_fixed <- yfv_params
m_fixed <- approxfun(times, rep(0, length(times)))
yfv_params_move_fixed$m <- m_fixed

# remove seasonally varying biting rate and movement
yfv_params_fixed <- yfv_params_bite_fixed
yfv_params_fixed$m <- m_fixed

# adjust YFV-induced monkey mortality rates
yfv_params_low_mu_v1 <- adjust_params(varName = 'mu_v1', varValue = 0.2)
yfv_params_high_mu_v1 <- adjust_params(varName = 'mu_v1', varValue = 0.8)

# adjust proportion Hm bites on NHP
yfv_params_low_p <- adjust_params(varName = 'p', varValue = 0.3)
yfv_params_mod_p <- adjust_params(varName = 'p', varValue = 0.5)

# adjust monkey movement rates
yfv_params_mod_move <- yfv_params
movement_mod <- movement_fn(x1 = 2)
mapprox_mod <- approxfun(times, movement_mod)
yfv_params_mod_move$m <- mapprox_mod

yfv_params_high_move <- yfv_params
movement_high <- movement_fn(x1 = 10)
mapprox_high <- approxfun(times, movement_high)
yfv_params_high_move$m <- mapprox_high

# interventions
intervention_date <- as.Date('2017-03-01')
intervention_date_id <- which(yfv_epidemic == intervention_date)

# reduce mosquitoes
int_params_reduce_mosq <- yfv_params
k_new <- c(k[1:(intervention_date_id-1)], k[intervention_date_id:length(k)]/2)
k_new <- approxfun(times, k_new)

int_params_reduce_mosq$K <- k_new

# reduce NHP movement
int_params_reduce_nhp_movement <- yfv_params
move_new <- c(movement[1:(intervention_date_id-1)], movement[intervention_date_id:length(movement)]/2)
move_new <- approxfun(times, move_new)
int_params_reduce_nhp_movement$m <- move_new

# increase vaccination: do we want to increase final prop vaccinated or start time of vaccination?
int_params_vax <- yfv_params
prevaxtimes_i <- as.numeric(difftime(intervention_date, start_date))
postvac_times_i <- length(times) - prevaxtimes
vaccination_rate_i <- (end_pop_vaccinated - start_pop_vaccinated)/postvac_times_i
vax_new <- c(rep(0, prevaxtimes), rep(vaccination_rate_i, postvac_times_i))

# create list of parameters lists
yfv_params_list <- list(
  yfv_params
  , yfv_params_bite_fixed
  , yfv_params_move_fixed
  , yfv_params_fixed
  , yfv_params_low_mu_v1
  , yfv_params_high_mu_v1
  , yfv_params_low_p
  , yfv_params_mod_p
  , yfv_params_mod_move
  , yfv_params_high_move
  , int_params_reduce_mosq
  , int_params_reduce_nhp_movement
  , int_params_vax
)

names(yfv_params_list) <- c(
  'base_model'
  , 'fixed_bite_rate'
  , 'fixed_movement'
  , 'fixed'
  , 'low_mu_v1'
  , 'high_mu_v1'
  , 'low_p'
  , 'mod_p'
  , 'mod_move'
  , 'high_move'
  , 'reduce_mosquitoes'
  , 'reduce_NHP_movement'
  , 'increase_vax'
)
