# load library
# library(lubridate)

# Load parameter and initial condition data
p <- read.csv('parameter_values.csv')
IC <- read.csv('init_conditions.csv')

# Load seasonal forcing function
source('function_seasonal_forcing.R')

# Set up simulation time vector
yfv_epidemic <- seq.Date(
  from = as.Date(IC$nonnumeric_value[IC$variable == 'start_date'], '%m/%d/%Y'),
  to = as.Date(IC$nonnumeric_value[IC$variable == 'end_date'], '%m/%d/%Y'),
  by = 'day'
)
times <- seq_along(yfv_epidemic)

# Vaccination timeline and rate
vax_start <- as.Date(IC$nonnumeric_value[IC$variable == 'start_date_vaccination'], '%m/%d/%Y') + 30
pre_vax_days <- as.numeric(difftime(vax_start, yfv_epidemic[1]))
post_vax_days <- length(times) - pre_vax_days
vaccination_rate <- with(p, value[variable == 'final_v'] - value[variable == 'initial_v']) / post_vax_days
v_ts <- c(rep(0, pre_vax_days), rep(vaccination_rate, post_vax_days))

# Seasonal forcings for K, birth rates, biting rates
k <- seasonal_forcing(times, high = p$high_value[p$variable == 'K'], low = p$low_value[p$variable == 'K'])
br_aa <- seasonal_forcing(times, high = p$high_value[p$variable == 'mu_aa'], low = p$low_value[p$variable == 'mu_aa'])
br_hg <- seasonal_forcing(times, high = p$high_value[p$variable == 'mu_hg'], low = p$low_value[p$variable == 'mu_hg'])
a_vals <- seasonal_forcing(times, high = p$high_value[p$variable == 'a'], low = p$low_value[p$variable == 'a'])

# External forest seeding (immigration)
# for long term, perhaps just keep the values low (maybe between 0.001 and 0.2)
trend <- c(seq(1.25, 0.3, length = 300)
           , seq(0.3, 4, length = 300)
           , rep(0.001, length = length(times) - (300+300)))

# plot(yfv_epidemic, rising_trend, type = 'l')
# Optional: add a sharp seasonal pattern
base_season <- seasonal_forcing(times, high = 1, low = 0, phase = -pi/4)
sharpened_season <- base_season^7  # adjust exponent for sharper peaks

# Combine to get final monkey seeding
monkey_seeding_vals <-  p$value[p$variable == 'monkey_seed'] * sharpened_season * trend
marmoset_seeding_vals <- p$value[p$variable == 'marmoset_seed'] * sharpened_season * trend
hg_seeding_vals <- p$value[p$variable == 'mosquito_seed'] * sharpened_season * trend

# Create approx functions for time-varying parameters
a_approx <- approxfun(times, a_vals)
Vapprox <- approxfun(times, v_ts)
Kapprox <- approxfun(times, k)
br1approx <- approxfun(times, br_hg)
br2approx <- approxfun(times, br_aa)
monkey_seed_approx <- approxfun(times, monkey_seeding_vals)
marmoset_seed_approx <- approxfun(times, marmoset_seeding_vals)
hg_seed_approx <- approxfun(times, hg_seeding_vals)

# Parameter list constructor
get_full_params <- function(alpha = a_approx
                            , monkey_seed = monkey_seed_approx
                            , marmoset_seed = marmoset_seed_approx
                            , hg_seed = hg_seed_approx) {
  list(alpha = alpha
    , sigma = p$value[p$variable == 'sigma']
    , b = p$value[p$variable == 'b']
    , pMI1 = p$value[p$variable == 'pMI1']
    , pMI2 = p$value[p$variable == 'pMI2']
    , pMI3 = p$value[p$variable == 'pMI3']
    , pMI4 = p$value[p$variable == 'pMI4']
    , PDR_hg = p$value[p$variable == 'PDR_hg']
    , PDR_aa = p$value[p$variable == 'PDR_aa']
    , mu_hg = p$value[p$variable == 'mu_hg']
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
    , V = Vapprox
    , K = Kapprox
    , birth_hg = br1approx
    , birth_aa = br2approx
    , monkey_seed = monkey_seed
    , marmoset_seed = marmoset_seed
    , hg_seed = hg_seed
  )
}

# Model comparison ---------------------------------------
# Full model parameter list
yfv_params <- get_full_params()

# Scenarios with fixed parameters

# constant biting rate
a_fixed <- approxfun(times, rep(p$value[p$variable == 'a'], length(times)))
constant_alpha <- get_full_params(alpha = a_fixed)

# no movement
no_seeding <- approxfun(times, rep(0, length(times)))
no_movement <- get_full_params(monkey_seed = no_seeding, marmoset_seed = no_seeding, hg_seed = no_seeding)

# constant biting rate and no movement
constant_alpha_and_movement <- get_full_params(alpha = a_fixed, monkey_seed = no_seeding, marmoset_seed = no_seeding, hg_seed = no_seeding)

# Parameter sensitivity ---------------------------------
adjust_params <- function(base, updates) {
  for (varName in names(updates)) {
    base[[varName]] <- updates[[varName]]
  }
  return(base)
}

# low/high YFV mortality rates for primates
low_mu_v1 <- adjust_params(yfv_params, list(mu_v1 = p$low_value[p$variable == 'mu_v1']))
high_mu_v1 <- adjust_params(yfv_params, list(mu_v1 = p$high_value[p$variable == 'mu_v1']))

# high pMI values
pMI_update <- list(
  pMI1 = p$high_value[p$variable == 'pMI1'],
  pMI2 = p$high_value[p$variable == 'pMI2'],
  pMI3 = p$high_value[p$variable == 'pMI3'],
  pMI4 = p$high_value[p$variable == 'pMI4']
)

high_pMI <- adjust_params(yfv_params, pMI_update)

# low/high movement rates for forest species to city
# low
monkey_seed_low <- approxfun(times, p$low_value[p$variable == 'monkey_seed'] * sharpened_season * trend)
marmoset_seed_low <- approxfun(times, p$low_value[p$variable == 'marmoset_seed'] * sharpened_season * trend)
mosquito_seed_low <- approxfun(times, p$low_value[p$variable == 'mosquito_seed'] * sharpened_season * trend)

low_movement_update <- list(
  monkey_seed = monkey_seed_low,
  marmoset_seed = marmoset_seed_low,
  hg_seed = mosquito_seed_low
)

low_move <- adjust_params(yfv_params, low_movement_update)

# high
monkey_seed_high <- approxfun(times, p$high_value[p$variable == 'monkey_seed'] * sharpened_season * trend)
marmoset_seed_high <- approxfun(times, p$high_value[p$variable == 'marmoset_seed'] * sharpened_season * trend)
mosquito_seed_high <- approxfun(times, p$high_value[p$variable == 'mosquito_seed'] * sharpened_season * trend)

high_movement_update <- list(
  monkey_seed = monkey_seed_high,
  marmoset_seed = marmoset_seed_high,
  hg_seed = mosquito_seed_high
)

high_move <- adjust_params(yfv_params, high_movement_update)

# Interventions --------------------------------------------

yfv_params_interventions <- get_full_params()

# extend simulatation to 8 years
int_times <- seq(1, 8 * 365)

# get seasonal forcings for extended time vector
a_vals_long <- seasonal_forcing(int_times, p$high_value[p$variable == 'a'], p$low_value[p$variable == 'a'])
br_aa_long <- seasonal_forcing(int_times, p$high_value[p$variable == 'mu_aa'], p$low_value[p$variable == 'mu_aa'])
br_hg_long <- seasonal_forcing(int_times, p$high_value[p$variable == 'mu_hg'], p$low_value[p$variable == 'mu_hg'])
k_long <- seasonal_forcing(int_times, p$high_value[p$variable == 'K'], p$low_value[p$variable == 'K'])
v_ts_long <- c(v_ts, rep(0.0001, length(int_times) - length(v_ts)))

base_season_long <- seasonal_forcing(seq(length(times)+1, length(int_times), 1), high = 0.2, low = 0.001, phase = -pi) 
trend_long <- c(seq(0.01, 0.2, length = 10), rep(1, length(int_times) - length(times) - 10))
monkey_seeding_long <- c(monkey_seeding_vals, p$value[p$variable == 'monkey_seed']/3 * base_season_long * trend_long)
marmoset_seeding_long <- c(marmoset_seeding_vals, p$value[p$variable == 'marmoset_seed']/2 * base_season_long * trend_long)
hg_seeding_long <- c(hg_seeding_vals, p$value[p$variable == 'mosquito_seed']/2 * base_season_long * trend_long)

# Approxfuns for interventions
a_approx_long <- approxfun(int_times, a_vals_long)
br_hg_approx_long <- approxfun(int_times, br_hg_long)
br_aa_approx_long <- approxfun(int_times, br_aa_long)
V_approx_long <- approxfun(int_times, v_ts_long)
K_approx_long <- approxfun(int_times, k_long)
monkey_seed_approx_long <- approxfun(int_times, monkey_seeding_long)
marmoset_seed_approx_long <- approxfun(int_times, marmoset_seeding_long)
hg_seed_approx_long <- approxfun(int_times, hg_seeding_long)

# Update parameters
interventions_update <- list(
  alpha = a_approx_long,
  birth_hg = br_hg_approx_long,
  birth_aa = br_aa_approx_long,
  V = V_approx_long,
  K = K_approx_long,
  monkey_seed = monkey_seed_approx_long,
  marmoset_seed = marmoset_seed_approx_long,
  hg_seed = hg_seed_approx_long
)

yfv_params_interventions <- adjust_params(yfv_params_interventions, interventions_update)

# Reduce mosquitoes: set up event to emulate vector control efforts every rainy season
# Set a threshold (you can adjust slightly based on signal sharpness)
threshold <- quantile(br_aa_long, 0.9)  

# Identify times when rainy season starts
rainy_start_days <- which(diff(br_aa_long > threshold) == 1)

# Now map rainy seasons
rainy_windows <- lapply(rainy_start_days, function(day) {
  list(
    start_day = day,
    end_day = day + 30  # intervention lasts 30 days
  )
})

event_function_reduce_mosquitoes <- function(t, state, parameters) {
  state <- unlist(state)
  
  # Always apply import
  if (abs(t - 15) < 1e-6) {
    state['Ic'] <- state['Ic'] + 10
  }
  
  # Only apply mosquito control if rainy_windows is defined
  if (!is.null(parameters$rainy_windows)) {
    for (window in parameters$rainy_windows) {
      if (t >= window$start_day && t <= window$end_day) {
        fraction <- (t - window$start_day) / (window$end_day - window$start_day)
        smooth_ramp <- 0.5 * (1 - cos(pi * fraction))
        final_reduction <- 0.5
        daily_multiplier <- 1 - (1 - final_reduction) * smooth_ramp
        
        Saa0 <- state['Saa']
        Shg0 <- state['Shg']
        
        state['Saa'] <- Saa0 * daily_multiplier
        state['Eaa'] <- state['Eaa'] * daily_multiplier
        state['Iaa'] <- state['Iaa'] * daily_multiplier
        state['Shg'] <- Shg0 * daily_multiplier
        state['Ehg'] <- state['Ehg'] * daily_multiplier
        state['Ihg'] <- state['Ihg'] * daily_multiplier
      }
    }
  }
  
  return(state)
}

  
# event_function_reduce_mosquitoes <- function(t, state, parameters) {
#   for (window in rainy_windows) {
#     if (t >= window$start_day && t <= window$end_day) {
#       
#       # apply cosine smoothing reduction
#       fraction <- (t - window$start_day) / (window$end_day - window$start_day)
#       smooth_ramp <- 0.5 * (1 - cos(pi * fraction))  # 0 → 1 smoothly over window
#       
#       final_reduction <- 0.5  # Target 50% final reduction
#       daily_multiplier <- 1 - (1 - final_reduction) * smooth_ramp
#       
#       # Apply the gradual decrease
#       state['Saa'] <- state['Saa'] * daily_multiplier
#       state['Eaa'] <- state['Saa'] * daily_multiplier
#       state['Iaa'] <- state['Saa'] * daily_multiplier
#       state['Shg'] <- state['Shg'] * daily_multiplier
#       state['Ehg'] <- state['Shg'] * daily_multiplier
#       state['Ihg'] <- state['Shg'] * daily_multiplier
#     }
#   }
#   return(state)
# }

vector_control <- yfv_params_interventions

# early vaccination
# Define intervention date and ID
intervention_date <- as.Date('2016-12-15')
intervention_date_id <- which(yfv_epidemic == intervention_date)

vax_start_i <- intervention_date_id + 30
vax_early <- c(v_ts[vax_start_i:length(v_ts)], rep(1/365, length(int_times) - (length(v_ts)-vax_start_i) - 1))
vax_early <- approxfun(int_times, vax_early)

early_vax <- yfv_params_interventions
early_vax$V <- vax_early

# limit movement of howler monkeys
monkey_consevation <- monkey_seeding_long
monkey_consevation[(length(times) + 1):length(int_times)] <- monkey_consevation[(length(times) + 1):length(int_times)]/2
monkey_consevation <- approxfun(int_times, monkey_consevation)

limit_monkey_movement <- yfv_params_interventions
limit_monkey_movement$monkey_seed <- monkey_consevation

# combine interventions
interventions_combined <- early_vax
interventions_combined$monkey_seed <- monkey_consevation

# Interventions - increased R0 ----------------------------------
scale_R0 <- function(params) {
  params[c('pMI1', 'pMI2', 'pMI3', 'pMI4')] <- lapply(
    params[c('pMI1', 'pMI2', 'pMI3', 'pMI4')], function(x) x * 2)
  return(params)
}

vector_control_Inc_R0 <- scale_R0(vector_control)
early_vax_Inc_R0 <- scale_R0(early_vax)
limit_monkey_movement_Inc_R0 <- scale_R0(limit_monkey_movement)
combined_interventions_Inc_R0 <- scale_R0(interventions_combined)

# Create list of parameter sets
yfv_params_list <- list(
  # Model comparison
  yfv_params
  , constant_alpha
  , no_movement
  , constant_alpha_and_movement
  # Sensitivity analysis
  , low_mu_v1
  , high_mu_v1
  , high_pMI
  , low_move
  , high_move
  # Interventions
  , vector_control
  , early_vax
  , limit_monkey_movement
  , interventions_combined
  # Interventions - increased R0
  , vector_control_Inc_R0
  , early_vax_Inc_R0
  , limit_monkey_movement_Inc_R0
  , combined_interventions_Inc_R0
  )

names(yfv_params_list) <- c(
  # Model comparison
  'full_model'
  , 'constant_bite_rate'
  , 'no_movement'
  , 'constant_bite_rate_and_no_movement'
  # Sensitivity analysis
  , 'low_mu_v1'
  , 'high_mu_v1'
  , 'high_pMI'
  , 'low_movement'
  , 'high_movement'
  # Interventions
  , 'reduce_mosquitoes'
  , 'shift_vax'
  , 'limit_monkey_movement'
  , 'combined_interventions'
  # Interventions - increased R0
  , 'reduce_mosquitoes_high_R0'
  , 'shift_vax_high_R0'
  , 'limit_monkey_movement_high_R0'
  , 'combined_interventions_high_R0'
)

# Helper functions for simulation event setup
# Return a list of time vectors matching the length of parameter sets
get_times_list <- function(times_short, times_long, n_short = 9, n_long = 8) {
  c(rep(list(times_short), n_short), rep(list(times_long), n_long))
}

# Generate corresponding time lists
times_list <- get_times_list(times_short = times, times_long = int_times)

# Return indices of parameter list names with partial matching
get_scenario_indices <- function(scenario_names, match_terms) {
  matches <- sapply(scenario_names, function(name) {
    any(sapply(match_terms, function(term) grepl(term, name)))
  })
  which(matches)
}

# Get scenario indices where events should apply
specific_idx <- get_scenario_indices(names(yfv_params_list), match_terms = c('reduce', 'combined'))

# Get event times
event_times <- unlist(rainy_windows)




