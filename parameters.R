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
k <- seasonal_forcing(times, high = 7e6, low = 4e5)
br_aa <- seasonal_forcing(times, high = p$value[p$variable == 'mu_aa_high'], low = p$value[p$variable == 'mu_aa_low'])
br_hg <- seasonal_forcing(times, high = p$value[p$variable == 'mu_hg_high'], low = p$value[p$variable == 'mu_hg_low'])
a_vals <- seasonal_forcing(times, high = p$value[p$variable == 'a_high'], low = p$value[p$variable == 'a_low'])

# External forest seeding (immigration)
trend <- c(seq(0.6, 0.3, length = 300)
           , seq(0.3, 1, length = 300)
           , rep(0.001, length = length(times) - (300+300)))
# plot(yfv_epidemic, rising_trend, type = 'l')
# Optional: add a sharp seasonal pattern
base_season <- seasonal_forcing(times, high = 1, low = 0, phase = -pi/6)
sharpened_season <- base_season^8  # adjust exponent for sharper peaks

# Combine to get final monkey seeding
monkey_seeding_vals <- 50 * sharpened_season * trend

# monkey_seeding_vals <- seasonal_forcing(times, high = 300, low = 1, phase = 5) * declining_trend
plot(yfv_epidemic, monkey_seeding_vals, type = 'l')

marmoset_seeding_vals <- 75 * sharpened_season * trend#seasonal_forcing(times, high = 5000, low = 1, phase = 5) * declining_trend
hg_seeding_vals <- 10 * sharpened_season * trend#seasonal_forcing(times, high = 5000, low = 1, phase = 5) * declining_trend

# Create approx functions for time-varying parameters
a_approx <- approxfun(times, a_vals)
Vapprox <- approxfun(times, v_ts)
Kapprox <- approxfun(times, k)
br1approx <- approxfun(times, br_hg)
br2approx <- approxfun(times, br_aa)
monkey_seed_approx <- approxfun(times, monkey_seeding_vals)
marmoset_seed_approx <- approxfun(times, marmoset_seeding_vals)
marmoset_seed_approx <- approxfun(times, marmoset_seeding_vals)
hg_seed_approx <- approxfun(times, hg_seeding_vals)

# Parameter list constructor
get_base_params <- function(alpha = a_approx) {
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
    , monkey_seed = monkey_seed_approx
    , marmoset_seed = marmoset_seed_approx
    , hg_seed = hg_seed_approx
  )
}

# Base model parameter list
yfv_params <- get_base_params()

# yfv_params_long <- get_base_params()
# 
# # Scenario variants derived from base params
# # constant biting rate
# a_fixed <- approxfun(times, rep(mean(a_vals), length(times)))
# params_fixed_alpha <- get_base_params(alpha = a_fixed)
# 
# # constant host weights
# wp_fixed <- approxfun(times, rep(median(wp_vals), length(times)))
# wh_fixed <- approxfun(times, rep(median(wh_vals), length(times)))
# wc_fixed <- approxfun(times, rep(median(wc_vals), length(times)))
# # wp_fixed <- approxfun(times, rep(1.2, length(times)))
# # wh_fixed <- approxfun(times, rep(0.6, length(times)))
# # wc_fixed <- approxfun(times, rep(0.2, length(times)))
# params_fixed_weights <- get_base_params(wp = wp_fixed, wh = wh_fixed, wc = wc_fixed)
# 
# # constant biting rates and host weights
# yfv_params_fixed <- get_base_params(alpha = a_fixed, wp = wp_fixed, wh = wh_fixed, wc = wc_fixed)
# 
# adjust_params <- function(base, varName, varValue) {
#   base[[varName]] <- varValue
#   return(base)
# }
# 
# yfv_params_low_mu_v1 <- adjust_params(yfv_params, 'mu_v1', p$value[p$variable == 'mu_v1_low'])
# yfv_params_high_mu_v1 <- adjust_params(yfv_params, 'mu_v1', p$value[p$variable == 'mu_v1_high'])
# 
# # Time-varying components for long simulations
# int_times <- seq(1, 8 * 365)
# k_long <- seasonal_forcing(int_times, p$value[p$variable == 'K_wet'], p$value[p$variable == 'K_dry'])
# br_aa_long <- seasonal_forcing(int_times, p$value[p$variable == 'mu_aa_high'], p$value[p$variable == 'mu_aa_low'])
# br_hg_long <- seasonal_forcing(int_times, p$value[p$variable == 'mu_hg_high'], p$value[p$variable == 'mu_hg_low'])
# a_vals_long <- seasonal_forcing(int_times, p$value[p$variable == 'a_high'], p$value[p$variable == 'a_low'])
# v_ts_long <- c(v_ts, rep(0, length(int_times) - length(v_ts)))
# wp_vals_long <- seasonal_forcing(int_times, high = 1, low = 0.2)  # howler monkeys
# wh_vals_long <- rep(1, length(int_times))                         # humans = always available
# wc_vals_long <- seasonal_forcing(int_times, high = 0.8, low = 0.4)  # marmosets
# 
# quant50 <- unname(quantile(k_long, 0.5))
# 
# # Approxfuns for extended period
# a_approx_long <- approxfun(int_times, a_vals_long)
# V_approx_long <- approxfun(int_times, v_ts_long)
# K_approx_long <- approxfun(int_times, k_long)
# br_hg_approx_long <- approxfun(int_times, br_hg_long)
# br_aa_approx_long <- approxfun(int_times, br_aa_long)
# wp_approx_long <- approxfun(int_times, wp_vals_long)
# wh_approx_long <- approxfun(int_times, wh_vals_long)
# wc_approx_long <- approxfun(int_times, wc_vals_long)
# 
# yfv_params_long$alpha <- a_approx_long
# yfv_params_long$V <- V_approx_long
# yfv_params_long$K <- K_approx_long
# yfv_params_long$birth_hg <- br_hg_approx_long
# yfv_params_long$birth_aa <- br_aa_approx_long
# 
# 
# # Helper functions for simulation event setup
# # Return a list of time vectors matching the length of parameter sets
# get_times_list <- function(times_short, times_long, n_short = 10, n_long = 8) {
#   c(rep(list(times_short), n_short), rep(list(times_long), n_long))
# }
# 
# # Identify event time indices where K is above a quantile threshold, starting after a date
# get_event_times <- function(K_ts, intervention_idx, threshold = 0.5) {
#   above_thresh <- which(K_ts > quantile(K_ts, threshold))
#   above_thresh[above_thresh >= intervention_idx]
# }
# 
# # Return indices of parameter list names matching specified scenario terms
# get_scenario_indices <- function(scenario_names, match_terms) {
#   which(scenario_names %in% match_terms)
# }
# 
# # Generate intervention parameter sets and associated time vectors
# # Requires yfv_params, yfv_params_long, times, int_times, etc. to be defined
# 
# # Define intervention date and ID
# intervention_date <- as.Date('2016-12-15')
# intervention_date_id <- which(yfv_epidemic == intervention_date)
# 
# # Define modified movement and vaccination schedules
# # move_new <- c(movement_long[1:(intervention_date_id - 1)], movement_long[intervention_date_id:length(movement_long)] / 4)
# # move_new <- approxfun(int_times, move_new)
# 
# vax_start_i <- intervention_date_id + 30
# vax_early <- c(rep(0, vax_start_i), rep(vaccination_rate, length(v_ts_long) - vax_start_i))
# vax_early <- approxfun(int_times, vax_early)
# 
# # Define parameter variants
# # int_params_reduce_nhp_movement <- yfv_params_long
# # int_params_reduce_nhp_movement$m <- move_new
# 
# int_params_vax <- yfv_params_long
# int_params_vax$V <- vax_early
# 
# int_params_combined <- int_params_vax
# # int_params_combined$m <- move_new
# 
# # Increased R0 scenarios
# scale_R0 <- function(params) {
#   params[c('pMI1', 'pMI2', 'pMI3', 'pMI4')] <- lapply(
#     params[c('pMI1', 'pMI2', 'pMI3', 'pMI4')], function(x) x * 2)
#   return(params)
# }
# 
# yfv_params_long_Inc_R0 <- scale_R0(yfv_params_long)
# # int_params_reduce_nhp_movement_Inc_R0 <- scale_R0(int_params_reduce_nhp_movement)
# int_params_vax_Inc_R0 <- scale_R0(int_params_vax)
# int_params_combined_Inc_R0 <- scale_R0(int_params_combined)
# 
# # Create list of parameter sets
# yfv_params_list <- list(
#   yfv_params,
#   params_fixed_alpha,
#   params_fixed_weights,
#   yfv_params_fixed,
#   yfv_params_low_mu_v1,
#   yfv_params_high_mu_v1,
#   # yfv_params_low_p,
#   # yfv_params_mod_p,
#   # yfv_params_mod_move,
#   # yfv_params_high_move,
#   yfv_params_long,
#   # int_params_reduce_nhp_movement,
#   int_params_vax,
#   int_params_combined,
#   yfv_params_long_Inc_R0,
#   # int_params_reduce_nhp_movement_Inc_R0,
#   int_params_vax_Inc_R0,
#   int_params_combined_Inc_R0
# )
# 
# names(yfv_params_list) <- c(
#   'base_model',
#   'fixed_bite_rate',
#   'fixed_biting_weights',
#   'fixed',
#   'low_mu_v1',
#   'high_mu_v1',
#   # 'low_p',
#   # 'mod_p',
#   # 'mod_move',
#   # 'high_move',
#   'reduce_mosquitoes',
#   # 'reduce_NHP_movement',
#   'shift_vax',
#   'combined_interventions',
#   'reduce_mosquitoes_high_R0',
#   # 'reduce_NHP_movement_high_R0',
#   'shift_vax_high_R0',
#   'combined_interventions_high_R0'
# )
# 
# yfv_params_list[["reduce_mosquitoes"]]$quant50 <- quant50
# yfv_params_list[["combined_interventions"]]$quant50 <- quant50
# yfv_params_list[["reduce_mosquitoes_high_R0"]]$quant50 <- quant50
# yfv_params_list[["combined_interventions_high_R0"]]$quant50 <- quant50
# 
# # Generate corresponding time lists
# times_list <- get_times_list(times_short = times, times_long = int_times, n_short = 10, n_long = 8)
# 
# # Get event times
# event_times <- get_event_times(K_ts = k_long, intervention_idx = intervention_date_id)
# 
# # Get scenario indices where events should apply
# specific_idx <- get_scenario_indices(names(yfv_params_list),
#                                      match_terms = c('reduce_mosquitoes', 'combined_interventions'))
# 
