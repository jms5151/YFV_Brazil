# population sizes
# primate_pop = S_p + I_p + R_p
# human_pop = S_h + E_h + I_h + R_h
# haem_pop = S_hm + E_hm + I_hm
# aedes_pop = S_aa + E_aa + I_aa

# set up time sequence for simulation
start_date <- as.Date('2016-12-15')
end_date <- as.Date('2018-12-15')
yfv_epidemic <- seq.Date(start_date, end_date, by = 'days')
times <- seq(from = 1, to = length(yfv_epidemic), by = 1)

# source other parameter values
p <- read.csv('parameter_values.csv')

# vaccination rate
start_vaccinated <- p$value[p$variable == 'initial_v']
final_vaccinated <- p$value[p$variable == 'final_v']
vaccination_rate <- (final_vaccinated - start_vaccinated)/length(times)

# carrying capacity
K_dry <- mosquitoes
# (Dec - April)
K_wet <- K_dry * 6

k <- ifelse(as.numeric(format(yfv_epidemic, '%m')) > 5 & as.numeric(format(yfv_epidemic, '%m')) < 11, K_dry, K_wet)


# Define parameters
# x <- 10  # Maximum value
# y <- 2   # Minimum value
# n_points <- 100  # Number of points in the wave
# frequency <- 2*pi/12  # Frequency for seasonal cycle (12 months in a year)
# 
# # Generate cosine wave
# time <- seq(0, 2*pi, length.out = length(yfv_epidemic))  # Time points
# cos_wave <- (K_wet + K_dry)/2 + (K_wet - K_dry)/2 * cos(time * frequency)
# 
# # Plot the wave
# plot(cos_wave, type = 'l', col = 'blue', xlab = 'Time', ylab = 'Value', main = 'Seasonal Cosine Wave')


# list parameters
yfv_params <- list(
  a1 = rnorm(n = length(times), mean = 0.5, sd = 0.4)#c(0.7, length(times))
  # , a1 = p$value[p$variable == 'a1']
  , a2 = rnorm(n = length(times), mean = 0.4, sd = 0.2)#c(0.35, length(times))
  # , a2 = p$value[p$variable == 'a2']
  , a3 = rnorm(n = length(times), mean = 0.4, sd = 0.2)#c(0.35, length(times))
  # , a3 = p$value[p$variable == 'a3']
  , b = p$value[p$variable == 'b']
  , pMI1 = p$value[p$variable == 'pMI1']
  , pMI2 = p$value[p$variable == 'pMI2']
  , pMI3 = p$value[p$variable == 'pMI3']
  , PDR_hm = p$value[p$variable == 'PDR_hm']
  , PDR_aa = p$value[p$variable == 'PDR_aa']
  , mu_hm = p$value[p$variable == 'mu_hm']
  , mu_aa = p$value[p$variable == 'mu_aa']
  , mu_p = p$value[p$variable == 'mu_p']
  , mu_h = p$value[p$variable == 'mu_h']
  , mu_v1 = p$value[p$variable == 'mu_v1']
  , mu_v2 = p$value[p$variable == 'mu_v2']
  , gamma_p = p$value[p$variable == 'gamma_p']
  , gamma_h = p$value[p$variable == 'gamma_h']
  , delta_h = p$value[p$variable == 'delta_h']
  , p = p$value[p$variable == 'p']
  , V = vaccination_rate#p$value[p$variable == 'V']
  , K = k
)
