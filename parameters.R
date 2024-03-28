# population sizes
# primate_pop = S_p + I_p + R_p
# human_pop = S_h + E_h + I_h + R_h
# haem_pop = S_hm + E_hm + I_hm
# aedes_pop = S_aa + E_aa + I_aa

# set up time sequence for simulation
start_date <- as.Date('2016-12-15') # may want to shift this to November
end_date <- as.Date('2018-12-15')
yfv_epidemic <- seq.Date(start_date, end_date, by = 'days')
times <- seq(from = 1, to = length(yfv_epidemic), by = 1)

# source other parameter values
p <- read.csv('parameter_values.csv')

# vaccination rate
vax_campaign_start <- as.Date('2018-02-25') + 30
prevaxtimes <- as.numeric(difftime(vax_campaign_start, start_date))
postvac_times <- length(times) - prevaxtimes

start_pop_vaccinated <- p$value[p$variable == 'initial_v']
end_pop_vaccinated <- p$value[p$variable == 'final_v']

vaccination_rate <- (end_pop_vaccinated - start_pop_vaccinated)/postvac_times

v_ts <- c(rep(0, prevaxtimes), rep(vaccination_rate, postvac_times))

# carrying capacity, model using cosine wave
yrs <- round(length(times)/365)
frequency <- 2*yrs*pi/length(times)
days <- seq(1, length(times), by = 1)
K_dry <- mosquitoes
K_wet <- K_dry * 6 # (Dec - April)
k <- (K_wet + K_dry)/2 + (K_wet - K_dry)/2 * cos(days * frequency)
# Plot the wave
# plot(yfv_epidemic, k, type = 'l', col = 'blue', xlab = 'Day', ylab = 'Value')

# birth rates as function of season
low_br_aa <- 1/35
high_br_aa <- 1/10
br_aa <- (high_br_aa + low_br_aa)/2 + (high_br_aa - low_br_aa)/2 * cos(days * frequency)
# plot(yfv_epidemic, br_aa, type = 'l', col = 'blue', xlab = 'Day', ylab = 'Value')

low_br_hm <- 1/27
high_br_hm <- 1/7
br_hm <- (high_br_hm + low_br_hm)/2 + (high_br_hm - low_br_hm)/2 * cos(days * frequency)
# plot(yfv_epidemic, br_hm, type = 'l', col = 'blue', xlab = 'Day', ylab = 'Value')

# # drought index
# source('format_drought_data.R')
# spei <- subset(spei, Date >= start_date & Date <= end_date) 
# 
# # biting rates
# source('biting_rate_drought_functions.R')
# p <- read.csv('parameter_values.csv')

# list parameters
yfv_params <- list(
  # a1 = sapply(spei$Drought, function(x) a1(x))
  a1 = rep(p$value[p$variable == 'a1'], length(times))
  # , a2 = sapply(spei$Drought, function(x) a2(x))
  , a2 = rep(p$value[p$variable == 'a2'], length(times))
  # , a3 = sapply(spei$Drought, function(x) a3(x))
  , a3 = rep(p$value[p$variable == 'a3'], length(times))
  , b = p$value[p$variable == 'b']
  , pMI1 = p$value[p$variable == 'pMI1']
  , pMI2 = p$value[p$variable == 'pMI2']
  , pMI3 = p$value[p$variable == 'pMI3']
  , pMI4 = p$value[p$variable == 'pMI4']
  , PDR_hm = p$value[p$variable == 'PDR_hm']
  , PDR_aa = p$value[p$variable == 'PDR_aa']
  , mu_hm = p$value[p$variable == 'mu_hm']
  , mu_aa = p$value[p$variable == 'mu_aa']
  , mu_p = p$value[p$variable == 'mu_p'] # multiplying by 5 or 10 helps
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
  , V = v_ts
  , K = k
  , br1 = br_hm
  , br2 = br_aa
  # , w = 1/(10*365) # waning immunity
  , m = br_hm/365 #spei2$Drought/365 # immigration
  )
