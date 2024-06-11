# set up time sequence for simulation
start_date <- as.Date('2016-11-01') # may want to shift this to November
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
plot(yfv_epidemic, k, type = 'l', col = 'blue', xlab = 'Day', ylab = 'Value')

# birth rates as function of season
low_br_aa <- 1/35
high_br_aa <- 1/10
br_aa <- (high_br_aa + low_br_aa)/2 + (high_br_aa - low_br_aa)/2 * cos(days * frequency)
# plot(yfv_epidemic, br_aa, type = 'l', col = 'blue', xlab = 'Day', ylab = 'Value')

low_br_hm <- 1/27
high_br_hm <- 1/7
br_hm <- (high_br_hm + low_br_hm)/2 + (high_br_hm - low_br_hm)/2 * cos(days * frequency)
# plot(yfv_epidemic, br_hm, type = 'l', col = 'blue', xlab = 'Day', ylab = 'Value')

# biting rates
# source('biting_rate_drought_functions.R')
# source('format_drought_data.R')
# spei <- subset(spei, Date >= start_date & Date <= end_date)
# spei$Drought <- spei$SPEI.1
# a1vals <-  sapply(spei$Drought, function(x) a1(x))
# a2vals <- sapply(spei$Drought, function(x) a2(x))
# a3vals <- sapply(spei$Drought, function(x) a3(x))

# 0.35 for Hm on people, 0.7 for Hm on monkeys and Aa on people, but could increase to 0.86
low_bite_rate <- 0.5
high_bite_rate <- 1
a1vals <- (high_bite_rate + low_bite_rate)/2 + (high_bite_rate - low_bite_rate)/2 * sin(days * frequency - 60)
# plot.ts(a1vals)
# divide wave in half for Hm on people

# Monkeys move from "R" compartment towards city at rate x
# look at range of m = from <1-15%/month? 
moveRate <- 0.5/(12*30)
movement <- (moveRate + 0)/2 + (moveRate - 0)/2 * sin(days * frequency - 60)
plot(yfv_epidemic, movement, type = 'l')
# lines(yfv_epidemic, movement, col = 'blue')
# movement <- ifelse(x$SPEI.1 > 1, movement, 0)

# create time varying functions
a1approx <- approxfun(times, a1vals)
a2approx <- approxfun(times, a1vals/2)
a3approx <- approxfun(times, a1vals)
Vapprox <- approxfun(times, v_ts)
Kapprox <- approxfun(times, k)
br1approx <- approxfun(times, br_hm)
br2approx <- approxfun(times, br_aa)
mapprox <- approxfun(times, movement)

# list parameters
yfv_params <- list(
  a1 = a1approx
  , a2 = a2approx
  , a3 = a3approx
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
  , mu_v1 = 0.5#p$value[p$variable == 'mu_v1'] # use 50 and 80% (or 85%)
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

# remove biting rate as function of drought
a1approx_fixed <- approxfun(times, rep(p$value[p$variable=='a1'], length(times)))
a2approx_fixed <- approxfun(times, rep(p$value[p$variable=='a2'], length(times)))
a3approx_fixed <- approxfun(times, rep(p$value[p$variable=='a3'], length(times)))

yfv_params_bite_fixed <- yfv_params
yfv_params_bite_fixed$a1 <- a1approx_fixed
yfv_params_bite_fixed$a2 <- a2approx_fixed
yfv_params_bite_fixed$a3 <- a3approx_fixed

# remove movement as function of season
m_fixed <- approxfun(times, rep(0, length(times)))

yfv_params_move_fixed <- yfv_params
yfv_params_move_fixed$m <- m_fixed

# remove biting rate and movement as function of season
yfv_params_fixed <- yfv_params_bite_fixed
yfv_params_fixed$m <- m_fixed

