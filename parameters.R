# population sizes
# primate_pop = S_p + I_p + R_p
# human_pop = S_h + E_h + I_h + R_h
# haem_pop = S_hm + E_hm + I_hm
# aedes_pop = S_aa + E_aa + I_aa

# source other parameter values
p <- read.csv('parameter_values.csv')

start_vaccinated <- 0.5
final_vaccinated <- 0.9
vaccination_rate <- (final_vaccinated - start_vaccinated)/length(times)


# list parameters
yfv_params <- list(
  # N_p = 10000
  # , N_h = 10000
  # , N_hm = 25000
  # , N_aa = 25000
   sigma_hm = p$value[p$variable == 'sigma_hm']
  , sigma_aa = p$value[p$variable == 'sigma_aa']
  , a1 = rnorm(n = length(times), mean = 0.5, sd = 0.4)#c(0.6, length(times))
  # , a1 = p$value[p$variable == 'a1']
  , a2 = rnorm(n = length(times), mean = 0.4, sd = 0.2)#c(0.6, length(times))
  # , a2 = p$value[p$variable == 'a2']
  , a3 = p$value[p$variable == 'a3']
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
  , K = mosquitoes
)
