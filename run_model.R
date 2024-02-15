# set up time sequence for simulation
years <- 3
times <- seq(from = 1, to = 360 * years, by = 1)

# source model and data
source('model.R')
source('state_variables.R')
source('parameters.R')

# run model
out <- as.data.frame(
  ode(
    state_start
    , times
    , yfv_model
    , yfv_params
    # , rtol = 1e-12
    # , hmax = 1 / 120
  )
)

# approximately 2200 confirmed cases
# 770 confirmed deaths

# total cases
sum(out$I_h)
# reporting rate
report = 0.10
# symptomatic cases
sum(out$I_h)*report
# deaths
sum(out$I_h)*p$value[p$variable=='V']
# primates infected
sum(out$I_p)


par(mfrow = c(2,1), mar = c(2.5,1,1,1))
plot(out$I_p, type = 'l')
plot(out$I_h, type = 'l')


plot(out$I_hm, type = 'l')
plot(out$I_aa, type = 'l')
plot(out$E_hm, type = 'l')
plot(out$S_hm, type = 'l')
plot(out$S_aa, type = 'l')
