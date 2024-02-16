
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
sum(out$I_h)*p$value[p$variable=='mu_v2']
# primates infected
sum(out$I_p)

par(mfrow = c(2,2), mar = c(2.5,4,1,1))
plot(yfv_epidemic, out$I_p, type = 'l', ylab = 'Inf primates')
plot(yfv_epidemic, out$I_h, type = 'l', ylab = 'Inf people')


plot(yfv_epidemic, out$I_hm, type = 'l', ylab = 'Inf Hm', ylim = c(0, max(out$I_hm)))
plot(yfv_epidemic, out$I_aa, type = 'l', ylab = 'Inf aa', ylim = c(0, max(out$I_hm)))
plot(out$E_hm, type = 'l')
plot(out$S_hm, type = 'l')
plot(out$S_aa, type = 'l')

plot(yfv_epidemic, out$S_p, type = 'l')

