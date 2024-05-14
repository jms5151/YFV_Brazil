# load validation data
realdata <- read.csv('../validation_data.csv')
realdata$Date <- as.Date(realdata$Month, '%m/%d/%Y')

# source model and data
source('model.R')
source('state_variables.R')
source('parameters.R')
source('biting_rate_drought_functions.R')

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

pdf('../Figures/model_v_data.pdf', height = 4, width = 10)
par(mfrow = c(1,2), mar = c(2.5,4,1,1))
plot(yfv_epidemic, out$I_h, type = 'l', ylab = 'Inf people')
lines(realdata$Date, realdata$MG_human, col = 'orange', type = 'b', pch = 16)
plot(yfv_epidemic, out$I_p, type = 'l', ylab = 'Inf primates', ylim =  c(0, 75))
lines(realdata$Date, realdata$MG_primate, col = 'orange', type = 'b', pch = 16)
dev.off()

# may want to sum to months and replot to match with observations
par(mfrow = c(2,2), mar = c(2.5,4,1,1))
plot(yfv_epidemic, out$I_p, type = 'l', ylab = 'Inf primates')
plot(yfv_epidemic, out$I_h, type = 'l', ylab = 'Inf people')
plot(yfv_epidemic, out$S_p, type = 'l', ylab = 'S primates')
plot(yfv_epidemic, out$R_p, type = 'l', ylab = 'R primates')



out$Date <- yfv_epidemic
out$YM <- format(out$Date, '%Y-%m')
library(tidyverse)
out2 <- out %>%
  group_by(YM) %>%
  summarise(InfPrimates = sum(I_p),
            InfHumans = sum(I_h))

plot.ts(out2$InfPrimates)
plot.ts(out2$InfHumans)

# plot(yfv_epidemic, out$R_p, type = 'l', ylab = 'R primates')
par(mfrow = c(3,3), mar = c(2.5,4,1,1))
plot(yfv_epidemic, out$I_p, type = 'l', ylab = 'Inf primates')
plot(yfv_epidemic, out$I_c, type = 'l', ylab = 'Inf marmosets')
plot(yfv_epidemic, out$I_h, type = 'l', ylab = 'Inf people')
plot(yfv_epidemic, out$S_p, type = 'l', ylab = 'Susceptible primates')
plot(yfv_epidemic, out$S_c, type = 'l', ylab = 'Susceptible marmosets')
plot(yfv_epidemic, out$S_h, type = 'l', ylab = 'Susceptible people')
plot(yfv_epidemic, out$I_hm, type = 'l', ylab = 'Inf Hm') # , ylim = c(0, max(out$I_hm))
plot(yfv_epidemic, out$S_hm, type = 'l', ylab = 'Susceptible Hm') # , ylim = c(0, max(out$I_hm))
plot(yfv_epidemic, out$I_aa, type = 'l', ylab = 'Inf aa') # , ylim = c(0, max(out$I_hm))

plot(yfv_epidemic, out$I_h, type = 'l', ylab = 'Inf people', ylim = c(0,2000))

plot(yfv_epidemic, out$a1, ylim = c(0,1), type = 'l', ylab = 'Hm biting rate on NHP')
plot(yfv_epidemic, out$a2, ylim = c(0,1), type = 'l', ylab = 'Hm biting rate on humans')
# plot(yfv_epidemic, out$a3, ylim = c(0,1), type = 'l')
plot(yfv_epidemic, out$S_hm, type = 'l', ylab = 'Susceptible Hm')
plot(yfv_epidemic, out$S_aa, type = 'l', ylab = 'Susceptible Aa')

par(mfrow = c(1,3), mar = c(2.5,4,1,1))
plot(yfv_epidemic, out$S_p, type = 'l', ylab = 'Susceptible primates')
plot(yfv_epidemic, out$S_c, type = 'l', ylab = 'Susceptible marmosets')
plot(yfv_epidemic, out$S_h, type = 'l', ylab = 'Susceptible people')
# outd <- out
# lines(yfv_epidemic, outd$I_p, col = 'orange')
# lines(yfv_epidemic, outd$I_c, col = 'orange')
# lines(yfv_epidemic, outd$I_h, col = 'orange')
# lines(yfv_epidemic, outd$I_hm, col = 'orange')
# lines(yfv_epidemic, outd$I_aa, col = 'orange')
# lines(yfv_epidemic, outd$S_c, col = 'orange')
# lines(yfv_epidemic, outd$S_p, col = 'orange')



lines(yfv_epidemic, out$I_p, col = 'orange')
lines(yfv_epidemic, out$I_c, type = 'l', ylab = 'Inf marmosets')
lines(yfv_epidemic, out$I_h, type = 'l', ylab = 'Inf people')
lines(yfv_epidemic, out$I_hm, type = 'l', ylab = 'Inf Hm') # , ylim = c(0, max(out$I_hm))
lines(yfv_epidemic, out$I_aa, type = 'l', ylab = 'Inf aa') # , ylim = c(0, max(out$I_hm))
lines(yfv_epidemic, out$S_c, type = 'l', ylab = 'Susceptible marmosets')
lines(yfv_epidemic, out$S_p, type = 'l', ylab = 'Susceptible primates')



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

par(mfrow = c(3,2), mar = c(2.5,4,1,1))
plot(yfv_epidemic, out$I_p, type = 'l', ylab = 'Inf primates')
plot(yfv_epidemic, out$I_c, type = 'l', ylab = 'Inf marmosets')
plot(yfv_epidemic, out$I_h, type = 'l', ylab = 'Inf people')
plot(yfv_epidemic, out$I_hm, type = 'l', ylab = 'Inf Hm') # , ylim = c(0, max(out$I_hm))
plot(yfv_epidemic, out$I_aa, type = 'l', ylab = 'Inf aa') # , ylim = c(0, max(out$I_hm))

plot(yfv_epidemic, out$S_p, type = 'l')
plot(yfv_epidemic, out$S_c, type = 'l')
plot(yfv_epidemic, out$S_hm, type = 'l')


# plot(out$E_hm, type = 'l')
# plot(out$S_hm, type = 'l')
# plot(out$S_aa, type = 'l')

plot(yfv_epidemic, out$S_p, type = 'l')

out <- out_no_spei_effect
par(mfrow = c(2,2), mar = c(2.5,4,1,1))
plot(yfv_epidemic, out$I_p, type = 'l', ylab = 'Inf primates')
lines(yfv_epidemic, out_no_spei_effect$I_p, col = 'orange')

plot(yfv_epidemic, out$I_h, type = 'l', ylab = 'Inf people')
lines(yfv_epidemic, out_no_spei_effect$I_h, col = 'orange')

plot(yfv_epidemic, out$I_hm, type = 'l', ylab = 'Inf Hm')
lines(yfv_epidemic, out_no_spei_effect$I_hm, col = 'orange') 

plot(yfv_epidemic, out$I_aa, type = 'l', ylab = 'Inf aa')
lines(yfv_epidemic, out_no_spei_effect$I_aa, col = 'orange') 

