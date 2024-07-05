# load the required packages
library(foreach)
library(doParallel)
library(deSolve)

# model
source('model.R')

# lists of parameters
source('parameters.R')

# lists of initial conditions/state variables
source('state_variables.R')

# Set up parallel backend to use multiple processors
numCores <- detectCores() - 1
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Run the ODE model in parallel for each combination of state_start_list and yfv_params_list
resultsNew <- foreach(yfv_params_idx = 1:length(yfv_params_list), .packages = 'deSolve') %:%
  foreach(state_start_idx = 1:length(state_start_list), .packages = 'deSolve') %dopar% {
    state_start <- state_start_list[[state_start_idx]]
    yfv_params <- yfv_params_list[[yfv_params_idx]]
    result <- as.data.frame(
      ode(
        y = unlist(state_start),  # Ensure state_start is a numeric vector
        times = times,
        func = yfv_model,
        parms = yfv_params
      )
    )
    list(
      yfv_params_idx = yfv_params_idx,
      state_start_idx = state_start_idx,
      result = result
    )
  }

# Stop the parallel backend
stopCluster(cl)

# View the results as a list of lists with identifiers
# print(results)
# plot.ts(results[[1]][[2]]$result$I_h)
# lines(results[[2]][[2]]$result$I_h)

# Assess results

# # simulate across a range of movements (<1-15%/month), inForst percent (10-70%), 
# # and mortality rates (50 & 80%)
# 
# pdf('../Figures/model_v_data.pdf', height = 4, width = 10)
# par(mfrow = c(1,2), mar = c(2.5,4,1,1))
# plot(yfv_epidemic, out$I_h, type = 'l', ylab = 'Inf people')
# lines(realdata$Date, realdata$MG_human, col = 'orange', type = 'b', pch = 16)
# plot(yfv_epidemic, out$I_p, type = 'l', ylab = 'Inf primates')
# lines(realdata$Date, realdata$MG_primate, col = 'orange', type = 'b', pch = 16)
# dev.off()


out_fixed_bite <- as.data.frame(
  ode(
    state_start
    , times
    , yfv_model
    , yfv_params_bite_fixed
  )
)

out_fixed_move <- as.data.frame(
  ode(
    state_start
    , times
    , yfv_model
    , yfv_params_move_fixed
  )
)

out_fixed <- as.data.frame(
  ode(
    state_start
    , times
    , yfv_model
    , yfv_params_fixed
  )
)

rho = 0.1

# overlaid plot
## NOTES: seasonally driven movement creates small second peak but at wrong time
# seasonally driven biting rate doesn't do anything, BUT
# seasonally driven movement and biting rate has synergistic response producing second and larger peak than
# either in isolation
pdf('../Figures/model_comparison_overlaid_v2.pdf', width = 12, height = 5)
par(mfrow = c(1,2), mar = c(2.5,4,1,1))
plot(yfv_epidemic, out$I_h*rho, type = 'l', ylab = 'Infected people', col = 'darkred', lwd = 2, ylim = c(0, 600))
lines(yfv_epidemic, out_fixed_move$I_h*rho, col = '#003285', lwd = 2, lty = 2)
lines(yfv_epidemic, out_fixed_bite$I_h*rho, col = '#2A629A', lwd = 2)
lines(yfv_epidemic, out_fixed$I_h*rho, col = '#FF7F3E', lwd = 2, lty = 2)
lines(realdata$Date, realdata$MG_human, type = 'b', pch = 16)
legend('topleft'
       , legend = c("Seasonal primate movement & seasonal biting rate", "No primate movement, seasonally-driven biting rate", "Seasonally-driven primate movement, fixed biting rate", "No primate movement & fixed biting rate", 'Data')
       , bty = 'n'
       , col = c("darkred", "#003285", "#2A629A", "#FF7F3E", 'black')
       , lwd = 2
       , lty = c(1, 2, 1, 2, 2)
       , pch=c(26,26,26,26,21))
plot(yfv_epidemic, out$I_p, type = 'l', ylab = 'Infected primates', col = 'darkred', lwd = 2, ylim = c(0, 600))
lines(yfv_epidemic, out_fixed_move$I_p, col = '#003285', lwd = 2, lty = 2)
lines(yfv_epidemic, out_fixed_bite$I_p, col = '#2A629A', lwd = 2)
lines(yfv_epidemic, out_fixed$I_p, col = '#FF7F3E', lwd = 2, lty = 2)
lines(realdata$Date, realdata$MG_primate, type = 'b', pch = 16)
dev.off()


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

