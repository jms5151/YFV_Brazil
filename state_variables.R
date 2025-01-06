# read in initial conditions
IC <- read.csv('init_conditions.csv')

# proportions for state variable initial conditions
sprop = 0.999
Iprop = 0.001
none = 0

# state variables that do not change
humans = IC$mean_value[IC$variable=='humans']
vaccinated = IC$mean_value[IC$variable=='vaccinated']

# set up N initial conditions
N = 500
mosquito_n = rnorm(N, IC$mean_value[IC$variable=='mosquitoes'], IC$sd[IC$variable=='mosquitoes'])
hm_prop_n = runif(N, IC$minimum[IC$variable == 'hm_prop'], IC$maximum[IC$variable == 'hm_prop'])
monkey_n = rnorm(N, IC$mean_value[IC$variable=='monkeys'], IC$sd[IC$variable=='monkeys'])
marmoset_n = monkey_n * 2
inForest_n = rnorm(N, IC$mean_value[IC$variable=='inForest'], IC$sd[IC$variable=='inForest']) # make sure value doesn't exceed 1
inForest_n[inForest_n >= 1] <- 0.99


# function to create state_start list for each set of initial conditions
create_state_start <- function(mosquitoes, hm_prop, monkeys, marmosets, inForest) {
  list(
    S_p = monkeys * (1 - inForest),
    I_p = monkeys * Iprop,
    R_p = monkeys * inForest,
    S_c = marmosets * sprop,
    I_c = marmosets * Iprop,
    R_c = none,
    S_h = humans * (1 - vaccinated),
    E_h = none,
    I_h = none,
    R_h = humans * vaccinated,
    S_hm = mosquitoes * hm_prop * sprop,
    E_hm = none,
    I_hm = mosquitoes * hm_prop * Iprop,
    S_aa = mosquitoes,
    E_aa = none,
    I_aa = none
  )
}

# create the list of lists
state_start_list <- vector("list", N)
for (i in 1:N) {
  state_start_list[[i]] <- create_state_start(mosquito_n[i], hm_prop = hm_prop_n[i], monkey_n[i], marmoset_n[i], inForest_n[i])
}
