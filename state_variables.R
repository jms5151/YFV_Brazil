# read in initial conditions
IC <- read.csv('init_conditions.csv')

# proportions for state variable initial conditions
Sprop = 0.999
Iprop = 0.001
none = 0

# state variables that do not change
humans = IC$mean_value[IC$variable=='humans']
vaccinated = IC$mean_value[IC$variable=='vaccinated']

# set up N initial conditions
N = 500
mosquito_n = rnorm(N, IC$mean_value[IC$variable=='mosquitoes'], IC$sd[IC$variable=='mosquitoes'])
hg_prop_n = runif(N, IC$minimum[IC$variable == 'hg_prop'], IC$maximum[IC$variable == 'hg_prop'])
monkey_n = rnorm(N, IC$mean_value[IC$variable=='monkeys'], IC$sd[IC$variable=='monkeys'])
marmoset_n = monkey_n * 1.5

# function to create state_start list for each set of initial conditions
create_state_start <- function(mosquitoes, hg_prop, monkeys, marmosets) {
  list(
    Sp = monkeys,
    Ip = none,
    Rp = none,
    Sc = marmosets,
    Ic = none,
    Rc = none,
    Sh = humans * (1 - vaccinated),
    Eh = none,
    Ih = none,
    Rh = humans * vaccinated,
    Shg = mosquitoes * hg_prop * Sprop,
    Ehg = none,
    Ihg = none,
    Saa = mosquitoes,
    Eaa = none,
    Iaa = none
  )
}

# create the list of lists
state_start_list <- vector("list", N)
for (i in 1:N) {
  state_start_list[[i]] <- create_state_start(mosquito_n[i], hg_prop = hg_prop_n[i], monkey_n[i], marmoset_n[i])
}

# create event function for initial infected import
# event_function_importI <- function(t, state, parameters) {
#   if (abs(t - 15) < 1e-6) {
#     state['Ic'] <- state['Ic'] + 10  
#     # state['Ip'] <- state['Ip'] + 1  
#     # state['Ihg'] <- state['Ihg'] + 5  
#     
#   }
#   return(state)
# }

