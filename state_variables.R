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
N = 1
mosquito_n = 2500000#rnorm(N, IC$mean_value[IC$variable=='mosquitoes'], IC$sd[IC$variable=='mosquitoes'])
hg_prop_n = 0.4#runif(N, IC$minimum[IC$variable == 'hg_prop'], IC$maximum[IC$variable == 'hg_prop'])
monkey_n = 20000#rnorm(N, IC$mean_value[IC$variable=='monkeys'], IC$sd[IC$variable=='monkeys'])
marmoset_n = monkey_n * 2
# may need this to be MUCH lower, like 0.05
# inForest_n = rnorm(N, IC$mean_value[IC$variable=='inForest'], IC$sd[IC$variable=='inForest']) # make sure value doesn't exceed 1
# inForest_n[inForest_n >= 1] <- 0.99
# inForest_n = 0.1

humans = 500000
mosquitoes = 30000
hg_prop = 0.2
monkeys = 500
marmosets = monkeys * 1.5



# function to create state_start list for each set of initial conditions
create_state_start <- function(mosquitoes, hg_prop, monkeys, marmosets) {
  list(
    Sp = monkeys * Sprop,
    Ip = monkeys * Iprop,
    Rp = none,
    Sc = marmosets * Sprop,
    Ic = marmosets * Iprop,
    Rc = none,
    Sh = humans * (1 - vaccinated),
    Eh = none,
    Ih = none,
    Rh = humans * vaccinated,
    Shg = mosquitoes * hg_prop * Sprop,
    Ehg = none,
    Ihg = mosquitoes * hg_prop * Iprop,
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

