# starting values for state variables
sprop = 0.999
Iprop = 0.001
# eprop = rprop = 0
monkeys = 16500 #15-20,000
humans = 2700000#20000000
mosquitoes = 25000
none = 0
vaccinated = 0.5

state_start <- c(
  S_p = monkeys*sprop
  , I_p = monkeys*Iprop
  , R_p = none
  , S_h = humans*(1-vaccinated)
  , E_h = none
  , I_h = none
  , R_h = humans*vaccinated
  , S_hm = mosquitoes*sprop
  , E_hm = none
  , I_hm = mosquitoes*Iprop
  , S_aa = mosquitoes
  , E_aa = none
  , I_aa = none
)