# starting values for state variables
sprop = 0.999
Iprop = 0.001
# eprop = rprop = 0
humans = 2700000#20000000
monkeys = humans/3000 #16500 #15-20,000
mosquitoes = 25000
marmosets = monkeys * 2
none = 0
vaccinated = 0.5
natimm = 0.5

state_start <- c(
  S_p = monkeys*sprop #- (monkeys * natimm)
  , I_p = monkeys*Iprop
  , R_p = none #monkeys * natimm
  , S_c = marmosets*sprop
  , I_c = marmosets*Iprop
  , R_c = none
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