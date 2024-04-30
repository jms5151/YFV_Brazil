library(deSolve)

source('biting_rate_drought_functions.R')

yfv_model <- function(t, state, parameters) {
  with(as.list(c(state,parameters)), {
    # browser() # can use in if t > x
    # Howler monkeys
    # a1(t)
    dS_p <- -(a1[t] * 0.5 * p * b * (I_hm / (S_p + I_p + R_p)) * S_p) + (mu_p * (I_p + R_p)) + (m[t] * R_p)
    dI_p <- (a1[t] * 0.5 * p * b * (I_hm / (S_p + I_p + R_p)) * S_p) - (gamma_p * I_p) - (mu_p * I_p) - (I_p * mu_v1)
    dR_p <- ((gamma_p * I_p) * (1-mu_v1)) - (mu_p * R_p) - (m[t] * R_p)
    
    # Marmosets
    dS_c <- -(a1[t] * 0.5 * p * b * (I_hm / (S_c + I_c + R_c)) * S_c) + (mu_c * (I_c + R_c))
    dI_c <- (a1[t] * 0.5 * p * b * (I_hm / (S_c + I_c + R_c)) * S_c) - (gamma_c * I_c) - (mu_c * I_c) 
    dR_c <- ((gamma_c * I_c) * (1-mu_v3)) - (mu_c * R_c)
    
    # Humans
    dS_h <- -((a2[t] * (1-p) * b * (I_hm / (S_h + E_h + I_h + R_h))) * S_h) - ((a3[t] * b * (I_aa / (S_h + E_h + I_h + R_h))) * S_h) - (V[t] * S_h) + (mu_h * (E_h + I_h + R_h)) #+ (w * R_h)
    dE_h <- ((a2[t] * (1-p) * b * (I_hm / (S_h + E_h + I_h + R_h))) * S_h)  + ((a3[t] * b * (I_aa / (S_h + E_h + I_h + R_h))) * S_h) - (delta_h * E_h) - (mu_h * E_h) - (V[t] * E_h)
    dI_h <- (delta_h * E_h) - (gamma_h * I_h) - (mu_h * I_h) - (V[t] * I_h)
    dR_h <- ((gamma_h * I_h) * (1-mu_v2)) + (V[t] * S_h) - (mu_h * R_h) + (V[t] * (S_h + E_h + I_h)) #- (w * R_h)
    
    # Haemagogus mosquitoes
    dS_hm <- -((a1[t] * 0.5 * p * pMI1 * (I_p / (S_p + I_p + R_p))) +
                 (a2[t] * (1-p) * pMI2 * (I_h / (S_h + E_h + I_h + R_h))) +
                 (a1[t] * 0.5 * p * pMI4 * (I_c/(S_c + I_c + R_c)))) * S_hm - (mu_hm * S_hm) + max(0, min(br1[t] * (S_hm + E_hm + I_hm), ((1 - ((S_hm + E_hm + I_hm)/K[t])) * (S_hm + E_hm + I_hm))))
    dE_hm <- ((a1[t] * 0.5 * p * pMI1 * (I_p / (S_p + I_p + R_p))) +
                (a2[t] * (1-p) * pMI2 * (I_h / (S_h + E_h + I_h + R_h))) +
                (a1[t] * 0.5 * p * pMI4 * (I_c/(S_c + I_c + R_c)))) * S_hm - ((PDR_hm + mu_hm) * E_hm)
    dI_hm <- (PDR_hm * E_hm) - (mu_hm * I_hm)
    
    # Aedes aegypti mosquitoes
    dS_aa <- -((a3[t] * pMI3 * (I_h / (S_h + E_h + I_h + R_h))) * S_aa) - (mu_aa * S_aa) + max(0, min(br2[t] * (S_aa + E_aa + I_aa), ((1 - ((S_aa + E_aa + I_aa)/K[t])) * (S_aa + E_aa + I_aa))))
    dE_aa <- ((a3[t] * pMI3 * (I_h / (S_h + E_h + I_h + R_h))) * S_aa) - ((PDR_aa + mu_aa) * E_aa)
    dI_aa <- (PDR_aa * E_aa) - (mu_aa * I_aa)
    
    list(c(dS_p, dI_p, dR_p, dS_c, dI_c, dR_c, dS_h, dE_h, dI_h, dR_h, dS_hm, dE_hm, dI_hm, dS_aa, dE_aa, dI_aa), 'a1' = (a1[t]), 'a2' = (a2[t]), 'a3' = (a3[t]), 'k' = k[t], 'vaccinated' = V[t] * (S_h + E_h + I_h + R_h))
  })
} 

