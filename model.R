library(deSolve)

# add carrying capacity
# biting rate | dehydration

yfv_model <- function(t, state, parameters) {
  with(as.list(c(state,parameters)), {
    # Primates
    dS_p <- -(a1 * p *b * (I_hm / (S_p + I_p + R_p)) * S_p) + mu_p * (I_p + R_p)
    dI_p <- (a1 * p * b * (I_hm / (S_p + I_p + R_p)) * S_p) - (gamma_p * I_p) - (mu_v1 * I_p)
    dR_p <- gamma_p * I_p 
    
    # Humans
    dS_h <- -(a2 * (1-p) * b * (I_hm / (S_h + E_h + I_h + R_h))) * S_h - (a3 * b * (I_aa / (S_h + E_h + I_h + R_h))) * S_h - V * S_h + mu_h * (E_h + I_h)
    dE_h <- (a2 * (1-p) * b * (I_hm / (S_h + E_h + I_h + R_h))) * S_h  + (a3 * b * (I_aa / (S_h + E_h + I_h + R_h))) * S_h - (delta_h * E_h)
    dI_h <- (delta_h * E_h) - (gamma_h * I_h) - (mu_v2 * I_h)
    dR_h <- gamma_h * I_h + V * S_h
    
    # Haemagogus mosquitoes
    dS_hm <- (sigma_hm * (S_hm + E_hm + I_hm)) - (a1 * p * pMI1 * (I_p / (S_p + I_p + R_p)) + mu_hm) * S_hm - (a2 * pMI2 * (I_h / (S_h + E_h + I_h + R_h))) * S_hm * (1 - (S_hm + E_hm + I_hm)/K)
    dE_hm <- (a1 * p * pMI1 * (I_p / (S_p + I_p + R_p)) + mu_hm) * S_hm + (a2 * (1-p) * pMI2 * (I_h / (S_h + E_h + I_h + R_h))) * S_hm - (PDR_hm + mu_hm) * E_hm
    dI_hm <- PDR_hm * E_hm - mu_hm * I_hm
    
    # Aedes aegypti mosquitoes
    dS_aa <- sigma_aa * (S_aa + E_aa + I_aa) - (a3 * pMI3 * (I_h / (S_h + E_h + I_h + R_h)) + mu_aa) * S_aa  * (1 - (S_aa + E_aa + I_aa)/K)
    dE_aa <- (a3 * pMI3 * (I_h / (S_h + E_h + I_h + R_h)) + mu_aa) * S_aa - (PDR_aa + mu_aa) * E_aa
    dI_aa <- PDR_aa * E_aa - mu_aa * I_aa
    
    list(c(dS_p, dI_p, dR_p, dS_h, dE_h, dI_h, dR_h, dS_hm, dE_hm, dI_hm, dS_aa, dE_aa, dI_aa))
  })
} 

