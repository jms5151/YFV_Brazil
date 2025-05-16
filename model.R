yfv_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Time-varying parameters
    alpha_val <- if (is.function(alpha)) alpha(t) else alpha
    V_val     <- if (is.function(V)) V(t) else V
    K_val     <- if (is.function(K)) K(t) else K
    birth_hg_val <- if (is.function(birth_hg)) birth_hg(t) else birth_hg
    birth_aa_val <- if (is.function(birth_aa)) birth_aa(t) else birth_aa
    monkey_seed_val <- if (is.function(monkey_seed)) monkey_seed(t) else monkey_seed
    marmoset_seed_val <- if (is.function(marmoset_seed)) marmoset_seed(t) else marmoset_seed
    hg_seed_val <- if (is.function(hg_seed)) hg_seed(t) else hg_seed
    
    # Prevent mosquito collapse
    Shg <- pmax(Shg, 1)
    Saa <- pmax(Saa, 1)
    
    # Total population sizes
    Np  <- Sp + Ip + Rp
    Nc  <- Sc + Ic + Rc
    Nh  <- Sh + Eh + Ih + Rh
    Nhg <- Shg + Ehg + Ihg
    Naa <- Saa + Eaa + Iaa
    
    # Force of infection (lambda)
    lambda_hg <- alpha_val * (
      pMI1 * Ip / Np +
        pMI2 * Ih / Nh +
        pMI4 * Ic / Nc
    )
    
    lambda_h <- alpha_val * (
      b * Ihg / Nh + Iaa / Nh
    )
    
    lambda_aa <- alpha_val * pMI3 * Ih / Nh
    
    # Howler monkeys
    dSp <- -alpha_val * b * Ihg / Np * Sp + mu_p * (Np - Sp) + monkey_seed_val
    dIp <-  alpha_val * b * Ihg / Np * Sp - gamma_p * Ip - mu_p * Ip - mu_v1 * Ip
    dRp <-  gamma_p * Ip - mu_p * Rp
    
    # Marmosets
    dSc <- -alpha_val * b * Ihg / Nc * Sc + mu_c * (Nc - Sc) + marmoset_seed_val
    dIc <-  alpha_val * b * Ihg / Nc * Sc - gamma_c * Ic - mu_c * Ic - mu_v3 * Ic
    dRc <-  gamma_c * Ic - mu_c * Rc
    
    # Humans
    dSh <- -lambda_h * Sh + mu_h * (Nh - Sh) - V_val * Sh
    dEh <-  lambda_h * Sh - delta_h * Eh - mu_h * Eh - V_val * Eh
    dIh <-  delta_h * Eh - gamma_h * Ih - mu_h * Ih - mu_v2 * Ih - V_val * Ih
    dRh <-  gamma_h * Ih - mu_h * Rh + V_val * (Sh + Eh + Ih)
    
    # Haemagogus mosquitoes
    dShg <- -lambda_hg * Shg + birth_hg_val * (Nhg - Shg) * (1 - Nhg / K_val) + hg_seed_val
    dEhg <-  lambda_hg * Shg - (PDR_hg + mu_hg) * Ehg
    dIhg <-  PDR_hg * Ehg - mu_hg * Ihg
    
    # Aedes mosquitoes
    dSaa <- -lambda_aa * Saa + birth_aa_val * (Naa - Saa) * (1 - Naa / K_val)
    dEaa <-  lambda_aa * Saa - (PDR_aa + mu_aa) * Eaa
    dIaa <-  PDR_aa * Eaa - mu_aa * Iaa
    
    list(c(dSp, dIp, dRp, dSc, dIc, dRc,
           dSh, dEh, dIh, dRh,
           dShg, dEhg, dIhg,
           dSaa, dEaa, dIaa))
  })
}

