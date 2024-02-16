p <- read.csv('parameter_values.csv')

alpha_adj <- function(x, min_biting){
  spei <- c(1, 2.5)
  biting <- c(min_biting, 0.86)
  mod <- lm(biting ~ spei)
  
  if(x < 1){
    a <- min_biting
  } else {
    a0 <- predict(mod, newdata = data.frame(spei = x))
    a <- min(a0, 1)
  }
  return(a)
}

a1 <- function(x){
  alpha_adj(x, min_biting = p$value[p$variable=='a1'])
}

a2 <- function(x){
  alpha_adj(x, min_biting = p$value[p$variable=='a2'])
}

a3 <- function(x){
  alpha_adj(x, min_biting = p$value[p$variable=='a3'])
}
