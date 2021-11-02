# Expected data:
# beta_est[1:J], V[1:J], J

model { 
  for (j in 1:J) {
    beta_est[j] ~ dnorm(beta_x_0 + beta_x_j[j], 1/V[j]) 
    beta_x_j[j] ~ dnorm(0, prec_beta_x_j) 
  }   

  ## Priors
  prec_beta_x_j <- 1/tau2_beta_x_j
  tau_beta_x_j <- sqrt(tau2_beta_x_j)
  tau2_beta_x_j ~ dnorm(0.0, nu)T(0,)
  
  beta_x_0 ~ dnorm(0, precision)
  
  ##
  precision <- .1
  nu <- 5
}

