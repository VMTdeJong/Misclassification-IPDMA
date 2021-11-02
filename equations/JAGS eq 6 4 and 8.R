# Note: this is a BUGS / JAGS file. It's only saved as .R to get fancy colours in Rstudio.

### Data
## Description
# v = covariate, v = 1,...,V, where 
# V = is number of covariates z.
# z = covariate, not missing, no error assumed.
# i = individual, i = 1,.., n
# n = sample size
# y = outcome, binary, no error assumed. May be missing in some i.
# x = exposure of interest, possibly missing in some i.
# s = surrogate observation of x, i.e. with error. May be missing in some i.
# j = study/cluster, j = 1,...,J
# J = number of studies/clusters.
# Note that values of x, y, or z may be missing and are automatically imputed.

## Expected data format:
# var N, J, V, x[N], y[N], s[N], j[N], z[N, V];

### The model
## Coefficients
# beta   = coefficient in outcome model
# gamma  = coefficient in exposure model
# lambda = coefficient in measurement model, if x_ij = 1
# phi    = coefficient in measurement model, if x_ij = 0
# tau    = precision of random effect

## First subscript
# _0 = intercept
# _x = coefficient for x
# _z = coefficient for z

## Second subscript
# _0 = summary effect / fixed effect
# _j = random effect

## Tau's third subscript:
# tau has one additional subscript (the first) to indicate the coefficient it 
# is the variance of.
# tau is equal to 1/sqrt(prec)
# T(0,) truncates the distribution at zero, so that it is nonnegative
# The sampler assumes this is a-priori known to be nonnegative
# http://www.stat.columbia.edu/~gelman/research/published/taumain.pdf
# https://stackoverflow.com/questions/34935606/cauchy-prior-in-jags


model{
  ### The model
  for (i in 1:N) {
    # 1. Measurement model  
    #    Stratified by x.
    s[i] ~ dbern(s_p[i])
    logit(s_p[i]) <- ifelse(x[i] == 1,
                            lambda_0_0 + lambda_0_j[j[i]] + inprod(lambda_z_0, z[i, ]),
                            phi_0_0    + phi_0_j[j[i]]    + inprod(phi_z_0,    z[i, ]))
    
    # 2. Exposure model
    #    Note that x may be missing, in those cases it is latent and imputed.
    x[i] ~ dbern(x_p[i])
    logit(x_p[i]) <- gamma_0_0 + gamma_0_j[j[i]] + inprod(gamma_z_0, z[i, ])
    
    # 3. Outcome model
    y[i] ~ dbern(y_p[i])
    logit(y_p[i]) <-  beta_0_0 + beta_0_j[j[i]] + beta_x_0 * x[i] + beta_x_j[j[i]] * x[i] + inprod(beta_z_0, z[i, ])
  }
  
  # ### Monitoring
  # for (i in 1:N) {
  #   x_pos[i] <- ifelse(x[i] == 1, s[i], 0)
  #   x_neg[i] <- ifelse(x[i] == 0, s[i], 0)
  #   
  #   y_pos[i] <- ifelse(y[i] == 1, x[i], 0)
  #   y_neg[i] <- ifelse(y[i] == 0, x[i], 0)
  # }
  
  # min_sum <- .0001 # if x_sum is zero, we would otherwise be dividing by zero.
  # 
  # x_sum <- sum(x)
  # x_sum2 <- ifelse(x_sum < min_sum, min_sum, x_sum)
  # sens_x_s <- sum(x_pos) / x_sum2
  # spec_x_s <- 1 - sum(x_neg) / (N - x_sum2)
  # 
  # y_sum <- sum(y)
  # y_sum2 <- ifelse(y_sum < min_sum, min_sum, y_sum)
  # sens_y_x <- sum(y_pos) / y_sum2
  # spec_y_x <- 1 - sum(y_neg) / (N - y_sum2)
  
  ###  Priors
  # 1. Measurement model
  # x = 1
  lambda_0_0  ~ dnorm(0.0, precision)
  for (v in 1:V) {
    lambda_z_0[v] ~ dnorm(0.0, precision)
  }
  # Random intercept
  for (jj in 1:J) {
    lambda_0_j[jj] ~ dnorm(0.0, prec_lambda_0_j)
  }
  prec_lambda_0_j <- 1/tau2_lambda_0_j
  tau_lambda_0_j <- sqrt(tau2_lambda_0_j)
  tau2_lambda_0_j ~ dnorm(0.0, nu)T(0,)
  
  # x = 0
  phi_0_0  ~ dnorm(0.0, precision)
  for (v in 1:V) {
    phi_z_0[v] ~ dnorm(0.0, precision)
  }
  # Random intercept
  for (jj in 1:J) {
    phi_0_j[jj] ~ dnorm(0.0, prec_phi_0_j)
  }
  prec_phi_0_j <- 1/tau2_phi_0_j
  tau_phi_0_j <- sqrt(tau2_phi_0_j)
  tau2_phi_0_j ~ dnorm(0.0, nu)T(0,)
  
  
  # 2. Exposure model
  gamma_0_0  ~ dnorm(0.0, precision)
  for (v in 1:V) {
    gamma_z_0[v] ~ dnorm(0.0, precision)
  }
  # Random intercepts
  for (jj in 1:J) {
    gamma_0_j[jj] ~ dnorm(0.0, prec_gamma_0_j)
  }
  prec_gamma_0_j <- 1/tau2_gamma_0_j
  tau_gamma_0_j <- sqrt(tau2_gamma_0_j)
  tau2_gamma_0_j ~ dnorm(0.0, nu)T(0,)
  
  
  # 3. Outcome model
  beta_0_0  ~ dnorm(0.0, precision)
  beta_x_0  ~ dnorm(0.0, precision)
  
  for (v in 1:V) {
    beta_z_0[v] ~ dnorm(0.0, precision)
  }
  # Random intercepts
  for (jj in 1:J) {
    beta_0_j[jj] ~ dnorm(0.0, prec_beta_0_j)
  }
  prec_beta_0_j <- 1/tau2_beta_0_j
  tau_beta_0_j <- sqrt(tau2_beta_0_j)
  tau2_beta_0_j ~ dnorm(0.0, nu)T(0,)
  
  # Random effects
  for (jj in 1:J) {
    beta_x_j[jj] ~ dnorm(0.0, prec_beta_x_j)
  }
  prec_beta_x_j <- 1/tau2_beta_x_j
  tau_beta_x_j <- sqrt(tau2_beta_x_j)
  tau2_beta_x_j ~ dnorm(0.0, nu)T(0,)
  
  ### Hyper priors
  precision <- .1
  nu <- 5
}
