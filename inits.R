
make_inits <- function(data, n_chains, params, init_limit = 2, ...) {
  inits <- list()
  for (chain in seq_len(n_chains)) {
    inits[[chain]] <- make_inits_single_chain(data,
                                              params = params,
                                              init_limit = init_limit,
                                              ...)
  }
  inits
}


make_inits_single_chain <- function(data,
                                        params,
                                        init_limit = 2,
                                        ...) {
  nz <- length(grep("z", colnames(data)))
  J <- length(unique(data$j))
  
  inits <- list(
    ### Coefficients
    ## Scalar
    eta_0_0    = runif(1, min = -init_limit, max = init_limit),
    theta_0_0  = runif(1, min = -init_limit, max = init_limit),
    psi_0_0    = runif(1, min = -init_limit, max = init_limit),
    omega_0_0  = runif(1, min = -init_limit, max = init_limit),
    lambda_0_0 = runif(1, min = -init_limit, max = init_limit),
    lambda_y_0 = runif(1, min = -init_limit, max = init_limit),
    phi_0_0    = runif(1, min = -init_limit, max = init_limit),
    phi_y_0    = runif(1, min = -init_limit, max = init_limit),
    gamma_0_0  = runif(1, min = -init_limit, max = init_limit),
    beta_0_0   = runif(1, min = -init_limit, max = init_limit),
    beta_x_0   = runif(1, min = -init_limit, max = init_limit),
    gamma_s_0  = runif(1, min = -init_limit, max = init_limit),
    
    ## Vectors
    eta_z_0    = runif(nz, min = -init_limit, max = init_limit),
    theta_z_0  = runif(nz, min = -init_limit, max = init_limit),
    psi_z_0    = runif(nz, min = -init_limit, max = init_limit),
    omega_z_0  = runif(nz, min = -init_limit, max = init_limit),
    lambda_z_0 = runif(nz, min = -init_limit, max = init_limit),
    phi_z_0    = runif(nz, min = -init_limit, max = init_limit),
    gamma_z_0  = runif(nz, min = -init_limit, max = init_limit),
    beta_z_0   = runif(nz, min = -init_limit, max = init_limit),
    
    eta_0_j    = runif(J, min = -init_limit, max = init_limit),
    theta_0_j  = runif(J, min = -init_limit, max = init_limit),
    psi_0_j    = runif(J, min = -init_limit, max = init_limit),
    omega_0_j  = runif(J, min = -init_limit, max = init_limit),
    lambda_0_j = runif(J, min = -init_limit, max = init_limit),
    phi_0_j    = runif(J, min = -init_limit, max = init_limit),
    gamma_0_j  = runif(J, min = -init_limit, max = init_limit),
    beta_0_j   = runif(J, min = -init_limit, max = init_limit),
    beta_x_j   = runif(J, min = -init_limit, max = init_limit),
    
    ### Variances
    ## Scalar
    tau2_gamma_0_j  = runif(1, min = .0001, max = init_limit),
    tau2_beta_0_j   = runif(1, min = .0001, max = init_limit),
    tau2_beta_x_j   = runif(1, min = .0001, max = init_limit),
    
    tau2_eta_0_j    = runif(1, min = .0001, max = init_limit),
    tau2_theta_0_j  = runif(1, min = .0001, max = init_limit),
    tau2_psi_0_j    = runif(1, min = .0001, max = init_limit),
    tau2_omega_0_j  = runif(1, min = .0001, max = init_limit),
    
    tau2_lambda_0_j = runif(1, min = .0001, max = init_limit),
    tau2_phi_0_j    = runif(1, min = .0001, max = init_limit)
  )
  
  inits[names(inits) %in% params]
}
