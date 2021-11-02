eq648 <- function(data) {
  dat <- list(N = nrow(data),
              s = data$s,
              x = data$x,
              y = data$y,
              z = as.matrix(data[, grep("z", colnames(data)), drop = FALSE]),
              V = length(grep("z", colnames(data))),
              j = data$j,
              J = length(unique(data$j)))
  
  params <- Hmisc::Cs(lambda_0_0,
                      lambda_z_0,
                      
                      phi_0_0,
                      phi_z_0,
                      
                      gamma_0_0,
                      gamma_z_0,
                      
                      beta_0_0,
                      beta_x_0,
                      beta_z_0,
                      
                      # sens_x_s,
                      # spec_x_s,
                      # sens_y_x,
                      # spec_y_x,
                      
                      tau_gamma_0_j,
                      tau_lambda_0_j,
                      tau_phi_0_j,
                      tau_beta_0_j,
                      tau_beta_x_j
  )
  
  load.module("glm", quiet = FALSE)
  
  init_params <- c(params, 
                   Hmisc::Cs(
                     lambda_0_j,
                     phi_0_j,
                     gamma_0_j,
                     beta_x_0,
                     beta_x_j))
  
  inits <- make_inits(data = data, n_chains = n_chains, params = init_params)
  
 
  
  converged <- FALSE
  
  while (!converged) {
    mod <- jags.model(file = "equations/JAGS eq 6 4 and 8.R",
                      data = dat,
                      inits = inits,
                      n.chains = n_chains, 
                      n.adapt = n_adapt)
    
    update(mod, n.iter = n_warmup)
    samples <- coda.samples(model = mod,
                            variable.names = params,
                            n.iter = n_iter,
                            thin = n_thin,
    )
    n_iter <- n_iter * 2
    mpsrf <- gelman.diag(samples)$mpsrf
    print(mpsrf)
    converged <- mpsrf < 1.05
  }
  
  samples
}
