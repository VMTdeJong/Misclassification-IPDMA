bayes_re_logit <- function(data, cca = FALSE) {
  
  if (cca) {
    data <- data[!apply(data, 1, anyNA), ]
    data$j <- as.numeric(as.factor(data$j))
  }
  
  dat <- list(N = nrow(data),
              x = data$x,
              y = data$y,
              z = as.matrix(data[, grep("z", colnames(data)), drop = FALSE]),
              V = length(grep("z", colnames(data))),
              j = data$j,
              J = length(unique(data$j)))
  
  params <- Hmisc::Cs(beta_0_0,
                      beta_x_0,
                      beta_z_0,
                      
                      tau_beta_0_j,
                      tau_beta_x_j
  )
  
  load.module("glm", quiet = FALSE)
  
  init_params <- c(params, 
                   Hmisc::Cs(
                     beta_x_0,
                     beta_x_j))
  
  inits <- make_inits(data = data, n_chains = n_chains, params = init_params)
  
  mod <- jags.model(file = "equations/JAGS bayes re logit.R",
                    data = dat,
                    inits = inits,
                    n.chains = n_chains,
                    n.adapt = n_adapt)
  
  update(mod, n.iter = n_warmup)
  
  coda.samples(model = mod,
               variable.names = params,
               n.iter = n_iter,
               thin = n_thin,
  )
}
