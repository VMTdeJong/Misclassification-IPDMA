adma <- function(data) {
  f <- formula(paste(c("y ~ x", colnames(data)[grep("z", colnames(data))]), collapse = " + "))
  all_models <- list()
  i <- 0
  for (j in sort(unique(data$j))) {
    i <- i + 1
    all_models[[i]] <- logistf::logistf(f, data = data[data$j == j, ])
  }
  
  vars <- function(object) diag(vcov(object))
  coefs <- sapply(all_models, coef)["x", ]
  variances <- sapply(all_models, vars)["x", ]
  
  dat <- list(beta_est = coefs, 
              V = variances, 
              J = length(unique(data$j)))
          
   params <- Hmisc::Cs(beta_x_0,
                      beta_x_j,
                      
                      prec_beta_x_j,
                      tau_beta_x_j
  )
  
  load.module("glm", quiet = FALSE)
  
  inits <- make_inits(data = data, n_chains = n_chains, params = c(params))
  
  mod <- jags.model(file = "equations/JAGS AD-MA.R",
                    data = dat,
                    inits = inits,
                    n.chains = n_chains, 
                    n.adapt = n_adapt)
  
  update(mod, n.iter = n_warmup)
  
  coda.samples(model = mod,
               variable.names = params,
               n.iter = n_iter,
               n.burnin = n_warmup,
               thin = n_thin,
  )
}
