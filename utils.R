get_summary_rjags <- function(object) {
  s <- summary(object)
  ss <- as.data.frame(s$statistics)
  sq <- as.data.frame(s$quantiles)
  data.frame(mean   = ss["beta_x_0", "Mean"],
             se     = ss["beta_x_0", "SD"],
             ci.lb  = sq[ "beta_x_0", "2.5%"],
             median = sq[ "beta_x_0", "50%"],
             ci.ub  = sq[ "beta_x_0", "97.5%"],
             tau.mean   = ss["tau_beta_x_j", "Mean"],
             tau.se     = ss["tau_beta_x_j", "SD"],
             tau.ci.lb  = sq[ "tau_beta_x_j", "2.5%"],
             tau.median = sq[ "tau_beta_x_j", "50%"],
             tau.ci.ub  = sq[ "tau_beta_x_j", "97.5%"]
             )
}

# Creates and saves a new seed
#
# These seeds use the L'Ecuyer-CMRG RNG, which is designed for parallel usage.

# new_seed <- function(path, default.seed = 2020) { # replaced with make_seeds() which makes all necessary seeds in one go.
#   
#   library(parallel)
#   RNGkind("L'Ecuyer-CMRG")
#   
#   if (!file.exists(path)) {
#     set.seed(default.seed)
#     seeds <- as.data.frame(matrix(.Random.seed, nrow = 1))
#   } else {
#     load(path)
#     s <- nextRNGStream(unlist(seeds[nrow(seeds), ]))
#     seeds <- rbind(seeds, s)
#   }
#   save(seeds, file = path)
#   return(seeds[nrow(seeds), ])
# }

see_seeds <- function(path) {
  load(path)
  return(seeds)
}

make_seeds <- function(path, I = 11000, default.seed = 2021) {
  library(parallel)
  RNGkind("L'Ecuyer-CMRG")
  
  if (!file.exists(path)) {
    set.seed(default.seed)
    seeds <- as.data.frame(matrix(.Random.seed, nrow = 1))
  } else {
    load(path)
  }
  
  if (I > 1) 
    for (i in 2:I)
      seeds <- rbind(seeds, nextRNGStream(unlist(seeds[nrow(seeds), ])))
  
  save(seeds, file = path)
  return(seeds)
}

get_seed <- function(path, i) {
  load(path)
  return(seeds[i, ])
}

#' @usage 
#' new_seed()
#' see_seeds()

make_params <- function(i) {
  sce <- i %% 11 + 1
  re <- i %/% 11
  I <- 1
  
  # All possible parameter values
  # J_gold, 3, 5, 7, or 9
  # n: 300, 500, 700, 900 or 1100
  # nz: 1, 2, 3, 4.
  
  # Default parameter values (sce = 1)
  J_gold <- 5 
  n <- 500
  nz <- 1
  
  # and a fixed one
  J <- 10
  
  # Alternate parameter values
  # Number of gold standard studies J_gold
  if (sce == 2)
    J_gold <- 3
  if (sce == 3)
    J_gold <- 7
  if (sce == 4)
    J_gold <- 9
  
  # Number of participants per study
  if (sce == 5)
    n <- 100
  if (sce == 6)
    n <- 300
  if (sce == 7)
    n <- 700
  
  # Number of covariates
  if (sce == 8)
    nz <- 2
  if (sce == 9)
    nz <- 3
  if (sce == 10)
    nz <- 4
  
  # Covariate strength
  b.z <- rep(0, nz)
  
  if (sce == 11)
    b.z <- rep(c(1, -1), nz)[seq_len(nz)]
  
  return(list(sce = sce,
              re = re,
              I = I,
              J_gold = J_gold,
              n = n,
              nz = nz,
              b.z = b.z,
              J = J))
}

