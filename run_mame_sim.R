c_args <- commandArgs(trailingOnly = T) # These are the arguments given through LINUX 

if (length(c_args) == 0) {
  # These args are for testing only.
  i <- 11
  
} else {
  # These args are for the real sim.
  i <- as.numeric(c_args[1])

  setwd("/hpc/shared/julius_bs/VMTdeJong/MAME/MAME_sim_2021")
} 

# JAGS parameter values
n_chains <- 2
n_thin <- 2
n_warmup <- 5000
n_adapt <- 1000
n_iter <- 1e4

# simulation parameters
source("utils.R")
p <- make_params(i)
 
print(paste0("Starting a run with scenario = ", p$sce, ", I = ", p$I, ", re = ", p$re, 
             ", n = ", p$n, ", J = ", p$J, ", J_gold = ", p$J_gold, ", number of covariates nz = ",
             p$nz, "."))



source("sampling.R")
source("sim.R")
seed_path <- "seed/parallel_seed_stream.RData"

# Sampling
sim_sampler <- function(a.tau = 1/4, b.x.tau = 0.15, nz = p$nz, center.z = F, n = p$n, J = p$J, b.z = p$b.z){
  sample_misclass_re(n = n, 
                     J = J, 
                     a = -1/2, 
                     a.tau = a.tau,
                     b.x = 1, 
                     b.x.tau = b.x.tau,
                     b.z = b.z,
                     z.mean.tau = 1/2,
                     cov = 1/4,
                     gamma.0 = -1/4, 
                     gamma.0.tau = 1/4,
                     gamma.z = rep(c(1, -1), nz)[seq_len(nz)], 
                     lambda.0.0 = 3,     
                     lambda.0.j.tau = 1, 
                     lambda.z = rep(c(1, -1), nz)[seq_len(nz)],     
                     phi.0.0 = -3,     
                     phi.0.j.tau = 1, 
                     phi.z = rep(c(1, -1), nz)[seq_len(nz)],         
                     center.z = center.z
  )
}

FUN_sampler_sim <- function(seed = .Random.seed) {
  .Random.seed <<- seed
  sim_sampler
}

# Run the simulation
system.time(den_sim(I = p$I, get_seed(seed_path, i), n = p$n, J_gold = p$J_gold, J = p$J, nz = p$nz, save_first_fit = i <= 11))

