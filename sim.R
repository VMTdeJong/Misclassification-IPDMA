# Functions
library(rjags)
library(abind)
library(logistf)

# Naive Bayesian random effects logit models
source("equations/bayes re logit.R") # Full data
source("equations/bayes re logit xors.R") # Naive model, calls same JAGS function as above

# Bayesian (random effects) logit models that account for misclassification
source("equations/eq 6 4 and 8.R")

# Aggregate data meta-analysis, for comparison.
source("equations/AD-MA.R")

# Initial values
source("inits.R") # requires nz and J to be in the global environment


# I Number of repetitions for the simulation
# seed seed to start with
# n sample size per center/study
# nz number of covariates z
# J number of centers/studies
# J_gold number of centers/studies that have recorded X

den_sim <- function(I, seed, n = 100, nz = 1, J = 10, J_gold = 5, save_first_fit = T) {
  seed_num <- row.names(seed)
  this_sim <- paste(" n_", n, " nz_", nz, " J_", J, " JGold_", J_gold,  " seed_", seed_num, sep = "")
  sampler <- FUN_sampler_sim(unlist(seed))
  
  # Sim parameters and initialization
  model_seeds <- data_seeds <- me_fit_list <- list()
  comparison_methods <- c("adma_gold", "adma_xors",
                          "ipd_gold", "ipd_xors")
  me_methods <- c("eq648")
  statistics <- c("mean", "se", "ci.lb", "median", "ci.ub", 
                  "tau.mean", "tau.se", "tau.ci.lb", "tau.median", "tau.ci.ub")
  
  estimates <- array(dim = c(length(statistics), length(comparison_methods) + length(me_methods), I),
                     dimnames = list(statistic = statistics,
                                     method = c(comparison_methods, me_methods),
                                     I = seq_len(I)))
  
  for (i in seq_len(I)) {
    model_seeds_current <- list()
    data_seeds[[i]] <- .Random.seed
    save(data_seeds, file = paste("output/data-seeds", this_sim,".RData", sep = ""))
    data_full <- sampler(n = n, nz = nz, J = J)
    data_miss <- add_missing(data_full, J_gold = J_gold)
    data_miss$x <- data_miss$x_miss
    data_miss$x_miss <- NULL
    
    ### For comparison: Methods that cheat or are naive
    begin_comparison <- proc.time()

    # x or s / naive
    model_seeds_current[["ipd_xors"]] <- .Random.seed
    print(paste("Starting method: ipd_xors"))
    ipd_xors_fit <- bayes_re_logit_xors(data_miss)
    estimates[ , "ipd_xors", i] <- unlist(get_summary_rjags(ipd_xors_fit))
    if (i == 1)
      me_fit_list[["ipd_xors"]] <- ipd_xors_fit

    # x or s adma / naive adma
    model_seeds_current[["adma_xors"]] <- .Random.seed
    print(paste("Starting method: adma_xors"))
    data_xors <- xors_replace(data_miss, "x", "s")
    adma_xors_fit <- adma(data_xors)
    estimates[ , "adma_xors", i] <- unlist(get_summary_rjags(adma_xors_fit))
    if (i == 1)
      me_fit_list[["adma_xors"]] <- adma_xors_fit

    # Gold only
    data_gold_bayes <- data_miss[!is.na(data_miss$x), ]
    J_saved <- J
    J <<- length(unique(data_gold_bayes$j)) # Because it uses J from the global env for the inits!
    model_seeds_current[["ipd_gold"]] <- .Random.seed
    print(paste("Starting method: ipd_gold"))
    ipd_gold_fit <- bayes_re_logit(data_gold_bayes, cca = TRUE)
    estimates[ , "ipd_gold", i] <- unlist(get_summary_rjags(ipd_gold_fit))
    if (i == 1)
      me_fit_list[["ipd_gold"]] <- ipd_gold_fit

    # AD-MA with gold only
    data_gold_bayes$s <- NULL
    model_seeds_current[["adma_gold"]] <- .Random.seed
    print(paste("Starting method: adma_gold"))
    adma_gold_fit <- adma(data_gold_bayes)
    estimates[ , "adma_gold", i] <- unlist(get_summary_rjags(adma_gold_fit))
    if (i == 1)
      me_fit_list[["adma_gold"]] <- adma_gold_fit

    J <<- J_saved # Set back to J for full data, for other methods.

    print(proc.time() - begin_comparison)

    # # # ### Bayesian error correction models
    for (me in me_methods) {
      model_seeds_current[[me]] <- .Random.seed
      begin_me <- proc.time()
      print(paste("Starting method:", me))
      me_fit <- match.fun(me)(data_miss)
      me_fit_saved <<- me_fit
      estimates[ , me, i] <- unlist(get_summary_rjags(me_fit))
      if (i == 1)
        me_fit_list[[me]] <- me_fit
      remove(me_fit)
      print(proc.time() - begin_me)
    }
    
    if (i == 1 && save_first_fit)
      save(me_fit_list, file = paste("output/fits_", this_sim, ".RData", sep = ""))
    model_seeds[[i]] <- model_seeds_current
    save(model_seeds, file = paste("output/model-seeds_", this_sim, ".RData", sep = ""))
    save(estimates, file = paste("output/est_", this_sim, ".RData", sep = ""))
  }
}



