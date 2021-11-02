library(tidyverse)

#' @param x variable to draw bootstrap samples from
#' @param fun function to apply on each bootstrap sample. E.g. mean, median, rmse
#' @param bootstrap_fun function to estimate statistic to return. E.g. SD for SE, or quantile for percentile
#' @param ... args passed to bootstrap_fun
#' @param n number of bootstrap samples
#' @param seed seed passed to set.seed, so that multiple calls have the same seed (dplyr::summarise will
#' then have the same samples for each method). Set to NULL for random seed.

bootstrap <- function(x, fun, bootstrap_fun = sd, ..., n = 2000, seed = 12345) {
  fun <- match.fun(fun)
  bootstrap_fun <- match.fun(bootstrap_fun)
  out <- rep(NA, n)
  
  if (!is.null(seed)) {set.seed(seed)}
  
  for (i in seq_len(n)) {
    z <- sample(x, size = length(x), replace = TRUE)
    out[i] <- fun(z)
  }
  bootstrap_fun(out, ...)
}

summarise_results <- function(x, true_value = 1)
  x %>%
  group_by(sce, statistic, method) %>%
  summarise(mean = mean(estimate),
            SD = sd(estimate),
            mean_SE = SD/sqrt(n()),
            
            median = median(estimate),
            median_lb = bootstrap(estimate, median, quantile, .025),
            median_ub = bootstrap(estimate, median, quantile, .975),
            
            bias = mean(estimate - true_value),
            bias_lb = bootstrap(estimate - true_value, mean, quantile, .025),
            bias_ub = bootstrap(estimate - true_value, mean, quantile, .975),
            
            median_bias = median(estimate - true_value),
            median_bias_lb = bootstrap(estimate - true_value, median, quantile, .025),
            median_bias_ub = bootstrap(estimate - true_value, median, quantile, .975),
            
            RMSE = sqrt(mean((estimate - true_value)^2)),
            RMSE_lb = bootstrap((estimate - true_value)^2, function(z) sqrt(mean(z)), quantile, .025),
            RMSE_ub = bootstrap((estimate - true_value)^2, function(z) sqrt(mean(z)), quantile, .975),
            
            median_RMSE = sqrt(median((estimate - true_value)^2)),
            median_RMSE_lb = bootstrap((estimate - true_value)^2, function(z) sqrt(median(z)), quantile, .025),
            median_RMSE_ub = bootstrap((estimate - true_value)^2, function(z) sqrt(median(z)), quantile, .975),
            
            n = n())

summarise_coverage <- function(x, true_value = 1) {
  spread(x[x$statistic %in% c("ci.lb", "ci.ub"), ], statistic, estimate) %>% 
    group_by(sce, method) %>%
    summarise(coverage = mean(true_value >= ci.lb & true_value <= ci.ub),
              coverage_lb = bootstrap(true_value >= ci.lb & true_value <= ci.ub, function(z) mean(z), quantile, .025),
              coverage_ub = bootstrap(true_value >= ci.lb & true_value <= ci.ub, function(z) mean(z), quantile, .975))
}

summarise_coverage_tau <- function(x, true_value = .15) {
  spread(x[x$statistic %in% c("tau.ci.lb", "tau.ci.ub"), ], statistic, estimate) %>% 
    group_by(sce, method) %>%
    summarise(coverage = mean(true_value >= tau.ci.lb & true_value <= tau.ci.ub),
              coverage_lb = bootstrap(true_value >= tau.ci.lb & true_value <= tau.ci.ub, function(z) mean(z), quantile, .025),
              coverage_ub = bootstrap(true_value >= tau.ci.lb & true_value <= tau.ci.ub, function(z) mean(z), quantile, .975))
}

summarise_ci <- function(x) {
  spread(x[x$statistic %in% c("ci.lb", "ci.ub", "se"), ], statistic, estimate) %>% 
    group_by(sce, method) %>%
    summarise(mean_ci_width = mean(ci.ub - ci.lb),
              median_ci_width = median(ci.ub - ci.lb),
              mean_se = mean(se),
              median_se = median(se))
}

summarise_ci_tau <- function(x) {
  spread(x[x$statistic %in% c("tau.ci.lb", "tau.ci.ub", "tau.se"), ], statistic, estimate) %>% 
    group_by(sce, method) %>%
    summarise(mean_ci_width = mean(tau.ci.ub - tau.ci.lb),
              median_ci_width = median(tau.ci.ub - tau.ci.lb),
              mean_se = mean(tau.se),
              median_se = median(tau.se))
}

summarise_overall_results <- function(x, true_value = 1)
  x %>%
  group_by(statistic, method) %>%
  summarise(mean = mean(estimate),
            SD = sd(estimate),
            mean_SE = SD/sqrt(n()),
            
            median = median(estimate),
            median_lb = bootstrap(estimate, median, quantile, .025),
            median_ub = bootstrap(estimate, median, quantile, .975),
            
            bias = mean(estimate - true_value),
            bias_lb = bootstrap(estimate - true_value, mean, quantile, .025),
            bias_ub = bootstrap(estimate - true_value, mean, quantile, .975),
            
            median_bias = median(estimate - true_value),
            median_bias_lb = bootstrap(estimate - true_value, median, quantile, .025),
            median_bias_ub = bootstrap(estimate - true_value, median, quantile, .975),
            
            RMSE = sqrt(mean((estimate - true_value)^2)),
            RMSE_lb = bootstrap((estimate - true_value)^2, function(z) sqrt(mean(z)), quantile, .025),
            RMSE_ub = bootstrap((estimate - true_value)^2, function(z) sqrt(mean(z)), quantile, .975),
            
            median_RMSE = sqrt(median((estimate - true_value)^2)),
            median_RMSE_lb = bootstrap((estimate - true_value)^2, function(z) sqrt(median(z)), quantile, .025),
            median_RMSE_ub = bootstrap((estimate - true_value)^2, function(z) sqrt(median(z)), quantile, .975),
            
            n = n())

plot_bias_rmse_coverage <- function(df, 
                                    methods = c("adma_full", "adma_gold","adma_xors",
                                                "eq648", "ipd_gold", "ipd_full", "ipd_xors")[c(2,3,4,5,7)],
                                    level = .95,
                                    sce_levels = unique(df$sce),
                                    ...) {
  library(ggplot2)
  
  df <- as.data.frame(df)
  df <- df[df$method %in% methods, ]
  
  bias <- df[ , c("sce", "method", "bias", "bias_lb", "bias_ub")]
  colnames(bias) <- c("Scenario", "Method", "Estimate", "lb", "ub")
  bias$Measure <- "Bias"
  
  rmse <- df[ , c("sce", "method", "RMSE", "RMSE_lb", "RMSE_ub")]
  colnames(rmse) <- c("Scenario", "Method", "Estimate", "lb", "ub")
  rmse$Measure <- "RMSE"
  
  cove <- df[ , c("sce", "method", "coverage", "coverage_lb", "coverage_ub")]
  colnames(cove) <- c("Scenario", "Method", "Estimate", "lb", "ub")
  cove$Measure <- "Coverage"
  
  df <- rbind(bias, rmse, cove)
  df$Scenario <- as.factor(df$Scenario)
  
  dodge <- 1/2
  
  df$Method[df$Method == "eq648"] <- "IPD, Misclassification model"
  df$Method[df$Method == "ipd_gold"] <- "IPD, gold"
  df$Method[df$Method == "ipd_xors"] <- "IPD, naive"
  df$Method[df$Method == "ipd_full"] <- "IPD, reference"
  
  df$Method[df$Method == "adma_gold"] <- "AD, gold"
  df$Method[df$Method == "adma_xors"] <- "AD, naive"
  df$Method[df$Method == "adma_full"] <- "AD, reference"
  
  df$Method <- factor(df$Method, levels = c("IPD, gold",
                                            "AD, gold",
                                            "IPD, naive",
                                            "AD, naive",
                                            "IPD, reference",
                                            "AD, reference",
                                            "IPD, Misclassification model"))
  
  ggplot(data = df,
         aes_string(x = "Scenario", y = "Estimate", group = "Method", ...)) +
    geom_point(aes(color = Method, shape = Method),
               position = position_dodge(dodge)) +
    geom_errorbar(width = 1,
                  aes(ymin = lb, ymax = ub, color = Method),
                  position = position_dodge(dodge)) +
    theme(legend.position = "bottom") +
    facet_wrap( . ~ Measure, scales = "free", ncol = 1) +
    scale_shape_manual(values = c(20, 24, 15, 3, 7, 11, 8)[seq_along(unique(df$Method))]) 
}

