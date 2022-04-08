source("utils_summarize.R")
load("aggregated/tidy estimates.RData")
# load("aggregated dt/tidy estimates.RData")

# Split into beta and tau
est_beta <- est_all[est_all$statistic %in% c("mean","se", "ci.lb", "median", "ci.ub"), ]
est_tau <- est_all[est_all$statistic %in% c("tau.mean", "tau.se", "tau.ci.lb", "tau.median", "tau.ci.ub"), ]

# average sims succeeded per sce:
nrow(est_all) / length(unique(est_all$sce)) / length(unique(est_all$method)) / length(unique(est_all$statistic))

# Beta
beta_summary <- summarise_results(est_beta)
beta_summary_overall <- summarise_overall_results(est_beta)

beta_coverage <- summarise_coverage(est_beta)
beta_ci <- summarise_ci(est_beta)

print(beta_summary[beta_summary$statistic == "median", c("sce", "method", "mean", "RMSE", "bias", "bias_lb", "bias_ub", "n")], n = 100)
print(beta_summary_overall[beta_summary_overall$statistic == "median", c("method", "mean", "RMSE", "bias", "bias_lb", "bias_ub")], n = 100)
print(beta_coverage, n = 100)

beta_median <- beta_summary[beta_summary$statistic == "median", ]


if (all(beta_median$sce == beta_coverage$sce) && all(beta_median$method == beta_coverage$method))
  beta_to_plot <- cbind(beta_median, beta_coverage[ , 3:5])

pdf("summarized/plot_beta.pdf")
plot_bias_rmse_coverage(beta_to_plot) 
dev.off()

pdf("summarized/plot_beta_sce1.pdf")
plot_bias_rmse_coverage(beta_to_plot[beta_to_plot$sce == 1, ], ncol = 3) 
dev.off()

# Tau
tau_summary <- summarise_results(est_tau, true_value = .15)
tau_summary_overall <- summarise_overall_results(est_tau, true_value = .15)
tau_coverage <- summarise_coverage_tau(est_tau, true_value = .15)
tau_ci <- summarise_ci_tau(est_tau)

tau_median <- tau_summary[tau_summary$statistic == "tau.median", ]

if (all(tau_median$sce == tau_coverage$sce) && all(tau_median$method == tau_coverage$method) )
  tau_to_plot <- cbind(tau_median, tau_coverage[ , 3:5])

pdf("summarized/plot_tau.pdf")
plot_bias_rmse_coverage(tau_to_plot) 
dev.off()

pdf("summarized/plot_tau.pdf")
plot_bias_rmse_coverage(tau_to_plot[tau_to_plot$sce == 1, ], ncol = 3) 
dev.off()

print(tau_summary[tau_summary$statistic == "tau.median", c("sce", "method", "mean", "RMSE", "bias", "bias_lb", "bias_ub", "n")], n = 100)
print(tau_summary_overall[tau_summary_overall$statistic == "tau.median", c("method", "mean", "RMSE", "bias", "bias_lb", "bias_ub")], n = 100)
print(tau_coverage, n = 100)

# save(beta_summary,         file = "summarized/beta_summary.RData")
# save(beta_summary_overall, file = "summarized/beta_summary_overall.RData")
# save(beta_coverage, file = "summarized/beta_coverage.RData")
# save(beta_ci, file = "summarized/tau_ci.RData")
# save(beta_median, file = "summarized/beta_median.RData")
# save(beta_to_plot, file = "summarized/beta_to_plot.RData")
# 
# save(tau_summary,          file = "summarized/tau_summary.RData")
# save(tau_summary_overall,  file = "summarized/tau_summary_overall.RData")
# save(tau_coverage, file = "summarized/tau_coverage.RData")
# save(tau_ci, file = "summarized/tau_ci.RData")
# save(tau_median, file = "summarized/tau_median.RData")
# save(tau_to_plot, file = "summarized/tau_to_plot.RData")
