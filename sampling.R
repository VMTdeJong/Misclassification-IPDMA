#' Sample clustered binary data with misclassified exposure
#' 
#' Sample clustered binary data with misclassified exposure, correctly classified exposure, covariates
#' and outcome. Random intercepts for outcome, measurement and exposure model. Random effects for
#' outcome model.
#' 
#' @param n sample size per cluster
#' @param J number of clusters
#' @param b.x coefficient for x in the outcome model
#' @param b.x.tau variance of random effect for x in the outcome model
#' @param b.z log odds ratios for covariates.
#' @param a intercept in the outcome model. if NULL, the incidence is fixed so that mean(y) = 1/2
#' @param a.tau variance of the random intercept in the outcome model
#' @param z.mean.tau SD of the mean of the z covariates across studies
#' @param gamma.z coefficients for z in the exposure model
#' @param gamma.0 intercept in the exposure model
#' @param gamma.0.tau variance of the random intercept in the exposure model
#' @param lambda.0.0 intercept in the measurement model, if x = 1
#' @param lambda.0.j.tau variance of the reandom intercept in the measurement model, if x = 1
#' @param lambda.z coefficients in the measurement model, if x = 1
#' @param phi.0.0 intercept in the measurement model, if x = 0
#' @param phi.0.j.tau variance of the random intercept in the measurement model, if x = 0
#' @param phi.z coefficients in the measurement model, if x = 0
#' @param col.names optional column names
#' @param center.z should covariates be centered? (post hoc)
#' 
#' @return A data.frame, including the mismeasured and correctly measured variables.

sample_misclass_re <- function(n = 50, J = 1,
                               a = NULL, a.tau = 0,
                               b.x = log(2), b.x.tau = 0, b.z = log(2), 
                               z.mean.tau = 1/2, cov = 1/4,
                               gamma.0 = 0, gamma.0.tau = 0, gamma.z = 1/2, 
                               lambda.0.0 = 3, lambda.0.j.tau = 1, lambda.z = -2,
                               phi.0.0 = -3, phi.0.j.tau = 1, phi.z = 2,
                               col.names = NULL, center.z = TRUE ) {
  ### Dimensions
  N <- J * n
  if (length(gamma.z) != length(lambda.z) || length(gamma.z) != length(phi.z) || length(gamma.z) != length(phi.z))
    stop("The number of coefficients should be the same for each submodel.")
  
  # Sample parameters for exposure model
  z_mean.j <- rnorm(J, mean = 0, sd = z.mean.tau)
  gamma.0.j <- rnorm(J, mean = gamma.0, sd = gamma.0.tau)
  
  # Sample exposure and covariate
  dat <- list()
  
  for (Jj in seq_len(J))
    dat[[Jj]] <- sample_x_z_cont_to_bin(n = n, gamma.0 = gamma.0.j[Jj], gamma.z = gamma.z, z_mean = z_mean.j[Jj], 
                                        cov = cov, j = Jj, outcome = NA)
  
  dat <- Reduce(rbind, dat)
  
  x <- dat[ , "x"]
  z <- dat[ , grep("z", colnames(dat))]
  j <- dat[ , "j"]
  
  ### Outcome model
  # Optionally set mean(y) = .5
  if (is.null(a)) 
    a <- -log(sqrt(prod(exp(c(b.x))))) 
  
  # Generate random intercepts and coefficients.
  a.j <- stats::rnorm(J, a, a.tau)[j]
  b.x.j <- stats::rnorm(J, b.x, b.x.tau)[j]
  coefficients <- cbind(matrix(a.j, nrow = N, ncol = 1), 
                        matrix(b.x.j, nrow = N, ncol = 1), 
                        matrix(rep(b.z, each = N), ncol = length(b.z), nrow = N))
  
  # Draw outcome
  y.p <- inv_logit(rowSums(coefficients * cbind(1, x, z)))
  y   <- stats::rbinom(N, size = 1, prob = y.p)
  
  ### Measurement model
  s <- rep(0, N)
  
  lambda.0.j <- stats::rnorm(J, lambda.0.0, lambda.0.j.tau)[j]
  s.x1.p <- inv_logit(lambda.0.j + lambda.z * z) 
  s[x == 1] <- stats::rbinom(N, 1, s.x1.p)[x == 1]
  
  phi.0.j <- stats::rnorm(J, phi.0.0, phi.0.j.tau)[j]
  s.x0.p <- inv_logit(phi.0.j + phi.z * z) 
  s[x == 0] <- stats::rbinom(N, 1, s.x0.p)[x == 0]
  
  # Return everything except the intercept
  out <- data.frame(y, j, s, x, z)
  if (!is.null(col.names))
    colnames(out) <- col.names
   
  if (center.z)
    out <- metamisc:::centerCovs(out, y.name = c("y", "x", "s"), cluster.name = "j")
  out
}

#'@usage
# sample_misclass_re(gamma.z = c(1/2, 1/2), lambda.z = c(-1, -1), phi.z = c(1, 1), b.z= c(1, 1))


inv_logit <- function(x) 
  1 / (1 + exp(-x))

logit <- function(x)
  log(x/(1-x))

sample_x_z_cont_to_bin <- function(n, gamma.0, gamma.z, z_mean, cov, j = NA, outcome = NA) {
  z <- sample_z_cont(n = n, z_n = length(gamma.z), z_mean = z_mean, cov = cov)
  x_p <- inv_logit(gamma.0 + z %*% gamma.z)
  x <- rbinom(n, 1, x_p)  
  cbind(x = x, z = z, j = rep(j, n), outcome = rep(outcome, n))
}

sample_z_cont <- function(n, z_n, z_mean, cov) {
  sigma <- matrix(cov, ncol = z_n, nrow = z_n)
  diag(sigma) <- 1
  out <- mvtnorm::rmvnorm(n = n, mean = rep(z_mean, z_n), sigma = sigma)
  colnames(out) <- paste("z", seq(from = 1, to = z_n), sep = "")
  out
}


#' Remove data from object
#' 
#' @param data data set.
#' @param cluster name of cluster variable
#' @param j vector containing cluster id's for each observation
#' @param columns Names or indices of columns for which data is to be removed
make_NA <- function(data, cluster = "j", j, columns, ...) {
  data[data[, cluster] %in% j, columns] <- NA
  data
}

divide_main_val_studies <- function(data) {
  j <- sort(unique(data$j))
  j_split <- split(j, sort(j%%2))
  data <- make_NA(data, j_split[[1]], "x")
  data <- make_NA(data, j_split[[2]], "y")
  data
}

add_missing <- function(data, J_gold = 5) {
  data$x_miss <- data$x
  data$x_miss[data$j > J_gold] <- NA
  data
}

#' Replace x with s, where necessary
#' @param data data.frame
#' @param x partially missing variable
#' @param s fully observed variable
xors <- function(data, x, s) {
  data$xors <- data[ , x]
  data$xors[is.na(data$xors)] <- data[is.na(data$xors), s]
  data
}

xors_replace <- function(data, x, s) {
  data[is.na(data[ , x]), x] <- data[is.na(data[ , x]), s]
  data
}
