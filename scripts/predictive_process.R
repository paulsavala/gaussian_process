library(tidyverse)
library(ggplot2)
library(reshape2)
library(scales)
library(microbenchmark)


# Full covariance matrix --------------------------------------------------
# Parameters
sigma = 1
ell = 0.5

# Data prep - 1d spatial grid
N = 25
x = seq(from=-1.5, to=1.5, length.out=N)

# Construct a full covariance matrix on a grid
gaussian_cov = function(xi, xj, sigma, ell) {
  sigma^2 * exp(-norm(xi - xj, type='2') / (2*ell^2))
}

# Vectorized version
cov_exp_quad = function(x, sigma, ell, y=NULL) {
  if (is.null(y)) {
    y = x
  }
  
  outer(x, y, gaussian_cov, sigma=sigma, ell=ell)
}

# Full covariance matrix
cov.full = cov_exp_quad(x, sigma, ell)

# Predictive process ------------------------------------------------------
make_pp_matrix = function(N.knots, x, sigma, ell, adjust_bias=FALSE) {
  # 1d knot grid
  N.knots = 5
  knots = seq(from=-1.25, to=1.25, length.out=N.knots)
  
  
  # Covariance sub-matrices -------------------------------------------------
  cov.space_knot = cov_exp_quad(x, sigma, ell, knots)
  cov.knot_knot = cov_exp_quad(knots, sigma, ell)
  cov.knot_knot.inv = solve(cov.knot_knot)
  
  # Low-rank covariance matrix
  cov.pp = cov.space_knot %*% cov.knot_knot.inv %*% t(cov.space_knot)
  
  if (adjust_bias) {
    # Bias adjustment (nugget term)
    cov.pp.nugget = diag(cov.full) - diag(cov.pp)
    
    # Bias-adjusted low-rank covariance matrix
    diag(cov.pp) = diag(cov.pp) + cov.pp.nugget
  }
  
  return(cov.pp)
}


# Spectral approximation of covariance matrix -----------------------------
# Square root of approximation function
approx_L = function(x, sigma, ell, M, scale) {
  epsilon = sqrt(1 / (2 * ell^2))
  alpha = 1 / scale
  beta = (1 + (2 * epsilon / alpha)^2)^0.25
  delta = sqrt(alpha^2 * (beta^2 - 1) / 2)
  
  N = length(x)
  Ht = matrix(0, nrow = N, ncol = M)
  xp = alpha * beta * x
  f = sqrt(epsilon^2 / (alpha^2 + delta^2 + epsilon^2))
  Ht[, 1] = sqrt(sqrt(alpha^2 / (alpha^2 + delta^2 + epsilon^2))) * sqrt(beta) * exp(-delta^2 * x * x)
  Ht[, 2] = f * sqrt(2) * xp * Ht[, 1]
  if (M > 2) {
    for(n in 3:M) {
      Ht[, n] = f * sqrt(2 / (n - 1)) * xp * Ht[, n - 1] - f^2 * sqrt((n - 2) / (n - 1)) * Ht[, n - 2]
    }
  }
  
  sigma * Ht
}

# Not giving good approximations to the covariance matrix, not sure why
approx_L.savala = function(M, alpha, sigma, x, z=NULL, sqrt_lambda=FALSE) {
  # If z is not supplied, assume it's centered at zero
  if (is.null(z)) {
    z = 0 * x
  }
  
  # Sanity check on x and z
  stopifnot(length(x) == length(z))
  
  # Calculate parameters used later
  beta = (1 + (2*sigma / alpha)**2)**(0.25)
  delta2 = alpha**2/2 * (beta**2 - 1)
  
  # Phi function
  phi = function(n, t) {
    if (n == 0) {
      return(0 * t) # Multiply by x to keep proper shape
    } else if (n == 1) {
      return(sqrt(beta) * exp(-delta2 * t**2))
    } else {
      return(sqrt(2/(n-1)) * phi(n-1, t) - sqrt((n-2) / (n-1))*phi(n-2, t))
    }
  }
  
  # Lambda function
  lambda_n = function(n, sqrt_lambda=sqrt_lambda) {
    term1 = sqrt(alpha**2 / (alpha**2 + sigma**2 + delta2))
    term2 = (sigma**2 / (alpha**2 + sigma**2 + delta2))**(n-1)
    
    if (sqrt_lambda) {
      term1 = sqrt(term1)
      term2 = sqrt(term2)
    }
    
    return(term1 * term2)
  }
  
  lambda = lambda_n(1:M, sqrt_lambda=sqrt_lambda)
  
  # Number of terms in x
  N = length(x)
  
  # Holds values of phi at all eigenfunctions (rows) and x/z (columns)
  phi_x_mat = matrix(nrow=M, ncol=N)
  phi_z_mat = matrix(nrow=M, ncol=N)
  
  # Generate phi_n(x) * lambda_n and phi_n(z) for all n and x/z
  # Only need to multiply lambda_n times one of them because doing both would
  # give lambda_n^2 in the summand
  for (m in 1:M) {
    phi_x_mat[m,] = phi(m, x) * lambda[m]
    phi_z_mat[m,] = phi(m, z)
  }
  
  # Make the (approximate) covariance matrix
  if (!sqrt_lambda) {
    K = t(phi_x_mat) %*% phi_z_mat
  } else {
    K = t(phi_x_mat) %*% phi_x_mat
  }
  
  return(K)
}

make_spectral_matrix = function(x, sigma, ell, M, scale, adjust_bias=FALSE, diag_cov_full=NULL) {
  cov.spectral.M_scale.sqrt = approx_L(x, sigma, ell, M, scale)
  cov.spectral.M_scale = cov.spectral.M_scale.sqrt %*% t(cov.spectral.M_scale.sqrt)
  
  if (adjust_bias & !is.null(diag_cov_full)) {
    diag(cov.spectral.M_scale) = diag(cov.spectral.M_scale) + (diag_cov_full - diag(cov.spectral.M_scale))
  }
  
  return(cov.spectral.M_scale)
}


# Error evaluation --------------------------------------------------------
# Plot the difference
plot_error_heatmap = function(cov.full, cov.low_rank, low_rank_title, title_l2_norm=FALSE) {
  cov.error = cov.full - cov.low_rank
  cov.error.df = melt(cov.error)
  names(cov.error.df) = c('x', 'y', 'value')
  
  plot_title = paste0("(", low_rank_title, ")")
  if (title_l2_norm) {
    plot_title = paste0(plot_title, "L2 norm = ", round(norm(cov.full - cov.low_rank, type="2"), 3))
  }
  
  g = ggplot(cov.error.df) +
    geom_tile(aes(x=x, y=y, fill=value)) +
    scale_fill_gradient2(low='blue', high='red', midpoint=0) +
    coord_fixed() +
    ggtitle(plot_title) 
  
  show(g)
}

plot_error_heatmap(cov.full, make_spectral_matrix(x, sigma, ell, M=5, scale=0.1), "Spectral (M=5, scale=0.1)")
plot_error_heatmap(cov.full, make_spectral_matrix(x, sigma, ell, M=5, scale=0.25), "Spectral (M=5, scale=0.25)")
plot_error_heatmap(cov.full, make_spectral_matrix(x, sigma, ell, M=5, scale=0.5), "Spectral (M=5, scale=0.5)")
plot_error_heatmap(cov.full, make_spectral_matrix(x, sigma, ell, M=5, scale=0.75), "Spectral (M=5, scale=0.75)")
plot_error_heatmap(cov.full, make_spectral_matrix(x, sigma, ell, M=5, scale=1), "Spectral (M=5, scale=1)")
plot_error_heatmap(cov.full, make_spectral_matrix(x, sigma, ell, M=5, scale=2), "Spectral (M=5, scale=2)")

plot_error_heatmap(cov.full, make_spectral_matrix(x, sigma, ell, M=10, scale=1), "Spectral (M=10, scale=1)")

plot_error_heatmap(cov.full, cov.pp, "Predictive process (5 knots)")

# Show how the error changes for different values of alpha and M
summarize_error = function(cov.full, alpha, M, cov.pp) {
  error_pp = norm(cov.full - cov.pp, type='2')
  
  error_mat = matrix(nrow=length(alpha), ncol=length(M))
  for (a in alpha) {
    a.idx = match(a, alpha)
    for (m in M) {
      m.idx = match(m, M)
      cov.low_rank = make_spectral_matrix(x, sigma, ell, M=m, scale=a)
      error_mat[a.idx, m.idx] = norm(cov.full - cov.low_rank, type='2') / error_pp
    }
  }
  
  # Convert to a long dataframe for plotting
  error_mat.df = melt(error_mat)
  names(error_mat.df) = c('scale', 'M', 'error')
  error_mat.df$log_error = log(error_mat.df$error)
  
  # Plot error
  ggplot(error_mat.df) +
    geom_tile(aes(x=M, y=scale, fill=log_error)) +
    scale_fill_gradient2(low='blue', mid='white', high='red', midpoint=0) +
    ggtitle('Comparison between error of predictive process and spectral covariance approximations')
}

summarize_error(cov.full, 
                alpha=seq(from=0.1, to=1, by=0.1),
                M=seq(from=2, to=10, by=1),
                cov.pp=cov.pp)


# Comparing timings -------------------------------------------------------







