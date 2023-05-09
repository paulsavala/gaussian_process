library(tidyverse)
library(rstan)
library(ggplot2)
library(reshape2)


# Paper 1 -----------------------------------------------------------------
# 1-d ---------------------------------------------------------------------
# Compute approximate decomposition
# HE'S DOING THE SQUARE ROOT OF THE COVARIANCE MATRIX, HENCE THE EXTRA SQUARE
# ROOT ON THOSE TERMS!!!
approx_L.1d = function(M, scale, x, sigma, l) {
  epsilon = 1 / sqrt(l) # (Multiplicative) length scale factor
  alpha = 1 / scale
  beta = (1 + (2 * epsilon / alpha)^2)^0.25 # Equation 3.4
  delta = sqrt(alpha^2 * (beta^2 - 1) / 2) # Equation 3.4
  
  N = length(x)
  Ht = matrix(0, nrow = N, ncol = M)
  
  # Term used inside Hermite polynomials
  xp = alpha * beta * x
  
  # Terms for the (square root of the) eigenvalue. Square root is because
  # we're working with the square root of the covariance matrix, so the eigenvalues
  # are also square rooted (eigenfunctions _do not_ change).
  lambda.term1 = sqrt(sqrt(alpha^2 / (alpha^2 + delta^2 + epsilon^2)))
  # Term that gets the power n-1 inside the eigenvalues
  lambda.term2 = sqrt(epsilon^2 / (alpha^2 + delta^2 + epsilon^2))
  
  # The eigenvalues satisfy lambda_n = lambda_{n-1} * lambda.term2 (easy to check this)
  
  # When n=1, lambda.term2 = 1
  Ht[, 1] = lambda.term1 * sqrt(beta) * exp(-delta^2 * x * x)
  # Use the recurrence relation on Hermite polynomials to get this second one
  Ht[, 2] = lambda.term2 * sqrt(2) * xp * Ht[, 1]
  # Use the recurrence relation on Hermite polynomials to get the rest
  for(n in 3:M) {
    Ht[, n] = lambda.term2 * sqrt(2 / (n - 1)) * xp * Ht[, n - 1] - lambda.term2^2 * sqrt((n - 2) / (n - 1)) * Ht[, n - 2]
  }
  
  sigma * Ht
}

approx_L = function(M, scale, x, sigma, l) {
  epsilon = sqrt(1 / (2 * l^2))
  alpha = 1 / scale
  beta = (1 + (2 * epsilon / alpha)^2)^0.25
  delta = sqrt(alpha^2 * (beta^2 - 1) / 2)
  
  N = length(x)
  Ht = matrix(0, nrow = N, ncol = M)
  xp = alpha * beta * x
  f = sqrt(epsilon^2 / (alpha^2 + delta^2 + epsilon^2))
  Ht[, 1] = sqrt(sqrt(alpha^2 / (alpha^2 + delta^2 + epsilon^2))) * sqrt(beta) * exp(-delta^2 * x * x)
  Ht[, 2] = f * sqrt(2) * xp * Ht[, 1]
  for(n in 3:M) {
    Ht[, n] = f * sqrt(2 / (n - 1)) * xp * Ht[, n - 1] - f^2 * sqrt((n - 2) / (n - 1)) * Ht[, n - 2]
  }
  
  sigma * Ht
}

# Compute exact covariance matrix
cov_exp_quad = function(x, sigma, l) {
  outer(x, x, function(xi, xj) { sigma^2 * exp(-(xi - xj)^2 / (2*l^2)) })
}

approx_L.1d.savala = function(M, alpha, sigma, x, z=NULL) {
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

# Evaluate norm of error for different values of alpha and M
num.M = 15
num.alpha = 20
norm_diff_mat = matrix(nrow=num.M, ncol=num.alpha)

N = 50
x = seq(-1.5, 1.5, length = N)
l = 0.25
sigma = 1.0

a.idx = 1
a.values = seq(0.01, 1, length.out=num.alpha)
a.values.df = data.frame(a.idx=1:length(a.values), a=a.values)
for (a in a.values) {
  for (m in 1:num.M) {
    K.actual = cov_exp_quad(x, sigma=1, l=l)
    K.spectral = approx_L(m, a, x, sigma=1, l=l)
    
    norm_diff_mat[m, a.idx] = norm(K.spectral %*% t(K.spectral) - K.actual, type="2")
  }
  a.idx = a.idx + 1
}

norm_diff_mat.df = melt(norm_diff_mat)
names(norm_diff_mat.df) = c('M', 'a.idx', 'norm')
norm_diff_mat.df = merge(norm_diff_mat.df, a.values.df, by='a.idx')
norm_diff_mat.df$a.idx = NULL

ggplot(norm_diff_mat.df) +
  geom_tile(aes(x=M, y=a, fill=norm)) +
  scale_fill_gradient2(low='blue', high='black') +
  xlab('M') + ylab('alpha')


# === Basketball code ===
# Or maybe just compare matrices
L = approx_L.1d(M, 0.25, x, sigma, l)
LLT = L %*% t(L)
dimnames(LLT) = list(x, x)
K = cov_exp_quad(x, sigma, l)
dimnames(K) = list(x, x)
bind_rows(LLT %>% melt %>% mutate(which = "approximate"),
          K %>% melt %>% mutate(which = "exact")) %>%
  rename(x1 = Var1, x2 = Var2) %>%
  ggplot(aes(x1, x2, fill = value)) +
  geom_tile() +
  facet_grid(. ~ which) +
  coord_equal()

# Or maybe just plot some eigenfunctions
L = approx_L.1d(M, 1, x, sigma, l)
bind_cols(as_tibble(L), as_tibble(x)) %>% rename(x = value) %>%
  gather(basis, value, 1:M) %>% ggplot(aes(x, value)) +
  geom_line(aes(group = basis, color = basis)) +
  ggtitle("Eigenfunctions (scaled by sqrt(eigenvalues))")


# 2-d ---------------------------------------------------------------------
# Compute approximate decomposition (use same parameters in both directions for simplicity)
approx_L.2d = function(M, scale, loc.df, sigma, l) {
  epsilon = sqrt(1 / (2 * l^2)) # Chosen freely (pg. 5), , but is probably just the scale parameter (different form for notational reasons?)
  alpha = 1 / scale # Chosen freely (pg. 5)
  beta = (1 + (2 * epsilon / alpha)^2)^0.25 # Equation 3.4
  delta = sqrt(alpha^2 * (beta^2 - 1) / 2) # Equation 3.4
  
  # Number of points
  N = nrow(loc.df) 
  # Holds evaluation of expansion at each point in data for each eigenfunction used in truncation
  Ht = matrix(0, nrow = N, ncol = M)
  # Used inside Hermite polynomials
  xp = alpha * beta * x
  # Part of equation 3.5b (pre-computed for simplicity)
  f = sqrt(alpha^2 / (alpha^2 + delta^2 + epsilon^2))
  
  # Compute summands
  # phi_n(x) = sqrt(beta / (2^(n-1) * (n-1)!)) * exp(-delta^2 * x^2) * H_n(xp)
  # lambda_n = sqrt(alpha^2 / (alpha^2 + delta^2 + epsilon^2)) * (epsilon^2 / (alpha^2 + delta^2 + epsilon^2))^(n-1)
  
  # = First summand = 
  # phi_1(x) = sqrt(beta) * exp(-delta^2 * x^2)
  # lamda_1 = f
  # Not sure where he gets the extra square root from. But empirically the fit _is_ much better with it!
  Ht[, 1] = sqrt(sqrt(alpha^2 / (alpha^2 + delta^2 + epsilon^2))) * sqrt(beta) * exp(-delta^2 * x * x)
  
  # = Second summand = 
  # phi_2(x) = sqrt(beta / 2) * exp(-delta^2 x^2) * H_2(xp)
  # Hermite polynomials recurrence relation: H_{n+1}(x) = x * H_n(x) - n*H_{n-1}(x)
  # phi_2(x) = sqrt(beta / 2) * (xp * H_1(xp))
  # lambda_2 = f * (epsilon^2 / (alpha^2 + delta^2 + epsilon^2))
  Ht[, 2] = f * sqrt(2) * xp * Ht[, 1]
  for(n in 3:M) {
    Ht[, n] = f * sqrt(2 / (n - 1)) * xp * Ht[, n - 1] - f^2 * sqrt((n - 2) / (n - 1)) * Ht[, n - 2]
  }
  
  return(sigma * Ht)
}


# Updated approach --------------------------------------------------------
alpha = 2
epsilon = 0.1
beta = (1 + 4*epsilon**2 / alpha**2)**(0.25)
delta2 = alpha**2 / 2 * (beta**2 - 1)

phi = function(n, x) {
  if (n == 0) {
    return(0 * x) # Multiply by x to keep proper shape
  } else if (n == 1) {
    return(sqrt(beta) * exp(-delta2 * x**2))
  } else {
    return(sqrt(2/n) * phi(n-1, x) - sqrt((n-1) / n)*phi(n-2, x))
  }
}

x = seq(from=-1.5, to=1.5, length.out=10)
N = length(x)
n = seq(from=1, to=N)

phi1 = matrix(NA, nrow=N, ncol=N)
for (i in 1:N) {
  for (j in 1:N) {
    phi1[i, j] = phi(n[j], x[i])
  }
}


