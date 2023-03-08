# https://bookdown.org/rbg/surrogates/chap5.html#chap5library
library(dplyr)
library(tidyr)
library(ggplot2)

library(mvtnorm)
library(plgp)
library(lhs)


# Generating from 1d GP ------------------------------------------------------
n = 100
# Column vector of 100 equally spaced values from 0 to 10
X = matrix(seq(0, 10, length=n), ncol=1)

# Distance matrix
D = plgp::distance(X)

# Construct the covariance matrix as exponentiated negative distance,
# and slightly increase the diagonal to ensure positive definite
eps = sqrt(.Machine$double.eps)
Sigma = exp(-D) + diag(eps, n)

# Generate a multivariate random normal (MVN) with mean zero and covariance Sigma
Y = mvtnorm::rmvnorm(1, sigma=Sigma) %>% t

# Plot the samples
ggplot() +
  geom_line(aes(x=X, y=Y)) +
  ggtitle('Random function under a GP prior')

# Correlation by distance (X)
ggplot() +
  geom_line(aes(x=X, y=exp(-X**2))) +
  ggtitle('Correlation by disatnce')

# Multiple realizations
Y = rmvnorm(3, sigma=Sigma) %>% t
df = data.frame(X=X)
df = cbind(df, data.frame(Y))
names(df) = c('X', 'Y1', 'Y2', 'Y3')
df.long = tidyr::pivot_longer(df, cols=starts_with('Y'))
df.long$name = as.factor(df.long$name)

ggplot(df.long) +
  geom_line(aes(x=X, y=value, color=name)) +
  ggtitle('Multiple realizations')

# Updating a 1d GP --------------------------------------------------
n = 8
X = matrix(seq(0, 2*pi, length=n), ncol=1)
y = sin(X)
D = distance(X)
Sigma = exp(-D) + diag(eps, ncol(D))

# New realizations XX
XX = matrix(seq(-0.5, 2*pi + 0.5, length=100), ncol=1)
DXX = distance(XX)
SXX = exp(-DXX) + diag(eps, ncol(DXX)) # Covariance matrix of XX

# Covariance between X and XX
DX = distance(XX, X)
SX = exp(-DX)

# Compute conditional mean and covariance of X, conditioned on observations XX
mup = SX %*% solve(Sigma) %*% y
Sigmap = SXX - SX %*% solve(Sigma) %*% t(SX)

# Generate observations using posterior
YY = rmvnorm(100, mup, Sigmap)

# Calculate 90% CI for posterior
q1 = mup + qnorm(0.05, mean=0, sd=sqrt(diag(Sigmap)))
q2 = mup + qnorm(0.95, mean=0, sd=sqrt(diag(Sigmap)))

# Prepare observations for plotting
df.XX = data.frame(XX)
df.YY = data.frame(YY %>% t)
names(df.YY) = sapply('Y', paste0, seq(1, ncol(YY)))
df = cbind(df.XX, df.YY)
df.long = pivot_longer(df, !XX)
df.long$name = as.factor(df.long$name)

# Plot observations
ggplot() +
  geom_line(data=df.long, aes(x=XX, y=value, group=name), color='gray', alpha=0.4) +
  geom_point(aes(x=X, y=y), size=5) +              # Original points
  geom_line(aes(x=XX, y=mup)) +                    # Posterior mean
  geom_line(aes(x=XX, y=sin(XX)), color='blue') +  # Posterior actual
  geom_line(aes(x=XX, y=q1), color='red', linetype='dashed') + # Lower CI bound
  geom_line(aes(x=XX, y=q2), color='red', linetype='dashed') + # Upper CI bound
  ggtitle('Posterior predictive distribution')


# Updating a 2d GP --------------------------------------------------------
# Load EPA data
epa.data = readRDS('data/df_data_12list.RDS')
# Single snapshot
epa.df = epa.data[[1]]
epa.df$case_cntl = NULL
epa.df$year = NULL
names(epa.df) = c('station_id', 'x', 'y', 'pm2_5')

# Convert to mean zero and standard deviation 1
epa.df$pm2_5 = epa.df$pm2_5 %>% scale

ggplot(epa.df) +
  geom_point(aes(x=x, y=y, color=pm2_5)) +
  scale_color_gradient2(low='blue', high='red', midpoint=0)

# Sample a small number of stations
epa.df.sample.idx = sample(1:nrow(epa.df), size=25)
epa.df.sample = epa.df[epa.df.sample.idx, ]

# Make a lattice of points across the US for predictions
x1 = seq(min(epa.df$x), max(epa.df$x), length=25)
x2 = seq(min(epa.df$y), max(epa.df$y), length=25)
X.lattice = expand.grid(x1, x2)

# Get X and y values for sample of stations
X.actual = epa.df.sample[, c('x', 'y')] %>% as.matrix
y.actual = epa.df.sample$pm2_5 %>% as.matrix(ncol=1)

# Compute covariance matrices
Sigma.actual = exp(-distance(X.actual)) + diag(eps, nrow(X.actual))
Sigma.lattice = exp(-distance(X.lattice)) + diag(eps, nrow(X.lattice))
Sigma.lattice_actual = exp(-distance(X.lattice, X.actual))

# Compute posterior mean and covariance matrix
mu.posterior = Sigma.lattice_actual %*% chol2inv(chol(Sigma.actual)) %*% y.actual
Sigma.posterior = Sigma.lattice - Sigma.lattice_actual %*% chol2inv(chol(Sigma.actual)) %*% t(Sigma.lattice_actual)
# Force to be symmetric (not symmetric due to rounding)
Sigma.posterior[lower.tri(Sigma.posterior)] = t(Sigma.posterior)[lower.tri(Sigma.posterior)]
# Force to be positive semi-definite
Sigma.posterior = Sigma.posterior + diag(0.01 + min(eigen(Sigma.posterior)$values), ncol(Sigma.posterior))

# Sample from MVN using posterior mean and covariance
YY = rmvnorm(n=100, mean=mu.posterior, sigma=Sigma.posterior)
YY.means = colMeans(YY)

# Prepare data for plotting
lattice.df = data.frame(X.lattice)
names(lattice.df) = c('x', 'y')
lattice.df$mean_pred = YY.means

# Plot actual and predicted
ggplot() +
  geom_point(data=epa.df, aes(x=x, y=y, color=pm2_5), size=3) +   # All measurements
  geom_point(data=lattice.df %>% filter(abs(mean_pred) < 4), aes(x=x, y=y, color=mean_pred), shape=18, size=2) + # Prediction
  scale_color_gradient2(low='blue', high='red', midpoint=0) +
  geom_point(data=epa.df.sample, aes(x=x, y=y), size=5) + # Stations
  ggtitle('Posterior predictive distribution')


# 1d GP regression --------------------------------------------------------
# Generate a noisy sine curve with amplitude > 1
n = 8
X = matrix(seq(0, 2*pi, length=n), ncol=1)
# Double up X
X = rbind(X, X)
n = nrow(X)
y = 5 * sin(X) + rnorm(n, sd=1)
D = distance(X)

# Negative log-likelihood of posterior, assuming scale is known
nlg <- function(g, D, Y) {
  n <- length(Y)
  K <- exp(-D) + diag(g, n)
  Ki <- solve(K)
  ldetK <- determinant(K, logarithm=TRUE)$modulus
  ll <- - (n/2)*log(t(Y) %*% Ki %*% Y) - (1/2)*ldetK
  counter <<- counter + 1
  return(-ll)
}

# Minimize the NLL
counter = 0
g = optimize(nlg, interval=c(eps, var(y)), D=D, Y=y)$minimum

# Compute tau and sigma
# Add the nugget to the covariance function
K = exp(-D) + diag(g, n)
# Compute \hat{tau}^2
tau2hat = (t(y) %*% chol2inv(chol(K)) %*% y) / n
# Convert to a scalar (from a matrix with a single entry)
tau2hat = drop(tau2hat)
tau = sqrt(tau2hat)
sigma = sqrt(tau2hat * g)

# Use these hyperparameters to compute kriging equations
XX = matrix(seq(-0.5, 2*pi + 0.5, length=100), ncol=1)
DXX = distance(XX)
DX = distance(XX, X)
SXX = exp(-DXX) + diag(g, ncol(DXX))
S = exp(-D) + diag(g, ncol(D))
SX = exp(-DX)

mu = SX %*% chol2inv(chol(S)) %*% y
Sigma = tau2hat * (SXX - SX %*% chol2inv(chol(S)) %*% t(SX))
# Force Sigma to be symmeric positive-definite
Sigma[lower.tri(Sigma)] = t(Sigma)[lower.tri(Sigma)]
Sigma = Sigma + diag(0.01 + min(eigen(Sigma)$values), ncol(Sigma))

# Simulate from th posterior
YY = rmvnorm(100, mu, Sigma)

# Calculate 90% CI for posterior
q1 = mu + qnorm(0.05, mean=0, sd=sqrt(diag(Sigma)))
q2 = mu + qnorm(0.95, mean=0, sd=sqrt(diag(Sigma)))

# Prepare observations for plotting
df.XX = data.frame(XX)
df.YY = data.frame(YY %>% t)
names(df.YY) = sapply('Y', paste0, seq(1, ncol(YY)))
df = cbind(df.XX, df.YY)
df.long = pivot_longer(df, !XX)
df.long$name = as.factor(df.long$name)

# Plot observations
ggplot() +
  geom_line(data=df.long, aes(x=XX, y=value, group=name), color='gray', alpha=0.4) +
  geom_point(aes(x=X, y=y), size=5) +              # Original points
  geom_line(aes(x=XX, y=mu)) +                    # Posterior mean
  geom_line(aes(x=XX, y=5*sin(XX)), color='blue') +  # Posterior actual
  geom_line(aes(x=XX, y=q1), color='red', linetype='dashed') + # Lower CI bound
  geom_line(aes(x=XX, y=q2), color='red', linetype='dashed') + # Upper CI bound
  ggtitle('Posterior predictive distribution')


# Lengthscale parameter inference using 2d EPA data -----------------------------------------
# Load EPA data
epa.data = readRDS('data/df_data_12list.RDS')
# Single snapshot
epa.df = epa.data[[1]]
epa.df$case_cntl = NULL
epa.df$year = NULL
names(epa.df) = c('station_id', 'x', 'y', 'pm2_5')

# Convert to mean zero and standard deviation 1
epa.df$pm2_5 = epa.df$pm2_5 %>% scale

# Sample a small number of stations
epa.df.sample.idx = sample(1:nrow(epa.df), size=25)
epa.df.sample = epa.df[epa.df.sample.idx, ]

# Make a lattice of points across the US for predictions
x1 = seq(min(epa.df$x), max(epa.df$x), length=25)
x2 = seq(min(epa.df$y), max(epa.df$y), length=25)
X.lattice = expand.grid(x1, x2)

# Get X and y values for sample of stations
X.actual = epa.df.sample[, c('x', 'y')] %>% as.matrix
y.actual = epa.df.sample$pm2_5 %>% as.matrix(ncol=1)

# Function to minmize the log likelihood, now incorporating the lengthscale and nugget
nl <- function(par, D, Y) {
  theta <- par[1]                                       # lengthscale parameter
  g <- par[2]                                           # nugget parameter
  n <- length(Y)
  K <- exp(-D/theta) + diag(g, n)                      
  Ki <- solve(K)
  ldetK <- determinant(K, logarithm=TRUE)$modulus
  ll <- - (n/2)*log(t(Y) %*% Ki %*% Y) - (1/2)*ldetK
  counter <<- counter + 1
  return(-ll)
}

gradnl <- function(par, D, Y) {
  ## extract parameters
  theta <- par[1]
  g <- par[2]
  
  ## calculate covariance quantities from data and parameters
  n <- length(Y)
  K <- exp(-D/theta) + diag(g, n)
  Ki <- solve(K)
  dotK <- K*D/theta^2
  KiY <- Ki %*% Y
  
  ## theta component
  dlltheta <- (n/2) * t(KiY) %*% dotK %*% KiY / (t(Y) %*% KiY) - 
    (1/2)*sum(diag(Ki %*% dotK))
  
  ## g component
  dllg <- (n/2) * t(KiY) %*% KiY / (t(Y) %*% KiY) - (1/2)*sum(diag(Ki))
  
  ## combine the components into a gradient vector
  return(-c(dlltheta, dllg))
}

# Minimize the log likelihood
D = distance(X.actual)
out = optim(c(0.1, 0.1*var(y.actual)),
            nl, gradnl,
            method="L-BFGS-B",
            lower=eps,
            upper=c(10, var(y.actual)),
            D=D, Y=y.actual)

theta = out$par[1]
g = out$par[2]

# Compute \hat{tau}^2
tau2hat = (t(y) %*% chol2inv(chol(K)) %*% y) / n
# Convert to a scalar (from a matrix with a single entry)
tau2hat = drop(tau2hat)

DXX = distance(X.lattice)
DX = distance(X.lattice, X.actual)

C = exp(-D / theta) + diag(g, ncol(C))
CXX = exp(-DXX / theta) + diag(g, ncol(CXX))
CX = exp(-DX / theta)

mu = CX %*% chol2inv(chol(C)) %*% y.actual
Sigma = tau2hat * (CXX - CX %*% chol2inv(chol(C)) %*% t(CX))

# Sample from the MVN using posterior
YY = rmvnorm(100, mu, Sigma)
YY.means = colMeans(YY)

# Prepare data for plotting
lattice.df = data.frame(X.lattice)
names(lattice.df) = c('x', 'y')
lattice.df$mean_pred = YY.means

# Plot actual and predicted
ggplot() +
  geom_tile(data=lattice.df, aes(x=x, y=y, fill=mean_pred)) +
  geom_point(data=epa.df, aes(x=x, y=y, fill=pm2_5), size=3, pch=21) +
  scale_fill_gradient2(low='blue', high='red', midpoint=0) +
  geom_point(data=epa.df.sample, aes(x=x, y=y), size=5, shape=3) + # Stations
  ggtitle('Posterior predictive distribution')


# Anisotropic modeling ----------------------------------------------------
# Friedman function
fried <- function(n=50, m=6) {
  if(m < 5) stop("must have at least 5 cols")
  X <- randomLHS(n, m)
  Ytrue <- 10*sin(pi*X[,1]*X[,2]) + 20*(X[,3] - 0.5)^2 + 10*X[,4] + 5*X[,5]
  Y <- Ytrue + rnorm(n, 0, 1)
  return(data.frame(X, Y, Ytrue))
}

# Generate data from Friedman function
m = 7 # Number of columns (5 used for Friedman, two irrelevant)
n = 200 # Training set
nprime = 1000 # Testing set
data = fried(n + nprime, m)
X = as.matrix(data[1:n, 1:m])
y = drop(data$Y[1:n])
XX = as.matrix(data[(n+1):(n+nprime), 1:m])
yy = drop(data$Y[(n+1):(n+nprime)])
yytrue = drop(data$Ytrue[(n+1):(n+nprime)])











