# Z(x) = Z_0(x) + sqrt(alpha) * sum[K(x-r_i)*Z_i(x)]
# Z_i = stationary GP with Matern covariance
# alpha = real, positive
# r_i = knot
# K = Epanechnikov kernel

# Parameters
# Matern:
# - tau2 = diagonal variance (nugget)
# - sigma = scale
# - rho = lengthscale
# - eta = variance (smoothness)
#
# Kernel:
# - h = bandwidth parameter

library(rSPDE)
library(fields)
library(ggforce)
library(dplyr)

# Data loading/prep -------------------------------------------------------
# Data parameters
N.KNOTS = 10
RADIUS = 0.2
# N.LATTICE.X = 25
# N.LATTICE.Y = 25

# Load EPA data
epa.data = readRDS('data/df_data_12list.RDS')
# Single snapshot
epa.df = epa.data[[1]]
epa.df$case_cntl = NULL
epa.df$year = NULL
names(epa.df) = c('station_id', 'x', 'y', 'pm2_5')
epa.df$station_id = NULL

# Filter to remove duplicates
epa.df = epa.df[!duplicated(epa.df[, c('x', 'y')]), ]

# # Filter to just western US
# epa.df = epa.df %>%
#   filter(x < 0.25)

# Convert to mean zero and standard deviation 1
epa.df$pm2_5 = epa.df$pm2_5 %>% scale

# Train-test split data
train_size = as.integer(0.7 * nrow(epa.df))
epa.train.idx = sample(1:nrow(epa.df), size=train_size)
epa.test.idx = setdiff(1:nrow(epa.df), epa.train.idx)
epa.train.df = epa.df[epa.train.idx, ]
epa.test.df = epa.df[epa.test.idx, ]

# Create X and y matrices
X.train = epa.train.df[, names(epa.train.df) != 'pm2_5']
X.test = epa.test.df[, names(epa.test.df) != 'pm2_5']
y.train = epa.train.df$pm2_5
y.test = epa.test.df$pm2_5

# Package into training and testing data frames
df.train = cbind(X.train, y.train)
df.test = cbind(X.test, y.test)
names(df.train) = c('x', 'y', 'pm2_5')
names(df.test) = c('x', 'y', 'pm2_5')


# Helper functions --------------------------------------------------------
# Epanechnikov kernel (distance function)
epan.kernel = function(u, h) {
  (2 / pi) * (1 / h**2) * (1-(u/h)**2)
}

generate.gp = function(X.knots, X, y, only.background=FALSE) {
  N.KNOTS = nrow(X.knots)
  RADIUS = 0.2
  
  # Empty list to hold all local spatial processes
  gp.list = list()
  
  if (!only.background) {
    for (knot.id in 1:N.KNOTS) {
      # Get knot coordinates
      r = X.knots[knot.id, ] %>% matrix(ncol=ncol(X.knots))
      
      # Find distance from all points to current knot
      D = rdist(r, X)
      
      # Get points within circle of radius RADIUS
      r.idx = which(D < RADIUS)
      r.nbhd = X[r.idx, ]
      r.nbhd.y = y[r.idx]
      
      # Fit spatial processes
      gp.fit = spatialProcess(x=r.nbhd, y=r.nbhd.y)
      
      # Store fitted spatial process
      gp.list[[knot.id]] = gp.fit
    }
  }
  
  # Fit background GP -------------------------------------------------------
  gp.background = spatialProcess(x=X, y=y)
  
  return(list("gp.background"=gp.background,
              "gp.list"=gp.list))
}

# Generate simulations from all spatial processes (local + background)
simulate.gp = function(gp.list,
                       gp.background,
                       X) {
  
  # Number of samples (used for estimating variability/calculating entropy)
  num_samples = 100
  
  # Simulate from background spatial process
  sim.background = sim.spatialProcess(gp.background, X, M=num_samples)
  
  # Simulate from local spatial processes
  sim.local = list()
  if (!is.null(gp.list)) {
    for (i in 1:length(gp.list)) {
      sim.gp = sim.spatialProcess(gp.list[[i]], X, M=num_samples)
      sim.local[[i]] = sim.gp
    }
  }
  
  return(list("sim.background"=sim.background,
              "sim.local"=sim.local))
}

gp.preds = function(sim.background, sim.local.list, X.knots, X.obs, y.obs, sqrt_alpha, ell) {
  # Dimensions:
  # sim.background: (num locations) x (num samples generated)
  # sim.local: same
  # epan.kernel(X.knots, X.obs): (num knots) x (num locations)
  
  # Holds predictions for each sample
  sim.total = matrix(NA, nrow=nrow(sim.background), ncol=ncol(sim.background))
  
  # Go sample-by-sample
  for (j in 1:ncol(sim.background)) {
    # Just the j'th sample from background spatial process
    sim.background.j = sim.background[, j]
    # Empty vector to hold the sum of all local spatial processes
    sim.local.sum.j = 0 * sim.background.j
    # Loop over spatial processes
    for (i in 1:length(sim.local.list)) {
      # Just the j'th sample of the i'th local spatial process
      sim.local.i.j = sim.local.list[[i]][, j]
      sim.local.sum.j = sim.local.sum.j + 
        epan.kernel(u=rdist(matrix(X.knots[i, ], ncol=2), X.obs), h=ell) * sim.local.i.j
    }
    # Combine local and background predictions into one row (all locations simultaneously)
    sim.total[, j] = sim.background.j + sqrt_alpha * sim.local.sum.j
  }
  
  return(sim.total)
}

# Function passed to optim to fit sqrt_alpha and ell
fit.params = function(sim.background, sim.local.list, X.knots, X.train, y.train, par) {
  # Extract parameters of interest
  sqrt_alpha = par[1]
  ell = par[2]
  
  sim.total = gp.preds(sim.background=sim.background,
                       sim.local.list=sim.local.list,
                       X.knots=X.knots,
                       X.obs=X.train,
                       y.obs=y.train,
                       sqrt_alpha=sqrt_alpha,
                       ell=ell)
  
  # Replicate actual values along each column to compute MSE
  y.train.matrix = matrix(y.train, nrow=length(y.train), ncol=ncol(sim.background))
  
  # Compute MSE 
  mse = sum((sim.total - y.train.matrix)**2 / length(y.train.matrix))
  return(mse)
}

# Fitting local GPs using training data ------------------------------------------------------------
fit.fuentes.gp = function(X.knots, X.train, y.train) {
  N.KNOTS = nrow(X.knots)
  
  gp.all = generate.gp(X.knots, X.train, y.train)
  gp.background = gp.all$gp.background
  gp.list = gp.all$gp.list
  
  # Generate predictions from all processes
  sim.train.values = simulate.gp(gp.list, gp.background, X.train)
  # Save background and local processes separately for easy reference
  sim.train.background = sim.train.values$sim.background
  sim.train.local.list = sim.train.values$sim.local
  
  params = optim(par=c(1, 1), fn=fit.params, 
                 method="L-BFGS",
                 lower=c(0, 0),
                 upper=c(10, 1),
                 sim.background=sim.train.background,
                 sim.local.list=sim.train.local.list,
                 X.knots=X.knots,
                 X.train=X.train,
                 y.train=y.train)
  
  sqrt_alpha = params$par[1]
  ell = params$par[2]
  
  return(list("sqrt_alpha"=sqrt_alpha, 
              "ell"=ell))
}

predict.fuentes.gp = function(X.knots, X, y, params) {
  gp.all = generate.gp(X.knots, X, y)
  gp.background = gp.all$gp.background
  gp.list = gp.all$gp.list
  
  # Generate predictions from all processes
  sim.values = simulate.gp(gp.list, gp.background, X)
  # Save background and local processes separately for easy reference
  sim.background = sim.values$sim.background
  sim.local.list = sim.values$sim.local
  
  best_fit_preds = gp.preds(sim.background=sim.background,
                            sim.local.list=sim.local.list,
                            X.knots=X.knots,
                            X.obs=X,
                            y.obs=y,
                            sqrt_alpha=params$sqrt_alpha,
                            ell=params$ell)
  
  best_fit_preds.mean = apply(best_fit_preds, 1, mean)
  best_fit_preds.median = apply(best_fit_preds, 1, median)
  # Bootstrap standard error
  best_fit_preds.se = apply(best_fit_preds, 1, sd) / sqrt(ncol(sim.background))
  
  best_fit_preds.df = data.frame(mean_pred=best_fit_preds.mean,
                                 median_pred=best_fit_preds.median,
                                 se_pred=best_fit_preds.se)
  return(best_fit_preds.df)
}


# Fit Fuentes GP on a number of different knots ---------------------------
N_KNOT_SETS = 25

knot_set.list.X = list()
knot_set.list.y = list()
knot_set.info.df = data.frame(id=NA, se=NA, utility=NA)

for (knot_set.idx in 1:N_KNOT_SETS) {
  print(paste0("Fitting knot set ", knot_set.idx, "/", N_KNOT_SETS, "..."))
  
  # Create a lattice of x and y coordinates from which to sample
  lattice.x = seq(from=min(epa.df$x), to=max(epa.df$x), length.out=100)
  lattice.y = seq(from=min(epa.df$y), to=max(epa.df$y), length.out=100)
  
  # Sample N.KNOTS random x and y locations
  X.knots = data.frame(x=sample(lattice.x, size=N.KNOTS),
                       y=sample(lattice.y, size=N.KNOTS))
  
  # Fit Fuentes model hyperparameters
  params = fit.fuentes.gp(X.knots, X.train, y.train)
  preds.df = predict.fuentes.gp(X.knots, X.test, y.test, params)
  
  # === Calculate utility ===
  # Calculate median PM 2.5 in neighborhood of knot
  knot_set.median_pm2_5 = c()
  for (i in 1:N.KNOTS) {
    # Get knot coordinates
    r = X.knots[i, ] %>% matrix(ncol=ncol(X.knots))
    
    # Find distance from all points to current knot
    D = rdist(r, X.train)
    
    # Get points within circle of radius RADIUS
    r.idx = which(D < RADIUS)
    r.nbhd.y = y.train[r.idx]
    
    # Calculate median PM 2.5 in this radius
    r.pm2_5 = median(r.nbhd.y)
    
    # Record this median
    knot_set.median_pm2_5 = c(knot_set.median_pm2_5, r.pm2_5)
  }
  
  # Calculate background GP using training data
  gp.background = generate.gp(X.knots, X.train, y.train, only.background=TRUE)$gp.background
  
  # Predict at knots using background GP
  knot_set.background.preds = simulate.gp(gp.list=NULL, gp.background, X.knots)$sim.background
  
  # Store median PM 2.5 and background preds together
  knot_set.df = data.frame(median_pm2_5=knot_set.median_pm2_5,
                           median_preds=apply(knot_set.background.preds, 1, median))
  
  # Record difference from median in z-score
  knot_set.df$pm2_5_z = scale(knot_set.df$median_preds - knot_set.df$median_pm2_5)
  
  # Record knots in knot_set.list
  knot_set.list.X[[knot_set.idx]] = X.knots
  
  # Record knot set info in knot_set.info.df
  knot_set.info.df[knot_set.idx, 'id'] = knot_set.idx
  knot_set.info.df[knot_set.idx, 'se'] = sum(preds.df$se_pred)
  knot_set.info.df[knot_set.idx, 'utility'] = sum(knot_set.df$pm2_5_z)
}

# Find the knots using the entropy (SE) and utility defined in her paper
max_se.idx = knot_set.info.df[knot_set.info.df$se == max(knot_set.info.df$se), 'id']
max_utility.idx = knot_set.info.df[knot_set.info.df$utility == max(knot_set.info.df$utility), 'id']

if (max_se.idx == max_utility.idx) {
  best_knots.idx = max_se.idx
} else {
  se_diff = knot_set.info.df[max_se.idx, 'se'] - knot_set.info.df[max_utility.idx, 'se']
  se_sum = knot_set.info.df[max_se.idx, 'se'] + knot_set.info.df[max_utility.idx, 'se']
  se_relative_gain = se_diff / se_sum
  
  utility_diff = knot_set.info.df[max_utility.idx, 'utility'] - knot_set.info.df[max_se.idx, 'utility']
  utility_sum = knot_set.info.df[max_utility.idx, 'utility'] + knot_set.info.df[max_se.idx, 'utility']
  utility_relative_gain = utility_diff / utility_sum
  
  if (se_relative_gain > utility_relative_gain) {
    best_knots.idx = max_se.idx
  } else {
    best_knots.idx = max_utility.idx
  }
}

# Select knots for this index
X.knots = knot_set.list.X[[best_knots.idx]]

# Fit model using these knots
params = fit.fuentes.gp(X.knots, X.train, y.train)

# Make predictions using this model
preds.df = predict.fuentes.gp(X.knots, X.test, y.test, params)












