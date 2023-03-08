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

# Sample a small number of stations as knots
epa.sample.idx = sample(1:nrow(epa.train.df), size=N.KNOTS)
epa.sample.df = epa.train.df[epa.sample.idx, ]

# # Make a lattice of points across the US for predictions
# x1 = seq(min(epa.df$x), max(epa.df$x), length=N.LATTICE.X)
# x2 = seq(min(epa.df$y), max(epa.df$y), length=N.LATTICE.X)
# X.lattice = expand.grid(x1, x2)
# df.lattice = data.frame(X.lattice)
# names(df.lattice) = c('x', 'y')

# Get X and y values for sample of stations
X.knots = epa.sample.df[, c('x', 'y')] %>% as.matrix
y.knots = epa.sample.df$pm2_5 %>% as.matrix(ncol=1)
df.knots = data.frame(X=X.knots, y=y.knots)
names(df.knots) = c('x', 'y', 'pm2_5')
rownames(df.knots) = 1:nrow(df.knots)

# Plot locations
ggplot(epa.df) +
  geom_point(aes(x=x, y=y, color=pm2_5)) +
  scale_color_gradient2(low='blue', high='red', midpoint=0) +
  geom_point(data=epa.sample.df, aes(x=x, y=y), shape=4, size=5) +
  geom_point(data=df.knots, aes(x=x, y=y), size=5, shape=10) +
  geom_circle(data=df.knots, aes(x0=x, y0=y, r=RADIUS)) +
  coord_fixed()

# ggplot(df.knots, aes(x=x, y=y, label=rownames(df.knots))) + 
#   geom_text() +
#   xlim(c(0, 1)) +
#   ylim(c(0, 1))


# Fitting local GPs using training data ------------------------------------------------------------
gp.list = list()

for (knot.id in 1:N.KNOTS) {
  # Get knot coordinates
  r = X.knots[knot.id, ] %>% matrix(ncol=ncol(X.knots))
  
  # Find distance from all points to current knot
  D = distance(r, X.train)
  
  # Get points within circle of radius RADIUS
  r.idx = which(D < RADIUS)
  r.nbhd = X.train[r.idx, ]
  r.nbhd.y = y.train[r.idx]
  
  # Fit spatial processes
  gp.fit = spatialProcess(x=r.nbhd, y=r.nbhd.y)
  
  # Store fitted spatial process
  gp.list[[knot.id]] = gp.fit
}


# Fit background GP -------------------------------------------------------
gp.background = spatialProcess(x=X.train, y=y.train)


# Combining GPs -----------------------------------------------------------
# Epanechnikov kernel (distance function)
epan.kernel = function(u, h) {
  (2 / pi) * (1 / h**2) * (1-(u/h)**2)
}

# Generate simulations from all spatial processes (local + background)
simulate.gp = function(gp.list,
                      gp.background,
                      X,
                      y) {
  
  # Number of samples (used for estimating variability/calculating entropy)
  num_samples = 100
  
  # Simulate from background spatial process
  sim.background = sim.spatialProcess(gp.background, X, M=num_samples)
  
  # Simulate from local spatial processes
  sim.local = list()
  for (i in 1:length(gp.list)) {
    sim.gp = sim.spatialProcess(gp.list[[i]], X, M=num_samples)
    sim.local[[i]] = sim.gp
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
        epan.kernel(u=distance(matrix(X.knots[i, ], ncol=2), X.obs), h=ell) * sim.local.i.j
    }
    # Combine local and background predictions into one row (all locations simultaneously)
    sim.total[, j] = sim.background.j + sqrt_alpha * sim.local.sum.j
  }
  
  return(sim.total)
}

# Function passed to optim to fit sqrt_alpha and ell
fit.params = function(sim.background, sim.local.list, X.knots, X.obs, y.obs, par) {
  # Extract parameters of interest
  sqrt_alpha = par[1]
  ell = par[2]
  
  sim.total = gp.preds(sim.background=sim.train.background,
                       sim.local.list=sim.train.local.list,
                       X.knots=X.knots,
                       X.obs=X.train,
                       y.obs=y.train,
                       sqrt_alpha=sqrt_alpha,
                       ell=ell)
  
  # Replicate actual values along each column to compute MSE
  y.obs.matrix = matrix(y.obs, nrow=length(y.obs), ncol=ncol(sim.background))
  
  # Compute MSE 
  mse = sum((sim.total - y.obs.matrix)**2 / length(y.obs.matrix))
  return(mse)
}


# Fit spatial process model -----------------------------------------------
# Generate predictions from all processes
sim.train.values = simulate.gp(gp.list, gp.background, X.train, y.train)
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
               X.obs=X.train,
               y.obs=y.train)


# Predict on training and testing sets ------------------------------------
# === Training set ===
best_fit_preds.train = gp.preds(sim.background=sim.background,
                          sim.local.list=sim.local.list,
                          X.knots=X.knots,
                          X.obs=X.train,
                          y.obs=y.train,
                          sqrt_alpha=params$par[1],
                          ell=params$par[2])

best_fit_preds.train.mean = apply(best_fit_preds.train, 1, mean)
best_fit_preds.train.median = apply(best_fit_preds.train, 1, median)
# Bootstrap standard error
best_fit_preds.train.se = apply(best_fit_preds.train, 1, sd) / sqrt(ncol(sim.background))

# Attach to training set
df.train$mean_pred = best_fit_preds.train.mean
df.train$median_pred = best_fit_preds.train.median
df.train$se_pred = best_fit_preds.train.se

# === Testing set ===
# Generate predictions from all processes
sim.test.values = simulate.gp(gp.list, gp.background, X.test, y.test)
# Save background and local processes separately for easy reference
sim.test.background = sim.test.values$sim.background
sim.test.local.list = sim.test.values$sim.local

best_fit_preds.test = gp.preds(sim.background=sim.test.background,
                                sim.local.list=sim.test.local.list,
                                X.knots=X.knots,
                                X.obs=X.test,
                                y.obs=y.test,
                                sqrt_alpha=params$par[1],
                                ell=params$par[2])

best_fit_preds.test.mean = apply(best_fit_preds.test, 1, mean)
best_fit_preds.test.median = apply(best_fit_preds.test, 1, median)
# Bootstrap standard error
best_fit_preds.test.se = apply(best_fit_preds.test, 1, sd) / sqrt(ncol(sim.background))

# Attach to testing set
df.test$mean_pred = best_fit_preds.test.mean
df.test$median_pred = best_fit_preds.test.median
df.test$se_pred = best_fit_preds.test.se


# Visualize predictions ---------------------------------------------------
# Mean prediction
ggplot() +
  geom_point(data=df.train, aes(x=pm2_5, y=mean_pred, color='train'), alpha=0.75) +
  geom_point(data=df.test, aes(x=pm2_5, y=mean_pred, color='test'), alpha=0.75) +
  geom_abline(aes(slope=1, intercept=0), color='black', linetype='dashed') +
  xlab('Actual PM2.5') + ylab('Predicted PM 2.5') + 
  ggtitle('Actual versus predicted (mean) PM 2.5')

# Median prediction
ggplot() +
  geom_point(data=df.train, aes(x=pm2_5, y=median_pred, color='train'), alpha=0.75) +
  geom_point(data=df.test, aes(x=pm2_5, y=median_pred, color='test'), alpha=0.75) +
  geom_abline(aes(slope=1, intercept=0), color='black', linetype='dashed') +
  xlab('Actual PM2.5') + ylab('Predicted PM 2.5') + 
  ggtitle('Actual versus predicted (median) PM 2.5')

# Standard error prediction (training set)
ggplot(df.train) +
  geom_point(aes(x=x, y=y, color=se_pred)) +
  scale_color_gradient(low='white', high='red') +
  geom_point(data=df.knots, aes(x=x, y=y), size=5, shape=10) +
  geom_circle(data=df.knots, aes(x0=x, y0=y, r=RADIUS)) +
  coord_fixed()

# Standard error prediction (testing set)
ggplot(df.test) +
  geom_point(aes(x=x, y=y, color=se_pred)) +
  scale_color_gradient(low='white', high='red') +
  geom_point(data=df.knots, aes(x=x, y=y), size=5, shape=10) +
  geom_circle(data=df.knots, aes(x0=x, y0=y, r=RADIUS)) +
  coord_fixed()







