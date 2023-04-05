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

source("scripts/fuentes_GP_model.R")

# # Helper functions --------------------------------------------------------
# # Epanechnikov kernel (distance function)
# epan.kernel = function(u, h) {
#   (2 / pi) * (1 / h**2) * (1-(u/h)**2)
# }
# 
# generate.gp = function(X.knots, X, y, only.background=FALSE) {
#   N.KNOTS = nrow(X.knots)
#   RADIUS = 0.2
#   
#   # Empty list to hold all local spatial processes
#   gp.list = list()
#   
#   if (!only.background) {
#     for (knot.id in 1:N.KNOTS) {
#       # Get knot coordinates
#       r = X.knots[knot.id, ] %>% matrix(ncol=ncol(X.knots))
#       
#       # Find distance from all points to current knot
#       D = rdist(r, X)
#       
#       # Get points within circle of radius RADIUS
#       r.idx = which(D < RADIUS)
#       r.nbhd = X[r.idx, ]
#       r.nbhd.y = y[r.idx]
#       
#       # Fit spatial processes
#       gp.fit = spatialProcess(x=r.nbhd, y=r.nbhd.y)
#       
#       # Store fitted spatial process
#       gp.list[[knot.id]] = gp.fit
#     }
#   }
#   
#   # Fit background GP -------------------------------------------------------
#   gp.background = spatialProcess(x=X, y=y)
#   
#   return(list("gp.background"=gp.background,
#               "gp.list"=gp.list))
# }
# 
# # Generate simulations from all spatial processes (local + background)
# simulate.gp = function(gp.list,
#                        gp.background,
#                        X) {
#   
#   # Number of samples (used for estimating variability/calculating entropy)
#   num_samples = 100
#   
#   # Simulate from background spatial process
#   sim.background = sim.spatialProcess(gp.background, X, M=num_samples)
#   
#   # Simulate from local spatial processes
#   sim.local = list()
#   if (!is.null(gp.list)) {
#     for (i in 1:length(gp.list)) {
#       sim.gp = sim.spatialProcess(gp.list[[i]], X, M=num_samples)
#       sim.local[[i]] = sim.gp
#     }
#   }
#   
#   return(list("sim.background"=sim.background,
#               "sim.local"=sim.local))
# }
# 
# gp.preds = function(sim.background, sim.local.list, X.knots, X.obs, y.obs, sqrt_alpha, ell) {
#   # Dimensions:
#   # sim.background: (num locations) x (num samples generated)
#   # sim.local: same
#   # epan.kernel(X.knots, X.obs): (num knots) x (num locations)
#   
#   # Holds predictions for each sample
#   sim.total = matrix(NA, nrow=nrow(sim.background), ncol=ncol(sim.background))
#   
#   # Go sample-by-sample
#   for (j in 1:ncol(sim.background)) {
#     # Just the j'th sample from background spatial process
#     sim.background.j = sim.background[, j]
#     # Empty vector to hold the sum of all local spatial processes
#     sim.local.sum.j = 0 * sim.background.j
#     # Loop over spatial processes
#     for (i in 1:length(sim.local.list)) {
#       # Just the j'th sample of the i'th local spatial process
#       sim.local.i.j = sim.local.list[[i]][, j]
#       sim.local.sum.j = sim.local.sum.j + 
#         epan.kernel(u=rdist(matrix(X.knots[i, ], ncol=2), X.obs), h=ell) * sim.local.i.j
#     }
#     # Combine local and background predictions into one row (all locations simultaneously)
#     sim.total[, j] = sim.background.j + sqrt_alpha * sim.local.sum.j
#   }
#   
#   return(sim.total)
# }
# 
# # Function passed to optim to fit sqrt_alpha and ell
# fit.params = function(sim.background, sim.local.list, X.knots, X.train, y.train, par) {
#   # Extract parameters of interest
#   sqrt_alpha = par[1]
#   ell = par[2]
#   
#   sim.total = gp.preds(sim.background=sim.background,
#                        sim.local.list=sim.local.list,
#                        X.knots=X.knots,
#                        X.obs=X.train,
#                        y.obs=y.train,
#                        sqrt_alpha=sqrt_alpha,
#                        ell=ell)
#   
#   # Replicate actual values along each column to compute MSE
#   y.train.matrix = matrix(y.train, nrow=length(y.train), ncol=ncol(sim.background))
#   
#   # Compute MSE 
#   mse = sum((sim.total - y.train.matrix)**2 / length(y.train.matrix))
#   return(mse)
# }

# Fitting local GPs using training data ------------------------------------------------------------
# fit.fuentes.gp = function(X.knots, X.train, y.train) {
#   N.KNOTS = nrow(X.knots)
#   
#   gp.all = generate.gp(X.knots, X.train, y.train)
#   gp.background = gp.all$gp.background
#   gp.list = gp.all$gp.list
#   
#   # Generate predictions from all processes
#   sim.train.values = simulate.gp(gp.list, gp.background, X.train)
#   # Save background and local processes separately for easy reference
#   sim.train.background = sim.train.values$sim.background
#   sim.train.local.list = sim.train.values$sim.local
#   
#   params = optim(par=c(1, 1), fn=fit.params, 
#                  method="L-BFGS",
#                  lower=c(0, 0),
#                  upper=c(10, 1),
#                  sim.background=sim.train.background,
#                  sim.local.list=sim.train.local.list,
#                  X.knots=X.knots,
#                  X.train=X.train,
#                  y.train=y.train)
#   
#   sqrt_alpha = params$par[1]
#   ell = params$par[2]
#   
#   return(list("sqrt_alpha"=sqrt_alpha, 
#               "ell"=ell))
# }
# 
# predict.fuentes.gp = function(X.knots, X, y, params) {
#   gp.all = generate.gp(X.knots, X, y)
#   gp.background = gp.all$gp.background
#   gp.list = gp.all$gp.list
#   
#   # Generate predictions from all processes
#   sim.values = simulate.gp(gp.list, gp.background, X)
#   # Save background and local processes separately for easy reference
#   sim.background = sim.values$sim.background
#   sim.local.list = sim.values$sim.local
#   
#   best_fit_preds = gp.preds(sim.background=sim.background,
#                             sim.local.list=sim.local.list,
#                             X.knots=X.knots,
#                             X.obs=X,
#                             y.obs=y,
#                             sqrt_alpha=params$sqrt_alpha,
#                             ell=params$ell)
#   
#   best_fit_preds.mean = apply(best_fit_preds, 1, mean)
#   best_fit_preds.median = apply(best_fit_preds, 1, median)
#   # Bootstrap standard error
#   best_fit_preds.se = apply(best_fit_preds, 1, sd) / sqrt(ncol(sim.background))
#   
#   best_fit_preds.df = data.frame(mean_pred=best_fit_preds.mean,
#                                  median_pred=best_fit_preds.median,
#                                  se_pred=best_fit_preds.se)
#   return(best_fit_preds.df)
# }


# Fit Fuentes GP on a number of different knots ---------------------------
fit.fuentes.knots = function(n.knots, X.train, y.train, n.knot_sets=25, radius=0.2, save_all_steps=FALSE) {
  knot_set.list.X = list()
  knot_set.list.y = list()
  knot_set.info.df = data.frame(id=NA, se=NA, utility=NA)
  
  for (knot_set.idx in 1:n.knot_sets) {
    print(paste0("Fitting knot set ", knot_set.idx, "/", n.knot_sets, "..."))
    # 
    # # Create a rough lattice
    # x.lattice = seq(min(X.train$x), max(X.train$x), length.out=25)
    # y.lattice = seq(min(X.train$y), max(X.train$y), length.out=25)
    # 
    # grid.df = make.surface.grid(list(x=x.lattice, y=y.lattice)) %>% as.data.frame
    # 
    # # Remove lattice points that have no training points to the left/right or above/below
    # rows_to_delete = c()
    # for (i in 1:nrow(grid.df)) {
    #   # Get lattice point coordinates
    #   r = grid.df[i, ]
    #   
    #   # Find distance to all training points
    #   D = rdist(r, X.train)
    #   
    #   # Find how many points it's within half the radius from
    #   r.idx = which(D < radius / 2)
    #   
    #   # Remove points that don't have at least 10 neighbors
    #   if (length(r.idx) < 10) {
    #     rows_to_delete = c(rows_to_delete, i)
    #   }
    # }
    # 
    # # Remove grid points too far from the training data
    # grid.df = grid.df[-rows_to_delete, ]
    
    # Sample n.knots random knots
    X.knots.idx = sample(1:nrow(X.train), size=n.knots, replace=FALSE)
    X.knots = X.train[X.knots.idx, ]
    
    # Fit Fuentes model hyperparameters
    params = fit.fuentes.gp(X.knots, X.train, y.train)
    preds.df = predict.fuentes.gp(X.knots, X.train, params)
    
    # === Calculate utility ===
    # Calculate median PM 2.5 in neighborhood of knot
    knot_set.median_pm2_5 = c()
    for (i in 1:n.knots) {
      # Get knot coordinates
      r = X.knots[i, ] %>% matrix(ncol=ncol(X.knots))
      
      # Find distance from all points to current knot
      D = rdist(r, X.train)
      
      # Get points within circle of radius
      r.idx = which(D < radius)
      r.nbhd.y = y.train[r.idx]
      
      # Calculate median PM 2.5 in this radius
      r.pm2_5 = median(r.nbhd.y)
      
      # Record this median
      knot_set.median_pm2_5 = c(knot_set.median_pm2_5, r.pm2_5)
    }
    
    # Calculate background GP using training data
    gp.background = generate.gp(X.knots, X.train, y.train)$gp.background
    
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
    
    if (save_all_steps) {
      saveRDS(list("knot_set.info.df"=knot_set.info.df,
                   "knot_set.list.X"=knot_set.list.X,
                   "knot_set.list.y"=knot_set.list.y),
              paste0("step_", knot_set.idx, ".rds"))
    }
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
  
  return(X.knots)
}












