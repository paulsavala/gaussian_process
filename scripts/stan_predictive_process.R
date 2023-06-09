# Setup -----------------------------------------------------------------------------------------------------------
set.seed(1)
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Imports (Graphing)
library(ggplot2)  # graphing; ggplot() + ...
library(ggforce)  # graphing circles; geom_circle()

# Imports (Modeling)
library(cmdstanr)
library(loo)
library(posterior)

# Imports
library(dplyr)    # piping and QoL functions; %>%
library(scales)   # normalizing lists between 0 and 1; rescale()
library(stringr)  # manipulating strings
library(mvtnorm)

# Imports (Knot Selection)
library(fields)   # cover.design; cover.design()

# Knot Evaluation
source("load_epa_data.R")
source("kaf_v6_sphere.R")
source("kaf_v6_ellipsoid.R")
source("stan_data_prep.R")

# TODO: [x] save summary.stats to file
# TODO: [x] add graphing
# TODO: [x] also save the files to the directory
# TODO: [ ] get hyperparameter testing
# TODO: [ ] maybe parallelize hyperparameter testing
# TODO: [ ]

# Data Loading ----------------------------------------------------------------------------------------------------

epa.df = load_epa_data()  # entire US

# Take a sample of the data for runtime purposes (California-ish)
# epa.df = epa.df[(epa.df$x < 0.15) & (epa.df$y < 0.75), ]
# epa.df = epa.df[epa.df$x < 0.35, ]

# show data
ggplot(epa.df) +
  geom_point(aes(x=x, y=y, fill=pm2_5), size=3, pch=21, color='black') +
  scale_fill_gradient2(low='blue', high='red', midpoint=0) +
  ggtitle('United States')

epa.sample = data.frame(epa.df)
# Get a sample of the rows as data
# epa.sample.idx = sample(1:nrow(epa.df), 250)
# epa.sample = epa.df[epa.sample.idx, ]
# rownames(epa.sample) = NULL


# Functions -------------------------------------------------------------------------------------------------------

# draws is the fit$draws() aka posterior samples
# returns a vector with a size matching the nrow of the original dataset
get_chain_preds_from_draws = function(draws, chain.number) {
  # get the draws and turn it into a data frame
  draws.df = draws %>% as_draws_df
  
  # find columns that match the pattern "R_space_knot[*,n]"
  cols = grep(paste0("R_space_knot\\[\\d+,", chain.number, "\\]"), colnames(draws.df))
  df_filtered = draws.df[, cols] %>% as.data.frame
  
  # get the mean from each column
  chain_means = colMeans(df_filtered) %>% as.vector
  
  return(chain_means)
}

# df.points contains 3 columns: `x`, `y`, `pm2_5`
# df.knots contains 3 columns: `x`, `y`, `pm2_5`
# gives you the fit object
fit_model = function(df.points, df.knots) {
  # Fitting Predictive Process Model
  model = cmdstan_model('../stan/predictive_process.stan')
  
  data = list(
    N_knots = nrow(df.knots),
    knot_locs = df.knots[, c('x', 'y')],
    N_spatial = nrow(df.points),
    spatial_locs = df.points[, c('x', 'y')],
    y_spatial = df.points$pm2_5 %>% as.vector,
    return_predictions = 1,
    return_log_likelihoods = 1
  )
  
  fit = model$sample(
    data = data,
    parallel_chains = 4,
    iter_warmup = 2000,
    max_treedepth = 10,
    init = function() list(
      sigma_z=0.1,
      ell_z=0.1,
      sigma_interp=0.1,
      ell_interp=0.1,
      lambda_y=0.1
    )
  )
  
  return(fit)
}

# also try closeAllConnections() if console output is not working
save_fit = function(fit.knots, knots, preds.knots, summary.stats, timestamp, name) {
  file.path = paste0("results/", timestamp, "/")
  if (!dir.exists(file.path)) { dir.create(file.path, recursive=T) }
  sink(paste0(file.path, name, "_results.txt"))
  print(fit.knots$loo())
  print(summary.stats)
  sink()
  fit.knots$save_object(paste0(file.path, name, "_fit.RDS"))
  knots$type = name
  write.csv(knots, paste0(file.path, name, "_knots.csv"), row.names=F)
  write.csv(preds.knots, paste0(file.path, name, "_preds.csv"), row.names=F)
}

get_preds = function(df.points, fit) {
  # Extract the predictions (these are slow to extract, so just be patient)
  fit.preds = fit$summary(variables = "y_spatial_sim")
  
  # Store the locations, true values, and predictions
  fit.preds.df = df.points[, c('x', 'y', 'pm2_5')]
  fit.preds.df$median_pred = fit.preds$median
  fit.preds.df$mean_pred = fit.preds$mean
  fit.preds.df$sd_pred = fit.preds$sd
  
  # Store the difference between the actual and predicted values
  fit.preds.df$diff_median = fit.preds.df$pm2_5 - fit.preds.df$median_pred
  fit.preds.df$diff_mean = fit.preds.df$pm2_5 - fit.preds.df$mean_pred
  
  # Also store the percentage difference
  fit.preds.df$pct_diff_median = (fit.preds.df$pm2_5 - fit.preds.df$median_pred) / fit.preds.df$pm2_5
  fit.preds.df$pct_diff_mean = (fit.preds.df$pm2_5 - fit.preds.df$mean_pred) / fit.preds.df$pm2_5
  
  return(fit.preds.df)
}

get_summary_statistics = function(preds.knots) {
  # Mean squared error
  mse = sum((preds.knots$pm2_5 - preds.knots$median_pred)**2) / nrow(preds.knots)
  # Mean absolute error
  mae = sum(abs(preds.knots$pm2_5 - preds.knots$median_pred)) / nrow(preds.knots)
  # R^2
  r2 = cor(preds.knots$pm2_5, preds.knots$median_pred)
  
  return(list(mse=mse, mae=mae, r2=r2[1, 1]))
}

# Create Different Kinds of Knots ---------------------------------------------------------------------------------

n.knots = 25

# # Random Knot Selection
# epa.sample.knots.rand.idx = sample(1:nrow(epa.sample), n.knots)
# epa.sample.knots.rand = epa.sample[epa.sample.knots.rand.idx, ]

# Uniform Grid Knot Selection
generate_grid = function(num.x, num.y) { expand.grid(x=seq(from=min(epa.sample$x), to=max(epa.sample$x), length.out=num.x), y=seq(from=min(epa.sample$y), to=max(epa.sample$y), length.out=num.y)) }
epa.sample.knots.grid = generate_grid(num.x=as.integer(sqrt(n.knots)), num.y=as.integer(sqrt(n.knots)))

# # Cover Design Knot Selection
# epa.sample.knots.cd.locs = epa.sample[, c("x", "y")]
# epa.sample.knots.cd = cover.design(epa.sample.knots.cd.locs, n.knots)$design %>% as.data.frame
# epa.sample.knots.cd = merge(epa.sample.knots.cd, epa.sample, by=c("x", "y"), sort=F)

# Isotropic Entropy Maximization Knot Selection
epa.sample.prepped = data.frame(epa.sample)
names(epa.sample.prepped) = c("x", "y", "signal")
epa.sample.prepped$z = 0
epa.sample.knots.sphere_entropy = entropy_max.sphere(epa.sample.prepped, 15, 1, n.knots)
epa.sample.knots.sphere_entropy = epa.sample.knots.sphere_entropy %>% select(-z, -radius, -entropy) %>% rename(pm2_5=signal)

ggplot(epa.sample.prepped) +
  geom_point(aes(x=x, y=y, fill=signal), size=3, pch=21, color='black') +
  geom_point(data=epa.sample.knots.sphere_entropy, aes(x=x, y=y), size=5, pch=18) +
  scale_fill_gradient2(low='blue', high='red') +
  coord_fixed()

# Anisotropic Entropy Maximization Knot Selection
get_ellipsoid_knots = function(df.points, max_knots.sphere, radius_mult.sphere, num_neighbors.sphere, spectral_components, max_knots.ellipsoid, radius_mult.ellipsoid) {
  # Circle algorithm to select knot candidates ------------------------------
  entropy.sphere.df = entropy_max.sphere(df.points, num_neighbors.sphere, radius_mult.sphere, max_knots.sphere)
  
  # Prepare data for Stan model ---------------------------------------------
  # Rename columns as needed
  knots.sphere.df = entropy.sphere.df[, c('x', 'y', 'signal')]
  knots.sphere.df$pm2_5 = knots.sphere.df$signal
  df.points$pm2_5 = df.points$signal
  
  # Prepare data for stan fitting
  data = prepare_data(df.points, M=spectral_components, plot=F, df.knots=knots.sphere.df)
  
  # Stan model fitting (fit ellipses) -----------------------------------------
  model.ellipse = cmdstan_model('../stan/aniso_process_axes_latent_knots_spectral_profile.stan')
  
  fit.ellipse = model.ellipse$sample(
    data=data,
    parallel_chains=4,
    iter_warmup=2000,
    max_treedepth=10,
    init=function() list(
      sigma_z=0.1,
      ell_z=0.1,
      sigma_interp=0.1,
      ell_interp=0.1,
      lambda_y=0.1
    )
  )
  
  # Extract the locations and PM2.5 data
  spatial.df = cbind(data$spatial_locs, data$y_spatial)
  names(spatial.df) = c('x', 'y', 'pm2_5')
  knots.df = data$knot_locs
  
  # Extract ellipses --------------------------------------------------------
  # Get the length and rotation of the ellipses
  ellipse.df = extract_ellipses(
    spatial.df,
    fit=fit.ellipse,
    plot=T,
    scale_ellipses=10, # Turn on if plotting
    psi_suffix="_all",
    return_df=TRUE
  )
  
  # Rename columns for ellipse algo
  ellipse.df$signal = ellipse.df$pm2_5
  ellipse.df$x0 = ellipse.df$x
  ellipse.df$y0 = ellipse.df$y
  ellipse.df$z0 = 0
  ellipse.df$c = 0
  
  ellipse.df$x = NULL
  ellipse.df$y = NULL
  ellipse.df$z = NULL
  ellipse.df$pm2_5 = NULL
  
  ellipse.df$alpha = atan(ellipse.df$b / ellipse.df$a)
  ellipse.df$beta = 0
  ellipse.df$gamma = 0
  
  ellipse.df = ellipse.df[, c('a', 'b', 'c', 'x0', 'y0', 'z0', 'alpha', 'beta', 'gamma', 'signal')]
  
  # Ellipse algo ------------------------------------------------------------
  knots.ellipsoid.df = entropy_max(ellipse.df, radius_mult.ellipsoid, max_knots.ellipsoid)
  return(list(knots.df=knots.ellipsoid.df, ellipse.df=ellipse.df))
}

epa.sample.knots.ellipsoid_entropy.list = get_ellipsoid_knots(
  df.points=epa.sample.prepped,
  max_knots.sphere=25,
  radius_mult.sphere=1,
  num_neighbors.sphere=15,
  spectral_components=5,
  max_knots.ellipsoid=n.knots,
  radius_mult.ellipsoid=1/50
)

# Test example
epa.sample.knots.ellipsoid_entropy.ellipses = epa.sample.knots.ellipsoid_entropy.list$ellipse.df
epa.sample.knots.ellipsoid_entropy = epa.sample.knots.ellipsoid_entropy.list$knots.df

knots.ellipsoid.df = entropy_max(epa.sample.knots.ellipsoid_entropy.ellipses, 1/5, 25)

ggplot(epa.sample.knots.ellipsoid_entropy.ellipses) +
  geom_point(aes(x=x0, y=y0, fill=signal), size=3, pch=21, color='black') +
  geom_ellipse(data=knots.ellipsoid.df, aes(x0=x0, y0=y0, a=a/7, b=b/7, angle=alpha)) +
  scale_fill_gradient2(low='blue', high='red') +
  coord_fixed()

knots.ellipsoid.parsed.df = knots.ellipsoid.df[, c('x0', 'y0', 'signal')] %>% rename(x=x0, y=y0, pm2_5=signal)
fit.knots.ellipsoid_entropy = fit_model(epa.sample, knots.ellipsoid.parsed.df)
preds.knots.ellipsoid_entropy = get_preds(epa.sample, fit.knots.ellipsoid_entropy)
preds.knots.ellipsoid_entropy$type = paste0("ellipsoid_entropy_nnsphere_", n.neighbors, "_rmellipsoid_", radius.mult)
summary.stats.ellipsoid_entropy = get_summary_statistics(preds.knots.ellipsoid_entropy)

ggplot(preds.knots.ellipsoid_entropy) +
  geom_point(aes(x=x, y=y, fill=log(abs((pm2_5-median_pred)/pm2_5))), color='black', size=3, pch=21) +
  scale_fill_gradient(low='white', high='red') +
  coord_fixed() +
  ggtitle('Anisotropic knot predictions') +
  # guides(fill=guide_legend(title="Log error as % of PM2.5")) +
  theme(text = element_text(size=20))
# End test example

ggplot(epa.sample.knots.ellipsoid_entropy.ellipses) +
  geom_point(aes(x=x0, y=y0, fill=signal), size=3, pch=21, color='black') +
  geom_ellipse(data=epa.sample.knots.ellipsoid_entropy, aes(x0=x0, y0=y0, a=a/10, b=b/10, angle=alpha)) +
  scale_fill_gradient2(low='blue', high='red') +
  coord_fixed() +
  theme(text = element_text(size=20))

epa.sample.knots.ellipsoid_entropy = epa.sample.knots.ellipsoid_entropy %>% select(x0, y0) %>% rename(x = x0, y = y0)

# Trying It Out ---------------------------------------------------------------------------------------------------

run.timestamp = Sys.time() %>% format("%Y-%m-%d_%H_%M_%S")

# # Random Knot Selection
# fit.knots.rand = fit_model(epa.sample, epa.sample.knots.rand)
# preds.knots.rand = get_preds(epa.sample, fit.knots.rand)
# preds.knots.rand$type = "random"
# summary.stats.rand = get_summary_statistics(preds.knots.rand)
# save_fit(fit.knots=fit.knots.rand, preds.knots=preds.knots.rand, summary.stats=summary.stats.rand, timestamp=run.timestamp, name="random")

# Uniform Grid Knot Selection
fit.knots.grid = fit_model(epa.sample, epa.sample.knots.grid)
preds.knots.grid = get_preds(epa.sample, fit.knots.grid)
preds.knots.grid$type = "grid"
summary.stats.grid = get_summary_statistics(preds.knots.grid)
save_fit(fit.knots=fit.knots.grid, knots=epa.sample.knots.grid, preds.knots=preds.knots.grid, summary.stats=summary.stats.grid, timestamp=run.timestamp, name="grid")

# Cover Design Knot Selection
fit.knots.cd = fit_model(epa.sample, epa.sample.knots.cd)
preds.knots.cd = get_preds(epa.sample, fit.knots.cd)
preds.knots.cd$type = "cover_design"
summary.stats.cd = get_summary_statistics(preds.knots.cd)
save_fit(fit.knots=fit.knots.cd, knots=epa.sample.knots.cd, preds.knots=preds.knots.cd, summary.stats=summary.stats.cd, timestamp=run.timestamp, name="cover_design")

epa.sample.prepped = data.frame(epa.sample)
# Isotropic Entropy Maximization Knot Selection
# Hyperparameter evaluation
for (n.neighbors in c(5, 7, 10, 15, 20)) {
  for (radius.mult in c(0.75, 1, 1.25)) {
    # Get knots using hyperparameters
    epa.sample.prepped = data.frame(epa.sample)
    names(epa.sample.prepped) = c("x", "y", "signal")
    epa.sample.prepped$z = 0
    epa.sample.knots.sphere_entropy = entropy_max.sphere(epa.sample.prepped, n.neighbors, radius.mult, n.knots)
    epa.sample.knots.sphere_entropy = epa.sample.knots.sphere_entropy %>% select(-z, -radius, -entropy) %>% rename(pm2_5=signal)
    
    # Fit model
    fit.knots.sphere_entropy = fit_model(epa.sample, epa.sample.knots.sphere_entropy)
    preds.knots.sphere_entropy = get_preds(epa.sample, fit.knots.sphere_entropy)
    preds.knots.sphere_entropy$type = paste0("sphere_entropy_nn_", n.neighbors, "_rm_", radius.mult)
    summary.stats.sphere_entropy = get_summary_statistics(preds.knots.sphere_entropy)
    save_fit(fit.knots=fit.knots.sphere_entropy, knots=epa.sample.knots.sphere_entropy, preds.knots=preds.knots.sphere_entropy, summary.stats=summary.stats.sphere_entropy, timestamp=run.timestamp, name=paste0("sphere_entropy_nn_", n.neighbors, "_rm_", radius.mult))
  }
}


for (n.neighbors in c(5, 7, 10, 15, 20)) {
  for (radius.mult in c(1/50, 1/25, 1/15, 1/10)) {
    epa.sample.knots.ellipsoid_entropy = get_ellipsoid_knots(
      df.points=epa.sample.prepped,
      max_knots.sphere=n.knots,
      radius_mult.sphere=1,
      num_neighbors.sphere=n.neighbors,
      spectral_components=5,
      max_knots.ellipsoid=n.knots,
      radius_mult.ellipsoid=radius.mult
    )
    
    epa.sample.knots.ellipsoid_entropy = epa.sample.knots.ellipsoid_entropy %>% select(x0, y0) %>% rename(x = x0, y = y0)
    
    # Anisotropic Entropy Maximization Knot Selection
    fit.knots.ellipsoid_entropy = fit_model(epa.sample, epa.sample.knots.ellipsoid_entropy)
    preds.knots.ellipsoid_entropy = get_preds(epa.sample, fit.knots.ellipsoid_entropy)
    preds.knots.ellipsoid_entropy$type = paste0("ellipsoid_entropy_nnsphere_", n.neighbors, "_rmellipsoid_", radius.mult)
    summary.stats.ellipsoid_entropy = get_summary_statistics(preds.knots.ellipsoid_entropy)
    save_fit(fit.knots=fit.knots.ellipsoid_entropy, knots=epa.sample.knots.ellipsoid_entropy, preds.knots=preds.knots.ellipsoid_entropy, summary.stats=summary.stats.ellipsoid_entropy, timestamp=run.timestamp, name=paste0("ellipsoid_entropy_nnsphere_", n.neighbors, "_rmellipsoid_", radius.mult))
  }
}


# Viewing Results -------------------------------------------------------------------------------------------------

# summary.stats.rand
summary.stats.grid
summary.stats.cd
summary.stats.sphere_entropy
# summary.stats.ellipsoid_entropy

# compare predictions between all sets of knots
# preds.knots.all.df = rbind(preds.knots.rand, preds.knots.grid, preds.knots.cd, preds.knots.sphere_entropy, preds.knots.ellipsoid_entropy)
preds.knots.all.df = rbind(preds.knots.grid, preds.knots.cd, preds.knots.sphere_entropy)

ggplot(preds.knots.all.df, aes(x=pm2_5, y=median_pred, color=type)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  geom_abline(slope=1, intercept=0, color='blue')


# Ellipses vs grid graph --------------------------------------------------
preds.knots.ellipsoid_entropy$type = 'ellipse'
preds.knots.ellipse_grid = rbind(preds.knots.ellipsoid_entropy, preds.knots.grid)

ggplot(preds.knots.ellipse_grid %>% filter(x<0.2, y<0.7)) +
  geom_point(aes(x=x, y=y, fill=log(abs((pm2_5-median_pred)/pm2_5))), color='black', size=3, pch=21) +
  scale_fill_gradient(low='blue', high='white') +
  facet_wrap(vars(type)) +
  coord_fixed() +
  ggtitle('Anisotropic knot predictions') +
  # guides(fill=guide_legend(title="Log error as % of PM2.5")) +
  theme(text = element_text(size=20))

