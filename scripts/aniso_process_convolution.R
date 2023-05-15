# Setup -----------------------------------------------------------------------------------------------------------
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

# Knot Evaluation
source("scripts/load_epa_data.R")


# Data loading ------------------------------------------------------------
epa.df = load_epa_data()

ggplot(epa.df) +
  geom_point(aes(x=x, y=y, fill=pm2_5), size=3, pch=21, color='black') +
  scale_fill_gradient2(low='blue', high='red', midpoint=0)


# Stan model fitting ------------------------------------------------------
# Prepare data for the stan model
prepare_data = function(df, M, num_knots=NULL, sample_pct=1, plot=FALSE, return_ellipses=TRUE, return_predictions=FALSE, return_log_likelihoods=FALSE, df.sample=NULL, df.knots=NULL) {
  # Reduce sample size if needed
  if (is.null(df.sample)) {
    if (sample_pct < 1) {
      df.sample.idx = sample(1:nrow(df), as.integer(sample_pct*nrow(df)))
      df.sample = df[df.sample.idx, ]
      rownames(df.sample) = 1:nrow(df.sample)
    } else {
      df.sample = df
    }
  }
  
  # Spatial process sites and data
  spatial_locs = df.sample[, c('x', 'y')]
  y_spatial = df.sample$pm2_5 %>% as.vector
  N_spatial = nrow(spatial_locs)
  
  # Package data
  data = list(N_spatial=N_spatial,
              spatial_locs=spatial_locs,
              y_spatial=y_spatial,
              M=M,
              return_ellipses=return_ellipses,
              return_predictions=return_predictions,
              return_log_likelihoods=return_log_likelihoods)
  
  # Include knots if needed
  if (!is.null(num_knots)) {
    # Randomly select num_knots
    df.knots.idx = sample(1:nrow(df.sample), size=num_knots)
    df.knots = df.sample[df.knots.idx, ]
  }
  
  data$N_knots = nrow(df.knots)
  data$knot_locs = df.knots[, c('x', 'y')]
  
  if (plot) {
    g = ggplot(df) +
      geom_point(aes(x=x, y=y, fill=pm2_5), size=3, pch=21, color='black', alpha=0.1) +
      geom_point(data=df.sample, aes(x=x, y=y, fill=pm2_5), size=3, pch=21, color='black') +
      scale_fill_gradient2(low='blue', high='red', midpoint=0) +
      ggtitle('US (sample points highlighted)')
    
    if (!is.null(num_knots)) {
      g = g + geom_point(data=df.knots, aes(x=x, y=y), shape=23, size=5, fill='black') +
        coord_fixed()
    }
    show(g)
  }
  
  return(data)
}


# Extract ellipses for plotting ------------------------------------------
# After the model is fit, extract (and optionally plot) the ellipses
extract_ellipses = function(df, df.all=NULL, fit, plot=TRUE, scale_ellipses=1, num_ellipses=NULL, knots=NULL, psi_suffix='', return_df=FALSE) {
  psi_x_col = paste0('psi_x', psi_suffix)
  psi_y_col = paste0('psi_y', psi_suffix)
  
  a = fit$summary(variables = psi_x_col)$median
  b = fit$summary(variables = psi_y_col)$median
  rotation = atan(b / a)
  
  # Append rotations and ellipse major/minor axes to original data
  df$a = a
  df$b = b
  df$alpha = rotation
  
  # If knots are supplied, extract their axes/rotation for plotting
  if (!is.null(knots)) {
    knots.df = merge(knots, df, by=c('x', 'y'), all.x=TRUE)
  }
  
  # Plot ellipses on map
  if (plot) {
    g = ggplot(df) +
          geom_point(aes(x=x, y=y, fill=pm2_5), size=3, pch=21, color='black') +
          scale_fill_gradient2(low='blue', high='red', midpoint=0) +
          ggtitle('US (sample points highlighted)') +
          coord_fixed()
    
    if (!is.null(df.all)) {
      g = g + geom_point(aes(x=x, y=y, fill=pm2_5), size=3, pch=21, color='black', alpha=0.1)
    }
    
    if (!is.null(num_ellipses)) {
      ellipse.idx = sample(1:nrow(df), size=num_ellipses)
      g = g + geom_ellipse(data=df[ellipse.idx,], aes(x0=x, y0=y, a=a/scale_ellipses, b=b/scale_ellipses, angle=alpha)) +
              geom_point(data=df[ellipse.idx,], aes(x=x, y=y), shape=13)
    } else if (!is.null(knots)) {
      g = g + geom_ellipse(data=knots.df, aes(x0=x, y0=y, a=a/scale_ellipses, b=b/scale_ellipses, angle=alpha)) +
        geom_point(data=knots.df, aes(x=x, y=y), shape=13)
    } else {
      g = g + geom_ellipse(data=df, aes(x0=x, y0=y, a=a/scale_ellipses, b=b/scale_ellipses, angle=alpha))
    }
    show(g)
  }
  
  if (return_df) {
    return(df)
  }
}


# # Fitting knot model ------------------------------------------------------
model.knots = cmdstan_model('stan/aniso_process_axes_latent_knots_spectral_profile.stan')

# Data
data.knots = prepare_data(epa.df, M=5, num_knots=15, sample_pct=0.32, plot=T, df.sample=data.knots.spatial.df)

fit.knots = model.knots$sample(data=data.knots,
                               parallel_chains=4,
                               iter_warmup=2000,
                               max_treedepth=10,
                               init=function() list(sigma_z=0.1,
                                                    ell_z=0.1,
                                                    sigma_interp=0.1,
                                                    ell_interp=0.1,
                                                    lambda_y=0.1))
#   
#   sink(paste0("M_", data.knots$M, "_numknots_", data.knots$N_knots, ".txt"))
#   paste0("N_spatial = ", data.knots$N_spatial, "\n") %>% cat
#   paste0("N_knots = ", data.knots$N_knots, "\n\n") %>% cat
#   paste0("M = ", data.knots$M) %>% cat
#   fit.knots$loo() %>% print
#   sink()
# }

# Extract the locations and data
data.knots.spatial.df = cbind(data.knots$spatial_locs, data.knots$y_spatial)
names(data.knots.spatial.df) = c('x', 'y', 'pm2_5')
data.knots.knots.df = data.knots$knot_locs

# # Plot the ellipses
# extract_ellipses(data.knots.spatial.df, 
#                  fit=fit.knots, 
#                  scale_ellipses=7, 
#                  # knots=data.knots.knots.df,
#                  psi_suffix="_all")
# 
# Extract the data frame with ellipse data
ellipse_df = extract_ellipses(data.knots.spatial.df,
                 fit=fit.knots,
                 scale_ellipses=10,
                 psi_suffix="_all",
                 return_df=TRUE)


# Fitting predictive process model ----------------------------------------
# Fitting knot model ------------------------------------------------------
# Generate initial data to re-use for all fits
data.raw = prepare_data(epa.df, M=0, num_knots=20, plot=T, sample_pct=0.32)

model.pp = cmdstan_model('stan/predictive_process.stan')

# Data
data.knots.spatial.df = cbind(data.raw$spatial_locs, data.raw$y_spatial)
names(data.knots.spatial.df) = c('x', 'y', 'pm2_5')

for (sample_pct in c(0.3, 0.4, 0.5, 0.6, 0.7)) {
  # data.pp = prepare_data(epa.df, M=0, num_knots=k, plot=T, df.sample=data.knots.spatial.df)
  data.pp = prepare_data(epa.df, M=0, num_knots=as.integer(sqrt(sample_pct*nrow(epa.df))), sample_pct=sample_pct)
  
  # ======== Fit PP model ===========
  data.pp$M = NULL
  data.pp$return_ellipses = NULL
  data.pp$return_predictions = 0
  
  fit.pp.start = Sys.time()
  print(paste("PP", "# knots =", k))
  fit.pp = model.pp$sample(data=data.pp,
                              parallel_chains=4,
                              iter_warmup=2000,
                              max_treedepth=10,
                              refresh=3000,
                              init=function() list(lambda_y=0.1,
                                                   ell_psi=0.1,
                                                   sigma_psi=0.1,
                                                   nugget_psi=0.1))
  
  # lp_summary = fit.pp$summary(variables="lp__")
  sink(paste0("PP_", "_numknots_", data.pp$N_knots, ".txt"))
  fit.pp.end = Sys.time()
  paste(fit.pp.end - fit.pp.start, "seconds elapsed \n") %>% cat
  paste0("N_spatial = ", data.pp$N_spatial, "\n") %>% cat
  paste0("N_knots = ", data.pp$N_knots, "\n\n") %>% cat
  # paste0("r_hat = ", lp_summary$rhat, "\n\n") %>% cat
  # paste0("ess_bulk = ", lp_summary$ess_bulk, "\n\n") %>% cat
  # paste0("ess_tail = ", lp_summary$ess_tail, "\n\n") %>% cat
  fit.pp$loo() %>% print
  sink()
  
  # ======== Fit spectral model ===========
  for (m in c(3, 5, 7, 10)) {
    data.pp$M = m
    data.pp$return_predictions = 0
    data.pp$return_ellipses = 0
    print(paste("SPEC", "# knots =", k, "M =", m))
    
    fit.start = Sys.time()
    fit.knots = model.knots$sample(data=data.pp,
                                   parallel_chains=4,
                                   iter_warmup=2000,
                                   max_treedepth=10,
                                   refresh=3000,
                                   init=function() list(sigma_z=0.1,
                                                        ell_z=0.1,
                                                        sigma_interp=0.1,
                                                        ell_interp=0.1,
                                                        lambda_y=0.1))
    fit.end = Sys.time()
    
    lp_summary = fit.knots$summary(variables="lp__")
    sink(paste0("SPEC_M_", data.pp$M, "_numknots_", data.pp$N_knots, ".txt"))
    paste(fit.end - fit.start, "seconds elapsed \n") %>% cat
    paste0("N_spatial = ", data.pp$N_spatial, "\n") %>% cat
    paste0("N_knots = ", data.pp$N_knots, "\n\n") %>% cat
    paste0("r_hat = ", lp_summary$rhat, "\n\n") %>% cat
    paste0("ess_bulk = ", lp_summary$ess_bulk, "\n\n") %>% cat
    paste0("ess_tail = ", lp_summary$ess_tail, "\n\n") %>% cat
    fit.knots$loo() %>% print
    sink()
  }
}


# Runtime for increasing knots --------------------------------------------
spec.results.df = data.frame(
  model='SPEC',
  m=c(rep(3, 6), 
      rep(5, 6), 
      rep(7, 6), 
      rep(10, 6)),
  k=c(rep(c(5, 10, 15, 20, 25, 30), 4)),
  t=c(24.0, 38.5, 47.6, 71.4, 108.6, 102,
      52.3, 60.0, 94.8, 163.8, 103.2, 183.0,
      81.0, 94.2, 119.4, 100.8, 288, 244.8,
      123.0, 113.4, 117.6, 192.0, 225.0, 288.0),
  elpd_loo=c(-344.6, -336.0, -338.0, -336.7, -331.9, -343.2,
             -345.9, -332.7, -338.9, -338.5, -327.3, -341.1,
             -341.9, -332.5, -337.1, -332.4, -334.4, -341.7,
             -340.4, -328.8, -337.8, -338.4, -331.2, -341.8),
  elpd_loo_se=c(14.7, 15.1, 15.1, 15.8, 15.6, 15.0,
                14.7, 15.0, 15.2, 15.6, 16.1, 14.9,
                14.8, 14.8, 14.6, 16.1, 15.7, 14.9,
                14.9, 15.7, 15.5, 16.6, 16.0, 14.8)
)

gpp.results.df = data.frame(
  model='GPP',
  k=c(5, 10, 15, 20, 25, 30),
  t=c(45.6, 92.4, 135.6, 196.2, 363.6, 435.0),
  elpd_loo=c(-337.4, -342.6, -320.9, -322.2, -322.9, -322.3),
  elpd_loo_se=c(13.5, 13.9, 16.6, 17.0, 16.3, 16.3)
)

ggplot() +
  geom_line(data=spec.results.df, aes(x=k, y=elpd_loo, color=as.factor(m))) +
  geom_line(data=gpp.results.df, aes(x=k, y=elpd_loo)) +
  theme(text = element_text(size=20)) +
  xlab('Number of knots') +
  ylab('ELPD-LOO') +
  guides(color=guide_legend(title="Spectral components"))

ggplot() +
  geom_line(data=spec.results.df, aes(x=k, y=t, color=as.factor(m))) +
  geom_line(data=gpp.results.df, aes(x=k, y=t), linetype='dashed') +
  theme(text = element_text(size=20)) +
  xlab('Number of knots') +
  ylab('Runtime (s)') +
  guides(color=guide_legend(title="Spectral components"))


# Runtime for increasing spacial points -----------------------------------
spec.results.df = data.frame(
  model='SPEC',
  m=c(rep(3, 5), 
      rep(5, 5), 
      rep(7, 5), 
      rep(10, 5)),
  n=c(rep(c(233, 311, 389, 466, 544), 4)),
  t=c(48.1, 684, 918, 1302, 1158, 
      66, 294, NA, 2058, 3780,
      102, 1314, 2850, 4968, 8388, 
      144, 840, 462, 2958, 9000),
  elpd_loo=c(-280.2, -400.6, -496.3, -590.0, -681.0, 
             -275.3, -399.5, -502.9, -582.4, -664.3,
             -273.1, -407.2, -500.9, -566.7, -667.2, 
             -279.3, -403.9, -461.7, -576.4, -693.1),
  elpd_loo_se=c(10.9, 13.0, 18.9, 21.6, 21.4, 
                11.3, 12.6, 20.6, 22.2, 22.7,
                11.5, 13.3, 17.9, 23.5, 23.8, 
                11.6, 13.6, 21.7, 21.1, 23.7)
)

gpp.results.df = data.frame(
  model='GPP',
  n=c(233, 311, 389, 466, 544),
  t=c(96, 234, 2094, 5400, 4320),
  elpd_loo=c(-263.5, -376.5, -449.2, -542.4, -618.7),
  elpd_loo_se=c(10.8, 15.9, 19.2, 23.2, 23.2)
)

ggplot() +
  geom_line(data=spec.results.df, aes(x=n, y=elpd_loo, color=as.factor(m))) +
  geom_line(data=gpp.results.df, aes(x=n, y=elpd_loo)) +
  theme(text = element_text(size=20)) +
  xlab('Number of knots') +
  ylab('ELPD-LOO') +
  guides(color=guide_legend(title="Spectral components"))

ggplot() +
  geom_line(data=spec.results.df, aes(x=n, y=t, color=as.factor(m))) +
  geom_line(data=gpp.results.df, aes(x=n, y=t), linetype='dashed') +
  theme(text = element_text(size=20)) +
  xlab('Number of knots') +
  ylab('Runtime (s)') +
  guides(color=guide_legend(title="Spectral components"))
