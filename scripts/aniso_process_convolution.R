# Setup -----------------------------------------------------------------------------------------------------------
# Imports (Graphing)
library(ggplot2)  # graphing; ggplot() + ...
library(ggforce)  # graphing circles; geom_circle()

# Imports (Modeling)
library(cmdstanr)

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

# # Sample of US ------------------------------------------------------------
# epa.sample.idx = sample(1:nrow(epa.df), 100)
# epa.sample.df = epa.df[epa.sample.idx, ]
# 
# ggplot(epa.df) +
#   geom_point(aes(x=x, y=y, fill=pm2_5), size=3, pch=21, color='black', alpha=0.1) +
#   geom_point(data=epa.sample.df, aes(x=x, y=y, fill=pm2_5), size=3, pch=21, color='black') +
#   scale_fill_gradient2(low='blue', high='red', midpoint=0) +
#   ggtitle('US (sample points highlighted)')
#
# # California (roughly) ----------------------------------------------------
# # Filter to approximately California
# epa.ca.df = epa.df %>%
#   dplyr::filter(x < 0.2, y > 0.25, y < 0.75) %>%
#   dplyr::filter((x < 0.125) | (y < 0.4))
# 
# epa.ca.df = epa.ca.df[, c('x', 'y', 'pm2_5')]
# 
# # Sample for easy testing
# epa.ca.sample.idx = sample(1:nrow(epa.ca.df), 100)
# epa.ca.sample.df = epa.ca.df[epa.ca.sample.idx, ]
# 
# # Plot CA and sample
# ggplot(epa.ca.df) +
#   geom_point(aes(x=x, y=y, fill=pm2_5), size=3, pch=21, color='black', alpha=0.1) +
#   geom_point(data=epa.ca.sample.df, aes(x=x, y=y, fill=pm2_5), size=3, pch=21, color='black') +
#   scale_fill_gradient2(low='blue', high='red', midpoint=0) +
#   ggtitle('California (sample points highlighted)')


# Stan model fitting ------------------------------------------------------
prepare_data = function(df, num_knots=NULL, sample_pct=1, plot=FALSE) {
  # Reduce sample size if needed
  if (sample_pct < 1) {
    df.sample.idx = sample(1:nrow(df), as.integer(sample_pct*nrow(df)))
    df.sample = df[df.sample.idx, ]
    rownames(df.sample) = 1:nrow(df.sample)
  } else {
    df.sample = df
  }
  
  # Spatial process sites and data
  spatial_locs = df.sample[, c('x', 'y')]
  y_spatial = df.sample$pm2_5 %>% as.vector
  N_spatial = nrow(spatial_locs)
  
  # Package data
  data = list(N_spatial=N_spatial,
              spatial_locs=spatial_locs,
              y_spatial=y_spatial)
  
  # Include knots if needed
  if (!is.null(num_knots)) {
    # Randomly select num_knots
    df.knots.idx = sample(1:nrow(df.sample), size=num_knots)
    df.knots = df.sample[df.knots.idx, ]
    
    data$N_knots = nrow(df.knots)
    data$knot_locs = df.knots[, c('x', 'y')]
  }
  
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
extract_ellipses = function(df, df.all=NULL, fit, foci=TRUE, plot=TRUE, scale_ellipses=1, num_ellipses=NULL, knots=NULL, psi_suffix='') {
  psi_x_col = paste0('psi_x', psi_suffix)
  psi_y_col = paste0('psi_y', psi_suffix)
  
  if (foci) {
    foci.x = fit$summary(variables = psi_x_col)$median
    foci.y = fit$summary(variables = psi_y_col)$median
    rotation = atan(foci.y / foci.x)
    
    focus = c()
    for (i in 1:length(rotation)) {
      # Rotate foci backwards to get focus
      neg_rot_mat = matrix(c(cos(-rotation[i]), -sin(-rotation[i]), sin(-rotation[i]), cos(-rotation[i])),
                           byrow=TRUE, nrow=2)
      
      focus.vec = neg_rot_mat %*% t(cbind(foci.x[i], foci.y[i]))
      # Extract just the x variable, since the vector has been rotated to lie on the x-axis
      focus = c(focus, focus.vec[1, 1])
    }
    
    # Major and minor axes (equations from Higdon's paper, I _can_ duplicate those equations)
    a = sqrt(sqrt(4*data$A**2 + focus**4*pi**2) / (2*pi) + focus**2/2)
    b = sqrt(sqrt(4*data$A**2 + focus**4*pi**2) / (2*pi) - focus**2/2)
  } else {
    a = fit$summary(variables = psi_x_col)$median
    b = fit$summary(variables = psi_y_col)$median
    rotation = atan(b / a)
  }
  
  # Append rotations and ellipse major/minor axes to original data
  df$a = a
  df$b = b
  df$alpha = rotation
  
  # If knots are supplied, extract their axes/rotation for plotting
  if (!is.null(knots)) {
    knots.df = merge(knots, df, by=c('x', 'y'), all.x=TRUE)
    # stopifnot(nrow(knots.df) == nrow(knots))
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
}


# Fitting model ---------------------------------------
# model = cmdstan_model('stan/aniso_process_axes_latent.stan')
# # Data
# data = prepare_data(epa.sample.df)
# fit = model$sample(data=data,
#                    parallel_chains=4,
#                    iter_warmup=1000,
#                    max_treedepth=15)
# 
# extract_ellipses(epa.sample.df, epa.df, fit, foci=F, scale_ellipses=7, num_ellipses=10)


# Fitting knot model ------------------------------------------------------
model.knots = cmdstan_model('stan/aniso_process_axes_latent_knots.stan')

# Data
data.knots = prepare_data(epa.df, num_knots=10, sample_pct=0.25, plot=T)

fit.knots = model.knots$sample(data=data.knots,
                                parallel_chains=4,
                                iter_warmup=1000,
                                max_treedepth=10)

data.knots.spatial.df = cbind(data.knots$spatial_locs, data.knots$y_spatial)
names(data.knots.spatial.df) = c('x', 'y', 'pm2_5')

data.knots.knots.df = data.knots$knot_locs

extract_ellipses(data.knots.spatial.df, 
                 fit=fit.knots, 
                 foci=F, 
                 scale_ellipses=15, 
                 knots=data.knots.knots.df, 
                 psi_suffix="_all")
extract_ellipses(epa.df, 
                 fit=fit.knots, 
                 foci=F, 
                 scale_ellipses=10, 
                 psi_suffix="_all")

