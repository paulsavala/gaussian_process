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



# Sample data (2d) -------------------------------------------------------------
# Get a sample of the rows as data
epa.sample.idx = sample(1:nrow(epa.df), 250)
epa.sample.df = epa.df[epa.sample.idx, ]

# Pick some random knots (just for testing purposes)
epa.knots.idx = sample(1:nrow(epa.sample.df), 15)
epa.sample.knots.df = epa.sample.df[epa.knots.idx, ]


# Fitting predictive process model ------------------------------------------------------
model = cmdstan_model('stan/predictive_process.stan')

data = list(N_knots=nrow(epa.sample.knots.df),
            knot_locs=epa.sample.knots.df[, c('x', 'y')],
            N_spatial=nrow(epa.sample.df),
            spatial_locs=epa.sample.df[, c('x', 'y')],
            y_spatial=epa.sample.df$pm2_5 %>% as.vector,
            return_predictions=1,
            return_log_likelihoods=1)

fit = model$sample(data=data,
                   parallel_chains=4,
                   iter_warmup=2000,
                   max_treedepth=10,
                   init=function() list(sigma_z=0.1,
                                        ell_z=0.1,
                                        sigma_interp=0.1,
                                        ell_interp=0.1,
                                        lambda_y=0.1))

# LOO
fit.knots$loo()
