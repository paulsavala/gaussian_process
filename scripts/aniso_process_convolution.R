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

# Knot Evaluation
source("scripts/load_epa_data.R")



# Parameters --------------------------------------------------------------
A = 1 # Area of each ellipse (fixed)


# Data loading ------------------------------------------------------------
epa.df = load_epa_data()

# California (roughly) ----------------------------------------------------
# Filter to approximately California
epa.ca.df = epa.df %>%
  dplyr::filter(x < 0.2, y > 0.25, y < 0.75) %>%
  dplyr::filter((x < 0.125) | (y < 0.4))

epa.ca.df = epa.ca.df[, c('x', 'y', 'pm2_5')]

# Sample for easy testing
epa.ca.sample.idx = sample(1:nrow(epa.ca.df), 20)
epa.ca.sample.df = epa.ca.df[epa.ca.sample.idx, ]

# Plot CA and sample
ggplot(epa.ca.df) +
  geom_point(aes(x=x, y=y, color=pm2_5), size=3) +
  geom_point(data=epa.ca.sample.df, aes(x=x, y=y), size=3, shape=4) +
  scale_color_gradient2(low='blue', high='red', midpoint=0) +
  ggtitle('California')



# Stan model fitting ------------------------------------------------------
# Spatial process sites and data
spatial_locs = epa.ca.sample.df[, c('x', 'y')]
y_spatial = epa.ca.sample.df$pm2_5 %>% as.vector
N_spatial = nrow(spatial_locs)

# Latent process sites and data
latent.ca.idx = sample(1:nrow(epa.ca.sample.df), 5)
latent.ca.df = epa.ca.sample.df[latent.ca.idx, ]
latent_locs = latent.ca.df[, c('x', 'y')]
N_latent = nrow(latent_locs)

# Data
data = list(N_spatial=N_spatial,
            spatial_locs=spatial_locs,
            y_spatial=y_spatial,
            N_latent=N_latent,
            latent_locs=latent_locs)

# Model
model = cmdstan_model('stan/aniso_process_convolution.stan')
fit = model$sample(data=data,
                   parallel_chains=4,
                   iter_warmup=4000,
                   max_treedepth=15)

# Data
# int<lower=0> N_spatial; // Number of sites at which spatial process is measured
# array[2] vector[N_spatial] spatial_locs; // x-y coordinates of spatial process
# vector[N_spatial] y_spatial; // Measured value of spatial process at each site
# int<lower=0> N_latent; // Number of sites at which latent process is modeled
# array[2] vector[N_latent] latent_locs; // x-y coordinates of latent process

