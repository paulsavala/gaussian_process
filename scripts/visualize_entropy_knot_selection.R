# Setup -----------------------------------------------------------------------------------------------------------
# Imports (Graphing)
library(ggplot2)  # graphing; ggplot() + ...
library(ggforce)  # graphing circles; geom_circle()
# library(plotly)   # 3d scatter plots; plot_ly()
# library(ggdark)   # dark theme; dark_theme_gray()

# Imports (Modeling)
library(mgcv)     # evaluating knots via GAMs; gam(), vis.gam()
library(fields)

# Imports
library(dplyr)    # piping and QoL functions; %>%
library(scales)   # normalizing lists between 0 and 1; rescale()
library(stringr)  # manipulating strings

# Knot Evaluation
source("scripts/fuentes_GP_model.R")
source("scripts/fuentes_entropy_knot_selection.R")
source("scripts/load_epa_data.R")
source("scripts/visualize_gp_predictions.R")
source("scripts/ka_v5_functions.R")



# Parameters --------------------------------------------------------------
N.KNOTS = 20


# Data loading ------------------------------------------------------------
epa.df = load_epa_data()

# # Knot selection: Entropy --------------------------------------------------
# knots.entropy = vkr_base(epa.df, list(n_neighbors=5, radius_mult=0.5, max_knots=N.KNOTS, cols_to_sort=c("entropy")))
# knots.entropy = knots.entropy[, c('x', 'y')]
# 
# ggplot(epa.df) +
#   geom_point(aes(x=x, y=y, color=pm2_5), size=3) +
#   scale_color_gradient2(low='blue', high='red', midpoint=0) +
#   geom_point(data=knots.entropy, aes(x=x, y=y), shape=4, size=5) +
#   ggtitle('Entropy knots (10 neighbors, radius multiplier = 0.5)')


# California (roughly) ----------------------------------------------------
# Filter to approximately California
epa.ca.df = epa.df %>%
  filter(x < 0.2, y > 0.25, y < 0.75) %>%
  filter((x<0.125) | (y < 0.4))

epa.ca.df = epa.ca.df[, c('x', 'y', 'pm2_5')]

ggplot(epa.ca.df) +
  geom_point(aes(x=x, y=y, color=pm2_5), size=3) +
  scale_color_gradient2(low='blue', high='red', midpoint=0) +
  ggtitle('California')

# Fit knots on just California
n_neighbors = 10
radius_mult = 2

knots.ca.entropy = vkr_base(epa.ca.df, list(n_neighbors=n_neighbors, 
                                            radius_mult=radius_mult, 
                                            max_knots=N.KNOTS, 
                                            cols_to_sort=c("entropy")))
knots.ca.entropy = knots.ca.entropy[, c('x', 'y', 'radius')]

ggplot(epa.ca.df) +
  geom_point(aes(x=x, y=y, color=pm2_5), size=3) +
  scale_color_gradient2(low='blue', high='red', midpoint=0) +
  geom_point(data=knots.ca.entropy, aes(x=x, y=y), shape=4, size=5, stroke=2) +
  # geom_circle(data=knots.ca.entropy, aes(x0=x, y0=y, r=radius*radius_mult)) +
  ggtitle(paste0('Entropy knots (', n_neighbors, ' neighbors, radius multiplier = ', radius_mult, ')'))









