# Setup -----------------------------------------------------------------------------------------------------------
# Imports (Graphing)
library(ggplot2)  # graphing; ggplot() + ...
library(ggforce)  # graphing circles; geom_circle()

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
RADIUS = 0.2

# Data loading ------------------------------------------------------------
epa.df = load_epa_data()
tts = train_test_split(epa.df)

X.train = tts$X.train
y.train = tts$y.train

X.test = tts$X.test
y.test = tts$y.test

df.train = tts$df.train
df.test = tts$df.test

grid.df = form_grid(df.train, RADIUS)

ggplot(epa.df) +
  geom_point(aes(x=x, y=y, color=pm2_5, size=2)) +
  scale_color_gradient2(low='blue', high='red', midpoint=0) +
  ggtitle('Training data')

for (N.KNOTS in c(25, 50, 75)) {
  # Knot selection methods --------------------------------------------------
  #### Knot Selection: Cover.Design ####
  make_knots = function(df, N_k, ...) {
    # remove duplicate locations
    df = df[!duplicated(df), ]
    design = fields::cover.design(R=df, nd=N_k, ...)
    knots_df = as.data.frame(design$design)
    return(knots_df)
  }
  
  knots.cd = make_knots(X.train, N.KNOTS)
  
  # ggplot(epa.df) +
  #   geom_point(aes(x=x, y=y, color=pm2_5, size=2)) +
  #   scale_color_gradient2(low='blue', high='red', midpoint=0) +
  #   geom_point(data=grid.df, aes(x=x, y=y), alpha=0.25, size=1) +
  #   geom_point(data=knots.cd, aes(x=x, y=y), shape=4, size=5) +
  #   ggtitle('Training data, grid and CD knots')
  
  
  #### Knot Selection: Entropy Maximization ####
  n_neighbors = 10
  radius_mult = 1
  knots.entropy = vkr_base(df.train, 
                           list(n_neighbors=n_neighbors, 
                                radius_mult=radius_mult, 
                                max_knots=N.KNOTS, 
                                cols_to_sort=c("entropy")))
  knots.entropy = knots.entropy[, c('x', 'y')]
  # vkr.gs = vkr_gs(df.points=df.train, seq.nn=1:10, seq.rm=seq(from=2, to=3, by=.25), n.knots=25, cols.to.sort=c("entropy"), gam.k=3)
  
  # ggplot(epa.df) +
  #   geom_point(aes(x=x, y=y, color=pm2_5, size=2)) +
  #   scale_color_gradient2(low='blue', high='red', midpoint=0) +
  #   geom_point(data=grid.df, aes(x=x, y=y), alpha=0.25, size=1) +
  #   geom_point(data=knots.entropy, aes(x=x, y=y), shape=4, size=5) +
  #   ggtitle('Training data, grid and Entropy knots')
  
  
  #### Knot Selection: Fuentes spatial process ####
  knots.fuentes = fit.fuentes.knots(N.KNOTS, X.train, y.train, n.knot_sets=10)
  # 
  # ggplot(epa.df) +
  #   geom_point(aes(x=x, y=y, color=pm2_5, size=2)) +
  #   scale_color_gradient2(low='blue', high='red', midpoint=0) +
  #   geom_point(data=grid.df, aes(x=x, y=y), alpha=0.25, size=1) +
  #   geom_point(data=knots.fuentes, aes(x=x, y=y), shape=4, size=5) +
  #   ggtitle('Training data, grid and Fuentes knots')
  
  # Knot Evaluation ------------------------------------------------------------------------------------
  # Fit Fuentes spatial process model using each set of knots
  gp.cd.params = fit.fuentes.gp(knots.cd, X.train, y.train)
  gp.entropy.params = fit.fuentes.gp(knots.entropy, X.train, y.train)
  # gp.fuentes.params = fit.fuentes.gp(knots.fuentes, X.train, y.train)
  
  # Make predictions on the test set
  gp.cd.preds = predict.fuentes.gp(knots.cd, X.test, gp.cd.params)
  gp.entropy.preds = predict.fuentes.gp(knots.entropy, X.test, gp.entropy.params)
  # gp.fuentes.preds = predict.fuentes.gp(knots.fuentes, X.test, gp.fuentes.params)
  
  # Compute MSE
  gp.cd.mse = sum((gp.cd.preds$median_pred - y.test)**2) / length(y.test)
  gp.entropy.mse = sum((gp.entropy.preds$median_pred - y.test)**2) / length(y.test)
  # gp.fuentes.mse = sum((gp.fuentes.preds$median_pred - y.test)**2) / length(y.test)
  
  # Compute median standard error
  gp.cd.se = median(gp.cd.preds$se)
  gp.entropy.se = median(gp.entropy.preds$se)
  # gp.fuentes.se = median(gp.fuentes.preds$se)
  
  # Print results in a nicely formatted table
  line1=paste0(N.KNOTS, " knots")
  line2=paste0("\nCover Design:\n\tMSE = ", gp.cd.mse %>% round(4), " --- Median standard error = ", gp.cd.se %>% round(4))
  line3=paste0("\nEntropy:\n\tMSE = ", gp.entropy.mse %>% round(4), " --- Median standard error = ", gp.entropy.se %>% round(4))
  # line4=paste0("\nFuentes:\n\tMSE = ", gp.fuentes.mse %>% round(4), " --- Median standard error = ", gp.fuentes.se %>% round(4))
  
  cat(line1, line2, line3) 
  # cat(line1, line2, line3, line4) 
}

