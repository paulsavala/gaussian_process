# Setup -----------------------------------------------------------------------------------------------------------
                                                            # ctrl+L (clear console output)

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
RADIUS = 0.2
N.KNOTS = 25


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
  geom_point(data=grid.df, aes(x=x, y=y), alpha=0.25, size=1) +
  ggtitle('Training data and grid')


# Knot selection methods --------------------------------------------------
#### Knot Selection: Cover.Design ####
make_knots = function(df, N_k, ...) {
  # remove duplicate locations
  df = df[!duplicated(df), ]
  design = fields::cover.design(R=df, nd=N_k, ...)
  knots_df = as.data.frame(design$design)
  return(knots_df)
}

knots.cd = make_knots(grid.df, N.KNOTS)

ggplot(epa.df) +
  geom_point(aes(x=x, y=y, color=pm2_5, size=2)) +
  scale_color_gradient2(low='blue', high='red', midpoint=0) +
  geom_point(data=knots.cd, aes(x=x, y=y), shape=4, size=5) +
  ggtitle('Training data and CD knots')


#### Knot Selection: Entropy Maximization ####

knots.entropy = vkr_base(df.train, list(n_neighbors=10, radius_mult=1, max_knots=N.KNOTS, cols_to_sort=c("entropy")))
knots.entropy = knots.entropy[, c('x', 'y')]
# vkr.gs = vkr_gs(df.points=df.train, seq.nn=1:10, seq.rm=seq(from=2, to=3, by=.25), n.knots=25, cols.to.sort=c("entropy"), gam.k=3)

ggplot(epa.df) +
  geom_point(aes(x=x, y=y, color=pm2_5, size=2)) +
  scale_color_gradient2(low='blue', high='red', midpoint=0) +
  geom_point(data=knots.entropy, aes(x=x, y=y), shape=4, size=5) +
  ggtitle('Training data and Entropy knots')


#### Knot Selection: Fuentes spatial process ####
knots.fuentes = fit.fuentes.knots(N.KNOTS, X.train, y.train, n.knot_sets=15)

ggplot(epa.df) +
  geom_point(aes(x=x, y=y, color=pm2_5, size=2)) +
  scale_color_gradient2(low='blue', high='red', midpoint=0) +
  geom_point(data=knots.fuentes, aes(x=x, y=y), shape=4, size=5) +
  ggtitle('Training data and Fuentes knots')

# Knot Evaluation ------------------------------------------------------------------------------------
# Fit Fuentes spatial process model using each set of knots
gp.cd.params = fit.fuentes.gp(knots.cd, X.train, y.train)
gp.entropy.params = fit.fuentes.gp(knots.entropy, X.train, y.train)
gp.fuentes.params = fit.fuentes.gp(knots.fuentes, X.train, y.train)

# Make predictions on the test set
gp.cd.preds = predict.fuentes.gp(knots.cd, X.test, gp.cd.params)
gp.entropy.preds = predict.fuentes.gp(knots.entropy, X.test, gp.entropy.params)
gp.fuentes.preds = predict.fuentes.gp(knots.fuentes, X.test, gp.fuentes.params)

# Compute MSE
gp.cd.mse = sum((gp.cd.preds$median_pred - y.test)**2) / length(y.test)
gp.entropy.mse = sum((gp.entropy.preds$median_pred - y.test)**2) / length(y.test)
gp.fuentes.mse = sum((gp.fuentes.preds$median_pred - y.test)**2) / length(y.test)

# Compute median standard error
gp.cd.se = median(gp.cd.preds$se)
gp.entropy.se = median(gp.entropy.preds$se)
gp.fuentes.se = median(gp.fuentes.preds$se)

# Print results in a nicely formatted table
line1=paste0(N.KNOTS, " knots")
line2=paste0("\nCover Design:\n\tMSE = ", gp.cd.mse %>% round(4), " --- Median standard error = ", gp.cd.se %>% round(4))
line3=paste0("\nEntropy:\n\tMSE = ", gp.entropy.mse %>% round(4), " --- Median standard error = ", gp.entropy.se %>% round(4))
line4=paste0("\nFuentes:\n\tMSE = ", gp.fuentes.mse %>% round(4), " --- Median standard error = ", gp.fuentes.se %>% round(4))

cat(line1, line2, line3, line4)

# eval_knots = function(df.points, df.knots, pred.knots.name="default_pred_knot_name", k=3, vis=F) {
#   # gam fitting
#   gam.eval = gam(
#     pm2_5 ~ s(x, y, k=k, bs="gp"),  # , k=12, bs="so"
#     data=df.knots,
#     method="REML",
#     family=gaussian
#   )
#   # gam visualization
#   if (vis == T) { vis.gam(gam.eval, type="response", plot.type="contour", color="gray") }
#   # predictions
#   df.points[, pred.knots.name] = predict(gam.eval, newdata=df.points)
#   # mean squared error; lower is better
#   # mse = (df.points$signal - df.points[, pred.knots.name])^2 %>% mean
#   # print(paste0("Mean Squared Error (lower better): ", mse))
#   # correlation squared; higher is better
#   # cor2 = cor(df.points$signal, df.points[, pred.knots.name])^2
#   # print(paste0("Correlation Squared (higher better): ", cor2))
#   df.points
# }
# 
# eval_knots_mse = function(df.points, df.knots, gam.k=3) {
#   gam.eval = gam(signal ~ te(x, y, k=gam.k, bs="gp"),data=df.knots,method="REML", family=gaussian)
#   df.points[, "pkn"] = predict(gam.eval, newdata=df.points)
#   (df.points$signal - df.points[, "pkn"])^2 %>% mean  # mse
# }
# 
# eval_knots_metrics = function(df.points.evaled) {
#   result.metrics = list()
#   n.knot.types = (names(df.points.evaled) %>% length) - 3
#   for(i.knot.type in 1:n.knot.types) {
#     i.name = names(df.points.evaled)[3+i.knot.type]
#     result.metrics[paste0(i.name, "_mse")] = (df.points.evaled$signal - df.points.evaled[, i.name])^2 %>% mean
#     result.metrics[paste0(i.name, "_cor2")] = cor(df.points.evaled$signal, df.points.evaled[, i.name])^2
#   }
#   result.metrics
# }

# Testing ---------------------------------------------------------------------------------------------------------

gam_k = 3
test.data = eval_knots(test.data, test.data, "all_data_preds", gam_k, vis=T)
test.data = eval_knots(test.data, knots.entropy, "entropy_knot_preds", gam_k, vis=T)
test.data = eval_knots(test.data, knots.cover.design, "cd_knot_preds", gam_k, vis=T)
test.knot.metrics = eval_knots_metrics(test.data)
tkm_results = function(kmetrics) {
  print(paste0("all:          ", kmetrics$all_data_preds_mse))
  print(paste0("entropy:      ", kmetrics$entropy_knot_preds_mse))
  print(paste0("cover.design: ", kmetrics$cd_knot_preds_mse))
}
tkm_results(test.knot.metrics)

alpha = .3
size = 3
ggplot(test.data) + geom_abline(slope=1, intercept=0) + theme_dark() +
  geom_point(aes(x=signal, y=all_data_preds, color="all"), alpha=alpha, size=size) +
  geom_point(aes(x=signal, y=entropy_knot_preds, color="entropy"), alpha=alpha, size=size)
ggplot(test.data) + geom_abline(slope=1, intercept=0) + theme_dark() +
  geom_point(aes(x=signal, y=all_data_preds, color="all"), alpha=alpha, size=size) +
  geom_point(aes(x=signal, y=cd_knot_preds, color="cd"), alpha=alpha, size=size)

# Fuentes Testing -------------------------------------------------------------------------------------------------
source("fuentes_GP_model.R")

# Training data with X columns x and y, and y column pm2_5
td = test.data[c("x", "y", "signal")]
names(td) = c("x", "y", "pm2_5")
td.train = td %>% sample_frac(.7)
td.test = td %>% anti_join(td.train)

td.X.train = td.train[, c("x", "y")]
td.y.train = td.train[, c("pm2_5")]

td.X.test = td.test[, c("x", "y")]
td.y.test = td.test[, c("pm2_5")]

# entropy fitting
params.entropy = fit.fuentes.gp(knots.entropy, td.X.train, td.y.train)
preds.df.entropy = predict.fuentes.gp(test.knots, td.X.test, params.entropy)

# Automation ------------------------------------------------------------------------------------------------------

cd_mse = function(n.knots, n.iter) {
  mses = c()
  for (i in 1:n.iter) {
    knots = generate_cd_knots(test.data, n.knots)
    # paste0("Currently Doing: i=", i) %>% print
    mses = c(mses, eval_knots_mse(test.data, knots))
  }
  mses
}

# TODO: create train-test split for more accurate MSEs
# TODO: add random knot selection and Fuentes knot seleciton
# TODO: use Fuentes model and compare to GAMs

auto_mse = function(seq.n.knots) {
  results = data.frame(nk=integer(), entropy=double(), cd=double())
  for (i.n.knots in seq.n.knots) {
    paste0("Currently Doing: n_knots=", i.n.knots) %>% print
    # gs.ent = gs_entropy(seq.nn=1:5, seq.rm=seq(from=1, to=3, by=.25), i.n.knots)
    gs.ent = vkr_gs(df.points=test.data, seq.nn=1:10, seq.rm=seq(from=2, to=3, by=.25), n.knots=n.knots, cols.to.sort=c("entropy"), gam.k=3)
    mse.entropy = gs.ent$grid$mse %>% min
    mse.cd = cd_mse(n.knots=i.n.knots, n.iter=10) %>% mean
    results = results %>% rbind(data.frame(n_k=i.n.knots, mse_entropy=mse.entropy, mse_cd=mse.cd))
    print(results)
  }
  results
}
test.mse = auto_mse(seq(from=10, to=20, by=10))
test.mse

ggplot(test.mse) + theme_dark() +
  labs(title="MSE for each knot selection method", subtitle="subtitle", caption="caption") +
  geom_line(aes(x=n_k, y=mse_entropy, color="mse_entropy"), size=2) + geom_line(aes(x=n_k, y=mse_cd, color="mse_cd"), size=2)

ggplot() + theme_dark() +
  labs(title="MSE for each knot selection method", subtitle="subtitle", caption="caption") +
  geom_line(data=gs.entropy, aes(x=nn, y=mse, color=as.factor(rm)), size=2) + geom_line(data=test.mse, aes(x=n_k, y=mse_cd, color="mse_cd"), size=2)

