library(ggplot2)
library(dplyr)



# Load data ---------------------------------------------------------------
folder.5 = '2023-06-01_14_52_25'
folder.10 = '2023-06-01_15_30_49'
folder.15 = '2023-06-01_20_14_56'
folder.20 = '2023-06-02_07_25_49'

load_preds = function(folder) {
  # Hacky way to read in lots of predictions. Read the grid predictions first...
  preds.df = read.csv(file.path('results', folder, 'grid_preds.csv'))
  
  # Then find all other predictions...
  for (f in list.files(file.path('results', folder))[list.files(file.path('results', folder)) %>% endsWith('_preds.csv')]) {
    # But make sure not to read in the grid predictions again
    if (!startsWith(f, 'grid_preds')) {
      preds.df = rbind(preds.df, read.csv(file.path('results', folder, f)))
    }
  }
  
  return(preds.df)
}

load_knots = function(folder) {
  # Hacky way to read in lots of predictions. Read the grid predictions first...
  knots.df = read.csv(file.path('results', folder, 'grid_knots.csv'))
  
  # Then find all other predictions...
  for (f in list.files(file.path('results', folder))[list.files(file.path('results', folder)) %>% endsWith('_knots.csv')]) {
    # But make sure not to read in the grid predictions again
    if (!startsWith(f, 'grid_knots')) {
      new_knots.df = read.csv(file.path('results', folder, f))
      new_knots.df$pm2_5 = NULL
      knots.df = rbind(knots.df, new_knots.df)
    }
  }
  
  return(knots.df)
}

preds.5 = load_preds(folder.5)
preds.10 = load_preds(folder.10)
preds.15 = load_preds(folder.15)
preds.20 = load_preds(folder.20)

knots.5 = load_knots(folder.5)
knots.10 = load_knots(folder.10)
knots.15 = load_knots(folder.15)
knots.20 = load_knots(folder.20)

# Set levels in order to make graphs easier to read
set_levels = function(df, include_ellipsoid=FALSE) {
  l = c('grid', 'cover_design', 
        'sphere_entropy_nn_5_rm_0.75', 'sphere_entropy_nn_5_rm_1', 'sphere_entropy_nn_5_rm_1.25',
        'sphere_entropy_nn_7_rm_0.75', 'sphere_entropy_nn_7_rm_1', 'sphere_entropy_nn_7_rm_1.25',
        'sphere_entropy_nn_10_rm_0.75', 'sphere_entropy_nn_10_rm_1', 'sphere_entropy_nn_10_rm_1.25',
        'sphere_entropy_nn_15_rm_0.75', 'sphere_entropy_nn_15_rm_1', 'sphere_entropy_nn_15_rm_1.25',
        'sphere_entropy_nn_20_rm_0.75', 'sphere_entropy_nn_20_rm_1', 'sphere_entropy_nn_20_rm_1.25')
  if (include_ellipsoid) {
    l = c(l, c('ellipsoid_entropy_nnsphere_5_rmellipsoid_0.02', 'ellipsoid_entropy_nnsphere_5_rmellipsoid_0.04', 'ellipsoid_entropy_nnsphere_5_rmellipsoid_0.0666666666666667', 'ellipsoid_entropy_nnsphere_5_rmellipsoid_0.1',
                'ellipsoid_entropy_nnsphere_7_rmellipsoid_0.02', 'ellipsoid_entropy_nnsphere_7_rmellipsoid_0.04', 'ellipsoid_entropy_nnsphere_7_rmellipsoid_0.0666666666666667', 'ellipsoid_entropy_nnsphere_7_rmellipsoid_0.1',
                'ellipsoid_entropy_nnsphere_10_rmellipsoid_0.02', 'ellipsoid_entropy_nnsphere_10_rmellipsoid_0.04', 'ellipsoid_entropy_nnsphere_10_rmellipsoid_0.0666666666666667', 'ellipsoid_entropy_nnsphere_10_rmellipsoid_0.1',
                'ellipsoid_entropy_nnsphere_15_rmellipsoid_0.02', 'ellipsoid_entropy_nnsphere_15_rmellipsoid_0.04', 'ellipsoid_entropy_nnsphere_15_rmellipsoid_0.0666666666666667', 'ellipsoid_entropy_nnsphere_15_rmellipsoid_0.1'))
  }
  factor(df$type, levels = l)
  return(df)
}
preds.5 = set_levels(preds.5, include_ellipsoid = TRUE)
knots.5 = set_levels(knots.5, include_ellipsoid = TRUE)

preds.10 = set_levels(preds.10, include_ellipsoid = TRUE)
knots.10 = set_levels(knots.10, include_ellipsoid = TRUE)

preds.15 = set_levels(preds.15)
knots.15 = set_levels(knots.15)

preds.20 = set_levels(preds.20)
knots.20 = set_levels(knots.20)


# Plot knots --------------------------------------------------------------
plot_knots = function(preds, num_knots, knots) {
  g1 = ggplot(preds) +
    geom_point(aes(x=x, y=y, fill=pm2_5), color='black', pch=21, size=3) +
    scale_fill_gradient2() +
    geom_point(data=knots, aes(x=x, y=y), size=5, pch=18) +
    facet_wrap(vars(type)) +
    ggtitle(paste0(' (', num_knots, ' knots) - SD of predictions'))
  
  show(g1)
}


# Plot prediction errors --------------------------------------------------
plot_pred_error = function(preds, num_knots, knots) {
  g1 = ggplot(preds) +
    geom_point(aes(x=x, y=y, fill=sd_pred), color='black', pch=21, size=3) +
    geom_point(data=knots, aes(x=x, y=y), size=5, pch=18) +
    facet_wrap(vars(type)) +
    scale_fill_gradient(limits=c(0.75, 0.92), low='white', high='red') +
    ggtitle(paste0(' (', num_knots, ' knots) - SD of predictions'))
  
  g2 = ggplot(preds) +
    geom_point(aes(x=x, y=y, fill=pct_diff_median), color='black', pch=21, size=3) +
    geom_point(data=knots, aes(x=x, y=y), size=5, pch=18) +
    facet_wrap(vars(type)) +
    scale_fill_gradient2(low='blue', high='red') +
    ggtitle(paste0(' (', num_knots, ' knots) - Log(abs(% median diff))'))
  
  show(g1)
  show(g2)
}


# Summarize predictions ---------------------------------------------------
summarise_preds = function(preds) {
  print("Prediction summaries")
  preds %>%
    group_by(type) %>%
    summarise(r2=cor(pm2_5, median_pred)**2,
              mse=mean((pm2_5-median_pred)**2),
              median_sd=median(sd_pred)) %>%
    as.data.frame %>%
    arrange(mse) %>%
    print
  
  g = ggplot(preds, aes(x=pm2_5, y=median_pred, color=as.factor(type))) +
    geom_point() +
    geom_smooth(method='lm', se=F) +
    geom_abline(slope=1, intercept=0, linetype='dashed')
  
  show(g)
}


# Results -----------------------------------------------------------------
plot_knots(preds.5, 5, knots.5)
plot_pred_error(preds.5, 5, knots.5)
summarise_preds(preds.5)

plot_knots(preds.10, 10, knots.10)
plot_pred_error(preds.10, 10, knots.10)
summarise_preds(preds.10)

plot_knots(preds.15, 15, knots.15)
plot_pred_error(preds.15, 15, knots.15)
summarise_preds(preds.15)
