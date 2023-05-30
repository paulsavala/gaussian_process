library(ggplot2)
library(dplyr)

# Spatial analysis (variograms)
library(sp)
library(gstat)

# Distance
library(proxy)

# Spectral analysis (periodogram)
library(TSA)

source('scripts/load_epa_data.R')


# Load data ---------------------------------------------------------------
df = load_epa_data()

ggplot(df) + 
  geom_point(aes(x=x, y=y, fill=pm2_5), size=3, color='black', pch=21) +
  scale_fill_gradient2(low='blue', high='red') +
  coord_fixed() + 
  ggtitle('Raw data')


# Fit polynomial to de-trend ----------------------------------------------
cubic.model = lm(data=df, formula=pm2_5 ~ I(x^3) + I(y^3) + I(x^2) + I(y^2) + I(x) + I(y) + 1)

cubic.preds = predict(cubic.model, newdata=df)

ggplot(df) + 
  geom_point(aes(x=x, y=y, fill=cubic.preds), size=3, color='black', pch=21) +
  scale_fill_gradient2(low='blue', high='red') +
  coord_fixed() + 
  ggtitle('Cubic trend data')

# De-trend
df$pm2_5_adj = df$pm2_5 - cubic.preds

ggplot(df) + 
  geom_point(aes(x=x, y=y, fill=pm2_5_adj), size=3, color='black', pch=21) +
  scale_fill_gradient2(low='blue', high='red') +
  coord_fixed() +
  ggtitle('Cubic trend data residuals')


# Convert to lattice ------------------------------------------------------
bins.x = seq(from=0, to=1, length.out=12)
bins.y = seq(from=0, to=1, length.out=12)
bins.grid = expand.grid(bins.x, bins.y)
names(bins.grid) = c('x', 'y')
bins.grid$idx = rownames(bins.grid) %>% as.numeric

# Find which lattice point is closest to each data point
dists = dist(df[, c('x', 'y')], bins.grid[, c('x', 'y')])
df$grid_idx = apply(dists, 1, which.min)
df.merged = merge(df, bins.grid, by.x='grid_idx', by.y='idx')
df.merged$grid_idx = NULL
names(df.merged) = c('x', 'y', 'pm2_5', 'pm2_5_adj', 'x_grid', 'y_grid')

df.grid = df.merged %>%
  group_by(x_grid, y_grid) %>%
  summarise(pm2_5 = median(pm2_5),
            pm2_5_adj = median(pm2_5_adj),
            n = n()) %>%
  as.data.frame

ggplot(df.grid) +
  geom_tile(aes(x=x_grid, y=y_grid, fill=pm2_5)) +
  scale_fill_gradient2(low='blue', high='red') +
  coord_fixed() +
  ggtitle('Binned PM2.5')

ggplot(df.grid) +
  geom_tile(aes(x=x_grid, y=y_grid, fill=pm2_5_adj)) +
  scale_fill_gradient2(low='blue', high='red') +
  coord_fixed() +
  ggtitle('Binned Adjusted PM2.5')


# Periodogram -------------------------------------------------------------
# Seems like I should compute this by hand. However, the observations need to
# occur over a regular grid. So I'll need to bin observations to make this work.
# To-do.


# Semivariogram -----------------------------------------------------------
# Convert to spatial points
coordinates(df) = ~x+y

# Variogram model for raw data (overall, and by direction)
v.raw = variogram(pm2_5 ~ 1, df)
v.raw.dir = variogram(pm2_5 ~ 1, df, alpha=c(0,45,90,135))

# Plot empirical variogram for raw data
plot(v.raw)
plot(v.raw.dir)

# Variogram model for adjusted data (overall, and by direction)
v.adj = variogram(pm2_5_adj ~ 1, df)
v.adj.dir = variogram(pm2_5_adj ~ 1, df, alpha=c(0,45,90,135))

# Plot empirical variogram for adjusted data
plot(v.adj)
plot(v.adj.dir)



