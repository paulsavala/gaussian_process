library(ggplot2)
library(dplyr)

library(fields)


# Data parameters ---------------------------------------------------------
N.KNOTS = 25
N.LATTICE.X = 25
N.LATTICE.Y = 25


# Load EPA data -----------------------------------------------------------
epa.data = readRDS('data/df_data_12list.RDS')
# Single snapshot
epa.df = epa.data[[1]]
epa.df$case_cntl = NULL
epa.df$year = NULL
names(epa.df) = c('station_id', 'x', 'y', 'pm2_5')
epa.df$station_id = NULL

# Remove duplicate locations
epa.df = epa.df[!duplicated(epa.df[, c('x', 'y')]), ]

# Convert to mean zero and standard deviation 1
epa.df$pm2_5 = epa.df$pm2_5 %>% scale

# Sample a small number of stations
epa.df.sample.idx = sample(1:nrow(epa.df), size=N.KNOTS)
epa.df.sample = epa.df[epa.df.sample.idx, ]

# Train-test split data
train_size = as.integer(0.7 * nrow(epa.df))
epa.train.idx = sample(1:nrow(epa.df), size=train_size)
epa.test.idx = setdiff(1:nrow(epa.df), epa.train.idx)
epa.train.df = epa.df[epa.train.idx, ]
epa.test.df = epa.df[epa.test.idx, ]

# Sample a small number of stations
epa.sample.idx = sample(1:nrow(epa.train.df), size=N.KNOTS)
epa.sample.df = epa.train.df[epa.sample.idx, ]

# Get X and y matrices
X.obs = epa.train.df[, names(epa.train.df) != 'pm2_5']
X.oos = epa.test.df[, names(epa.test.df) != 'pm2_5']
y.obs = epa.train.df$pm2_5

# Make a lattice of points across the US for predictions
x1 = seq(min(epa.df$x), max(epa.df$x), length=N.LATTICE.X)
x2 = seq(min(epa.df$y), max(epa.df$y), length=N.LATTICE.X)
X.lattice = expand.grid(x1, x2)
df.lattice = data.frame(X.lattice)
names(df.lattice) = c('x', 'y')

# Get X and y values for sample of stations
X.stations = epa.df.sample[, c('x', 'y')] %>% as.matrix
y.stations = epa.df.sample$pm2_5 %>% as.matrix(ncol=1)
df.stations = data.frame(X=X.stations, y=y.stations)

# Plot locations
# ggplot(epa.df) +
#   geom_point(aes(x=x, y=y, color=pm2_5)) +
#   scale_color_gradient2(low='blue', high='red', midpoint=0) +
#   geom_point(data=epa.df.sample, aes(x=x, y=y), shape=4, size=5)


# Fit stationary GP -------------------------------------------------------
gp1 = spatialProcess(x=X.obs, y=y.obs)

