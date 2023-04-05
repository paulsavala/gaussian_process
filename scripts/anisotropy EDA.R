# Setup -----------------------------------------------------------------------------------------------------------
# Imports (Graphing)
library(ggplot2)  # graphing; ggplot() + ...
library(ggforce)  # graphing circles; geom_circle()
# library(plotly)   # 3d scatter plots; plot_ly()
# library(ggdark)   # dark theme; dark_theme_gray()

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


# Data loading ------------------------------------------------------------
epa.df = load_epa_data()

# Take a sample of the data for easier computation
epa.df = epa.df[sample(1:nrow(epa.df), size=500), ]

# Filter data to a smaller set for easier computations
# epa.df = epa.df %>%
#   filter(x > 0.5)

ggplot(epa.df) + 
  geom_point(aes(x=x, y=y, fill=pm2_5), pch=21, color='black', size=2) +
  scale_fill_gradient2(low='blue', high='red', midpoint=0)



# Empirical semivariogram contour (ESC) plots -----------------------------
# Define rectangular bins to hold points based on displacement
bins.x.width = 0.1
bins.y.width = 0.1
# x displaments can be positive or negative, with a maximum possible value
# of +- max(epa.df$x) (assuming all coordinates are positive)
bins.x = seq(-max(epa.df$x), max(epa.df$x), by=bins.x.width)
# y displacements are restricted to be positive, so no such issue
bins.y = seq(0, max(epa.df$y), by=bins.y.width)
# Create data frames for bin cutoffs for easy merging later
bins.x.df = data.frame(cutoff=bins.x, midpoint = bins.x - bins.x.width / 2, idx=1:length(bins.x))
bins.y.df = data.frame(cutoff=bins.y, midpoint = bins.y - bins.y.width / 2, idx=1:length(bins.y))

# Create a long data frame to house bin indices for each pair of points
bins.epa.df = expand.grid(1:nrow(epa.df), 1:nrow(epa.df)) %>% as.data.frame
names(bins.epa.df) = c('row.epa1', 'row.epa2')
bins.epa.df$x.idx = NA
bins.epa.df$y.idx = NA
# Filter to just points with second coordinate after the first
bins.epa.df = bins.epa.df[bins.epa.df$row.epa2 > bins.epa.df$row.epa1, ]
rownames(bins.epa.df) = 1:nrow(bins.epa.df)

# Go through all pairs of points and...
for (i in 1:nrow(epa.df)) {
  for (j in i:nrow(epa.df)) {
    row.i = epa.df[i, ]
    row.j = epa.df[j, ]
    
    # 1. Calculate the displacement in the x and y axes
    h.x = row.i$x - row.j$x
    h.y = row.i$y - row.j$y
    
    if (h.y < 0) {
      h.x = -h.x
      h.y = -h.y
    }
    
    # 2. Put into rectangular bins of distances
    x.idx = ((h.x < bins.x) == TRUE) %>% which %>% min
    y.idx = ((h.y < bins.y) == TRUE) %>% which %>% min
    if (x.idx == Inf) { x.idx = 1 }
    if (y.idx == Inf) { y.idx = 1 }
    bins.epa.df[(bins.epa.df$row.epa1 == i) & (bins.epa.df$row.epa2 == j), c('x.idx', 'y.idx')] = c(x.idx, y.idx)
    }
}

# 3. Compute the within-bin variance
# Merge in the pm2.5 values for each point
epa.df$idx = rownames(epa.df) %>% as.numeric
bins.epa.df = merge(bins.epa.df, epa.df[, c('idx', 'pm2_5')], by.x='row.epa1', by.y='idx')
names(bins.epa.df)[names(bins.epa.df) == 'pm2_5'] = 'pm2_5.1'
bins.epa.df = merge(bins.epa.df, epa.df[, c('idx', 'pm2_5')], by.x='row.epa2', by.y='idx')
names(bins.epa.df)[names(bins.epa.df) == 'pm2_5'] = 'pm2_5.2'

bins.var.df = bins.epa.df %>%
  group_by(x.idx, y.idx) %>%
  summarise(bucket_var=sum((pm2_5.1 - pm2_5.2)**2) / (2*n())) %>%
  as.data.frame

# 4. Plot a heatmap of the within-bin variance
bins.epa.df = merge(bins.var.df, bins.x.df[, c('idx', 'midpoint')], 
                    by.x='x.idx',
                    by.y='idx')
names(bins.epa.df)[names(bins.epa.df) == 'midpoint'] = 'x.midpoint'

bins.epa.df = merge(bins.epa.df, bins.x.df[, c('idx', 'midpoint')], 
                    by.x='y.idx',
                    by.y='idx')
names(bins.epa.df)[names(bins.epa.df) == 'midpoint'] = 'y.midpoint'

ggplot(bins.epa.df) +
  geom_tile(aes(x=x.midpoint, y=y.midpoint, fill=bucket_var)) +
  scale_fill_gradient(low='white', high='blue') +
  coord_fixed() +
  ggtitle('Empirical semivariogram') +
  xlab('Displacement in x') +
  ylab('Displacement in y')
