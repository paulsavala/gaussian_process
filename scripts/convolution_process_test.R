library(gsignal) # conv2 - 2d convolution function
library(cmdstanr)
library(mvtnorm)
library(plgp) # distance function


# Generate grid (grid process) --------------------------------------------
# Define a fine lattice of points
grid.num_x = 20
grid.num_y = 20
grid.mat = matrix(rnorm(grid.num_x*grid.num_y), grid.num_x, grid.num_y)
grid.df = melt(grid.mat)
names(grid.df) = c('x', 'y', 'value')


# Generate spatial process locations --------------------------------------
# In general the spatial process will be unrelated to the grid process. But
# for the sake for experimentation, I'm defining the spatial process to be
# the grid process with some jitter.

# Jitter the spatial process locations to get more realistic locations
spatial.jitter.dist = (grid.df[2, 'x'] - grid.df[1, 'x']) / 3
spatial.jitter.x = rnorm(nrow(grid.df), mean=0, sd=spatial.jitter.dist)
spatial.jitter.y = rnorm(nrow(grid.df), mean=0, sd=spatial.jitter.dist)

spatial.df = as.data.frame(grid.df)
# spatial.df$x = spatial.df$x + spatial.jitter.x
# spatial.df$y = spatial.df$y + spatial.jitter.y

# Plot locations and grid
ggplot() +
  geom_point(data=grid.df, aes(x=x, y=y), color='black') +
  geom_point(data=spatial.df, aes(x=x, y=y), color='blue', size=3) +
  coord_fixed()


# Generate spatial process data for testing (GP) -------------------------
# Scale and rotate points to make it anisotropic
rot.mat = function(a, b, theta, plot_density=FALSE) {
  stretch = matrix(c(a, 0, 0, b), nrow=2, ncol=2, byrow=T)
  rotate = matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow=2, ncol=2, byrow=T)
  out = rotate %*% stretch
  
  if (plot_density) {
    g = dmvnorm(grid.df[, c('x', 'y')], mean = c(10, 10), sigma = solve(out %*% t(out))) %>%
      matrix(nrow=sqrt(nrow(grid.df)), ncol=sqrt(nrow(grid.df)), byrow=T) %>%
      melt %>%
      ggplot() +
      geom_tile(aes(x=Var1, y=Var2, fill=value)) +
      coord_fixed()
    
    show(g)
  }
  
  return(out)
}

# (Optional) Scale and rotate points to induce direction. Set to 1, 1, 0 to do nothing.
scale_rotate.mat = rot.mat(1, 1, 0, plot_density = F)

xy_transform = scale_rotate.mat %*% t(spatial.df[, c('x', 'y')])
spatial.transform.df = xy_transform %>% t %>% as.data.frame
names(spatial.transform.df) = c('x', 'y')

# Compute distance between points
D = distance(spatial.transform.df[, c('x', 'y')])

# Gaussian covariance matrix
eps = sqrt(.Machine$double.eps)
n = nrow(D)
Sigma = exp(-D) + diag(eps, n)

# Generate from GP with covariance Sigma and mean zero
y = rmvnorm(1, sigma=Sigma)

# Matrix is used for model fitting
y.mat = matrix(y, nrow=grid.num_x, ncol=grid.num_y, byrow=T)

# Dataframe is used for plotting and analysis
spatial.df$value = t(y)

# Visualize GP-based data
ggplot(spatial.df) +
  # geom_point(data=grid.df, aes(x=x, y=y), color='black') +
  # geom_point(aes(x=x, y=y, fill=value), size=3, pch=21, color='black') +
  geom_tile(aes(x=x, y=y, fill=value)) +
  scale_fill_gradient2(low='blue', high='red', midpoint=0) +
  coord_fixed()


# Define convolution kernels -----------------------------------------------

# Basis locations (toy example)
basis.locs = data.frame(x=c(5, 10, 15), 
                        y=c(5, 10, 15))

n_basis = nrow(basis.locs)

# Basis inverse square root covariance matrices (test values)
# Function for easily making inverse square root covariance matrices
rot.mat = function(a, b, theta) {
  stretch = matrix(c(a, 0, 0, b), nrow=2, ncol=2, byrow=T)
  rotate = matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow=2, ncol=2, byrow=T)
  out = rotate %*% stretch
  return(out)
}

conv.kernels.basis = list(rot.mat(0.5, 1, 0),
                    rot.mat(1, 0.5, 0),
                    rot.mat(0.5, 1, -pi/4))

# Plot kernels and associated densities (sanity check)
view.densities = matrix(0, nrow=grid.num_x, ncol=grid.num_y)
for (i in 1:n_basis) {
  ck = conv.kernels.basis[[i]]
  bl = basis.locs[i, c('x', 'y')]
  view.densities = view.densities + 
    dmvnorm(grid.df[, c('x', 'y')], mean = c(bl$x, bl$y), sigma = solve(ck %*% t(ck)))
}

view.densities %>%
  melt %>%
  ggplot() +
    geom_tile(aes(x=Var1, y=Var2, fill=value)) +
    coord_fixed()

# Interpolate inverse square root covariance matrices to all locations
# Holds interpolated kernels at all locations
conv.kernels.all = list()

# Distance weight function
w_k = function(s, basis.loc) {
  exp(-0.5 * norm(s - basis.loc, type='2'))
}

# Interpolate density kernels to all other locations using equation on pg 179
for (s.idx in 1:nrow(spatial.df)) {
  # Get spatial process location
  s = spatial.df[s.idx, c('x', 'y')]
  
  # Keep track of sum of w (for normalizing)
  w.sum = 0
  
  # Interpolated kernel at location s
  ck = matrix(0, nrow=2, ncol=2)
  for (basis.loc.idx in 1:nrow(basis.locs)) {
    # Get basis location
    basis.loc = basis.locs[basis.loc.idx, c('x', 'y')]
    
    # Calculate basis location weight
    w = w_k(s, basis.loc)
    
    # Update sum of w
    w.sum = w.sum + w
    
    # Add weighted value of convolution from basis location
    ck = ck + w * conv.kernels.basis[[basis.loc.idx]]
    
    # Null out to avoid floating variables
    basis.loc = NULL
  }
  
  # Store interpolated kernel density
  conv.kernels.all[[s.idx]] = ck / w.sum
}

# Convert inverse square root covariance matrices to actual covariance matrices
conv.kernels.all.cov = list()

for (i in 1:nrow(spatial.df)) {
  # Extract values
  mu = spatial.df[i, c('x', 'y')]
  kernel = conv.kernels.all[[i]]
  
  # Convert kernel matrix to covariance matrix (inverse square root)
  conv.kernels.all.cov[[i]] = solve(kernel %*% t(kernel))
}

# Generate kernel densities using covariance matrices
conv.kernel.densities = list()

for (i in 1:nrow(spatial.df)) {
  # Basis location
  s = spatial.df[i, c('x', 'y')]
  
  # Basis kernel density (setting mean to basis location is equivalent to mean (0, 0)
  # and evaluating at displacement. Doing it this one to save a computational step).
  density.s = dmvnorm(grid.df[, c('x', 'y')], mean = c(s$x, s$y), sigma = conv.kernels.all.cov[[i]])
  
  # Reshape the density vector into a matrix
  conv.kernel.densities[[i]] = density.s %>%
    matrix(nrow=grid.num_x, ncol=grid.num_y, byrow=TRUE)
}

# # Plot densities (sanity check)
# conv.kernel.densities.df = data.frame(matrix(nrow=0, ncol=3))
# names(conv.kernel.densities.df) = c('x', 'y', 'value')
# 
# for (i in 1:n_basis) {
#   conv.kernel.densities.i.df = melt(conv.kernel.densities[[i]])
#   names(conv.kernel.densities.i.df) = c('x', 'y', 'value')
#   conv.kernel.densities.i.df$idx = i
#   conv.kernel.densities.df = rbind(conv.kernel.densities.df, conv.kernel.densities.i.df)
#   conv.kernel.densities.i.df = NULL
# }
# 
# ggplot(conv.kernel.densities.df) +
#   geom_tile(aes(x=x, y=y, fill=value)) +
#   facet_wrap(vars(idx)) +
#   coord_fixed()


# # Interpolate kernel densities ---------------------------------------------
# # Holds interpolated kernels at all locations
# C_s.list = list()
# 
# # Interpolate density kernels to all other locations using equation on pg 179
# for (s.idx in 1:nrow(spatial.df)) {
#   # Get spatial process location
#   s = spatial.df[s.idx, c('x', 'y')]
#   
#   # Keep track of sum of w (for normalizing)
#   w.sum = 0
#   
#   # Interpolated kernel at location s
#   C_s = matrix(0, nrow=nrow(grid.mat), ncol=ncol(grid.mat))
#   for (basis.loc.idx in 1:nrow(basis.locs)) {
#     # Get basis location
#     basis.loc = basis.locs[basis.loc.idx, c('x', 'y')]
#     
#     # Calculate basis location weight
#     w = w_k(s, basis.loc)
#     
#     # Update sum of w
#     w.sum = w.sum + w
#     
#     # Add weighted value of convolution from basis location
#     C_s = C_s + w * conv.kernel.densities[[basis.loc.idx]]
#     
#     # Null out to avoid floating variables
#     basis.loc = NULL
#   }
#   
#   # Store interpolated kernel density
#   C_s.list[[s.idx]] = C_s / w.sum
# }

# Plot densities at intermediate locations (sanity check)
# First check that "interpolated" densities at basis locations are the same as 
# densities defined at each basis location

# Merge with spatial data (keeping track of indices)
spatial.df$spatial_idx = 1:nrow(spatial.df)
basis.locs.idx = merge(spatial.df, basis.locs, by=c('x', 'y'))

# Merge with basis locations (keeping track of indices)
basis.locs$basis_idx = 1:nrow(basis.locs)
merge(basis.locs.idx, basis.locs, by=c('x', 'y'))

# Sanity check that number of basis locations didn't change
stopifnot(nrow(basis.locs.idx) == n_basis)

# Compare "interpolated" and defined densities at basis locations
for (i in 1:n_basis) {
  spatial.idx = basis.locs.idx[i, 'spatial_idx']
  C_s.interp = C_s.list[[spatial.idx]]
  
  basis.idx = basis.locs.idx[i, 'basis_idx']
  C_s.defined = conv.kernel.densities[[basis.idx]]
  
  C_s.interp.df = melt(C_s.interp)
  names(C_s.interp.df) = c('x', 'y', 'value')
  C_s.interp.df$type = 'interp'
  
  C_s.defined.df = melt(C_s.defined)
  names(C_s.defined.df) = c('x', 'y', 'value')
  C_s.defined.df$type = 'defined'
  
  C_s.df = rbind(C_s.interp.df, C_s.defined.df)
  
  g = ggplot(C_s.df) +
    geom_tile(aes(x=x, y=y, fill=value)) +
    facet_wrap(vars(type)) +
    coord_fixed() +
    ggtitle(paste('Densities for basis location', basis.idx))
  
  show(g)
}

# Plot basis densities and several interpolated densities

# Start by combining all basis densities into one matrix
densities.basis.mat = matrix(0, nrow=grid.num_x, ncol=grid.num_y)

for (i in 1:n_basis) {
  densities.basis.mat = densities.basis.mat + conv.kernel.densities[[i]]
}

densities.basis.df = melt(densities.basis.mat)
names(densities.basis.df) = c('x', 'y', 'value')
densities.basis.df$type = 'basis'

# Then hand select a few locations to show interpolated densities at
locs.to_show = data.frame(x=c(5, 10, 10, 15),
                          y=c(15, 5, 15, 5))

locs.to_show.idx = merge(spatial.df, locs.to_show, by=c('x', 'y')) %>%
  dplyr::select(spatial_idx) %>%
  as.vector %>%
  unname %>%
  unlist

densities.interp.mat = matrix(0, nrow=grid.num_x, ncol=grid.num_y)
for (i in 1:nrow(locs.to_show)) {
  locs.to_show.idx.i = locs.to_show.idx[i]
  densities.interp.mat = densities.interp.mat + conv.kernel.densities[[locs.to_show.idx.i]]
  
  g = melt(conv.kernel.densities[[locs.to_show.idx.i]]) %>%
    ggplot() +
      geom_tile(aes(x=Var1, y=Var2, fill=value)) +
      geom_point(data=locs.to_show[i, ], aes(x=x, y=y), size=3, color='red') +
      coord_fixed() +
      ggtitle(locs.to_show.idx)
  
  show(g)
}

densities.interp.df = melt(densities.interp.mat)
names(densities.interp.df) = c('x', 'y', 'value')
densities.interp.df$type = 'interp'

densities.df = rbind(densities.basis.df, densities.interp.df)

ggplot(densities.interp.df) +
  geom_tile(aes(x=x, y=y, fill=value)) +
  geom_point(data=locs.to_show, aes(x=x, y=y), size=3, color='red') +
  facet_wrap(vars(type)) +
  coord_fixed()

# Calculate covariance matrix K (pg. 180) -------------------------------------
K = matrix(NA, nrow=nrow(spatial.df), ncol=nrow(grid.df))
for (i in 1:nrow(spatial.df)) {
  for (j in 1:nrow(grid.df)) {
    K[i, j] = C_s.list[[i]][j]
  }
}

K.df = melt(K)
names(K.df) = c('x', 'y', 'value')

ggplot(K.df) +
  geom_tile(aes(x=x, y=y, fill=value)) +
  scale_fill_gradient2(low='blue', high='red') +
  coord_fixed()

# Sample from model priors ---------------------------------------------
model = cmdstan_model('stan/convolution_process.stan')

model.data = list(
  N = nrow(spatial.df),
  N_grid = nrow(grid.df),
  grid_locs = grid.df[, c('x', 'y')],
  basis_locs = basis.locs[, c('x', 'y')],
  K = K,
  y = spatial.df$value %>% as.vector
)

model.out = model$sample(
  data=model.data,
  chains=4,
  parallel_chains=4,
  iter_warmup=2000,
  iter_sampling=1000,
  fixed_param=T
)
