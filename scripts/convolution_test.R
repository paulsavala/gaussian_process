library(gsignal) # conv2 - 2d convolution
library(cmdstanr)
library(mvtnorm)
library(plgp) # distance function
library(tidyr)
library(dplyr)
library(reshape2)
library(ggplot2)

# ==== Testing convolution in R using toy data ====
rot.mat = function(a, b, theta) {
  stretch = matrix(c(a, 0, 0, b), nrow=2, ncol=2, byrow=T)
  rotate = matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow=2, ncol=2, byrow=T)
  out = rotate %*% stretch
  return(out)
}

density.mat = function(a, sigma) {
  out = dmvnorm(a, 0, sigma, log=FALSE)
  return(out)
}

mat.to_long = function(m, type) {
  m.df = melt(m)
  names(m.df) = c('x', 'y', 'value')
  m.df$type = type
  return(m.df)
}

conv_fft <- function(matrix1, matrix2) {
  # Get the dimensions of the input matrices
  nrow1 <- nrow(matrix1)
  ncol1 <- ncol(matrix1)
  nrow2 <- nrow(matrix2)
  ncol2 <- ncol(matrix2)
  
  # Compute the size of the output matrix
  output_rows <- nrow1 + nrow2 - 1
  output_cols <- ncol1 + ncol2 - 1
  
  # Pad the matrices to the size of the output matrix
  padded_matrix1 <- cbind(matrix1, matrix(0, nrow = nrow1, ncol = output_cols - ncol1))
  padded_matrix1 <- rbind(padded_matrix1, matrix(0, nrow = output_rows - nrow1, ncol = output_cols))
  padded_matrix2 <- cbind(matrix2, matrix(0, nrow = nrow2, ncol = output_cols - ncol2))
  padded_matrix2 <- rbind(padded_matrix2, matrix(0, nrow = output_rows - nrow2, ncol = output_cols))
  
  # Perform FFT on the padded matrices
  fft1 <- fft(padded_matrix1)
  fft2 <- fft(padded_matrix2)
  
  # Compute element-wise multiplication in the frequency domain
  fft_multiply <- fft1 * fft2
  
  # Perform inverse FFT to obtain the convolution result
  convolution_result <- Re(fft(fft_multiply, inverse = TRUE) / (output_rows * output_cols))
  
  # Return the convolution result
  return(convolution_result)
}


# Define the mean and covariance matrix
mean_vec <- c(5, 5)

# Generate a grid of spatial locations
grid_x <- seq(0, 10, length.out = 10)
grid_y <- seq(0, 10, length.out = 10)
grid_points <- expand.grid(x = grid_x, y = grid_y)

k = rot.mat(1, 2, pi/4)
k.cov = solve(k %*% t(k))

# Calculate the density for each spatial location
density <- dmvnorm(grid_points, mean = mean_vec, sigma = k.cov)

# Reshape the density vector into a matrix
density_matrix <- matrix(density, nrow = 10, ncol = 10, byrow = TRUE)

# For plotting later
density_matrix.df = mat.to_long(density_matrix, 'density')

ggplot(density_matrix.df) +
  geom_tile(aes(x=x, y=y, fill=value)) +
  coord_fixed()

# Space
a = matrix(rnorm(10000), 100, 100)
a.df = mat.to_long(a, 'orig')

# Convolution
# a.conv_k = conv_fft(a, density_matrix)
a.conv_k = conv2(a, density_matrix, 'same')
a.conv_k.df = mat.to_long(a.conv_k, 'conv')
a.conv_k.df$value = a.conv_k.df$value

# Stack together for graphing
a_k = rbind(a.conv_k.df, density_matrix.df)

ggplot(a_k) +
  geom_tile(aes(x=x, y=y, fill=value)) +
  facet_wrap(vars(type)) +
  scale_fill_gradient2(low='blue', high='red') +
  coord_fixed()

# Testing convolution in Stan
model = cmdstan_model('stan/convolution_test.stan')
model.out = model$sample(
  data=list(
    N = nrow(a),
    N_kernel = nrow(k),
    a = a,
    kernel = k
  ),
  fixed_param = TRUE
)

out.stan = matrix(model.out$summary()$mean, 4, 4)


# === Test on generated data ===
# = Generate background data from a GP =
# Define a fine lattice of points
lattice.x = seq(from=0, to=10, length.out=20) # 20 x-values
lattice.y = seq(from=0, to=10, length.out=20) # 20 y-values
lattice.df = expand.grid(lattice.x, lattice.y) %>% as.data.frame
names(lattice.df) = c('x', 'y')

# Compute distance between points
D = distance(lattice.df[, c('x', 'y')])

# Gaussian covariance matrix
eps = sqrt(.Machine$double.eps)
n = nrow(D)
Sigma = exp(-D) + diag(eps, n)

# Generate from GP with covariance Sigma and mean zero
y = rmvnorm(1, sigma=Sigma)
y.mat = matrix(y, nrow=length(lattice.x), ncol=length(lattice.y), byrow=T)
lattice.df$value = t(y)

# Visualize GP
ggplot(lattice.df) +
  geom_tile(aes(x=x, y=y, fill=value)) +
  scale_fill_gradient2(low='blue', high='red', midpoint=0) +
  coord_fixed()

# = Define identity convolution kernel (sanity check) =
# Identity kernel
kernel.identity = matrix(c(1, 0, 0, 1), byrow=T, nrow=2)

# = Define white noise process =
whitenoise.vals = rnorm(nrow(lattice.df), mean=0, sd=1)
whitenoise.df = expand.grid(lattice.x, lattice.y) %>% as.data.frame
names(whitenoise.df) = c('x', 'y')
whitenoise.df$value = whitenoise.vals
whitenoise.mat = matrix(whitenoise.vals, 
                        nrow=length(lattice.x), ncol=length(lattice.y))

# = Convolution =
model.toy.out = model$sample(
  data=list(
    N = nrow(y.mat),
    N_kernel = nrow(kernel.identity),
    a = y.mat %>% t,
    kernel = kernel.identity
  ),
  fixed_param = TRUE
)

whitenoise.df$conv_value = model.toy.out$summary()$mean
out.stan.mat = matrix(model.toy.out$summary()$mean, nrow(whitenoise.mat), nrow(whitenoise.mat))


# Graphing
ggplot(whitenoise.df) +
  geom_tile(aes(x=x, y=y, fill=value)) +
  scale_fill_gradient2() + 
  coord_fixed() +
  ggtitle('White noise')

ggplot(whitenoise.df) +
  geom_tile(aes(x=x, y=y, fill=conv_value)) +
  scale_fill_gradient2() + 
  coord_fixed() +
  ggtitle('GP (kernel convolved with identity matrix)')
