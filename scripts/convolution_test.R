library(gsignal)
library(cmdstanr)
library(mvtnorm)
library(plgp) # distance function

# Testing convolution in R
a = matrix(1:16, 4, 4)
b = matrix(1:4, 2, 2) / 10

# Testing convolution in Stan
model = cmdstan_model('stan/convolution_test.stan')
model.out = model$sample(
  data=list(
    N = nrow(a),
    N_kernel = nrow(b),
    a = a,
    kernel = b
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
# y.mat = matrix(y, nrow=length(lattice.x), ncol=length(lattice.y), byrow=T)
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
ggplot(lattice.df) +
  geom_tile(aes(x=x, y=y, fill=value)) +
  scale_fill_gradient2() + 
  coord_fixed() +
  ggtitle('White noise')

ggplot(whitenoise.df) +
  geom_tile(aes(x=x, y=y, fill=conv_value)) +
  scale_fill_gradient2() + 
  coord_fixed() +
  ggtitle('GP (kernel convolved with white noise)')