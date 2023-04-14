# Setup -----------------------------------------------------------------------------------------------------------
# Imports (Graphing)
library(ggplot2)  # graphing; ggplot() + ...
library(ggforce)  # graphing circles; geom_circle()

# Imports (Modeling)
library(cmdstanr)

# Imports
library(dplyr)    # piping and QoL functions; %>%
library(scales)   # normalizing lists between 0 and 1; rescale()
library(stringr)  # manipulating strings

# Knot Evaluation
source("scripts/load_epa_data.R")



# Parameters --------------------------------------------------------------
A = 1 # Area of each ellipse (fixed)


# Data loading ------------------------------------------------------------
epa.df = load_epa_data()

# California (roughly) ----------------------------------------------------
# Filter to approximately California
epa.ca.df = epa.df %>%
  dplyr::filter(x < 0.2, y > 0.25, y < 0.75) %>%
  dplyr::filter((x < 0.125) | (y < 0.4))

epa.ca.df = epa.ca.df[, c('x', 'y', 'pm2_5')]

# Sample for easy testing
epa.ca.sample.idx = sample(1:nrow(epa.ca.df), 20)
epa.ca.sample.df = epa.ca.df[epa.ca.sample.idx, ]


# Stan model fitting ------------------------------------------------------
# Plot CA and sample
ggplot(epa.ca.df) +
  geom_point(aes(x=x, y=y, fill=pm2_5), size=3, pch=21, color='black', alpha=0.1) +
  geom_point(data=epa.ca.sample.df, aes(x=x, y=y, fill=pm2_5), size=3, pch=21, color='black') +
  scale_fill_gradient2(low='blue', high='red', midpoint=0) +
  ggtitle('California (sample points highlighted)')

# Spatial process sites and data
spatial_locs = epa.ca.sample.df[, c('x', 'y')]
y_spatial = epa.ca.sample.df$pm2_5 %>% as.vector
N_spatial = nrow(spatial_locs)

# Data
data = list(N_spatial=N_spatial,
            spatial_locs=spatial_locs,
            y_spatial=y_spatial)

# Model
model = cmdstan_model('stan/aniso_process_convolution.stan')
fit = model$sample(data=data,
                   parallel_chains=4,
                   iter_warmup=1000,
                   max_treedepth=15)


# Visualizing transformation ----------------------------------------------
# === Function to plot 2d normal given covariance matrix
normal_pdf = function(sigma=matrix(c(1, 0, 0, 1), byrow=T, nrow=2), plot=TRUE) {
  x = seq(from=-3, to=3, length.out=20)
  y = seq(from=-3, to=3, length.out=20)
  x_y = expand.grid(x, y)
  
  if (plot) {
    x_y$density = dmvnorm(x_y, mean = rep(0, 2), sigma = sigma)
    ggplot(x_y) + 
      geom_tile(aes(x=Var1, y=Var2, fill=density)) +
      coord_fixed()
  }
}

# === Transformation matrix using foci ===
A = pi
psi = c(0, 0)
alpha = 0

transform_mat.foci = function(A, psi, alpha) {
  psi.norm = norm(psi, type='2')
  term1 = sqrt(4*A**2 + psi.norm**4*pi**2) / (2*pi)
  term2 = psi.norm**2 / 2
  
  scale_mat = matrix(c(sqrt(term1 + term2), 0, 0, sqrt(term1 - term2)), 
                     byrow=TRUE, nrow=2)
  rotation_mat = matrix(c(cos(alpha), sin(alpha), -sin(alpha), cos(alpha)),
                        byrow=TRUE, nrow=2)
  
  # return(scale_mat %*% rotation_mat)
  return(rotation_mat %*% scale_mat)
}

# === Correlation function rho ===
# Centers of each ellipse (covariance matrix)
s.1 = matrix(c(0, 0), nrow=2)
s.2 = matrix(c(1, 1), nrow=2)

# Covariance (ellipse) matrices
m.foci.1 = transform_mat.foci(A=pi, psi=c(1, 0), alpha=0)
m.foci.2 = transform_mat.foci(A=pi, psi=c(1, 0), alpha=pi/4)

m.foci.1.cov = m.foci.1 %*% t(m.foci.1)
m.foci.2.cov = m.foci.2 %*% t(m.foci.2)

# Plot covariance matrices (sanity check)
normal_pdf(m.foci.1.cov)
normal_pdf(m.foci.2.cov)

# Terms used to compute correlation matrix
a.1 = sqrt(m.foci.1.cov[1, 1])
b.1 = sqrt(m.foci.1.cov[2, 2])
p.1 = m.foci.1.cov[1, 2] / (a.1 * b.1)

a.2 = sqrt(m.foci.2.cov[1, 1])
b.2 = sqrt(m.foci.2.cov[2, 2])
p.2 = m.foci.2.cov[1, 2] / (a.2 * b.2)

# W
W = m.foci.1.cov + m.foci.2.cov
W[1, 2] = -W[1, 2]
W[2, 1] = -W[2, 1]

# q1
q1.term1 = 2*pi*a.1*a.2*b.1*b.2*sqrt((1-p.1**2)*(1-p.2**2))
q1.term2 = sqrt(-((p.1**2-1)*b.1**2 + (p.2**2-1)*b.2**2) / (p.1**2-1)*(p.2**2-1)*b.1**2*b.2**2)
q1.term3.num = 2*p.1*p.2*a.1*a.2*b.1*b.2 + a.1**2*((p.1**2-1)*b.1**2-b.2**2) + a.2**2*((p.2**2-1)*b.2**2-b.1**2)
q1.term3.den = a.1**2*a.2**2*((p.1**2-1)*b.1**2 + (p.2**2-1)*b.2**2)
q1.term3 = sqrt(q1.term3.num /q1.term3.den)
q1 = q1.term1 * q1.term2 * q1.term3

# q2
q2 = 2*(2*p.1*p.2*a.1*a.2*b.1*b.2 + a.1**2*((p.1**2-1)*b.1**2-b.2**2) - a.2**2*((p.2**2-1)*b.2**2-b.1**2))

# rho
rho = (1 / q1) * exp(-(1 / q2) * t(s.1-s.2) %*% W %*% (s.1 - s.2))



