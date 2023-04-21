data {
  int<lower=0> N_spatial; // Number of sites at which spatial process is measured
  array[N_spatial] vector[2] spatial_locs; // x-y coordinates of spatial process
  vector[N_spatial] y_spatial; // Measured value of spatial process at each site
}

transformed data {
  // Compute the maximum distance between any two points. This is used
  // for spatial dependence priors.
  real max_dist = 0;
  real temp_dist;
  for (i in 1:N_spatial) {
    for (j in (i+1):N_spatial) {
      temp_dist = distance(spatial_locs[i], spatial_locs[j]);
      max_dist = max([max_dist, temp_dist]);
    }
  }
}

parameters {
  real mu; // Intercept
  
  // Kernel (ellipse) process
  real<lower=0> nugget_psi;
  real<lower=0> sigma_psi; // Variance of foci process
  real<lower=0> tau_psi; // Range of ellipse dependence
  vector<lower=0>[N_spatial] psi_x; // x component of foci
  vector<lower=0>[N_spatial] psi_y; // y component of foci
  
  // (Latent) spatial process
  real<lower=0> nugget_z; // Precision of spatial error term
  real<lower=0> sigma_z; // Variance of spatial process
  real<lower=0> tau_z; // Range of spatial process
  
  // Overall spatial process
  vector[N_spatial] eta;
  real<lower=0> lambda_y; // Precision of model error term
}

model {
  // Priors
  mu ~ std_normal();
  
  nugget_psi ~ std_normal();
  sigma_psi ~ normal(0, max_dist);
  ell_psi ~ normal(0, max_dist);
  tau_psi ~ normal(0, max_dist);
  
  nugget_z ~ std_normal();
  sigma_z ~ normal(0, max_dist);
  ell_z ~ normal(0, max_dist);
  
  eta ~ std_normal();
  lambda_y ~ cauchy(0, 2.5);

  // Covariance matrix for foci
  // Used in multi_normal priors for foci
  matrix[N_spatial, N_spatial] R_psi = gp_exp_quad_cov(spatial_locs, sigma_psi, ell_psi);
  for (i in 1:N_spatial) {
    R_psi[i, i] = R_psi[i, i] + nugget_psi;
  }

  psi_x ~ multi_normal(rep_vector(0., N_spatial), R_psi);
  psi_y ~ multi_normal(rep_vector(0., N_spatial), R_psi);

  // Hold covariance matrices at each site
  array[N_spatial] matrix[2, 2] Sigma_array;

  // Construct kernel covariance matrices at each site (using foci)
  array[N_spatial] vector[2] spatial_locs_transformed;
  for (i in 1:N_spatial) {
    // Rotation (give in terms of foci)
    real alpha = atan(psi_y[i] / psi_x[i]);
    // Component of covariance matrix that handles ellipse scaling
    matrix[2, 2] ellipse_scale = [[sqrt(psi_x[i]), 0], [0, sqrt(psi_y[i])]];
    // Component of covariance matrix that handles ellipse rotation
    matrix[2, 2] rotation = [[cos(alpha), sin(alpha)], [-sin(alpha), cos(alpha)]];
    // Full kernel covariance matrix (square root)
    matrix[2, 2] Sigma_sqrt = tau_psi * rotation * ellipse_scale;
    // Full kernel covariance matrix (matrix times its transpose)
    // Equivalent to Sigma_sqrt * Sigma_sqrt'
    Sigma_array[i] = tcrossprod(Sigma_sqrt);
    // Transform each location according to its elliptical covariance matrix
    spatial_locs_transformed[i] = Sigma_array[i] * spatial_locs[i];
  }

  // Latent variable formulation
  // Reference: https://mc-stan.org/docs/stan-users-guide/fit-gp.html
  vector[N_spatial] f;
  {
    // Compute covariance matrix between transformed coordinates (simplified version)
    matrix[N_spatial, N_spatial] R_z = gp_exp_quad_cov(spatial_locs_transformed, sigma_z, ell_z);

    // diagonal elements
    for (n in 1:N_spatial) {
      R_z[n, n] = R_z[n, n] + nugget_z;
    }

    matrix[N_spatial, N_spatial] R_z_chol = cholesky_decompose(R_z);
    f = R_z_chol * eta;
  }

  // === Sampling from linear model ===
  y_locs ~ normal(mu + f, lambda_y);
}
