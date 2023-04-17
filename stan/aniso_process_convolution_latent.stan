data {
  int<lower=0> N_spatial; // Number of sites at which spatial process is measured
  array[N_spatial] vector[2] spatial_locs; // x-y coordinates of spatial process
  vector[N_spatial] y_spatial; // Measured value of spatial process at each site
}

transformed data {
  real ell = 0; // Lower bound of range of spatial process
  real upp = 1; // Upper bound of range of spatial process
  real a = 2; // Shape parameter of diffuse gamma distribution
  real b = 0.1; // Shape parameter of diffuse gamma distribution
  real A = 0.1; // (Fixed) area of each ellipse
}

parameters {
  real mu; // Intercept
  // real<lower=0> lambda_psi; // Precision of ellipse process
  real<lower=0> lambda_y; // Precision of error term
  real<lower=0> sigma_z; // Variance of spatial process
  real<lower=0> nugget_z;
  real<lower=0> tau_z; // Range of spatial process
  vector<lower=0>[N_spatial] psi_x; // x component of foci
  vector<lower=0>[N_spatial] psi_y; // y component of foci
  real<lower=0> tau_psi; // Range of ellipse dependence
  vector[N_spatial] eta;
}

transformed parameters {
}

model {
  // Priors
  lambda_y ~ cauchy(0, 2.5);
  mu ~ std_normal();
  // lambda_psi ~ gamma(a, b);
  sigma_z ~ cauchy(0, 2.5);
  tau_z ~ cauchy(0, 2.5);
  tau_psi ~ cauchy(0, 2.5);
  nugget_z ~ cauchy(0, 2.5);
  
  // Covariance matrix for foci
  // Used in multi_normal priors for foci
  matrix[N_spatial, N_spatial] R_psi = gp_exp_quad_cov(spatial_locs, sigma_z, tau_psi);
  for (i in 1:N_spatial) {
    R_psi[i, i] = R_psi[i, i] + nugget_z;
  }
  
  psi_x ~ multi_normal(rep_vector(0., N_spatial), R_psi);
  psi_y ~ multi_normal(rep_vector(0., N_spatial), R_psi);
  
  // Hold covariance matrices at each site
  array[N_spatial] matrix[2, 2] Sigma_array;
  
  // Construct kernel covariance matrices at each site (using foci)
  array[N_spatial] vector[2] spatial_locs_transformed;
  for (i in 1:N_spatial) {
    // Compute norm of foci vector
    real psi_norm = norm2([psi_x[i], psi_y[i]]);
    // Term 1 on diagonal of kernel covariance matrix
    real term1 = sqrt(4*square(A) + pow(psi_norm, 4) * square(pi())) / 2*pi();
    // Term 2 on diagonal of kernel covariance matrix
    real term2 = square(psi_norm) / 2;
    // Rotation (give in terms of foci)
    real alpha = atan(psi_y[i] / psi_x[i]);
    // Component of covariance matrix that handles ellipse scaling
    matrix[2, 2] ellipse_scale = [[sqrt(term1 + term2), 0], [0, sqrt(term1 - term2)]];
    // Component of covariance matrix that handles ellipse rotation
    matrix[2, 2] rotation = [[cos(alpha), sin(alpha)], [-sin(alpha), cos(alpha)]];
    // Full kernel covariance matrix (square root)
    matrix[2, 2] Sigma_sqrt = tau_z * rotation * ellipse_scale ;
    // Full kernel covariance matrix (matrix times its transpose)
    // Equivalent to Sigma_sqrt * Sigma_sqrt'
    Sigma_array[i] = tcrossprod(Sigma_sqrt);
    // if (is_inf(Sigma_array[i][1, 1]) || is_inf(Sigma_array[i][2, 2])) {
    //   print(Sigma_array[i]);
    //   print(psi_x[i]);
    //   print(psi_y[i]);
    //   print("======");
    // }
    // Transform each location according to its elliptical covariance matrix
    spatial_locs_transformed[i] = Sigma_array[i] * spatial_locs[i];
  }
  
  // Latent variable formulation
  // Reference: https://mc-stan.org/docs/stan-users-guide/fit-gp.html
  vector[N_spatial] f;
  {
    // Compute covariance matrix between transformed coordinates (simplified version)
    matrix[N_spatial, N_spatial] R_z = gp_exp_quad_cov(spatial_locs_transformed, 1, tau_z);

    // diagonal elements
    for (n in 1:N_spatial) {
      R_z[n, n] = R_z[n, n] + 0.1;
    }

    matrix[N_spatial, N_spatial] R_z_chol = cholesky_decompose(R_z);
    f = R_z_chol * eta;
  }
  
  eta ~ std_normal();
  
  // === Sampling from linear model ===
  y_spatial ~ normal(mu + f, lambda_y);
}
