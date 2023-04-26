functions {
  real generalized_inverse_gaussian_lpdf(real x, int p,
                                        real a, real b) {
    return p * 0.5 * log(a / b)
      - log(2 * modified_bessel_second_kind(p, sqrt(a * b)))
      + (p - 1) * log(x)
      - (a * x + b / x) * 0.5;
 }
 
 real partial_sum_lm_lpdf(array[] real y_spatial_slice,
                  int start, int end,
                  real mu,
                  vector f,
                  real lambda_y) {
    return normal_lupdf(y_spatial_slice | mu + f[start:end], lambda_y);
  }
  
  real partial_sum_spatial_transform()
}

data {
  int<lower=0> N_knots;
  array[N_knots] vector[2] knot_locs;
  
  int<lower=0> fixed_rank_dim; // Dimension to project covariance matrix into
  
  int<lower=0> N_spatial; // Number of sites at which spatial process is measured
  array[N_spatial] vector[2] spatial_locs; // x-y coordinates of spatial process
  array[N_spatial] real y_spatial; // Measured value of spatial process at each site
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
  
  // fixed_rank_dim-dimension projection matrix
  matrix[N_spatial, fixed_rank_dim] Omega = to_matrix(normal_rng(rep_vector(0, N_spatial*fixed_rank_dim), 1./fixed_rank_dim), N_spatial, fixed_rank_dim);

  // Left spectral factor of projection matrix
  matrix[fixed_rank_dim, N_spatial] Omega_Phi = svd_U(Omega)';
}

parameters {
  // Intercept
  real mu;
  
  // Kernel (ellipse) process
  // real<lower=0> nugget_psi;
  real<lower=0> sigma_psi; // Variance of foci process
  real<lower=0> tau_psi; // Range of ellipse dependence
  real<lower=0> ell_psi; // Range of ellipse dependence
  real<lower=0> nugget_psi;
  vector<lower=0>[N_knots] psi_x; // x component of foci
  vector<lower=0>[N_knots] psi_y; // y component of foci
  
  // Kernel interpolation process
  real<lower=0> sigma_interp;
  real<lower=0> ell_interp;
  
  // (Latent) spatial process
  // real<lower=0> nugget_z; // Precision of spatial error term
  real<lower=0> sigma_z; // Variance of spatial process
  real<lower=0> ell_z; // Range of spatial process
  
  // Overall spatial process
  vector[N_spatial] eta;
  real<lower=0> lambda_y; // Precision of model error term
}

model {
  // Priors
  mu ~ std_normal();
  
  // Reference for priors: https://mc-stan.org/docs/stan-users-guide/fit-gp.html
  // nugget_psi ~ std_normal();
  nugget_psi ~ lognormal(0, 1);
  // target += generalized_inverse_gaussian_lpdf(nugget_psi | -2, 1, 0.1);
  // nugget_psi ~ inv_gamma(3, 0.1);
  sigma_psi ~ std_normal();
  // ell_psi ~ inv_gamma(3, 0.1);
  // target += generalized_inverse_gaussian_lpdf(ell_psi | -2, 1, 0.1);
  ell_psi ~ normal(0, max_dist/3);
  tau_psi ~ std_normal();
  
  sigma_interp ~ std_normal();
  ell_interp ~ inv_gamma(3, 0.1);
  // target += generalized_inverse_gaussian_lpdf(ell_interp | -2, 1, 0.1);
  // ell_interp ~ normal(0, max_dist/3);
  
  // nugget_z ~ normal(0, 0.1);
  // nugget_z ~ inv_gamma(3, 0.1);
  sigma_z ~ std_normal();
  ell_z ~ inv_gamma(3, 0.1);
  // target += generalized_inverse_gaussian_lpdf(ell_z | -2, 1, 0.1);
  // ell_z ~ normal(0, max_dist/3);
  
  eta ~ std_normal();
  lambda_y ~ std_normal();

  // Covariance matrix for foci
  // Used in multi_normal priors for foci
  matrix[N_knots, N_knots] R_psi = gp_exp_quad_cov(knot_locs, sigma_psi, ell_psi);
  for (i in 1:N_knots) {
    R_psi[i, i] = R_psi[i, i] + nugget_psi;
  }
  
  matrix[N_knots, N_knots] R_psi_chol = cholesky_decompose(R_psi);

  psi_x ~ multi_normal_cholesky(rep_vector(0., N_knots), R_psi_chol);
  psi_y ~ multi_normal_cholesky(rep_vector(0., N_knots), R_psi_chol);
  
  matrix[N_spatial, N_knots] W_interp = gp_exp_quad_cov(spatial_locs, knot_locs, sigma_interp, ell_interp);

  vector[N_spatial] psi_x_all = W_interp * psi_x;
  vector[N_spatial] psi_y_all = W_interp * psi_y;

  // Compute rotations (sampling both psi_x and psi_y near zero is problematic)
  vector[N_spatial] alpha = atan(psi_y_all ./ psi_x_all);
  for (i in 1:N_spatial) {
    if (is_nan(alpha[i])) {
      alpha[i] = 0;
    }
  }
  
  // Construct kernel covariance matrices at each site
  array[N_spatial] vector[2] spatial_locs_transformed;
  for (i in 1:N_spatial) {
    // // Component of covariance matrix that handles ellipse scaling
    // matrix[2, 2] ellipse_scale = [[sqrt(psi_x_all[i]), 0], [0, sqrt(psi_y_all[i])]];
    // // Component of covariance matrix that handles ellipse rotation
    // matrix[2, 2] rotation = [[cos(alpha[i]), sin(alpha[i])], [-sin(alpha[i]), cos(alpha[i])]];
    // Full kernel covariance matrix (square root)
    // matrix[2, 2] Sigma_sqrt = tau_psi * rotation * ellipse_scale;
    // matrix[2, 2] Sigma_sqrt = tau_psi * [[sqrt(psi_x_all[i])*cos(alpha[i]), sqrt(psi_y_all[i])*sin(alpha[i])], [-sqrt(psi_x_all[i])*sin(alpha[i]), sqrt(psi_y_all[i])*cos(alpha[i])]];
    // Full kernel covariance matrix (matrix times its transpose)
    // Equivalent to Sigma_sqrt * Sigma_sqrt'
    // Sigma_array[i] = tcrossprod(Sigma_sqrt);
    // Transform each location according to its elliptical covariance matrix
    spatial_locs_transformed[i] = tcrossprod(tau_psi * [[sqrt(psi_x_all[i])*cos(alpha[i]), sqrt(psi_y_all[i])*sin(alpha[i])], [-sqrt(psi_x_all[i])*sin(alpha[i]), sqrt(psi_y_all[i])*cos(alpha[i])]]) * spatial_locs[i];
    
    // if (is_nan(spatial_locs_transformed[i][1]) || is_nan(spatial_locs_transformed[i][2])) {
    //   print(psi_x_all[i]);
    //   print(psi_y_all[i]);
    //   print(alpha[i]);
    //   print(tau_psi);
    //   print("=========");
    // }
  }

  // Latent variable formulation
  // Reference: https://mc-stan.org/docs/stan-users-guide/fit-gp.html
  vector[N_spatial] f;
  {
    // Compute covariance matrix between transformed coordinates (simplified version)
    matrix[N_spatial, N_spatial] R_z = gp_exp_quad_cov(spatial_locs_transformed, sigma_z, ell_z);
    
    // === Banerjee linear projection method ===
    // Project down to lower dimension
    matrix[N_spatial, N_spatial] C_lp = (Omega_Phi * R_z)' / quad_form(R_z, Omega_Phi') * (Omega_Phi * R_z);
    // matrix[N_spatial, N_spatial] C_lp = quad_form(Omega_Phi * R_z, solve(Omega_Phi * R_z * Omega_Phi'));
    // matrix[N_spatial, fixed_rank_dim] C_lp = (Omega_Phi * R_z)' * solve(Omega_Phi * R_z * Omega_Phi') * Omega_Phi * R_z;
    
    // Add a nugget term (correct for variance)
    C_lp = add_diag(C_lp, diagonal(R_z) - diagonal(C_lp));
    for (n in 1:N_spatial) {
      C_lp[n, n] = C_lp[n, n] + 0.01;
    }
    
    // Cholesky decomposition for later efficiency
    matrix[N_spatial, N_spatial] C_lp_chol = cholesky_decompose(C_lp);
    // === End linear projection method ===

    // Latent variable formulation uses Cholesky decomposition (use this _or_ Banerjee linear projection)
    // matrix[N_spatial, N_spatial] R_z_chol = cholesky_decompose(R_z);
    f = C_lp_chol * eta;
  }

  // // === Sampling from linear model ===
  // y_spatial ~ normal(mu + f, lambda_y);
  int grainsize = 1;
  target += reduce_sum(partial_sum_lm_lpdf, y_spatial, 
                        grainsize,
                        mu, f, lambda_y);
}

generated quantities {
  matrix[N_spatial, N_knots] W_interp = gp_exp_quad_cov(spatial_locs, knot_locs, sigma_interp, ell_interp);

  vector[N_spatial] psi_x_all = W_interp * psi_x;
  vector[N_spatial] psi_y_all = W_interp * psi_y;
}
