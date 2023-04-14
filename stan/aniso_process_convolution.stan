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
  // real<lower=0> lambda_z; // Precision of spatial process
  real<lower=0> tau_z; // Range of spatial process
  vector<lower=0>[N_spatial] psi_x; // x component of foci
  vector<lower=0>[N_spatial] psi_y; // y component of foci
  // array[N_spatial] unit_vector[2] xy; // vector representing rotation of ellipse
  real<lower=0> tau_psi; // Range of ellipse dependence
  vector[N_spatial] z; // Spatial process
  vector[N_spatial] eta;
}

transformed parameters {
}

model {
  // Priors
  lambda_y ~ gamma(a, b);
  mu ~ std_normal();
  // lambda_psi ~ gamma(a, b);
  // lambda_z ~ gamma(a, b);
  tau_z ~ gamma(a, b);
  tau_psi ~ gamma(a, b);
  
  // Covariance matrix for foci
  // Used in multi_normal priors for foci
  matrix[N_spatial, N_spatial] R_psi = gp_exp_quad_cov(spatial_locs, 1, tau_psi / 2);
  for (i in 1:N_spatial) {
    R_psi[i, i] = R_psi[i, i] + 0.1;
  }
  
  psi_x ~ multi_normal(rep_vector(0., N_spatial), R_psi);
  psi_y ~ multi_normal(rep_vector(0., N_spatial), R_psi);
  
  // // Rotation of ellipse
  // array[N_spatial] real<lower=-pi(), upper=pi()> alpha;
  // for (i in 1:N_spatial) {
  //   alpha[i] = atan2(xy[i][2], xy[i][1]);
  // }
  
  // Hold covariance matrices at each site
  array[N_spatial] matrix[2, 2] Sigma_array;
  
  // // Construct kernel covariance matrices at each site (using major/minor axes)
  // for (i in 1:N_spatial) {
  //   real s1 = square(cos(alpha[i]))/(2*square(psi_x[i])) + square(sin(alpha[i]))/(2*square(psi_y[i]));
  //   real s12 = sin(2*alpha[i])/(4*square(psi_y[i])) - sin(2*alpha[i])/(4*square(psi_x[i]));
  //   real s2 = square(sin(alpha[i]))/(2*square(psi_x[i])) + square(cos(alpha[i]))/(2*square(psi_y[i]));
  //   Sigma_array[i] = [[s1, s12], [s12, s2]];
  // }
  
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
    // matrix[2, 2] Sigma_sqrt = tau_z * ellipse_scale;
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
  
  // Compute covariance matrix between transformed coordinates (simplified version)
  matrix[N_spatial, N_spatial] R_z = gp_exp_quad_cov(spatial_locs_transformed, 1, 1);
  
  // // === Compute spatial process correlation matrix (pg. 3, eq. 2)
  // matrix[N_spatial, N_spatial] R_z;
  // for (i in 1:N_spatial) {
  //   // Access components used for later computations
  //   matrix[2, 2] Sigma_i = Sigma_array[i];
  //   real a_i = sqrt(Sigma_i[1, 1]);
  //   real b_i = sqrt(Sigma_i[2, 2]);
  //   real p_i = Sigma_i[1, 2] / (a_i * b_i);
  //   vector[2] s_i = spatial_locs[i];
  // 
  //   // Sub-terms used regularly
  //   real a_ib_i = a_i * b_i;
  //   real a_i2 = Sigma_i[1, 1];
  //   real b_i2 = Sigma_i[2, 2];
  //   real p_i2 = square(p_i);
  // 
  //   for (j in (i+1):N_spatial) {
  //     // Access components used for later computations
  //     matrix[2, 2] Sigma_j = Sigma_array[j];
  //     real a_j = sqrt(Sigma_j[1, 1]);
  //     real b_j = sqrt(Sigma_j[2, 2]);
  //     real p_j = Sigma_j[1, 2] / (a_j * b_j);
  //     vector[2] s_j = spatial_locs[j];
  // 
  //     // Sub-terms used regularly
  //     real a_jb_j = a_j * b_j;
  //     real a_j2 = Sigma_j[1, 1];
  //     real b_j2 = Sigma_j[2, 2];
  //     real p_j2 = square(p_j);
  // 
  //     // Compute W
  //     matrix[2, 2] W = [[b_i2 + b_j2, -(p_i*a_ib_i + p_j*a_jb_j)], [-(p_i*a_ib_i + p_j*a_jb_j), a_i2 + a_j2]];
  //     // if (is_inf(W[1, 1]) || is_inf(W[1, 2]) || is_inf(W[2, 2])) {
  //     //   print(Sigma_i);
  //     //   print(Sigma_j);
  //     //   print("==========");
  //     // }
  // 
  //     // Compute q1
  //     real q1_term1 = 2*pi()*a_ib_i*a_jb_j * sqrt((1-p_i2)*(1-p_j2));
  //     real q1_term2 = sqrt(-((p_i2-1)*b_i2 + (p_j2-1)*b_j2) / ((p_i2-1)*(p_j2-1)*b_i2*b_j2));
  //     // Higdon may have a typo here (q1_term3_num). Second and third terms seem like they should be symmetric, but are not (different parantheses).
  //     // I wrote the symmetric version, but it's probably worth working it out at some point to verify.
  //     real q1_term3_num = 2*a_ib_i*a_jb_j*p_i*p_j + a_i2*((p_i2-1)*b_i2-b_j2) + a_j2*((p_j2-1)*b_j2-b_i2);
  //     real q1_term3_den = a_i2*a_j2*((p_i2 - 1)*b_i2 + (p_j2 - 1)*b_j2);
  //     real q1_term3 = sqrt(q1_term3_num / q1_term3_den);
  //     real q1 = q1_term1 * q1_term2 * q1_term3;
  //     // if (abs(q1) < 1e-10) {
  //     //   print("q1 = 0");
  //     //   print("Sigma_i = ", Sigma_i);
  //     //   print("Sigma_j = ", Sigma_j);
  //     //   print("=====");
  //     // }
  // 
  //     // Compute q2
  //     real q2 = 2*(2*p_i*p_j*a_ib_i*a_jb_j + a_i2*((p_i2-1)*b_i2-b_j2) - a_j2*((p_j2-1)*b_j2-b_i2));
  //     
  //     // === Override values for testing ===
  //     q1 = 1;
  //     q2 = 1;
  //     
  //     real rho;
  //     if (abs(q2) < 1e-10) {
  //       rho = 0;
  //     } else {
  //       // Compute correlation function rho(s, s')
  //       // real rho = (1 / q1) * exp(-(1 / q2) * (spatial_locs[i] - spatial_locs[j])' * W * (spatial_locs[i] - spatial_locs[j]));
  //       // rho = exp(-quad_form_sym(W, spatial_locs[i] - spatial_locs[j]) / q2) / q1;
  //     }
  // 
  //     // if (is_inf(rho) || is_nan(rho)) {
  //     //   print("Sigma_i = ", Sigma_i);
  //     //   print("Sigma_j = ", Sigma_j);
  //     //   print("s_i = ", s_i);
  //     //   print("s_j = ", s_j);
  //     //   print("=====");
  //     // }
  // 
  //     // Store correlation value in covariance matrix R_z
  //     R_z[i, j] = rho;
  //     R_z[j, i] = rho; // Symmetric matrix
  //   }
  //   R_z[i, i] = 1; // Variance of spatial process
  // }
  
  // Convert correlation matrix R_z to a covariance matrix
  // https://stats.stackexchange.com/questions/62850/obtaining-covariance-matrix-from-correlation-matrix
  // vector[N_spatial] z_var = rep_vector(lambda_z, N_spatial);
  // matrix[N_spatial, N_spatial] R_z_cov = quad_form_diag(R_z, z_var);
  // for (i in 1:N_spatial) {
  //   R_z_cov[i, i] = R_z_cov[i, i] + 0.1;
  // }
  // if (min(eigenvalues_sym(R_z_cov)) < 0) {
  //   print(min(eigenvalues_sym(R_z_cov)));
  // }
  // print("=====");
  // matrix[N_spatial, N_spatial] R_z_chol = cholesky_decompose(R_z);
  
  // === Sample from spatial process ===
  matrix[N_spatial, N_spatial] R_z_chol = cholesky_decompose(R_z);
  z ~ multi_normal_cholesky(rep_vector(0, N_spatial), R_z_chol);
  
  // === Sampling from linear model ===
  y_spatial ~ normal(mu + z, lambda_y);
  
  // // === Spatial process k at each distance s-omega_j (see Higdon "Space and 
  // space-time modeling using process convolutions" pg. 6 eq. 4) ===
  // matrix[N_spatial, N_latent] spatial_at_dist;
  // for (i in 1:N_spatial) {
  //   for (j in 1:N_latent) {
  //     spatial_at_dist[i, j] = multi_normal_lpdf(spatial_locs[i] - latent_locs[j] | rep_vector(0, 2), R_z_array[j]);  
  //   }
  // }
  // 
  // // === Convolve latent process x(omega_j) and spatial process k(s-omega_j) ===
  // vector[N_spatial] z;
  // z = spatial_at_dist * x_j;
}
