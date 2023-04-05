functions {
  // Two-dimensional convolution
  matrix conv2(matrix a, matrix kernel) {
    int N = rows(a);
    int N_kernel = rows(kernel);
    
    matrix[N, N] conv_out;
    matrix[N, N] kernel_pad = rep_matrix(0, N, N);
  
    // Pad kernel if necessary
    if (N_kernel < N) {
      kernel_pad[1:N_kernel, 1:N_kernel] = kernel;
    } else {
      kernel_pad = kernel;
    }
    conv_out = get_real(inv_fft2(fft2(a) .* fft2(kernel_pad)));
    
    return conv_out;
  }
}

data {
  int<lower=0> N_spatial; // Number of sites at which spatial process is measured
  array[N_spatial] vector[2] spatial_locs; // x-y coordinates of spatial process
  vector[N_spatial] y_spatial; // Measured value of spatial process at each site
  int<lower=0> N_latent; // Number of sites at which latent process is modeled
  array[N_latent] vector[2] latent_locs; // x-y coordinates of latent process
}

transformed data {
  real ell = 3; // Lower bound of range of spatial process
  real u = 200; // Upper bound of range of spatial process
  real a = 1; // Shape parameter of diffuse gamma distribution
  real b = 10; // Shape parameter of diffuse gamma distribution
  real A = 0.5; // (Fixed) area of each ellipse
  
  // // Calculate distance between each site for the spatial process and each site for the latent process
  // matrix[N_spatial, N_latent] spatial_latent_distances = rep_matrix(0, N_spatial, N_latent);
  // for (i in 1:N_spatial) {
  //   for (j in 1:N_latent) {
  //     spatial_latent_distances[i, j] = distance(spatial_locs[i], latent_locs[j]);
  //   }
  // }
}

parameters {
  real mu; // Intercept
  real<lower=0> lambda_x; // Precision (1 / standard deviation) of latent process
  real<lower=0> lambda_y; // Precision of error term
  real<lower=0> lambda_z; // Precision of spatial process
  real<lower=ell, upper=u> tau_z; // Range of spatial process
  vector<lower=0>[N_latent] psi_x; // x component of foci
  vector<lower=0>[N_latent] psi_y; // y component of foci
  real<lower=ell> tau_psi; // Range of ellipse dependence
  vector[N_latent] x_j; // Latent process
}

transformed parameters {
  // Covariance matrix for foci
  matrix[N_latent, N_latent] R_psi = gp_exp_quad_cov(latent_locs, 1, tau_psi / 2);
}

model {
  // Priors
  lambda_x ~ inv_gamma(a, b);
  lambda_y ~ inv_gamma(a, b);
  mu ~ std_normal();
  lambda_z ~ inv_gamma(a, b);
  tau_z ~ uniform(ell, u);
  tau_psi ~ uniform(ell, u);
  psi_x ~ multi_normal(rep_vector(0., N_latent), R_psi);
  psi_y ~ multi_normal(rep_vector(0., N_latent), R_psi);
  x_j ~ normal(0, lambda_y);
  
  // Hold covariance matrices at each latent site
  array[N_latent] matrix[2, 2] R_z_array;
  
  // Construct kernel covariance matrices
  for (j in 1:N_latent) {
    // Compute norm of foci vector
    real psi_norm = norm2([psi_x[j], psi_y[j]]);
    // Term 1 on diagonal of kernel covariance matrix
    real term1 = sqrt(4*square(A) + pow(psi_norm, 4) * square(pi())) / 2*pi();
    // Term 2 on diagonal of kernel covariance matrix
    real term2 = square(psi_norm) / 2;
    // Rotation (give in terms of foci)
    real alpha = atan(psi_y[j] / psi_x[j]);
    // Component of covariance matrix that handles ellipse scaling
    matrix[2, 2] ellipse_scale = [[sqrt(term1 + term2), 0], [0, sqrt(term1 - term2)]];
    // Component of covariance matrix that handles ellipse rotation
    matrix[2, 2] rotation = [[cos(alpha), sin(alpha)], [-sin(alpha), cos(alpha)]];
    // Full kernel covariance matrix (square root)
    matrix[2, 2] R_z_sqrt = tau_z * ellipse_scale * rotation;
    // Full kernel covariance matrix
    matrix[2, 2] R_z = R_z_sqrt * R_z_sqrt' / lambda_z;
    // Store covariance matrix
    R_z_array[j] = R_z;
  }
  
  // === Spatial process k at each distance s-omega_j (see Higdon "Space and 
  // space-time modeling using process convolutions" pg. 6 eq. 4) ===
  matrix[N_spatial, N_latent] spatial_at_dist;
  for (i in 1:N_spatial) {
    for (j in 1:N_latent) {
      spatial_at_dist[i, j] = multi_normal_lpdf(spatial_locs[i] - latent_locs[j] | rep_vector(0, 2), R_z_array[j]);  
    }
  }
  
  // === Convolve latent process x(omega_j) and spatial process k(s-omega_j) ===
  vector[N_spatial] z;
  z = spatial_at_dist * x_j;
  
  // === Sampling statement ===
  y_spatial ~ normal(mu + z, 1 / lambda_y);
}

