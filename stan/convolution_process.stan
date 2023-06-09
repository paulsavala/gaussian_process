// functions {
//   matrix conv2(matrix a, matrix kernel) {
//     int N = rows(a);
//     int N_kernel = rows(kernel);
//     
//     matrix[N, N] conv_out;
//     matrix[N, N] kernel_pad;
//   
//     // Pad kernel if necessary
//     if (N_kernel < N) {
//       kernel_pad = rep_matrix(0, N, N);
//       kernel_pad[1:N_kernel, 1:N_kernel] = kernel;
//     } else {
//       kernel_pad = kernel;
//     }
//     conv_out = get_real(inv_fft2(fft2(a) .* fft2(kernel_pad)));
//     
//     return conv_out;
//   }
// }

data {
  int<lower=0> N;
  int<lower=0> N_grid;
  
  matrix[N, N_grid] K;
  
  vector[N] y;
}

parameters {
  // Grid process parameters
  real mu_x;
  real<lower=0> sigma_x2_inv;
  vector[N_grid] grid_process; 
  
  // Error variance
  real<lower=0> sigma_eps2_inv;
}

transformed parameters {
  real sigma_x2 = 1 / sigma_x2_inv;
  real sigma_eps2 = 1 / sigma_eps2_inv;
}

model {
  // Grid process
  mu_x ~ std_normal();
  sigma_x2_inv ~ gamma(1, 0.005);
  
  // Error variance
  sigma_eps2_inv ~ gamma(1, 0.005);
  
  grid_process ~ normal(mu_x, sigma_x2);
  
  target += normal_lpdf(y | K*grid_process, sigma_eps2);
}

generated quantities {
  // Return predictions
  array[N] real y_sim;
  y_sim = normal_rng(K*grid_process, sigma_eps2);
  
  // // Return log likelihood (for LOOIC computations)
  // vector[N] log_lik;
  // for (i in 1:N) {
  //   log_lik[i] = normal_lpdf(y[i] | (K*grid_process)[i], sigma_eps2)
  // }
  
}
