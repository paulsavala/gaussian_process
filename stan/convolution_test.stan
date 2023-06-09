functions {
  matrix conv2(matrix a, matrix kernel) {
    int N = rows(a);
    int N_kernel = rows(kernel);
    
    matrix[N, N] conv_out;
    matrix[N, N] kernel_pad;
  
    // Pad kernel if necessary
    if (N_kernel < N) {
      kernel_pad = rep_matrix(0, N, N);
      kernel_pad[1:N_kernel, 1:N_kernel] = kernel;
    } else {
      kernel_pad = kernel;
    }
    
    // Convolution
    conv_out = get_real(inv_fft2(fft2(a) .* fft2(kernel_pad)));
    
    return conv_out;
  }
}

data {
  int<lower=0> N;
  int<lower=0> N_kernel;
  matrix[N, N] a;
  matrix[N_kernel, N_kernel] kernel;
}

generated quantities {
  matrix[N, N] out = conv2(a, kernel);
}

