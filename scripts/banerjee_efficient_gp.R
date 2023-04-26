# https://academic.oup.com/biomet/article-pdf/100/1/75/481197/ass068.pdf

# ======= Data prep =========
# Test data
x = seq(from=-3, to=3, length.out=100)

xx = expand.grid(x, x)

# Full covariance matrix
C = matrix(exp(-abs(xx$x-xx$y)**2), nrow=length(x), ncol=length(x))

# "Nice to have" values
n = nrow(C)
m = as.integer(sqrt(n))

# Condition number of original covariance matrix
condition.C = eigen(C, symmetric=T, only.values=T)$values[1] / eigen(C, symmetric=T, only.values=T)$values[n]
condition.C

# ======= Algorithm =========
algo_1 = function(C, m) {
  # JL matrix (low-rank projection matrix)
  Omega = matrix(rnorm(n*m, mean=0, sd=1/m), nrow=n, ncol=m)
  
  # SVD of matrix (only need the left factor, svd_U in Stan)
  Phi = svd(Omega)$u %>% t # Phi is the transpose of the left factor
  
  # # Low-rank (m x m) embedding of C
  # C.1 = Phi %*% C %*% t(Phi)
  # 
  # # Cholesky decomposition of low-rank embedding
  # B = chol(C.1) %>% t # R returns Cholesky factor as transposed from paper
  # 
  # # Nystrom factor (used later)
  # C.2 = C %*% t(Phi) %*% solve(t(B))
  # 
  # # Left factor and singular values of Nystrom factor
  # U.2 = svd(C.2)$u
  # D.2 = svd(C.2)$d %>% diag
  # 
  # # Approximate spectral decomposition
  # C.fr = U.2 %*% (D.2**2) %*% t(U.2)
  
  # Low-rank embedding of covariance matrix C (not used in algorithm)
  C.lp = t(Phi %*% C) %*% solve(Phi %*% C %*% t(Phi)) %*% Phi %*% C
  
  return(C.lp)
}

# ======= Evaluation =======
norms = c()

for (m in seq(from=2, to=25)) {
  norms = c(norms, norm(C - algo_1(C, m), type="2"))
}
