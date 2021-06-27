#@title Variance for fuzzy RD at the boundary
#
#@description This function computes various quantities needed to construct the
#  variance and adjusted variance for a fuzzy RD.
#
#@param dat A list with the following components: \describe{ \item{y}{A vector
#  treatment outcomes. In a fuzzy RD, this variable can also be a vector of
#  treatment assignment variables.} \item{x}{A vector of covariates.}
#  \item{z}{A vector of treatment assignment for different observations. 1 for
#  receiving the treatment and 0 otherwise.} \item{q}{Number of quantiles to be
#  used in estimation. Defaults to 5. It has to be an odd number.} \item{h}{A
#  scalar bandwidth.} \item{tau}{Quantile positions that correspond to the q
#  quantiles. They are obtained by \code{tau = (1:q)/(q+1)}.} \item{p}{The
#  polynomial degree.} \item{n_all}{Total number of observations in the RD
#  study.} \item{f0_hat}{Estimated density of the covariate at x = 0.}
#  \item{fd1}{Estimated first derivative of the density of the covariate at x =
#  0.} }
#@param t0 Treatment effect under the null. Defaults to 0.
#@param kernID Kernel id number. Defaults to 0. \enumerate{ \item \code{kernID
#  = 0}: triangular kernel. \item \code{kernID = 1}:   biweight kernel. \item
#  \code{kernID = 2}:  Epanechnikov kernel. \item \code{kernID = 3}:  Gaussian
#  kernel. \item \code{kernID = 4}:  tricube kernel. \item \code{kernID = 5}:
#  triweight kernel. \item \code{kernID = 6}:  uniform kernel. }
#@param left A logical variable that takes the value \code{TRUE} for data to
#  the left of (below) the cutoff. Defaults to \code{TRUE}.
#@param maxit maximum iteration number in the MM algorithm for quantile
#  estimation. Defaults to 100.
#@param tol  Convergence criterion in the MM algorithm. Defaults to 1.0e-4.
#@param para A 0/1 variable specifying whether to use parallel computing.
#  Defaults to 1.
#@param grainsize Minimum chunk size for parallelization. Defaults to 1.
#@param llr.residuals Whether to use the residuals from a local linear
#  regression as the input to compute the LCQR standard errors and the
#  corresponding bandwidths. Defaults to \code{FALSE}. If this option is set to
#  \code{TRUE}, the treatment estimate and the bias-correction is still done in
#  LCQR. We use the same kernel function used in LCQR in the local linear
#  regression to obtain the residuals and use them to compute the unadjusted
#  and adjusted asymptotic standard errors and the bandwidths. This option will
#  improve the speed but its accuracy is unclear. One can use this option to
#  get a quick estimate of the standard errors when the sample size is large.
#
#@return \code{boundary_var_fuzzy} returns a list with the following
#  components: \item{var_my}{Estimated variance for the conditional mean of the
#  treatment outcome variable.} \item{var_mt}{Estimated variance for the
#  conditional mean of the treatment assignment variable.}
#  \item{cov_mymt}{Estimated covariance between the conditional means of the
#  treatment outcome variable and the treatment assignment variable.}
#  \item{var_adj_y}{Adjusted estimate for the variance of the conditional mean
#  of the treatment outcome variable.} \item{var_adj_t}{Adjusted estimate for
#  the variance of the conditional mean of the treatment assignment variable.}
#  \item{cov_mybmt}{Estimated covariance between the conditional mean of the
#  treatment outcome variable and the bias of the treatment assignment
#  variable.} \item{cov_bmymt}{Estimated covariance between the bias of the
#  conditonal mean of the treatment variable and the conditional mean of the
#  treatment assignment variable.} \item{cov_bmybmt}{Estimated covariance
#  between the two biases of the conditional mean of the treatment outcome and
#  treatment assignment variables.} \item{nr_var}{Variance in the
#  null-restricted hypothesis testing.} \item{nr_var_adj}{Adjusted variance in
#  the null-restricted hypothesis testing.}
#
#@usage boundary_var_fuzzy(dat, t0 = 0, kernID = 0, left = TRUE, maxit = 100,
#  tol = 1.0e-4, para = 1, grainsize = 1, llr.residuals = TRUE,
#  ls.derivative = TRUE)

boundary_var_fuzzy <- function(dat, t0 = 0, kernID = 0, left = TRUE, maxit = 100, 
                               tol = 1.0e-4, para = 1, grainsize = 1, llr.residuals = TRUE,
                               ls.derivative = TRUE){
  
  n   = dim(dat$y)[1]
  x   = dat$x
  y   = dat$y
  z   = dat$z    
  tau = dat$tau
  h_y = dat$h_y
  h_t = dat$h_z   # t = z
  p   = dat$p
  q   = dat$q
  n_all = dat$n_all
  f0_hat = dat$f0_hat
  fd1 = dat$fd1
  
  vu0 = kmoments_yt(kernID = kernID, left = left, h_y = h_y, h_t = h_t)$vu0
  vu1 = kmoments_yt(kernID = kernID, left = left, h_y = h_y, h_t = h_t)$vu1
  vu2 = kmoments_yt(kernID = kernID, left = left, h_y = h_y, h_t = h_t)$vu2
  vu3 = kmoments_yt(kernID = kernID, left = left, h_y = h_y, h_t = h_t)$vu3
  vu4 = kmoments_yt(kernID = kernID, left = left, h_y = h_y, h_t = h_t)$vu4
  vu5 = kmoments_yt(kernID = kernID, left = left, h_y = h_y, h_t = h_t)$vu5
  vu6 = kmoments_yt(kernID = kernID, left = left, h_y = h_y, h_t = h_t)$vu6
  vu7 = kmoments_yt(kernID = kernID, left = left, h_y = h_y, h_t = h_t)$vu7
  vu8 = kmoments_yt(kernID = kernID, left = left, h_y = h_y, h_t = h_t)$vu8
  
  # Construct dat_y and dat_t
  dat_y = list("x" = as.matrix(x),
               "y" = as.matrix(y),
               "h" = h_y,
               "q" = q,
               "tau" = tau,
               "p" = p,
               "n_all" = n_all,
               "f0_hat" = f0_hat,
               "fd1" = fd1)
  
  dat_t = list("x" = as.matrix(x),
               "y" = as.matrix(z),
               "h" = h_t,
               "q" = q,
               "tau" = tau,
               "p" = p,
               "n_all" = n_all,
               "f0_hat" = f0_hat,
               "fd1" = fd1)
  
  var_mY = boundary_var(dat_y, kernID = kernID, left = left, maxit = maxit, tol = tol, para = para, grainsize = grainsize, llr.residuals = llr.residuals)
  var_mT = boundary_var(dat_t, kernID = kernID, left = left, maxit = maxit, tol = tol, para = para, grainsize = grainsize, llr.residuals = llr.residuals)
  
  e_y_hat = var_mY$e_hat
  ss_0_y  = var_mY$ss_0
  sig0_hat_y = var_mY$sig0_hat
  
  e_t_hat = var_mT$e_hat
  ss_0_t  = var_mT$ss_0
  sig0_hat_t = var_mT$sig0_hat
  
  c_k_y = stats::quantile(e_y_hat,tau)
  c_k_t = stats::quantile(e_t_hat,tau)
  
  eq = matrix(1,nrow = q,ncol = 1)
  
  e_yt_hat = cbind(e_y_hat,e_t_hat)
  cov_eig  = eigen(stats::cov(e_yt_hat))
  sqrtm    = cov_eig$vectors %*% diag(sqrt(cov_eig$values)) %*% solve(cov_eig$vectors)
  e_yt_hat = e_yt_hat %*% solve(sqrtm)
  
  phi_mat = matrix(0,nrow = q, ncol = q)
  for(i in 1:q) {
    for (j in 1:q) {
      phi_mat[i,j] = mean(e_yt_hat[,1] < c_k_y[i] & e_yt_hat[,2] < c_k_t[j]) - tau[i]*tau[j] 
    }
  }
  
  sg11  = phi_mat * vu0
  sg12  = as.matrix(rowSums(phi_mat)) * vu1 * h_y / h_t # The vu term for 12 block has an extra term h_y/h_t.
  sg21  = t(colSums(phi_mat)) * vu1
  sg22  = sum(phi_mat) * vu2
  sg_yt = rbind(cbind(sg11,sg12), cbind(sg21,sg22))
  byt   = t(eq) %*% (solve(ss_0_y) %*% sg_yt %*% solve(ss_0_t))[1:q,1:q] %*% eq / q^2
  cov_mymt = byt*sig0_hat_y*sig0_hat_t / (n_all*sqrt(h_y)*sqrt(h_t)*f0_hat)
  
  var_Y_tilde = var_mY$var + t0^2*var_mT$var - 2*t0*cov_mymt
  
  # Compute Var(Bias_Y)
  var_bias_mY = var_mY$var_bias
  ss_1_Y      = var_mY$ss_1
  ac_Y        = var_mY$ac
  var_bias_mT = var_mT$var_bias
  ss_1_T      = var_mT$ss_1
  ac_T        = var_mT$ac
  
  # Reconstruct Sigma_YT matrix with p = 2.
  p = 2
  sg11_1 = phi_mat * vu0 
  sg12_1 = cbind(rowSums(phi_mat)*vu1, rowSums(phi_mat)*vu2, rowSums(phi_mat)*vu3)[,1:p] * h_y / h_t
  sg21_1 = t(cbind(colSums(phi_mat)*vu1, colSums(phi_mat)*vu2, colSums(phi_mat)*vu3)[,1:p])
  vu_vec = c(vu1,vu2,vu3,vu4,vu5,vu6)
  sg22_1 = sum(phi_mat)*cbind(vu_vec[2:(p+1)], vu_vec[3:(p+2)], vu_vec[4:(p+3)])[,1:p]
  
  if(q == 1){
    sg_1 = rbind(c(sg11_1, sg12_1), cbind(t(sg21_1), sg22_1))
  }else{
    sg_1 = rbind(cbind(sg11_1,sg12_1), cbind(sg21_1,sg22_1))
  }
  
  sand   = solve(ss_1_Y) %*% sg_1 %*% solve(ss_1_T)
  e_2    = matrix(0,nrow = p, ncol = 1)
  e_2[2] = 1
  cov_bmYbmT = ac_Y*ac_T*sig0_hat_y*sig0_hat_t* t(e_2) %*% sand[(q+1):(q+p),(q+1):(q+p)] %*% e_2 /
    (n_all*sqrt(h_y)*sqrt(h_t)*f0_hat)
  
  var_bias = var_bias_mY + t0^2*var_bias_mT - 2*t0*cov_bmYbmT
  
  # Compute cov(Y_tilde, Bias).
  var_cov_Y = var_mY$var_cov
  var_cov_T = var_mT$var_cov
  
  cov_mYbmT = ac_T*sig0_hat_y*sig0_hat_t*sum(sand[1:q,(q+2)]) /
    (q*n_all*sqrt(h_y)*sqrt(h_t))
  
  sand = solve(ss_1_T) %*% sg_1 %*% solve(ss_1_Y)
  
  cov_mTbmY = ac_Y*sig0_hat_y*sig0_hat_t*sum(sand[1:q,(q+2)]) /
    (q*n_all*sqrt(h_y)*sqrt(h_t))
  
  cov_Y_tilde_bias = var_cov_Y - t0*cov_mYbmT - 
    t0*cov_mTbmY + t0^2*var_cov_T
  
  # The null-restricted version.
  nr_var_adj = var_Y_tilde + var_bias - 2*cov_Y_tilde_bias
  nr_var0 = var_Y_tilde
  
  # Add the fixed-n results, starting from line 126 "ss_0_y."
  ss_0_y_fn = var_mY$ss_0_fixedn
  ss_0_t_fn = var_mT$ss_0_fixedn
  ss_1_y_fn = var_mY$ss_1_fixedn
  ss_1_t_fn = var_mT$ss_1_fixedn
  
  # Compute cov_mymt_fn.
  x0    = 0
  u0_y  = (x - x0)/h_y
  w0_y  = kernel_fun(u0_y,kernID = kernID)
  u0_t  = (x - x0)/h_t
  w0_t  = kernel_fun(u0_t,kernID = kernID)
  mkykt = sum(w0_y * w0_t) / (n * sqrt(h_y * h_t))
  mkyktx = sum(w0_y * w0_t * u0_t) / (n * sqrt(h_y * h_t))
  mkyktx2 = sum(w0_y * w0_t * u0_t^2) / (n * sqrt(h_y * h_t))
  mkyktx3 = sum(w0_y * w0_t * u0_t^3) / (n * sqrt(h_y * h_t))
  mkyxkt = sum(w0_y * w0_t * u0_y) / (n * sqrt(h_y * h_t))
  mkyx2kt = sum(w0_y * w0_t * u0_y^2) / (n * sqrt(h_y * h_t))
  mkyx3kt = sum(w0_y * w0_t * u0_y^3) / (n * sqrt(h_y * h_t))
  mkyxktx = sum(w0_y * w0_t * u0_y * u0_t) / (n * sqrt(h_y * h_t))
  mkyxktx2 = sum(w0_y * w0_t * u0_y * u0_t^2) / (n * sqrt(h_y * h_t))
  mkyx2ktx = sum(w0_y * w0_t * u0_y^2 * u0_t) / (n * sqrt(h_y * h_t))
  mkyx2ktx2 = sum(w0_y * w0_t * u0_y^2 * u0_t^2) / (n * sqrt(h_y * h_t))
  
  sg11_fn = phi_mat * mkykt
  sg12_fn = as.matrix(rowSums(phi_mat)) * mkyktx
  sg21_fn = t(as.matrix(colSums(phi_mat))) * mkyxkt
  sg22_fn = sum(phi_mat) * mkyxktx
  sg_yt_fn = rbind(cbind(sg11_fn,sg12_fn), cbind(sg21_fn,sg22_fn))
  
  
  if (1 / kappa(ss_0_y_fn) < 1e-14){
    message("Reciprical of matrix condition is less than 1e-14. Generalized inverse is used.")
    ss_0_y_inv_fn = ginv(ss_0_y_fn)
  } else {
    ss_0_y_inv_fn = solve(ss_0_y_fn)
  }
  
  if (1 / kappa(ss_0_t_fn) < 1e-14){
    message("Reciprical of matrix condition is less than 1e-14. Generalized inverse is used.")
    ss_0_t_inv_fn = ginv(ss_0_t_fn)
  } else {
    ss_0_t_inv_fn = solve(ss_0_t_fn)
  }
  
  cov_mymt_fn = t(eq) %*% ss_0_y_inv_fn[1:q,] %*% sg_yt_fn %*% t(ss_0_t_inv_fn[1:q,]) %*% eq /
    (q^2 * n * sqrt(h_y * h_t))
  
  var_Y_tilde_fn = var_mY$var_fixedn + t0^2*var_mT$var_fixedn - 2*t0*cov_mymt_fn
  
  # Compute Var(Bias_Y_fn)
  var_bias_mY_fn = var_mY$var_bias_fixedn
  var_bias_mT_fn = var_mT$var_bias_fixedn
  
  # Reconstruct Sigma_YT_fn matrix with p = 2.
  p = 2
  sg11_1_fn = phi_mat * mkykt
  sg12_1_fn = cbind(rowSums(phi_mat)*mkyktx, rowSums(phi_mat)*mkyktx2, rowSums(phi_mat)*mkyktx3)[,1:p]
  sg21_1_fn = t(cbind(colSums(phi_mat)*mkyxkt, colSums(phi_mat)*mkyx2kt, colSums(phi_mat)*mkyx3kt)[,1:p])
  sg22_1_fn = sum(phi_mat)*cbind(c(mkyxktx,mkyx2ktx),c(mkyxktx2,mkyx2ktx2))
  
  if(q == 1){
    sg_1_fn = rbind(c(sg11_1_fn, sg12_1_fn), cbind(t(sg21_1_fn), sg22_1_fn))
  }else{
    sg_1_fn = rbind(cbind(sg11_1_fn,sg12_1_fn), cbind(sg21_1_fn,sg22_1_fn))
  }
  
  D_1_y_fn = var_mY$D_1_fn
  D_1_t_fn = var_mT$D_1_fn
  if (1 / kappa(ss_1_y_fn) < 1e-14){
    message("Reciprical of matrix condition is less than 1e-14. Generalized inverse is used.")
    ss_1_y_inv_fn = ginv(ss_1_y_fn)
  } else {
    ss_1_y_inv_fn = solve(ss_1_y_fn)
  }
  if (1 / kappa(ss_1_t_fn) < 1e-14){
    message("Reciprical of matrix condition is less than 1e-14. Generalized inverse is used.")
    ss_1_t_inv_fn = ginv(ss_1_t_fn)
  } else {
    ss_1_t_inv_fn = solve(ss_1_t_fn)
  }
  cov_bmYbmT_fn = 4 * D_1_y_fn * D_1_t_fn * t(e_2) %*% 
    (ss_1_y_inv_fn %*% sg_1_fn %*% ss_1_t_inv_fn)[(q+1):(q+p),(q+1):(q+p)] %*% e_2 / 
    (n * sqrt(h_y * h_t))
  
  var_bias_fn = var_bias_mY_fn + t0^2 * var_bias_mT_fn - 2 * t0 * cov_bmYbmT_fn
  
  # Compute cov(Y_tilde_fn, Bias)
  var_cov_Y_fn = var_mY$var_cov_fixedn
  var_cov_T_fn = var_mT$var_cov_fixedn
  
  cov_mYbmT_fn = 2 * D_1_t_fn * 
    sum((ss_1_y_inv_fn %*% sg_1_fn %*% ss_1_t_inv_fn)[1:q,(q+2)]) / (q * n * sqrt(h_y * h_t))
  cov_mTbmY_fn = 2 * D_1_y_fn *
    sum((ss_1_t_inv_fn %*% sg_1_fn %*% ss_1_y_inv_fn)[1:q,(q+2)]) / (q * n * sqrt(h_y * h_t))
  
  cov_Y_tilde_bias_fn = var_cov_Y_fn - t0*cov_mYbmT_fn - 
    t0*cov_mTbmY_fn + t0^2*var_cov_T_fn
  
  # The null-restricted version.
  nr_var_adj_fn = var_Y_tilde_fn + var_bias_fn - 2*cov_Y_tilde_bias_fn
  nr_var0_fn = var_Y_tilde_fn

  return(list(var_my     = var_mY$var, 
              var_mt     = var_mT$var,
              cov_mymt   = cov_mymt,
              var_adj_y  = var_mY$var_adj, 
              var_adj_t  = var_mT$var_adj,
              cov_mybmt  = cov_mYbmT, 
              cov_bmymt  = cov_mTbmY, 
              cov_bmybmt = cov_bmYbmT,
              nr_var     = nr_var0, 
              nr_var_adj = nr_var_adj,
              var_my_fixedn = var_mY$var_fixedn,
              var_mt_fixedn = var_mT$var_fixedn,
              cov_mymt_fixedn = cov_mymt_fn,
              var_adj_y_fixedn  = var_mY$var_adj_fixedn, 
              var_adj_t_fixedn  = var_mT$var_adj_fixedn,
              bias_y_fixedn     = var_mY$bias_fixedn,
              bias_t_fixedn     = var_mT$bias_fixedn,
              cov_mybmt_fixedn  = cov_mYbmT_fn, 
              cov_bmymt_fixedn  = cov_mTbmY_fn, 
              cov_bmybmt_fixedn = cov_bmYbmT_fn,
              nr_var_fixedn     = nr_var0_fn, 
              nr_var_adj_fixedn = nr_var_adj_fn)
              )
  }