#'@title Variance on the boundary
#'
#'@description This function computes variance and adjusted variance on the
#'  boundary for a sharp RD using the local composite quantile regression
#'  method. It also returns several other quantities that will be used in
#'  computing the variance in a fuzzy RD.
#'
#'@param dat A list with the following components: \describe{ \item{y}{A vector
#'  of treatment outcomes, In a fuzzy RD, this variable can also be a vector of
#'  treatment assignment variables.} \item{x}{A vector of covariates.}
#'  \item{q}{Number of quantiles to be used in estimation. Defaults to 5. It
#'  needs to be an odd number.} \item{h}{A scalar bandwidth.}
#'  \item{tau}{Quantile positions that correspond to the q quantiles. They are
#'  obtained by \code{tau = (1:q)/(q+1)}.} \item{p}{Degree of polynomial in LCQR
#'  estimation. Set it to 1 when estimating the residuals.}
#'  \item{n_all}{Total number of observations in the data set.}
#'  \item{f0_hat}{Estimated density of the covariate at x = 0.}
#'  \item{fd1}{Estimated first derivative of the density of the covariate at x =
#'  0.} }
#'@param kernID Kernel id number. Defaults to 0. \enumerate{ \item \code{kernID
#'  = 0}: triangular kernel. \item \code{kernID = 1}:   biweight kernel. \item
#'  \code{kernID = 2}:  Epanechnikov kernel. \item \code{kernID = 3}:  Gaussian
#'  kernel. \item \code{kernID = 4}:  tricube kernel. \item \code{kernID = 5}:
#'  triweight kernel. \item \code{kernID = 6}:  uniform kernel. }
#'@param left A logical variable that takes the value \code{TRUE} for data to
#'  the left of (below) the cutoff. Defaults to \code{TRUE}.
#'@param maxit Maximum iteration number in the MM algorithm for quantile
#'  estimation. Defaults to 100.
#'@param tol  Convergence criterion in the MM algorithm. Defaults to 1.0e-4.
#'@param para A 0/1 variable specifying whether to use parallel computing.
#'  Defaults to 1.
#'@param grainsize Minimum chunk size for parallelization. Defaults to 1.
#'@param llr.residuals Whether to use residuals from the local linear
#'   regression as the input to compute the LCQR standard errors and the
#'   corresponding bandwidths. Defaults to \code{TRUE}. If this option is set to
#'   \code{TRUE}, the treatment effect estimate and the bias-correction is still done
#'   in LCQR. We use the same kernel function used in LCQR in the local linear
#'   regression to obtain the residuals and use them to compute the unadjusted
#'   and adjusted asymptotic standard errors and the bandwidths. This option
#'   will improve the speed. One can use this opition to get a quick estimate of
#'   the standard errors when the sample size is large. To use residuals from
#'   the LCQR method, set \code{llr.residuals = FALSE}.
#'@param ls.derivative Whether to use a global quartic and quintic polynomial to
#'  estimate the second and third derivatives of the conditional mean function.
#'  Defaults to \code{TRUE}.
#'
#'@return \code{boundary_var} returns a list with the following components:
#'  \item{var}{Estimated variance on the boundary.} \item{var_adj}{Adjusted
#'  variance estimate on the boundary.} \item{e_hat}{Residuals adjusted by the
#'  estimated standard deviations.} \item{ss_0}{The S matrix in Huang and Zhan
#'  (2021) when the degree of the polynomial in estimation is equal to 1.}
#'  \item{ss_1}{The S matrix in Huang and Zhan (2021) when the degree of the
#'  polynomial in estimation is equal to 2.} \item{sig0_hat}{Estimated
#'  conditional standard deviation at the cutoff 0.} \item{var_bias}{Estimated
#'  variance of the bias.} \item{var_cov}{Estimated covariance between the
#'  conditional mean and its bias.} \item{ac}{The constant in the bias
#'  expression.} \item{bias_fixedn}{The bias based on fixed-n approximation.}
#'  There are also several returned values with suffix "fixedn", which are the 
#'  counterparts of the above values based on the fixed-n approximation. 
#'
#'@export
#'
#'@usage boundary_var(dat, kernID = 0, left = TRUE, maxit = 100, tol = 1.0e-4,
#'  para = 1, grainsize = 1, llr.residuals = TRUE, ls.derivative = TRUE)
#'
#' @examples
#' \dontrun{
#' # Use the headstart data.
#' data(headstart)
#' data_p = subset(headstart, headstart$poverty > 0)
#' p = 1
#' q = 5
#' tau = (1:q) / (q + 1)
#' h_d0   = ks::hns(x, deriv.order = 0)
#' f0_hat = ks::kdde(x, h = h_d0, deriv.order = 0, eval.points = c(0))$estimate
#' h_d1   = ks::hns(x, deriv.order = 1)
#' fd1    = ks::kdde(x, h = h_d1, deriv.order = 1, eval.points = c(0))$estimate
#'
#' # Set up the list to be passed to the boundary_var function.
#' # Supply a bandwidth equal to 2.0.
#' dat_p = list("x"      = data_p$poverty,
#'              "y"      = data_p$mortality,
#'              "q"      = q,
#'              "h"      = 2.0,
#'              "tau"    = tau,
#'              "p"      = p,
#'              "n_all"  = n_all,
#'              "f0_hat" = f0_hat,
#'              "fd1"    = fd1)
#'
#' # Use the residuals from local linear regression for a quick try.
#' boundary_var(dat = dat_p, left = FALSE, llr.residuals = TRUE, ls.derivative = TRUE)
#' }
#'
#'@references{
#'
#'\cite{Huang and Zhan (2021) "Local Composite Quantile Regression for
#'Regression Discontinuity," working paper.}
#'
#'}
boundary_var <- function(dat, kernID = 0, left = TRUE, maxit = 100, tol = 1.0e-4, para = 1, 
                         grainsize = 1, llr.residuals = TRUE, ls.derivative = TRUE){
  
  n   = length(dat$y)
  x   = dat$x
  y   = dat$y
  tau = dat$tau
  h   = dat$h
  p   = dat$p
  q   = dat$q
  n_all = dat$n_all
  f0_hat = dat$f0_hat
  fd1 = dat$fd1
  
  mu0 = kmoments(kernID = kernID, left = left)$mu0
  mu1 = kmoments(kernID = kernID, left = left)$mu1
  mu2 = kmoments(kernID = kernID, left = left)$mu2
  mu3 = kmoments(kernID = kernID, left = left)$mu3
  mu4 = kmoments(kernID = kernID, left = left)$mu4
  mu5 = kmoments(kernID = kernID, left = left)$mu5
  mu6 = kmoments(kernID = kernID, left = left)$mu6
  mu7 = kmoments(kernID = kernID, left = left)$mu7
  mu8 = kmoments(kernID = kernID, left = left)$mu8
  vu0 = kmoments(kernID = kernID, left = left)$vu0
  vu1 = kmoments(kernID = kernID, left = left)$vu1
  vu2 = kmoments(kernID = kernID, left = left)$vu2
  vu3 = kmoments(kernID = kernID, left = left)$vu3
  vu4 = kmoments(kernID = kernID, left = left)$vu4
  vu5 = kmoments(kernID = kernID, left = left)$vu5
  vu6 = kmoments(kernID = kernID, left = left)$vu6
  vu7 = kmoments(kernID = kernID, left = left)$vu7
  vu8 = kmoments(kernID = kernID, left = left)$vu8
  
  switch(kernID + 1,
         {kern_llr = "tria" },
         {kern_llr = "bisq" },
         {kern_llr = "epan" },
         {kern_llr = "gauss"},
         {kern_llr = "tcub" },
         {kern_llr = "trwt" },
         {kern_llr = "rect" }
  )
  
  if(llr.residuals){
    y_hat   = stats::fitted(locfit::locfit(y ~ locfit::lp(x, h = h, deg = 1), kern = kern_llr))
    u_hat   = y - y_hat
    sig2_hat = stats::fitted(locfit::locfit(u_hat^2 ~ locfit::lp(x, h = h, deg = 0), kern = kern_llr))
    sig2_hat[sig2_hat <= 0] = .Machine$double.eps
    sig_hat = sqrt(sig2_hat)
    e_hat   = u_hat / sig_hat
  }else{
    est_hat = est_cqr(x, y, kernID = kernID, tau, h, p = 1, maxit = maxit, tol = tol,
                      parallel = para, grainsize = grainsize)
    y_hat   = est_hat$y_hat
    u_hat   = est_hat$u_hat
    sig_hat = est_hat$sig_hat
    e_hat   = est_hat$e_hat
  }
  
  x0       = 0
  u0       = (x - x0)/h
  w0       = kernel_fun(u0,kernID = kernID)
  sig0_hat = sqrt( sum(t(w0) %*% u_hat^2)/sum(w0) )
  
  c_k = stats::quantile(e_hat,tau)
  f_e = matrix(0,q,1)
  h_e = stats::density(e_hat)$bw
  for(i in 1:q){
    u      = (e_hat - c_k[i])/h_e
    w_ck   = kernel_fun(u, kernID = kernID)
    f_e[i] = sum(w_ck)/(n*h_e)
  }
  
  tau_mat = matrix(0,nrow = q, ncol = q)
  for (i in 1:q) {
    for (j in 1:q) {
      tau_mat[i,j] = min(tau[i],tau[j]) - tau[i]*tau[j]
    }
  }
  
  f_mat = matrix(0,nrow = q, ncol = q)
  for(i in 1:q){
    for(j in 1:q){
      f_mat[i,j] = 1/(f_e[i]*f_e[j])
    }
  }
  
  p = 1
  if(q == 1){
    s11_0 = diag(f_e) * mu0
  }else{
    s11_0 = diag(c(f_e))*mu0
  }
  s12_0  = cbind(f_e*mu1,f_e*mu2,f_e*mu3)[,1:p]
  s21_0  = t(s12_0)
  mu_vec = c(mu1,mu2,mu3,mu4,mu5,mu6)
  s22_0  = sum(f_e)*cbind(mu_vec[2:(p+1)], mu_vec[3:(p+2)], mu_vec[4:(p+3)])[,1:p]
  ss_0   = rbind(cbind(s11_0,s12_0),cbind(s21_0,s22_0))
  sg11_0 = tau_mat * vu0
  sg12_0 = cbind(rowSums(tau_mat)*vu1, rowSums(tau_mat)*vu2, rowSums(tau_mat)*vu3)[,1:p]
  sg21_0 = t(sg12_0)
  vu_vec = c(vu1,vu2,vu3,vu4,vu5,vu6)
  sg22_0 = sum(tau_mat)*cbind(vu_vec[2:(p+1)], vu_vec[3:(p+2)], vu_vec[4:(p+3)])[,1:p]
  sg_0   = rbind(cbind(sg11_0,sg12_0), cbind(sg21_0,sg22_0))
  eq     = matrix(1,nrow = q,ncol = 1)
  sand_0 = solve(ss_0) %*% sg_0 %*% solve(ss_0)
  bY     = t(eq) %*% sand_0[1:q,1:q] %*% eq/(q^2)
  var0   = bY*sig0_hat^2/(n_all*h*f0_hat)
  
  # Add fixed-n computation.
  if(ls.derivative){
    est_ls  = stats::lm(y ~ I(x) + I(x^2) + I(x^3) + I(x^4))
    md2     = est_ls$coefficients[3][[1]] * 2
    est_ls  = stats::lm(y ~ I(x) + I(x^2) + I(x^3) + I(x^4) + I(x^5))
    md3     = est_ls$coefficients[4][[1]] * 6
  }else{
    est_cqr = cqrMMcpp(x0 = 0, x, y, kernID = kernID, tau, h = h, p = 2, maxit = maxit, tol = tol)
    md2     = est_cqr$beta1[2]*2
    est_cqr = cqrMMcpp(x0 = 0, x, y, kernID = kernID, tau, h = h, p = 3, maxit = maxit, tol = tol)
    md3     = est_cqr$beta1[3]*6
  }
  
  p = 2
  if(q == 1){
    s11_1 = diag(f_e) * mu0
  }else{
    s11_1 = diag(c(f_e))*mu0
  }
  s12_1  = cbind(f_e*mu1,f_e*mu2,f_e*mu3)[,1:p]
  s21_1  = t(s12_1)
  mu_vec = c(mu1,mu2,mu3,mu4,mu5,mu6)
  s22_1  = sum(f_e)*cbind(mu_vec[2:(p+1)], mu_vec[3:(p+2)], mu_vec[4:(p+3)])[,1:p]
  if(q == 1){
    ss_1 = rbind(c(s11_1,s12_1), cbind(t(s21_1), s22_1))
  }else{
    ss_1 = rbind(cbind(s11_1,s12_1),cbind(s21_1,s22_1))
  }
  sg11_1 = tau_mat * vu0
  sg12_1 = cbind(rowSums(tau_mat)*vu1, rowSums(tau_mat)*vu2, rowSums(tau_mat)*vu3)[,1:p]
  sg21_1 = t(sg12_1)
  vu_vec = c(vu1,vu2,vu3,vu4,vu5,vu6)
  sg22_1 = sum(tau_mat)*cbind(vu_vec[2:(p+1)], vu_vec[3:(p+2)], vu_vec[4:(p+3)])[,1:p]
  if(q == 1){
    sg_1 = rbind(c(sg11_1, sg12_1), cbind(t(sg21_1), sg22_1))
  }else{
    sg_1 = rbind(cbind(sg11_1,sg12_1), cbind(sg21_1,sg22_1))
  }
  eq     = matrix(1,nrow = q,ncol = 1)
  sand_1 = solve(ss_1) %*% sg_1 %*% solve(ss_1)
  e_2    = matrix(0,nrow = p, ncol = 1)
  e_2[2] = 1
  ac     = (mu2^2 - mu1*mu3)/(mu0*mu2 - mu1^2)
  var_bias = ac^2*sig0_hat^2*
    (t(e_2) %*% sand_1[(q+1):(q+p),(q+1):(q+p)] %*% e_2) /
    (n_all*h*f0_hat)
  var_cov  = ac*sig0_hat^2* sum(sand_1[1:q,(q+2)]) / (q*n_all*h*f0_hat)
  var_adj  = var0 + var_bias - 2*var_cov
  
  # Compute the fixed-n variance.
  
  mks    = sum(w0 / sig_hat) / (n * h)
  mkx1s = sum(w0 * u0 / sig_hat) / (n * h)
  mkx2s = sum(w0 * u0^2 / sig_hat) / (n * h)
  mkx3s = sum(w0 * u0^3 / sig_hat) / (n * h)
  mkx4s = sum(w0 * u0^4 / sig_hat) / (n * h)
  mkx5s = sum(w0 * u0^5 / sig_hat) / (n * h)
  mkx6s = sum(w0 * u0^6 / sig_hat) / (n * h)
  
  mk2    = sum(w0^2) / (n * h)
  mk2x1 = sum(w0^2 * u0) / (n * h)
  mk2x2 = sum(w0^2 * u0^2) / (n * h)
  mk2x3 = sum(w0^2 * u0^3) / (n * h)
  mk2x4 = sum(w0^2 * u0^4) / (n * h)
  mk2x5 = sum(w0^2 * u0^5) / (n * h)
  mk2x6 = sum(w0^2 * u0^6) / (n * h) 
  
  mksig = sum(w0 * (sig_hat - sig0_hat) / sig_hat) / (n * h)
  mksigx1 = sum(w0 * u0 * (sig_hat - sig0_hat) / sig_hat) / (n * h)
  mksigx2 = sum(w0 * u0^2 * (sig_hat - sig0_hat) / sig_hat) / (n * h)
  
  p = 1   # Set p = 1 when calculating var0. Repeat the earlier calculation steps.
  if(q == 1){
    s11_0_fn = diag(f_e) * mks
  }else{
    s11_0_fn = diag(c(f_e)) * mks
  }
  s12_0_fn = cbind(f_e * mkx1s, f_e * mkx2s, f_e * mkx3s)[,1:p]
  s21_0_fn = t(s12_0_fn)
  mu_vec_fn = c(mkx1s, mkx2s, mkx3s, mkx4s, mkx5s, mkx6s)
  s22_0_fn = sum(f_e) * cbind(mu_vec_fn[2:(p+1)], mu_vec_fn[3:(p+2)], mu_vec_fn[4:(p+3)])[,1:p]
  ss_0_fn  = rbind(cbind(s11_0_fn,s12_0_fn),cbind(s21_0_fn,s22_0_fn))
  
  sg11_0_fn = tau_mat * mk2
  sg12_0_fn = cbind(rowSums(tau_mat)*mk2x1, rowSums(tau_mat)*mk2x2, rowSums(tau_mat)*mk2x3)[,1:p]
  sg21_0_fn = t(sg12_0_fn)
  vu_vec_fn = c(mk2x1, mk2x2, mk2x3, mk2x4, mk2x5, mk2x6)
  sg22_0_fn = sum(tau_mat)*cbind(vu_vec_fn[2:(p+1)], vu_vec_fn[3:(p+2)], vu_vec_fn[4:(p+3)])[,1:p]
  sg_0_fn   = rbind(cbind(sg11_0_fn,sg12_0_fn), cbind(sg21_0_fn,sg22_0_fn))
  
  if( 1 / kappa(ss_0_fn) < 1e-14) {
    message("Reciprical of matrix condition is less than 1e-14. Generalized inverse is used.")
    ss_0_inv_fn = ginv(ss_0_fn)
  } else {
    ss_0_inv_fn = solve(ss_0_fn)
  }
  
  A_0_fn = t(eq) %*% (ss_0_inv_fn[1:q,1:q] %*% f_e * md2 * 0.5 * mkx2s +
                        ss_0_inv_fn[1:q,(q+1):(q+p)] * sum(f_e) * md2 * 0.5 * mkx3s) / q

  B_0_fn = t(eq) %*% (ss_0_inv_fn %*% sg_0_fn %*% ss_0_inv_fn)[1:q,1:q] %*% eq / q^2
  
  h_mse_fn = (B_0_fn / (4 * A_0_fn^2))^(0.2) * n^(-0.2)
  
  var0_fn = B_0_fn / (n * h)    # for computing variance(m - bias(m)).
  
  # Compute A_1_fn
  A_1_fn = t(eq) %*% (ss_0_inv_fn[1:q,1:q] %*% f_e * md3 * (1/6) * mkx3s +
                        ss_0_inv_fn[1:q,(q+1):(q+p)] * sum(f_e) * md3 * (1/6) * mkx4s) / q

  p = 2  # Set p = 2 when calculating var_adj.
  if(q == 1){
    s11_1_fn = diag(f_e) * mks
  }else{
    s11_1_fn  = diag(c(f_e)) * mks
  }
  s12_1_fn = cbind(f_e * mkx1s, f_e * mkx2s, f_e * mkx3s)[,1:p]
  s21_1_fn = t(s12_1_fn)
  s22_1_fn  = sum(f_e)*cbind(mu_vec_fn[2:(p+1)], mu_vec_fn[3:(p+2)], mu_vec_fn[4:(p+3)])[,1:p]
  if(q == 1){
    ss_1_fn = rbind(c(s11_1_fn,s12_1_fn), cbind(t(s21_1_fn), s22_1_fn))
  }else{
    ss_1_fn = rbind(cbind(s11_1_fn,s12_1_fn),cbind(s21_1_fn,s22_1_fn))
  }
  
  sg11_1_fn = tau_mat * mk2
  sg12_1_fn = cbind(rowSums(tau_mat)*mk2x1, rowSums(tau_mat)*mk2x2, rowSums(tau_mat)*mk2x3)[,1:p]
  sg21_1_fn = t(sg12_1_fn)
  sg22_1_fn = sum(tau_mat)*cbind(vu_vec_fn[2:(p+1)], vu_vec_fn[3:(p+2)], vu_vec_fn[4:(p+3)])[,1:p]
  
  if(q == 1){
    sg_1_fn = rbind(c(sg11_1_fn, sg12_1_fn), cbind(t(sg21_1_fn), sg22_1_fn))
  }else{
    sg_1_fn = rbind(cbind(sg11_1_fn,sg12_1_fn), cbind(sg21_1_fn,sg22_1_fn))
  }
  
  # Next compute var_bias_fn and var_cov_fn
  D_1_fn = A_0_fn / md2
  
  if (1 / kappa(ss_1_fn) < 1e-14){
    message("Reciprical of matrix condition is less than 1e-14. Generalized inverse is used.")
    ss_1_inv_fn = ginv(ss_1_fn)
  } else {
    ss_1_inv_fn = solve(ss_1_fn)
  }
  
  sand_1_fn = ss_1_inv_fn %*% sg_1_fn %*% ss_1_inv_fn
  var_bias_fn = 4 * D_1_fn^2 * (t(e_2) %*% sand_1_fn[(q+1):(q+p),(q+1):(q+p)] %*% e_2) /
    (n*h)
  var_cov_fn = 2 * D_1_fn * sum(sand_1_fn[1:q,(q+2)]) / (q * n *h)
  
  var_adj_fn = var0_fn + var_bias_fn - 2 * var_cov_fn
  
  # Compute and export the fixed-n bias up to O(h^2).
  bias_fixedn = A_0_fn * h^2
  
  return(list(var      = var0, 
              var_adj  = var_adj, 
              e_hat    = e_hat, 
              ss_0     = ss_0, 
              ss_1     = ss_1,
              sig0_hat = sig0_hat, 
              var_bias = var_bias, 
              var_cov  = var_cov, 
              ac       = ac,
              var_fixedn = var0_fn,
              var_adj_fixedn = var_adj_fn,
              bias_fixedn = bias_fixedn,
              ss_0_fixedn = ss_0_fn,
              ss_1_fixedn = ss_1_fn,
              var_bias_fixedn = var_bias_fn,
              D_1_fn = D_1_fn,
              var_cov_fixedn = var_cov_fn))
}
