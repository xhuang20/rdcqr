#'@title Bandwidth computation
#'
#'@description This function computes several bandwidths for local composite
#'  quantile estimation on the boundary.
#'
#'@param dat A list with the following components: \describe{ \item{y}{A vector
#'  treatment outcomes. In a fuzzy RD, this variable can also be a vector of
#'  treatment assignment variables.} \item{x}{A vector of covariates.}
#'  \item{q}{Number of quantiles to be used in estimation. Defaults to 5. It
#'  needs to be an odd number.} \item{h}{A scalar bandwidth.}
#'  \item{tau}{Quantile positions that correspond to the q quantiles. They are
#'  obtained by \code{tau = (1:q)/(q+1)}.} \item{p}{The degree of polynomial in
#'  LCQR estimation. Set it to 1 when estimating the residuals.}
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
#'  estimation. Defaults to 20.
#'@param tol  Convergence criterion in the MM algorithm. Defaults to 1.0e-3.
#'@param para A 0/1 variable specifying whether to use parallel computing.
#'  Defaults to 1.
#'@param grainsize Minimum chunk size for parallelization. Defaults to 1.
#'@param llr.residuals Whether to use residuals from the local linear regression
#'  as the input to compute the LCQR standard errors and the corresponding
#'  bandwidths. Defaults to \code{TRUE}. If this option is set to \code{TRUE},
#'  the treatment effect estimate and the bias-correction is still done in LCQR.
#'  We use the same kernel function used in LCQR in the local linear regression
#'  to obtain the residuals and use them to compute the unadjusted and adjusted
#'  asymptotic standard errors and the bandwidths. This option will improve the
#'  speed. One can use this option to get a quick estimate of the standard
#'  errors when the sample size is large. To use residuals from the LCQR method,
#'  set \code{llr.residuals = FALSE}.
#'@param ls.derivative Whether to use a global quartic and quintic polynomial to
#'  estimate the second and third derivatives of the conditional mean function.
#'  Defaults to \code{TRUE}.
#'
#'@return \code{h_cqr} returns a list with the following components:
#'  \item{h_mse}{MSE-optimal bandwidth on the boundary. This bandwidth has order
#'  \eqn{O(n^{-1/5})}, where n is the sample size for data either below or above
#'  the cutoff.} \item{h_opt}{Bandwidth based on the adjusted MSE on the
#'  boundary. See Huang and Zhan (2020) for details about the adjusted MSE. This
#'  bandwidth has order \eqn{O(n^{-1/7})}, where n is the sample size for data
#'  either below or above the cutoff.} \item{h_rot}{A transform of the
#'  rule-of-thumb bandwidth for the local linear regression. The rule-of-thumb
#'  bandwidth is close to the Mean Integrated Squared Error optimal
#'  (MISE-optimal) bandwidth in a local linear regression. This bandwidth has
#'  order \eqn{O(n^{-1/5})}, where n is the sample size for data either below or
#'  above the cutoff.}
#'
#'@export
#'
#'@usage h_cqr(dat, kernID = 0, left = TRUE, maxit = 20, tol = 1.0e-3, para = 1,
#'  grainsize = 1, llr.residuals = TRUE, ls.derivative = TRUE)
#'
#'@examples
#' \dontrun{
#' # Use the headstart data.
#' data(headstart)
#' data_n = subset(headstart, headstart$poverty < 0)
#' p = 1
#' q = 5
#' tau = (1:q) / (q + 1)
#' x = headstart$poverty
#' h_d0   = ks::hns(x, deriv.order = 0)
#' f0_hat = ks::kdde(x, h = h_d0, deriv.order = 0, eval.points = c(0))$estimate
#' h_d1   = ks::hns(x, deriv.order = 1)
#' fd1    = ks::kdde(x, h = h_d1, deriv.order = 1, eval.points = c(0))$estimate
#'
#' # Set up the list to be passed to the h_cqr function.
#' # Supply a bandwidth that is equal to 3.5.
#' dat_n = list("x"      = data_n$poverty,
#'              "y"      = data_n$mortality,
#'              "q"      = q,
#'              "h"      = 3.5,
#'              "tau"    = tau,
#'              "p"      = p,
#'              "n_all"  = n_all,
#'              "f0_hat" = f0_hat,
#'              "fd1"    = fd1)
#'
#' # Use the residuals from local linear regression for a quick try.
#' h_cqr(dat = dat_n, llr.residuals = TRUE)
#' }
#'
#'@references{
#'
#'\cite{Huang and Zhan (2020) "Local Composite Quantile Regression for
#'Regression Discontinuity," working paper.}
#'
#'}
h_cqr <- function(dat, kernID = 0, left = TRUE, maxit = 20, tol = 1.0e-3, para = 1, grainsize = 1, llr.residuals = TRUE, ls.derivative = TRUE){
  
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
    sig_hat = sqrt(stats::fitted(locfit::locfit(u_hat^2 ~ locfit::lp(x, h = h, deg = 0), kern = kern_llr))) + .Machine$double.eps
    e_hat   = u_hat / sig_hat
  }else{
    est_hat = est_cqr(x, y, kernID = kernID, tau, h, p = p, maxit = maxit, tol = tol,
                      parallel = para, grainsize = grainsize)
    y_hat   = est_hat$y_hat
    u_hat   = est_hat$u_hat
    sig_hat = est_hat$sig_hat
    e_hat   = est_hat$e_hat
  }
  
  x0       = 0
  u0       = (x - x0)/h
  w0       = kernel_fun(u0, kernID = kernID)
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
  
  p = 1   # Set p = 1 when calculating var0.
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
  
  ac = (mu2^2 - mu1*mu3)/(mu0*mu2 - mu1^2)
  R1 = sum(tau_mat*f_mat)/ (q^2)
  bc = (mu2^2*vu0 - 2*mu1*mu2*vu1 + mu1^2*vu2)/(mu0*mu2 - mu1^2)^2
  bY2 = bc*R1
  
  a2c = (mu2^2 - mu1*mu4)/(mu0*mu2 - mu1^2)
  
  p = 2   # Set p = 2 when calculating var_adj.
  if(q == 1){
    s11_1 = diag(f_e) * mu0
  }else{
    s11_1  = diag(c(f_e))*mu0
  }
  #s11_1  = diag(c(f_e))*mu0
  s12_1  = cbind(f_e*mu1,f_e*mu2,f_e*mu3)[,1:p]
  s21_1  = t(s12_1)
  mu_vec = c(mu1,mu2,mu3,mu4,mu5,mu6)
  s22_1  = sum(f_e)*cbind(mu_vec[2:(p+1)], mu_vec[3:(p+2)], mu_vec[4:(p+3)])[,1:p]
  if(q == 1){
    ss_1 = rbind(c(s11_1,s12_1), cbind(t(s21_1), s22_1))
  }else{
    ss_1   = rbind(cbind(s11_1,s12_1),cbind(s21_1,s22_1))
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
  
  e_2      = matrix(0,nrow = p, ncol = 1)
  e_2[2]   = 1
  ac       = (mu2^2 - mu1*mu3)/(mu0*mu2 - mu1^2)
  var_bias = ac^2*sig0_hat^2*
    (t(e_2) %*% sand_1[(q+1):(q+p),(q+1):(q+p)] %*% e_2) /
    (n_all*h*f0_hat)
  var_cov  = ac*sig0_hat^2* sum(sand_1[1:q,(q+2)]) / (q*n_all*h*f0_hat)
  
  var_adj  = var0 + var_bias - 2*var_cov
  
  v1 = var0 * n_all * h
  v2 = var_bias * n_all * h
  v3 = var_cov * n_all * h
  B  = v1 + v2 - 2 * v3
  a3c = (mu2 * mu3 - mu1 * mu4) / (mu0 * mu2 - mu1^2)
  A  = 0.5 * a2c * fd1 / f0_hat * md2 + a3c * md3 / 6
  A2 = A^2
  h_opt = (B/6/A2)^(1/7)*(n_all)^(-1/7)
  h_mse = (sig0_hat^2 * bY / (ac^2 * md2^2 * f0_hat))^(0.2) * n_all^(-0.2)
  h_rot = h*R1^0.2
  
  return(list(h_mse = c(h_mse), 
              h_opt = c(h_opt), 
              h_rot = c(h_rot),
              A     = c(A),
              B     = c(B)))
}
