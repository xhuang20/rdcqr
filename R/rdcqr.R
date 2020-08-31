#' @title Local composite quantile estimation in regression discontinuity
#'
#' @description This function computes the local composite quantile regression
#'   (LCQR) estimator of treatment effect for both sharp and fuzzy regression
#'   discontinuity (RD) designs. It also computes the bias- corrected estimator
#'   and adjusts its standard error by incorporating the variability due to
#'   bias-correction.
#'
#' @param y A vector of treatment outcomes.
#' @param x A vector of covariates.
#' @param fuzzy A vector of treatment assignments in a fuzzy RD. Defaults to
#'   \code{NULL} in a sharp RD. 1 for receiving the treatment and 0 otherwise.
#' @param t0 Treatment effect under the null. Defaults to 0.
#' @param cutoff Cutoff for treatment assignment. Defaults to 0.
#' @param q Number of quantiles to be used in estimation. Defaults to 5. It
#'   needs to be an odd number.
#' @param bandwidth In a sharp RD, if the supplied bandwidth is a numeric vector
#'   of length two, the first element is the bandwidth for data below the cutoff
#'   and the second element is the bandwidth for data above the cutoff. In a
#'   fuzzy RD, the supplied bandwidth vector needs to have four elements in it:
#'   the first two bandwidths are for treatment outcomes below and above the
#'   cutoff and the last two bandwidths are for treatment assignments below and
#'   above the cutoff. If it is a string, the following types of bandwidth
#'   selector are provided: \enumerate{ \item \code{rot}: Rule-of-thumb
#'   bandwidth selector. Two bandwidths are used on each side of the cutoff,
#'   each of which is a transform of the rule-of-thumb bandwidth for the local
#'   linear regression.\item \code{adj.mseone}: One bandwidth based on the
#'   adjusted MSE function. \item \code{adj.msetwo}: Two bandwidths based on the
#'   adjusted MSE function. \item \code{msetwo}: Two bandwidths based on the MSE
#'   function, each of which is the MSE-optimal bandwidth.  }
#' @param kernel.type Kernel type that includes \enumerate{ \item
#'   \code{triangular}: triangular kernel. \code{kernID = 0}. \item
#'   \code{biweight}:   biweight kernel. \code{kernID = 1}. \item
#'   \code{epanechnikov}:  Epanechnikov kernel. \code{kernID = 2}. \item
#'   \code{gaussian}:  Gaussian kernel. \code{kernID = 3}. \item \code{tricube}:
#'   tricube kernel. \code{kernID = 4}. \item \code{triweight}:  triweight
#'   kernel. \code{kernID = 5}. \item \code{uniform}:  uniform kernel
#'   \code{kernID = 6}. }
#' @param maxit Maximum iteration number in the MM algorithm for quantile
#'   estimation. Defaults to 20.
#' @param tol  Convergence criterion in the MM algorithm. Defaults to 1.0e-3.
#' @param parallel A logical value specifying whether to use parallel computing.
#'   Defaults to  \code{TRUE}.
#' @param numThreads Number of threads used in parallel computation. The option
#'   \code{auto} uses all threads. Defaults to the number of cores minus 1.
#' @param grainsize Minimum chunk size for parallelization. Defaults to 1.
#' @param llr.residuals Whether to use residuals from the local linear
#'   regression as the input to compute the LCQR standard errors and the
#'   corresponding bandwidths. Defaults to \code{FALSE}. If this option is set
#'   to \code{TRUE}, the treatment estimate and the bias-correction is still
#'   done in LCQR. We use the same kernel function used in LCQR in the local
#'   linear regression to obtain the residuals and use them to compute the
#'   unadjusted and adjusted asymptotic standard errors and the bandwidths. This
#'   option will improve the speed but its accuracy is unclear. One can use this
#'   opition to get a quick estimate of the standard errors when the sample size
#'   is large.
#' @param ls.derivative Whether to use a global quartic polynomial to estimate
#'   the second derivative of the conditional mean function. Defaults to FALSE.
#'
#' @details This is the main function of the package and it estimates the
#'   treatment effect for both sharp and fuzzy RD designs. The LCQR estimate is
#'   obtained from an iterative algorithm and the estimation speed is much
#'   slower compared to that of the lcoal linear regression. Most computation
#'   time is spend on the calculation of the standard errors. By default, the
#'   code to compute the standard error and bandwidth is paralleled and the
#'   argument \code{numThreads} is set to the number of physical cores on a
#'   computer. The option \code{numThreads = "auto"} uses all available threads
#'   on a computer.
#'
#'   To further speed up things, the option \code{llr.residuals = TRUE}. This is
#'   particularly suitable when the sample size is large.
#'
#'   The two arguments \code{maxit} and \code{tol} have a large impact on the
#'   computation speed. For example, using \code{maxit = 500} and \code{tol =
#'   1e-6} will take much longer to complete compared to the default setting,
#'   though the results are more precise. Our limited experience with some of
#'   the popular RD data suggests that the treatment effect can usually be
#'   estimated precisely with low computation cost while the standard errors may
#'   have non-negligible change when one changes \code{maxit} and \code{tol}.
#'   This certainly depends on the data. One should experiment with different
#'   settings during estimation.
#'
#'   In estimating the bandwidths \code{adj.mseone}, \code{adj.msetwo}, and
#'   \code{msetwo}, we need an estimate for the second derivative of the
#'   conditional mean function. By default, the second derivative is estimated
#'   based on a quadratic polynomical in LCQR. The option \code{ls.derivative},
#'   when set to \code{TRUE}, estimates the second derivative using a simple,
#'   global quartic polynomial at the boundary point. Sometimes it can be
#'   difficult for a nonparametric method such as LCQR to estimate the second
#'   derivative accurately, and this option provides an alternative way for
#'   estimation.
#'
#' @return \code{rdcqr} returns a list with the following components:
#'   \item{estimate}{Treatment effect estimate and the bias-corrected treatment
#'   estimate} \item{se}{Asymptotic standard error and the adjusted asymptotic
#'   standard error} \item{bws}{Bandwidths used in estimation. There are two and
#'   four bandwidths for sharp and fuzzy RD, respectively. In a fuzzy RD, the
#'   first two bandwidths are associated with the treatment outcome variable
#'   below and above the cutoff. The last two bandwidths are associated with the
#'   treatment assignment variable below and above the curoff.} \item{nr_t}{The
#'   null-restricted t statistic to mitigate weak identification in a fuzzy RD.
#'   The second element is the bias-corrected and s.e.-adjusted version of this
#'   test.}
#'
#' @export
#'
#' @usage rdcqr(y, x, fuzzy = NULL, t0 = 0, cutoff = 0, q = 5, bandwidth =
#'   "rot", kernel.type = "triangular", maxit = 20, tol = 1.0e-3, parallel =
#'   TRUE, numThreads = "default", grainsize = 1, llr.residuals = FALSE,
#'   ls.derivative = FALSE)
#'
#' @examples
#' \dontrun{
#' # An example of using the Head Start data.
#' data(headstart)
#' y = headstart$mortality
#' x = headstart$poverty
#'
#' # Use the defaul rule-of-thumb bandwidth in estimation.
#' # Also use the residuals from a local linear regression to estimate the
#' # standard errors of LCQR.
#' rdcqr(y, x, bandwidth = "rot", llr.residuals = TRUE)
#'
#' # Supply bandwidths to data below and above the cutoff 0.
#' # The \code{poverty} variable is preprocessed to have a cutoff equal to 0.
#' rdcqr(y, x, bandwidth = c(10, 3), llr.residuals = TRUE)
#'
#' # Try the MSE-optimal bandwidths for data below and above the cutoff.
#' rdcqr(y, x, bandwidth = "msetwo", llr.residuals = TRUE)
#'
#' # Use the default settings.
#' # This approach uses residuals from a LCQR estimation when computing the
#' # standard errors. It is slow for large data sets. By default, the option
#' \code{parallel = TRUE} is used.
#' rdcqr(y, x)
#' }
rdcqr <- function(y, x, fuzzy = NULL, t0 = 0, cutoff = 0, q = 5, bandwidth = "rot", 
                  kernel.type = "triangular", maxit = 20, tol = 1.0e-3, 
                  parallel = TRUE, numThreads = "default", grainsize = 1,
                  llr.residuals = FALSE, ls.derivative = FALSE){
  
  switch(kernel.type,
         triangular   = {kernID = 0},
         biweight     = {kernID = 1},
         epanechnikov = {kernID = 2},
         gaussian     = {kernID = 3},
         tricube      = {kernID = 4},
         triweight    = {kernID = 5},
         uniform      = {kernID = 6},
         {stop("Invalid kernel type.")}
  )
  
  if(parallel) {
    para = 1
    if(is.character(numThreads)){
      if(numThreads == "default"){
        RcppParallel::setThreadOptions(numThreads = RcppParallel::defaultNumThreads()/2 - 1)
      } else if(numThreads == "auto") {
        RcppParallel::setThreadOptions(numThreads = "auto")
      } else {
        stop("Invalid numThreads type.")
      }
    }
    if(is.numeric(numThreads)){
      if (numThreads > RcppParallel::defaultNumThreads()) 
        message("The specified number of threads is larger than the number of logical processors on this computer.")
    }
  } else {
    para = 0
  }
  
  if(is.null(fuzzy)){
    data_all = data.frame("y" = y, "x" = x - cutoff)
    data_all = data_all[order(data_all$x),]
    x        = data_all$x
    y        = data_all$y
    data_p   = subset(data_all, data_all$x >= 0)
    data_n   = subset(data_all, data_all$x < 0 )
    tau      = (1:q)/(q+1)
    n_all    = dim(data_all)[1]
    # h_f0     = 1.06*stats::sd(x)*n_all^(-0.2)
    # f0_hat   = sum(exp(-0.5*(x/h_f0)^2)/(sqrt(2*pi)))/(n_all*h_f0)
    h_d0   = ks::hns(x, deriv.order = 0)
    fd0    = ks::kdde(x,h = h_d0, deriv.order = 0,eval.points = c(0))$estimate
    f0_hat = fd0
    h_d1   = ks::hns(x, deriv.order = 1)
    fd1    = ks::kdde(x,h = h_d1, deriv.order = 1,eval.points = c(0))$estimate
    
    p     = 1
    dat_n = list("x"      = as.matrix(data_n$x),
                 "y"      = as.matrix(data_n$y),
                 "q"      = q,
                 "tau"    = tau,
                 "p"      = p,
                 "n_all"  = n_all,
                 "f0_hat" = f0_hat,
                 "fd1"    = fd1)
    
    dat_p = list("x"      = as.matrix(data_p$x),
                 "y"      = as.matrix(data_p$y),
                 "q"      = q,
                 "tau"    = tau,
                 "p"      = p,
                 "n_all"  = n_all,
                 "f0_hat" = f0_hat,
                 "fd1"    = fd1)
    
    if(is.numeric(bandwidth)){
      if(length(bandwidth) == 1){
        h_n = bandwidth
        h_p = bandwidth
      }else if(length(bandwidth) == 2){
        h_n = bandwidth[1]
        h_p = bandwidth[2]
      }else{
        stop("Plese supply only two bandwidths.")
      }
    }else {
      if(!(bandwidth %in% c("rot", "adj.mseone", "adj.msetwo", "msetwo"))) {
        stop("The requested bandwidth selector is not programmed.")
      }
      h_pilot      = hrotlls(data_n$y,data_n$x, p = 1, kernID = kernID)
      dat_n[["h"]] = h_pilot
      h_n_set      = h_cqr(dat_n, kernID = kernID, left = TRUE, maxit = maxit, tol = tol, para = para, grainsize = grainsize, llr.residuals = llr.residuals, ls.derivative = ls.derivative) 
      h_pilot      = hrotlls(data_p$y, data_p$x, p = 1, kernID = kernID)
      dat_p[["h"]] = h_pilot
      h_p_set      = h_cqr(dat_p, kernID = kernID, left = FALSE, maxit = maxit, tol = tol, para = para, grainsize = grainsize, llr.residuals = llr.residuals, ls.derivative = ls.derivative) 
      A_n = h_n_set$A
      B_n = h_n_set$B
      A_p = h_p_set$A
      B_p = h_p_set$B
      h   = ((B_n + B_p)/(A_p - A_n)^2/6)^(1/7) * n_all^(-1/7)
      if(bandwidth == "rot"){
        h_n = h_n_set$h_rot
        h_p = h_p_set$h_rot
      }
      if(bandwidth == "adj.mseone"){
        h_n = h
        h_p = h
      }
      if(bandwidth == "adj.msetwo"){
        h_n = h_n_set$h_opt
        h_p = h_p_set$h_opt
      }
      if(bandwidth == "msetwo"){
        h_n = h_n_set$h_mse
        h_p = h_p_set$h_mse
      }
    }
    
    dat_n[["h"]] = h_n
    dat_p[["h"]] = h_p
    est_n  = cqrMMcpp(x0 = 0, dat_n$x, dat_n$y, kernID = kernID, tau, h = h_n, p = 1, maxit = maxit, tol = tol)
    est_p  = cqrMMcpp(x0 = 0, dat_p$x, dat_p$y, kernID = kernID, tau, h = h_p, p = 1, maxit = maxit, tol = tol)
    bias_n = boundary_bias(x0 = 0, dat_n, kernID = kernID, left = TRUE, maxit = maxit, tol = tol)
    bias_p = boundary_bias(x0 = 0, dat_p, kernID = kernID, left = FALSE, maxit = maxit, tol = tol)
    var_n  = boundary_var(dat_n, kernID = kernID, left = TRUE, maxit = maxit, tol = tol, para = para, grainsize = grainsize, llr.residuals = llr.residuals)
    var_p  = boundary_var(dat_p, kernID = kernID, left = FALSE, maxit = maxit, tol = tol, para = para, grainsize = grainsize, llr.residuals = llr.residuals)
    est    = mean(est_p$beta0) - mean(est_n$beta0)
    est_bc = mean(est_p$beta0) - bias_p - mean(est_n$beta0) + bias_n
    est_var     = var_n$var + var_p$var
    est_var_adj = var_n$var_adj + var_p$var_adj
    
    
    result_est           = cbind(est, est_bc)
    colnames(result_est) = c("est", "est_bc")
    result_se            = cbind(sqrt(est_var), sqrt(est_var_adj))
    colnames(result_se)  = c("se", "se_adj")
    result_bws           = cbind(h_n, h_p)
    colnames(result_bws) = c("h_y_n", "h_y_p")
    
    return(list(estimate = result_est, 
                se       = result_se,
                bws      = result_bws))
  }
  
  if(!is.null(fuzzy)){
    
    data_all = data.frame("y" = y, "x" = x - cutoff, "z" = fuzzy)
    data_all = data_all[order(data_all$x),]
    x        = data_all$x
    y        = data_all$y
    z        = data_all$z
    data_p   = subset(data_all, data_all$x >= 0)
    data_n   = subset(data_all, data_all$x <  0)
    tau      = (1:q)/(q+1)
    n_all    = dim(data_all)[1]
    h_f0     = 1.06*stats::sd(x)*n_all^(-0.2)
    f0_hat   = sum(exp(-0.5*(x/h_f0)^2)/(sqrt(2*pi)))/(n_all*h_f0)
    
    h_d0   = ks::hns(x, deriv.order = 0)
    fd0    = ks::kdde(x,h = h_d0, deriv.order = 0,eval.points = c(0))$estimate
    f0_hat = fd0
    h_d1   = ks::hns(x, deriv.order = 1)
    fd1    = ks::kdde(x,h = h_d1, deriv.order = 1,eval.points = c(0))$estimate
    
    p     = 1
    dat_n = list("x"      = as.matrix(data_n$x),
                 "y"      = as.matrix(data_n$y),
                 "z"      = as.matrix(data_n$z),
                 "q"      = q,
                 "tau"    = tau,
                 "p"      = p,
                 "n_all"  = n_all,
                 "f0_hat" = f0_hat,
                 "fd1"    = fd1)
    
    dat_p = list("x"      = as.matrix(data_p$x),
                 "y"      = as.matrix(data_p$y),
                 "z"      = as.matrix(data_p$z),
                 "q"      = q,
                 "tau"    = tau,
                 "p"      = p,
                 "n_all"  = n_all,
                 "f0_hat" = f0_hat,
                 "fd1"    = fd1)
    
    if(is.numeric(bandwidth)){
      if(length(bandwidth) == 4){
        h_y_n = bandwidth[1]
        h_y_p = bandwidth[2]
        h_z_n = bandwidth[3]
        h_z_p = bandwidth[4]
      }else if(length(bandwidth) < 4){
        stop("Please supply 4 bandwidths.")
      }else{
        stop("Plese supply no more than 4 bandwidths.")
      }
    }else {
      if(!(bandwidth %in% c("rot", "adj.mseone", "adj.msetwo", "msetwo"))) {
        stop("The requested bandwidth selector is not programmed.")
      }
      # Compute the bandwidth for the outcome variable y.
      h_pilot      = hrotlls(data_n$y,data_n$x, p = 1, kernID = kernID)
      dat_n[["h"]] = h_pilot
      h_n_set_y    = h_cqr(dat_n, kernID = kernID, left = TRUE, maxit = maxit, tol = tol, para = para, grainsize = grainsize, llr.residuals = llr.residuals, ls.derivative = ls.derivative) 
      h_pilot      = hrotlls(data_p$y, data_p$x, p = 1, kernID = kernID)
      dat_p[["h"]] = h_pilot
      h_p_set_y    = h_cqr(dat_p, kernID = kernID, left = FALSE, maxit = maxit, tol = tol, para = para, grainsize = grainsize, llr.residuals = llr.residuals, ls.derivative = ls.derivative) 
      A_n = h_n_set_y$A
      B_n = h_n_set_y$B
      A_p = h_p_set_y$A
      B_p = h_p_set_y$B
      h_y = ((B_n + B_p)/(A_p - A_n)^2/6)^(1/7) * n_all^(-1/7)
      
      # Repeat the above process to compute the bandwidth for the treatment variable z.
      h_pilot      = hrotlls(data_n$z,data_n$x, p = 1, kernID = kernID)
      dat_n_temp   = dat_n
      dat_n_temp[["y"]] = as.matrix(data_n$z)
      dat_n_temp[["h"]] = h_pilot
      h_n_set_z      = h_cqr(dat_n_temp, kernID = kernID, left = TRUE, maxit = maxit, tol = tol, para = para, grainsize = grainsize, llr.residuals = llr.residuals, ls.derivative = ls.derivative) 
      
      h_pilot      = hrotlls(data_p$y, data_p$x, p = 1, kernID = kernID)
      dat_p_temp   = dat_p
      dat_p_temp[["y"]] = as.matrix(data_p$z)
      dat_p_temp[["h"]] = h_pilot
      h_p_set_z    = h_cqr(dat_p_temp, kernID = kernID, left = FALSE, maxit = maxit, tol = tol, para = para, grainsize = grainsize, llr.residuals = llr.residuals, ls.derivative = ls.derivative) 
      A_n = h_n_set_z$A
      B_n = h_n_set_z$B
      A_p = h_p_set_z$A
      B_p = h_p_set_z$B
      h_z = ((B_n + B_p)/(A_p - A_n)^2/6)^(1/7) * n_all^(-1/7)
      
      if(bandwidth == "rot"){
        h_y_n = h_n_set_y$h_rot
        h_y_p = h_p_set_y$h_rot
        h_z_n = h_n_set_z$h_rot
        h_z_p = h_p_set_z$h_rot
      }
      if(bandwidth == "adj.mseone"){
        h_y_n = h_y
        h_y_p = h_y
        h_z_n = h_z
        h_z_p = h_z
      }
      if(bandwidth == "adj.msetwo"){
        h_y_n = h_n_set_y$h_opt
        h_y_p = h_p_set_y$h_opt
        h_z_n = h_n_set_z$h_opt
        h_z_p = h_p_set_z$h_opt
      }
      if(bandwidth == "msetwo"){
        h_y_n = h_n_set_y$h_mse
        h_y_p = h_p_set_y$h_mse
        h_z_n = h_n_set_z$h_mse
        h_z_p = h_p_set_z$h_mse
      }
    }
    
    dat_n[["h_y"]] = h_y_n
    dat_p[["h_y"]] = h_y_p
    dat_n[["h_z"]] = h_z_n
    dat_p[["h_z"]] = h_z_p
    est_n_y     = cqrMMcpp(x0 = 0, dat_n$x, dat_n$y, kernID = kernID, tau, h = dat_n$h_y, p = 1, maxit = maxit, tol = tol)
    est_p_y     = cqrMMcpp(x0 = 0, dat_p$x, dat_p$y, kernID = kernID, tau, h = dat_p$h_y, p = 1, maxit = maxit, tol = tol)
    bias_n_y    = boundary_bias(x0 = 0, dat_n, kernID = kernID, left = TRUE,  maxit = maxit, tol = tol)
    bias_p_y    = boundary_bias(x0 = 0, dat_p, kernID = kernID, left = FALSE, maxit = maxit, tol = tol)
    
    dat_n_temp        = dat_n
    dat_n_temp[["y"]] = dat_n$z
    dat_n_temp[["h"]] = dat_n$h_z
    dat_p_temp        = dat_p
    dat_p_temp[["y"]] = dat_p$z
    dat_p_temp[["h"]] = dat_n$h_z
    
    est_n_z  = cqrMMcpp(x0 = 0, dat_n$x, dat_n$z, kernID = kernID, tau, h = dat_n$h_z, p = 1, maxit = maxit, tol = tol)
    est_p_z  = cqrMMcpp(x0 = 0, dat_p$x, dat_p$z, kernID = kernID, tau, h = dat_p$h_z, p = 1, maxit = maxit, tol = tol)
    bias_n_z = boundary_bias(x0 = 0, dat_n_temp, kernID = kernID, left = TRUE,  maxit = maxit, tol = tol)
    bias_p_z = boundary_bias(x0 = 0, dat_p_temp, kernID = kernID, left = FALSE, maxit = maxit, tol = tol)
    
    est     = (mean(est_p_y$beta0) - mean(est_n_y$beta0)) / 
              (mean(est_p_z$beta0) - mean(est_n_z$beta0))
    est_adj = (mean(est_p_y$beta0) - bias_p_y - mean(est_n_y$beta0) + bias_n_y) / 
              (mean(est_p_z$beta0) - bias_p_z - mean(est_n_z$beta0) + bias_n_z)
    
    var_n  = boundary_var_fuzzy(dat_n, t0 = t0, kernID = kernID, left = TRUE,  maxit = maxit, tol = tol, para = para, grainsize = grainsize, llr.residuals = llr.residuals)
    var_p  = boundary_var_fuzzy(dat_p, t0 = t0, kernID = kernID, left = FALSE, maxit = maxit, tol = tol, para = para, grainsize = grainsize, llr.residuals = llr.residuals)
    
    dif_y = mean(est_p_y$beta0) - mean(est_n_y$beta0)
    dif_z = mean(est_p_z$beta0) - mean(est_n_z$beta0)
    
    est_var = (var_p$var_my + var_n$var_my) / dif_z^2 +
              (var_p$var_mt + var_n$var_mt) * dif_y^2 / dif_z^4 -
              2 * (var_p$cov_mymt + var_n$cov_mymt) * dif_y / dif_z^3
    est_var_adj = (var_p$var_adj_y + var_n$var_adj_y) / dif_z^2 +
                  (var_p$var_adj_t + var_n$var_adj_t) * dif_y^2 / dif_z^4 -
              2 * (var_p$cov_mymt  - var_p$cov_mybmt - var_p$cov_bmymt + var_p$cov_bmybmt +
                   var_n$cov_mymt  - var_n$cov_mybmt - var_n$cov_bmymt + var_n$cov_bmybmt) * dif_y / dif_z^3
    # Compute the null-restricted t stat, with and without bias-correction and s.e.-adjustment.
    nr_nume     = mean(est_p_y$beta0) - t0 * mean(est_p_z$beta0) -
                 (mean(est_n_y$beta0) - t0 * mean(est_n_z$beta0)) 
    
    nr_nume_adj =  mean(est_p_y$beta0) - bias_p_y - 
                   t0 * mean(est_p_z$beta0) + t0 * bias_p_y -
                   (mean(est_n_y$beta0) - bias_n_y - 
                   t0 * mean(est_n_z$beta0) + t0 * bias_n_z) 
    
    nr_deno     = sqrt(var_n$nr_var + var_p$nr_var)
    nr_deno_adj = sqrt(var_n$nr_var_adj + var_p$nr_var_adj)
    
    result_est           = cbind(est, est_adj)
    colnames(result_est) = c("est", "est_bc")
    result_se            = cbind(sqrt(est_var), sqrt(est_var_adj))
    colnames(result_se)  = c("se", "se_adj")
    result_bws           = cbind(h_y_n, h_y_p, h_z_n, h_z_p)
    colnames(result_bws) = c("h_y_n", "h_y_p", "h_z_n", "h_z_p")
    result_nrt           = cbind(nr_nume/nr_deno, nr_nume_adj/nr_deno_adj)
    colnames(result_nrt) = c("nr_t_stat", "nr_t_stat_adj")
    
    return(list(estimate = result_est,
                se       = result_se,
                bws      = result_bws,
                nr_t     = result_nrt))
  } 
}