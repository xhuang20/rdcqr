#' @title Bias on the boundary
#'
#' @description This function computes the bias on the boundary, typicall at x =
#'   0, in a local composite quantile regression.
#'
#' @param x0 A scalar, boundary point in estimation.
#' @param dat A list with the following components: \describe{ \item{y}{A vector
#'   treatment outcomes. In a fuzzy RD, this variable can also be a vector of
#'   treatment assignment variables.} \item{x}{A vector of covariates.}
#'   \item{h}{A scalar bandwidth.} \item{p}{The polynomial degree. p = 2 is used
#'   to estimate the second derivative of the conditional mean function for bias
#'   estimation. LCQR is used to estimate the second derivative.}
#'   \item{q}{Number of quantiles used in estimation. It needs to be an odd
#'   integer such as 5 or 9.} }
#'
#' @param kernID Kernel id number. \enumerate{ \item \code{kernID = 0}:
#'   triangular kernel. \item \code{kernID = 1}: biweight kernel. \item
#'   \code{kernID = 2}: Epanechnikov kernel. \item \code{kernID = 3}: Gaussian
#'   kernel. \item \code{kernID = 4}: tricube kernel. \item \code{kernID = 5}:
#'   triweight kernel. \item \code{kernID = 6}: uniform kernel. }
#' @param left A logical variable that takes the value \code{TRUE} for data to
#'   the left of (below) the cutoff. Defaults to \code{TRUE}.
#' @param maxit Maximum iteration number in the MM algorithm for quantile
#'   estimation. Defaults to 20.
#' @param tol  Convergence criterion in the MM algorithm. Defaults to 1.0e-3.
#' @return Estimated bias on the boundary.
#' @usage boundary_bias(x0 = 0, dat, kernID = 0, left = TRUE, maxit = 20, tol =
#'   1e-3)
boundary_bias <- function(x0 = 0, dat, kernID = 0, left = TRUE, maxit = 20, tol = 1e-3){
  
  n   = length(dat$y)
  x   = dat$x
  y   = dat$y
  h   = dat$h
  p   = dat$p
  q   = dat$q
  tau = (1:q)/(q+1)
  mu0 = kmoments(kernID = kernID, left = left)$mu0
  mu1 = kmoments(kernID = kernID, left = left)$mu1
  mu2 = kmoments(kernID = kernID, left = left)$mu2
  mu3 = kmoments(kernID = kernID, left = left)$mu3

  est_cqr = cqrMMcpp(x0 = x0, x, y, kernID = kernID, tau, h = h, p = 2, maxit = maxit, tol = tol)
  md2     = est_cqr$beta1[2]*2
  bias    = 0.5*(mu2^2 - mu1*mu3)*md2*(h^2)/(mu0*mu2 - mu1^2)
  return(bias)
}