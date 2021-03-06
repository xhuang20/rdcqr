% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boundary_var.R
\name{boundary_var}
\alias{boundary_var}
\title{Variance on the boundary}
\usage{
boundary_var(dat, kernID = 0, left = TRUE, maxit = 100, tol = 1.0e-4,
 para = 1, grainsize = 1, llr.residuals = TRUE, ls.derivative = TRUE)
}
\arguments{
\item{dat}{A list with the following components: \describe{ \item{y}{A vector
of treatment outcomes, In a fuzzy RD, this variable can also be a vector of
treatment assignment variables.} \item{x}{A vector of covariates.}
\item{q}{Number of quantiles to be used in estimation. Defaults to 5. It
needs to be an odd number.} \item{h}{A scalar bandwidth.}
\item{tau}{Quantile positions that correspond to the q quantiles. They are
obtained by \code{tau = (1:q)/(q+1)}.} \item{p}{Degree of polynomial in LCQR
estimation. Set it to 1 when estimating the residuals.}
\item{n_all}{Total number of observations in the data set.}
\item{f0_hat}{Estimated density of the covariate at x = 0.}
\item{fd1}{Estimated first derivative of the density of the covariate at x =
0.} }}

\item{kernID}{Kernel id number. Defaults to 0. \enumerate{ \item \code{kernID
 = 0}: triangular kernel. \item \code{kernID = 1}:   biweight kernel. \item
\code{kernID = 2}:  Epanechnikov kernel. \item \code{kernID = 3}:  Gaussian
kernel. \item \code{kernID = 4}:  tricube kernel. \item \code{kernID = 5}:
triweight kernel. \item \code{kernID = 6}:  uniform kernel. }}

\item{left}{A logical variable that takes the value \code{TRUE} for data to
the left of (below) the cutoff. Defaults to \code{TRUE}.}

\item{maxit}{Maximum iteration number in the MM algorithm for quantile
estimation. Defaults to 100.}

\item{tol}{Convergence criterion in the MM algorithm. Defaults to 1.0e-4.}

\item{para}{A 0/1 variable specifying whether to use parallel computing.
Defaults to 1.}

\item{grainsize}{Minimum chunk size for parallelization. Defaults to 1.}

\item{llr.residuals}{Whether to use residuals from the local linear
regression as the input to compute the LCQR standard errors and the
corresponding bandwidths. Defaults to \code{TRUE}. If this option is set to
\code{TRUE}, the treatment effect estimate and the bias-correction is still done
in LCQR. We use the same kernel function used in LCQR in the local linear
regression to obtain the residuals and use them to compute the unadjusted
and adjusted asymptotic standard errors and the bandwidths. This option
will improve the speed. One can use this opition to get a quick estimate of
the standard errors when the sample size is large. To use residuals from
the LCQR method, set \code{llr.residuals = FALSE}.}

\item{ls.derivative}{Whether to use a global quartic and quintic polynomial to
estimate the second and third derivatives of the conditional mean function.
Defaults to \code{TRUE}.}
}
\value{
\code{boundary_var} returns a list with the following components:
\item{var}{Estimated variance on the boundary.} \item{var_adj}{Adjusted
variance estimate on the boundary.} \item{e_hat}{Residuals adjusted by the
estimated standard deviations.} \item{ss_0}{The S matrix in Huang and Zhan
(2021) when the degree of the polynomial in estimation is equal to 1.}
\item{ss_1}{The S matrix in Huang and Zhan (2021) when the degree of the
polynomial in estimation is equal to 2.} \item{sig0_hat}{Estimated
conditional standard deviation at the cutoff 0.} \item{var_bias}{Estimated
variance of the bias.} \item{var_cov}{Estimated covariance between the
conditional mean and its bias.} \item{ac}{The constant in the bias
expression.} \item{bias_fixedn}{The bias based on fixed-n approximation.}
There are also several returned values with suffix "fixedn", which are the
counterparts of the above values based on the fixed-n approximation.
}
\description{
This function computes variance and adjusted variance on the
boundary for a sharp RD using the local composite quantile regression
method. It also returns several other quantities that will be used in
computing the variance in a fuzzy RD.
}
\examples{
\dontrun{
# Use the headstart data.
data(headstart)
data_p = subset(headstart, headstart$poverty > 0)
p = 1
q = 5
tau = (1:q) / (q + 1)
h_d0   = ks::hns(x, deriv.order = 0)
f0_hat = ks::kdde(x, h = h_d0, deriv.order = 0, eval.points = c(0))$estimate
h_d1   = ks::hns(x, deriv.order = 1)
fd1    = ks::kdde(x, h = h_d1, deriv.order = 1, eval.points = c(0))$estimate

# Set up the list to be passed to the boundary_var function.
# Supply a bandwidth equal to 2.0.
dat_p = list("x"      = data_p$poverty,
             "y"      = data_p$mortality,
             "q"      = q,
             "h"      = 2.0,
             "tau"    = tau,
             "p"      = p,
             "n_all"  = n_all,
             "f0_hat" = f0_hat,
             "fd1"    = fd1)

# Use the residuals from local linear regression for a quick try.
boundary_var(dat = dat_p, left = FALSE, llr.residuals = TRUE, ls.derivative = TRUE)
}

}
\references{
{

\cite{Huang and Zhan (2021) "Local Composite Quantile Regression for
Regression Discontinuity," working paper.}

}
}
