
<!-- README.md is generated from README.Rmd. Please edit that file -->
# rdcqr

<!-- badges: start -->
<!-- badges: end -->
Regression discontinuity (RD) design is a popular quasi-experiment method to identify the local average treatment effect. **rdcqr** applies the local composite quantile regression method to estimate the treatment effect in both sharp and fuzzy RD designs. The package can be used to perform the following analysis.

-   Estimation of local treatment effect for both sharp and fuzzy RD designs, with and without bias correction. Results based on asympotitic and fixed-n approximations are both available.

-   Hypothesis testing on the treatment effect. Adjusted standard errors are computed to incorporate the additional variability due to bias-correction.

-   A modified *t* test that eliminates the asymptotic size distortion due to weak identification in a fuzzy RD.

In addition, the package also offers functions to perform sharp RD analysis with covariates and to estimate the treatment effect in a sharp kink RD.

## Installation

You can install the development version from [GitHub](https://github.com/xhuang20/rdcqr.git) with:

``` r
# install.packages("devtools")
devtools::install_github("xhuang20/rdcqr")
```

## Examples

Use the data on elections to the U.S. House of Representatives in Lee (2008) as an example.

``` r
library(rdcqr)
data(lee)
x = lee$margin
y = lee$voteshare

# Example 1 
# Use residuals from a local linear regression to compute the standard errors.
rdcqr(y, x, llr.residuals = TRUE)

# Example 2
# Use residuals from a local composite regression to compute the standard errors.
rdcqr(y, x, llr.residuals = FALSE)

# Example 3
# Use additional covariates in a sharp RD.
n   = length(x)
p   = 1                              # Select a polynomial order.
q   = 7                              # Select the number of quantile positions.
tau = (1:q) / (q + 1)                # quantile positions
h   = 0.3                            # Try a bandwidth.
treat = as.numeric(x >= 0)           # Create the treatment dummy variable.
z   = matrix(rnorm(2*n), nrow = n)   # Simulate an ad hoc covariate matrix.
xcv = cbind(treat, treat*x, z)       # Make sure the "treat" variable is the first column.
est = cqrMMxcv(x0     = 0,           # The cutoff.
               x_vec  = x,        
               y      = y,  
               xcv    = xcv,  
               kernID = 0,           # Use the triangular kernel.
               tau    = tau,
               h      = h,
               p      = p)
tau_hat = est$beta1[p+1]             # This gives the sharp RD estimate.

# Example 4
# Estimate sharp kink RD treatment effect.
p = 2      # Use a quadratic polynomial to estimate the first derivative.
est_n = cqrMMcpp(x0     = 0,
                 x_vec  = x[x < 0],
                 y      = y[x < 0],
                 kernID = 0,
                 tau    = tau,
                 h      = h,
                 p      = p)
est_p = cqrMMcpp(x0     = 0,
                 x_vec  = x[x >= 0],
                 y      = y[x >= 0],
                 kernID = 0,
                 tau    = tau,
                 h      = h,
                 p      = p)
kink_hat = est_p$beta1[1] - est_n$beta1[1] # This gives the sharp kink RD estimate.
```
