
<!-- README.md is generated from README.Rmd. Please edit that file -->
rdcqr
=====

<!-- badges: start -->
<!-- badges: end -->
Regression discontinuity (RD) design is a popular quasi-experiment method to identify the local average treatment effect. **rdcqr** applies the local composite quantile regression method to estimate the treatment effect in both sharp and fuzzy RD designs. The package can be used to perform the following analysis.

-   Estimation of local treatment effect for both sharp and fuzzy RD designs, with and without bias correction.

-   Hypothesis testing on the treatment effect. Adjusted standard errors are computed to incorporate the additional variability due to bias-correction.

-   A modified *t* test that eliminates the asymptotic size distortion due to weak identification in a fuzzy RD.

Installation
------------

You can install the development version from [GitHub](https://github.com/xhuang20/rdcqr.git) with:

``` r
# install.packages("devtools")
devtools::install_github("xhuang20/rdcqr")
```

Examples
--------

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
rdcqr(y, x)
```
