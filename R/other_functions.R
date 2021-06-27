# kernel function
kernel_fun <- function(x, kernID = 0){
  switch(kernID + 1,
         {(1 - abs(x))*(abs(x) <= 1)          },
         {15/16*(1 - x^2)^2*(abs(x) <= 1)     },
         {0.75*(1 - x^2)*(abs(x) <= 1)        },
         {stats::dnorm(x,0,1)                 },
         {70/81*(1 - abs(x)^3)^3*(abs(x) <= 1)},
         {35/32*(1 - x^2)^3*(abs(x) <= 1)     },
         {0.5*(abs(x) <= 1)                   },
         {stop("Invalid kernel type.")        }
  )
}

# kernel moments
kmoments <- function(kernID = 0, left = TRUE){
  
  kern_fun_mu_0 <- function(x) x^0*kernel_fun(x, kernID = kernID)
  kern_fun_mu_1 <- function(x) x^1*kernel_fun(x, kernID = kernID)
  kern_fun_mu_2 <- function(x) x^2*kernel_fun(x, kernID = kernID)
  kern_fun_mu_3 <- function(x) x^3*kernel_fun(x, kernID = kernID)
  kern_fun_mu_4 <- function(x) x^4*kernel_fun(x, kernID = kernID)
  kern_fun_mu_5 <- function(x) x^5*kernel_fun(x, kernID = kernID)
  kern_fun_mu_6 <- function(x) x^6*kernel_fun(x, kernID = kernID)
  kern_fun_mu_7 <- function(x) x^7*kernel_fun(x, kernID = kernID)
  kern_fun_mu_8 <- function(x) x^8*kernel_fun(x, kernID = kernID)
  kern_fun_vu_0 <- function(x) x^0*kernel_fun(x, kernID = kernID)*kernel_fun(x, kernID = kernID)
  kern_fun_vu_1 <- function(x) x^1*kernel_fun(x, kernID = kernID)*kernel_fun(x, kernID = kernID)
  kern_fun_vu_2 <- function(x) x^2*kernel_fun(x, kernID = kernID)*kernel_fun(x, kernID = kernID)
  kern_fun_vu_3 <- function(x) x^3*kernel_fun(x, kernID = kernID)*kernel_fun(x, kernID = kernID)
  kern_fun_vu_4 <- function(x) x^4*kernel_fun(x, kernID = kernID)*kernel_fun(x, kernID = kernID)
  kern_fun_vu_5 <- function(x) x^5*kernel_fun(x, kernID = kernID)*kernel_fun(x, kernID = kernID)
  kern_fun_vu_6 <- function(x) x^6*kernel_fun(x, kernID = kernID)*kernel_fun(x, kernID = kernID)
  kern_fun_vu_7 <- function(x) x^7*kernel_fun(x, kernID = kernID)*kernel_fun(x, kernID = kernID)
  kern_fun_vu_8 <- function(x) x^8*kernel_fun(x, kernID = kernID)*kernel_fun(x, kernID = kernID)
  
  int_bound = 0.00001
  if(left){
    mu0 = stats::integrate(kern_fun_mu_0,-1,int_bound)$value
    mu1 = stats::integrate(kern_fun_mu_1,-1,int_bound)$value
    mu2 = stats::integrate(kern_fun_mu_2,-1,int_bound)$value
    mu3 = stats::integrate(kern_fun_mu_3,-1,int_bound)$value
    mu4 = stats::integrate(kern_fun_mu_4,-1,int_bound)$value
    mu5 = stats::integrate(kern_fun_mu_5,-1,int_bound)$value
    mu6 = stats::integrate(kern_fun_mu_6,-1,int_bound)$value
    mu7 = stats::integrate(kern_fun_mu_7,-1,int_bound)$value
    mu8 = stats::integrate(kern_fun_mu_8,-1,int_bound)$value
    vu0 = stats::integrate(kern_fun_vu_0,-1,int_bound)$value
    vu1 = stats::integrate(kern_fun_vu_1,-1,int_bound)$value
    vu2 = stats::integrate(kern_fun_vu_2,-1,int_bound)$value
    vu3 = stats::integrate(kern_fun_vu_3,-1,int_bound)$value
    vu4 = stats::integrate(kern_fun_vu_4,-1,int_bound)$value
    vu5 = stats::integrate(kern_fun_vu_5,-1,int_bound)$value
    vu6 = stats::integrate(kern_fun_vu_6,-1,int_bound)$value
    vu7 = stats::integrate(kern_fun_vu_7,-1,int_bound)$value
    vu8 = stats::integrate(kern_fun_vu_8,-1,int_bound)$value
  }else{
    mu0 = stats::integrate(kern_fun_mu_0,-int_bound,1)$value
    mu1 = stats::integrate(kern_fun_mu_1,-int_bound,1)$value
    mu2 = stats::integrate(kern_fun_mu_2,-int_bound,1)$value
    mu3 = stats::integrate(kern_fun_mu_3,-int_bound,1)$value
    mu4 = stats::integrate(kern_fun_mu_4,-int_bound,1)$value
    mu5 = stats::integrate(kern_fun_mu_5,-int_bound,1)$value
    mu6 = stats::integrate(kern_fun_mu_6,-int_bound,1)$value
    mu7 = stats::integrate(kern_fun_mu_7,-int_bound,1)$value
    mu8 = stats::integrate(kern_fun_mu_8,-int_bound,1)$value
    vu0 = stats::integrate(kern_fun_vu_0,-int_bound,1)$value
    vu1 = stats::integrate(kern_fun_vu_1,-int_bound,1)$value
    vu2 = stats::integrate(kern_fun_vu_2,-int_bound,1)$value
    vu3 = stats::integrate(kern_fun_vu_3,-int_bound,1)$value
    vu4 = stats::integrate(kern_fun_vu_4,-int_bound,1)$value
    vu5 = stats::integrate(kern_fun_vu_5,-int_bound,1)$value
    vu6 = stats::integrate(kern_fun_vu_6,-int_bound,1)$value
    vu7 = stats::integrate(kern_fun_vu_7,-int_bound,1)$value
    vu8 = stats::integrate(kern_fun_vu_8,-int_bound,1)$value
  }
  return(list(mu0 = mu0,
              mu1 = mu1,
              mu2 = mu2,
              mu3 = mu3,
              mu4 = mu4,
              mu5 = mu5,
              mu6 = mu6,
              mu7 = mu7,
              mu8 = mu8,
              vu0 = vu0,
              vu1 = vu1,
              vu2 = vu2,
              vu3 = vu3,
              vu4 = vu4,
              vu5 = vu5,
              vu6 = vu6,
              vu7 = vu7,
              vu8 = vu8))
}


# Compute the equivalent kernel for an interior point x.
equivalent.kernel <- function(kernID = 0, p = 1, vu = 0){
  lower = ifelse(kernID == 3, -4.0, -1.0)
  upper = ifelse(kernID == 3,  4.0,  1.0)
  s_mat = matrix(0, nrow = p + 1, ncol = p + 1)
  for (i in 1:(p+1)) {
    for (j in 1:(p+1)) {
      inte_fun <- function(x) x^(i + j - 2) * kernel_fun(x, kernID = kernID)
      s_mat[i,j] = stats::integrate(inte_fun, lower, upper)$value
    }
  }
  e_vu = matrix(0, nrow = p + 1, ncol = 1)
  e_vu[vu + 1] = 1
  ek_fun <- function(x){
    d = length(x)
    f_vec = matrix(0,nrow = d, ncol = 1)
    for (i in 1:d) {
      f_vec[i] = c(t(e_vu) %*% solve(s_mat) %*% x[i]^seq(0,p) * kernel_fun(x[i], kernID = kernID))
    }
    return(f_vec)
  }
  return(ek_fun)
}

# Compute the constant C_vp in equation (4.3) in Fan and Gijbels (1996)
c_vp <- function(kernID = 0, p = 1, vu = 0){
  lower   = ifelse(kernID == 3, -4.0, -1.0)
  upper   = ifelse(kernID == 3,  4.0,  1.0)
  eqk_fun <- equivalent.kernel(kernID, p = p, vu = vu)
  eqk_num <- function(x) eqk_fun(x)*eqk_fun(x)
  eqk_den <- function(x) x^(p + 1) * eqk_fun(x)
  numerator   = factorial(p + 1)^2 * (2 * vu + 1) * stats::integrate(eqk_num, lower, upper)$value
  denominator = 2 * (p + 1 - vu) * (stats::integrate(eqk_den, lower, upper)$value)^2
  c_value = (numerator/denominator)^(1/(2*p + 3))
  return(c_value)
}

# Compute the global polynomial regression to estimate the second derivative when p = 1.
polyreg <- function(y, x, p = 1, weight = rep(1, length(y))){
  xname    = paste("x^", 1:(p+3), sep = "")
  xname    = paste("I(",xname, ")")
  equation = stats::as.formula(paste("y ~ ", paste(xname, collapse= "+")))
  fit      = stats::lm(equation, data = data.frame(x,y), weights = weight)
  cf       = cumprod(1:(p + 3))
  deri     = stats::coef(fit)[[p + 2]] * cf[p + 1] + stats::coef(fit)[[p + 3]] * cf[p + 2] * x +
    stats::coef(fit)[[p + 4]] * cf[p + 3] * x^2/2
  return(list(residuals = stats::residuals(fit), derivatives = deri))
}

# Compute the rule-of-thumb bandwidth for local linear regression.
hrotlls <- function(y, x, p = 1, kernID = 0, weight = rep(1, length(y))){
  res_polyreg = polyreg(y, x, p = p, weight = weight)
  residual    = res_polyreg$residuals
  derivative  = res_polyreg$derivatives
  numerator   = mean(residual^2)
  denominator = sum(derivative^2)
  c_value     = c_vp(kernID = kernID, p = 1, vu = 0)
  h_rot_lls   = c_value * (numerator/denominator)^(1/(2*p + 3))
  return(h_rot_lls)
} 

# Generalized inverse function.
# Copied from the R MASS package.
ginv <- function(X, tol = sqrt(.Machine$double.eps))
{
  #
  # based on suggestions of R. M. Heiberger, T. M. Hesterberg and WNV
  #
  if(length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X)))
    stop("'X' must be a numeric or complex matrix")
  if(!is.matrix(X)) X <- as.matrix(X)
  Xsvd <- svd(X)
  if(is.complex(X)) Xsvd$u <- Conj(Xsvd$u)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  if (all(Positive)) Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
  else if(!any(Positive)) array(0, dim(X)[2L:1L])
  else Xsvd$v[, Positive, drop=FALSE] %*% ((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop=FALSE]))
}

# Kernel moments involving both Y and T.
kmoments_yt <- function(kernID = 0, left = TRUE, h_y = 1, h_t = 1){
  
  kern_fun_mu_0 <- function(x) x^0*kernel_fun(x, kernID = kernID)
  kern_fun_mu_1 <- function(x) x^1*kernel_fun(x, kernID = kernID)
  kern_fun_mu_2 <- function(x) x^2*kernel_fun(x, kernID = kernID)
  kern_fun_mu_3 <- function(x) x^3*kernel_fun(x, kernID = kernID)
  kern_fun_mu_4 <- function(x) x^4*kernel_fun(x, kernID = kernID)
  kern_fun_mu_5 <- function(x) x^5*kernel_fun(x, kernID = kernID)
  kern_fun_mu_6 <- function(x) x^6*kernel_fun(x, kernID = kernID)
  kern_fun_mu_7 <- function(x) x^7*kernel_fun(x, kernID = kernID)
  kern_fun_mu_8 <- function(x) x^8*kernel_fun(x, kernID = kernID)
  kern_fun_vu_0 <- function(x) x^0*kernel_fun(x, kernID = kernID)*kernel_fun(x * h_y / h_t, kernID = kernID) * sqrt(h_y / h_t)
  kern_fun_vu_1 <- function(x) x^1*kernel_fun(x, kernID = kernID)*kernel_fun(x * h_y / h_t, kernID = kernID) * sqrt(h_y / h_t)
  kern_fun_vu_2 <- function(x) x^2*kernel_fun(x, kernID = kernID)*kernel_fun(x * h_y / h_t, kernID = kernID) * sqrt(h_y / h_t)
  kern_fun_vu_3 <- function(x) x^3*kernel_fun(x, kernID = kernID)*kernel_fun(x * h_y / h_t, kernID = kernID) * sqrt(h_y / h_t)
  kern_fun_vu_4 <- function(x) x^4*kernel_fun(x, kernID = kernID)*kernel_fun(x * h_y / h_t, kernID = kernID) * sqrt(h_y / h_t)
  kern_fun_vu_5 <- function(x) x^5*kernel_fun(x, kernID = kernID)*kernel_fun(x * h_y / h_t, kernID = kernID) * sqrt(h_y / h_t)
  kern_fun_vu_6 <- function(x) x^6*kernel_fun(x, kernID = kernID)*kernel_fun(x * h_y / h_t, kernID = kernID) * sqrt(h_y / h_t)
  kern_fun_vu_7 <- function(x) x^7*kernel_fun(x, kernID = kernID)*kernel_fun(x * h_y / h_t, kernID = kernID) * sqrt(h_y / h_t)
  kern_fun_vu_8 <- function(x) x^8*kernel_fun(x, kernID = kernID)*kernel_fun(x * h_y / h_t, kernID = kernID) * sqrt(h_y / h_t)
  
  int_bound = 0.00001
  if(left){
    mu0 = stats::integrate(kern_fun_mu_0,-1,int_bound)$value
    mu1 = stats::integrate(kern_fun_mu_1,-1,int_bound)$value
    mu2 = stats::integrate(kern_fun_mu_2,-1,int_bound)$value
    mu3 = stats::integrate(kern_fun_mu_3,-1,int_bound)$value
    mu4 = stats::integrate(kern_fun_mu_4,-1,int_bound)$value
    mu5 = stats::integrate(kern_fun_mu_5,-1,int_bound)$value
    mu6 = stats::integrate(kern_fun_mu_6,-1,int_bound)$value
    mu7 = stats::integrate(kern_fun_mu_7,-1,int_bound)$value
    mu8 = stats::integrate(kern_fun_mu_8,-1,int_bound)$value
    vu0 = stats::integrate(kern_fun_vu_0,-1,int_bound)$value
    vu1 = stats::integrate(kern_fun_vu_1,-1,int_bound)$value
    vu2 = stats::integrate(kern_fun_vu_2,-1,int_bound)$value
    vu3 = stats::integrate(kern_fun_vu_3,-1,int_bound)$value
    vu4 = stats::integrate(kern_fun_vu_4,-1,int_bound)$value
    vu5 = stats::integrate(kern_fun_vu_5,-1,int_bound)$value
    vu6 = stats::integrate(kern_fun_vu_6,-1,int_bound)$value
    vu7 = stats::integrate(kern_fun_vu_7,-1,int_bound)$value
    vu8 = stats::integrate(kern_fun_vu_8,-1,int_bound)$value
  }else{
    mu0 = stats::integrate(kern_fun_mu_0,-int_bound,1)$value
    mu1 = stats::integrate(kern_fun_mu_1,-int_bound,1)$value
    mu2 = stats::integrate(kern_fun_mu_2,-int_bound,1)$value
    mu3 = stats::integrate(kern_fun_mu_3,-int_bound,1)$value
    mu4 = stats::integrate(kern_fun_mu_4,-int_bound,1)$value
    mu5 = stats::integrate(kern_fun_mu_5,-int_bound,1)$value
    mu6 = stats::integrate(kern_fun_mu_6,-int_bound,1)$value
    mu7 = stats::integrate(kern_fun_mu_7,-int_bound,1)$value
    mu8 = stats::integrate(kern_fun_mu_8,-int_bound,1)$value
    vu0 = stats::integrate(kern_fun_vu_0,-int_bound,1)$value
    vu1 = stats::integrate(kern_fun_vu_1,-int_bound,1)$value
    vu2 = stats::integrate(kern_fun_vu_2,-int_bound,1)$value
    vu3 = stats::integrate(kern_fun_vu_3,-int_bound,1)$value
    vu4 = stats::integrate(kern_fun_vu_4,-int_bound,1)$value
    vu5 = stats::integrate(kern_fun_vu_5,-int_bound,1)$value
    vu6 = stats::integrate(kern_fun_vu_6,-int_bound,1)$value
    vu7 = stats::integrate(kern_fun_vu_7,-int_bound,1)$value
    vu8 = stats::integrate(kern_fun_vu_8,-int_bound,1)$value
  }
  return(list(mu0 = mu0,
              mu1 = mu1,
              mu2 = mu2,
              mu3 = mu3,
              mu4 = mu4,
              mu5 = mu5,
              mu6 = mu6,
              mu7 = mu7,
              mu8 = mu8,
              vu0 = vu0,
              vu1 = vu1,
              vu2 = vu2,
              vu3 = vu3,
              vu4 = vu4,
              vu5 = vu5,
              vu6 = vu6,
              vu7 = vu7,
              vu8 = vu8))
}

