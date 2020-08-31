// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include "RcppArmadillo.h"
#include <RcppParallel.h>

// Function to compute kernel weights
// 
arma::mat kernel_triangular(const arma::colvec & x) {
  arma::colvec m = (1 - arma::abs(x))%(arma::abs(x) <= 1);
  return m;
}

arma::mat kernel_biweight(const arma::colvec & x) {
  arma::colvec m = arma::square(1 - arma::square(x)) % (arma::abs(x) <= 1) * 15/16;
  return m;
}

arma::mat kernel_epanechnikov(const arma::colvec & x) {
  arma::colvec m = 0.75 * (1 - arma::square(x)) % (arma::abs(x) <= 1);
  return m;
}

arma::mat kernel_gaussian(const arma::colvec & x) {
  arma::colvec m =  ( 1 / ( 1.0 * sqrt(2*M_PI) ) ) * exp( -0.5 * arma::pow( (x - 0.0)/1.0, 2.0 ) );
  return m;
}

arma::mat kernel_tricube(const arma::colvec & x) {
  arma::colvec m = arma::pow((1 - arma::pow(arma::abs(x),3)),3) % (arma::abs(x) <= 1) * 70/81;
  return m;
}

arma::mat kernel_triweight(const arma::colvec & x) {
  arma::colvec m = arma::pow((1 - arma::square(x)),3) % (arma::abs(x) <= 1) * 35/32;
  return m;
}

arma::mat kernel_uniform(const arma::colvec & x) {
  arma::colvec m =  0.5 * arma::conv_to<arma::mat>::from(arma::abs(x) <= 1);
  return m;
}

typedef struct{
  arma::mat (*kern)(const arma::colvec &);
}kernels;

kernels kernPtr[] = {kernel_triangular, 
                     kernel_biweight, 
                     kernel_epanechnikov,
                     kernel_gaussian,
                     kernel_tricube,
                     kernel_triweight,
                     kernel_uniform};

//' @title Function to compute the local composite quantile regression estimate
//' 
//' @description This function computes the local composite quantile regression 
//' estimate. The point of interest can be either interior or boundary. 
//' 
//' @param x0 point of interest
//' @param x_vec a vector of covariates
//' @param y a vector of dependent variable, the treatment outcome variable 
//' in the case of regression discontinuity.
//' @param kernID kernel ID for different kernels.
//' \enumerate{
//'   \item \code{kernID = 0}: triangular kernel.
//'   \item \code{kernID = 1}: biweight kernel.
//'   \item \code{kernID = 2}: Epanechnikov kernel. 
//'   \item \code{kernID = 3}: Gaussian kernel.
//'   \item \code{kernID = 4}: tricube kernel.
//'   \item \code{kernID = 5}: triweight kernel. 
//'   \item \code{kernID = 6}: uniform kernel. 
//' }
//' @param tau A vector of quantile positions. They are obtained by 
//' \code{tau = (1:q)/(q+1)}.
//' @param h A scalar bandwidth.
//' @param p The polynomial degree.
//' @param maxit Maximum iteration number in the MM algorithm for quantile
//'  estimation.
//' @param tol Convergence criterion in the MM algorithm.
//' @export
//' @return \code{cqrMMcpp} returns a list with the following components:
//' \item{beta0}{A q by 1 vector of estimates for q quantiles.}
//' \item{beta1}{A p by 1 vector of estimates for the first p derivatives of
//' the conditional mean function.}
//' @examples
//' # Use the Head Start data as an example.
//' data(headstart)
//' data_n = subset(headstart, headstart$poverty < 0)
//' q      = 5
//' tau    = (1:q) / (q + 1)
//' 
//' # Compute the local composite quantile estimate
//' est = cqrMMcpp(x0     = 0,
//'                x_vec  = data_n$poverty,
//'                y      = data_n$mortality,
//'                kernID = 2,
//'                tau    = tau,
//'                h      = 4.0,
//'                p      = 1,
//'                maxit  = 10,
//'                tol    = 1.0e-3)
//'                
//' # Estimate of the conditional mean on the boundary
//' est_mean = mean(est$beta0)
//' 
//' # Estimate of the first derivative of the conditional mean function
//' est_d1   = est$beta1[1]
// [[Rcpp::export]]
Rcpp::List cqrMMcpp(double x0, const arma::colvec & x_vec, const arma::colvec & y, int kernID,
                    const arma::colvec & tau, double h, int p, 
                    int maxit = 20, double tol = 1e-3) {
  int n = arma::size(y)(0);
  int q = arma::size(tau)(0);
  int it;
  double delta, tn, e0, eps, sumw, sig;
  arma::mat x(n,p), x_mat(n,p+1), y_mat(n,1), p1(p,q), p5;
  arma::colvec beta(p + 1), res(n), b(p), w(n), a(q), c;
  arma::cube p6(p,p,q);
  arma::colvec xd = x_vec - x0;

  for (int j = 1; j < (p + 1); j++) {
    x.col(j - 1) = pow(xd,j);
  }
  w = kernPtr[kernID].kern(xd/h);
  if (arma::accu(w != 0) == 1) {
    Rcpp::stop("There is no data around the point of interest. Increase the bandwidth.");
  }
  
  /* d = p;*/
  x_mat = arma::join_rows(arma::ones(n,1),x) % arma::repmat(arma::sqrt(w), 1, (p+1));
  y_mat = y % sqrt(w);
  beta  = arma::solve(x_mat, y_mat); 
  res   = y_mat - x_mat*beta;
  sig   = sqrt(std::inner_product(res.begin(), res.end(), res.begin(), 0.0)/(n));
  b     = beta.rows(1,p);
  for (int j = 1; j <= q; j++) {
    a(j-1) = beta(0) + R::qnorm5(tau(j - 1), 0.0, sig, 1, 0);
  }
  c     = 2*tau - 1;
  delta = 1;
  it    = 0;
  tn    = tol/n;
  e0    = -tn/log(tn);
  eps   = (e0-tn)/(1+log(e0));
  sumw  = arma::as_scalar(arma::accu(w));
  p5    = trans(w) * x;

  arma::mat ar, A, p2, p3, p4, tp, denor, numer, aa;
  arma::colvec a0, b0;
  while(delta > tol && it < maxit){
    it++;
    b0 = b;
    a0 = a;
    ar = arma::abs(arma::repmat(y-x*b0,1,q) - arma::repmat(trans(a0),n,1));
    A  = 1 / (ar + eps);
    p2 = trans(y % w) * A;
    p3 = trans(x % arma::repmat(w,1,p)) * A;
    p4 = trans(w) * A;
    denor = arma::zeros(p,p);
    numer = arma::zeros(p,1);
    for(int k = 0; k < q; k++){
      tp          = arma::repmat(A.col(k) % w, 1, p) % x;
      p1.col(k)   = trans(tp) * y;
      p6.slice(k) = trans(tp) * x;
      denor      += p6.slice(k) - p3.col(k) * trans(p3.col(k)) / p4(k);
      numer      += p1.col(k) - p2(k) * p3.col(k) / p4(k) -
        sumw * c(k) * p3.col(k) / p4(k) + c(k) * trans(p5);
    }
    b = arma::solve(denor,numer);
    a = trans( (trans((y - x * b) % w) * A + sumw * trans(c) ) / p4 );
    delta = arma::as_scalar(arma::abs( arma::mean(a-a0,0)) + arma::sum(arma::abs(b-b0) ));
  }
  
  return Rcpp::List::create(Rcpp::Named("beta0")     = a,
                            Rcpp::Named("beta1")     = b,
                            Rcpp::Named("iteration") = it);
}

struct estimate {
  arma::colvec beta0;
  arma::colvec beta1;
};

estimate cqrMM(double x0, const arma::colvec & x_vec, const arma::colvec & y, int kernID,
                  const arma::colvec & tau, double h, int p, 
                  int maxit = 500, double tol = 1e-4) {
  int n = arma::size(y)(0);
  int q = arma::size(tau)(0);
  int it;
  double delta, tn, e0, eps, sumw, sig;
  arma::mat x(n,p), x_mat(n,p+1), y_mat(n,1), p1(p,q), p5;
  arma::colvec beta(p + 1), res(n), b(p), w(n), a(q), c;
  arma::cube p6(p,p,q);
  arma::colvec xd = x_vec - x0;
  
  for (int j = 1; j < (p + 1); j++) {
    x.col(j - 1) = pow(xd,j);
  }
  w = kernPtr[kernID].kern(xd/h);
  if (arma::accu(w != 0) == 1) {
    w = kernPtr[kernID].kern(xd/(3*h));
  }
  
  /* d = p;*/
  x_mat = arma::join_rows(arma::ones(n,1),x) % arma::repmat(arma::sqrt(w), 1, (p+1));
  y_mat = y % sqrt(w);
  beta  = arma::solve(x_mat, y_mat); 
  res   = y_mat - x_mat*beta;
  sig   = sqrt(std::inner_product(res.begin(), res.end(), res.begin(), 0.0)/(n));
  b     = beta.rows(1,p);
  for (int j = 1; j <= q; j++) {
    a(j-1) = beta(0) + R::qnorm5(tau(j - 1), 0.0, sig, 1, 0);
  }
  c     = 2*tau - 1;
  delta = 1;
  it    = 0;
  tn    = tol/n;
  e0    = -tn/log(tn);
  eps   = (e0-tn)/(1+log(e0));
  sumw  = arma::as_scalar(arma::accu(w));
  p5    = trans(w) * x;
  
  arma::mat ar, A, p2, p3, p4, tp, denor, numer, aa;
  arma::colvec a0, b0;
  while(delta > tol && it < maxit){
    it++;
    b0 = b;
    a0 = a;
    ar = arma::abs(arma::repmat(y-x*b0,1,q) - arma::repmat(trans(a0),n,1));
    A  = 1 / (ar + eps);
    p2 = trans(y % w) * A;
    p3 = trans(x % arma::repmat(w,1,p)) * A;
    p4 = trans(w) * A;
    denor = arma::zeros(p,p);
    numer = arma::zeros(p,1);
    for(int k = 0; k < q; k++){
      tp          = arma::repmat(A.col(k) % w, 1, p) % x;
      p1.col(k)   = trans(tp) * y;
      p6.slice(k) = trans(tp) * x;
      denor      += p6.slice(k) - p3.col(k) * trans(p3.col(k)) / p4(k);
      numer      += p1.col(k) - p2(k) * p3.col(k) / p4(k) -
        sumw * c(k) * p3.col(k) / p4(k) + c(k) * trans(p5);
    }
    b = arma::solve(denor,numer);
    a = trans( (trans((y - x * b) % w) * A + sumw * trans(c) ) / p4 );
    delta = arma::as_scalar(arma::abs( arma::mean(a-a0,0)) + arma::sum(arma::abs(b-b0) ));
  }
  
  estimate result;
  result.beta0 = a;
  result.beta1 = b;
  
  return result;
}

struct Allest1 : public RcppParallel::Worker {
  
  const  RcppParallel::RVector<double> x_vec;
  const  RcppParallel::RVector<double> y;
  int    kernID;
  const  RcppParallel::RVector<double> tau;
  double h;
  int    p;
  int    maxit;
  double tol;
  
  arma::colvec convert(const RcppParallel::RVector<double> input_mat)
  {
    arma::colvec output_mat(input_mat.begin(), input_mat.length());
    return output_mat;
  }
  
  arma::colvec x_vec2 = convert(x_vec);
  arma::colvec y2     = convert(y);
  arma::colvec tau2   = convert(tau);
  
  RcppParallel::RMatrix<double> rmat1;
  
  Allest1(const Rcpp::NumericVector x_vec, const Rcpp::NumericVector y, int kernID, 
         const Rcpp::NumericVector tau, double h, int p, int maxit, double tol, Rcpp::NumericMatrix rmat1)
    : x_vec(x_vec), y(y), kernID(kernID), tau(tau), h(h), p(p), maxit(maxit), tol(tol), rmat1(rmat1) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      estimate result = cqrMM(x_vec2(i), x_vec2, y2, kernID, tau2, h, p, maxit, tol);
      rmat1(i,0) = arma::as_scalar(arma::mean(result.beta0));
      rmat1(i,1) = arma::as_scalar(y2(i) - rmat1(i,0));
    }
  }
};

struct Allest2 : public RcppParallel::Worker {

  const  RcppParallel::RVector<double> x_vec;
  const  RcppParallel::RVector<double> y;
  const  RcppParallel::RVector<double> y_hat;
  const  RcppParallel::RVector<double> u_hat;
  int    kernID;
  double h;
  int    n = x_vec.length();
  double epsilon = 2.220446e-16;
  
  arma::colvec convert(const RcppParallel::RVector<double> input_mat)
  {
    arma::colvec output_mat(input_mat.begin(), input_mat.length());
    return output_mat;
  }
  
  arma::colvec x_vec2 = convert(x_vec);
  arma::colvec y2     = convert(y);
  arma::colvec y_hat2 = convert(y_hat);
  arma::colvec u_hat2 = convert(u_hat);
  
  RcppParallel::RMatrix<double> rmat2;
  
  Allest2(const Rcpp::NumericVector x_vec, const Rcpp::NumericVector y,  
          const Rcpp::NumericVector y_hat, const Rcpp::NumericVector u_hat,
          int kernID, double h, Rcpp::NumericMatrix rmat2)
    : x_vec(x_vec), y(y), y_hat(y_hat), u_hat(u_hat), kernID(kernID), h(h), rmat2(rmat2) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      arma::colvec u(n), w(n);
      u = (x_vec2 - x_vec2(i))/h;
      w = kernPtr[kernID].kern(u);
      rmat2(i,0) = std::sqrt( arma::as_scalar((trans(w) * arma::pow(u_hat2,2)) / arma::sum(w)) ) + epsilon;
      rmat2(i,1) = arma::as_scalar((y2(i) - y_hat2(i) + epsilon) / rmat2(i,0));
    }
  }
};



//' @title Function to compute the fitted values and residuals
//' 
//' @description This function computes the fitted values and residuals 
//' in a local composite quantile regression.
//' 
//' @param x_vec A vector of covariates.
//' @param y a vector of dependent variable, the treatment outcome variable
//' in the case of regression discontinuity.
//' @param kernID The kernel id that includes
//' \enumerate{
//'   \item \code{kernID = 0}: triangular kernel.
//'   \item \code{kernID = 1}: biweight kernel.
//'   \item \code{kernID = 2}: Epanechnikov kernel. 
//'   \item \code{kernID = 3}: Gaussian kernel.
//'   \item \code{kernID = 4}: tricube kernel.
//'   \item \code{kernID = 5}: triweight kernel. 
//'   \item \code{kernID = 6}: uniform kernel. 
//' }
//' @param tau A vector of quantile positions. They are obtained by 
//' \code{tau = (1:q)/(q+1)}.
//' @param h A scalar bandwidth.
//' @param p The polynomial degree. Defaults to 1.
//' @param maxit Maximum number of iterations in the MM algorithm. Defaults to 20.
//' @param tol The convergence criterion. Defaults to 1.0e-3.
//' @param parallel Set it to 1 if using parallel computing. Default is 1.
//' @param grainsize The minimum chunk size for parallelization. Defaults to 1.
//' @export
//' @return \code{est_cqr} returns a list with the following components:
//' \item{y_hat}{The fitted value at each point of the input vector \code{x_vec}.}
//' \item{u_hat}{The residual vector.}
//' \item{sig_hat}{The estimated standard deviation at each point of the input 
//' vector \code{x_vec}.}
//' \item{e_hat}{The scaled residual vector. It is the residual vector divided 
//' by the \code{sig_hat} vector.}
//' @examples
//' # Use the Head Start data as an example.
//' data(headstart)
//' data_n = subset(headstart, headstart$poverty < 0)
//' q      = 5
//' tau    = (1:q) / (q + 1)
//' est_cqr(x_vec  = data_n$poverty,
//'         y      = data_n$mortality,
//'         kernID = 2,
//'         tau    = tau,
//'         h      = 4.0,
//'         p      = 1) 
// [[Rcpp::export]]
Rcpp::List est_cqr(const arma::colvec & x_vec, const arma::colvec & y, int kernID, 
                   const arma::colvec & tau, double h, int p = 1, 
                   int maxit = 20, double tol = 1e-3, int parallel = 0, int grainsize = 1) {
  
  int n = arma::size(y)(0);
  
  if (parallel == 0) {
    
    arma::colvec y_hat(n), u_hat(n), e_hat(n), sig_hat(n), u(n),w(n);
    double epsilon = 2.220446e-16;
    Rcpp::NumericVector beta0;
    Rcpp::List est;
    
    for(int i = 0; i < n; i++){
      est   = cqrMMcpp(x_vec(i), x_vec, y, kernID = kernID, tau, h, p, maxit = maxit, tol = tol);
      beta0 = est[0];
      y_hat(i) = Rcpp::mean(beta0);
      u_hat(i) = y(i) - y_hat(i) + epsilon;
    }
    
    for(int i = 0; i < n; i++){
      u = (x_vec - x_vec(i))/h;
      w = kernPtr[kernID].kern(u);
      sig_hat(i) = std::sqrt( arma::as_scalar((trans(w) * arma::pow(u_hat,2)) / arma::sum(w)) ) + epsilon;
      e_hat(i)   = (y(i) - y_hat(i) + epsilon) / sig_hat(i);
    }
    
    return Rcpp::List::create(Rcpp::Named("y_hat")   = y_hat,
                              Rcpp::Named("u_hat")   = u_hat,
                              Rcpp::Named("sig_hat") = sig_hat,
                              Rcpp::Named("e_hat")   = e_hat);
    
  } else {

    Rcpp::NumericMatrix rmat1(n, 2);
    Rcpp::NumericMatrix rmat2(n, 2);

    Allest1 est1(Rcpp::NumericVector(x_vec.begin(), x_vec.end()), 
                 Rcpp::NumericVector(y.begin(), y.end()), kernID, 
                 Rcpp::NumericVector(tau.begin(), tau.end()), 
                 h, p, maxit, tol, rmat1);

    RcppParallel::parallelFor(0, n, est1, grainsize);
    
    Rcpp::NumericVector y_hat = rmat1( Rcpp::_ , 0);
    Rcpp::NumericVector u_hat = rmat1( Rcpp::_ , 1);
    Allest2 est2(Rcpp::NumericVector(x_vec.begin(), x_vec.end()), 
                 Rcpp::NumericVector(y.begin(), y.end()), 
                 y_hat, u_hat, kernID, h, rmat2);
    
    RcppParallel::parallelFor(0, n, est2, grainsize);
    
    return Rcpp::List::create(Rcpp::Named("y_hat")   = rmat1( Rcpp::_ , 0),
                              Rcpp::Named("u_hat")   = rmat1( Rcpp::_ , 1),
                              Rcpp::Named("sig_hat") = rmat2( Rcpp::_ , 0),
                              Rcpp::Named("e_hat")   = rmat2( Rcpp::_ , 1));
  }
}
