#include <Rmath.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <roptim.h>
// [[Rcpp::depends(roptim)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
// [[Rcpp::depends(RcppProgress)]]

/* Library functions of the paper:
Estimation of C*_pm under Progressive 
 Multiple Interval Type-I Censoring: Bayesian and Classical 
 Approaches with Log-Logistic Distribution
 */

#include<iostream>
#include<algorithm>

using namespace Rcpp;
using namespace std;
using namespace Numer;
using namespace arma;
using namespace roptim;

// gamse codes for gam progress bar inside gam loop
#include <progress.hpp>
#include <progress_bar.hpp>


// LogLogistic CDF
// [[Rcpp::export]]
double cdf_cpp(double x, NumericVector par) {
  double bet = par[0];
  double gam = par[1];
  return pow(x,gam)/(pow(x,gam)+pow(bet,gam));
}
// LogLogistic PDF
// [[Rcpp::export]]
double pdf_cpp(double x, NumericVector par) {
  double bet = par[0];
  double gam = par[1];
  return gam*pow(bet,gam)*pow(x,gam-1)/pow(pow(x,gam)+pow(bet,gam),2.0);
}
// LogLogistic Survival
// [[Rcpp::export]]
double sur_cpp(double x, NumericVector par) {
  return 1.0-cdf_cpp(x, par);
}
// LogLogistic Quantile
// [[Rcpp::export]]
double qnt_cpp(double x, NumericVector par) {
  double bet = par[0];
  double gam = par[1];
  return bet*pow(x/(1.0-x), 1.0/gam);
}

// The moments of LogLogistic Dist
// [[Rcpp::export]]
List moments_cpp(NumericVector par) {
  
  double bet = par[0];
  double gam = par[1];
  
  // ---- Beta terms used repeatedly ----
  double B1 = R::beta(1.0 + 1.0/gam, 1.0 - 1.0/gam);
  double B2 = R::beta(1.0 + 2.0/gam, 1.0 - 2.0/gam);
  double B3 = R::beta(1.0 + 3.0/gam, 1.0 - 3.0/gam);
  
  // ---- Mean μ ----
  double mu = bet*B1;
  
  // ---- Variance σ² ----
  double var = bet*bet*(B2-B1*B1);
  
  // ---- Third central moment μ₃ ----
  double mu3 = bet*bet*bet*(B3-3.0*B1*B2+2.0*std::pow(B1, 3.0));
  
  // Return results
  return List::create(
    _["mu"]    = mu,
    _["var"]   = var,
    _["mu3"]   = mu3
  );
}

// definition of C*pm function
// [[Rcpp::export]]
double Cpm_star_fun_cpp(NumericVector par,
                        NumericVector Cpm_par) {
  
  double bet = par[0];
  double gam = par[1];
  
  double l = Cpm_par[0];
  double u = Cpm_par[1];
  double t = Cpm_par[2];
  double lambda = Cpm_par[3];
  
  if (bet <= 0 || gam <= 0 || u <= l) {
    return NA_REAL;
  }
  
  // Quartiles and median
  double Q1 = qnt_cpp(0.25, par);
  double Q2 = qnt_cpp(0.50, par);
  double Q3 = qnt_cpp(0.75, par);
  
  // Equation (1.2)
  double Omega = (Q3 - Q1) / 2.0;
  
  double z = Q2 - t;
  
  double Lambda;
  
  // Numerical stability when lambda ≈ 0
  if (std::abs(lambda) < 1e-10) {
    Lambda = z * z;
  } else {
    Lambda = (2.0 * std::exp(lambda * z)
                - lambda * z
                - 1.0) / (lambda * lambda);
  }
  
  // Equation (1.1)
  double Cpm_star =
    (u - l) /
      (6.0 * std::sqrt(Omega * Omega + Lambda));
  
  return Cpm_star;
}


//------------------------------------------------
// Numerical gradient of C*_pm
//------------------------------------------------

// [[Rcpp::export]]
NumericVector Cpm_star_grad_cpp(NumericVector par,
                                NumericVector Cpm_par,
                                double eps = 1e-6) {
  
  NumericVector grad(2);
  
  //------------------------------------------------
  // derivative wrt beta
  //------------------------------------------------
  NumericVector par_plus  = clone(par);
  NumericVector par_minus = clone(par);
  
  par_plus[0]  += eps;
  par_minus[0] -= eps;
  
  grad[0] =
    (Cpm_star_fun_cpp(par_plus, Cpm_par) -
    Cpm_star_fun_cpp(par_minus, Cpm_par))/ (2.0 * eps);
  
  //------------------------------------------------
  // derivative wrt gamma
  //------------------------------------------------
  par_plus  = clone(par);
  par_minus = clone(par);
  
  par_plus[1]  += eps;
  par_minus[1] -= eps;
  
  grad[1] =
    (Cpm_star_fun_cpp(par_plus, Cpm_par) -
    Cpm_star_fun_cpp(par_minus, Cpm_par)) / (2.0 * eps);
  
  grad.names() = Rcpp::CharacterVector::create("beta", "gamma");
  
  return grad;
}


//------------------------------------------------
// Numerical Hessian matrix of C*_pm
//------------------------------------------------

// [[Rcpp::export]]
NumericMatrix Cpm_star_hessian_cpp(NumericVector par,
                                   NumericVector Cpm_par,
                                   double eps = 1e-5) {
  
  NumericMatrix H(2,2);

  //------------------------------------------------
  // Second derivative wrt beta
  //------------------------------------------------
  NumericVector p_plus  = clone(par);
  NumericVector p_minus = clone(par);
  
  p_plus[0]  += eps;
  p_minus[0] -= eps;
  
  double f0 = Cpm_star_fun_cpp(par,    Cpm_par);
  double fp = Cpm_star_fun_cpp(p_plus, Cpm_par);
  double fm = Cpm_star_fun_cpp(p_minus,Cpm_par);
  
  H(0,0) = (fp - 2.0*f0 + fm) / (eps*eps);
  
  //------------------------------------------------
  // Second derivative wrt gamma
  //------------------------------------------------
  p_plus  = clone(par);
  p_minus = clone(par);
  
  p_plus[1]  += eps;
  p_minus[1] -= eps;
  
  fp = Cpm_star_fun_cpp(p_plus, Cpm_par);
  fm = Cpm_star_fun_cpp(p_minus,Cpm_par);
  
  H(1,1) = (fp - 2.0*f0 + fm) / (eps*eps);
  
  //------------------------------------------------
  // Mixed derivative
  //------------------------------------------------
  NumericVector pp = clone(par);
  NumericVector pm = clone(par);
  NumericVector mp = clone(par);
  NumericVector mm = clone(par);
  
  pp[0] += eps; pp[1] += eps;
  pm[0] += eps; pm[1] -= eps;
  mp[0] -= eps; mp[1] += eps;
  mm[0] -= eps; mm[1] -= eps;
  
  double fpp = Cpm_star_fun_cpp(pp, Cpm_par);
  double fpm = Cpm_star_fun_cpp(pm, Cpm_par);
  double fmp = Cpm_star_fun_cpp(mp, Cpm_par);
  double fmm = Cpm_star_fun_cpp(mm, Cpm_par);
  
  double mixed =
    (fpp - fpm - fmp + fmm) /
      (4.0 * eps * eps);
  
  H(0,1) = mixed;
  H(1,0) = mixed;
  
  rownames(H) =
    CharacterVector::create("beta","gamma");
  
  colnames(H) =
    CharacterVector::create("beta","gamma");
  
  return H;
}


//------------------------------------------------
// Generated data from PMITICS 
//------------------------------------------------

// [[Rcpp::export]]
List simulate_PMITICS_cpp(NumericVector par,
                          int n,
                          int m,
                          double censor_prop,
                          double l_prop,
                          NumericVector Rprog) {
  
  if(Rprog.size() != 2 * m) {
    stop("Rprog must have length 2*m.");
  }
  
  // Step 1: Generate lifetimes
  NumericVector lifetimes(n);
  NumericVector u = runif(n);
  
  for(int i = 0; i < n; i++) {
    lifetimes[i] = qnt_cpp(u[i], par);
  }
  
  std::sort(lifetimes.begin(), lifetimes.end());
  
  // Step 2: Choose terminal time T_{2m}
  int idx = std::floor((1.0 - censor_prop) * n);
  
  if(idx < 0) idx = 0;
  if(idx >= n) idx = n - 1;
  
  double T2m = lifetimes[idx];
  
  // Step 3: Construct inspection times
  NumericVector T(2 * m + 1);
  
  T[0] = 0.0;
  
  double total_odd_length = l_prop * T2m;
  double odd_length  = total_odd_length / m;
  double even_length = (T2m - total_odd_length) / m;
  
  for(int j = 1; j <= m; j++) {
    T[2*j - 1] = T[2*j - 2] + odd_length;
    T[2*j]     = T[2*j - 1] + even_length;
  }
  
  // Step 4: Storage
  std::vector<int> alive(n, 1);
  
  std::vector<double> X;        // exact failures only
  std::vector<double> L;        // lower limits of interval failures
  std::vector<double> R;        // upper limits of interval failures
  
  IntegerVector r(m);           // exact counts
  IntegerVector d(m);           // interval counts
  IntegerVector Ractual(2 * m); // actual progressive removals
  
  // Step 5: PMITICS algorithm
  for(int j = 1; j <= m; j++) {
    
    // Exact failures in [T_{2j-2}, T_{2j-1}]
    double a1 = T[2*j - 2];
    double b1 = T[2*j - 1];
    
    for(int i = 0; i < n; i++) {
      if(alive[i] == 1 &&
         lifetimes[i] >= a1 &&
         lifetimes[i] <= b1) {
        
        X.push_back(lifetimes[i]);
        alive[i] = 0;
        r[j - 1]++;
      }
    }
    
    // Progressive removal at T_{2j-1}
    int k1 = 2*j - 2;
    int rem1 = Rprog[k1];
    
    std::vector<int> surv1;
    
    for(int i = 0; i < n; i++) {
      if(alive[i] == 1 && lifetimes[i] > b1) {
        surv1.push_back(i);
      }
    }
    
    if(rem1 > (int)surv1.size()) {
      rem1 = surv1.size();
    }
    
    if(rem1 > 0) {
      IntegerVector id1 = sample(surv1.size(), rem1, false) - 1;
      
      for(int h = 0; h < rem1; h++) {
        int id = surv1[id1[h]];
        alive[id] = 0;
        Ractual[k1]++;
      }
    }
    
    // Interval failures in (T_{2j-1}, T_{2j}]
    double a2 = T[2*j - 1];
    double b2 = T[2*j];
    
    for(int i = 0; i < n; i++) {
      if(alive[i] == 1 &&
         lifetimes[i] > a2 &&
         lifetimes[i] <= b2) {
        
        L.push_back(a2);
        R.push_back(b2);
        
        alive[i] = 0;
        d[j - 1]++;
      }
    }
    
    // Progressive removal at T_{2j}
    int k2 = 2*j - 1;
    int rem2 = Rprog[k2];
    
    std::vector<int> surv2;
    
    for(int i = 0; i < n; i++) {
      if(alive[i] == 1 && lifetimes[i] > b2) {
        surv2.push_back(i);
      }
    }
    
    if(rem2 > (int)surv2.size()) {
      rem2 = surv2.size();
    }
    
    if(rem2 > 0) {
      IntegerVector id2 = sample(surv2.size(), rem2, false) - 1;
      
      for(int h = 0; h < rem2; h++) {
        int id = surv2[id2[h]];
        alive[id] = 0;
        Ractual[k2]++;
      }
    }
  }
  
  // Remaining units beyond T_{2m}
  int final_right_censored = 0;
  
  for(int i = 0; i < n; i++) {
    if(alive[i] == 1 && lifetimes[i] > T2m) {
      final_right_censored++;
    }
  }
  
  std::sort(X.begin(), X.end());
  
  return List::create(
    _["X"] = wrap(X),
    _["L"] = wrap(L),
    _["R"] = wrap(R),
    _["r"] = r,
    _["d"] = d,
    _["T"] = T,
    _["Rprog"] = Ractual,
    _["final_right_censored"] = final_right_censored,
    _["T2m"] = T2m,
    _["lifetimes"] = lifetimes
  );
}

//------------------------------------------------
// Log-likelihood for PMITICS data from Algorithm 1
//------------------------------------------------

// [[Rcpp::export]]
double loglik_cpp(NumericVector par,
                  NumericVector X,
                  NumericVector L,
                  NumericVector R,
                  NumericVector T,
                  NumericVector Rprog,
                  int final_right_censored = 0) {
  
  double bet = par[0];
  double gam = par[1];
  
  if(bet <= 0 || gam <= 0) return R_NegInf;
  
  double logL = 0.0;
  double eps = 1e-12;
  
  // Exact failures: product f(x_i)
  for(int i = 0; i < X.size(); i++) {
    
    double dens = pdf_cpp(X[i], par);
    
    if(!R_finite(dens) || dens < eps) dens = eps;
    
    logL += std::log(dens);
  }
  
  // Interval failures: product [F(R_i) - F(L_i)]
  for(int i = 0; i < L.size(); i++) {
    
    double diff = cdf_cpp(R[i], par) - cdf_cpp(L[i], par);
    
    if(!R_finite(diff) || diff < eps) diff = eps;
    
    logL += std::log(diff);
  }
  
  // Progressive removals: product [1 - F(T_k)]^{R_k}
  if(Rprog.size() != T.size() - 1) {
    stop("Rprog must have length T.size() - 1.");
  }
  
  for(int k = 1; k < T.size(); k++) {
    
    double surv = 1.0 - cdf_cpp(T[k], par);
    
    if(!R_finite(surv) || surv < eps) surv = eps;
    
    logL += Rprog[k - 1] * std::log(surv);
  }
  
  // Remaining survivors beyond T_{2m}
  if(final_right_censored > 0) {
    
    double surv_T2m = 1.0 - cdf_cpp(T[T.size() - 1], par);
    
    if(!R_finite(surv_T2m) || surv_T2m < eps) {
      surv_T2m = eps;
    }
    
    logL += final_right_censored * std::log(surv_T2m);
  }
  
  return logL;
}

//------------------------------------------------
// Likelihood for PMITICS data from Algorithm 1
//------------------------------------------------

// [[Rcpp::export]]
double lik_cpp(NumericVector par,
               NumericVector X,
               NumericVector L,
               NumericVector R,
               NumericVector T,
               NumericVector Rprog,
               int final_right_censored = 0) {
  
  double val =
    loglik_cpp(par, X, L, R, T, Rprog, final_right_censored);
  
  if(!R_finite(val)) return 0.0;
  
  return std::exp(val);
}

//------------------------------------------------
// Log-product spacing function for PMITICS
//------------------------------------------------

// [[Rcpp::export]]
double logmps_cpp(NumericVector par,
                  NumericVector X,
                  NumericVector L,
                  NumericVector R,
                  NumericVector T,
                  NumericVector Rprog,
                  int final_right_censored = 0) {
  
  double bet = par[0];
  double gam = par[1];
  
  if(bet <= 0 || gam <= 0) return R_NegInf;
  
  double logS = 0.0;
  double eps = 1e-12;
  
  //------------------------------------------------
  // Exact spacings
  //------------------------------------------------
  NumericVector x = clone(X);
  
  std::sort(x.begin(), x.end());
  
  int r = x.size();
  
  if(r > 0) {
    
    NumericVector xx(r + 2);
    
    xx[0] = 0.0;
    
    for(int i = 0; i < r; i++) {
      xx[i + 1] = x[i];
    }
    
    xx[r + 1] = R_PosInf;
    
    for(int i = 1; i <= r + 1; i++) {
      
      double Fi   = cdf_cpp(xx[i], par);
      double Fim1 = cdf_cpp(xx[i - 1], par);
      
      // last spacing uses F(inf)=1
      if(std::isinf(xx[i])) Fi = 1.0;
      
      double spacing = Fi - Fim1;
      
      if(!R_finite(spacing) || spacing < eps) {
        spacing = eps;
      }
      
      logS += std::log(spacing);
    }
  }
  
  //------------------------------------------------
  // Interval spacings
  //------------------------------------------------
  for(int i = 0; i < L.size(); i++) {
    
    double diff =
      cdf_cpp(R[i], par) -
      cdf_cpp(L[i], par);
    
    if(!R_finite(diff) || diff < eps) {
      diff = eps;
    }
    
    logS += std::log(diff);
  }
  
  //------------------------------------------------
  // Progressive removals
  //------------------------------------------------
  if(Rprog.size() != T.size() - 1) {
    stop("Rprog must have length T.size()-1.");
  }
  
  for(int k = 1; k < T.size(); k++) {
    
    double surv =
      1.0 - cdf_cpp(T[k], par);
    
    if(!R_finite(surv) || surv < eps) {
      surv = eps;
    }
    
    logS +=
      Rprog[k - 1] * std::log(surv);
  }
  
  //------------------------------------------------
  // Final right censored observations
  //------------------------------------------------
  if(final_right_censored > 0) {
    
    double surv_T2m =
      1.0 - cdf_cpp(T[T.size() - 1], par);
    
    if(!R_finite(surv_T2m) || surv_T2m < eps) {
      surv_T2m = eps;
    }
    
    logS +=
      final_right_censored *
      std::log(surv_T2m);
  }
  
  return logS;
}

//------------------------------------------------
// Product spacing function for PMITICS
//------------------------------------------------

// [[Rcpp::export]]
double mps_cpp(NumericVector par,
               NumericVector X,
               NumericVector L,
               NumericVector R,
               NumericVector T,
               NumericVector Rprog,
               int final_right_censored = 0) {
  
  double val =
    logmps_cpp(
      par,
      X,
      L,
      R,
      T,
      Rprog,
      final_right_censored
    );
  
  if(!R_finite(val)) return 0.0;
  
  return std::exp(val);
}

// Returns -log-likelihood or -log-MPS function for PMITICS
class LogObjective : public Functor {
private:
  arma::vec X;
  arma::vec L;
  arma::vec R;
  arma::vec T;
  arma::vec Rprog;
  int final_right_censored;
  std::string type;
  
public:
  LogObjective(arma::vec XX_,
               arma::vec LL_,
               arma::vec RR_,
               arma::vec TT_,
               arma::vec Rprog_,
               int final_right_censored_,
               std::string Type_)
    : X(XX_),
      L(LL_),
      R(RR_),
      T(TT_),
      Rprog(Rprog_),
      final_right_censored(final_right_censored_),
      type(Type_) {}
  
  double operator()(const arma::vec &y) override {
    
    NumericVector par =
      Rcpp::as<NumericVector>(Rcpp::wrap(y));
    
    if(type == "LK") {
      return -loglik_cpp(
          par,
          Rcpp::as<NumericVector>(Rcpp::wrap(X)),
          Rcpp::as<NumericVector>(Rcpp::wrap(L)),
          Rcpp::as<NumericVector>(Rcpp::wrap(R)),
          Rcpp::as<NumericVector>(Rcpp::wrap(T)),
          Rcpp::as<NumericVector>(Rcpp::wrap(Rprog)),
          final_right_censored
      );
    }
    
    if(type == "PS") {
      return -logmps_cpp(
          par,
          Rcpp::as<NumericVector>(Rcpp::wrap(X)),
          Rcpp::as<NumericVector>(Rcpp::wrap(L)),
          Rcpp::as<NumericVector>(Rcpp::wrap(R)),
          Rcpp::as<NumericVector>(Rcpp::wrap(T)),
          Rcpp::as<NumericVector>(Rcpp::wrap(Rprog)),
          final_right_censored
      );
    }
    
    return R_PosInf;
  }
};

   
// Gamma prior distribution
// h = (a1, b1, a2, b2)
   
double prior_cpp(Rcpp::NumericVector par,
                 Rcpp::NumericVector h) {
     
     double bet = par[0];
     double gam = par[1];
     
     double a1 = h[0];
     double b1 = h[1];
     
     double a2 = h[2];
     double b2 = h[3];
     
     // parameter constraints
     if(bet <= 0 || gam <= 0)
       return 0.0;
     
     // independent gamma priors
     double prior_bet =
       std::pow(bet, a1 - 1.0) *
       std::exp(-b1 * bet);
     
     double prior_gam =
       std::pow(gam, a2 - 1.0) *
       std::exp(-b2 * gam);
     
     return prior_bet * prior_gam;
}
   
   
// Function to compute HPD interval
// -------------------------------------------------------
// HPD interval function
// -------------------------------------------------------
NumericVector HPD_cpp(NumericVector x, double prob = 0.95) {
  
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  
  int n = y.size();
  int m = std::floor(prob * n);
  
  if (m < 1) {
    stop("Not enough MCMC samples for HPD interval.");
  }
  
  double min_width = R_PosInf;
  int best = 0;
  
  for (int i = 0; i < n - m; i++) {
    double width = y[i + m] - y[i];
    if (width < min_width) {
      min_width = width;
      best = i;
    }
  }
  
  return NumericVector::create(y[best], y[best + m]);
}

arma::mat num_hess_cpp(const arma::vec& par,
                       const arma::vec& X,
                       const arma::vec& L,
                       const arma::vec& R,
                       const arma::vec& T,
                       const arma::vec& Rprog,
                       int final_right_censored,
                       const arma::vec& lower,
                       const arma::vec& upper,
                       std::string type,
                       double h = 1e-5) {
  
  int p = par.n_elem;
  arma::mat H(p, p, arma::fill::zeros);
  
  // create objective inside the function
  LogObjective obj(X, L, R, T, Rprog,
                   final_right_censored, type);
  
  auto neg_obj = [&](arma::vec z) {
    
    if(arma::any(z <= lower) || arma::any(z >= upper)) {
      return 1e20;
    }
    
    return obj(z);
  };
  
  double f0 = neg_obj(par);
  
  for(int i = 0; i < p; i++) {
    
    arma::vec ei(p, arma::fill::zeros);
    ei(i) = 1.0;
    
    double hi =
      std::min(h, 0.25 * std::max(1e-8, par(i) - lower(i)));
    
    arma::vec pp = par + hi * ei;
    arma::vec pm = par - hi * ei;
    
    H(i,i) =
      (neg_obj(pp) - 2.0 * f0 + neg_obj(pm)) /
        (hi * hi);
    
    for(int j = i + 1; j < p; j++) {
      
      arma::vec ej(p, arma::fill::zeros);
      ej(j) = 1.0;
      
      double hj =
        std::min(h, 0.25 * std::max(1e-8, par(j) - lower(j)));
      
      arma::vec ppp = par + hi * ei + hj * ej;
      arma::vec ppm = par + hi * ei - hj * ej;
      arma::vec pmp = par - hi * ei + hj * ej;
      arma::vec pmm = par - hi * ei - hj * ej;
      
      H(i,j) =
        (neg_obj(ppp) - neg_obj(ppm) -
        neg_obj(pmp) + neg_obj(pmm)) /
          (4.0 * hi * hj);
      
      H(j,i) = H(i,j);
    }
  }
  
  return H;
}


// ---------- 1D golden-section minimization ----------
double optimize_1d_cpp(std::function<double(double)> f,
                       double lower,
                       double upper,
                       double tol = 1e-8,
                       int max_iter = 1000) {
  
  const double gr = (std::sqrt(5.0) - 1.0) / 2.0;
  
  double a = lower;
  double b = upper;
  
  double c = b - gr * (b - a);
  double d = a + gr * (b - a);
  
  double fc = f(c);
  double fd = f(d);
  
  for (int i = 0; i < max_iter; i++) {
    if (std::abs(b - a) < tol) break;
    
    if (fc < fd) {
      b = d;
      d = c;
      fd = fc;
      c = b - gr * (b - a);
      fc = f(c);
    } else {
      a = c;
      c = d;
      fc = fd;
      d = a + gr * (b - a);
      fd = f(d);
    }
  }
  
  return 0.5 * (a + b);
}

//Computing the profile maximum likelihood function

// [[Rcpp::export]]
Rcpp::List Estim_profile_cpp(arma::vec par,
                             arma::vec par_init,
                             arma::vec X,
                             arma::vec L,
                             arma::vec R,
                             arma::vec T,
                             arma::vec Rprog,
                             int final_right_censored,
                             arma::vec Cpm_par,
                             arma::vec lower_lim,
                             arma::vec upper_lim,
                             std::string type = "LK") {
  
  bool Err = false;
  
  arma::vec Est(3, arma::fill::value(NA_REAL));
  arma::vec ESE(3, arma::fill::value(NA_REAL));
  arma::vec MSE(3, arma::fill::value(NA_REAL));
  arma::mat ACI(3, 4, arma::fill::value(NA_REAL));
  
  double opt_betha = par_init(0);
  double opt_gamta = par_init(1);
  
  try {
    
    LogObjective obj(X, L, R, T, Rprog,
                     final_right_censored, type);
    
    opt_gamta = optimize_1d_cpp(
      [&](double gamta) {
        arma::vec p = {opt_betha, gamta};
        return obj(p);
      },
      lower_lim(1),
      upper_lim(1)
    );
    
    opt_betha = optimize_1d_cpp(
      [&](double betha) {
        arma::vec p = {betha, opt_gamta};
        return obj(p);
      },
      lower_lim(0),
      upper_lim(0)
    );
    
  } catch (...) {
    Err = true;
  }
  
  if(!Err) {
    
    arma::vec par_hat = {opt_betha, opt_gamta};
    
    arma::mat H = num_hess_cpp(
      par_hat,
      X, L, R, T, Rprog,
      final_right_censored,
      lower_lim,
      upper_lim,
      type
    );
    
    // force symmetry
    H = 0.5 * (H + H.t());
    
    arma::mat VarCov;
    bool inv_ok = arma::inv_sympd(VarCov, H);
    
    if(!inv_ok) {
      inv_ok = arma::inv(VarCov, H);
    }
    
    if(!inv_ok || !VarCov.is_finite()) {
      Err = true;
    }
    
    if(!Err && (VarCov(0,0) <= 0.0 || VarCov(1,1) <= 0.0)) {
      Err = true;
    }
    
    if(!Err) {
      
      double Cpm_hat =
        Cpm_star_fun_cpp(wrap(par_hat), wrap(Cpm_par));
      
      arma::vec Cpm_grad =
        Cpm_star_grad_cpp(wrap(par_hat), wrap(Cpm_par));
      
      double var_Cpm =
        arma::as_scalar(Cpm_grad.t() * VarCov * Cpm_grad);
      
      if(!std::isfinite(var_Cpm) || var_Cpm <= 0.0) {
        Err = true;
      }
      
      if(!Err) {
        
        Est(0) = par_hat(0);
        Est(1) = par_hat(1);
        Est(2) = Cpm_hat;
        
        ESE(0) = std::sqrt(VarCov(0,0));
        ESE(1) = std::sqrt(VarCov(1,1));
        ESE(2) = std::sqrt(var_Cpm);
        
        double true_Cpm =
          Cpm_star_fun_cpp(wrap(par.subvec(0,1)), wrap(Cpm_par));
        
        arma::vec true_par = {par(0), par(1), true_Cpm};
        
        MSE = arma::square(Est - true_par);
        
        for(int j = 0; j < 3; j++) {
          ACI(j,0) = std::max(0.0, Est(j) - 1.96 * ESE(j));
          ACI(j,1) = Est(j) + 1.96 * ESE(j);
          ACI(j,2) = ACI(j,1) - ACI(j,0);
          ACI(j,3) =
            (ACI(j,0) <= true_par(j) &&
            true_par(j) <= ACI(j,1)) ? 1.0 : 0.0;
        }
      }
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("Err") = Err,
    Rcpp::Named("Est") = Est,
    Rcpp::Named("ESE") = ESE,
    Rcpp::Named("MSE") = MSE,
    Rcpp::Named("ACI") = ACI
  );
}

// ---------------------------------------------------
// Percentile bootstrap interval
// ---------------------------------------------------

// [[Rcpp::export]]
Rcpp::NumericVector make_boot_p_cpp(arma::vec x,
                                    double true_val,
                                    double alpha) {
  
  x = x.elem(find_finite(x));
  
  NumericVector xr = wrap(x);
  
  Function quantile("quantile");
  
  NumericVector probs =
    NumericVector::create(alpha / 2.0,
                          1.0 - alpha / 2.0);
  
  NumericVector q = quantile(
    xr,
    Named("probs") = probs,
    Named("na.rm") = true,
    Named("names") = false
  );
  
  double lower = std::max(0.0, q[0]);
  double upper = q[1];
  
  return NumericVector::create(
    lower,
    upper,
    upper - lower,
    (lower < true_val && true_val < upper) ? 1.0 : 0.0
  );
}

// ---------------------------------------------------
// Bootstrap-t interval
// ---------------------------------------------------

// [[Rcpp::export]]
Rcpp::NumericVector make_boot_t_cpp(arma::vec tstar,
                                    double est,
                                    double se0,
                                    double true_val,
                                    double alpha) {
  
  tstar = tstar.elem(find_finite(tstar));
  
  NumericVector tr = wrap(tstar);
  
  Function quantile("quantile");
  
  NumericVector probs =
    NumericVector::create(1.0 - alpha / 2.0,
                          alpha / 2.0);
  
  NumericVector q = quantile(
    tr,
    Named("probs") = probs,
    Named("na.rm") = true,
    Named("names") = false
  );
  
  double lower = est - q[0] * se0;
  double upper = est - q[1] * se0;
  
  lower = std::max(0.0, lower);
  
  return NumericVector::create(
    lower,
    upper,
    upper - lower,
    (lower < true_val && true_val < upper) ? 1.0 : 0.0
  );
}

// ---------------------------------------------------
// Log bootstrap-t interval
// ---------------------------------------------------

// [[Rcpp::export]]
Rcpp::NumericVector make_boot_t_log_cpp(arma::vec tlog,
                                        double est,
                                        double se0,
                                        double true_val,
                                        double alpha) {
  
  tlog = tlog.elem(find_finite(tlog));
  
  NumericVector tr = wrap(tlog);
  
  Function quantile("quantile");
  
  NumericVector probs =
    NumericVector::create(1.0 - alpha / 2.0,
                          alpha / 2.0);
  
  NumericVector q = quantile(
    tr,
    Named("probs") = probs,
    Named("na.rm") = true,
    Named("names") = false
  );
  
  double lower =
    est * std::exp(-q[0] * se0 / est);
  
  double upper =
    est * std::exp(-q[1] * se0 / est);
  
  lower = std::max(0.0, lower);
  
  return NumericVector::create(
    lower,
    upper,
    upper - lower,
    (lower < true_val && true_val < upper) ? 1.0 : 0.0
  );
}



// ---------------------------------------------------
// Main bootstrap function for PMITICS
// ---------------------------------------------------

// [[Rcpp::export]]
Rcpp::List Boot_fun_cpp(arma::vec par,
                        arma::vec init,
                        arma::vec se,
                        int n,
                        int m,
                        double censor_prop,
                        double l_prop,
                        arma::vec Rprog,
                        arma::vec lower,
                        arma::vec upper,
                        arma::vec Cpm_par,
                        int BOOT = 1000,
                        double alpha = 0.05,
                        std::string type = "LK") {
  
  arma::mat Boot_chain(BOOT, 3, arma::fill::value(NA_REAL));
  arma::mat Tstar_chain(BOOT, 3, arma::fill::value(NA_REAL));
  arma::mat Tlog_chain(BOOT, 3, arma::fill::value(NA_REAL));
  
  Function txtProgressBar("txtProgressBar");
  Function setTxtProgressBar("setTxtProgressBar");
  Function closePB("close");
  
  Rcpp::RObject pb = txtProgressBar(
    Named("min") = 0,
    Named("max") = BOOT,
    Named("style") = 3
  );
  
  for(int u = 0; u < BOOT; u++) {
    
    List dat_boot = simulate_PMITICS_cpp(
      wrap(init.subvec(0,1)),
      n,
      m,
      censor_prop,
      l_prop,
      wrap(Rprog)
    );
    
    arma::vec X_boot =
      as<arma::vec>(dat_boot["X"]);
    
    arma::vec L_boot =
      as<arma::vec>(dat_boot["L"]);
    
    arma::vec R_boot =
      as<arma::vec>(dat_boot["R"]);
    
    arma::vec T_boot =
      as<arma::vec>(dat_boot["T"]);
    
    arma::vec Rprog_boot =
      as<arma::vec>(dat_boot["Rprog"]);
    
    int final_right_censored_boot =
      as<int>(dat_boot["final_right_censored"]);
    
    List out;
    
    try {
      
      out = Estim_profile_cpp(
        par,
        init.subvec(0,1),
        X_boot,
        L_boot,
        R_boot,
        T_boot,
        Rprog_boot,
        final_right_censored_boot,
        Cpm_par,
        lower,
        upper,
        type
      );
      
    } catch (...) {
      
      setTxtProgressBar(pb, u + 1);
      continue;
    }
    
    bool Err = as<bool>(out["Err"]);
    
    if(Err) {
      setTxtProgressBar(pb, u + 1);
      continue;
    }
    
    arma::vec boot_est =
      as<arma::vec>(out["Est"]);
    
    arma::vec se_boot =
      as<arma::vec>(out["ESE"]);
    
    if(!boot_est.is_finite() ||
       !se_boot.is_finite() ||
       arma::any(se_boot <= 0.0) ||
       arma::any(boot_est <= 0.0)) {
      
      setTxtProgressBar(pb, u + 1);
      continue;
    }
    
    Boot_chain.row(u) = boot_est.t();
    
    Tstar_chain.row(u) =
      ((boot_est - init) / se_boot).t();
    
    Tlog_chain.row(u) =
      ((arma::log(boot_est) - arma::log(init)) /
        (se_boot / boot_est)).t();
    
    setTxtProgressBar(pb, u + 1);
  }
  
  closePB(pb);
  
  // ------------------------------------
  // Remove bad rows
  // ------------------------------------
  
  std::vector<int> keep_id;
  
  for(int i = 0; i < BOOT; i++) {
    
    bool keep = true;
    
    for(int j = 0; j < 3; j++) {
      
      if(!std::isfinite(Boot_chain(i,j))) keep = false;
      if(!std::isfinite(Tstar_chain(i,j))) keep = false;
      if(!std::isfinite(Tlog_chain(i,j))) keep = false;
      
      if(std::abs(Tstar_chain(i,j)) >= 10.0) keep = false;
      if(std::abs(Tlog_chain(i,j)) >= 10.0) keep = false;
    }
    
    if(keep) keep_id.push_back(i);
  }
  
  int Bgood = keep_id.size();
  
  arma::mat Boot_good(Bgood, 3);
  arma::mat Tstar_good(Bgood, 3);
  arma::mat Tlog_good(Bgood, 3);
  
  for(int i = 0; i < Bgood; i++) {
    Boot_good.row(i)   = Boot_chain.row(keep_id[i]);
    Tstar_good.row(i)  = Tstar_chain.row(keep_id[i]);
    Tlog_good.row(i)   = Tlog_chain.row(keep_id[i]);
  }
  
  // ------------------------------------
  // Percentile intervals
  // ------------------------------------
  
  NumericVector BOT_beta =
    make_boot_p_cpp(Boot_good.col(0), par[0], alpha);
  
  NumericVector BOT_gamma =
    make_boot_p_cpp(Boot_good.col(1), par[1], alpha);
  
  NumericVector BOT_Cpm =
    make_boot_p_cpp(Boot_good.col(2), par[2], alpha);
  
  // ------------------------------------
  // Bootstrap-t intervals
  // ------------------------------------
  
  NumericVector BTT_beta =
    make_boot_t_cpp(Tstar_good.col(0), init[0], se[0], par[0], alpha);
  
  NumericVector BTT_gamma =
    make_boot_t_cpp(Tstar_good.col(1), init[1], se[1], par[1], alpha);
  
  NumericVector BTT_Cpm =
    make_boot_t_cpp(Tstar_good.col(2), init[2], se[2], par[2], alpha);
  
  // ------------------------------------
  // Log bootstrap-t intervals
  // ------------------------------------
  
  NumericVector BLT_beta =
    make_boot_t_log_cpp(Tlog_good.col(0), init[0], se[0], par[0], alpha);
  
  NumericVector BLT_gamma =
    make_boot_t_log_cpp(Tlog_good.col(1), init[1], se[1], par[1], alpha);
  
  NumericVector BLT_Cpm =
    make_boot_t_log_cpp(Tlog_good.col(2), init[2], se[2], par[2], alpha);
  
  return List::create(
    
    Named("Boot_p.beta")  = BOT_beta,
    Named("Boot_p.gamma") = BOT_gamma,
    Named("Boot_p.Cpm")    = BOT_Cpm,
    
    Named("Boot_t.beta")  = BTT_beta,
    Named("Boot_t.gamma") = BTT_gamma,
    Named("Boot_t.Cpm")    = BTT_Cpm,
    
    Named("Boot_t_log.beta")  = BLT_beta,
    Named("Boot_t_log.gamma") = BLT_gamma,
    Named("Boot_t_log.Cpm")    = BLT_Cpm,
    
    Named("Boot_chain")  = Boot_good,
    Named("Tstar_chain") = Tstar_good,
    Named("Tlog_chain")  = Tlog_good,
    
    Named("n_success") = Bgood
  );
}

// Function to compute Bayesian MCMC sample for PMITICS data
// [[Rcpp::export]]
Rcpp::List MH_sample_cpp(NumericVector par,
                         std::string type,
                         NumericVector init,
                         NumericVector X,
                         NumericVector L,
                         NumericVector R,
                         NumericVector T,
                         NumericVector Rprog,
                         int final_right_censored,
                         NumericVector se,
                         NumericVector h,
                         NumericVector Cpm_par,
                         int MC_size,
                         NumericVector lower,
                         NumericVector upper) {
  
  NumericMatrix MH(MC_size, 3);
  
  MH(0,0) = init[0];
  MH(0,1) = init[1];
  MH(0,2) = Cpm_star_fun_cpp(init, Cpm_par);
  
  auto logpost_eval = [&](NumericVector pp) -> double {
    
    if(pp[0] <= 0.0 || pp[1] <= 0.0)
      return R_NegInf;
    
    double lp = std::log(prior_cpp(pp, h));
    
    double ll;
    
    if(type == "LK") {
      ll = loglik_cpp(pp, X, L, R, T, Rprog, final_right_censored);
    } else if(type == "PS") {
      ll = logmps_cpp(pp, X, L, R, T, Rprog, final_right_censored);
    } else {
      Rcpp::stop("Unknown type. Use 'LK' or 'PS'.");
    }
    
    if(!R_finite(lp) || !R_finite(ll))
      return R_NegInf;
    
    return lp + ll;
  };
  
  Function txtProgressBar("txtProgressBar");
  Function setTxtProgressBar("setTxtProgressBar");
  Function closePB("close");
  
  Rcpp::RObject pb = txtProgressBar(
    Named("min") = 0,
    Named("max") = MC_size,
    Named("style") = 3
  );
  
  NumericVector prop(2);
  NumericVector curr(2);
  
  curr[0] = init[0];
  curr[1] = init[1];
  
  int i = 1;
  int accepts = 0;
  
  while(i < MC_size) {
    
    prop[0] = std::exp(rnorm(1, std::log(curr[0]), se[0])[0]);
    prop[1] = std::exp(rnorm(1, std::log(curr[1]), se[1])[0]);
    
    if(NumericVector::is_na(prop[0]) ||
       NumericVector::is_na(prop[1])) {
      continue;
    }
    
    if(prop[0] < lower[0] || prop[0] > upper[0]) continue;
    if(prop[1] < lower[1] || prop[1] > upper[1]) continue;
    
    double logratio = logpost_eval(prop) - logpost_eval(curr);
    
    if(!R_finite(logratio)) continue;
    
    double del = std::min(1.0, std::exp(logratio));
    double dad = runif(1)[0];
    
    if(dad < del) {
      curr[0] = prop[0];
      curr[1] = prop[1];
      accepts++;
    }
    
    MH(i,0) = curr[0];
    MH(i,1) = curr[1];
    MH(i,2) = Cpm_star_fun_cpp(curr, Cpm_par);
    
    setTxtProgressBar(pb, i + 1);
    
    i++;
  }
  
  closePB(pb);
  
  colnames(MH) =
    CharacterVector::create("beta", "gamma", "Cpm");
  
  return Rcpp::List::create(
    Rcpp::Named("MH_sample") = MH,
    Rcpp::Named("Accpts") =
      (double)accepts / (MC_size - 1)
  );
}
// [[Rcpp::export]]
Rcpp::List MCMC_cpp(NumericVector par,
                    std::string type,
                    NumericVector init,
                    NumericVector X,
                    NumericVector L,
                    NumericVector R,
                    NumericVector T,
                    NumericVector Rprog,
                    int final_right_censored,
                    NumericVector se,
                    NumericVector h,
                    NumericVector Cpm_par,
                    int MC_size,
                    int MC_burn,
                    NumericVector c,
                    NumericVector lower,
                    NumericVector upper) {
  
  if(MC_burn >= MC_size) {
    stop("MC_burn must be smaller than MC_size.");
  }
  
  Rcpp::List MH_out = MH_sample_cpp(
    par, type, init, X, L, R, T, Rprog,
    final_right_censored, se, h, Cpm_par,
    MC_size, lower, upper
  );
  
  NumericMatrix MH_all = MH_out["MH_sample"];
  
  int n_keep = MC_size - MC_burn;
  int p = MH_all.ncol();
  
  NumericMatrix MH_sel(n_keep, p);
  
  for(int i = 0; i < n_keep; i++) {
    for(int j = 0; j < p; j++) {
      MH_sel(i,j) = MH_all(i + MC_burn, j);
    }
  }

  
  NumericVector est_SEL(p);
  
  for(int j = 0; j < p; j++) {
    double sum_j = 0.0;
    
    for(int i = 0; i < n_keep; i++) {
      sum_j += MH_sel(i,j);
    }
    
    est_SEL[j] = sum_j / n_keep;
  }
  
  NumericVector est_LIN1(p);
  NumericVector est_LIN2(p);
  
  for(int j = 0; j < p; j++) {
    
    double sum_lin1 = 0.0;
    double sum_lin2 = 0.0;
    
    for(int i = 0; i < n_keep; i++) {
      sum_lin1 += std::exp(-c[0] * MH_sel(i,j));
      sum_lin2 += std::exp(-c[1] * MH_sel(i,j));
    }
    
    est_LIN1[j] =
      -(1.0 / c[0]) * std::log(sum_lin1 / n_keep);
      
      est_LIN2[j] =
      -(1.0 / c[1]) * std::log(sum_lin2 / n_keep);
  }
  
  NumericMatrix HPD(3, 4);
  
  for(int j = 0; j < 3; j++) {
    
    NumericVector sample_j(n_keep);
    
    for(int i = 0; i < n_keep; i++) {
      sample_j[i] = MH_sel(i,j);
    }
    
    NumericVector hpd_j =
      HPD_cpp(sample_j, 0.95);
    
    HPD(j,0) = hpd_j[0];
    HPD(j,1) = hpd_j[1];
    HPD(j,2) = HPD(j,1) - HPD(j,0);
    
    if(j < par.size()) {
      HPD(j,3) =
        (HPD(j,0) <= par[j] &&
        par[j] <= HPD(j,1)) ? 1.0 : 0.0;
    } else {
      HPD(j,3) = NA_REAL;
    }
  }
  
  colnames(HPD) =
    CharacterVector::create("lower", "upper", "AL", "CP");
  
  colnames(MH_sel) =
    CharacterVector::create("beta", "gamma", "Cpm");
  
  est_SEL.names() =
    CharacterVector::create("beta", "gamma", "Cpm");
  
  est_LIN1.names() =
    CharacterVector::create("beta", "gamma", "Cpm");
  
  est_LIN2.names() =
    CharacterVector::create("beta", "gamma", "Cpm");
  
  return Rcpp::List::create(
    Rcpp::Named("est_SEL")    = est_SEL,
    Rcpp::Named("est_LIN1")   = est_LIN1,
    Rcpp::Named("est_LIN2")   = est_LIN2,
    Rcpp::Named("HPD")        = HPD,
    Rcpp::Named("MH_sample")  = MH_sel,
    Rcpp::Named("acceptance") = MH_out["Accpts"]
  );
}

// -------------------------------------------------------
// Hyperparameters based on PMITICS samples
// -------------------------------------------------------

class HyperObj_PMITICS : public Functor {
  
private:
  Rcpp::NumericVector X, L, R, T, Rprog;
  int final_right_censored;
  
public:
  
  HyperObj_PMITICS(Rcpp::NumericVector X_,
                   Rcpp::NumericVector L_,
                   Rcpp::NumericVector R_,
                   Rcpp::NumericVector T_,
                   Rcpp::NumericVector Rprog_,
                   int final_right_censored_)
    : X(X_),
      L(L_),
      R(R_),
      T(T_),
      Rprog(Rprog_),
      final_right_censored(final_right_censored_) {}
  
  double operator()(const arma::vec& x) override {
    
    if(x[0] <= 0.0 || x[1] <= 0.0)
      return 1e100;
    
    Rcpp::NumericVector par =
      Rcpp::NumericVector::create(x[0], x[1]);
    
    double ll = loglik_cpp(
      par,
      X,
      L,
      R,
      T,
      Rprog,
      final_right_censored
    );
    
    if(!R_finite(ll))
      return 1e100;
    
    return -ll;   // minimization
  }
};

// [[Rcpp::export]]
Rcpp::NumericVector HyperPara_PMITICS_cpp(Rcpp::NumericVector para,
                                          int n,
                                          int m,
                                          double censor_prop,
                                          double l_prop,
                                          Rcpp::NumericVector Rprog,
                                          Rcpp::NumericVector lower,
                                          Rcpp::NumericVector upper,
                                          int Hyper_It = 100,
                                          std::string method = "Nelder-Mead") {
  
  Rcpp::NumericMatrix MLE(Hyper_It, 2);
  
  int it = 0;
  
  while(it < Hyper_It) {
    
    Rcpp::List dat;
    
    try {
      dat = simulate_PMITICS_cpp(
       para, n, m, censor_prop, l_prop, Rprog
      );
    } catch(...) {
      continue;
    }
    
    Rcpp::NumericVector X_boot = dat["X"];
    Rcpp::NumericVector L_boot = dat["L"];
    Rcpp::NumericVector R_boot = dat["R"];
    Rcpp::NumericVector T_boot = dat["T"];
    Rcpp::NumericVector Rprog_boot = dat["Rprog"];
    
    int final_right_censored_boot =
      Rcpp::as<int>(dat["final_right_censored"]);
    
    HyperObj_PMITICS obj(
        X_boot,
        L_boot,
        R_boot,
        T_boot,
        Rprog_boot,
        final_right_censored_boot
    );
    
    arma::vec init(2);
    init[0] = para[0];
    init[1] = para[1];
    
    Roptim<HyperObj_PMITICS> opt(method);
    opt.control.trace = 0;
    
    if(method == "L-BFGS-B") {
      opt.set_lower(Rcpp::as<arma::vec>(lower));
      opt.set_upper(Rcpp::as<arma::vec>(upper));
    }
    
    try {
      opt.minimize(obj, init);
    } catch(...) {
      continue;
    }
    
    arma::vec est = opt.par();
    
    if(!std::isfinite(est[0]) || !std::isfinite(est[1]) ||
       est[0] <= lower[0] || est[0] >= upper[0] ||
       est[1] <= lower[1] || est[1] >= upper[1]) {
      continue;
    }
    
    MLE(it,0) = est[0];
    MLE(it,1) = est[1];
    
    it++;
  }
  
  if(it < 2) {
    return Rcpp::NumericVector::create(
      Rcpp::Named("a1") = NA_REAL,
      Rcpp::Named("b1") = NA_REAL,
      Rcpp::Named("a2") = NA_REAL,
      Rcpp::Named("b2") = NA_REAL
    );
  }
  
  Rcpp::NumericVector bet(it);
  Rcpp::NumericVector gam(it);
  
  for(int i = 0; i < it; i++) {
    bet[i] = MLE(i,0);
    gam[i] = MLE(i,1);
  }
  
  double mean_bet = Rcpp::mean(bet);
  double mean_gam = Rcpp::mean(gam);
  
  double se_bet = Rcpp::sd(bet);
  double se_gam = Rcpp::sd(gam);
  
  if(se_bet <= 0.0 || se_gam <= 0.0 ||
     !R_finite(se_bet) || !R_finite(se_gam)) {
     
     return Rcpp::NumericVector::create(
       Rcpp::Named("a1") = NA_REAL,
       Rcpp::Named("b1") = NA_REAL,
       Rcpp::Named("a2") = NA_REAL,
       Rcpp::Named("b2") = NA_REAL
     );
  }
  
  double a1 = std::pow(mean_bet, 2.0) / std::pow(se_bet, 2.0);
  double b1 = mean_bet / std::pow(se_bet, 2.0);
  
  double a2 = std::pow(mean_gam, 2.0) / std::pow(se_gam, 2.0);
  double b2 = mean_gam / std::pow(se_gam, 2.0);
  
  return Rcpp::NumericVector::create(
    Rcpp::Named("a1") = a1,
    Rcpp::Named("b1") = b1,
    Rcpp::Named("a2") = a2,
    Rcpp::Named("b2") = b2
  );
}

// logartihm of indep. gamma prior distribution

double logprior_cpp(Rcpp::NumericVector par,
                    Rcpp::NumericVector h) {
  
  double beta  = par[0];
  double gamma = par[1];
  
  double a1 = h[0];
  double b1 = h[1];
  double a2 = h[2];
  double b2 = h[3];
  
  if(beta <= 0.0 || gamma <= 0.0)
    return R_NegInf;
  
  return
  (a1 - 1.0) * std::log(beta)  - b1 * beta +
    (a2 - 1.0) * std::log(gamma) - b2 * gamma;
}

// Bayes using TK under LINEX

class TK_Obj_GEL : public Functor {
  
private:
  arma::vec X, L, R, T, Rprog, Cpm_par, h;
  int Final_right_censored, N, NoPara;
  double Q;
  std::string Type;
  
public:
  
  TK_Obj_GEL(arma::vec xx_,
             arma::vec ll_,
             arma::vec rr_,
             arma::vec tt_,
             arma::vec rprog_,
             int final_right_censored_,
             int nn_,
             arma::vec cpm_par_,
             arma::vec h_,
             double q_,
             std::string type_,
             int nopara_)
    : X(xx_), L(ll_), R(rr_), T(tt_), Rprog(rprog_),
      Cpm_par(cpm_par_), h(h_),
      Final_right_censored(final_right_censored_),
      N(nn_), NoPara(nopara_), Q(q_), Type(type_) {}
  
  double operator()(const arma::vec &y) override {
    
    if(y[0] <= 0.0 || y[1] <= 0.0)
      return 1e100;
    
    Rcpp::NumericVector par = Rcpp::wrap(y);
    
    Rcpp::NumericVector x   = Rcpp::wrap(X);
    Rcpp::NumericVector l   = Rcpp::wrap(L);
    Rcpp::NumericVector r   = Rcpp::wrap(R);
    Rcpp::NumericVector tt  = Rcpp::wrap(T);
    Rcpp::NumericVector rp  = Rcpp::wrap(Rprog);
    Rcpp::NumericVector hh  = Rcpp::wrap(h);
    Rcpp::NumericVector Cpm = Rcpp::wrap(Cpm_par);
    
    double log_obj;
    
    if(Type == "LK") {
      log_obj = loglik_cpp(par, x, l, r, tt, rp,
                           Final_right_censored);
    } else if(Type == "PS") {
      log_obj = logmps_cpp(par, x, l, r, tt, rp,
                           Final_right_censored);
    } else {
      return 1e100;
    }
    
    double lp = logprior_cpp(par, hh);
    
    if(!R_finite(log_obj) || !R_finite(lp))
      return 1e100;
    
    log_obj += lp;
    
    double gval = 1.0;
    
    if(NoPara == 1) gval = y[0];
    if(NoPara == 2) gval = y[1];
    if(NoPara == 3) gval = Cpm_star_fun_cpp(par, Cpm);
    
    if(!R_finite(gval) || gval <= 0.0)
      return 1e100;
    
    log_obj += -Q * std::log(gval);
    
    return -log_obj / N;
  }
};

// [[Rcpp::export]]
Rcpp::List TK_GEL_cpp(arma::vec True_para,
                      arma::vec Init_para,
                      int n,
                      arma::vec X,
                      arma::vec L,
                      arma::vec R,
                      arma::vec T,
                      arma::vec Rprog,
                      int final_right_censored,
                      arma::vec Cpm_par,
                      arma::vec h,
                      double q,
                      std::string type) {
  
  Rcpp::NumericVector est_GEL(3), bias_GEL(3), mse_GEL(3);
  
  if(q == 0.0) {
    Rcpp::stop("q must be nonzero for GEL.");
  }
  
  auto run_TK = [&](double qval, int nopara) {
    
    TK_Obj_GEL obj(
        X, L, R, T, Rprog,
        final_right_censored,
        n,
        Cpm_par,
        h,
        qval,
        type,
        nopara
    );
    
    Roptim<TK_Obj_GEL> opt("Nelder-Mead");
    opt.control.trace = 0;
    opt.set_hessian(true);
    
    arma::vec init = Init_para.subvec(0, 1);
    
    try {
      opt.minimize(obj, init);
    } catch(...) {
      return Rcpp::NumericVector::create(NA_REAL, NA_REAL);
    }
    
    double val = -opt.value();
    
    arma::mat H =
      Rcpp::as<arma::mat>(Rcpp::wrap(opt.hessian()));
    
    H = 0.5 * (H + H.t());
    
    arma::mat I;
    bool ok = arma::inv_sympd(I, H);
    
    if(!ok) ok = arma::inv(I, H);
    
    if(!ok || !I.is_finite())
      return Rcpp::NumericVector::create(NA_REAL, NA_REAL);
    
    double detI = arma::det(I);
    
    if(!R_finite(detI) || detI <= 0.0)
      return Rcpp::NumericVector::create(NA_REAL, NA_REAL);
    
    return Rcpp::NumericVector::create(val, detI);
  };
  
  Rcpp::NumericVector base = run_TK(0.0, 1);
  
  double obj0 = base[0];
  double det0 = base[1];
  
  if(!R_finite(obj0) || !R_finite(det0) || det0 <= 0.0) {
    return Rcpp::List::create(Rcpp::Named("Err") = true);
  }
  
  for(int j = 0; j < 3; j++) {
    
    Rcpp::NumericVector tmp = run_TK(q, j + 1);
    
    double objj = tmp[0];
    double detj = tmp[1];
    
    if(!R_finite(objj) || !R_finite(detj) || detj <= 0.0) {
      est_GEL[j] = NA_REAL;
    } else {
      
      double ratio =
        std::sqrt(detj / det0) *
        std::exp(n * (objj - obj0));
      
      est_GEL[j] =
        std::pow(ratio, -1.0 / q);
    }
    
    bias_GEL[j] = est_GEL[j] - True_para[j];
    mse_GEL[j]  = std::pow(bias_GEL[j], 2.0);
  }
  
  est_GEL.names() =
    Rcpp::CharacterVector::create("beta", "gamma", "Cpm");
  
  return Rcpp::List::create(
    Rcpp::Named("Err") = false,
    Rcpp::Named("est_GEL") = est_GEL,
    Rcpp::Named("bias_GEL") = bias_GEL,
    Rcpp::Named("mse_GEL") = mse_GEL
  );
}

// Bayes using TK under LINEX

class TK_Obj_LINEX : public Functor {
  
private:
  arma::vec X, L, R, T, Rprog, Cpm_par, h;
  int Final_right_censored, N, NoPara;
  double C;
  std::string Type;
  
public:
  
  TK_Obj_LINEX(arma::vec xx_,
               arma::vec ll_,
               arma::vec rr_,
               arma::vec tt_,
               arma::vec rprog_,
               int final_right_censored_,
               int nn_,
               arma::vec cpm_par_,
               arma::vec h_,
               double c_,
               std::string type_,
               int nopara_)
    : X(xx_), L(ll_), R(rr_), T(tt_), Rprog(rprog_),
      Cpm_par(cpm_par_), h(h_),
      Final_right_censored(final_right_censored_),
      N(nn_), NoPara(nopara_), C(c_), Type(type_) {}
  
  double operator()(const arma::vec &y) override {
    
    if(y[0] <= 0.0 || y[1] <= 0.0)
      return 1e100;
    
    Rcpp::NumericVector par = Rcpp::wrap(y);
    
    Rcpp::NumericVector x   = Rcpp::wrap(X);
    Rcpp::NumericVector l   = Rcpp::wrap(L);
    Rcpp::NumericVector r   = Rcpp::wrap(R);
    Rcpp::NumericVector tt  = Rcpp::wrap(T);
    Rcpp::NumericVector rp  = Rcpp::wrap(Rprog);
    Rcpp::NumericVector hh  = Rcpp::wrap(h);
    Rcpp::NumericVector Cpm = Rcpp::wrap(Cpm_par);
    
    double log_obj;
    
    if(Type == "LK") {
      log_obj = loglik_cpp(par, x, l, r, tt, rp,
                           Final_right_censored);
    } else if(Type == "PS") {
      log_obj = logmps_cpp(par, x, l, r, tt, rp,
                           Final_right_censored);
    } else {
      return 1e100;
    }
    
    double lp = logprior_cpp(par, hh);
    
    if(!R_finite(log_obj) || !R_finite(lp))
      return 1e100;
    
    log_obj += lp;
    
    double gval = 1.0;
    
    if(NoPara == 1) gval = y[0];
    if(NoPara == 2) gval = y[1];
    if(NoPara == 3) gval = Cpm_star_fun_cpp(par, Cpm);
    
    if(!R_finite(gval) || gval <= 0.0)
      return 1e100;
    
    log_obj += -C * gval;
    
    return -log_obj / N;
  }
};
// [[Rcpp::export]]
Rcpp::List TK_LINEX_cpp(arma::vec True_para,
                        arma::vec Init_para,
                        int n,
                        arma::vec X,
                        arma::vec L,
                        arma::vec R,
                        arma::vec T,
                        arma::vec Rprog,
                        int final_right_censored,
                        arma::vec Cpm_par,
                        arma::vec h,
                        double c,
                        std::string type) {
  
  Rcpp::NumericVector est_LINEX(3), bias_LINEX(3), mse_LINEX(3);
  
  if(c == 0.0) {
    Rcpp::stop("c must be nonzero for LINEX.");
  }
  
  auto run_TK = [&](double cval, int nopara) {
    
    TK_Obj_LINEX obj(
        X, L, R, T, Rprog,
        final_right_censored,
        n,
        Cpm_par,
        h,
        cval,
        type,
        nopara
    );
    
    Roptim<TK_Obj_LINEX> opt("Nelder-Mead");
    opt.control.trace = 0;
    opt.set_hessian(true);
    
    arma::vec init = Init_para.subvec(0, 1);
    
    try {
      opt.minimize(obj, init);
    } catch(...) {
      return Rcpp::NumericVector::create(NA_REAL, NA_REAL);
    }
    
    double val = -opt.value();
    
    arma::mat H =
      Rcpp::as<arma::mat>(Rcpp::wrap(opt.hessian()));
    
    H = 0.5 * (H + H.t());
    
    arma::mat I;
    bool ok = arma::inv_sympd(I, H);
    
    if(!ok) ok = arma::inv(I, H);
    
    if(!ok || !I.is_finite())
      return Rcpp::NumericVector::create(NA_REAL, NA_REAL);
    
    double detI = arma::det(I);
    
    if(!R_finite(detI) || detI <= 0.0)
      return Rcpp::NumericVector::create(NA_REAL, NA_REAL);
    
    return Rcpp::NumericVector::create(val, detI);
  };
  
  Rcpp::NumericVector base = run_TK(0.0, 1);
  
  double obj0 = base[0];
  double det0 = base[1];
  
  if(!R_finite(obj0) || !R_finite(det0) || det0 <= 0.0) {
    return Rcpp::List::create(Rcpp::Named("Err") = true);
  }
  
  for(int j = 0; j < 3; j++) {
    
    Rcpp::NumericVector tmp = run_TK(c, j + 1);
    
    double objj = tmp[0];
    double detj = tmp[1];
    
    if(!R_finite(objj) || !R_finite(detj) || detj <= 0.0) {
      est_LINEX[j] = NA_REAL;
    } else {
      
      double ratio =
        std::sqrt(detj / det0) *
        std::exp(n * (objj - obj0));
      
      est_LINEX[j] =
        -(1.0 / c) * std::log(ratio);
    }
    
    bias_LINEX[j] = est_LINEX[j] - True_para[j];
    mse_LINEX[j]  = std::pow(bias_LINEX[j], 2.0);
  }
  
  est_LINEX.names() =
    Rcpp::CharacterVector::create("beta", "gamma", "Cpm");
  
  return Rcpp::List::create(
    Rcpp::Named("Err") = false,
    Rcpp::Named("est_LINEX") = est_LINEX,
    Rcpp::Named("bias_LINEX") = bias_LINEX,
    Rcpp::Named("mse_LINEX") = mse_LINEX
  );
}

// Computing Bayes estimaties using Lindley's approximation


double Lindley_core_cpp(arma::vec th,
                        arma::vec X, 
                        arma::vec L, 
                        arma::vec R,
                        arma::vec T, 
                        arma::vec Rprog,
                        int final_right_censored,
                        arma::vec h,
                        arma::vec Cpm_par,
                        std::string type,
                        int NoPara,
                        std::string loss,
                        double q,
                        double eps = 1e-5) {

  auto loglik_fun = [&](arma::vec z) {
    NumericVector zz = wrap(z);
    if(type == "LK")
      return loglik_cpp(zz, wrap(X), wrap(L), wrap(R), wrap(T), wrap(Rprog), final_right_censored);
    else
      return logmps_cpp(zz, wrap(X), wrap(L), wrap(R), wrap(T), wrap(Rprog), final_right_censored);
  };

  auto basic_g = [&](arma::vec z) {
    NumericVector zz = wrap(z);
    NumericVector cp = wrap(Cpm_par);

    if(NoPara == 1) return z[0];
    if(NoPara == 2) return z[1];
    return Cpm_star_fun_cpp(zz, cp);
  };

  auto trans_g = [&](arma::vec z) {
    double g0 = basic_g(z);
    if(g0 <= 0.0 || !R_finite(g0)) return NA_REAL;

    if(loss == "GEL")   return std::pow(g0, -q);
    if(loss == "LINEX") return std::exp(-q * g0);

    return g0;
  };

  double bet = th[0];
  double gam = th[1];

  double r1 = (h[0] - 1.0) / bet - h[1];
  double r2 = (h[2] - 1.0) / gam - h[3];

  arma::mat H(2,2,arma::fill::zeros);
  double f0 = loglik_fun(th);

  for(int i=0;i<2;i++) {
    arma::vec ei(2,arma::fill::zeros);
    ei[i] = 1.0;

    H(i,i) =
      (loglik_fun(th + eps*ei) - 2.0*f0 + loglik_fun(th - eps*ei)) /
      (eps*eps);

    for(int j=i+1;j<2;j++) {
      arma::vec ej(2,arma::fill::zeros);
      ej[j] = 1.0;

      H(i,j) =
        (loglik_fun(th + eps*ei + eps*ej)
       - loglik_fun(th + eps*ei - eps*ej)
       - loglik_fun(th - eps*ei + eps*ej)
       + loglik_fun(th - eps*ei - eps*ej)) /
        (4.0*eps*eps);

      H(j,i) = H(i,j);
    }
  }

  H = 0.5 * (H + H.t());

  arma::mat Sigma;
  bool ok = arma::inv_sympd(Sigma, -H);
  if(!ok) ok = arma::inv(Sigma, -H);
  if(!ok || !Sigma.is_finite()) return NA_REAL;

  double s11 = Sigma(0,0);
  double s12 = Sigma(0,1);
  double s22 = Sigma(1,1);

  double g = trans_g(th);

  arma::vec g1(2);
  arma::mat g2(2,2,arma::fill::zeros);

  for(int i=0;i<2;i++) {
    arma::vec ei(2,arma::fill::zeros);
    ei[i] = 1.0;

    g1[i] = (trans_g(th + eps*ei) - trans_g(th - eps*ei)) / (2.0*eps);

    g2(i,i) =
      (trans_g(th + eps*ei) - 2.0*g + trans_g(th - eps*ei)) /
      (eps*eps);

    for(int j=i+1;j<2;j++) {
      arma::vec ej(2,arma::fill::zeros);
      ej[j] = 1.0;

      g2(i,j) =
        (trans_g(th + eps*ei + eps*ej)
       - trans_g(th + eps*ei - eps*ej)
       - trans_g(th - eps*ei + eps*ej)
       + trans_g(th - eps*ei - eps*ej)) /
        (4.0*eps*eps);

      g2(j,i) = g2(i,j);
    }
  }

  auto Hess_loglik = [&](arma::vec z) {
    arma::mat HH(2,2,arma::fill::zeros);
    double ff0 = loglik_fun(z);

    for(int i=0;i<2;i++) {
      arma::vec ei(2,arma::fill::zeros);
      ei[i] = 1.0;

      HH(i,i) =
        (loglik_fun(z + eps*ei) - 2.0*ff0 + loglik_fun(z - eps*ei)) /
        (eps*eps);
    }

    arma::vec e1(2,arma::fill::zeros), e2(2,arma::fill::zeros);
    e1[0] = 1.0; e2[1] = 1.0;

    HH(0,1) =
      (loglik_fun(z + eps*e1 + eps*e2)
     - loglik_fun(z + eps*e1 - eps*e2)
     - loglik_fun(z - eps*e1 + eps*e2)
     + loglik_fun(z - eps*e1 - eps*e2)) /
      (4.0*eps*eps);

    HH(1,0) = HH(0,1);
    return HH;
  };

  arma::vec e1(2,arma::fill::zeros), e2(2,arma::fill::zeros);
  e1[0] = 1.0; e2[1] = 1.0;

  arma::mat H_p1 = Hess_loglik(th + eps*e1);
  arma::mat H_m1 = Hess_loglik(th - eps*e1);
  arma::mat H_p2 = Hess_loglik(th + eps*e2);
  arma::mat H_m2 = Hess_loglik(th - eps*e2);

  double l111 = (H_p1(0,0) - H_m1(0,0)) / (2.0*eps);
  double l112 = (H_p2(0,0) - H_m2(0,0)) / (2.0*eps);
  double l122 = (H_p2(0,1) - H_m2(0,1)) / (2.0*eps);
  double l222 = (H_p2(1,1) - H_m2(1,1)) / (2.0*eps);

  double g_alp = g1[0];
  double g_gam = g1[1];

  double g_aa = g2(0,0);
  double g_ag = g2(0,1);
  double g_gg = g2(1,1);

  double out =
    g +
    0.5 * (
      (g_aa + 2.0*g_alp*r1)*s11 +
      (g_ag + 2.0*g_gam*r1)*s12 +
      (g_ag + 2.0*g_alp*r2)*s12 +
      (g_gg + 2.0*g_gam*r2)*s22
    ) +
    0.5 * (
      (g_alp*s11 + g_gam*s12) *
      (l111*s11 + 2.0*l112*s12 + l122*s22) +
      (g_alp*s12 + g_gam*s22) *
      (l112*s11 + 2.0*l122*s12 + l222*s22)
    );

  if(!R_finite(out) || out <= 0.0) return NA_REAL;

  return out;
}

// [[Rcpp::export]]
Rcpp::List Lindley_GEL_cpp(arma::vec True_para,
                           arma::vec Init_para,
                           int n,
                           arma::vec X,
                           arma::vec L,
                           arma::vec R,
                           arma::vec T,
                           arma::vec Rprog,
                           int final_right_censored,
                           arma::vec Cpm_par,
                           arma::vec h,
                           double q,
                           std::string type) {
  
  arma::vec th = Init_para.subvec(0,1);
  
  NumericVector est(3), bias(3), mse(3);
  
  for(int j=0;j<3;j++) {
    double val = Lindley_core_cpp(
      th, X,L,R,T,Rprog, final_right_censored,
      h, Cpm_par, type, j+1, "GEL", q
    );
    
    est[j] = std::pow(val, -1.0/q);
    bias[j] = est[j] - True_para[j];
    mse[j] = std::pow(bias[j], 2.0);
  }
  
  est.names() = CharacterVector::create("beta","gamma","Cpm");
  
  return List::create(
    Named("est_GEL") = est,
    Named("bias_GEL") = bias,
    Named("mse_GEL") = mse
  );
}

// [[Rcpp::export]]
Rcpp::List Lindley_LINEX_cpp(arma::vec True_para,
                             arma::vec Init_para,
                             int n,
                             arma::vec X,
                             arma::vec L,
                             arma::vec R,
                             arma::vec T,
                             arma::vec Rprog,
                             int final_right_censored,
                             arma::vec Cpm_par,
                             arma::vec h,
                             double c,
                             std::string type) {
  
  arma::vec th = Init_para.subvec(0,1);
  
  NumericVector est(3), bias(3), mse(3);
  
  for(int j=0;j<3;j++) {
    double val = Lindley_core_cpp(
      th, X,L,R,T,Rprog, final_right_censored,
      h, Cpm_par, type, j+1, "LINEX", c
    );
    
    est[j] = -(1.0/c) * std::log(val);
    bias[j] = est[j] - True_para[j];
    mse[j] = std::pow(bias[j], 2.0);
  }
  
  est.names() = CharacterVector::create("beta","gamma","Cpm");
  
  return List::create(
    Named("est_LINEX") = est,
    Named("bias_LINEX") = bias,
    Named("mse_LINEX") = mse
  );
}
