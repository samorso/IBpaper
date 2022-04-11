// Negative binomial regression with random Poisson left censoring process
#include "common.h"

// ------------------
// Implementation
// ------------------
// ------------------
// Log-likelihood
// ------------------
double trim_pois_pmf(
    double y,
    double lambda,
    double cens
){
  double z;
  boost::math::poisson_distribution<>  pois_distr(lambda);
  if(y < cens){
    z = pdf(pois_distr, y);
  } else {
    z = cdf(complement(pois_distr, cens - 1));
  }
  return z;
}

double trim_pois_cdf(
    double y,
    double lambda,
    double cens
){
  if(y<0){
    return 0.0;
  }
  double z;
  if(y < cens){
    boost::math::poisson_distribution<>  pois_distr(lambda);
    z = cdf(pois_distr, y);
  } else {
    z = 1.0;
  }
  return z;
}

double pois_cdf(
    double y,
    double lambda
){
  if(y<0){
    return 0.0;
  }
  double z;
  boost::math::poisson_distribution<>  pois_distr(lambda);
  z = cdf(pois_distr, y);
  return z;
}

double negbin_pmf(
    double y,
    double mu,
    double alpha
){
  double r = 1.0 / alpha;
  double p = r / (r + mu);
  boost::math::negative_binomial_distribution<>  negbin_distr(r,p);
  return pdf(negbin_distr, y);
}

double negbin_cdf(
    double y,
    double mu,
    double alpha
){
  if(y<0){
    return 0.0;
  }
  double r = 1.0 / alpha;
  double p = r / (r + mu);
  boost::math::negative_binomial_distribution<>  negbin_distr(r,p);
  return cdf(negbin_distr, y);
}

double nll_min(
    const Eigen::VectorXd& theta,
    const Eigen::VectorXd& y,
    const Eigen::MatrixXd& x,
    const double& lambda,
    const double& cens
){
  unsigned int n = y.size();
  unsigned int p = theta.size();
  Eigen::VectorXd eta(n);
  double mu, alpha, of(0.0);

  eta = x * theta.head(p-1);
  alpha = std::exp(theta(p-1));

  for(unsigned int i(0); i<n; ++i){
    mu = std::exp(eta(i));
    of -= std::log(negbin_pmf(y(i),mu,alpha) + trim_pois_pmf(y(i),lambda,cens) +
      negbin_cdf(y(i)-1,mu,alpha) * trim_pois_cdf(y(i)-1,lambda,cens) -
      negbin_cdf(y(i),mu,alpha) * trim_pois_cdf(y(i),lambda,cens));
  }

  return of / n;
}

double nll_max(
    const Eigen::VectorXd& theta,
    const Eigen::VectorXd& y,
    const Eigen::MatrixXd& x,
    const double& lambda
){
  unsigned int n = y.size();
  unsigned int p = theta.size();
  Eigen::VectorXd eta(n);
  double mu, alpha, of(0.0);

  eta = x * theta.head(p-1);
  alpha = std::exp(theta(p-1));

  for(unsigned int i(0); i<n; ++i){
    mu = std::exp(eta(i));
    of -= std::log(negbin_cdf(y(i),mu,alpha) * pois_cdf(y(i),lambda) -
      negbin_cdf(y(i)-1,mu,alpha) * pois_cdf(y(i)-1,lambda));
  }

  return of / n;
}

//' Negative log-likelihood for beta coefficients (negative binomial regression)
//'
//' @param beta a p-vector of coefficients
//' @param alpha parameter of negative binomial
//' @param y a n-vector of response
//' @param x a n x p matrix of design
//' @param lambda mean of Poisson censoring process
//' @export
// [[Rcpp::export]]
double nll_max_beta(
    const Eigen::VectorXd& beta,
    const double& alpha,
    const Eigen::VectorXd& y,
    const Eigen::MatrixXd& x,
    const double& lambda
){
  unsigned int n = y.size();
  Eigen::VectorXd eta(n);
  double mu, of(0.0);

  eta = x * beta;

  for(unsigned int i(0); i<n; ++i){
    mu = std::exp(eta(i));
    of -= std::log(negbin_cdf(y(i),mu,alpha) * pois_cdf(y(i),lambda) -
      negbin_cdf(y(i)-1,mu,alpha) * pois_cdf(y(i)-1,lambda));
  }

  return of / n;
}

//' Negative log-likelihood for the overdispersion parameter (negative binomial regression)
//'
//' @param alpha parameter of negative binomial
//' @param beta a p-vector of coefficients
//' @param y a n-vector of response
//' @param x a n x p matrix of design
//' @param lambda mean of Poisson censoring process
//' @export
// [[Rcpp::export]]
double nll_max_alpha(
    const double& alpha,
    const Eigen::VectorXd& beta,
    const Eigen::VectorXd& y,
    const Eigen::MatrixXd& x,
    const double& lambda
){
  unsigned int n = y.size();
  Eigen::VectorXd eta(n);
  double mu, of(0.0);

  eta = x * beta;

  for(unsigned int i(0); i<n; ++i){
    mu = std::exp(eta(i));
    of -= std::log(negbin_cdf(y(i),mu,alpha) * pois_cdf(y(i),lambda) -
      negbin_cdf(y(i)-1,mu,alpha) * pois_cdf(y(i)-1,lambda));
  }

  return of / n;
}

//' Log-likelihood for negative binomial with interfered responses
//'
//' @param beta a p-vector of coefficients
//' @param alpha overdispersion parameter of negative binomial
//' @param y a n-vector of response
//' @param x a n x p matrix of design
//' @param lambda mean of Poisson censoring process
//' @export
// [[Rcpp::export]]
double logLike_negbin(
    const Eigen::VectorXd& beta,
    const double& alpha,
    const Eigen::VectorXd& y,
    const Eigen::MatrixXd& x,
    const double& lambda
){
  unsigned int n = y.size();
  Eigen::VectorXd eta(n);
  double mu, of(0.0);

  eta = x * beta;

  for(unsigned int i(0); i<n; ++i){
    mu = std::exp(eta(i));
    of += std::log(negbin_cdf(y(i),mu,alpha) * pois_cdf(y(i),lambda) -
      negbin_cdf(y(i)-1,mu,alpha) * pois_cdf(y(i)-1,lambda));
  }

  return of;
}
