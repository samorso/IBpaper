// Beta regression with rounding (to the nearest tenth)
#include "common.h"

// ------------------
// Implementation
// ------------------
// ------------------
// Link function
// ------------------
double logit(double x){return std::log(x / (0.1e1 - x));}
double logistic(double x){return 0.1e1 / (0.1e1 + std::exp(-x));}

// ------------------
// Log-likelihood
// ------------------
double ll(
    const Eigen::VectorXd& theta,
    const Eigen::MatrixXd& y,
    const Eigen::MatrixXd& x
){
  unsigned int n = y.rows();
  unsigned int p = theta.size();
  Eigen::VectorXd eta(n);
  double of(0.0),phi,mu,t1,t2;

  phi = std::exp(theta(p-1));
  eta = x * theta.head(p-1);

  for(unsigned int i(0);i<n;++i){
    mu = logistic(eta(i));
    t1 = mu * phi;
    t2 = phi - t1;
    of += std::lgamma(t1) - std::lgamma(phi) + std::lgamma(t2) - (t1 - 0.1e1) * y(i,0) - (t2 - 0.1e1) * y(i,1);
  }

  return of / n;
}

Eigen::VectorXd d_ll(
    const Eigen::VectorXd& theta,
    const Eigen::MatrixXd& y,
    const Eigen::MatrixXd& x
){
  unsigned int n = y.rows();
  unsigned int p = theta.size();
  Eigen::VectorXd eta(n);
  double phi,mu,t1,t2,t3,t4,delta;
  Eigen::VectorXd grad = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd d_mu(p-1);

  phi = std::exp(theta(p-1));
  eta = x * theta.head(p-1);
  t4 = boost::math::digamma(phi, my_policy());

  for(unsigned int i(0);i<n;++i){
    mu = logistic(eta(i));
    t1 = mu * phi;
    t2 = phi - t1;
    t3 = mu - mu * mu;
    delta = boost::math::digamma(t1, my_policy()) - boost::math::digamma(t2, my_policy()) - y(i,0) + y(i,1);
    d_mu = t3 * x.row(i);
    grad.head(p-1) += phi * delta * d_mu;
    grad(p-1) += mu * delta + boost::math::digamma(t2, my_policy()) - y(i,1) - t4;
  }

  grad(p-1) *= phi;

  return grad / n;
}

// ------------------
// M-Step : Maximum likelihood
// ------------------
class mle_betareg: public Numer::MFuncGrad
{
private:
  const Eigen::MatrixXd y;
  const Eigen::MatrixXd x;

public:
  mle_betareg(const Eigen::MatrixXd y_, const Eigen::MatrixXd x_) : y(y_), x(x_) {}
  double f_grad(Numer::Constvec& theta, Numer::Refvec gr);
};

double mle_betareg::f_grad(
    Numer::Constvec& theta,
    Numer::Refvec gr
){
  const double of = ll(theta,y,x);
  gr = d_ll(theta,y,x);
  return of;
}

Rcpp::List optim_mle_betareg(
    Eigen::VectorXd& theta,
    Eigen::MatrixXd& y,
    Eigen::MatrixXd& x,
    int maxit = 300,
    double eps_f = 1e-6,
    double eps_g = 1e-6
){
  double fopt;
  mle_betareg f(y,x);
  int status = Numer::optim_lbfgs(f, theta, fopt, maxit, eps_f, eps_g);
  return Rcpp::List::create(
    Rcpp::Named("par") = theta,
    Rcpp::Named("value") = fopt,
    Rcpp::Named("conv") = status
  );
}

Rcpp::List of_mle_betareg(
    Eigen::VectorXd& theta,
    Eigen::MatrixXd& y,
    Eigen::MatrixXd& x
){
  double fopt;
  unsigned int p = theta.size();
  Eigen::VectorXd gr(p);
  mle_betareg f(y,x);
  fopt = f.f_grad(theta,gr);
  return Rcpp::List::create(
    Rcpp::Named("f") = fopt,
    Rcpp::Named("grad") = gr
  );
}

Eigen::VectorXd m_step(
    Eigen::VectorXd& theta,
    Eigen::MatrixXd& y,
    Eigen::MatrixXd& x
){
  double fopt;
  mle_betareg f(y,x);
  Numer::optim_lbfgs(f, theta, fopt);
  return theta;
}

// ------------------
// E-Step
// ------------------
class logx_beta_pdf: public Numer::Func
{
private:
  double mu;
  double phi;
  double t1 = mu * phi;
  double t2 = phi - t1;

public:
  logx_beta_pdf(const double mu_, const double phi_) : mu(mu_), phi(phi_) {}
  double operator()(const double& x) const
  {
    return std::log(x) * std::tgamma(phi) * std::pow(x, t1 - 1.0) * std::pow(1.0 - x, t2 - 1.0) / std::tgamma(t1) / std::tgamma(t2);
  }
};

class log1x_beta_pdf: public Numer::Func
{
private:
  double mu;
  double phi;
  double t1 = mu * phi;
  double t2 = phi - t1;

public:
  log1x_beta_pdf(const double mu_, const double phi_) : mu(mu_), phi(phi_) {}
  double operator()(const double& x) const
  {
    return std::log(1.0 - x) * std::tgamma(phi) * std::pow(x, t1 - 1.0) * std::pow(1.0 - x, t2 - 1.0) / std::tgamma(t1) / std::tgamma(t2);
  }
};

Rcpp::List integrate_test(
    const double mu,
    const double phi,
    const double lower,
    const double upper
)
{
  double alpha = mu * phi;
  double beta = phi - alpha;
  boost::math::beta_distribution<>  beta_distr(alpha, beta);
  double t1 = cdf(beta_distr, upper) - cdf(beta_distr, lower);
  logx_beta_pdf f(mu, phi);
  double err_est;
  int err_code;
  const double res = Numer::integrate(f, lower, upper, err_est, err_code) / t1;
  return Rcpp::List::create(
    Rcpp::Named("approximate") = res,
    Rcpp::Named("error_estimate") = err_est,
    Rcpp::Named("error_code") = err_code
  );
}

Eigen::MatrixXd e_step(
    const Eigen::VectorXd& theta,
    const Eigen::VectorXd& y,
    const Eigen::MatrixXd& x
){
  unsigned int n = x.rows();
  unsigned int p = theta.size();
  double alpha, beta, mu, phi, t1, upper, lower, t2, t3;
  double err_est;
  int err_code;
  Eigen::VectorXd eta(n);
  Eigen::MatrixXd logs(n,2);
  double mygrid[] = {0.00,0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.00};
  std::vector<double> grid(mygrid,mygrid+12);
  std::vector<double>::iterator it;

  phi = std::exp(theta(p-1));
  eta = x * theta.head(p-1);

  for(unsigned int i=0;i<n;++i){
    it = std::lower_bound(grid.begin(), grid.end(), y(i));
    lower = *(it-1);
    if(lower<0){
      lower = 0.0;
    }
    upper = *it;
    if(upper>=1.0){
      upper = 1. - 1e-6;
    }
    mu = logistic(eta(i));
    alpha = mu * phi;
    beta = phi - alpha;
    boost::math::beta_distribution<>  beta_distr(alpha, beta);
    if(lower == 0.){
      t2 = 0.;
    } else {
      t2 = cdf(beta_distr, lower);
    }
    if(upper == 1.){
      t3 = 1.;
    } else {
      t3 = cdf(beta_distr, upper);
    }
    t1 = t3 - t2;
    logx_beta_pdf f(mu, phi);
    logs(i,0) = Numer::integrate(f, lower, upper, err_est, err_code) / t1;
    log1x_beta_pdf f1(mu, phi);
    logs(i,1) = Numer::integrate(f1, lower, upper, err_est, err_code) / t1;
  }

  return logs;
}

// ------------------
// EM algorithm
// ------------------
//' EM algorithm to compute the MLE for beta regression with rounded responses
//'
//' @param y a n-vector of response
//' @param x a n x p matrix of design
//' @param theta vector of parameter
//' @param maxit maximum number of iteration
//' @param eps tolerance
//' @param verbose boolean for printing information
//' @export
// [[Rcpp::export]]
Rcpp::List em(
    Eigen::VectorXd& y,
    Eigen::MatrixXd& x,
    Eigen::VectorXd& theta,
    unsigned int maxit=30,
    double eps=1e-7,
    bool verbose = true
){
  unsigned int n = y.size();
  unsigned int p = theta.size();
  double test = eps + 1.;
  Eigen::VectorXd t0 = theta;
  Eigen::VectorXd t1(p);
  Eigen::MatrixXd y2(n,2);
  y2.col(0) = y.array().log();
  y2.col(1) = (1.0 - y.array()).log();

  unsigned int it(0);

  do{
    // M-step
    t1 = m_step(t0,y2,x);

    // E-step
    y2 = e_step(t1,y,x);

    // test
    if(it > 0){
      test = (t0-t1).squaredNorm();
    }

    // update
    t0 = t1;
    ++it;
  } while (it < maxit && test > eps);

  return Rcpp::List::create(
    Rcpp::Named("par") = t0,
    Rcpp::Named("iteration") = it,
    Rcpp::Named("error") = test
  );
}
