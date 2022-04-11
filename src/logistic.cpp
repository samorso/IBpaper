// Logistic regression with misclassification (given FP and FN rates)
#include "common.h"

// ------------------
// Implementation
// ------------------
// --------------
// Misc
// --------------
struct Sigmoid {
  Sigmoid(){}
  const double operator()(const double& x) const {return 0.1e1 / (0.1e1 + std::exp(-x));}
};

struct Indicator {
  Indicator(){}
  const double operator()(const double &t) const {return (t <= 0.5) ? 0.0 : 1.0;}
};

// --------------
// Logistic regression
// --------------
class mle_logistic: public Numer::MFuncGrad
{
private:
  const Eigen::ArrayXd y;
  const Eigen::MatrixXd x;
  const unsigned int n = y.size();
  const unsigned int p = x.cols();
  double fp;
  double fn;

public:
  mle_logistic(const Eigen::ArrayXd& y_,const Eigen::MatrixXd& x_,
               const double& fp_, const double& fn_) :
  y(y_), x(x_), fp(fp_), fn(fn_){}
  double f_grad(Numer::Constvec& beta, Numer::Refvec grad);
};

double mle_logistic::f_grad(
    Numer::Constvec& beta,
    Numer::Refvec grad
){
  // data storage
  Eigen::ArrayXd sig(n);
  Eigen::ArrayXd sig_star(n);
  Eigen::VectorXd mu(n);
  Eigen::MatrixXd x1(n,p+1);
  double t1;

  // pre-computation
  x1.rightCols(p) = x;
  x1.col(0) = Eigen::VectorXd::Constant(n,1.0);
  mu = x1 * beta;
  sig = mu.unaryExpr(Sigmoid());
  t1 = 1.0 - fp - fn;
  sig_star = sig * t1 + fp;


  // computation
  // objective function
  const double f = -((0.1e1 - y) * (0.1e1 - sig_star).log()).sum() / n - (y * sig_star.log()).sum() / n;

  // gradient
  grad = -t1 * x1.transpose() * ((y - sig_star) * sig * (1.0 - sig) / sig_star / (1.0 - sig_star)).matrix() / n;
  return f;
}

//' MLE for logistic regression with misclassified responses
//'
//' @param x a n x p matrix of design
//' @param y a n-vector of response
//' @param fp false positive rate
//' @param fn false negative rate
//' @export
// [[Rcpp::export]]
Eigen::VectorXd logistic_misclassification_mle(
    Eigen::MatrixXd& x,
    Eigen::ArrayXd& y,
    double fp,
    double fn
){
  // Regress
  unsigned int p = x.cols() + 1;
  double fopt;
  Eigen::VectorXd beta(p);
  beta.setZero();
  double prop = y.mean();
  beta(0) = std::log(prop) - std::log(1.0 - prop);
  mle_logistic f(y,x,fp,fn);
  Numer::optim_lbfgs(f,beta,fopt);
  return beta;
}

Rcpp::List check(
    Eigen::ArrayXd& y,
    Eigen::MatrixXd& x,
    Eigen::VectorXd& beta,
    double fp,
    double fn
){
  double fopt;
  Eigen::VectorXd grad(beta.size());
  mle_logistic f(y,x,fp,fn);
  fopt = f.f_grad(beta,grad);

  return Rcpp::List::create(
    Rcpp::Named("of") = fopt,
    Rcpp::Named("grad") = grad
  );
}


//' Inverse Fisher information matrix for logistic regression with misclassified
//' responses
//'
//' @param x a n x p matrix of design
//' @param beta a p-vector of parameter
//' @param fp false positive rate
//' @param fn false negative rate
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd inverse_FIM(
    Eigen::MatrixXd& x,
    Eigen::VectorXd& beta,
    double fp,
    double fn
){
  unsigned int n = x.rows();
  unsigned int p = x.cols();

  // data storage
  Eigen::ArrayXd sig(n);
  Eigen::ArrayXd sig_star(n);
  Eigen::VectorXd mu(n);
  Eigen::MatrixXd x1(n,p+1);
  Eigen::MatrixXd fish(p,p);
  double t1;

  // pre-computation
  x1.rightCols(p) = x;
  x1.col(0) = Eigen::VectorXd::Constant(n,1.0);
  mu = x1 * beta;
  sig = mu.unaryExpr(Sigmoid());
  t1 = 1.0 - fp - fn;
  sig_star = sig * t1 + fp;

  // Fisher information matrix
  fish = x1.transpose() * (sig * sig * (1.0 - sig) * (1.0 - sig) / sig_star / (1.0 - sig_star)).matrix().asDiagonal() * x1;
  return fish.inverse() / t1 / t1;
}

// ------------------
// Parametric Bootstrap
// ------------------
class mle_logistic2 : public Numer::MFuncGrad
{
private:
  const Eigen::ArrayXd y;
  const Eigen::MatrixXd x;
  const unsigned int n = y.size();
  const unsigned int p = x.cols();

public:
  mle_logistic2(const Eigen::ArrayXd& y_,const Eigen::MatrixXd& x_) : y(y_), x(x_) {}
  double f_grad(Numer::Constvec& beta, Numer::Refvec grad);
};

double mle_logistic2::f_grad(
    Numer::Constvec& beta,
    Numer::Refvec grad
){
  // data storage
  Eigen::ArrayXd sig(n);
  Eigen::VectorXd v(n);
  Eigen::MatrixXd x1(n,p+1);

  // pre-computation
  x1.rightCols(p) = x;
  x1.col(0) = Eigen::VectorXd::Constant(n,1.0);
  v = x1 * beta;
  sig = v.unaryExpr(Sigmoid());

  // computation
  // objective function
  const double f = -(0.1e1 - sig).log().sum() / n - v.dot(y.matrix()) / n;

  // gradient
  grad = x1.transpose() * (sig - y).matrix() / n;
  return f;
}

Eigen::ArrayXd r_logistic(
    Eigen::VectorXd& beta,
    Eigen::MatrixXd& x,
    unsigned int seed,
    double fp=0,
    double fn=0
){
  // data storage
  unsigned int p = x.cols();
  unsigned int n = x.rows();
  Eigen::ArrayXd sig(n);
  Eigen::ArrayXd sig2(n);
  Eigen::VectorXd v(n);
  Eigen::MatrixXd x1(n,p+1);
  Eigen::ArrayXd y(n);
  double u;

  // pre-computation
  x1.rightCols(p) = x;
  x1.col(0) = Eigen::VectorXd::Constant(n,1.0);
  v = x1 * beta;
  sig = v.unaryExpr(Sigmoid());
  sig2 = (1.0 - fp - fn) * sig + fp;
  std::mt19937_64 engine(seed);  // Mersenne twister random number engine
  std::uniform_real_distribution<double> unif(0.0, 1.0);

  for(unsigned int i(0);i<n;++i){
    u = unif(engine);
    y(i) = (u > sig2(i)) ? 0.0 : 1.0;
  }
  return y;
}

//' Parametric bootstrap for logistic regression with misclassified responses
//'
//' @param beta a p-vector of parameter
//' @param x a n x p design matrix
//' @param B the number of bootstrap replicates
//' @param seed for random number generator
//' @param ncores number of cores for parallelisation
//' @param fp false positive rate
//' @param fn false negative rate
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd par_bootstrap_mle(
    Eigen::VectorXd& beta,
    Eigen::MatrixXd& x,
    unsigned int B,
    unsigned int seed,
    unsigned int ncores,
    double fp=0,
    double fn=0
){
  unsigned int p = beta.size();
  Eigen::MatrixXd boot(B,p);

#pragma omp parallel for num_threads(ncores)
  for(unsigned int i=0; i<B; ++i){
    Eigen::VectorXd theta = beta;
    unsigned int se = seed + i;
    Eigen::ArrayXd y = r_logistic(beta,x,se,fp,fn);
    double fopt;
    mle_logistic2 f(y,x);
    int maxit = 300;
    double eps_f = 1e-8;
    double eps_g = 1e-7;
    Numer::optim_lbfgs(f, theta, fopt, maxit, eps_f, eps_g);
    boot.row(i) = theta;
  }

  return boot;
}
