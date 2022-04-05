#ifndef _COMMON
#define _COMMON
// ------------------
// Headers
// ------------------
// [[Rcpp::depends(RcppEigen,RcppNumerical,BH)]]

// Enable C++14 via the Rcpp plugin
// [[Rcpp::plugins("cpp14")]]

// Libraries
#include <RcppEigen.h>
#include <random>
#include <math.h>
#include <numeric>
#include <algorithm>
#include <vector>
#include <RcppNumerical.h>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/distributions/beta.hpp>

// using namespace Eigen;
// using namespace std;

// Add a flag to enable OpenMP at compile time
// [[Rcpp::plugins(openmp)]]

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

// Change error handling for Boost functions
// Define a specific policy:
typedef boost::math::policies::policy<
  boost::math::policies::digits10<5>,
  boost::math::policies::overflow_error<boost::math::policies::ignore_error>
> my_policy;

#endif
