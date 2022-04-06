# IBpaper

  [![R-CMD-check](https://github.com/samorso/IBpaper/workflows/R-CMD-check/badge.svg)](https://github.com/samorso/IBpaper/actions)

<img src="https://img.shields.io/badge/C%2B%2B-00599C?style=for-the-badge&logo=c%2B%2B&logoColor=white"> + <img src="https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white"> inside

## Overview 

The purpose of this package is to share the code my co-authors and myself used for our paper [Guerrier et
al (2020)](https://arxiv.org/pdf/2002.08757.pdf).

In order to install the package:

``` r
## if not installed
## install.packages("remotes")
remotes::install_github("samorso/IBpaper")
```

We highly rely on the `ib` package (v0.2.0): 

``` r
install.packages("ib")
```

You also need the following packages:
``` r
install.packages(c("betareg", "BH", "MASS", "Rcpp", "RcppEigen", "RcppNumerical"))
```

## Usage 
We essentially provides code for running simulation for logistic,
negative binomial and beta regressions with different assumptions on the data generating process, and data for a case study.

Here is a simple way to run a simulation for a logistic regression:
``` r
library(IBpaper)
library(ib)

x <- matrix(rnorm(300), ncol = 3) # design matrix
beta <- 1:4 # regression coefficients
logistic_object <- make_logistic(x, beta) # see ?`make_logistic`
y <- simulation(logistic_object) # from `ib` package

fit_mle <- glm(y ~ x, family = binomial(link = "logit")) # fit logistic regression
fit_jini <- ib(fit_mle, control=list(H=200, verbose=TRUE)) # iterative bootstrap procedure from `ib` package
results <- cbind(coef(fit_mle),coef(fit_jini))
colnames(results) <- c("MLE","JINI")
results
```

``` r
                 MLE     JINI
(Intercept) 1.294089 1.062380
x1          2.309151 1.952350
x2          3.793111 3.045767
x3          5.045143 4.091039
```

