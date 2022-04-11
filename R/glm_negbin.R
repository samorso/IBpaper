#' @title MLE for negative binomial with interfered responses
#' @description
#' Implementation of the MLE for negative binomial regression
#' with interfered responses: one observes \eqn{\max(y_i,z_i)},
#' where \eqn{y_i} is the negative binomial response and \eqn{z_i}
#' is a Poisson random variable with mean \eqn{\lambda}.
#' @param y response vector
#' @param x design matrix
#' @param lambda mean parameter for Poisson censoring process
#' @param maxit maximum number of iteration
#' @param epsilon tolerance parameter
#' @param trace boolean for printing information
#' @export
#' @importFrom stats coef optim optimize
#' @importFrom MASS glm.nb
glm_negbin <- function(y, x, lambda, maxit=50L, epsilon=1e-07, trace=FALSE){
  if(is.null(lambda) || lambda <= 0 || !is.numeric(lambda)) stop("lambda should be a positive value")

  fit <- glm.nb(y ~ x)
  bh <- coef(fit)
  th <-  1.0 / fit$theta

  d1 <- sqrt(2.0 * max(1.0, fit$df.residual))
  del <- 1
  Lm <- logLike_negbin(beta = bh, alpha = th, y = y, x = cbind(1,x), lambda = lambda)
  Lm0 <- Lm + 2.0 * d1
  iter <- 0L

  # alternate fitting beta and alpha
  while((iter <- iter + 1L) <= maxit && (abs(Lm0 - Lm)/d1 + abs(del)) > epsilon){
    b0 <- bh
    t0 <- th

    # fit beta
    fit <- optim(par = bh, fn = nll_max_beta, method = "BFGS", alpha = th, y = y, x = cbind(1,x), lambda = lambda)
    bh <- fit$par

    # fit alpha
    fit <- optimize(f = nll_max_alpha, interval = c(1e-03,1e02), beta = bh, y = y, x = cbind(1,x), lambda = lambda)
    th <- fit$minimum

    # update values
    del <- t0 - th
    Lm0 <- Lm
    Lm <- logLike_negbin(beta = bh, alpha = th, y = y, x = cbind(1,x), lambda = lambda)

    # info
    if(trace) message(sprintf("Theta(%d) = %f, error = %f", iter, signif(th), signif(abs(Lm0 - Lm)/d1 + abs(del))), domain = NA)
  }

  if(iter > maxit) warning("alternation limit reached")

  list(
    par = c(bh, th),
    iter = iter,
    of = Lm
  )
}
