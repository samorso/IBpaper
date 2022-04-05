## https://github.com/wch/r-source/blob/79298c499218846d14500255efd622b5021c10ec/src/library/stats/R/glm.R

#' @importFrom stats .getXlevels as.formula binomial cor dbeta delete.response gaussian glm.control is.empty.model lm.wfit make.link model.extract model.matrix model.offset model.response model.weights qlogis terms var
#' @importFrom Formula Formula as.Formula model.part
#' @importFrom betareg betareg.control
#' @importFrom utils getFromNamespace
makeglm <- function(formula, family = gaussian, data, coefficients,
                    weights, subset, na.action, start = NULL,
                    etastart, mustart, offset, control = list(...),
                    model = TRUE, method = "glm.nofit", x = FALSE, y = TRUE,
                    singular.ok = TRUE, contrasts = NULL, ...)
{
  cal <- match.call()
  ## coefficients
  if(is.null(coefficients))
    stop("'coefficients' must be provided")


  ## family
  if(is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if(is.function(family)) family <- family()
  if(is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  ## extract x, y, etc from the model formula and frame
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "etastart", "mustart", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  ## need stats:: for non-standard evaluation
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if(identical(method, "model.frame")) return(mf)

  if (!is.character(method) && !is.function(method))
    stop("invalid 'method' argument")
  ## for back-compatibility in return result
  if (identical(method, "glm.fit"))
    control <- do.call("glm.control", control)

  mt <- attr(mf, "terms") # allow model.frame to have updated it

  Y <- model.response(mf, "any") # e.g. factors are allowed
  ## avoid problems with 1D arrays, but keep names
  if(length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if(!is.null(nm)) names(Y) <- nm
  }
  ## null model support
  X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts) else matrix(,NROW(Y), 0L)
  ## avoid any problems with 1D or nx1 arrays by as.vector.
  weights <- as.vector(model.weights(mf))
  if(!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
  ## check weights and offset
  if( !is.null(weights) && any(weights < 0) )
    stop("negative weights not allowed")

  offset <- as.vector(model.offset(mf))
  if(!is.null(offset)) {
    if(length(offset) != NROW(Y))
      stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                    length(offset), NROW(Y)), domain = NA)
  }
  ## these allow starting values to be expressed in terms of other vars.
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")

  ## We want to set the name on this call and the one below for the
  ## sake of messages from the fitter function
  fit <- eval(call(if(is.function(method)) "method" else method,
                   x = X, y = Y, coefficients = coefficients, weights = weights, start = start,
                   etastart = etastart, mustart = mustart,
                   offset = offset, family = family, control = control,
                   intercept = attr(mt, "intercept") > 0L, singular.ok = singular.ok))

  ## This calculated the null deviance from the intercept-only model
  ## if there is one, otherwise from the offset-only model.
  ## We need to recalculate by a proper fit if there is intercept and
  ## offset.
  ##
  ## The glm.fit calculation could be wrong if the link depends on the
  ## observations, so we allow the null deviance to be forced to be
  ## re-calculated by setting an offset (provided there is an intercept).
  ## Prior to 2.4.0 this was only done for non-zero offsets.
  if(length(offset) && attr(mt, "intercept") > 0L) {
    fit2 <-
      eval(call(if(is.function(method)) "method" else method,
                x = X[, "(Intercept)", drop=FALSE], y = Y,
                coefficients = coefficients,
                ## starting values potentially required (PR#16877):
                mustart = fit$fitted.values,
                weights = weights, offset = offset, family = family,
                control = control, intercept = TRUE))
    ## That fit might not have converged ....
    if(!fit2$converged)
      warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?")
    fit$null.deviance <- fit2$deviance
  }
  if(model) fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  if(x) fit$x <- X
  if(!y) fit$y <- NULL
  structure(c(fit,
              list(call = cal, formula = formula,
                   terms = mt, data = data,
                   offset = offset, control = control, method = method,
                   contrasts = attr(X, "contrasts"),
                   xlevels = .getXlevels(mt, mf))),
            class = c(fit$class, c("glm", "lm")))
}

glm.nofit <- function (x, y, coefficients, weights = rep.int(1, nobs), start = NULL,
                       etastart = NULL, mustart = NULL, offset = rep.int(0, nobs),
                       family = gaussian(), control = list(), intercept = TRUE,
                       singular.ok = TRUE)
{
  control <- do.call("glm.control", control)
  x <- as.matrix(x)
  xnames <- dimnames(x)[[2L]]
  ynames <- if(is.matrix(y)) rownames(y) else names(y)
  conv <- FALSE
  nobs <- NROW(y)
  nvars <- ncol(x)
  EMPTY <- nvars == 0
  ## define weights and offset if needed
  if (is.null(weights))
    weights <- rep.int(1, nobs)
  if (is.null(offset))
    offset <- rep.int(0, nobs)

  ## get family functions:
  variance <- family$variance
  linkinv  <- family$linkinv
  if (!is.function(variance) || !is.function(linkinv) )
    stop("'family' argument seems not to be a valid family object", call. = FALSE)
  dev.resids <- family$dev.resids
  aic <- family$aic
  mu.eta <- family$mu.eta
  unless.null <- function(x, if.null) if (is.null(x)) if.null else x
  valideta <- unless.null(family$valideta, function(eta) TRUE)
  validmu  <- unless.null(family$validmu, function(mu) TRUE)
  if(is.null(mustart)) {
    ## calculates mustart and may change y and weights and set n (!)
    eval(family$initialize)
  } else {
    mukeep <- mustart
    eval(family$initialize)
    mustart <- mukeep
  }
  if(EMPTY) {
    eta <- rep.int(0, nobs) + offset
    if (!valideta(eta))
      stop("invalid linear predictor values in empty model", call. = FALSE)
    mu <- linkinv(eta)
    ## calculate initial deviance and coefficient
    if (!validmu(mu))
      stop("invalid fitted means in empty model", call. = FALSE)
    dev <- sum(dev.resids(y, mu, weights))
    w <- sqrt((weights * mu.eta(eta)^2)/variance(mu))
    residuals <- (y - mu)/mu.eta(eta)
    good <- rep_len(TRUE, length(residuals))
    boundary <- conv <- TRUE
    coef <- numeric()
    iter <- 0L
  } else {
    coefold <- NULL
    eta <- if(!is.null(etastart)){ etastart} else {
      if(!is.null(start))
        if (length(start) != nvars)
          stop(gettextf(
            "length of 'start' should equal %d and correspond to initial coefs for %s",
            nvars, paste(deparse(xnames), collapse=", ")),
            domain = NA)
      else {
        coefold <- start
        offset + as.vector(if (NCOL(x) == 1L) x * start else x %*% start)
      }
      else family$linkfun(mustart)
    }
    mu <- linkinv(eta)
    if (!(validmu(mu) && valideta(eta)))
      stop("cannot find valid starting values: please specify some", call. = FALSE)
    ## calculate initial deviance and coefficient
    devold <- sum(dev.resids(y, mu, weights))
    boundary <- conv <- FALSE

    ##------------- THE Iteratively Reweighting L.S. iteration -----------
    ##------------- No fit modification!
    for (iter in 1L) {
      good <- weights > 0
      varmu <- variance(mu)[good]
      if (anyNA(varmu))
        stop("NAs in V(mu)")
      if (any(varmu == 0))
        stop("0s in V(mu)")
      mu.eta.val <- mu.eta(eta)
      if (any(is.na(mu.eta.val[good])))
        stop("NAs in d(mu)/d(eta)")
      ## drop observations for which w will be zero
      good <- (weights > 0) & (mu.eta.val != 0)

      if (all(!good)) {
        conv <- FALSE
        warning(gettextf("no observations informative at iteration %d",
                         iter), domain = NA)
        break
      }
      z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
      w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
      ## call Fortran code via C wrapper
      # fit <- .Call(C_Cdqrls, x[good, , drop = FALSE] * w, z * w,
      #              min(1e-7, control$epsilon/1000), check=FALSE)
      fit <- list()
      ## replace by user-defined coefficients
      fit$coefficients <- coefficients
      ## stop if not enough parameters
      ## calculate updated values of eta and mu with the new coef:
      start <- fit$coefficients
      eta <- drop(x %*% start)
      mu <- linkinv(eta <- eta + offset)
      dev <- sum(dev.resids(y, mu, weights))
      if (control$trace)
        cat("Deviance = ", dev, " Iterations - ", iter, "\n", sep = "")
      ## check for divergence
      boundary <- FALSE
      if (!is.finite(dev)) {
        if(is.null(coefold))
          stop("no valid set of coefficients has been found: please supply starting values", call. = FALSE)
        warning("step size truncated due to divergence", call. = FALSE)
        ii <- 1
        while (!is.finite(dev)) {
          if (ii > control$maxit)
            stop("inner loop 1; cannot correct step size", call. = FALSE)
          ii <- ii + 1
          start <- (start + coefold)/2
          eta <- drop(x %*% start)
          mu <- linkinv(eta <- eta + offset)
          dev <- sum(dev.resids(y, mu, weights))
        }
        boundary <- TRUE
        if (control$trace)
          cat("Step halved: new deviance = ", dev, "\n", sep = "")
      }
      ## check for fitted values outside domain.
      if (!(valideta(eta) && validmu(mu))) {
        if(is.null(coefold))
          stop("no valid set of coefficients has been found: please supply starting values", call. = FALSE)
        warning("step size truncated: out of bounds", call. = FALSE)
        ii <- 1
        while (!(valideta(eta) && validmu(mu))) {
          if (ii > control$maxit)
            stop("inner loop 2; cannot correct step size", call. = FALSE)
          ii <- ii + 1
          start <- (start + coefold)/2
          eta <- drop(x %*% start)
          mu <- linkinv(eta <- eta + offset)
        }
        boundary <- TRUE
        dev <- sum(dev.resids(y, mu, weights))
        if (control$trace)
          cat("Step halved: new deviance = ", dev, "\n", sep = "")
      }
      ## check for convergence
      if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
        conv <- TRUE
        coef <- start
        break
      } else {
        devold <- dev
        coef <- coefold <- start
      }
    } ##-------------- end IRLS iteration -------------------------------

    # if (!conv)
    #   warning("glm.fit: algorithm did not converge", call. = FALSE)
    # if (boundary)
    #   warning("glm.fit: algorithm stopped at boundary value", call. = FALSE)
    # eps <- 10*.Machine$double.eps
    # if (family$family == "binomial") {
    #   if (any(mu > 1 - eps) || any(mu < eps))
    #     warning("glm.fit: fitted probabilities numerically 0 or 1 occurred", call. = FALSE)
    # }
    # if (family$family == "poisson") {
    #   if (any(mu < eps))
    #     warning("glm.fit: fitted rates numerically 0 occurred", call. = FALSE)
    # }
    ## If X matrix was not full rank then columns were pivoted,
    ## hence we need to re-label the names ...
    ## Original code changed as suggested by BDR---give NA rather
    ## than 0 for non-estimable parameters
    # if (fit$rank < nvars) coef[fit$pivot][seq.int(fit$rank+1, nvars)] <- NA
    # xxnames <- xnames[fit$pivot]
    ## update by accurate calculation, including 0-weight cases.
    residuals <-  (y - mu)/mu.eta(eta)
    ##        residuals <- rep.int(NA, nobs)
    ##        residuals[good] <- z - (eta - offset)[good] # z does not have offset in.
    # fit$qr <- as.matrix(fit$qr)
    nr <- min(sum(good), nvars)
    # if (nr < nvars) {
    Rmat <- diag(nvars)
    # Rmat[1L:nr, 1L:nvars] <- fit$qr[1L:nr, 1L:nvars]
    # }
    # else Rmat <- fit$qr[1L:nvars, 1L:nvars]
    Rmat <- as.matrix(Rmat)
    Rmat[row(Rmat) > col(Rmat)] <- 0
    names(coef) <- xnames
    # colnames(fit$qr) <- xxnames
    # dimnames(Rmat) <- list(xxnames, xxnames)
  }
  names(residuals) <- ynames
  names(mu) <- ynames
  names(eta) <- ynames
  # for compatibility with lm, which has a full-length weights vector
  wt <- rep.int(0, nobs)
  wt[good] <- w^2
  names(wt) <- ynames
  names(weights) <- ynames
  names(y) <- ynames
  # if(!EMPTY)
  #   names(fit$effects) <-
  #   c(xxnames[seq_len(fit$rank)], rep.int("", sum(good) - fit$rank))
  # ## calculate null deviance -- corrected in glm() if offset and intercept
  wtdmu <- if (intercept){sum(weights * y)/sum(weights)} else {linkinv(offset)}
  nulldev <- sum(dev.resids(y, wtdmu, weights))
  ## calculate df
  n.ok <- nobs - sum(weights==0)
  nulldf <- n.ok - as.integer(intercept)
  rank <- 0# if(EMPTY){ 0 }else {fit$rank}
  resdf  <- n.ok - rank
  ## calculate AIC
  aic.model <- aic(y, n.ok, mu, weights, dev) + 2*rank
  ##     ^^ is only initialize()d for "binomial" [yuck!]
  list(coefficients = coef, residuals = residuals, fitted.values = mu,
       # effects = if(!EMPTY) fit$effects, R = if(!EMPTY) Rmat,
       rank = rank,
       # qr = if(!EMPTY) structure(fit[c("qr", "rank", "qraux", "pivot", "tol")], class = "qr"),
       family = family,
       linear.predictors = eta, deviance = dev, aic = aic.model,
       null.deviance = nulldev, iter = iter, weights = wt,
       prior.weights = weights, df.residual = resdf, df.null = nulldf,
       y = y, converged = conv, boundary = boundary)
}

# https://github.com/cran/MASS/blob/master/R/negbin.R
makeglm.nb <- function(formula, data, weights, coefficients=NULL,
                       subset, na.action, start = NULL, etastart, mustart,
                       control = glm.control(...), method = "glm.fit",
                       model = TRUE, x = FALSE, y = TRUE, contrasts = NULL, ...,
                       init.theta=NULL, link = log)
{
  ## must provide parameters
  ## coefficients
  if(is.null(coefficients))
    stop("'coefficients' must be provided")

  ## theta
  if(is.null(init.theta)) stop("'theta' must be provided")

  loglik <- function(n, th, mu, y, w)
    sum(w*(lgamma(th + y) - lgamma(th) - lgamma(y + 1) + th * log(th) +
             y * log(mu + (y == 0)) - (th + y) * log(th + mu)))

  link <- substitute(link)
  fam0 <- if(missing(init.theta))
    do.call("poisson", list(link = link))
  else
    do.call("negative.binomial", list(theta = init.theta, link = link))

  mf <- Call <- match.call()
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "etastart", "mustart", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval.parent(mf)
  Terms <- attr(mf, "terms")
  if(method == "model.frame") return(mf)
  Y <- model.response(mf, "numeric")
  ## null model support
  X <- if (!is.empty.model(Terms)) model.matrix(Terms, mf, contrasts) else matrix(,NROW(Y),0)
  w <- model.weights(mf)
  if(!length(w)) w <- rep(1, nrow(mf))
  else if(any(w < 0)) stop("negative weights not allowed")
  offset <- model.offset(mf)
  ## these allow starting values to be expressed in terms of other vars.
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")
  n <- length(Y)
  # if(!missing(method)) {
  #   if(!exists(method, mode = "function"))
  #     stop(gettextf("unimplemented method: %s", sQuote(method)),
  #          domain = NA)
  #   glm.fitter <- get(method)
  # } else {
  #   method <- "glm.fit"
  #   glm.fitter <- stats::glm.fit
  # }
  # if(control$trace > 1) message("Initial fit:")
  fit <- glm.nofit(x = X, y = Y, coefficients, weights = w, start = start,
                   etastart = etastart, mustart = mustart,
                   offset = offset, family = fam0,
                   control = list(maxit=control$maxit,
                                  epsilon = control$epsilon,
                                  trace = control$trace > 1),
                   intercept = attr(Terms, "intercept") > 0)
  class(fit) <- c("glm", "lm")
  mu <- fit$fitted.values
  # th <- as.vector(theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
  #                            control$trace> 2)) # drop attributes
  th <- init.theta
  if(control$trace > 1)
    message(gettextf("Initial value for 'theta': %f", signif(th)),
            domain = NA)
  fam <- do.call("negative.binomial", list(theta = th, link = link))
  iter <- 0
  d1 <- sqrt(2 * max(1, fit$df.residual))
  d2 <- del <- 1
  g <- fam$linkfun
  Lm <- loglik(n, th, mu, Y, w)
  Lm0 <- Lm + 2 * d1
  # while((iter <- iter + 1) <= control$maxit &&
  #       (abs(Lm0 - Lm)/d1 + abs(del)/d2) > control$epsilon) {
  #   eta <- g(mu)
  #   fit <- glm.fitter(x = X, y = Y, w = w, etastart =
  #                       eta, offset = offset, family = fam,
  #                     control = list(maxit=control$maxit,
  #                                    epsilon = control$epsilon,
  #                                    trace = control$trace > 1),
  #                     intercept = attr(Terms, "intercept") > 0)
  #   t0 <- th
  #   th <- theta.ml(Y, mu, sum(w), w, limit=control$maxit,
  #                  trace = control$trace > 2)
  #   fam <- do.call("negative.binomial", list(theta = th, link = link))
  #   mu <- fit$fitted.values
  #   del <- t0 - th
  #   Lm0 <- Lm
  #   Lm <- loglik(n, th, mu, Y, w)
  #   if(control$trace) {
  #     Ls <- loglik(n, th, Y, Y, w)
  #     Dev <- 2 * (Ls - Lm)
  #     message(sprintf("Theta(%d) = %f, 2(Ls - Lm) = %f",
  #                     iter, signif(th), signif(Dev)), domain = NA)
  #   }
  # }
  # if(!is.null(attr(th, "warn"))) fit$th.warn <- attr(th, "warn")
  # if(iter > control$maxit) {
  #   warning("alternation limit reached")
  #   fit$th.warn <- gettext("alternation limit reached")
  # }

  # If an offset and intercept are present, iterations are needed to
  # compute the Null deviance; these are done here, unless the model
  # is NULL, in which case the computations have been done already
  #
  # if(length(offset) && attr(Terms, "intercept")) {
  #   null.deviance <-
  #     if(length(Terms))
  #       glm.fitter(X[, "(Intercept)", drop = FALSE], Y, w,
  #                  offset = offset, family = fam,
  #                  control = list(maxit=control$maxit,
  #                                 epsilon = control$epsilon,
  #                                 trace = control$trace > 1),
  #                  intercept = TRUE
  #       )$deviance
  #   else fit$deviance
  #   fit$null.deviance <- null.deviance
  # }
  class(fit) <- c("negbin", "glm", "lm")
  fit$terms <- Terms
  fit$formula <- as.vector(attr(Terms, "formula"))
  ## make result somewhat reproducible
  Call$init.theta <- signif(as.vector(th), 10)
  Call$link <- link
  fit$call <- Call
  if(model) fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  if(x) fit$x <- X
  if(!y) fit$y <- NULL
  fit$theta <- as.vector(th)
  fit$SE.theta <- attr(th, "SE")
  fit$twologlik <- as.vector(2 * Lm)
  fit$aic <- -fit$twologlik + 2*fit$rank + 2
  fit$contrasts <- attr(X, "contrasts")
  fit$xlevels <- .getXlevels(Terms, mf)
  fit$method <- method
  fit$control <- control
  fit$offset <- offset
  fit
}

library(Formula)
# https://github.com/cran/betareg/blob/master/R/betareg.R
makebetareg <- function(formula, data, coefficients, subset, na.action, weights, offset,
                        link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog"),
                        link.phi = NULL, type = c("ML", "BC", "BR"),
                        control = betareg.control(...),
                        model = TRUE, y = TRUE, x = FALSE, ...)
{
  ## call
  cl <- match.call()
  ## coefficients
  if(is.null(coefficients))
    stop("'coefficients' must be provided")

  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  ## formula
  oformula <- as.formula(formula)
  formula <- as.Formula(formula)
  if(length(formula)[2L] < 2L) {
    formula <- as.Formula(formula(formula), ~ 1)
    simple_formula <- TRUE
  } else {
    if(length(formula)[2L] > 2L) {
      formula <- Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts")
    }
    simple_formula <- FALSE
  }
  mf$formula <- formula

  ## evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## extract terms, model matrix, response
  mt <- terms(formula, data = data)
  mtX <- terms(formula, data = data, rhs = 1L)
  mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
  Y <- model.response(mf, "numeric")
  X <- model.matrix(mtX, mf)
  Z <- model.matrix(mtZ, mf)

  ## sanity checks
  if(length(Y) < 1) stop("empty model")
  if(!(min(Y) > 0 & max(Y) < 1)) stop("invalid dependent variable, all observations must be in (0, 1)")

  ## convenience variables
  n <- length(Y)

  ## type of estimator
  type <- match.arg(type)

  ## links
  if(is.character(link)) link <- match.arg(link)
  if(is.null(link.phi)) link.phi <- if(simple_formula) "identity" else "log"
  if(is.character(link.phi)) link.phi <- match.arg(link.phi, c("identity", "log", "sqrt"))

  ## weights
  weights <- model.weights(mf)
  if(is.null(weights)) weights <- 1
  if(length(weights) == 1) weights <- rep.int(weights, n)
  weights <- as.vector(weights)
  names(weights) <- rownames(mf)

  ## offsets
  expand_offset <- function(offset) {
    if(is.null(offset)) offset <- 0
    if(length(offset) == 1) offset <- rep.int(offset, n)
    as.vector(offset)
  }
  ## in mean part of formula
  offsetX <- expand_offset(model.offset(model.part(formula, data = mf, rhs = 1L, terms = TRUE)))
  ## in precision part of formula
  offsetZ <- expand_offset(model.offset(model.part(formula, data = mf, rhs = 2L, terms = TRUE)))
  ## in offset argument (used for mean)
  if(!is.null(cl$offset)) offsetX <- offsetX + expand_offset(mf[, "(offset)"])
  ## collect
  offset <- list(mean = offsetX, precision = offsetZ)

  ## call the actual workhorse: betareg.fit()
  rval <- betareg.nofit(X, Y, Z, coefficients, weights, offset, link, link.phi, type, control)

  ## further model information
  rval$call <- cl
  rval$formula <- oformula
  rval$terms <- list(mean = mtX, precision = mtZ, full = mt)
  rval$levels <- list(mean = .getXlevels(mtX, mf), precision = .getXlevels(mtZ, mf), full = .getXlevels(mt, mf))
  rval$contrasts <- list(mean = attr(X, "contrasts"), precision = attr(Z, "contrasts"))
  if(model) rval$model <- mf
  if(y) rval$y <- Y
  if(x) rval$x <- list(mean = X, precision = Z)

  class(rval) <- "betareg"
  return(rval)
}

betareg.nofit <- function(x, y, z = NULL, coefficients, weights = NULL, offset = NULL,
                          link = "logit", link.phi = "log", type = "ML", control = betareg.control())
{
  ## response and regressor matrix
  n <- NROW(x)
  k <- NCOL(x)
  if(is.null(weights)) weights <- rep.int(1, n)
  nobs <- sum(weights > 0)
  if(is.null(offset)) offset <- rep.int(0, n)
  if(!is.list(offset)) offset <- list(mean = offset, precision = rep.int(0, n))
  if(is.null(z)) {
    m <- 1L
    z <- matrix(1, ncol = m, nrow = n)
    colnames(z) <- "(Intercept)"
    rownames(z) <- rownames(x)
    phi_const <- TRUE
  } else {
    m <- NCOL(z)
    if(m < 1L) stop("dispersion regression needs to have at least one parameter")
    phi_const <- (m == 1L) && isTRUE(all.equal(as.vector(z[, 1L]), rep.int(1, n)))
  }

  ## link processing
  if(is.character(link)) {
    linkstr <- link
    if(linkstr != "loglog") {
      linkobj <- make.link(linkstr)
      ## add dmu.deta potentially needed for BC/BR
      make.dmu.deta <- getFromNamespace("make.dmu.deta", ns = "betareg")
      linkobj$dmu.deta <- make.dmu.deta(linkstr)
    } else {
      linkobj <- structure(list(
        linkfun = function(mu) -log(-log(mu)),
        linkinv = function(eta) pmax(pmin(exp(-exp(-eta)), 1 - .Machine$double.eps), .Machine$double.eps),
        mu.eta = function(eta) {
          eta <- pmin(eta, 700)
          pmax(exp(-eta - exp(-eta)), .Machine$double.eps)
        },
        dmu.deta = function(eta) pmax(exp(-exp(-eta) - eta) * expm1(-eta), .Machine$double.eps),
        valideta = function(eta) TRUE,
        name = "loglog"
      ), class = "link-glm")
    }
  } else {
    linkobj <- link
    linkstr <- link$name
    if(type != "ML" && is.null(linkobj$dmu.deta)) warning("link needs to provide dmu.deta component for BC/BR")
  }
  linkfun <- linkobj$linkfun
  linkinv <- linkobj$linkinv
  mu.eta <- linkobj$mu.eta
  dmu.deta <- linkobj$dmu.deta
  if(is.character(link.phi)) {
    phi_linkstr <- link.phi
    phi_linkobj <- make.link(phi_linkstr)
    phi_linkobj$dmu.deta <- make.dmu.deta(phi_linkstr)
  } else {
    phi_linkobj <- link.phi
    phi_linkstr <- link.phi$name
    if(type != "ML" && is.null(phi_linkobj$dmu.deta)) warning("link.phi needs to provide dmu.deta component for BC/BR")
  }
  phi_linkfun <- phi_linkobj$linkfun
  phi_linkinv <- phi_linkobj$linkinv
  phi_mu.eta <- phi_linkobj$mu.eta
  phi_dmu.deta <- phi_linkobj$dmu.deta
  ## y* transformation
  ystar <- qlogis(y)

  ## control parameters
  ocontrol <- control
  phi_full <- control$phi
  method <- control$method
  hessian <- control$hessian
  start <- control$start
  fsmaxit <- control$fsmaxit
  fstol <- control$fstol
  control$phi <- control$method <- control$hessian <- control$start <- control$fsmaxit <- control$fstol <- NULL

  ## starting values
  if(is.null(start)) {
    auxreg <- lm.wfit(x, linkfun(y), weights, offset = offset[[1L]])
    beta <- auxreg$coefficients
    yhat <- linkinv(auxreg$fitted.values)
    dlink <- 1/mu.eta(linkfun(yhat))
    res <- auxreg$residuals
    sigma2 <- sum(weights * res^2)/((sum(weights) - k) * (dlink)^2)
    phi_y <- weights * yhat * (1 - yhat)/(sum(weights) * sigma2) - 1/n
    phi <- rep.int(0, ncol(z))
    phi[1L] <- suppressWarnings(phi_linkfun(sum(phi_y)))
    ## i.e., start out from the fixed dispersion model as described
    ## in Ferrari & Cribari-Neto (2004) (and differing from Simas et al. 2009)
    ## An alternative would be
    ##   phi <- lm.wfit(z, phi_linkfun(phi_y), weights)$coefficients
    ## but that only works in general if all(phi_y > 0) which is not necessarily
    ## the case.
    ##
    ## Additionally, sum(phi_y) might not even be > 0 which should be caught.
    if(!isTRUE(phi_linkinv(phi[1L]) > 0)) {
      warning("no valid starting value for precision parameter found, using 1 instead")
      phi[1L] <- 1
    }
    start <- list(mean = beta, precision = phi)
  }
  if(is.list(start)) start <- do.call("c", start)

  ## various fitted quantities (parameters, linear predictors, etc.)
  fitfun <- function(par, deriv = 0L) {
    beta <- par[seq.int(length.out = k)]
    gamma <- par[seq.int(length.out = m) + k]
    eta <- as.vector(x %*% beta + offset[[1L]])
    phi_eta <- as.vector(z %*% gamma + offset[[2L]])
    mu <- linkinv(eta)
    phi <- phi_linkinv(phi_eta)
    mustar <- if(deriv >= 1L) digamma(mu * phi) - digamma((1 - mu) * phi) else NULL
    psi1 <- if(deriv >= 2L) trigamma(mu * phi) else NULL
    psi2 <- if(deriv >= 2L) trigamma((1 - mu) * phi) else NULL
    list(
      beta = beta,
      gamma = gamma,
      eta = eta,
      phi_eta = phi_eta,
      mu = mu,
      phi = phi,
      mustar = mustar,
      psi1 = psi1,
      psi2 = psi2
    )
  }

  ## objective function
  loglikfun <- function(par, fit = NULL) {
    ## extract fitted parameters
    if(is.null(fit)) fit <- fitfun(par)
    alpha <- fit$mu * fit$phi
    beta <- (1 - fit$mu) * fit$phi

    ## compute log-likelihood
    if(any(!is.finite(fit$phi)) | any(alpha > 1e300) | any(beta > 1e300)) NaN else { ## catch extreme cases without warning
      ll <- suppressWarnings(dbeta(y, alpha, beta, log = TRUE))
      if(any(!is.finite(ll))) NaN else sum(weights * ll) ## again: catch extreme cases without warning
    }
  }

  ## gradient (by default) or gradient contributions (sum = FALSE)
  gradfun <- function(par, sum = TRUE, fit = NULL) {
    ## extract fitted means/precisions
    if(is.null(fit)) fit <- fitfun(par, deriv = 1L)
    mu <- fit$mu
    phi <- fit$phi
    eta <- fit$eta
    phi_eta <- fit$phi_eta
    mustar <- fit$mustar

    ## compute gradient contributions
    rval <- cbind(
      phi * (ystar - mustar) * mu.eta(eta) * weights * x,
      (mu * (ystar - mustar) + log(1-y) - digamma((1-mu)*phi) + digamma(phi)) *
        phi_mu.eta(phi_eta) * weights * z
    )
    if(sum) colSums(rval) else rval
  }

  ## analytical Hessian (expected information) or covariance matrix (inverse of Hessian)
  hessfun <- function(par, inverse = FALSE, fit = NULL) {
    ## extract fitted means/precisions
    if(is.null(fit)) fit <- fitfun(par, deriv = 2L)
    mu <- fit$mu
    phi <- fit$phi
    eta <- fit$eta
    phi_eta <- fit$phi_eta
    mustar <- fit$mustar
    psi1 <- fit$psi1
    psi2 <- fit$psi2

    ## auxiliary transformations
    a <- psi1 + psi2
    b <- psi1 * mu^2 + psi2 * (1-mu)^2 - trigamma(phi)
    ## compute elements of W
    wbb <- phi^2 * a * mu.eta(eta)^2
    wpp <- b * phi_mu.eta(phi_eta)^2
    wbp <- phi * (mu * a - psi2) * mu.eta(eta) * phi_mu.eta(phi_eta)
    ## compute elements of K
    kbb <- if(k > 0L) crossprod(sqrt(weights) * sqrt(wbb) * x) else crossprod(x)
    kpp <- if(m > 0L) crossprod(sqrt(weights) * sqrt(wpp) * z) else crossprod(z)
    kbp <- if(k > 0L & m > 0L) crossprod(weights * wbp * x, z) else crossprod(x, z)

    ## put together K (= expected information)
    K <- cbind(rbind(kbb, t(kbp)), rbind(kbp, kpp))
    if(!inverse) K else chol2inv(chol(K))
    ## previously computed K^(-1) via partitioned matrices, but this appears to be
    ## slower - even for moderately sized problems
    ##   kbb1 <- if(k > 0L) chol2inv(qr.R(qr(sqrt(weights) * sqrt(wbb) * x))) else kbb
    ##   kpp1 <- if(m > 0L) solve(kpp - t(kbp) %*% kbb1 %*% kbp) else kpp
    ##   vcov <- cbind(rbind(kbb1 + kbb1 %*% kbp %*% kpp1 %*% t(kbp) %*% kbb1,
    ##     -kpp1 %*% t(kbp) %*% kbb1), rbind(-kbb1 %*% kbp %*% kpp1, kpp1))
  }

  ## compute biases and adjustment for bias correction/reduction
  # biasfun <- function(par, fit = NULL, vcov = NULL) {
  #   if (is.null(fit)) fit <- fitfun(par, deriv = 2L)
  #   InfoInv <- if(is.null(vcov)) try(hessfun(par, inverse = TRUE), silent = TRUE) else vcov
  #   mu <- fit$mu
  #   phi <- fit$phi
  #   eta <- fit$eta
  #   phi_eta <- fit$phi_eta
  #   D1 <- mu.eta(eta)
  #   D2 <- phi_mu.eta(phi_eta)
  #   D1dash <- dmu.deta(eta)
  #   D2dash <- phi_dmu.deta(phi_eta)
  #   Psi2 <- fit$psi2
  #   dPsi1 <-  psigamma(mu * phi, 2)       ## potentially move to fitfun() when we add support for
  #   dPsi2 <-  psigamma((1 - mu) * phi, 2) ## observed information (as opposed to expected)
  #   kappa2 <- fit$psi1 + Psi2
  #   kappa3 <- dPsi1 - dPsi2
  #   Psi3 <- psigamma(phi, 1)
  #   dPsi3 <- psigamma(phi, 2)
  #   ## PQsum produces the the adustments to the score functions and is suggested for iteration
  #   PQsum <- function(t) {
  #     if (t <= k)  {
  #       Xt <- x[,t]
  #       bb <- if (k > 0L)
  #         crossprod(x, weights * phi^2 * D1 * (phi * D1^2 * kappa3 + D1dash * kappa2) * Xt * x)
  #       else
  #         crossprod(x)
  #       bg <- if ((k > 0L) & (m > 0L))
  #         crossprod(x, weights * phi * D1^2 * D2 * (mu * phi * kappa3 + phi * dPsi2 + kappa2) * Xt * z)
  #       else
  #         crossprod(x, z)
  #       gg <- if (m > 0L)
  #         crossprod(z, weights * phi * D1 * D2^2 * (mu^2 * kappa3 - dPsi2 + 2 * mu * dPsi2) * Xt * z) +
  #         crossprod(z, weights * phi * D1 * D2dash * (mu * kappa2 - Psi2) * Xt * z)
  #       else
  #         crossprod(z)
  #     } else {
  #       Zt <- z[, t - k]
  #       bb <- if (k > 0L)
  #         crossprod(x, weights * phi * D2 * (phi * D1^2 * mu * kappa3 + phi * D1^2 * dPsi2 + D1dash * mu * kappa2 - D1dash * Psi2) * Zt * x)
  #       else
  #         crossprod(x)
  #       bg <- if ((k > 0L) & (m > 0L))
  #         crossprod(x, weights * D1 * D2^2 * (phi * mu^2 * kappa3 + phi * (2 * mu - 1) * dPsi2 + mu * kappa2 - Psi2) * Zt * z)
  #       else
  #         crossprod(x, z)
  #       gg <- if (m > 0L)
  #         crossprod(z, weights * D2^3 * (mu^3 * kappa3 + (3 * mu^2 - 3 * mu + 1) * dPsi2 - dPsi3) * Zt * z) +
  #         crossprod(z, weights * D2dash * D2 * (mu^2 * kappa2 + (1 - 2 * mu) * Psi2 - Psi3) * Zt * z)
  #       else
  #         crossprod(z)
  #     }
  #     pq <- rbind(cbind(bb, bg), cbind(t(bg), gg))
  #     sum(diag(InfoInv %*% pq))/2
  #   }
  #   if (inherits(InfoInv, "try-error")) {
  #     bias <- adjustment <- rep.int(NA_real_, k + m)
  #   }
  #   else {
  #     adjustment <- sapply(1:(k + m), PQsum)
  #     bias <- - InfoInv %*% adjustment
  #   }
  #   list(bias = bias, adjustment = adjustment)
  # }
  #

  ## optimize likelihood
  # opt <- optim(par = start, fn = loglikfun, gr = gradfun,
  #              method = method, hessian = hessian, control = control)
  # par <- opt$par
  par <- coefficients

  ## conduct further (quasi) Fisher scoring to move ML derivatives
  ## even further to zero or conduct bias reduction
  ## (suppressed if fsmaxit = 0 or if only numerical optim result desired)
  # if(type == "BR" & fsmaxit <= 0) warning("BR cannot be performed with fsmaxit <= 0")
  # step <- .Machine$integer.max
  # iter <- 0
  # if(fsmaxit > 0 & !(hessian & type == "ML"))
  # {
  #   for (iter in 1:fsmaxit) {
  #     stepPrev <- step
  #     stepFactor <- 0
  #     testhalf <- TRUE
  #     while (testhalf & stepFactor < 11) {
  #       fit <- fitfun(par, deriv = 2L)
  #       scores <- gradfun(par, fit = fit)
  #       InfoInv <- try(hessfun(par, fit = fit, inverse = TRUE))
  #       if(failedInv <- inherits(InfoInv, "try-error")) {
  #         warning("failed to invert the information matrix: iteration stopped prematurely")
  #         break
  #       }
  #       bias <- if(type == "BR") biasfun(par, fit = fit, vcov = InfoInv)$bias else 0
  #       par <- par + 2^(-stepFactor) * (step <- InfoInv %*% scores - bias)
  #       stepFactor <- stepFactor + 1
  #       testhalf <- drop(crossprod(stepPrev) < crossprod(step))
  #     }
  #     if (failedInv | (all(abs(step) < fstol))) {
  #       break
  #     }
  #   }
  # }

  ## check whether both optim() and manual iteration converged IK:
  ## modified the condition a bit... optim might fail to converge but
  ## if additional iteration are requested Fisher scoring might get
  ## there
  # if((fsmaxit == 0 & opt$convergence > 0) | iter >= fsmaxit) {
  #   converged <- FALSE
  #   warning("optimization failed to converge")
  # } else {
  converged <- TRUE
  # }

  ## conduct single bias correction (if BC selected) else do not
  ## estimate the first order biases
  # if(type == "BC") {
  #   bias <- as.vector(biasfun(par)$bias)
  #   par <- par - bias
  # }
  # else {
  #   bias <- rep.int(NA_real_, k + m)
  # }

  bias <- rep.int(NA_real_, k + m)

  ## extract fitted values/parameters
  fit <- fitfun(par, deriv = 2L)
  beta <- fit$beta
  gamma <- fit$gamma
  eta <- fit$eta
  mu <- fit$mu
  phi <- fit$phi

  ## log-likelihood/gradients/covariance matrix at optimized parameters
  ll <- loglikfun(par, fit = fit)
  ## No need to evaluate ef below.
  ef <- gradfun(par, fit = fit, sum = FALSE)
  # vcov <- if (hessian & (type == "ML")) solve(-as.matrix(opt$hessian)) else hessfun(fit = fit, inverse = TRUE)
  vcov <- hessfun(fit = fit, inverse = TRUE)

  ## R-squared
  pseudor2 <- if(var(eta) * var(ystar) <= 0) NA else cor(eta, linkfun(y))^2

  ## names
  names(beta) <- colnames(x)
  names(gamma) <- if(phi_const & phi_linkstr == "identity") "(phi)" else colnames(z)
  rownames(vcov) <- colnames(vcov) <- names(bias) <- c(colnames(x),
                                                       if(phi_const & phi_linkstr == "identity") "(phi)" else paste("(phi)", colnames(z), sep = "_"))

  ## set up return value
  rval <- list(
    coefficients = list(mean = beta, precision = gamma),
    residuals = y - mu,
    fitted.values = structure(mu, .Names = names(y)),
    type = type,
    optim = NULL,
    method = method,
    control = ocontrol,
    scoring = NULL,
    start = start,
    weights = if(identical(as.vector(weights), rep.int(1, n))) NULL else weights,
    offset = list(mean = if(identical(offset[[1L]], rep.int(0, n))) NULL else offset[[1L]],
                  precision = if(identical(offset[[2L]], rep.int(0, n))) NULL else offset[[2L]]),
    n = n,
    nobs = nobs,
    df.null = nobs - 2,
    df.residual = nobs - k - m,
    phi = phi_full,
    loglik = ll,
    vcov = vcov,
    bias = bias,
    pseudo.r.squared = pseudor2,
    link = list(mean = linkobj, precision = phi_linkobj),
    converged = converged
  )
  return(rval)
}

#' @title
#' Construct a `glm` object for logistic regression with predefined coefficients
#' @description
#' It creates an object of class \code{glm} with \code{family=binomial(link = "logit")},
#' as one would obtain when fitting a logistic regression, but here there is no
#' fit, one simply provides a design matrix and coefficients. This function is
#' particularly useful in conjunction with \code{simulation} method from the
#' \code{ib} package to simulate from this model.
#' @param x n x p design matrix (no intercept)
#' @param coefficients p+1 vector of coefficients, the first coefficient is the intercept.
#' @importFrom stats rbinom
#' @seealso \code{\link[stats]{glm}}, \code{\link[ib]{simulation}}
#' @export
make_logistic <- function(x, coefficients){
  if(ncol(x) != length(coefficients) - 1L) stop("Design matrix mistmaches coefficients size")
  y_tmp <- rbinom(nrow(x),1,0.5)
  makeglm(y_tmp ~ x, family = binomial(link = "logit"), coefficients = coefficients)
}

#' @title
#' Construct a `negbin` object for negative binomial regression with predefined coefficients
#' @description
#' It creates an object of class \code{negbin} from \code{MASS} package,
#' as one would obtain when fitting a negative binomial regression, but here there is no
#' fit, one simply provides a design matrix and coefficients. This function is
#' particularly useful in conjunction with \code{simulation} method from the
#' \code{ib} package to simulate from this model.
#' @param x n x p design matrix (no intercept)
#' @param coefficients p+1 vector of coefficients, the first coefficient is the intercept.
#' @param theta additional parameter for negative binomial regression
#' @importFrom MASS rnegbin
#' @seealso \code{\link[MASS]{glm.nb}}, \code{\link[ib]{simulation}}
#' @export
make_negbin <- function(x, coefficients, theta){
  if(ncol(x) != length(coefficients) - 1L) stop("Design matrix mistmaches coefficients size")
  y_tmp <- rnegbin(nrow(x),1,1)
  makeglm.nb(y_tmp ~ x, coefficients = coefficients, init.theta = theta)
}

#' @title
#' Construct a `betareg` object for beta regression with predefined coefficients
#' @description
#' It creates an object of class \code{betareg} from \code{betareg} package,
#' as one would obtain when fitting a beta regression, but here there is no
#' fit, one simply provides a design matrix and coefficients. This function is
#' particularly useful in conjunction with \code{simulation} method from the
#' \code{ib} package to simulate from this model.
#' @param x n x p design matrix (no intercept)
#' @param coefficients p+2 vector of coefficients, the first coefficient is the intercept,
#' the last coefficient is the precision parameter.
#' @importFrom stats rbeta
#' @seealso \code{\link[betareg]{betareg}}, \code{\link[ib]{simulation}}
#' @export
make_betareg <- function(x, coefficients){
  if(ncol(x) != length(coefficients) - 2L) stop("Design matrix mistmaches coefficients size")
  y_tmp <- rbeta(nrow(x),1,1)
  makebetareg(y_tmp ~ x, link.phi = "log", coefficients = coefficients)
}
