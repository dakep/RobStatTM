#' M-scale estimator
#'
#' This function computes an M-scale, which is a robust
#' scale (spread) estimator.
#' M-estimators of scale are a robust alternative to
#' the sample standard deviation. Given a vector of
#' residuals \code{r}, the M-scale estimator \code{s}
#' solves the non-linear equation \code{mean(rho(r/s, cc))=b},
#' where \code{b} and \code{cc} are user-chosen tuning constants.
#' In this package the function \code{rho} is one of
#' Tukey's bisquare family.
#' The breakdown point of the estimator is \code{min(b, 1-b)},
#' so the optimal choice for \code{b} is 0.5. To obtain a
#' consistent estimator the constant
#' \code{cc} should be chosen such that E(rho(Z, cc)) = b, where
#' Z is a standard normal random variable.
#'
#' The iterative algorithm starts from the scaled median of
#' the absolute values of the input vector, and then
#' cycles through the equation s^2 = s^2 * mean(rho(r/s, cc)) / b.
#'
#' @export mscale scaleM
#' @aliases mscale scaleM
#' @rdname mscale
#'
#' @param u vector of residuals
#' @param delta the right hand side of the M-scale equation
#' @param family string specifying the name of the family of loss function to be used (current valid
#' options are "bisquare", "opt" and "mopt").
#' @param tuning.chi the tuning object for the rho function as returned
#' by \code{\link{lmrobdet.control}}, \link{bisquare}, \link{mopt} or \link{opt}.
#' It should correspond to the family of rho functions specified in the argument \code{family}.
#' @param tol relative tolerance for convergence
#' @param max.it maximum number of iterations allowed
#' @param tolerancezero smallest (in absolute value) non-zero value accepted as a scale. Defaults to \code{.Machine$double.eps}
#'
#' @return The scale estimate value at the last iteration or at convergence.
#'
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}
#'
#' @examples
#' set.seed(123)
#' r <- rnorm(150, sd=1.5)
#' mscale(r)
#' sd(r)
#' # 10% of outliers, sd of good points is 1.5
#' set.seed(123)
#' r2 <- c(rnorm(135, sd=1.5), rnorm(15, mean=-5, sd=.5))
#' mscale(r2)
#' sd(r2)
#'
scaleM <- mscale <- function(u, delta=0.5, tuning.chi=1.547645, family ="bisquare",
                             max.it=100, tol=1e-6, tolerancezero=.Machine$double.eps) {
  # M-scale of a sample u
  # tol: accuracy
  # delta: breakdown point (right side)
  # Initial
  s0 <- median(abs(u))/.6745
  if(s0 < tolerancezero) return(0)
  err <- tol + 1
  it <- 0
  while( (err > tol) && ( it < max.it) ) {
    it <- it+1
    s1 <- sqrt( s0^2 * mean(rho(u/s0, family = family, cc = tuning.chi)) / delta )
    err <- abs(s1-s0)/s0
    s0 <- s1
  }
  return(s0)
}

#' MM regression estimator using Pen~a-Yohai candidates
#'
#' This function computes MM-regression estimator using Pen~a-Yohai
#' candidates for the initial S-estimator. This function is used
#' internally by \code{\link{lmrobdetMM}}, and not meant to be used
#' directly.
#'
#' @param X design matrix
#' @param y response vector
#' @param control a list of control parameters as returned by \code{\link{lmrobdet.control}}
#' @param mf model frame
#'
#' @return an \code{\link{lmrob}} object witht the M-estimator
#' obtained starting from the S-estimator computed with the
#' Pen~a-Yohai initial candidates. The properties of the final
#' estimator (efficiency, etc.) are determined by the tuning constants in
#' the argument \code{control}.
#'
#' @rdname MMPY
#' @author Victor Yohai, \email{victoryohai@gmail.com}, Matias Salibian-Barrera, \email{matias@stat.ubc.ca}
#' @references \url{http://www.wiley.com/go/maronna/robust}
#' @seealso \code{\link{DCML}}, \code{\link{MMPY}}, \code{\link{SMPY}}
#'
#' @export
MMPY <- function(X, y, control, mf) {
  # INPUT
  # X nxp matrix, where n is the number of observations and p the number of  columns
  # y vector of dimension  n with the responses
  #
  # OUTPUT
  # outMM output of the MM estimator (lmrob) with 85% of efficiency and PY as initial
  n <- nrow(X)
  p <- ncol(X)
  dee <- control$bb
  if(control$corr.b) dee <- dee * (1-(p/n))
  a <- pyinit::pyinit(x=X, y=y, intercept=FALSE, delta=dee,
                      cc=1.54764,
                      psc_keep=control$psc_keep*(1-(p/n)), resid_keep_method=control$resid_keep_method,
                      resid_keep_thresh = control$resid_keep_thresh, resid_keep_prop=control$resid_keep_prop,
                      maxit = control$py_maxit, eps=control$py_eps,
                      mscale_maxit = control$mscale_maxit, mscale_tol = control$mscale_tol,
                      mscale_rho_fun="bisquare")
  # refine the PY candidates to get something closer to an S-estimator for y ~ X1
  kk <- dim(a$coefficients)[2]
  best.ss <- +Inf

  for(i in 1:kk) {
    tmp <- refine.sm(x=X, y=y, initial.beta=a$coefficients[,i], initial.scale=a$objective[i],
                     k=control$refine.PY, conv=1, b=dee, family=control$family, cc=control$tuning.chi, step='S')
    if(tmp$scale.rw < best.ss) {
      best.ss <- tmp$scale.rw # initial$objF[1]
      betapy <- tmp$beta.rw # initial$initCoef[,1]
    }
  }
  S.init <- list(coef=betapy, scale=best.ss)
  orig.control <- control

  control$psi <- control$family # tuning.psi$name
  # control$tuning.psi <- control$tuning.psi # $cc

  control$method <- 'M'
  control$cov <- ".vcov.w"
  control$subsampling <- 'simple'

  # # lmrob() does the above when is.list(init)==TRUE, in particular:
  outMM <- lmrob.fit(X, y, control, init=S.init, mf=mf)
  outMM$control <- orig.control
  outMM$init <- S.init
  coefnames <- names(coef(outMM))
  residnames <- names(resid(outMM))
  # if weights were all zero, then return the S estimator
  # S.resid <- as.vector(y - X %*% S.init$coef) / S.init$scale
  # ws <- rhoprime(S.resid, family=orig.control$family, cc=orig.control$tuning.psi)
  if(all( outMM$rweights == 0 ) ) {
    outMM$coefficients <- as.vector(S.init$coef)
    outMM$scale <- S.init$scale
    names(outMM$coefficients) <- coefnames
    outMM$residuals <- as.vector(y - X %*% S.init$coef)
    names(outMM$residuals) <- residnames
  }
  return(outMM)
}
