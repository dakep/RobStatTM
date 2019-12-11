lmrob.fit <- function(x, y, control, init=NULL, mf=NULL) {
  if(!is.matrix(x)) x <- as.matrix(x)
  ## old notation: MM -> SM
  if (control$method == "MM") control$method <- "SM"
  ## Assumption:  if(is.null(init))  method = "S..."   else  method = "..."
  ## ---------    where "..." consists of letters {"M", "D"}
  est <- if (is.null(init)) {
    ## --- initial S estimator
    if ((M1 <- substr(control$method,1,1)) != 'S') {
      warning(gettextf("Initial estimator '%s' not supported; using S-estimator instead",
                       M1), domain = NA)
      substr(control$method,1,1) <- 'S'
    }
    init <- lmrob.S(x, y, control = control, mf = mf)
    'S'
  } else {
    stopifnot(is.list(init))
    if (is.null(init$converged)) init$converged <- TRUE
    if (is.null(init$control)) {
      init$control <- control
      M <- init$control$method <- 'l'
    } else if(!length(M <- init$control$method) || !nzchar(M))
      M <- "l"
    M
  }
  stopifnot(is.numeric(init$coef), length(init$coef) == ncol(x),
            is.numeric(init$scale), init$scale >= 0)

  trace.lev <- control$trace.lev
  if (init$converged) {
    ## --- loop through the other estimators; build up 'est' string
    method <- sub(paste0("^", est), '', control$method)
    if(trace.lev) {
      cat(sprintf("init converged (remaining method = \"%s\") -> coef=\n", method))
      print(init$coef) }
    for (step in strsplit(method,'')[[1]]) {
      ## now we have either M or D steps
      est <- paste0(est, step)
      init <- switch(step, ## 'control' may differ from 'init$control' when both (init, control) are spec.
                     ## D(AS)-Step
                     D = lmrob..D..fit(init, x, mf = mf,
                                       control=control, method = init$control$method),
                     ## M-Step
                     M = lmrob..M..fit(x = x, y = y, obj = init, mf = mf,
                                       control=control, method = init$control$method),
                     stop('only M and D are steps supported after "init" computation'))
      if(trace.lev) { cat(sprintf("step \"%s\" -> new coef=\n", step)); print(init$coef) }
      ## break if an estimator did not converge
      if (!init$converged) {
        warning(gettextf(
          "%s-step did NOT converge. Returning unconverged %s-estimate",
          step, est),
          domain = NA)
        # break # now we return unconverged estimator with all associated inference
      }
    }
  }
  ## << FIXME? qr(.)  should be available from earlier
  if (is.null(init$qr)) init$qr <- qr(x * sqrt(init$rweights))
  if (is.null(init$rank)) init$rank <- init$qr$rank
  control$method <- est ## ~= original 'method', but only with the steps executed.
  init$control <- control#c
  df <- NROW(y) - init$rank ## sum(init$r?weights)-init$rank
  init$degree.freedom <- init$df.residual <- df
  init
}

lmrob..M..fit <- function (x = obj$x, y = obj$y, beta.initial = obj$coef,
                           scale = obj$scale, control = obj$control,
                           obj,
                           mf = obj$model,
                           ##   ^^^^^^^^^ not model.frame(obj) to avoid errors.
                           method = obj$control$method) #<- also when 'control' is not obj$control
{
  c.psi <- .psi.conv.cc(control$psi, control$tuning.psi)
  ipsi <- .psi2ipsi(control$psi)
  stopifnot(is.matrix(x))
  n <- nrow(x)
  p <- ncol(x)
  if (is.null(y) && !is.null(obj$model))
    y <- model.response(obj$model, "numeric")
  # Only optimal and modified.optimal have more than 4 tuning constants
  stopifnot(length(y) == n,
            length(c.psi) > 0, c.psi[-5] >= 0,
            scale >= 0, length(beta.initial) == p)

  ret <- .C(R_lmrob_MM,
            x = as.double(x),
            y = as.double(y),
            n = as.integer(n),
            p = as.integer(p),
            beta.initial = as.double(beta.initial),
            scale = as.double(scale),
            coefficients = double(p),
            residuals = double(n),
            iter = as.integer(control$max.it),
            c.psi = as.double(c.psi),
            ipsi = as.integer(ipsi),
            loss = double(1),
            rel.tol = as.double(control$rel.tol),
            converged = logical(1),
            trace.lev = as.integer(control$trace.lev),
            mts = as.integer(control$mts),
            ss =  .convSs(control$subsampling)
  )[c("coefficients",  "scale", "residuals", "loss", "converged", "iter")]
  ## FIXME?: Should rather warn *here* in case of non-convergence
  ret$fitted.values <- drop(x %*% ret$coefficients)
  names(ret$coefficients) <- colnames(x)
  names(ret$residuals) <- rownames(x)
  ret$rweights <- lmrob.rweights(ret$residuals, scale, control$tuning.psi, control$psi)
  ret$control <- control
  if (!missing(obj)) {
    if (!grepl('M$', method)) {
      ## update method if it's not there already
      method <- paste0(method, 'M')
    }
    if (!is.null(obj$call)) {
      ret$call <- obj$call
      ret$call$method <- method
    }
    if (method %in% c('SM', 'MM')) {
      ret$init.S <- obj
    } else {
      ret$init <-
        obj[intersect(names(obj),
                      c("coefficients", "scale", "residuals", "loss", "converged",
                        "iter", "rweights", "fitted.values", "control", "ostats",
                        "init.S", "init", "kappa", "tau"))]
      class(ret$init) <- 'lmrob'
      ret <- c(ret,
               obj[intersect(names(obj),
                             c("df.residual", "degree.freedom",
                               "xlevels", "terms", "model", "x", "y",
                               "na.action", "contrasts", "MD"))])
    }
    ret$qr <- qr(x * sqrt(ret$rweights))
    ret$rank <- ret$qr$rank
    if (!is.null(obj$assign)) ret$assign <- obj$assign
    if (method %in% control$compute.outlier.stats)
      ret$ostats <- outlierStats(ret, x, control)
  }
  class(ret) <- "lmrob"
  ret
}


lmrob.S <- function (x, y, control, trace.lev = control$trace.lev, mf = NULL)
{
  if (!is.matrix(x)) x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  nResample <- as.integer(control$nResample)
  groups <- as.integer(control$groups)
  nGr <- as.integer(control$n.group)
  large_n <- (n > control$fast.s.large.n)
  if (large_n) {
    if (nGr <= p)
      stop("'control$n.group' must be larger than 'p' for 'large_n' algorithm")
    if (nGr * groups > n)
      stop("'groups * n.group' must be smaller than 'n' for 'large_n' algorithm")
    if (nGr <= p + 10) ## FIXME (be smarter ..)
      warning("'control$n.group' is not much larger than 'p', probably too small")
  }
  if (length(seed <- control$seed) > 0) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      seed.keep <- get(".Random.seed", envir = .GlobalEnv,
                       inherits = FALSE)
      on.exit(assign(".Random.seed", seed.keep, envir = .GlobalEnv))
    }
    assign(".Random.seed", seed, envir = .GlobalEnv) ## why not set.seed(seed)
  }

  bb <- as.double(control$bb)
  c.chi <- .psi.conv.cc(control$psi, control$tuning.chi)
  best.r <- as.integer(control$best.r.s)
  #Only optimal and modified.optimal have more than 4 tuning constants
  stopifnot(length(c.chi) > 0, c.chi[-5] >= 0, length(bb) > 0,
            length(best.r) > 0, best.r >= 1, length(y) == n, n > 0)

  b <- .C(R_lmrob_S,
          x = as.double(x),
          y = as.double(y),
          n = as.integer(n),
          p = as.integer(p),
          nResample = nResample,
          scale = double(1),
          coefficients = double(p),
          as.double(c.chi),
          .psi2ipsi(control$psi),
          bb,
          best_r = best.r,
          groups = groups,
          n.group = nGr,
          k.fast.s = as.integer(control$k.fast.s),
          k.iter = as.integer(control$k.max),
          maxit.scale = as.integer(control$maxit.scale),
          refine.tol = as.double(control$refine.tol),
          inv.tol = as.double(control$solve.tol),
          converged = logical(1),
          trace.lev = as.integer(trace.lev),
          mts = as.integer(control$mts),
          ss = .convSs(control$subsampling),
          fast.s.large.n = as.integer(if (large_n) control$fast.s.large.n else n+1)
          ## avoids the use of NAOK = TRUE for control$fast.s.large.n == Inf
  )[c("coefficients", "scale", "k.iter", "converged")]
  scale <- b$scale
  if (scale < 0)
    stop("C function R_lmrob_S() exited prematurely")
  if (scale == 0)
    warning("S-estimated scale == 0:  Probably exact fit; check your data")
  ## FIXME: get 'res'iduals from C
  if(trace.lev) {
    cat(sprintf("lmrob.S(): scale = %g; coeff.=\n", scale)); print(b$coefficients) }
  b$fitted.values <- x %*% b$coefficients
  b$residuals <- setNames(drop(y - b$fitted.values), rownames(x))
  names(b$coefficients) <- colnames(x)
  ## robustness weights
  b$rweights <- lmrob.rweights(b$residuals, scale, control$tuning.chi, control$psi)
  ## set method argument in control
  control$method <- 'S'
  b$control <- control
  ## add call if called from toplevel
  if (identical(parent.frame(), .GlobalEnv))
    b$call <- match.call()
  class(b) <- 'lmrob.S'
  if ("S" %in% control$compute.outlier.stats)
    b$ostats <- outlierStats(b, x, control)
  b
}



.Mpsi.R.names <- c('bisquare', 'lqq', 'welsh', 'opt', 'hampel', 'ggw', 'mopt')

.Mpsi.M.names <- c('huber') ## .M: the monotone ones:

.Mpsi.names <- c(R= .Mpsi.R.names, M= .Mpsi.M.names)

.regularize.Mpsi <- function(psi, redescending = TRUE) {
    stopifnot(is.character(psi), length(psi) == 1)
    psi <- tolower(psi)
    psi <- switch(psi,
		  'tukey'= , 'biweight'= "bisquare",
		  ## otherwise keep
		  psi)
    nms <- if(redescending) .Mpsi.R.names else .Mpsi.names
    if (is.na(i <- pmatch(psi, nms)))
	stop(gettextf("'psi' should be one of %s", paste(dQuote(nms), collapse=", ")),
	     domain = NA)
    nms[i]
}

.Mpsi.tuning.defaults <- list(
    'huber' = 1.345
    , 'bisquare' = 4.685061
    , 'welsh' = 2.11
    , 'ggw' = c(-0.5, 1.5, .95, NA) ## (min{slope}, b ,  eff, bp)
    , 'lqq' = c(-0.5, 1.5, .95, NA) ## (min{slope}, b/c, eff, bp)
    , 'opt' = c(a = 0.01317965,
                    lower = 0.03305454,
                    upper = 3.003281,
                    c = 1.0,
                    "Psi_Opt(lower)" = -0.0005459033,
                    "rho(Inf)" = 3.331370)
    , 'hampel' = c(1.5, 3.5, 8) * 0.9016085 ## a, b, r
    , 'mopt' = c(a = 0.01316352,
                             normConst = 1.05753107,
                             upper = 3.00373940,
                             c = 1.0,
                             "Psi_Opt(1)" = 0.46057111,
                             "rho(Inf)" = 3.53690811)
    )

.Mpsi.tuning.default <- function(psi) {
    if(is.null(p <- .Mpsi.tuning.defaults[[psi]]))
        stop(gettextf("invalid 'psi'=%s; possibly use .regularize.Mpsi(%s)",
                      psi, "psi, redescending=FALSE"), domain=NA)
    p
}

.Mchi.tuning.defaults <- list(
    ## Here, psi must be redescending! -> 'huber' not possible
    'bisquare' = 1.54764
    , 'welsh' = 0.5773502
    , 'ggw' = c(-0.5, 1.5, NA, .50) ## (min{slope}, b ,  eff, bp)
    , 'lqq' = c(-0.5, 1.5, NA, .50) ## (min{slope}, b/c, eff, bp)
    , 'opt' = c(a = 0.01317965,
                    lower = 0.03305454,
                    upper = 3.003281,
                    c = 0.2618571,
                    "Psi_Opt(lower)" = -0.0005459033,
                    "rho(Inf)" = 3.331370)
    , 'hampel' = c(1.5, 3.5, 8) * 0.2119163 ## a, b, r
    , 'mopt' = c(a = 0.01316352,
                             normConst = 1.05753107,
                             upper = 3.00373940,
                             c = 0.38124404,
                             "Psi_Opt(1)" = 0.46057111,
                             "rho(Inf)" = 3.53690811)
    )
.Mchi.tuning.default <- function(psi) {
    if(is.null(p <- .Mchi.tuning.defaults[[psi]]))
	stop(gettextf("invalid 'psi'=%s; possibly use .regularize.Mpsi(%s)",
		      psi, "psi"), domain=NA)
    p
}


lmrob.control <-
    function(setting, seed = NULL, nResample = 500,
	     tuning.chi = NULL, bb = 0.5,
	     tuning.psi = NULL, max.it = 50,
	     groups = 5, n.group = 400, k.fast.s = 1, best.r.s = 2,
	     k.max = 200, maxit.scale = 200, k.m_s = 20,
	     ##           ^^^^^^^^^^^ had MAX_ITER_FIND_SCALE 200 in ../src/lmrob.c
	     refine.tol = 1e-7, rel.tol = 1e-7,
	     solve.tol = 1e-7,
	     ## had  ^^^^^^^^  TOL_INVERSE 1e-7 in ../src/lmrob.c
	     trace.lev = 0, mts = 1000,
	     subsampling = c("nonsingular", "simple"),
	     compute.rd = FALSE,
	     method = 'MM',
	     psi = 'bisquare',
	     numpoints = 10, cov = NULL,
	     split.type = c("f", "fi", "fii"),
	     fast.s.large.n = 2000,
             eps.outlier = function(nobs) 0.1 / nobs,
             eps.x = function(maxx) .Machine$double.eps^(.75)*maxx,
             compute.outlier.stats = method,
             warn.limit.reject = 0.5,
             warn.limit.meanrw = 0.5,
             ...)
{
    p.ok <- missing(psi) # if(p.ok) psi does not need regularization
    if (!missing(setting)) {
        if (setting %in% c('KS2011', 'KS2014')) {
            if (missing(method)) method <- 'SMDM'
	    psi <- if(p.ok) 'lqq' else .regularize.Mpsi(psi) ; p.ok <- TRUE
            if (missing(max.it)) max.it <- 500
            if (missing(k.max)) k.max <- 2000
            if (setting == 'KS2014') {
                if (missing(best.r.s)) best.r.s <- 20
                if (missing(k.fast.s)) k.fast.s <- 2
                if (missing(nResample)) nResample <- 1000
            }
        } else {
            warning("Unknown setting '", setting, "'. Using defaults.")
        }
    } else {
	if(p.ok && grepl('D', method)) psi <- 'lqq'
    }
    if(!p.ok) psi <- .regularize.Mpsi(psi)
    subsampling <- match.arg(subsampling)

    ## in ggw, lqq:  if tuning.{psi|chi}  are non-standard, calculate coefficients:
    compute.const <- (psi %in% c('ggw', 'lqq'))

    if(is.null(tuning.chi))
	tuning.chi <- .Mchi.tuning.default(psi)
    else if(compute.const)
	tuning.chi <- .psi.const(tuning.chi, psi)

    if(is.null(tuning.psi))
	tuning.psi <- .Mpsi.tuning.default(psi)
    else if(compute.const)
	tuning.psi <- .psi.const(tuning.psi, psi)

    c(list(setting = if (missing(setting)) NULL else setting,
           seed = as.integer(seed), nResample=nResample, psi=psi,
           tuning.chi=tuning.chi, bb=bb, tuning.psi=tuning.psi,
           max.it=max.it, groups=groups, n.group=n.group,
           best.r.s=best.r.s, k.fast.s=k.fast.s,
           k.max=k.max, maxit.scale=maxit.scale, k.m_s=k.m_s, refine.tol=refine.tol,
           rel.tol=rel.tol, solve.tol=solve.tol, trace.lev=trace.lev, mts=mts,
           subsampling=subsampling,
           compute.rd=compute.rd, method=method, numpoints=numpoints,
           split.type = match.arg(split.type),
           fast.s.large.n=fast.s.large.n,
           eps.outlier = eps.outlier, eps.x = eps.x,
           compute.outlier.stats = sub("^MM$", "SM", compute.outlier.stats),
           warn.limit.reject = warn.limit.reject,
           warn.limit.meanrw = warn.limit.meanrw),
      list(...))
}

lmrob.control.neededOnly <- function(control) {
    if(is.null(control)) return(control)
    switch(sub("^(S|M-S).*", "\\1", control$method),
	   S = {                       # remove all M-S specific control pars
	       control$k.m_s <- NULL
	       control$split.type <- NULL
					# if large_n is not used, remove corresp control pars
	       if (length(residuals) <= control$fast.s.large.n) {
		   control$groups <- NULL
		   control$n.group <- NULL
	       }
	   },
	   `M-S` = {                # remove all fast S specific control pars
	       control$refine.tol <- NULL
	       control$groups <- NULL
	       control$n.group <- NULL
	       control$best.r.s <- NULL
	       control$k.fast.s <- NULL
	   }, { # else: do not keep parameters used by initial ests. only
	       control$tuning.chi <- NULL
	       control$bb <- NULL
	       control$refine.tol <- NULL
	       control$nResample <- NULL
	       control$groups <- NULL
	       control$n.group <- NULL
	       control$best.r.s <- NULL
	       control$k.fast.s <- NULL
	       control$k.max <- NULL
	       control$k.m_s <- NULL
	       control$split.type <- NULL
	       control$mts <- NULL
	       control$subsampling <- NULL
	   } )
    if (!grepl("D", control$method)) control$numpoints <- NULL
    if (control$method == 'SM') control$method <- 'MM'
    control
}

lmrob.fit.MM <- function(x, y, control) ## defunct
    .Defunct("lmrob.fit(*, control) with control$method = 'SM'")
## .Deprecated() till robustbase 0.92-6 (2016-05-28)


.psi2ipsi <- function(psi)
{
    psi <- .regularize.Mpsi(psi, redescending=FALSE)
    i <- match(psi, c(
	'huber', 'bisquare', 'welsh', 'opt',
	## 0	    1	        2	 3
	'hampel', 'ggw', 'lqq', 'mopt'
	## 4	    5	   6      7
	))
    if(is.na(i)) stop("internal logic error in psi() function name: ", psi,
		      "  Please report!")
    i - 1L
}

.psi.conv.cc <- function(psi, cc)
{
    if (!is.character(psi) || length(psi) != 1)
        stop("argument 'psi' must be a string (denoting a psi function)")
    if(!is.numeric(cc))
        stop("tuning constant 'cc' is not numeric")

    ## "FIXME": For (ggw, lqq) this is much related to  .psi.const() below
    switch(tolower(psi),
           'ggw' = {
               ## Input: 4 parameters, (minimal slope, b, efficiency, breakdown point) _or_ c(0, a,b,c, m.rho)
               ## Output 'k': either k in {1:6} or  k = c(0, k[2:5])

               ## prespecified 6 cases all treated in C ( ../src/lmrob.c ) :
               if (     isTRUE(all.equal(cc, c(-.5, 1  , 0.95, NA)))) return(1)
               else if (isTRUE(all.equal(cc, c(-.5, 1  , 0.85, NA)))) return(2)
               else if (isTRUE(all.equal(cc, c(-.5, 1. , NA, 0.5)))) return(3)
               else if (isTRUE(all.equal(cc, c(-.5, 1.5, 0.95, NA)))) return(4)
               else if (isTRUE(all.equal(cc, c(-.5, 1.5, 0.85, NA)))) return(5)
               else if (isTRUE(all.equal(cc, c(-.5, 1.5, NA, 0.5)))) return(6)
               else if (length(cc) == 5 && cc[1] == 0 ||
                        (length(cc <- attr(cc, 'constants')) == 5 && cc[1] == 0))
                   return(cc)
               else stop('Coefficients for ',psi,' function incorrectly specified.\n',
			 'Use c(minimal slope, b, efficiency, breakdown point) or c(0, a,b,c, max_rho).')
           },
           'lqq' = {
               ## Input: 4 parameters, (minimal slope, b/c, efficiency, breakdown point) _or_ (b, c, s) [length 3]
               ## Output: k[1:3] = (b, c, s)
               if (isTRUE(all.equal(cc, c(-.5, 1.5, 0.95, NA))))
                   return(c(1.4734061, 0.9822707, 1.5))
               else if (isTRUE(all.equal(cc, c(-.5, 1.5, NA, 0.5))))
                   return(c(0.4015457, 0.2676971, 1.5))
               else if (length(cc) == 3 || length(cc <- attr(cc, 'constants')) == 3)
                   return(cc)
               else stop('Coefficients for ',psi,' function incorrectly specified.\n',
                         'Use c(minimal slope, b, efficiency, breakdown point) [2 cases only] or  c(b, c, s)')
           },
           'hampel' = {
               ## just check length of coefficients
               if (length(cc) != 3)
                   stop('Coef. for Hampel psi function not of length 3')
           },
           'opt' = {
               ## just check length of coefficients
               if (length(cc) != 6)
                   stop('Coef. for Optimal psi function not of length 6')
           },
           'mopt' = {
               ## just check length of coefficients
               if (length(cc) != 6)
                   stop('Coef. for Modified Optimal psi function not of length 6')
           }, {
               ## otherwise: should have length 1
               if (length(cc) != 1)
                   stop('Coef. for psi function ', psi,' not of length 1')
           })

    return(cc)
}

.psi.ggw.mxs <- function(a, b, c, tol = .Machine$double.eps^0.25) {
    ipsi <- .psi2ipsi('ggw')
    ccc <- c(0, a, b, c, 1) ## == .psi.conv.cc('ggw', cc=c(0, a, b, c, 1))
    optimize(.Mpsi, c(c, max(a+b+2*c, 0.5)), ccc=ccc, ipsi=ipsi, deriv = 1,
             tol = tol)
}

.psi.ggw.ms <- function(a, b, c, tol = .Machine$double.eps^0.25) ## find minimal slope
    .psi.ggw.mxs(a, b, c, tol=tol)[["objective"]]

.psi.ggw.finda <- function(ms, b, c, tol = .Machine$double.eps^0.25, maxiter = 1000,
                            ms.tol = tol / 64,...) ## find constant 'a' (reparametrized to 1/o scale).
{
    val <- uniroot(function(a) .psi.ggw.ms(1/a, b, c, tol=ms.tol) - ms,
                   c(200, if (b > 1.4) 1/400 else if (b > 1.3) 1/50 else 1/20),
                   tol=tol, maxiter=maxiter)
    1/val$root
}

.psi.ggw.eff <- function(a, b, c) ## calculate asymptotic efficiency
{
    ipsi <- .psi2ipsi('ggw')
    ccc <- c(0, a, b, c, 1)
    lmrob.E(.Mpsi(r, ccc, ipsi, deriv=1), use.integrate = TRUE)^2 /
    lmrob.E(.Mpsi(r, ccc, ipsi) ^2,       use.integrate = TRUE)
}

.psi.ggw.bp <- function(a, b, c, ...) { ## calculate kappa
    ipsi <- .psi2ipsi('ggw')
    abc <- c(0, a, b, c)
    nc <- integrate(.Mpsi, 0, Inf, ccc = c(abc, 1), ipsi=ipsi, ...)$value
    lmrob.E(.Mchi(r, ccc = c(abc, nc), ipsi), use.integrate = TRUE)
}

.psi.ggw.findc <- function(ms, b, eff = NA, bp = NA,
                           subdivisions = 100L,
                           rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol,
                           tol = .Machine$double.eps^0.25, ms.tol = tol/64, maxiter = 1000) {
    ## find c by eff for bp
    c. <- if (!is.na(eff)) {
        if (!is.na(bp))
            warning('tuning constants for ggw psi: both eff and bp specified, ignoring bp')
        ## find c by eff
        tryCatch(uniroot(function(x) .psi.ggw.eff(.psi.ggw.finda(ms, b, x, ms.tol=ms.tol),
                                                  b, x) - eff,
                         c(0.15, if (b > 1.61) 1.4 else 1.9), tol=tol, maxiter=maxiter)$root,
		   error=function(e)e)
    } else {
        if (is.na(bp))
	    stop("neither breakdown point 'bp' nor efficiency 'eff' specified")
        ## find c by bp
        tryCatch(uniroot(function(x) .psi.ggw.bp(.psi.ggw.finda(ms, b, x, ms.tol=ms.tol),
                                                 b, x) - bp,
                c(0.08, if (ms < -0.4) 0.6 else 0.4), tol=tol, maxiter=maxiter)$root,
		   error=function(e)e)
    }
    if (inherits(c., 'error'))
        stop(gettextf('unable to find constants for "ggw" psi function: %s',
                      c.$message), domain=NA)
    a <- .psi.ggw.finda(ms, b, c., ms.tol=ms.tol)
    nc <- integrate(.Mpsi, 0, Inf, ccc= c(0, a, b, c., 1), ipsi = .psi2ipsi('ggw'))$value
    ## return
    c(0, a, b, c., nc)
}

lmrob.efficiency <-  function(psi, cc, ...) {
    ipsi <- .psi2ipsi(psi)
    ccc <- .psi.conv.cc(psi, cc=cc)

    integrate(function(x) .Mpsi(x, ccc=ccc, ipsi=ipsi, deriv=1)*dnorm(x),
	      -Inf, Inf, ...)$value^2 /
    integrate(function(x) .Mpsi(x, ccc=ccc, ipsi=ipsi)^2 *dnorm(x),
	      -Inf, Inf, ...)$value
}

lmrob.bp <- function(psi, cc, ...)
  integrate(function(x) Mchi(x, cc, psi)*dnorm(x), -Inf, Inf, ...)$value

.psi.lqq.findc <-
    function(ms, b.c, eff = NA, bp = NA,
             interval = c(0.1, 4), subdivisions = 100L,
             rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol,
             tol = .Machine$double.eps^0.25, maxiter = 1000)
{
    ## b.c == b/c
    bcs <- function(cc) c(b.c*cc, cc, 1-ms)
    t.fun <- if (!is.na(eff)) { ## 'eff' specified
	if (!is.na(bp))
	    warning("tuning constants for \"lqq\" psi: both 'eff' and 'bp' specified, ignoring 'bp'")
	## find c by b, s and eff
	function(c)
	    lmrob.efficiency('lqq', bcs(c), subdivisions=subdivisions,
                             rel.tol=rel.tol, abs.tol=abs.tol) - eff
    } else {
	if (is.na(bp))
	    stop('Error: neither breakdown point nor efficiency specified')
        ## breakdown point  'bp' specified
	function(c)
	    lmrob.bp('lqq', bcs(c), subdivisions=subdivisions,
                     rel.tol=rel.tol, abs.tol=abs.tol) - bp
    }
    c. <- tryCatch(uniroot(t.fun, interval=interval, tol=tol, maxiter=maxiter)$root,
		   error=function(e)e)
    if (inherits(c., 'error'))
        stop(gettextf('unable to find constants for "lqq" psi function: %s',
                      c.$message), domain=NA)
    else bcs(c.)
}

.psi.const <- function(cc, psi)
{
    switch(psi,
           "ggw" = { ## only calculate for non-standard coefficients
               if (!(isTRUE(all.equal(cc, c(-.5, 1,   0.95, NA))) ||
                     isTRUE(all.equal(cc, c(-.5, 1,   0.85, NA))) ||
                     isTRUE(all.equal(cc, c(-.5, 1,   NA,  0.5))) ||
                     isTRUE(all.equal(cc, c(-.5, 1.5, 0.95, NA))) ||
                     isTRUE(all.equal(cc, c(-.5, 1.5, 0.85, NA))) ||
                     isTRUE(all.equal(cc, c(-.5, 1.5, NA, 0.5))))) {
		   attr(cc, 'constants') <-
			.psi.ggw.findc(ms=cc[[1]], b=cc[[2]], eff=cc[[3]], bp=cc[[4]])
               }
           },
           "lqq" = { ## only calculate for non-standard coefficients
               if (!(isTRUE(all.equal(cc, c(-.5, 1.5, 0.95, NA))) ||
                     isTRUE(all.equal(cc, c(-.5, 1.5, NA, 0.5))))) {
		   attr(cc, 'constants') <-   ##  b.c :== b/c
		       .psi.lqq.findc(ms=cc[[1]], b.c=cc[[2]], eff=cc[[3]], bp=cc[[4]])
               }
           },
           stop("method for psi function ", psi, " not implemented"))
    cc
}

Mpsi <- function(x, cc, psi, deriv=0) {
    x[] <- .Call(R_psifun, x, .psi.conv.cc(psi, cc), .psi2ipsi(psi), deriv)
    x
}
.Mpsi <- function(x, ccc, ipsi, deriv=0) .Call(R_psifun, x, ccc, ipsi, deriv)

Mchi <- function(x, cc, psi, deriv=0) {
    x[] <- .Call(R_chifun, x, .psi.conv.cc(psi, cc), .psi2ipsi(psi), deriv)
    x
}
.Mchi <- function(x, ccc, ipsi, deriv=0) .Call(R_chifun, x, ccc, ipsi, deriv)

Mwgt <- function(x, cc, psi) {
    x[] <- .Call(R_wgtfun, x, .psi.conv.cc(psi, cc), .psi2ipsi(psi))
    x
}
.Mwgt <- function(x, ccc, ipsi) .Call(R_wgtfun, x, ccc, ipsi)

## only for nlrob() -- and to use instead of MASS:::psi.huber etc:
## returns a *function* a la  psi.huber() :
.Mwgt.psi1 <- function(psi, cc = .Mpsi.tuning.default(psi)) {
    ipsi <- .psi2ipsi(psi)
    ccc <- .psi.conv.cc(psi, cc)
    ## return function *closure* :
    function(x, deriv = 0)
    if(deriv) .Mpsi(x, ccc, ipsi, deriv=deriv) else .Mwgt(x, ccc, ipsi)
}

# The normalizing constant for  rho(.) <--> rho~(.)
MrhoInf <- function(cc, psi) {
    cc <- .psi.conv.cc(psi, cc)
    .Call(R_rho_inf, cc, .psi2ipsi(psi))
}

.MrhoInf <- function(ccc, ipsi) .Call(R_rho_inf, ccc, ipsi)


lmrob.rweights <- function(resid, scale, cc, psi, eps = 16 * .Machine$double.eps) {
    if (scale == 0) { ## exact fit
	m <- max(ar <- abs(resid))
	if(m == 0) numeric(seq_len(ar)) else
	as.numeric(ar < eps * m)# 1 iff res ~= 0
    } else
	Mwgt(resid / scale, cc, psi)
}

lmrob.E <- function(expr, control, dfun = dnorm, use.integrate = FALSE, obj, ...)
{
    expr <- substitute(expr)
    if (missing(control) && !missing(obj))
        control <- obj$control

    lenvir <-
      if (!missing(control)) {
        psi <- control$psi
        if (is.null(psi)) stop('parameter psi is not defined')
	c.psi <- control[[if (control$method %in% c('S', 'SD'))
			  "tuning.chi" else "tuning.psi"]]
        if (!is.numeric(c.psi)) stop('tuning parameter (chi/psi) is not numeric')

        list(psi = function(r, deriv = 0) Mpsi(r, c.psi, psi, deriv),
             chi = function(r, deriv = 0) Mchi(r, c.psi, psi, deriv), ## change?
             wgt = function(r) Mwgt(r, c.psi, psi)) ## change?

      } else list()

    pf <- parent.frame()
    FF <- function(r)
	eval(expr, envir = c(list(r = r), lenvir), enclos = pf) * dfun(r)
    if (isTRUE(use.integrate)) {
	integrate(FF, -Inf,Inf, ...)$value
    ## This would be a bit more accurate .. *AND* faster notably for larger 'numpoints':
    ## } else if(use.integrate == "GQr") {
    ##     require("Gqr")# from R-forge [part of lme4 project]
    ##     ## initialize Gauss-Hermite Integration
    ##     GH <- GaussQuad(if(is.null(control$numpoints)) 13 else control$numpoints,
    ##                     "Hermite")
    ##     ## integrate
    ##     F. <- function(r) eval(expr, envir = c(list(r = r), lenvir), enclos = pf)
    ##     sum(GH$weights * F.(GH$knots))
    } else {
	## initialize Gauss-Hermite Integration
	gh <- ghq(if(is.null(control$numpoints)) 13 else control$numpoints)
	## integrate
	sum(gh$weights * FF(gh$nodes))
    }
}

ghq <- function(n = 1, modify = TRUE) {
    ## Adapted from gauss.quad in statmod package
    ## which itself has been adapted from Netlib routine gaussq.f
    ## Gordon Smyth, Walter and Eliza Hall Institute

    n <- as.integer(n)
    if(n<0) stop("need non-negative number of nodes")
    if(n==0) return(list(nodes=numeric(0), weights=numeric(0)))
    ## i <- seq_len(n) # 1 .. n
    i1 <- seq_len(n-1L)

    muzero <- sqrt(pi)
    ## a <- numeric(n)
    b <- sqrt(i1/2)

    A <- numeric(n*n)
    ## A[(n+1)*(i-1)+1] <- a # already 0
    A[(n+1)*(i1-1)+2] <- b
    A[(n+1)*i1] <- b
    dim(A) <- c(n,n)
    vd <- eigen(A,symmetric=TRUE)
    n..1 <- n:1L
    w <- vd$vectors[1, n..1]
    w <- muzero * w^2
    x <- vd$values[n..1] # = rev(..)
    list(nodes=x, weights= if (modify) w*exp(x^2) else w)
}

.convSs <- function(ss)
    switch(ss,
           "simple"= 0L,
           "nonsingular"= 1L,
           stop("unknown setting for parameter ss"))

