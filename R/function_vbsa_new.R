#' Variance Based Sensitivity Analysis
#' 
#' Variance based Sensitivity analysis for the input factors of a 
#' hydrological model following Saltelli et al., 2010.
#' 
#' @param fn character with the name of a valid R function to be analised.
#' Its first argument MUST contain a K-dimensional variable representing the
#' input factors to be analised. Its output must always be a single value of
#' type numeric!
#' 
#' @param ... additional arguments for 'fn'
#' 
#' @param lower K-dimensional (optionally named) numeric vector with the minimum values for
#' each input factor.
#' 
#' @param upper K-dimensional (optionally named) numeric vector with the maximum values for
#' each input factor.
#' 
#' @param N number of input factor samples used for the low-discrepancy (Sobol')
#' or LHS design of the matrices A and B.
#' 
#' @param rtype character, indicating the type of (quasi-)random sampling design used
#' for matrices A and B. Valid values are: 'quasirandom' (quasi-random low-discrepancy
#' Sobol' sequences) and 'lhs' (latin hypercube sampling)
#' 
#' @param scrambling Shall scrambling be applied to the quasi-random sequences? 
#'  Default value: 'none'. See \code{\link[randtoolbox]{sobol}}.
#'  
#' @param subsamp If given, sensitivity indices will be calculated for
#' (and optionally bootstraping be applied to) sub-samples of N for monitoring of
#' convergence. Either a vector of numeric values <= N or a single numeric value
#' resulting in a sequence of subsamples \code{seq(subsamp, N, by = subsamp)}.
#' NOTE: Values of one make no sense and will be disregarded with a warning.
#'  
#' @param Nb Number of bootstrapping samples to assess the robustness (i.e. a confidence
#' interval) of calculated sensitivity indices. Default: 1 (i.e. bootstrapping is
#' effectively not conducted). Note: This requires no additional evaluations of \code{fn}.
#' 
#' @param Nb.sig Significance level defining upper and lower bounds of the confidence
#' intervals of sensitivity indices. Default: 0.05. Only used if \code{Nb > 1}.
#' 
#' @param ncores \code{integer} specifying the number of CPU cores to be used.
#' If > 1, packages \code{\link[doMC]{doMC}} (Linux only!) and \code{\link[parallel]{parallel}}
#' are needed. Values > 1 only useful if \code{fun} is very complex and computationally
#' demanding, otherwise multiple thread handling will cause function slowdown! Default: 1.
#' 
#' @param verbose logical, should progress messages be printed out to the screen?
#' 
#' @param full.output logical, should the sampling matrices 'A' and 'B' along with all
#' the input factor sets and its corresponding outputs of 'fn' be included in the output
#' of this function? Default: \code{FALSE}.
#' 
#' @param na.handle character. How to handle \code{NA} and non-finite values returned by \code{fn}?
#' One of: \code{stop} - function will stop with an error (default); \code{remove} - problematic
#' values are removed with a warning (affects N, fn.counts, A, B, Ab, fn.A, fn.B, fn.Ab)
#' and the calculation of indices proceeded.
#' 
#' @param debug logical. If \code{TRUE}, the function's complete internal environment
#' will be saved for closer inspection. The following files will be written to \code{getwd()},
#' unless an error occurred beforehand: \code{vbsa_backup1.RData} after initialisation
#' and generation of sampling matrices; \code{vbsa_backup2.RData} after evaluation of
#' \code{fn}; \code{vbsa_backup3.RData} after checking (and possibliy alterating if
#' \code{na.handle = TRUE}) the output of \code{fn}; \code{vbsa_backup4.RData} after
#' calculation of sensitivity indices, before compiling the output. Default: \code{FALSE}.
#' 
#' @return A list with the following elements:
#' 
#'  *) N                      : number of input factor samples employed in the SA
#'  
#'  *) fn.counts              : total number of model evaluations
#'  
#'  *) Si                     : First order sensitivity indices for each of the K input factors (if \code{Nb > 1} this is the mean of the bootstrapping distribution)
#'  
#'  *) St                     : Total effects sensitivity indices for each of the K input factors (if \code{Nb > 1} this is the mean of the bootstrapping distribution)
#'  
#'  *) ranking                : Importance ranking each of the K input factors according to Si
#'  
#'  If \code{Nb > 1} in addition:
#'  
#'  *) Si.sd                  : Standard deviation of the Si bootrstapping distribution
#'  
#'  *) Si.lo                  : Lower bound of the confidence interval of Si (defined by significance level \code{Nb.sig})
#'
#'  *) Si.up                  : Upper bound of the confidence interval of Si (defined by significance level \code{Nb.sig})
#'  
#'  *) St.sd                  : Standard deviation of the St bootrstapping distribution
#'  
#'  *) St.lo                  : Lower bound of the confidence interval of St (defined by significance level \code{Nb.sig})
#'
#'  *) St.up                  : Upper bound of the confidence interval of St (defined by significance level \code{Nb.sig})
#'  
#'  *) ranking.boot           : As ranking above, but for each bootstrapping sample
#'    
#'  OPTIONAL (if \code{full.output = TRUE}):
#' 
#'  *) Matrix.A               : Sampling matrix A of dimension (N, K)
#'  
#'  *) Matrix.B               : Sampling matrix B of dimension (N, K)
#'  
#'  *) Matrix.Ab              : Re-sampled matrix Ab (radial sampling) of dimension (N*K, K)
#'  
#'  *) fn.A                   : N-dimensional vector of fn evaluations for Matrix A
#'  
#'  *) fn.B                   : N-dimensional vector of fn evaluations for Matrix B
#'  
#'  *) fn.Ab                  : Matrix of fn evaluations for Matrix Ab of dimension (N, K)
#'  
#'  If \code{Nb > 1} in addition:
#'  
#'  *) Si.boot                : Matrix of Si bootstrapping distribution of dimension (Nb, K)
#'  
#'  *) St.boot                : Matrix of St bootstrapping distribution of dimension (Nb, K)
#'  
#'  *) B                      : Bootstrap matrix (indices for resampling of fn.A, fn.B, and fn.Ab) of dimension (N, Nb)
#'  
#'  NOTE: if argument \code{subsamp} was given, each sensitivity index vector or matrix
#'  as well as matrix B will be extended by a dimension giving the results for each the sub-sample.
#'  
#' @references Theory:
#'  
#'  Andrea Saltelli, Paola Annoni, Ivano Azzini, Francesca Campolongo,
#'  Marco Ratto, Stefano Tarantola, Variance based sensitivity
#'  analysis of model output. Design and estimator for the total
#'  sensitivity index, Computer Physics Communications, Volume 181,
#'  Issue 2, February 2010, Pages 259-270, DOI: 10.1016/j.cpc.2009.09.018.
#'  \url{http://www.sciencedirect.com/science/article/pii/S0010465509003087}
#'  
#'  Code based on R-script by Mauricio Zambrano-Bigiarini obtained from
#'  \url{https://ec.europa.eu/jrc/en/samo/simlab} (version from 20-May-2013),
#'  based on Matlab code by Stefano Tarantola
#'  
#' @author Tobias Pilz \email{tpilz@@uni-potsdam.de}
#'  
#' @export
#'  
#' @examples 
#'  if ( !require(randtoolbox) ) install.packages(randtoolbox)
#'  if ( !require(lhs) )         install.packages(lhs)
#'
#'  nparam <- 3
#'
#'  # Example: Ishigami test function 
#'  ishigami <- function(x, a=7, b=0.1) {
#'    sin(x[1]) + a*(sin(x[2]))^2 + b*(x[3]^4)*sin(x[1]) 
#'  }
#'  
#'  set.seed(123) # for reproducible results
#'  
#'  res <- vbsa(
#'    fn="ishigami", 
#'    rtype = "quasirandom",
#'    lower=rep(-pi, nparam), 
#'    upper=rep(pi, nparam),
#'    N=2000, scrambling="none"
#'  )
#'  
vbsa <- function(
  fn,  
  ...,
  lower=-Inf,
  upper=Inf,
  N=1000,
  rtype=c("quasirandom", "lhs"),   
  scrambling=c("Owen","Faure-Tezuka","Owen+Faure-Tezuka","none"),          # See ?randtoolbox::sobol
  subsamp = NULL,
  Nb = 1,
  Nb.sig = 0.05,
  ncores=1,
  verbose= TRUE,                  # logical, indicating if progress messages have to be printed          
  full.output=FALSE,
  na.handle = "stop",
  debug = FALSE
) {

  #### INPUT CHECKING AND INITIALISATIONS ####
  if(verbose) {
    message("                                                              ")
    message("[ START SENSITIVITY ANALYSIS! ]")
  }
  
  # Assigning the objective function
  fn <- match.fun(fn)
  
  # checking length of 'lower' and 'upper'
  if (length(lower) != length(upper) )
    stop( paste( "Invalid argument: 'length(lower) != length(upper) (", length(lower), "!=", length(upper), ")'", sep="" ) )        
  
  rtype      <- match.arg(rtype, c("quasirandom", "lhs"))  
  scrambling <- match.arg(scrambling, c("Owen","Faure-Tezuka","Owen+Faure-Tezuka","none"))  
  scrambling <- pmatch(scrambling, c("Owen","Faure-Tezuka","Owen+Faure-Tezuka","none"))  
  if (scrambling > 3) scrambling <- 0
  
  # parallel version
  if(ncores > 1)
    registerDoMC(ncores)
  
  # Checking that 'N' is integer
  if ( trunc(N) != N ) stop( "Invalid argument: 'N' must be integer !" )
  
  # check na.handle
  if(!any(na.handle == c("stop", "remove"))) stop("Invalid argument: 'na.handle' must be one of {'stop', 'remove'}!")
  
  # checking 'X.Boundaries' 
  if ( (lower[1L] == -Inf) || (upper[1L] == Inf) ) {
    #if (any(upper==Inf | lower==-Inf))
    stop( "Invalid argument: 'lower' and 'upper' boundaries must be finite !!'" )
  } else X.Boundaries <- cbind(lower, upper)              
  
  # Computing 'K', the Dimension of the Solution Space
  K <- nrow(X.Boundaries)
  
  # Meaningful name of each one of the parameters
  if (is.null(rownames(X.Boundaries))) {
    param.IDs <- paste("Param", 1:K, sep="")
  } else param.IDs <- rownames(X.Boundaries)   
  
  if (rtype=="quasirandom") {
    # Backing up the original boundaries
    lower.ini <- lower
    upper.ini <- upper
    X.Boundaries.ini <- X.Boundaries
    LOWER.ini <- matrix( rep(lower.ini, N), nrow=N, byrow=TRUE)
    UPPER.ini <- matrix( rep(upper.ini, N), nrow=N, byrow=TRUE)
    
    # normalising
    lower <- rep(0, K)
    upper <- rep(1, K)
    X.Boundaries <- cbind(lower, upper)
    rownames(X.Boundaries) <- param.IDs
  } # IF end 
  
  # check for sub-sampling
  if(!is.null(subsamp)) {
    if(!is.numeric(subsamp)) stop("Argument 'subsamp' must be of type numeric or NULL!")
    if(any(subsamp == 1)) {
      warning("Argument 'subsamp' should only contain values > 1. Values of 1 will be disregarded!")
      subsamp <- subsamp[-which(subsamp == 1)]
      if(length(subsamp) == 0) subsamp <- N
    }
    if(length(subsamp) == 1) {
      subsamp <- seq(subsamp, N, by = subsamp)
    }
    if(max(subsamp) > N) stop("Maximum value of 'subsamp' must be <= N!")
    if(max(subsamp) < N) subsamp <- c(subsamp, N)
  } else {
    subsamp <- N
  }
  
  # Total Number of parameter sets to be used in the GSA
  nparamsets  <- N*(K+2)
  
  if(verbose) message("[ Number of Parameter Sets to be run ( N*(K+2) ): ", nparamsets, " ]")
  

  #### MAIN BODY ####
  
  #### Matrix sampling ####
  if(verbose) message("[ Generate sampling matrices A, B, and Ab ]")
  
  # Matrices A and B
  if (rtype=="quasirandom") {
    rs <- randtoolbox::sobol(n=N, dim=2*K, scrambling=scrambling)  
    A <- rs[, 1:K] * (UPPER.ini - LOWER.ini) + LOWER.ini
    B <- rs[, (K+1):(2*K)] * (UPPER.ini - LOWER.ini) + LOWER.ini
  } else if (rtype=="lhs") {
    # tp corrected (rLHS which does not exist replaced by randomLHS)
    A <- randomLHS(n=N, k = K)
    A <- t(apply(A, 1, function(x) qunif(x, min = X.Boundaries[,1], max = X.Boundaries[,2])))
    B <- randomLHS(n=N, k = K)
    B <- t(apply(B, 1, function(x) qunif(x, min = X.Boundaries[,1], max = X.Boundaries[,2])))
    # tp end
  } # ELSE end
  
  # Matrix Ab dim(N*K,K) by radial sampling
  Ab <- A[rep(1:nrow(A), each = K),]
  for(k in 1:K) Ab[seq(k, N*K, K), k] <- B[,k]
  
  colnames(A) <- param.IDs
  colnames(B) <- param.IDs
  colnames(Ab) <- param.IDs
  
  if(debug) save(list = ls(all.names = TRUE), file = "vbsa_backup1.RData") 
  
  
  #### Function evaluation at sample points ####
  
  if(ncores > 1) {
    if(verbose) message(paste("[ Evaluation of matrix A in parallel mode, involves", N, "calls of fn ]"))
    yA <- foreach(j=1:nrow(A), .combine = "c", .errorhandling = "pass") %dopar% fn(A[j,], ...)
    if(verbose) message(paste("[ Evaluation of matrix B in parallel mode, involves", N, "calls of fn ]"))
    yB <- foreach(j=1:nrow(B), .combine = "c", .errorhandling = "pass") %dopar% fn(B[j,], ...)
    if(verbose) message(paste("[ Evaluation of matrix Ab in parallel mode, involves", N*K, "calls of fn ]"))
    yAb <- foreach(j=1:nrow(Ab), .combine = "c", .errorhandling = "pass") %dopar% fn(Ab[j,], ...)
  } else {
    if(verbose) message(paste("[ Evaluation of matrix A, involves", N, "calls of fn ]"))
    yA <- apply(A, 1, fn, ...)
    if(verbose) message(paste("[ Evaluation of matrix B, involves", N, "calls of fn ]"))
    yB <- apply(B, 1, fn, ...)
    if(verbose) message(paste("[ Evaluation of matrix Ab, involves", N*K, "calls of fn ]"))
    yAb <- apply(Ab, 1, fn, ...)
  }
  
  if(debug) save(list = ls(all.names = TRUE), file = "vbsa_backup2.RData") 
  
  # check results
  if(any(!is.numeric(c(yA, yB, yAb)))) {
    stop("Evaluation of 'fn' produced unexpected non-numeric results! Check your 'fn' to always return a numeric value!")
  }
  if(length(yA) != N) stop("Evaluation of 'fn' returned more results than expected for matrix A!")
  if(length(yB) != N) stop("Evaluation of 'fn' returned more results than expected for matrix B!")
  if(length(yAb) != N*K) stop("Evaluation of 'fn' returned more results than expected for matrix Ab!")
  if(any(!is.finite(c(yA, yB, yAb)))) {
    if(na.handle == "stop") stop("Evaluation of 'fn' produced non-finite results! Consider argument 'na.handle'.")
    
    if(na.handle == "remove") {
      n_na <- lapply(list(yA, yB, yAb), function(x) which(!is.finite(x)))
      n_na[[3]] <- ceiling(n_na[[3]] / K)
      n_na <- unique(unlist(n_na))
      if(length(n_na) > 0.5*N) {
        stop(paste("Evaluation of 'fn' produced many non-finite results! Stopping function as", length(n_na)/N*100, "% of N would have to be removed. Check your function and/or parameter ranges!"))
      } else {
        warning("Non-finite output of 'fn' detected! Removing problematic values and adapting parameters (check element N of output).")
      }
      N <- N - length(n_na)
      if(!is.null(subsamp)) subsamp[length(subsamp)] <- N
      nparamsets  <- N*(K+2)
      A <- A[-n_na,]
      yA <- yA[-n_na]
      B <- B[-n_na,]
      yB <- yB[-n_na]
      # can be helpful here...
      seq_vectorised <- Vectorize(seq.default, vectorize.args = c("from", "to"))
      Ab <- Ab[-c(seq_vectorised((n_na-1)*K+1, n_na*K)),]
      yAb <- yAb[-c(seq_vectorised((n_na-1)*K+1, n_na*K))]
    }
  }
  
  if(debug) save(list = ls(all.names = TRUE), file = "vbsa_backup3.RData") 
  
  
  #### Computation of sensitivity indices ####
  if(verbose) message("[ Calculate sensitivity indices ]")
  
  # yAb as a matrix from here on
  yAb <- matrix(yAb, ncol = K, byrow = T)
  
  # calculate indices for each sub-sampling and bootstrapping realisation
  ind_out <- lapply(subsamp, function(s) indices_vb_boot(yA[1:s], yB[1:s], yAb[1:s,], s, Nb, K))
  
  # convert output
  Si <- sapply(ind_out, function(x) t(sapply(x$ind_out, function(y) y$Si)), simplify = "array")
  if(dim(Si)[3] > 1) {
    Si <- aperm(Si, perm = c(3,1,2))
    dimnames(Si) <- list(subsamp, NULL, param.IDs)
  } else {
    Si <- matrix(Si, ncol=K)
    colnames(Si) <- param.IDs
  }
  
  St <- sapply(ind_out, function(x) t(sapply(x$ind_out, function(y) y$St)), simplify = "array")
  if(dim(St)[3] > 1) {
    St <- aperm(St, perm = c(3,1,2))
    dimnames(St) <- list(subsamp, NULL, param.IDs)
  } else {
    St <- matrix(St, ncol=K)
    colnames(St) <- param.IDs
  }
  
  Boot <- lapply(ind_out, function(x) x$B)
  names(Boot) <- subsamp
  
  if(debug) save(list = ls(all.names = TRUE), file = "vbsa_backup4.RData") 
  
  
  #### OUTPUT ####
  if(verbose) message("[ Compile output ]")
  
  # generic output
  if(length(dim(Si)) == 3) apply_dims <- c(1,3) else apply_dims <- 2
  out <- list(
    N = N,
    fn.counts = nparamsets,
    Si = apply(Si, apply_dims, mean),
    St = apply(St, apply_dims, mean)
  )
  
  # parameter ranking based on Si
  if(length(dim(Si)) > 2) {
    ranks <- t(apply(out$Si, 1, order, decreasing = T))
    colnames(ranks) <- colnames(out$Si)
  } else {
    ranks <- order(out$Si, decreasing = T)
    names(ranks) <- names(out$Si)
  }
  out <- c(out, list(ranking = ranks))
  
  # Bootstrapping output
  if(Nb > 1) {
    n_min <- max(1, round(Nb * Nb.sig/2))
    n_max <- round(Nb * (1 - Nb.sig/2) )
    if(length(dim(Si)) == 3) {
      Si_sorted <- aperm(apply(Si, c(1,3), sort), c(2,1,3))
      Si.lo <- Si_sorted[,n_min,]
      Si.up <- Si_sorted[,n_max,]
      St_sorted <- aperm(apply(St, c(1,3), sort), c(2,1,3))
      St.lo <- St_sorted[,n_min,]
      St.up <- St_sorted[,n_max,]
    } else {
      Si_sorted <- apply(Si, 2, sort)
      Si.lo <- Si_sorted[n_min,]
      Si.up <- Si_sorted[n_max,]
      St_sorted <- apply(St, 2, sort)
      St.lo <- St_sorted[n_min,]
      St.up <- St_sorted[n_max,]
    }
    out.distr <- list(
      Si.sd = apply(Si, apply_dims, sd),
      Si.lo = Si.lo,
      Si.up = Si.up,
      St.sd = apply(St, apply_dims, sd),
      St.lo = St.lo,
      St.up = St.up
    )
    out <- c(out, out.distr)
    
    # ranking
    if(length(apply_dims) > 1) {
      ranks.boot <- aperm(apply(Si, c(1,2), order, decreasing = T), c(3,1,2))
      dimnames(ranks.boot)[[2]] <- param.IDs
    } else {
      ranks.boot <- t(apply(Si, 1, order, decreasing = T))
      colnames(ranks.boot) <- param.IDs
    }
    out <- c(out, list(ranking.boot = ranks.boot))
  }
  
  # optional output
  if(full.output) {
    out.opt <- list(
      Matrix.A = A,
      Matrix.B = B,
      Matrix.Ab = Ab,
      fn.A = yA,
      fn.B = yB,
      fn.Ab = yAb
    )
    out <- c(out, out.opt)
    
    if(Nb > 1) {
      out.distr <- list(
        Si.boot = Si,
        St.boot = St,
        B = Boot
      )
      out <- c(out, out.distr)
    }
  }
  
  if(verbose) message("[ DONE! ]")
  return(out)
}
