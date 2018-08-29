#' Variance Based Sensitivity Analysis
#' 
#' Variance based Sensitivity analysis for the input factors of a 
#' hydrological model following Saltelli et al., 2010.
#' 
#' @param fn character with the name of a valid R function to be analised.
#' Its first argument MUST contain a K-dimensional variable representing the
#' input factors to be analised.
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
#' @return A list with the following elements:
#' 
#'  *) N                      : number of input factor samples employed in the SA
#'  
#'  *) fn.counts              : total number of model evaluations
#'  
#'  *) Si                     : First order sensitivity indices for each of the K input factors
#'  
#'  *) St                     : Total effects sensitivity indices for each of the K input factors
#'  
#'  OPTIONAL if \code{full.output = TRUE}
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
#'  *) fn.Ab                  : N*K-dimensional vector of fn evaluations for Matrix Ab
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
  ncores=1,
  verbose= TRUE,                  # logical, indicating if progress messages have to be printed          
  full.output=FALSE
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
  
  
  #### Function evaluation at sample points ####
  
  if(ncores > 1) {
    if(verbose) message(paste("[ Evaluation of matrix A in parallel mode, involves", N, "calls of fn ]"))
    yA <- mclapply(1:nrow(A), function(j) fn(A[j,], ...), mc.cores = ncores)
    yA <- unlist(yA)
    if(verbose) message(paste("[ Evaluation of matrix B in parallel mode, involves", N, "calls of fn ]"))
    yB <- mclapply(1:nrow(B), function(j) fn(B[j,], ...), mc.cores = ncores)
    yB <- unlist(yB)
    if(verbose) message(paste("[ Evaluation of matrix Ab in parallel mode, involves", N*K, "calls of fn ]"))
    yAb <- mclapply(1:nrow(Ab), function(j) fn(Ab[j,], ...), mc.cores = ncores)
    yAb <- unlist(yAb)
  } else {
    if(verbose) message(paste("[ Evaluation of matrix A, involves", N, "calls of fn ]"))
    yA <- apply(A, 1, fn, ...)
    if(verbose) message(paste("[ Evaluation of matrix B, involves", N, "calls of fn ]"))
    yB <- apply(B, 1, fn, ...)
    if(verbose) message(paste("[ Evaluation of matrix Ab, involves", N*K, "calls of fn ]"))
    yAb <- apply(Ab, 1, fn, ...)
  }
  
  
  #### Computation of sensitivity indices ####
  if(verbose) message("[ Calculate sensitivity indices ]")
  
  # total variance, TODO: not really sure about this, in the SAFER implemantation it's var(yA)
  Vtot <- var(c(yA, yB))
  
  # first order index, Table 2 (b)
  Si <- sapply(1:N, function(n) yB[n] * (yAb[((n-1)*K+1):(n*K)] - yA[n]) )
  Si <- apply(Si, 1, mean) / Vtot
  names(Si) <- param.IDs
  
  # total effects index, Table 2 (f)
  St <- sapply(1:N, function(n) (yA[n] - yAb[((n-1)*K+1):(n*K)])^2 )
  St <- apply(St, 1, mean) / (2*Vtot)
  names(St) <- param.IDs
  
  
  #### OUTPUT ####
  if(verbose) message("[ Compile output ]")
  
  # generic output
  out <- list(
    N = N,
    fn.counts = nparamsets,
    Si = Si,
    St = St
  )
  
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
  }
  
  if(verbose) message("[ DONE! ]")
  return(out)
}
