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
#' @param reltol numeric, relative tolerance. If a value different from zero,
#' the SA will finish when the maximum difference between the first
#' order sensitivity indices of all the input factors ('DSi') satisfy
#' the following relation:  Di = Max {|S[i] - S[i-1]} <= reltol
#' 
#' @param ncores \code{integer} specifying the number of CPU cores to be used.
#' If > 1, packages \code{\link[doMC]{doMC}} (Linux only!) and \code{\link[parallel]{parallel}}
#' are needed. Values > 1 only useful if \code{fun} is very complex and computationally
#' demanding, otherwise multiple thread handling will cause function slowdown! Default: 1.
#' 
#' @param verbose logical, should progress messages be printed out to the screen?
#' 
#' @param REPORT numeric, only used when 'verbose=TRUE'. A progress message is
#' printed every 'REPORT' number of model evaluations.
#' 
#' @param full.output logical, should the sampling matrices 'A' and 'B' along with all
#' the input factor sets and its corresponding outputs of 'fn' be included in the output
#' of this function? Default: \code{FALSE}.
#' 
#' @param monitor.ind logical, shall sensitivity indices Si and St be monitored?
#' If \code{FALSE} (default), output elements \code{First.Order.Indices.Si} and
#' \code{Total.Order.Indices.Sti} are merely K-dimensional numerical vectors considering
#' the full sample size. If \code{TRUE}, they are N-by-K dimensional matrices.
#' 
#' @return A list with the following elements:
#' 
#'  *) N                      : number of input factor samples employed in the SA
#'  
#'  *) counts                 : total number of model evaluations
#'  
#'  *) First.Order.Indices.Si : First order sensitivity indices for each of the K input factors
#'  
#'  *) Sum.Si                 : sum of all the K first order sensitivity indices Si
#'  
#'  *) Total.Order.Indices.Sti: Total effects sensitivity indices for each of the K input factors
#'  
#'  *) Total.Indices.of.Pairs.of.Factors: matrix with the total indices of pairs of factors
#'  for each possible combination of input factors
#'  
#'  *) Ranking                : matrix with a ranking of the most important input factors,
#'  sorted in decreasing order according to the first order sensitivity index Si
#' 
#'  *) Matrix.A               : (optional). A matrix used in the radial sampling
#'  
#'  *) Matrix.B               : (optional). B matrix used in the radial sampling
#'  
#'  *) ParameterSets          : (optional) matrix with all the parameter sets used in the analysis
#'  
#'  *) GoFs                   : (optional) numeric vector with the outputs of 'fn' for each
#'  input factor set in 'ParameterSets'
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
#'    N=2000, scrambling="none", REPORT=100
#'  )
#'  
vbsa_old <- function(
  fn,  
  ...,
  lower=-Inf,
  upper=Inf,
  N=1000,                         
  rtype=c("quasirandom", "lhs"),   
  scrambling=c("Owen","Faure-Tezuka","Owen+Faure-Tezuka","none"),          # See ?randtoolbox::sobol
  reltol=0,                       # 
  ncores=1,
  verbose= TRUE,                  # logical, indicating if progress messages have to be printed          
  REPORT=10, 
  full.output=FALSE,
  monitor.ind=FALSE
) {
  ##############################################################################
  #                            INPUT CHECKING                                  #
  ##############################################################################
  
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
  
  ########################################################################
  ##################### Dummy checkings ##################################
  
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
  
  # Checking report
  if (N < REPORT) {
    REPORT <- N
    warning("[ 'REPORT' is greater than 'N' => 'REPORT = N' ]")
  } # IF end
  
  message("                                                              ")
  message("[ Number of Parameter Sets to be run ( N*(K+2) ): ", nparamsets, " ]")     
  
  ##############################################################################
  #                            MAIN BODY
  ##############################################################################
  
  # sampling design
  if (rtype=="quasirandom") {
    rs <- randtoolbox::sobol(n=N, dim=2*K, scrambling=scrambling)  
    A <- rs[, 1:K]
    B <- rs[, (K+1):(2*K)]
  } else if (rtype=="lhs") {
    # tp corrected (rLHS which does not exist replaced by randomLHS)
    A <- randomLHS(n=N, k = K)
    A <- t(apply(A, 1, function(x) qunif(x, min = X.Boundaries[,1], max = X.Boundaries[,2])))
    B <- randomLHS(n=N, k = K)
    B <- t(apply(B, 1, function(x) qunif(x, min = X.Boundaries[,1], max = X.Boundaries[,2])))
    # tp end
  } # ELSE end
  
  # Output Matrix with First Order Indices
  S           <- matrix(NA, nrow=N, ncol=K)
  colnames(S) <- param.IDs
  
  # Output Matrix with Total Order Indices
  St          <- S
  
  # Matrixes for first and total order variances
  Vi <- matrix(NA, nrow=N, ncol=K)
  Vt <- Vi
  
  # Vector of model outputs used to compute the output variance
  y4var <- rep(NA, 2*N)
  
  # Preparation of the radial sample matrix X [K+2,K]
  X           <- matrix(NA, nrow=nparamsets, ncol=K)
  colnames(X) <- param.IDs
  
  # Numerators for total indices of pairs of factors
  Vtij <- array( NA, dim=c(N, K, K) )
  
  # Si values in the (i-1)-th point
  Si.previous <- rep(0, K)
  
  # Differences between the first order sensitivity index in the i-th point and 
  # the previous one
  Si.Delta <- rep(NA, K)
  
  # Vector with goodness-of-fit of model outputs
  y <- rep(NA, nparamsets)
  
  # Vector with model outputs
  ModelOut <- vector("list", nparamsets)
  
  # Counter of the points belonging to the initial LHS
  i <- 1    
  
  # relative convergence
  reltol.conv <- FALSE
  
  ##############################################################################
  #                4) Loop for each row of the A and B matrixes                #
  ##############################################################################
  while ( (i <= N) & !reltol.conv ) {   
    
    # Matrix with parameter sets
    ParamSets  <- matrix(NA, nrow=K+2, ncol=K)
    
    # Creating of the Ab matrix for radial sampling
    Ab <- matrix(NA, nrow=K+2, ncol=K)
    
    # Goodness-of-fit of each parameter set used in the GSA
    gof <- rep(NA, K+2)
    
    if (rtype == "quasirandom") {
      An <- A * (UPPER.ini - LOWER.ini) + LOWER.ini
      Bn <- B * (UPPER.ini - LOWER.ini) + LOWER.ini
    } else {
      An <- A
      Bn <- B
    } # ELSE end
    
    # Creating of the Ab matrix for radial sampling
    Ab      <- matrix(NA, nrow=K+2, ncol=K)
    Ab[1, ] <- An[i, ]
    Ab[2, ] <- Bn[i, ]
    for (j in 1:K){
      Ab[j+2, ]  <- An[i, ]
      Ab[j+2, j] <- Bn[i, j]
    } # FOR end    
    
    Xn <- Ab
    colnames(Xn) <- rownames(X.Boundaries)
    
    # Filling the X matrix
    X[((K+2)*(i-1) + 1):((K+2)*i), ]  <- Xn
    
    ############################################################################
    ############################################################################
    # Loop for each parameter of the i point of the initial LHS
    # 3.a) Evaluate the particles fitness
    
    if(ncores > 1) {
      GoF <- mclapply(1:nrow(Xn), function(j) fn(Xn[j,], ...), mc.cores = ncores)
      GoF <- unlist(GoF)
    } else {
      GoF <- apply(Xn, fn, MARGIN=1, ...)
    }
    
    gof[1:(K+2)]      <- GoF
    ModelOut[1:(K+2)] <- GoF  ###
    
    ##########################################################################
    # 9)                  Updating the sensitivity matrix                    #                                 
    ##########################################################################
    
    # Storing the GoF corresponding to the 'i' parameter set and the related ones
    y[((K+2)*(i-1) + 1):((K+2)*i)]  <- gof
    
    yA  <- gof[1]
    yB  <- gof[2]
    yAb <- gof[3:(K+2)]
    
    # Vector for computation of total variance
    y4var[(2*i-1):(2*i)] <- c(yA, yB)
    
    # Numerator for first and total Indices
    Vi[i, 1:K] <- yB*(yAb[1:K]-yA)
    Vt[i, 1:K] <- (yA - yAb[1:K])^2
    
    # Computation of total variance
    Vtot <- var(y4var, na.rm=TRUE)
    
    # Computation of first and total order sensitivity indices (Si and St, respectively)
    S[i, ]  <- apply(Vi, MARGIN=2, FUN=mean, na.rm=TRUE) / Vtot
    St[i, ] <- apply(Vt, MARGIN=2, FUN=mean, na.rm=TRUE) / (2*Vtot)   
    
    # Numerators for total indices of pairs of factors (no extra computational cost)
    # Vtij[N,k,k] is a set of triangular matrixes
    for (j in 1:K) {
      yAbi <- yAb[j]
      j2 <- j +1
      while(j2 <= K) {
        yAbj <- yAb[j2]
        Vtij[i, j, j2] <-  (yAbi - yAbj) * (yAbi-yAbj)
        j2 <- j2 + 1
      } # FOR 'j2' END
    } # FOR 'j' END
    
    # Updating values of the first order sensitivity index for the i-th point and
    # its difference with the previous point    
    Si.Delta    <- abs( S[i, ] - Si.previous )
    Si.Delta.mx <- max(Si.Delta)
    Si.previous <- S[i, ]
    
    if (reltol==0) {
      reltol.conv <- FALSE
    } else {
      if (Si.Delta.mx <= reltol) {
        reltol.conv <- TRUE
      } else reltol.conv <- FALSE
    } # ELSE end
    
    
    if ( (i/REPORT == floor(i/REPORT)) & (verbose) ) {
      if (K <= 7) {
        message( "[ Model Runs : ", 
                 format( i*(K+2), width=7, justify="left" ), "/", nparamsets, ". ",   
                 paste("dS", 1:K, ": ", format(Si.Delta, scientific=TRUE, digits=3, width=9), ".  ", collapse=" ", sep=""),  
                 "]" )
      } else {
        message( "[ Model Runs : ", 
                 format( i*(K+2), width=7, justify="left" ), "/", nparamsets, ". ",   
                 paste("Max(|dS|): ", 
                       formatC(max(Si.Delta), format="E", digits=3, flag=" "), 
                       " (i=", which.max(Si.Delta), 
                       "). Min(|dS|): ", 
                       formatC(min(Si.Delta), format="E", digits=3, flag=" "),
                       " (i=", which.min(Si.Delta), 
                       ")", sep=""),  "]" )
      } # ELSE end
    } # IF end
    
    i <- i + 1    
    
  } # WHILE i end
  
  ##############################################################################
  # 7)                    Sensitivity of each Parameter                       #                                 
  ##############################################################################
  
  # Computation of totals indices of pairs of factors (no extra computational cost)
  Stij <- matrix(NA, nrow=K, ncol=K)
  for (j in 1:K) {
    j2 <- 1
    while (j2 <= K) {
      Stij[j2, j] <-  mean(Vtij[, j, j2], na.rm=TRUE) / (2*Vtot)
      j2 <- j2 + 1
    } # FOR 'j2' END
  } # FOR 'j' END
  
  # Number of totals indices of pairs of factors
  #nStij <- 1 + sum(1:(K-1)) # including St11
  nStij <- sum(1:(K-1)) # not including St11
  
  # Creating a vector with names for Stij
  Stij.names <- rep(NA, nStij)
  
  # Creating a vector with values for Stij
  Stij.values <- rep(NA, nStij)
  
  #p <- 2 # including St11
  p <- 1 # not including St11
  for (i in 1:K) { 
    j <- i + 1
    while (j <= K) {
      Stij.names[p]  <- paste("St_", param.IDs[i], "_", param.IDs[j], sep="")
      Stij.values[p] <- Stij[j,i]
      j <- j + 1 
      p <- p + 1
    } # WHILE end
  } # FOR end
  colnames(Stij) <- param.IDs 
  rownames(Stij) <- param.IDs 
  
  ##############################################################################
  # 9)          Computing  the  Final  Ranking  of  Sensitivity                #
  #                sorted by First Order Sensitivity Indices                   #                                 
  ##############################################################################
  if(monitor.ind) {
    # First order sensitivity indices: final values
    Si.out        <- S
    
    # Total order sensitivity indices: final values
    St.out        <- St
  } else {
    # First order sensitivity indices: final values
    Si.out        <- S[N, ]
    names(Si.out) <- param.IDs
    
    # Total order sensitivity indices: final values
    St.out        <- St[N, ]
    names(St.out) <- param.IDs
  }
  
  # Sorting the parameters, from the most sensitive to the least one
  Ranking <- sort(S[N,], decreasing=TRUE, na.last=TRUE)
  
  # Parameter Order
  index <- pmatch(names(Ranking), param.IDs)
  
  # Adding a column with a sequential number for the ranking
  Ranking <- data.frame(Ranking.Nmbr=format(as.character(1:K), width=11, justify="left"), 
                        Parameter.Name=format(param.IDs[index], width=13, justify="left"), 
                        First.Order.Index=as.numeric(S[N,index]),
                        Total.Order.Index=as.numeric(St[N,index]) 
  )                        
  Ranking[, "Ranking.Nmbr"]   <- as.character(Ranking[, "Ranking.Nmbr"])
  Ranking[, "Parameter.Name"] <- as.character(Ranking[, "Parameter.Name"]) 
  
  ##############################################################################
  # 10)                    Creating the output                                 #                                 
  ##############################################################################
  if (full.output) {   
    nelements <- 11
    first     <- 5
    colnames(A) <- param.IDs 
    colnames(B) <- param.IDs 
  } else {
    nelements <- 7
    first     <- 1
  } # ELSE end
  
  # Creating the R output
  ## "pre-allocate" an empty list of length 'nelements'
  out <- vector("list", nelements)
  
  out[[first]]   <- N           # number of points used in the random design
  out[[first+1]] <- nparamsets  # total number of parameter sets   
  out[[first+2]] <- Si.out          # First.Order.Indices
  out[[first+3]] <- sum(S[N,])     # Sum of First.Order.Indices
  out[[first+4]] <- St.out          # Total.Order.Indices
  out[[first+5]] <- Stij        # Total Indices of Pairs of Factors
  out[[first+6]] <- Ranking     # Ranking
  names(out)[first:(first+6)] <- c("N", "counts", "First.Order.Indices.Si", 
                                   "Sum.Si", "Total.Order.Indices.Sti", 
                                   "Total.Indices.of.Pairs.of.Factors", 
                                   "Ranking")    
  
  if (full.output) {
    
    out[[1]] <- An # Matrix A
    out[[2]] <- Bn # Matrix B
    out[[3]] <- X  # All Parameter Sets
    out[[4]] <- y  # all goodness-of-fit values
    
    names(out)[1:4] <- c("Matrix.A", "Matrix.B", "ParameterSets", "GoFs")
    
  } # IF end
  
  return(out)
}
