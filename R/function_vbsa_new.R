#' Variance Based Sensitivity Analysis
#'
#' Variance based Sensitivity analysis for the input factors of a
#' hydrological model following Saltelli et al., 2010.
#'
#' @param fn character with the name of a valid R function to be analised.
#' Its first argument MUST contain a K-dimensional variable representing the
#' input factors to be analised. Its output must either be a single value of
#' type numeric or an optionally named numeric vector. If the latter, sensitivity
#' indices will be calculated and output be compiled for each of the fn outputs.
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
#' If \code{fn} returns multivariate output, there might be superfluous output messages!
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
#' @param restart Character list of two elements: 'backup': path to backup files
#' created by a previous run with \code{debug = TRUE}; 'log': path to logfiles
#' containing parameter values and \code{fn} outputs from a previous run (dee details).
#'
#' @return For each output of \code{fn} an optionally named (sub-)list with the following elements:
#'
#'  *) N                      : number of input factor samples employed in the SA
#'
#'  *) fn.counts              : total number of model evaluations
#'
#'  *) Si                     : First order sensitivity indices for each of the K input factors (if \code{Nb > 1} this is the mean of the bootstrapping distribution)
#'
#'  *) St                     : Total effects sensitivity indices for each of the K input factors (if \code{Nb > 1} this is the mean of the bootstrapping distribution)
#'
#'  *) ranking                : Importance ranking of input factors from 1 (highest) to K (lowest importance) according to Si
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
#' @note
#' If \code{debug = TRUE}, backup files \code{vbsa_backup3.RData} and \code{vbsa_backup4.RData}
#' will contain reasonable output only if \code{fn} returns a single value! Otherwise the output
#' in these two files relate to the last output variable in the output list/vector of \code{fn}!
#'
#' If \code{fn} returns (lots of) values equal to zero, no total variance can be computed
#' and \code{NA} will be returned for the respective cases (a warning will be given).
#' Note that NA values might cause trouble with some output metrics! Unless it is a problem
#' with your model, you could try to increase \code{N} and/or reduce subsampling intervalls
#' (i.e. increase N for each subsample) to avoid NAs.
#' 
#' @details 
#' There is an option to re-start previous function executions that terminated with
#' an error. To do so, use argument \code{restart}. However, currently this is still
#' rather unflexible. There must be logfiles in directory \code{restart$log} detectable
#' via regexp '^log[-_0-9a-zA-Z]*.dat$'. Each file has to correspond to a single
#' evaluation of \code{fn}. It must be a table with the elements 'group' (requires
#' entries 'pars' and 'output'), 'variable' (the parameter and output element names),
#' and 'value' (the actual parameter / output value). Internally, the parameter
#' realisations found in the logfiles will be compared to the parameter matrices stored
#' in debug file 'vbsa_backup1.RData'. Missing evaluations of \code{fn} will then be
#' conducted and sensitivity indices be calculated. If \code{restart$log = NULL},
#' only the indices will be calculated based on information stored in 'vbsa_backup2.RData',
#' which will be required. For a re-start, the same parametrisation will be used,
#' except for parameters \code{ncores}, \code{verbose}, and (obviously) \code{restart}.
#' Changes to other parameters will be ignored!
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
  debug = FALSE,
  restart = NULL
) {

  #### INPUT CHECKING AND INITIALISATIONS ####
  if(verbose) {
    message("                                                              ")
    message("[ START SENSITIVITY ANALYSIS! ]")
  }

  # check if this is a restart of a previous run
  if(!is.null(restart)) {
    
    if(verbose) message("[ Re-start of a previous run ]")
    
    # parallel version
    if(ncores > 1)
      registerDoMC(ncores)
    
    # save runtime parameters that might have changed for restart before they get overwritten when loading vbsa_backup1.RData
    ncores_save <- ncores
    verbose_save <- verbose
    restart_save <- restart
    
    # load vbsa_backup1.RData
    if(!file.exists(paste(restart$backup, "vbsa_backup1.RData", sep="/")))
      stop("Could not find backup file 'vbsa_backup1.RData' which is required!")
    
    load(paste(restart$backup, "vbsa_backup1.RData", sep="/"))
    
    # restore saved parameters
    ncores <- ncores_save
    verbose <- verbose_save
    restart <- restart_save
    
    if(verbose) message("[ Loaded vbsa_backup1.RData ]")
    
  } else { # no restart
    
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
  
  } # restart?


  #### Function evaluation at sample points ####

  # if restart, check for vbsa_backup2.RData (in this case no fn evaluations are necessary)
  do_eval <- T
  if(!is.null(restart)) {
    
    # load vbsa_backup2.RData
    if(is.null(restart$log) & file.exists(paste(restart$backup, "vbsa_backup2.RData", sep="/"))) {
      load(paste(restart$backup, "vbsa_backup2.RData", sep="/"))
      do_eval <- F
      if(verbose) message("[ Found and loaded vbsa_backup2.RData ]")
      
    } else { # try to get results from runs so far
      
      # param matrices as tidy df
      dat_pars <- map_dfr(list("A"=A, "B"=B, "Ab"=Ab), function(m) {
        as_tibble(round(m, 4)) %>%
          mutate(n = 1:n())
      }, .id = "matrix")
      
      # get data from logfiles
      logfiles <- dir(restart$log, pattern = "^log[-_0-9a-zA-Z]*.dat$", full.names = T)
      if(length(logfiles) == 0) stop("Could not find logfiles!")
      dat_log <- map_dfr(logfiles, function(f) {
        dat_t <- read_tsv(f, col_types = "ccc")
        if(nrow(dat_t) == 0) return(NULL)
        dat_t %>%
          filter(group %in% c("pars", "output")) %>%
          mutate(value = as.numeric(value)) %>%
          select(-group) %>%
          spread(key = variable, value = value)
      }, .id = "file")
      
      # identify rows (col 'n') of dat_pars from data in dat_log
      dat_key <- left_join(dat_log, dat_pars, by = param.IDs)
      
      if(!all(dat_key$matrix %in% c("A", "B", "Ab")) || !all(dat_key$n %in% dat_pars$n))
        stop("Could not merge log data with parameter matrices!")
      
      # get rows in A, B, Ab already processed and adapt matrices accordingly
      n_A <- sort(filter(dat_key, matrix == "A")$n)
      A_save <- A
      A <- A[-n_A,]
      n_B <- sort(filter(dat_key, matrix == "B")$n)
      B_save <- B
      B <- B[-n_B,]
      n_Ab <- sort(filter(dat_key, matrix == "Ab")$n)
      Ab_save <- Ab
      Ab <- Ab[-n_Ab,]
      
      # check for duplicates
      if(any( any(duplicated(n_A)) | any(duplicated(n_B)) | any(duplicated(n_Ab)) ))
        stop("The same parametrisation was detected in multiple log files! There must be something wrong, check your log files!")
      
      # matrices/vectors of results
      dat_res <- dat_key %>%
        select(-file, -!!param.IDs) %>%
        arrange(matrix, n)
      
      yA_t <- dat_res %>%
        filter(matrix == "A") %>%
        complete(n=1:nrow(A_save), matrix) %>%
        select(-n,-matrix) %>%
        as.matrix() %>% t()
      yB_t <- dat_res %>%
        filter(matrix == "B") %>%
        complete(n=1:nrow(B_save), matrix) %>%
        select(-n,-matrix) %>%
        as.matrix() %>% t()
      yAb_t <- dat_res %>%
        filter(matrix == "Ab") %>%
        complete(n=1:nrow(Ab_save), matrix) %>%
        select(-n,-matrix) %>%
        as.matrix() %>% t()
      
      if(verbose) message("[ Extracted data from log files ]")
    } # get results from log data
      
  } # restart?
  
  # do fn evaluations?
  if(do_eval) {
    
    if(ncores > 1) {
      if(verbose) {
        if(is.null(restart)) {
          message(paste("[ Evaluation of matrix A in parallel mode, involves", N, "calls of fn ]"))
        } else {
          message(paste("[ Evaluation of matrix A in parallel mode,", nrow(A), "calls of fn left ]"))
        }
      }
      if(nrow(A) > 0) {
        yA <- foreach(j=1:nrow(A), .combine = "cbind", .errorhandling = "pass") %dopar% unlist(fn(A[j,], ...))
      } else {
        yA <- NA
      }
      
      if(verbose) {
        if(is.null(restart)) {
          message(paste("[ Evaluation of matrix B in parallel mode, involves", N, "calls of fn ]"))
        } else {
          message(paste("[ Evaluation of matrix B in parallel mode,", nrow(B), "calls of fn left ]"))
        }
      }
      if(nrow(B) > 0) {
        yB <- foreach(j=1:nrow(B), .combine = "cbind", .errorhandling = "pass") %dopar% unlist(fn(B[j,], ...))
      } else {
        yB <- NA
      }
      
      if(verbose) {
        if(is.null(restart)) {
          message(paste("[ Evaluation of matrix Ab in parallel mode, involves", N, "calls of fn ]"))
        } else {
          message(paste("[ Evaluation of matrix Ab in parallel mode,", nrow(Ab), "calls of fn left ]"))
        }
      }
      if(nrow(Ab) > 0) {
        yAb <- foreach(j=1:nrow(Ab), .combine = "cbind", .errorhandling = "pass") %dopar% unlist(fn(Ab[j,], ...))
      } else {
        yAb <- NA
      }
    
    } else {
      
      if(verbose) {
        if(is.null(restart)) {
          message(paste("[ Evaluation of matrix A, involves", N, "calls of fn ]"))
        } else {
          message(paste("[ Evaluation of matrix A,", nrow(A), "calls of fn left ]"))
        }
      }
      if(nrow(A) > 0) {
        yA <- apply(A, 1, function(x) unlist(fn(x, ...)))
      } else {
        yA <- NA
      }
      
      if(verbose) {
        if(is.null(restart)) {
          message(paste("[ Evaluation of matrix B, involves", N, "calls of fn ]"))
        } else {
          message(paste("[ Evaluation of matrix B,", nrow(B), "calls of fn left ]"))
        }
      }
      if(nrow(B) > 0) {
        yB <- apply(B, 1, function(x) unlist(fn(x, ...)))
      } else {
        yB <- NA
      }
      
      if(verbose) {
        if(is.null(restart)) {
          message(paste("[ Evaluation of matrix Ab, involves", N, "calls of fn ]"))
        } else {
          message(paste("[ Evaluation of matrix Ab,", nrow(Ab), "calls of fn left ]"))
        }
      }
      if(nrow(Ab) > 0) {
        yAb <- apply(Ab, 1, function(x) unlist(fn(x, ...)))
      } else {
        yAb <- NA
      }
      
    }
    
    # if 'fn' returns single NAs, try to convert results to matrix
    y.to.mat <- function(x) {
      if(any(is.na(x))) {
        x <- rep(NA, length(param.IDs))
        names(x) <- param.IDs
      }
      return(x)
    }
    if(any(is.na(yA), is.na(yB), is.na(yAb)) & any(is.list(yA), is.list(yB), is.list(yAb))) {
      warning("Function 'fn' produced NAs! Try to calculate indices anyway but check the results!")
      yA <- sapply(yA, y.to.mat)
      yB <- sapply(yB, y.to.mat)
      yAb <- sapply(yAb, y.to.mat)
    }
  
  } # evaluations of fn necessary?
  
  # merge with data from logfiles (if restart)
  if(!is.null(restart) & do_eval) {
    if(nrow(A) > 0) yA_t[,-n_A] <- yA
    yA <- yA_t
    if(nrow(B) > 0) yB_t[,-n_B] <- yB
    yB <- yB_t
    if(nrow(Ab) > 0) yAb_t[,-n_Ab] <- yAb
    yAb <- yAb_t
    rm(yA_t, yB_t, yAb_t)
    
    # restore original parameter matrices
    A <- A_save
    B <- B_save
    Ab <- Ab_save
  }

  if(debug & do_eval) save(list = ls(all.names = TRUE), file = "vbsa_backup2.RData")

  # use current environment within internal function eval_fn()
  environment(eval_fn) <- environment()

  # evaluate and compile SA results for each output function of fn
  if(all(is.matrix(yA), is.matrix(yB), is.matrix(yAb)) & all(!is.list(yA), !is.list(yB), !is.list(yAb))) {
    if(verbose) message("[ NOTE: found multivariate output of 'fn' calls ]")
    if(ncores > 1 && nrow(yA) >= ncores) {
      out <- foreach(i=1:nrow(yA), .errorhandling = "stop") %dopar% eval_fn(yA[i,], yB[i,], yAb[i,])
    } else {
      out <- lapply(1:nrow(yA), function(i) eval_fn(yA[i,], yB[i,], yAb[i,]))
    }
    names(out) <- rownames(yA)
  } else if(all(is.vector(yA, mode = "numeric"), is.vector(yB, mode = "numeric"), is.vector(yAb, mode = "numeric"))){
    out <- eval_fn(yA, yB, yAb)
  } else stop("Outputs of 'fn' (objects yA, yB, yAb) are of unexpected types (expected all to be vectors or matrices)! Consider argument 'debug'.")

  if(anyNA(out, recursive = T)) warning("Output contains NA values! Read Notes section in the documentation for information regarding NAs!")

  if(verbose) message("[ DONE! ]")
  return(out)
}
