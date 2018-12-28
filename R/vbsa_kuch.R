
# test example
ishigami <- function(x, a=7, b=0.1) {
  sin(x[1]) + a*(sin(x[2]))^2 + b*(x[3]^4)*sin(x[1])
}
linearsum <- function(x) sum(x)
rho <- 1
sig=2
cov_mat <- matrix(c(1  ,0  ,  0,
                    0  ,1  , rho,
                    0  ,rho, sig^2), ncol=3)

nparam <- 3
fn="linearsum"
rtype = "quasirandom"
lower=rep(-pi, nparam)
upper=rep(pi, nparam)
N=100
scrambling=0

vbsa_kuch(fn="linearsum",
          lower = rep(-pi, nparam), upper = rep(pi, nparam),
          N = 10000, rtype = "quasirandom",
          cov_mat = cov_mat, mu = c(0,0,0))

Si_a <- c(1/(2+sig^2+2*sig*rho),
          (1+rho*sig)^2 / (2+sig^2+2*sig*rho),
          (sig+rho)^2 / (2+sig^2+2*sig*rho))
St_a <- c(1 / (2+sig^2+2*sig*rho),
          (1-rho^2) / (2+sig^2+2*rho*sig),
          (sig^2*(1-rho^2)) / (2+sig^2+2*sig*rho))


vbsa_kuch <- function(
  fn,
  ...,
  lower=-Inf,
  upper=Inf,
  N=1000,
  rtype=c("quasirandom", "lhs"),
  scrambling = 0,
  mu = NULL,
  cov_mat
) {
  
  # Assigning the objective function
  fn <- match.fun(fn)

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
  
  # Backing up the original boundaries
  lower.ini <- lower
  upper.ini <- upper
  X.Boundaries.ini <- X.Boundaries
  LOWER.ini <- matrix( rep(lower.ini, N), nrow=N, byrow=TRUE)
  UPPER.ini <- matrix( rep(upper.ini, N), nrow=N, byrow=TRUE)
  if(is.null(mu))
    mu <- lower.ini+0.5*(upper.ini-lower.ini)
  
  # normalising
  lower <- rep(0, K)
  upper <- rep(1, K)
  X.Boundaries <- cbind(lower, upper)
  rownames(X.Boundaries) <- param.IDs
  
  # parameter sampling from U[0,1] (u and u' in Kucherenko)
  if (rtype=="quasirandom") {
    samp <- randtoolbox::sobol(n=N, dim=2*K, scrambling=scrambling)
  } else if (rtype=="lhs") {
    samp <- rbind(randomLHS(n=N, k = K), randomLHS(n=N, k = K))
  } # ELSE end
  
  # transfer samples to standard normal space N[0,1]
  A_norm <- t(apply(samp[,1:K], 1, qnorm)) # u in Kucherenko
  B_norm <- t(apply(samp[,(K+1):(2*K)], 1, qnorm)) # u' in Kucherenko
  
  # Cholesky decomposition of correlation matrix
  cov_chol <- suppressWarnings(chol(cov_mat, pivot = T))
  
  # N[0,1] (u in Kucherenko) to N[mu,cov] (x)
  xA <- t(apply(A_norm, 1, function(x) crossprod(cov_chol, x) + mu)) # (y,z) in Kucherenko
  xB <- t(apply(B_norm, 1, function(x) crossprod(cov_chol, x) + mu)) # (y',z') in Kucherenko
  
  # get samples from conditional distributions (correlation among input factors)
  #xA_zy <- NULL
  xA_yz <- NULL
  xB_zy <- NULL
  for (j in 1:K) {
    # calculate (y, z_bar_prime)
    #xA_zy_t <- t(sapply(1:N, function(i) cond_distr(xA[i,], B_norm[i,], j, mu, cov_mat)))
    #xA_zy <- rbind(xA_zy, xA_zy_t)
    
    # calculate (y_bar_prime, z)
    j_c <- (1:K)[-j]
    xA_yz_t<- t(sapply(1:N, function(i) cond_distr(xA[i,], B_norm[i,], j_c, mu, cov_mat)))
    xA_yz <- rbind(xA_yz, xA_yz_t)
    
    # calculate (y_prime, z_hat)
    xB_zy_t<- t(sapply(1:N, function(i) cond_distr(xB[i,], A_norm[i,], j, mu, cov_mat)))
    xB_zy <- rbind(xB_zy, xB_zy_t)
  }
  
  # evaluation of fn
  fnA <- apply(xA, 1, function(x) fn(x, ...))
  fnB <- apply(xB, 1, function(x) fn(x, ...))
  #fnA_zy <- apply(xA_zy, 1, function(x) fn(x, ...))
  fnA_yz <- apply(xA_yz, 1, function(x) fn(x, ...))
  fnB_zy <- apply(xB_zy, 1, function(x) fn(x, ...))
  
  #fnA_zy <- matrix(fnA_zy, ncol=3)
  fnA_yz <- matrix(fnA_yz, ncol=3)
  fnB_zy <- matrix(fnB_zy, ncol=3)
  
  # calculate indices
  f0 <- mean(fnA)
  D <- mean(fnA^2) - f0^2
  
  #Si <- apply(fnA * (fnA_zy-fnB), 2, mean) / D
  Si <- apply(fnB * (fnB_zy-fnA), 2, mean) / D
  St <- apply( (fnA - fnA_yz)^2, 2, mean) / (2*D)
  
  return(list(Si = Si, St = St))
}
