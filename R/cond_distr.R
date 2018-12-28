#' Transfrom multivariate standard normal to conditional distribution function
#' 
#' @param z vector of length K; samples of K parameters in N(mu,cov)
#' @param q_n vector of length K; alternative parameter samples in N(0,1) to be transformed
#' @param s vector of length d defining the subset dividing z and q_n
#' @param mu mean values of the input parameter distributions
#' @param cov_mat covariance matrix describing dependencies of the parameters
#' 
#' @return Vector of length K; values of q_n from N(0,1) space transformed to
#' conditional normal distribution (derived from subsets s and K-s of z) respecting
#' the given mu and cov_mat.
#' 
#' @author Tobias Pilz \email{tpilz@@uni-potsdam.de}
#' 
cond_distr <- function(z, q_n, s, mu, cov_mat) {
  
  # subsets
  k <- length(z)
  s_c <- (1:k)[-s]
  x <- z[s]
  y <- z[s_c]
  mu_x <- mu[s]
  mu_y <- mu[s_c]
  c_x <- cov_mat[s,s, drop=F]
  tryCatch({c_xi <- solve(c_x)},
           error = function(e) {
             mat_str <- NULL
             for(i in 1:nrow(c_x))
               mat_str <- paste(mat_str, paste(c_x[i,], collapse="  "), "\n")
             stop(paste0("A problem occurred while inverting matrix\n",
                   mat_str, "\n",
                   "Error message: ", e))})
  c_y <- cov_mat[s_c,s_c, drop=F]
  c_xy <- cov_mat[s,s_c, drop=F]
  c_yx <- cov_mat[s_c,s, drop=F]
  
  # conditional mean E(y|x)
  mu_yc <- mu_y + c_yx %*% c_xi %*% (x - mu_x)
  
  # conditional covariance matrix var(y|x)
  c_yc <- c_y - c_yx %*% c_xi %*% c_xy
  
  # Cholesky decomposition
  c_yc_chol <- suppressWarnings(chol(c_yc, pivot = T))
  
  # realisation of (y|x) following N[mu_yc,c_yc] from sampled y in N[0,1]
  yc <- crossprod(c_yc_chol, q_n[s_c]) + mu_yc
  
  # merge with x
  out <- c(x,yc)
  
  return(out)
}
