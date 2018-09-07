
# compute bootstraping samples of variance-based sensitivity indices based on fn output of matrices A, B and Ab
indices_vb_boot <- function(yA, yB, yAb, N, Nb, K) {
  
  # get bootstrapping realisations (i.e. sampling matrix indices) as a matrix (N,Nb)
  if(Nb > 1) {
    B <- matrix(sample.int(N, N*Nb, replace = T), nrow = N, ncol = Nb, byrow = F)
  } else {
    B <- as.matrix(1:N)
  }
  
  # use external function for computation of indices for each bootstrapping realisation
  ind_out <- apply(B, 2, function(b) indices_vb(yA[b], yB[b], yAb[b,], K) )
  
  return(list(ind_out = ind_out, B = B))
}


# compute variance-based sensitivity indices based on fn output of matrices A, B and Ab
indices_vb <- function(yA, yB, yAb, K) {
  
  # total variance, TODO: not really sure about this, in the SAFER implemantation it's var(yA)
  Vtot <- var(c(yA, yB))
  
  # first order index, Table 2 (b)
  Si <- yB * (yAb - yA)
  Si <- apply(Si, 2, mean) / Vtot
  
  # total effects index, Table 2 (f)
  St <- (yA - yAb)^2
  St <- apply(St, 2, mean) / (2*Vtot)
  
  return(list(Si = Si, St = St))
}


# function checking output of fn for non-finite values and adjusting it if desired
check_output <- function(yA, yB, yAb, A, B, Ab, N, K, subsamp, nparamsets, na.handle) {
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
  
  return(list(yA = yA, yB=yB, yAb=yAb,
              A = A, B = B, Ab = Ab,
              N = N, nparamsets = nparamsets, subsamp = subsamp))
}
