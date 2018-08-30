
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
