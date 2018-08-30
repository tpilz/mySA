
# compute variance-based sensitivity indices based on fn output of matrices A, B and Ab
indices_vb <- function(yA, yB, yAb, K) {
  
  # yAb as a matrix for more efficient calculations
  #yAb_t <- matrix(yAb, ncol = K, byrow = T)
  
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
