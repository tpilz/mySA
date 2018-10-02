# function for checking of fn output, calculation of sensitivity indices, and (pre-)compilation of output
eval_fn <- function(yA, yB, yAb) {

  # check (and optionally adjust) output of fn
  checked <- check_output(yA, yB, yAb, A, B, Ab, N, K, subsamp, nparamsets, na.handle)
  yA <- checked$yA
  yB <- checked$yB
  yAb <- checked$yAb
  A <- checked$A
  B <- checked$B
  Ab <- checked$Ab
  N <- checked$N
  subsamp <- checked$subsamp
  nparamsets <- checked$nparamsets

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
    ranks <- t(apply(out$Si, 1, function(x)  match(param.IDs, param.IDs[order(x, decreasing = T)])))
    colnames(ranks) <- param.IDs
  } else {
    ranks <- match(param.IDs, param.IDs[order(out$Si, decreasing = T)])
    names(ranks) <- param.IDs
  }
  out <- c(out, list(ranking = ranks))

  # Bootstrapping output
  if(Nb > 1) {
    n_min <- max(1, round(Nb * Nb.sig/2))
    n_max <- round(Nb * (1 - Nb.sig/2) )
    if(length(dim(Si)) == 3) {
      Si_sorted <- aperm(apply(Si, c(1,3), sort, na.last = F), c(2,1,3))
      Si.lo <- Si_sorted[,n_min,]
      Si.up <- Si_sorted[,n_max,]
      St_sorted <- aperm(apply(St, c(1,3), sort, na.last = F), c(2,1,3))
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
      ranks.boot <- aperm(apply(Si, c(1,2), function(x)  match(param.IDs, param.IDs[order(x, decreasing = T)])), c(2,3,1))
      dimnames(ranks.boot)[[3]] <- param.IDs
    } else {
      ranks.boot <- t(apply(Si, 1, function(x)  match(param.IDs, param.IDs[order(x, decreasing = T)])))
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

  return(out)
}

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

  if(Vtot == 0) return(list(Si = rep(NA, K), St = rep(NA, K)))

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
      if(!is.null(subsamp) && any(subsamp > N)) subsamp <- subsamp[-which(subsamp > N)]
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
