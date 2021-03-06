#' Computing the publication probability in replication studies
#'
#'`mle_replication()` calculates the publication probability, its variance and robust standard errors
#' of replication studies by a GMM approach.  Check package vignette for further information.
#'
#' @param Z A `n x 2` matrix where the first (second) column contains the standardized original estimates (replication estimates), where `n` is the number of estimates.
#' @param sigmaZ2 A `n x 1` matrix containing the standard errors (se_replication divided by se_original) of the estimates, where `n` is the number of estimates.
#' @param symmetric If set to TRUE, the publication probability is assumed to be symmetric around zero. If set to FALSE, asymmetry is allowed.
#' @param cluster_ID A `n x 1` matrix containing IDs going from 1 to `n`, where `n` is the number of estimates.
#' @param cutoffs A matrix containing the thresholds for the steps of the publication probability. Should be strictly increasing column
#' vector of size `k x 1` where `k` is the number of cutoffs.
#' @param C A `n x 1` matrix with all values being 1. Controls.
#'
#' @return Returns a list object with the publication probability (`Psihat`), its variance (`Varhat`) and robust standard error (`se_robust`).
#' @export
#'
mle_replication <- function(Z, sigmaZ2, symmetric, cluster_ID, cutoffs, C) {

  # Stepsize
  stepsize <- 10 ^ (-6)
  # n
  n <- nrow(Z)

  # Run model
  LLH <- function(Psi) {
    f <- replication_analytic_llh(Psi[1],
                                  Psi[2],
                                  c(reshape(t(Psi[-c(1, 2)]), c(
                                    length(Psi[-c(1, 2)]) / length(cutoffs), length(cutoffs)
                                  )), 1),
                                  cutoffs,
                                  symmetric,
                                  Z,
                                  sigmaZ2,
                                  C)
    return(f)
  }


  nn = n
  # Necessary for optimization procedure
  LLH_only <- function (Psi) {
    f <- LLH(Psi)

    return(f$LLH)
  }

  #########################
  # Get maximum likelihood estimator using only LLH

  # Starting values
  Psihat0 <- c(1, 1, rep(1, length(cutoffs)))

  upper.b = rep(Inf, length(Psihat0))
  lower.b = c(-Inf, -Inf,-Inf)

  # Optimize based on starting values
  mini <-
    suppressWarnings(nlminb(
      objective = LLH_only,
      start = Psihat0,
      lower = lower.b,
      upper = upper.b
    )
    )

  # More accurate Optimization, aka. optimizing again
  Psihat1 <- mini$par
  mini <-
    suppressWarnings(nlminb(
      objective = LLH_only,
      start = Psihat1,
      lower = lower.b,
      upper = upper.b
    )
    )

  # Optimal values, Objval = max. likelihood
  Psihat <<- mini$par
  Objval <- mini$value

  # Variance and robust SEs
  Varhat <- robust_variance(stepsize, nn, Psihat, LLH, cluster_ID)
  se_robust <- sqrt(diag(Varhat))

  return(list(
    "Psihat" = Psihat,
    "Varhat" = Varhat,
    "se_robust" = se_robust
  ))
}
