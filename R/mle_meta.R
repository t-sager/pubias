#' Computing the publication probability in meta analyses
#'
#'`mle_meta()`  calculates the publication probability, its variance and robust standard errors
#' of replication studies by a Maximum Likelihood approach.  Check package vignette for further information.
#'
#' @param X A `n x 1` matrix containing the estimates, where `n` is the number of estimates.
#' @param sigma A `n x 1` matrix containing the standard errors of the estimates, where `n` is the number of estimates.
#' @param symmetric If set to TRUE, the publication probability is assumed to be symmetric around zero. If set to FALSE, asymmetry is allowed.
#' @param cluster_ID A `n x 1` matrix containing IDs going from 1 to `n`, where `n` is the number of estimates.
#' @param cutoffs A matrix containing the thresholds for the steps of the publication probability. Should be strictly increasing column
#' vector of size `k x 1` where `k` is the number of cutoffs.
#' @param C A `n x 1` matrix with all values being 1. Controls.
#'
#' @return Returns a list object with the publication probability (`Psihat`), its variance (`Varhat`) and robust standard error (`se_robust`).
#' @export
#'
mle_meta <- function(X, sigma, symmetric, cluster_ID, cutoffs, C) {

  # Stepsize
  stepsize <- 10 ^ (-6)
  # n
  n <- length(X)

  # Run model
  LLH <-
    function (Psi) {
      f <- variation_variance_llh(Psi[1],
                                  Psi[2],
                                  c(reshape(t(Psi[-c(1, 2)]), c(
                                    length(Psi[-c(1, 2)]) / length(cutoffs), length(cutoffs)
                                  )), 1),
                                  cutoffs,
                                  symmetric,
                                  X,
                                  sigma,
                                  C)
      return(f)
    }


  nn <- n

  # Necessary for optimization procedure
  LLH_only <- function (Psi) {
    f <- LLH(Psi)

    return(f$LLH)
  }

  #########################
  # Get maximum likelihood estimator using only LLH

  # Starting Values
  Psihat0 <- c(1, 1, rep(1, length(cutoffs)))

  # Optimize based on starting values
  mini <-
    optim(
      par = Psihat0,
      fn = LLH_only,
      method = "BFGS",
      control = list(abstol = 10 ^ -8, maxit = 10 ^ 5)
    )

  # More accurate Optimization, aka. optimizing again
  Psihat1 <- mini$par
  mini <-
    optim(
      par = Psihat1,
      fn = LLH_only,
      method = "BFGS",
      control = list(abstol = 10 ^ -8, maxit = 10 ^ 5)
    )

  # Optimal values, Objval = max. likelihood
  Psihat <<- mini$par
  Objval <- mini$value

  # Variance and robust SEs
  Varhat <- robust_variance(stepsize, nn, Psihat, LLH, cluster_ID)
  se_robust <- sqrt(diag(Varhat))

  # Return estimation results
  return(list(
    "Psihat" = Psihat,
    "Varhat" = Varhat,
    "se_robust" = se_robust
  ))
}
