#' Computing the publication probability in replication studies
#'
#'`gmm_replication()`  calculates the publication probability, its variance and robust standard errors
#' of meta-analyses by a GMM approach.  Check package vignette for further information..
#'
#' @param Z A `n x 2` matrix where the first (second) column contains the standardized original estimates (replication estimates), where `n` is the number of estimates.
#' @param sigmaZ2 A `n x 1` matrix containing the standard errors (se_replication divided by se_original) of the estimates, where `n` is the number of estimates.
#' @param symmetric If set to `1`, the publication probability is assumed to be symmetric around zero. If set to `0`, asymmetry is allowed.
#' @param cluster_ID A `n x 1` matrix containing IDs going from 1 to `n`, where `n` is the number of estimates.
#' @param cutoffs A matrix containing the thresholds for the steps of the publication probability. Should be strictly increasing column
#' vector of size `k x 1` where `k` is the number of cutoffs.
#' @param studynames A vector of type `character` containing all the Studynames of size `n` in the same order as the argument `data`.
#'
#' @return Returns a list object with the publication probability (`Psihat`), its variance (`Varhat`) and robust standard errors (`se_robust`).
#' @export
#'
gmm_replication <- function(Z, sigmaZ2, symmetric, cluster_ID, cutoffs, studynames) {

  # Starting Values
  Psihat0<-c(rep(0.1,length(cutoffs)))

  # GMM Objective Function
  GMM_obj <- function(Psi) {
    replication_gmm_objective(c(Psi, 1), cutoffs, symmetric, Z, sigmaZ2) +
      max(-min(Psi), 0) * 10 ^ 5
  }

  mini <- optim(par = Psihat0, fn = GMM_obj, method = "BFGS", control = list(abstol = 10 ^ -8, maxit = 10 ^ 5))

  # More accurate Optimization:
  Psihat1 <- mini$par

  mini<-optim(par=Psihat1,fn=GMM_obj)

  Psihat <- mini$par
  Objval <- mini$value

  moments <-replication_moments(c(Psihat, 1), cutoffs, symmetric, Z, sigmaZ2)
  Sigma_hat <- cov(moments)
  stepsize <- 10 ^ -3
  G <- zeros(ncol(moments), length(Psihat))

  for (n1 in 1:length(Psihat)) {
    beta_plus <- Psihat

    beta_plus[n1] <- beta_plus[n1] + stepsize

    moments_plus <-
      replication_moments(c(beta_plus, 1),
                         cutoffs,
                         symmetric,
                         Z,
                         sigmaZ2)

    beta_minus <- Psihat

    beta_minus[n1] <- beta_minus[n1] - stepsize

    moments_minus <-
      replication_moments(c(beta_minus, 1),
                         cutoffs,
                         symmetric,
                         Z,
                         sigmaZ2)

    G[, n1] = t(mean((moments_plus - moments_minus) / (2 * stepsize)))
  }


  Varhat <- solve(t(G)%*%solve(Sigma_hat)%*%G) / nrow(moments)
  se_robust <- sqrt(diag(Varhat))

  return(list("Psihat"= Psihat, "Varhat" = Varhat, "se_robust" = se_robust))

}
