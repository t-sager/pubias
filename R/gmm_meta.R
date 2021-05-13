#' Computing the publication probability in meta analyses
#'
#'`gmm_meta()` calculates the publication probability, its variance and robust standard errors
#' of meta-analyses by a Maximum Likelihood approach. Check package vignette for further information.
#'
#' @param X A `n x 1` matrix containing the estimates, where `n` is the number of estimates.
#' @param sigma A `n x 1` matrix containing the standard errors of the estimates, where `n` is the number of estimates.
#' @param symmetric If set to TRUE, the publication probability is assumed to be symmetric around zero. If set to FALSE, asymmetry is allowed.
#' @param cluster_ID A `n x 1` matrix containing IDs going from 1 to `n`, where `n` is the number of estimates.
#' @param cutoffs A matrix containing the thresholds for the steps of the publication probability. Should be strictly increasing column
#' vector of size `k x 1` where `k` is the number of cutoffs.
#' @param studynames A vector of type `character` containing all the Studynames of size `n` in the same order as the argument `data`.
#'
#' @return Returns a list object with the publication probability (`Psihat`), its variance (`Varhat`) and robust standard error (`se_robust`).
#' @export
#'
gmm_meta <- function(X, sigma, symmetric, cluster_ID, cutoffs, studynames) {

  # Starting Values
  Psihat0<-c(rep(1,length(cutoffs)))

  # GMM Objective Function
  GMM_obj <- function(Psi) {
    meta_gmm_objective(c(Psi, 1), cutoffs, symmetric, X, sigma, cluster_ID) + max(-min(Psi), 0) * 10 ^ 5
  }

  # Optimizing the objective function w.r.t starting values
  mini<-optim(par=Psihat0,fn=GMM_obj,method="BFGS",control = list(abstol=10^-8,maxit=10^5))

  # More accurate Optimization (optimizing again, based on first optimization)
  Psihat1 <- mini$par
  mini<-optim(par=Psihat1,fn=GMM_obj,method="BFGS",control = list(abstol=10^-8,maxit=10^5))

  # Optimal Values
  Psihat <- mini$par
  Objval <- mini$value

  # Calculate moments with optimal publication probability
  mom <- meta_moments(c(Psihat, 1), cutoffs, symmetric, X, sigma)
  moments <- mom$moment_mean
  Sigma_temp <- mom$moment_var
  rhat <- mom$raw_moments

  Sigma_hat <- clustered_covariance_estimate(rhat, cluster_ID)

  stepsize <- 10 ^ -3
  G <- zeros(ncol(moments), length(Psihat))

  for (n1 in 1:length(Psihat)) {
    beta_plus = Psihat
    beta_plus[n1] = beta_plus[n1] + stepsize
    mom = meta_moments(c(beta_plus, 1), cutoffs, symmetric, X, sigma)
    moments_plus = mom$moment_mean
    beta_minus = Psihat
    beta_minus[n1] = beta_minus[n1] - stepsize
    mom = meta_moments(c(beta_minus, 1), cutoffs, symmetric, X, sigma)
    moments_minus = mom$moment_mean
    G[, n1] = t(moments_plus - moments_minus) / (2 * stepsize)
  }

  # Variance and robust SEs
  Varhat <- solve(t(G) %*% solve(Sigma_hat) %*% G) / length(X)
  se_robust <- sqrt(diag(Varhat))
  dof <- ncol(moments)

  # Returning the results
  return(list("Psihat"= Psihat, "Varhat" = Varhat, "se_robust" = se_robust))
}
