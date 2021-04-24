#' Title
#'
#' @param X
#' @param sigma
#' @param symmetric
#' @param cluster_ID
#' @param cutoffs
#' @param Studynames
#'
#' @return
#' @export
#'
#' @examples
gmm_meta <- function(X, sigma, symmetric, cluster_ID, cutoffs, studynames) {

  # Starting Values
  Psihat0<-c(rep(1,length(cutoffs)))

  # GMM Objective Function
  GMM_obj <- function(Psi) {
    meta_gmm_objective(c(Psi, 1), cutoffs, symmetric, X, sigma, cluster_ID) +
      max(-min(Psi), 0) * 10 ^ 5
  }

  mini<-optim(par=Psihat0,fn=GMM_obj,method="BFGS",control = list(abstol=10^-8,maxit=10^5))

  Psihat <- mini$par
  Objval <- mini$value

  mom = meta_moments(c(Psihat, 1), cutoffs, symmetric, X, sigma)
  moments = mom$moment_mean
  Sigma_temp = mom$moment_var
  rhat = mom$raw_moments

  Sigma_hat <- clustered_covariance_estimate(rhat, cluster_ID)

  stepsize = 10 ^ -3
  G = matrix(0, dim(moments)[2], length(Psihat))

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


  Varhat <- solve(t(G) %*% solve(Sigma_hat) %*% G) / length(X)
  se_robust <- sqrt(diag(Varhat))
  dof <- ncol(moments)

  if (length(cutoffs) == 1) {
    Psi_grid = c(seq(0.001, 5, 0.001), 10, 10 ^ 5)
    S_store = matrix(0, length(Psi_grid), 1)
    for (m in c(1:length(Psi_grid))) {
      S_store[m, 1] = GMM_obj(Psi_grid[m])
    }

    CS_LB = min(Psi_grid[S_store < qchisq(0.95, 1)])
    CS_UB = max(Psi_grid[S_store < qchisq(0.95, 1)])
  } else if (length(cutoffs) == 2) {
    gridsteps = 200

    Psi_grid1 = seq(0.01, 0.25, (0.25 - 0.01) / (gridsteps - 1))

    Psi_grid2 = seq(0.05, 5, (5 - 0.05) / (gridsteps - 1))

    S_store = matrix(0, length(Psi_grid1), length(Psi_grid2))

    for (m1 in c(1:length(Psi_grid1))) {
      for (m2 in c(1:length(Psi_grid2))) {
        Psi_temp = c(Psi_grid1[m1], Psi_grid2[m2])

        S_store[m1, m2] = GMM_obj(Psi_temp)

      }

    }
    Psi_grid <- c(t(Psi_grid1), t(Psi_grid2))

  } else {
    Psi_grid <- c(seq(0.001, 5, 0.001), 10, 10 ^ 5)
    S_store <- zeros(length(Psi_grid), 1)
  }

  Psi_grid <- as.matrix(Psi_grid)

  return(list("Psihat"= Psihat, "Varhat" = Varhat, "se_robust" = se_robust))
}
