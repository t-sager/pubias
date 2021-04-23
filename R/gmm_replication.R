#' Title
#'
#' @param Z
#' @param sigmaZ2
#' @param symmetric
#' @param cluster_ID
#' @param cutoffs
#' @param Studynames
#'
#' @return
#' @export
#'
#' @examples
gmm_replication <- function(Z,sigmaZ2,symmetric,cluster_ID,cutoffs,Studynames) {

  # Starting Values
  Psihat0<-c(rep(1,length(cutoffs)))

  # GMM Objective Function
  GMM_obj <- function(Psi) {
    ReplicationGMMObjective(c(Psi, 1), cutoffs, symmetric, Z, sigmaZ2) +
      max(-min(Psi), 0) * 10 ^ 5
  }

  mini <- optim(par = Psihat0, fn = GMM_obj, method = "BFGS", control = list(abstol = 10 ^ -8, maxit = 10 ^ 5))

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
  se_robust=sqrt(diag(Varhat))
  dof=ncol(moments)


  if (length(cutoffs) == 1) {
    Psi_grid <- c(seq(from = 0.001, to = 5, by = 0.001), 10, 10 ^ 5)
    S_store <- zeros(length(Psi_grid), 1)

    for (m in 1:length(Psi_grid)) {
      S_store[m, 1] <- GMM_obj(Psi_grid[m])
    }

    CS_LB <- min(Psi_grid[S_store < qchisq(.95, 1)])
    CS_UB <- max(Psi_grid[S_store < qchisq(.95, 1)])

  } else if (length(cutoffs) == 2) {
    gridsteps <- 200
    Psi_grid1 <- seq(0.01,0.25,(0.25-0.01)/(gridsteps-1))
    Psi_grid2 <- seq(0.05,5,(5-0.05)/(gridsteps-1))
    S_store <- matrix(0,length(Psi_grid1),length(Psi_grid2))
    for (m1 in 1:length(Psi_grid1)) {
      for (m2 in 1:length(Psi_grid2)) {
        Psi_temp <- c(Psi_grid1[m1],Psi_grid2[m2])
        S_store[m1,m2]=GMM_obj(Psi_temp)
      }

    }

    Psi_grid <- c(t(Psi_grid1), t(Psi_grid2))

  }

  Psi_grid <- as.matrix(Psi_grid)

  # SelectionTableGMM(Psihat, se_robust, name, S_store, Psi_grid, dof)

  return(list("Psihat"= Psihat, "Varhat" = Varhat))

}
