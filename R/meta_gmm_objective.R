meta_gmm_objective <- function(betap,cutoffs,symmetric,X,sigma,cluster_ID) {
  mom <- meta_moments(betap,cutoffs,symmetric,X,sigma)
  moments_mean <- mom$moment_mean
  rhat <- mom$raw_moments
  Sigma_hat <- clustered_covariance_estimate(rhat,cluster_ID)
  objective=length(X) * moments_mean %*% solve(Sigma_hat) %*% t(moments_mean)
  return(objective)
}
