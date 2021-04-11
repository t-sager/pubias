MetastudyGMMObjective <- function(beta_theta,betap,cutoffs,symmetric,X,sigma,cluster_ID) {
  betap=c(beta_theta,betap)
  betap=as.matrix(betap)
  betap=t(betap)
  betap=cbind(betap,1)
  #Calculate continuously updating GMM objective
  mom=ComputingMetastudyMoments(betap,cutoffs,symmetric,X,sigma);
  moments_mean=mom$moment_mean;
  rhat=mom$raw_moments;
  Sigma_hat=Clustered_covariance_estimate(rhat,cluster_ID);
  objective=length(X)*moments_mean%*%solve(Sigma_hat)%*%t(moments_mean);
  return(objective)
}
