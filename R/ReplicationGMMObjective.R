ReplicationGMMObjective <- function(betap, cutoffs, symmetric, Z, sigmaZ2){
moments <- ReplicationMoments(betap, cutoffs, symmetric, Z, sigmaZ2)
moments_mean <- mean(moments)
moments_var <- cov(moments)
moments_var[is.na(moments_var)] <- 0
objective <- nrow(moments)*moments_mean%*%pracma::pinv(moments_var)%*%t(moments_mean)
}

