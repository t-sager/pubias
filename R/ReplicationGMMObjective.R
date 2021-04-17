ReplicationGMMObjective <- function(betap, cutoffs, symmetric, Z, sigmaZ2,momentchoice){
#Calculates continiously updating GMM objective
moments <- ReplicationMoments(betap, cutoffs, symmetric, Z, sigmaZ2,momentchoice)
moments_mean <- mean(moments)
moments_var <- cov(moments)
moments_var[is.na(moments_var)] <- 0
objective <- nrow(moments)%*%moments_mean%*%pinv(moments_var)%*%t(moments_mean)
}

