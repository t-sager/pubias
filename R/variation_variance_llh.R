variation_variance_llh <- function(lambdabar, tauhat, betap, cutoffs, symmetric, X, sigma, C) {
  n <- nrow(X)
  betap <- t(betap)
  TT <- X/sigma

  # Tpowers
  Tpowers <-  zeros(n,length(cutoffs)+1)

  if (symmetric == TRUE) {

    Tpowers[,1] <- abs(TT)<cutoffs[1]


     if (length(cutoffs)>1) {
         for (m in (2:length(cutoffs))) {
              Tpowers[,m]=(abs(TT)<cutoffs[m])*(abs(TT)>=cutoffs[m-1]);
          }
          Tpowers[,length(cutoffs)+1]=abs(TT)>=cutoffs[length(cutoffs)];
     }
    Tpowers[,ncol(Tpowers)]=abs(TT)>=cutoffs[length(cutoffs)];
    } else {

        Tpowers[,1] <- TT<cutoffs[1];


    if (length(cutoffs)>1) {
          for (m in 2:length(cutoffs)) {
              Tpowers[,m] <- (TT<cutoffs[m])*(TT>=cutoffs[m-1]);
          }


    }
        Tpowers[,length(cutoffs)+1] <- (TT)>=cutoffs[length(cutoffs)];
  }





# Calculate objective function LLH

# Calculating components of the likelihood
# vector of estimated publication probabilities

phat <- zeros(nrow(Tpowers), 1)

for (m in (1:ncol(Tpowers))) {
  Cmat <- repmat(Tpowers[ ,m],1,ncol(as.matrix(C)))*C
  phat <- phat + Cmat * betap[,m]

}

#vector of un-truncated likelihoods
if (symmetric == TRUE){
    #Monte-carlo integration
    set.seed(1)
    draw <- matrix(runif(10^5),1)
    theta_vec <- qgamma(draw, lambdabar,scale = tauhat)
    theta_mat <- matrix(theta_vec,n,length(theta_vec))
    X_mat <- repmat(X,1,length(theta_vec))
    sigma_mat <- repmat(sigma,1,length(theta_vec))
    g <-  (0.5*dnorm((X_mat-theta_mat)/sigma_mat)+0.5*dnorm((-X_mat-theta_mat)/sigma_mat))/sigma_mat # TSA: CHECK!
    fX <- rowMeans(g,2)
} else {
  fX <- dnorm(X,lambdabar, sqrt(sigma^2 + tauhat^2))
}


#normalizing constant
mu_vec <- lambdabar*ones(n,1)
sigma_tilde_vec <- sqrt(((tauhat)^2 +sigma^2))
prob_vec <- zeros(n,length(cutoffs)+1)
if (symmetric == TRUE){
  for (m in (1:length(cutoffs))){

      #Monte Carlo Integration
      g <- (pnorm(cutoffs[m]-theta_mat/sigma_mat)-pnorm(-cutoffs[m]-theta_mat/sigma_mat))
      prob_vec[,m+1] <- rowMeans(g)

  }
  prob_vec <- cbind(prob_vec,1)
  mean_Z1 <- prob_vec[,2:ncol(prob_vec)]-prob_vec[,1:ncol(prob_vec)-1]

} else {
  for (m in (1:length(cutoffs))){
    prob_vec[,m+1] <- pnorm((cutoffs[m]*sigma-mu_vec)/sigma_tilde_vec)
  }
  prob_vec <- cbind(prob_vec,1)
  mean_Z1 <- prob_vec[,2:ncol(prob_vec)]-prob_vec[,1:ncol(prob_vec)-1]
  }

normalizingconst <- zeros(n,1)
parameter_space_violation <- 0
for (m in (1:ncol(mean_Z1))){
  Cmat <- C*repmat(mean_Z1[,m],1,ncol(as.matrix(C)))
  normalizingconst <- normalizingconst+Cmat*betap[,m]
}

#vector of likelihoods
L <- phat*fX/normalizingconst
logL <- log(phat)+log(fX)-log(normalizingconst)

LLH <- -sum(log(L)) #objective function; note the sign flip, since minimization

return(list("LLH"= LLH, "logL" = logL))
}





