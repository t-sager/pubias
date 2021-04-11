VariationVarianceLogLikelihoodtdist <- function(lambdabar, tauhat, nuhat, betap, cutoffs, symmetric, X, sigma,C,numerical_integration) {
  
  n=length(X);
  
  #%regressors for step function p
  TT=X/sigma;
  
  
  # Tpowers
  Tpowers=matrix (0,n,length(cutoffs)+1);
  
  if (symmetric==1 ) {
    
    Tpowers[,1]=abs(TT)<cutoffs[1];
    
    
    if (length(cutoffs)>1) {
      for (m in (2:length(cutoffs))) {
        Tpowers[,m]=(abs(TT)<cutoffs[m])*(abs(TT)>=cutoffs[m-1]);
      }
      Tpowers[,length(cutoffs)+1]=abs(TT)>=cutoffs[length(cutoffs)];  
    }
    Tpowers[,length(cutoffs)+1]=abs(TT)>=cutoffs[length(cutoffs)];
  } else {
    
    Tpowers[,1]=TT<cutoffs[1];      
    
    
    if (length(cutoffs)>1) {
      for (m in 2:length(cutoffs)) {
        Tpowers[,m]=(TT<cutoffs[m])*(TT>=cutoffs[m-1]);  
      }
      
      
    }
    Tpowers[,length(cutoffs)+1]=(TT)>=cutoffs[length(cutoffs)]; 
  }
  
  
  
  
  
  # Calculate objective function LLH
  
  # Calculating components of the likelihood
  # vector of estimated publication probabilities
  
  phat <- matrix (0, dim(Tpowers)[1], 1);
  
  for (m in (1:dim(Tpowers)[2])) {
    Cmat <- matrix(rep(Tpowers[, m], dim(C)[2]), ncol=dim(C)[2], byrow = FALSE) * C  
    phat <- phat + Cmat %*% betap[m]
    
  }
  
  #vector of un-truncated likelihoods
    if (numerical_integration==1){
      #likelihoods calculated by numerical integration, t distribution
      g <- function(theta) {dnorm((X-theta)/sigma)/sigma*dt((theta-lambdabar)/tauhat,nuhat)/tauhat}
      fX <- integrate(g,-Inf,Inf, subdivisions=2000)$value # numerically solving the function
      
    } else {
      
      #Monte-carlo integration
      set.seed(1)
      draw <- matrix(runif(10^5),1)
      theta_vec=tinv(draw,nuhat)*tauhat+lambdabar;
      theta_mat=matrix(rep(theta_vec,n),ncol=length(theta_vec),byrow = FALSE);
      X_mat=matrix(rep(X,length(theta_vec)),ncol=length(theta_vec),byrow = FALSE);
      sigma_mat=matrix(rep(sigma,length(theta_vec)),ncol=length(theta_vec),byrow = FALSE);
      g <-  dnorm((X_mat-theta_mat)/sigma_mat)/sigma_mat
      fX <- mean(g,2)
      
    }
  
  #normalizing constant
  prob_vec <- zeros(n,length(cutoffs)+1)
  
  #Normalizing constant, t distribution
  for (m in (1:length(cutoffs))){ 
      if (numerical_integration==1){
        g <- function(theta) {dnorm(cutoffs(m)-theta/sigma)*dt((theta-lambdabar)/tauhat,nuhat)/tauhat}
        prob_vec[,m+1] <- integrate(g,-Inf,Inf, subdivisions=2000)$value
      } else {
        #Monte Carlo Integration
        g <- dnorm(cutoffs(m)-theta_mat/sigma_mat)
        prob_vec[,m+1] <- mean(g,2)
      }
  }
prob_vec <- cbind(prob_vec,1)
mean_Z1 <- prob_vec[,2:ncol(prob_vec)]-prob_vec[,1:ncol(prob_vec)-1]

  
normalizingconst <- zeros(n,1)
for (m in (1:dim(mean_Z1)[2])){
    normalizingconst <- normalizingconst+mean_Z1[,m]*C%*%betap[m] # @TSA: Korrekt mit betap?
  }
  
#vector of likelihoods
L <- phat*fX/normalizingconst
logL <- log(L)
  
LLH <- -sum(log(L)) #objective function; note the sign flip, since minimization
}

