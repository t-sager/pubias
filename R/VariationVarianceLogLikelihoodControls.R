VariationVarianceLogLikelihoodControls <- function(lambdabar, tauhat, betap, cutoffs, symmetric, X, sigma,C,numerical_integration) {

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
if (symmetric==1){
  
  if (numerical_integration==1){
    #likelihoods calculated by numerical integration, gamma distribution
    g <- function(theta) {(0.5*dnorm((X-theta)/sigma)+0.5*dnorm((-X-theta)/sigma))/sigma*dgamma(theta,lambdabar,tauhat)}
    fX <- integrate(g,-Inf,Inf, subdivisions=2000)$value # numerically solving the function
    
  } else {
    
    #Monte-carlo integration
    set.seed(1)
    draw <- matrix(runif(10^5),1)
    theta_vec=dinvgamma(draw,lambdabar,tauhat);
    theta_mat=matrix(rep(theta_vec,n),ncol=length(theta_vec),byrow = FALSE);
    X_mat=matrix(rep(X,length(theta_vec)),ncol=length(theta_vec),byrow = FALSE);
    sigma_mat=matrix(rep(sigma,length(theta_vec)),ncol=length(theta_vec),byrow = FALSE);
    g <-  (0.5*((X_mat-theta_mat)/sigma_mat)+0.5*dnorm((-X_mat-theta_mat)/sigma_mat))/sigma_mat
    fX <- mean(g,2)
  
    }
  
} else {
  fX <- dnorm(X,lambdabar, sqrt(sigma^2 + tauhat^2))
}


#normalizing constant
mu_vec <- lambdabar*ones(n,1)
sigma_tilde_vec <- sqrt(((tauhat)^2 +sigma^2))
prob_vec <- zeros(n,length(cutoffs)+1)
if (symmetric==1){
  for (m in (1:length(cutoffs))){

    if (numerical_integration==1){
      #Normalizing constant, gamma distribution
      g <- function(theta) {(dnorm(cutoffs[m]-theta/sigma)-dnorm(-cutoffs[m]-theta/sigma))*dgamma(theta,lambdabar,tauhat)}
      prob_vec[,m+1] <- integrate(g,-Inf,Inf, subdivisions=2000)$value
    } else {
      #Monte Carlo Integration
      g <- (dnorm(cutoffs[m]-theta_mat/sigma_mat)-dnorm(-cutoffs[m]-theta_mat/sigma_mat))
      prob_vec[,m+1] <- mean(g,2)
    }
  }
  prob_vec <- cbind(prob_vec,1)
  mean_Z1 <- prob_vec[,2:ncol(prob_vec)]-prob_vec[,1:ncol(prob_vec)-1]
} else {
  for (m in (1:length(cutoffs))){
    prob_vec[,m+1] <- dnorm((cutoffs[m]*sigma-mu_vec)/sigma_tilde_vec)
  }
  prob_vec <- cbind(prob_vec,1)
  mean_Z1 <- prob_vec[,2:ncol(prob_vec)]-prob_vec[,1:ncol(prob_vec)-1]
  }

normalizingconst <- zeros(n,1)
parameter_space_violation <- 0
for (m in (1:dim(mean_Z1)[2])){
  Cmat <- C*matrix(rep(mean_Z1[,m],dim(C)[2]), ncol = dim(C)[2],byrow = FALSE)
  normalizingconst <- normalizingconst+Cmat*betap[m]
  if (min(Cmat*betap[m] < 0)){
    parameter_space_violation <- 1
  }
}

#vector of likelihoods
L <- phat*fX/normalizingconst
logL <- log(phat)+log(fX)-log(normalizingconst)

LLH <- -sum(log(L)) #objective function; note the sign flip, since minimization
if (parameter_space_violation==1){
  LLH <- 10^5
}
}



# 
# test<-TRUE
# if (test) {
#   C <- matrix(1,10,1)
#   cutoffs=c(1.96);
#   betap = c(-1,0.099);
#   symmetric<-1
#   tauhat<-1
#   numerical_integration<-1
#   X<-1:10
#   sigma<-1:10
#   lambdabar<-1
#   cdf<-VariationVarianceLogLikelihoodControls(lambdabar, tauhat, betap, cutoffs, symmetric, X, sigma,C,numerical_integration)
# }


