ReplicationAnalyticLogLikelihoodControls <-function(nuhat,tauhat, betap, cutoffs, symmetric,Z,sigmaZ2,C,numerical_integration) {
  

  #arguments: mean and stdev of distribution pi of mu
  #coefficient vector for step function p, in increasing order
  #vector of thresholds for step function p, in increasing order
  #dummy for symmetric step function: if one, all cutoffs should be positive
  #n by 2 matrix of estimates
  #vector of stdevs of X2
  
  n<-dim(Z)[1];
  k<-length(betap);
  betap<-as.matrix(betap,length(betap),1); # @TSA: betap in matrix umgewandelt. Auch in anderen Skripten nötig?
  #betap<-t(betap);
  Z1_dummies<-matrix(0,n,length(cutoffs)+1);
  
  if (all(sort(cutoffs)!=cutoffs)) {
    stop ("Unsorted cutoffs!")
  }
  
  if (symmetric ==1) {
    # check that cutoffs are positive
    if (!all(cutoffs>=0)) {
      stop("Needs positive cutoffs!")
    }
    Z1_dummies[,1]<-abs(Z[,1])<cutoffs[1];
    if (length(cutoffs)>1) {
      for (m in (2:length(cutoffs))) {
        Z1_dummies[,m]<-(abs(Z[,1])<cutoffs[m])*(abs(Z[,1])>=cutoffs[m-1]);
      }
    }
    Z1_dummies[,length(cutoffs)+1]<-abs(Z[,1])>=cutoffs[length(cutoffs)];
  } else {
    Z1_dummies[,1]<-Z[,1]<cutoffs[1];
    if (length(cutoffs)>1) {
      for (m in 2:length(cutoffs)) {
        Z1_dummies[,m]=(Z[,1]<cutoffs[m])*(Z[,1]>=cutoffs[m-1]);
      }
    }
    Z1_dummies[,length(cutoffs)+1]<-Z[,1]>=cutoffs[length(cutoffs)];
    
  }

  
########################
# Calculate objective function LLH
phat <- matrix (0, dim(Z1_dummies)[1], 1);
  
for (m in (1:dim(Z1_dummies)[2])) {
    Cmat <- matrix(rep(Z1_dummies[, m], dim(C)[2]), ncol=dim(C)[2], byrow = FALSE) * C  
    phat <- phat + Cmat %*% betap[m]
  }  


#vector of un-truncated likelihoods
if (symmetric==1){
  
  if (numerical_integration==1){
    #likelihoods calculated by numerical integration, gamma distribution
    g <- function(theta) {(0.5*dnorm(Z[,1]-theta)+0.5*dnorm(-Z[,1]-theta))*dnorm((Z[,2]-theta)/sigmaZ2)/sigmaZ2*dgamma(theta,nuhat,tauhat)}
    piZ1Z2 <- integrate(g,-Inf,Inf, subdivisions=2000)$value # numerically solving the function
    
  } else {
    
    #Monte-carlo integration
    set.seed(1)
    draw <- matrix(runif(10^5),1)
    theta_vec=dinvgamma(draw,nuhat,tauhat);
    theta_mat=matrix(rep(theta_vec,n),ncol=length(theta_vec),byrow = FALSE);
    Z1mat=matrix(rep(Z[,1],length(theta_vec)),ncol=length(theta_vec),byrow = FALSE);
    Z2mat=matrix(rep(Z[,2],length(theta_vec)),ncol=length(theta_vec),byrow = FALSE);
    sigmaZ2mat=matrix(rep(sigmaZ2,length(theta_vec)),ncol=length(theta_vec),byrow = FALSE);
    g <-  (0.5*dorm(Z1mat-theta_mat)+0.5*dnorm(-Z1mat-theta_mat))*dnorm((Z2mat-theta_mat)/sigmaZ2mat)/sigmaZ2mat
    piZ1Z2 <- mean(g,2)
  }
  
} else {
  piZ1Z2=rep(0,n); 
  for (i in 1:n) {
    Omega<-c(1+tauhat^2,tauhat^2,tauhat^2,sigmaZ2[i]^2+tauhat^2);
    Omega<-matrix(Omega,2,2);
    piZ1Z2[i]<-dmvnorm(c(Z[i,1], Z[i,2]),nuhat*c(1,1),Omega,log=FALSE)*0.5 +
      dmvnorm(c(Z[i,1], Z[i,2]),-nuhat*c(1,1),Omega,log=FALSE)*0.5
      }
  }



  
  
 # %normalizing constant
  prob_vec<-rep(0,length(cutoffs)+2);
  
  if (symmetric==1){
    for (m in (1:length(cutoffs))){
      
      if (numerical_integration==1){
        #Normalizing constant, gamma distribution
        g <- function(theta) {(dnorm(cutoffs[m]-theta)-dnorm(-cutoffs[m]-theta))*dgamma(theta,nuhat,tauhat)}
        prob_vec[,m+1] <- integrate(g,-Inf,Inf, subdivisions=2000)$value
      } else {
        #Monte Carlo Integration
        g <- (dnorm(cutoffs[m]-theta_vec)-dnorm(-cutoffs[m]-theta_vec))
        prob_vec[,m+1] <- mean(g)
      }
    }
prob_vec[length(cutoffs)+2]<-1
mean_Z1<-prob_vec[2:length(prob_vec)]-prob_vec[1:(length(prob_vec)-1)]

  } else {
    mu=nuhat
    sigma=sqrt((tauhat^2)+1)
    
    for (m in (1:length(cutoffs))){
      prob_vec[,m+1] <- 0.5*dnorm((cutoffs[m]-mu)/sigma)+0.5*dnorm((cutoffs[m]+mu)/sigma)
    }
    prob_vec[length(cutoffs)+2]<-1
    mean_Z1<-prob_vec[2:length(prob_vec)]-prob_vec[1:(length(prob_vec)-1)]
  }
  
  normalizingconst <- zeros(size(Z1_dummies,1),1)
  
  for (m in (1:dim(Z1_dummies)[2])){
      Cmat<-mean_Z1[m,1]*C
      normalizingconst=normalizingconst+Cmat*betap[m];
  }
  
  
  
  
  

  L=as.vector(phat)*as.vector(piZ1Z2)/normalizingconst; 

  if (normalizingconst<0) {
  }
 
  logL=log(L);
  #  %objective function; note the sign flip, since we are doing minimization

  LLH=-sum(log(L));
  
  return(list(LLH=LLH,logL=logL))

}

### Testing 
# if (test) {
#   
# 
# nuhat<-5
# tauhat<-1
# betap<-c(0.1,0.2,0.3,0.4,0.5)
# cutoffs<-c(-0.1503876, -0.2202119, 0.6473313, 0.9394543)
# symmetric<-0
# set.seed(1)
# Z<-matrix(rnorm(2*100),100,2)
# Z<-cbind(Test$V1,Test$V2)
# sigmaZ2<-rep(1,100)
# 
# 
# LLH<-ReplicationAnalyticLogLikelihood(nuhat, 
#   tauhat, 
#   betap, 
#   cutoffs, 
#   symmetric, 
#   Z, 
#   sigmaZ2)
# 
# }