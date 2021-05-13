
replication_analytic_llh <-function(nuhat,tauhat, betap, cutoffs, symmetric,Z,sigmaZ2,C) {
  # nuhat <- Psihat0[1]
  # tauhat <- Psihat0[2]
  # betap <- c(reshape(t(Psihat0[-c(1,2)]), c(length(Psihat0[-c(1,2)]) / length(cutoffs), length(cutoffs))), 1)

  n <- nrow(Z)
  k <- length(betap)
  betap <- t(t(betap))
  Z1_dummies <- zeros(n,length(cutoffs)+1)

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
phat <- zeros(nrow(Z1_dummies), 1)

for (m in (1:ncol(Z1_dummies))) {
    Cmat <- repmat(Z1_dummies[ ,m],1,ncol(C)) * C
    phat <- phat + Cmat * betap[m]
  }


#vector of un-truncated likelihoods
if (symmetric==1){

    #Monte-carlo integration
    set.seed(1)
    draw <- matrix(runif(10^5),1)
    theta_vec <- qgamma(draw, nuhat,scale = tauhat)
    theta_mat <- matrix(theta_vec,n,length(theta_vec))
    Z1mat <- repmat(Z[,1],1,length(theta_vec))
    Z2mat <- repmat(Z[,2],1,length(theta_vec))
    sigmaZ2mat <- repmat(sigmaZ2,1,length(theta_vec))
    g <-  (0.5*dnorm(Z1mat-theta_mat)+0.5*dnorm(-Z1mat-theta_mat))*dnorm((Z2mat-theta_mat)/sigmaZ2mat)/sigmaZ2mat
    piZ1Z2 <- rowMeans(g,2)

  } else {

  piZ1Z2=rep(0,n)
  for (i in 1:n) {
    Omega<-c(1+tauhat^2,tauhat^2,tauhat^2,sigmaZ2[i]^2+tauhat^2);
    Omega<-matrix(Omega,2,2);
    piZ1Z2[i]<-dmvnorm(c(Z[i,1], Z[i,2]),nuhat*c(1,1),Omega,log=FALSE)*0.5 +
      dmvnorm(c(Z[i,1], Z[i,2]),-nuhat*c(1,1),Omega,log=FALSE)*0.5
      }
  }





 # normalizing constant
prob_vec<-rep(0,length(cutoffs)+2);

  if (symmetric==1){
    for (m in (1:length(cutoffs))){

        #Monte Carlo Integration
        g <- (pnorm(cutoffs[m]-theta_vec)-pnorm(-cutoffs[m]-theta_vec))
        prob_vec[m+1] <- mean(g)

        prob_vec[length(cutoffs)+2]<-1
        mean_Z1<-prob_vec[2:length(prob_vec)]-prob_vec[1:(length(prob_vec)-1)]
    }
  } else {
    mu=nuhat
    sigma=sqrt((tauhat^2)+1)

    for (m in (1:length(cutoffs))){
      prob_vec[m+1] <- 0.5*pnorm((cutoffs[m]-mu)/sigma)+0.5*pnorm((cutoffs[m]+mu)/sigma)
    }
    prob_vec[length(cutoffs)+2]<-1
    mean_Z1<-prob_vec[2:length(prob_vec)]-prob_vec[1:(length(prob_vec)-1)]
  }

  normalizingconst <- zeros(nrow(Z1_dummies),1)

  for (m in (1:ncol(Z1_dummies))){
      Cmat<-mean_Z1[m]*C
      normalizingconst=normalizingconst+Cmat*betap[m];
  }

  L <- as.vector(phat)*as.vector(piZ1Z2)/normalizingconst;

  logL=log(L)

  LLH=-sum(log(L))

  return(list("LLH"= LLH, "logL" = logL))

}


