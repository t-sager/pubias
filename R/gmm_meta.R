# This Code follows closley the Matlab Code 2019. But it only uses the GMM/Meta Approach.
# I can do a similar thing with the other 3 dimensions --> Translate the Matlab Code to R Code
# for each dimension!

Step_function_normal_cdf <- function(X,theta,sigma,betap,cutoffs,symmetric) {
  #Arguments:
  #X: point at which to evaluate cdf
  #theta: parameter value under which to evaluate cdf
  #sigma:standarad deviation of (untruncated) normal variable
  #sigma: stdev of distribution pi of mu
  #cutoffs: vector of thresholds for step function p, in increasing order
  #cutoffs are given in terms of X, not z statistics
  #Symmetric: dummy indicating whether publication probability is symmetric around zero.
  #In symmetric case, cutoffs should include only positive values
  #Note: publication probability for largest category (i.e. for those points beyond largest cutoff) normalized to one.

  if (length(betap)!=(length(cutoffs)+1)) {
    stop("length of betap must be one greater than length of cutoffs ")
  }

  #For symmetric case, create symmetrized version of cutoffs and coefficients
  if (symmetric==1) {
    cutoffs_u=matrix(0,2*length(cutoffs),1)
    betap_u=matrix(0,1,2*length(cutoffs))
    for (n in c(1:length(cutoffs))) {
      cutoffs_u[n]=-cutoffs[length(cutoffs)+1-n]
      betap_u[n]=betap[length(betap)+1-n]
    }
    for (n in c(1:length(cutoffs))) {
      cutoffs_u[length(cutoffs)+n]=cutoffs[n]
      betap_u[length(cutoffs)+n]=betap[n]
    }
    betap_u[length(betap_u)+1]=1
    betap_u=as.matrix(betap_u)
  } else {
    cutoffs_u=cutoffs
    betap_u=betap
  }
  betap_u=as.matrix(betap_u)

  #Calculate denominatore in cdf
  prob_vec=matrix(0,length(cutoffs_u)+2,1)
  for (m in c(1:length(cutoffs_u))) {
    prob_vec[m+1]=pnorm((cutoffs_u[m]-theta)/sigma)
  }
  prob_vec[length(prob_vec)]=1
  mean_Z1=as.matrix(prob_vec[2:length(prob_vec),1])-as.matrix(prob_vec[1:(length(prob_vec)-1),1])
  denominator=t(mean_Z1)%*%betap_u

  #Calculate numerator in cdf
  cutoffs_u[length(cutoffs_u)+1]=Inf
  if (X<=cutoffs_u[1]) {
    numerator=pnorm((X-theta)/sigma)*betap_u[1]
  } else {
    numerator=pnorm((cutoffs_u[1]-theta)/sigma)*betap_u[1]
    m=1
    while (X>cutoffs_u[m]) {
      Xcap=min(X,cutoffs_u[m+1])
      numerator=numerator+(pnorm((Xcap-theta)/sigma)-pnorm((cutoffs_u[m]-theta)/sigma))*betap_u[m+1]
      m=m+1
    }
  }
  #evaluate cdf
  cdf=numerator/denominator

  return(cdf)
}

MetastudyMoments <- function(betap,cutoffs,symmetric,X,sigma) {
  #Calculate GMM Moments for metastudy applications
  n=length(X);
  betap=as.matrix(betap);
  #betap=t(betap)
  betap=rbind(betap,1)

  #regressors for step function p
  T=X/sigma;

  Tpowers=matrix(0,n,length(cutoffs)+1)
  if (symmetric==1) {
    Tpowers[,1]=ifelse(abs(T)<cutoffs[1],1,0)
    if (length(cutoffs)>1) {
      for (m in c(2:length(cutoffs))) {
        Tpowers[,m]=(ifelse(abs(T)<cutoffs[m],1,0))*(ifelse(abs(T)>=cutoffs[m-1]));
      }
    }
    Tpowers[,length(cutoffs)+1]=ifelse(abs(T)>=cutoffs[length(cutoffs)],1,0);
  } else {
    Tpowers[,1]=ifelse(T<cutoffs[1],1,0);
    if (length(cutoffs)>1) {
      for (m in c(2:length(cutoffs))) {
        Tpowers[,m]=(ifelse(T<cutoffs[m],1,0))*(ifelse(T>=cutoffs[m-1],1,0));
      }
    }
    Tpowers[,length(cutoffs)+1]=ifelse(T>=cutoffs[length(cutoffs)],1,0)
  }

  phat=Tpowers%*%betap;

  Xmat1=matrix(X,length(X),length(X));
  Xmat2=t(Xmat1);
  sigmamat1=matrix(sigma,length(X),length(X));
  sigmamat2=t(sigmamat1);
  indicator=ifelse(sigmamat1>=sigmamat2,1,0);
  sigmadiff=sqrt(indicator*(sigmamat1^2-sigmamat2^2))+sqrt((1-indicator)*(sigmamat2^2-sigmamat1^2))
  sigmadiff=Re(sigmadiff)
  pmat1=matrix(phat,length(X),length(X))
  pmat2=t(pmat1)

  moment_mean=matrix(0,1,length(cutoffs))
  rhat=matrix(0,n,length(cutoffs))

  if (symmetric==1) {
    for (k in c(1:length(cutoffs))) {
      c=cutoffs[k]
      base_moments=(pmat2^-1)*(pmat1^-1)*indicator*((pnorm((c*sigmamat1-Xmat2)/sigmadiff)-pnorm((-c*sigmamat1-Xmat2)/sigmadiff))-(((ifelse(Xmat1<=c*sigmamat1,1,0))-(ifelse(Xmat1<=-c*sigmamat1,1,0)))))+(pmat2^-1)*(pmat1^-1)*(1-indicator)*((pnorm((c*sigmamat2-Xmat1)/sigmadiff)-pnorm((-c*sigmamat2-Xmat1)/sigmadiff))-(((ifelse(Xmat2<=c*sigmamat2,1,0))-(ifelse(Xmat2<=-c*sigmamat2,1,0)))))

      base_moments=base_moments-diag(diag(base_moments))
      base_moments[is.nan(base_moments)]=0
      moment_mean[1,k]=((choose(n,2)^-1)/2)*sum(base_moments)
      #Built normalization to variance into rhat to ensure correct answers with culstered variance estimators
      rhat[,k] = 2*((n-1)^-1)*as.matrix(rowSums(base_moments,na.rm=FALSE))
    }
  } else {
    for (k in c(1:length(cutoffs))) {
      c=cutoffs[k]
      base_moments = (pmat2^-1)*(pmat1^-1)*indicator*(pnorm((c*sigmamat1-Xmat2)/sigmadiff)-(ifelse(Xmat1<=c*sigmamat1,1,0)))+(pmat2^-1)*(pmat1^-1)*(1-indicator)*(pnorm((c*sigmamat2-Xmat1)/sigmadiff)-(ifelse(Xmat2<=c*sigmamat2,1,0)))
      base_moments=base_moments-diag(diag(base_moments))
      base_moments[is.nan(base_moments)]=0
      moment_mean[1,k]=((choose(n,2)^-1)/2)*sum(base_moments)
      #Built normalization to variance into rhat to ensure correct answers with culstered variance estimators
      rhat[,k] = 2*((n-1)^-1)*as.matrix(rowSums(base_moments,na.rm=FALSE))
    }
  }
  moment_var=cov(rhat)
  raw_moments=rhat
  return(list("moment_mean"= moment_mean, "moment_var" = moment_var, "raw_moments" = raw_moments))
}

ComputingMetastudyMoments <- function(betap,cutoffs,symmetric,X,sigma) {
  #Calculate GMM Moments for metastudy applications
  n=length(X);
  gamma_1=as.matrix(betap[1])
  gamma_2=as.matrix(betap[2])
  betap=as.matrix(betap[3:length(betap)]);

  #regressors for step function p
  T=X/sigma;

  Tpowers=matrix(0,n,length(cutoffs)+1)
  if (symmetric==1) {
    Tpowers[,1]=ifelse(abs(T)<cutoffs[1],1,0)
    if (length(cutoffs)>1) {
      for (m in c(2:length(cutoffs))) {
        Tpowers[,m]=(ifelse(abs(T)<cutoffs[m],1,0))*(ifelse(abs(T)>=cutoffs[m-1]));
      }
    }
    Tpowers[,length(cutoffs)+1]=ifelse(abs(T)>=cutoffs[length(cutoffs)],1,0);
  } else {
    Tpowers[,1]=ifelse(T<cutoffs[1],1,0);
    if (length(cutoffs)>1) {
      for (m in c(2:length(cutoffs))) {
        Tpowers[,m]=(ifelse(T<cutoffs[m],1,0))*(ifelse(T>=cutoffs[m-1],1,0));
      }
    }
    Tpowers[,length(cutoffs)+1]=ifelse(T>=cutoffs[length(cutoffs)],1,0)
  }

  phat=Tpowers%*%betap;

  Xmat1=matrix(X,length(X),length(X));
  Xmat2=t(Xmat1);
  sigmamat1=matrix(sigma,length(X),length(X));
  sigmamat2=t(sigmamat1);
  indicator=ifelse(sigmamat1>=sigmamat2,1,0);
  sigmadiff=sqrt(indicator*(sigmamat1^2-sigmamat2^2))+sqrt((1-indicator)*(sigmamat2^2-sigmamat1^2))
  sigmadiff=Re(sigmadiff)
  pmat1=matrix(phat,length(X),length(X))
  pmat2=t(pmat1)

  moment_mean=matrix(0,1,length(cutoffs))
  rhat=matrix(0,n,length(cutoffs))

  if (symmetric==1) {
    for (k in c(1:length(cutoffs))) {
      c=cutoffs[k]
      base_moments=(pmat2^-1)*(pmat1^-1)*indicator*((pnorm((c*sigmamat1-Xmat2)/sigmadiff)-pnorm((-c*sigmamat1-Xmat2)/sigmadiff))-(((ifelse(Xmat1<=c*sigmamat1,1,0))-(ifelse(Xmat1<=-c*sigmamat1,1,0)))))+(pmat2^-1)*(pmat1^-1)*(1-indicator)*((pnorm((c*sigmamat2-Xmat1)/sigmadiff)-pnorm((-c*sigmamat2-Xmat1)/sigmadiff))-(((ifelse(Xmat2<=c*sigmamat2,1,0))-(ifelse(Xmat2<=-c*sigmamat2,1,0)))))

      base_moments=base_moments-diag(diag(base_moments))
      base_moments[is.nan(base_moments)]=0
      moment_mean[1,k]=((choose(n,2)^-1)/2)*sum(base_moments)
      #Built normalization to variance into rhat to ensure correct answers with culstered variance estimators
      rhat[,k] = 2*((n-1)^-1)*as.matrix(rowSums(base_moments,na.rm=FALSE))
    }
  } else {
    for (k in c(1:length(cutoffs))) {
      c=cutoffs[k]
      base_moments = (pmat2^-1)*(pmat1^-1)*indicator*(pnorm((c*sigmamat1-Xmat2)/sigmadiff)-(ifelse(Xmat1<=c*sigmamat1,1,0)))+(pmat2^-1)*(pmat1^-1)*(1-indicator)*(pnorm((c*sigmamat2-Xmat1)/sigmadiff)-(ifelse(Xmat2<=c*sigmamat2,1,0)))
      base_moments=base_moments-diag(diag(base_moments))
      base_moments[is.nan(base_moments)]=0
      moment_mean[1,k]=((choose(n,2)^-1)/2)*sum(base_moments)
      #Built normalization to variance into rhat to ensure correct answers with culstered variance estimators
      rhat[,k] = 2*((n-1)^-1)*as.matrix(rowSums(base_moments,na.rm=FALSE))
    }
  }
  base_mean_moment=(phat^-1)*(X-as.vector(gamma_1))
  base_sd_moment=(phat^-1)*((X-as.vector(gamma_1))^2-sigma^2-as.vector(gamma_2)^2)
  mean_moment=(1/n)*sum(base_mean_moment,na.rm=FALSE)
  sd_moment=(1/n)*sum(base_sd_moment,na.rm=FALSE)
  moment_mean=cbind(mean_moment, sd_moment, moment_mean)
  rhat=cbind(base_mean_moment,base_sd_moment,rhat)
  moment_var=cov(rhat)
  raw_moments=rhat
  return(list("moment_mean"= moment_mean, "moment_var" = moment_var, "raw_moments" = raw_moments))
}

MetastudyGMMObjective <- function(betap,cutoffs,symmetric,X,sigma,cluster_ID) {
  #Calculate continuously updating GMM objective
  mom=MetastudyMoments(betap,cutoffs,symmetric,X,sigma);
  moments_mean=mom$moment_mean;
  rhat=mom$raw_moments;
  Sigma_hat=Clustered_covariance_estimate(rhat,cluster_ID);
  objective=length(X)*moments_mean%*%solve(Sigma_hat)%*%t(moments_mean);
  return(objective)
}

ComputingMetastudyGMMObjective <- function(beta_theta,betap,cutoffs,symmetric,X,sigma,cluster_ID) {
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

DescriptiveStats <- function(X, sigma) {
  #X has 1 column, initial estimate unscaled, sigma=sigma1

  critval=1.96;
  n=length(X);
  circlesize=75/sqrt(n);
  color<-brewer.pal(8,"Greys");
  GREY=color[3];
  SLIDE_GREY=color[5];
  color<-brewer.pal(8,"Purples")
  BLUE=color[6];

  circlesize=circlesize*4;

  significant=(ifelse(abs(X/sigma)>critval,1,0));
  nooutlier=(ifelse(sigma<50,1,0));

  #Scatter plot for significant vs. insignificant values
  rangeX=1.1*max(max(abs(X)),max(abs(sigma[nooutlier]))*critval);

  pdf(here("FiguresandTables","GMMresultsScatter.pdf"), width = 3, height = 3)
  par(mar = c(3, 3, 1, 1))
  plot(x = X[!significant&nooutlier], y = sigma[!significant&nooutlier],
       pch = 19, frame = FALSE,
       xlab = '' , ylab = '', xlim = c(-rangeX,rangeX), ylim = c(0,rangeX/critval), col = SLIDE_GREY,
       cex.axis=0.75);
  par(new=TRUE)
  lines(c(0,rangeX),c(0,rangeX/critval),col=SLIDE_GREY);
  par(new=TRUE)
  lines(c(-rangeX,0),c(rangeX/critval,0),col=SLIDE_GREY);
  par(new=TRUE)
  plot(x = X[!significant&nooutlier], y = sigma[!significant&nooutlier],
       pch = 19, frame = FALSE,
       xlab = '' , ylab = '', xlim = c(-rangeX,rangeX), ylim = c(0,rangeX/critval), col = SLIDE_GREY,
       cex.axis=0.75,add=TRUE);
  title(ylab=TeX('$\\Sigma$'), xlab=TeX('$X$'), line=1.75, cex.lab=0.8)
  par(new=TRUE)
  plot(x = X[significant&nooutlier], y = sigma[significant&nooutlier],
       pch = 19, frame = FALSE,
       xlab = '', ylab = '', xlim = c(-rangeX,rangeX), ylim = c(0,rangeX/critval), col = BLUE,
       cex.axis=0.75,add=TRUE);
  dev.off()

  #Histogram of Z statistics
  Zuse=X/sigma
  if (max(abs(Zuse))>10) {
    Zuse=Zuse[abs(Zuse)<=6]
  }

  ll=floor(min(Zuse))
  lleven=floor(ll/2)*2
  uu=ceiling(max(Zuse))

  if (n>=30) {
    uu2=ceiling((uu-0.36)/0.32)*0.32+0.36;
    ll2=floor((ll+0.36)/0.32)*0.32-0.36;
    edges=as.numeric(c(seq(ll2,-0.36,0.32),0,seq(0.36,uu2,0.32)))
  } else {
    uu2=ceiling((uu-0.68)/0.64)*0.64+0.68;
    ll2=floor((ll+0.68)/0.64)*0.64-0.68;
    edges=as.numeric(c(seq(ll2,-0.68,0.64),0,seq(0.68,uu2,0.64)))
  }

  pdf(here("FiguresandTables", "GMMresultsScatterHist.pdf"), width = 6, height = 2.5)
  h<-hist(Zuse,
          plot=FALSE,
          right=FALSE,
          breaks=edges
  ) #right=FALSE means the bins are right open left closed (,]

  h$counts=h$counts/sum(h$counts)

  par(mfrow = c(1, 2),mar = c(3, 3, 1, 1))
  plot(h,col=BLUE,xlab='',
       ylab='',
       main="",
       xlim=c(min(edges),max(edges)),
       ylim=c(0,max(h$counts)*1.1),freq=TRUE,cex.axis=0.75)
  title(ylab="Density", xlab=TeX('$X/\ \\Sigma$'), line=1.75, cex.lab=0.8)
  par(new=TRUE)
  lines(c(-1.96,-1.96),c(0,max(Zuse)*1.1),col=SLIDE_GREY,lty=1, lwd=1)
  par(new=TRUE)
  lines(c(1.96,1.96),c(0,max(Zuse)*1.1),col=SLIDE_GREY,lty=1, lwd=1)
  par(new=TRUE)
  lines(c(0,0),c(0,max(Zuse)*1.1),col=SLIDE_GREY,lty=2, lwd=1)

  plot(x = X[!significant&nooutlier], y = sigma[!significant&nooutlier],
       pch = 19, frame = FALSE,
       xlab = '', ylab = '', xlim = c(-rangeX,rangeX), ylim = c(0,rangeX/critval), col = SLIDE_GREY,cex.axis=0.75);
  par(new=TRUE)
  lines(c(0,rangeX),c(0,rangeX/critval),col=SLIDE_GREY);
  par(new=TRUE)
  lines(c(-rangeX,0),c(rangeX/critval,0),col=SLIDE_GREY);
  par(new=TRUE)
  plot(x = X[!significant&nooutlier], y = sigma[!significant&nooutlier],
       pch = 19, frame = FALSE,
       xlab = '', ylab = '', xlim = c(-rangeX,rangeX), ylim = c(0,rangeX/critval), col = SLIDE_GREY,cex.axis=0.75,add=TRUE);
  par(new=TRUE)
  plot(x = X[significant&nooutlier], y = sigma[significant&nooutlier],
       pch = 19, frame = FALSE,
       xlab = '', ylab = '', xlim = c(-rangeX,rangeX), ylim = c(0,rangeX/critval), col = BLUE,cex.axis=0.75,add=TRUE);
  title(ylab=TeX('$\\Sigma$'), xlab=TeX('$X$'), line=1.75, cex.lab=0.8)

  dev.off()
}

Clustered_covariance_estimate <- function(g,cluster_index) {
  sorted=sort(cluster_index,decreasing=FALSE,index.return=TRUE)
  cluster_index=as.matrix(sorted$x)
  I=as.matrix(sorted$ix)
  g=as.matrix(g[I,])
  g=g-matrix(rep(colMeans(g),dim(g)[1]),ncol=dim(g)[2],byrow=TRUE)
  gsum=matrix(cumsum(g),dim(g)[1],dim(g)[2])
  index_diff=as.matrix(ifelse(cluster_index[2:length(cluster_index),1]!=cluster_index[1:(length(cluster_index)-1),1],1,0))
  index_diff=rbind(index_diff,1)
  gsum=subset(gsum,index_diff==1)
  gsum=rbind(gsum[1,],diff(gsum))
  Sigma=(1/(dim(g)[1]-1))*(t(gsum)%*%gsum)
  return(Sigma)
}

EstimatingSelection <- function(X,sigma,symmetric,cluster_ID,cutoffs,Studynames) {
  GMM_obj <- function(Psi){
    MetastudyGMMObjective(Psi, cutoffs, symmetric, X, sigma, cluster_ID)+max(-min(Psi),0)*10^5
  }

  mini<-optim(par=Psihat0,fn=GMM_obj,method="BFGS",control = list(abstol=10^-8,maxit=10^5));

  #For more accurate optimization, use the following:
  #Psihat1=mini$par;
  #Objval=mini$value;

  #mini<-optim(par=Psihat1,fn=GMM_obj);

  Psihat=mini$par;
  Objval=mini$value;

  GMM_obj_new<-function(Psi){
    #all_in_one=cbind(Psi,Psihat)
    ComputingMetastudyGMMObjective(Psi,Psihat,cutoffs,symmetric,X,sigma,cluster_ID)
  }

  mini<-optim(par=Psihat0_theta,fn=GMM_obj_new,method="BFGS",control = list(abstol=10^-8,maxit=10^5));

  #For more accurate optimization, use the following:
  #Psihat1_theta=mini$par;
  #Objval_theta=mini$value;

  #mini<-optim(par=Psihat1_theta,fn=GMM_obj_new);

  Psihat_theta=mini$par;
  Objval_theta=mini$value;

  Psihat_theta[2]=abs(Psihat_theta[2])

  Psihat_all = cbind(Psihat_theta,Psihat);

  mom=ComputingMetastudyMoments(cbind(Psihat_all,1),cutoffs,symmetric,X,sigma)
  moments=mom$moment_mean;
  Sigma_temp=mom$moment_var;
  rhat=mom$raw_moments;

  Sigma_hat=Clustered_covariance_estimate(rhat,cluster_ID);

  stepsize=10^-3;
  G=matrix(0,dim(moments)[2],length(Psihat_all));
  for (n1 in c(1:length(Psihat_all))) {
    beta_plus=Psihat_all;
    beta_plus[n1]=beta_plus[n1]+stepsize
    mom=ComputingMetastudyMoments(c(beta_plus,1),cutoffs,symmetric,X,sigma);
    moments_plus=mom$moment_mean;
    beta_minus=Psihat_all;
    beta_minus[n1]=beta_minus[n1]-stepsize
    mom=ComputingMetastudyMoments(cbind(beta_minus,1),cutoffs,symmetric,X,sigma);
    moments_minus=mom$moment_mean;
    G[,n1]=t(moments_plus-moments_minus)/(2*stepsize);
  }

  #In any case you get an error on an almost numerically singular matrix,
  #try using ``MASS::ginv'' instead of ``solve''. This command will run a generalized
  #inverse matrix and usually deals with this error.

  Varhat_all=solve(t(G)%*%solve(Sigma_hat)%*%G)/length(X);
  se_robust=sqrt(diag(Varhat_all));
  dof=dim(moments)[2];
  Varhat=Varhat_all[3:dim(Varhat_all)[1],3:dim(Varhat_all)[2]];

  if (length(cutoffs)==1) {
    Psi_grid=c(seq(0.001,5,0.001),10,10^5)
    S_store=matrix(0,length(Psi_grid),1)
    for (m in c(1:length(Psi_grid))) {
      S_store[m,1]=GMM_obj(Psi_grid[m])
    }

    CS_LB=min(Psi_grid[S_store<qchisq(0.95,1)]);
    CS_UB=max(Psi_grid[S_store<qchisq(0.95,1)]);
  } else if (length(cutoffs)==2) {
    gridsteps=200;
    Psi_grid1=seq(0.01,0.25,(0.25-0.01)/(gridsteps-1));
    Psi_grid2=seq(0.05,5,(5-0.05)/(gridsteps-1));
    S_store=matrix(0,length(Psi_grid1),length(Psi_grid2));
    for (m1 in c(1:length(Psi_grid1))) {
      for (m2 in c(1:length(Psi_grid2))) {
        Psi_temp=c(Psi_grid1[m1],Psi_grid2[m2]);
        S_store[m1,m2]=GMM_obj(Psi_temp);
      }

    }
    Psi_grid=cbind(t(Psi_grid1),t(Psi_grid2));
  } else {
    Psi_grid=c(seq(0.001,5,0.001),10,10^5)
    S_store=matrix(0,length(Psi_grid),1)
  }

  Psi_grid=as.matrix(Psi_grid)

  # Ab hier SelectionTableGMM
  sink(here("FiguresandTables","GMMEstimatesAllSelectionModel.tex"))
  columns_names<-c("$\\Theta$","$\\Sigma$")
  if (length(Psihat)==1) {
    columns_names<-c(columns_names,"$\\beta_{p}$")
  } else {
    for (i in c(1:length(Psihat))) {
      columns_names<-c(columns_names,paste0("$\\beta_{p,",i,"}$"))
    }

  }
  results<-matrix(rbind(format(Psihat_all,digits=2),format(se_robust,digits=2)),nrow=2,dimnames = list(c(1:2),columns_names))
  results[2,]<-paste0("(",results[2,],")")
  results<-as.matrix(results)
  alignment <- c("c")
  for (i in c(1:(length(Psihat_all)))) {
    alignment<-c(alignment,"c")
  }
  print(xtable(results,align = alignment), sanitize.rownames.function = function(a) {''}, sanitize.colnames.function = function(a) {a})
  sink()

  sink(here("FiguresandTables","GMMresultsSelectionModel.tex"))
  if (length(Psihat)>1) {
    columns_names<-c("$\\beta_{p,1}$")
    for (i in c(2:length(Psihat))) {
      columns_names<-c(columns_names,paste0("$\\beta_{p,",i,"}$"))
    }
  } else {
    columns_names<-c("$\\beta_{p}$")
  }

  results<-matrix(rbind(format(Psihat,digits=2),format(se_robust[3:length(se_robust)],digits=2)),nrow=2,dimnames = list(c(1:2),columns_names))
  results[2,]<-paste0("(",results[2,],")")
  results<-as.matrix(results)
  alignment <- c("c")
  if (length(Psihat)>=1) {
    for (i in c(1:(length(Psihat)))) {
      alignment<-c(alignment,"c")
    }
  }

  print(xtable(results,align = alignment), sanitize.rownames.function = function(a) {''}, sanitize.colnames.function = function(a) {a})
  sink()

  dof=length(moments[3:length(moments)])
  len=length(Psihat)
  if (dorobust==1) {
    if (len==1) {
      CS_LB = min(Psi_grid[S_store<qchisq(0.95,dof)])
      CS_UB = max(Psi_grid[S_store<qchisq(0.95,dof)])

      if (CS_LB==min(Psi_grid)) {
        CS_LB=0
      }
      if (CS_UB==max(Psi_grid)) {
        CS_UB=c("$\\infty$")
      }
      sink(here("FiguresandTables","GMMresultsSelectionModelCS.tex"))
      columns_CS_names<-c(paste0("$\\beta_p$", " Lower Bound"),paste0("$\\beta_p$"," Upper Bound"))
      CS_bounds<-matrix(c(CS_LB,CS_UB),nrow=1,dimnames = list(c(1),columns_CS_names))
      alignment <- c("c","c","c")
      print(xtable(CS_bounds,align = alignment), sanitize.text.function = function(a) {a}, sanitize.rownames.function = function(a) {''}, sanitize.colnames.function = function(a) {a})
      sink()
    } else if (len==2) {
      color<-brewer.pal(8,"Greys");
      GREY=color[3];
      S_CS = ifelse(S_store<=qchisq(0.95,dof),1,0)
      S_CS=as.matrix(S_CS)
      S_CS=t(S_CS)
      pdf(here("FiguresandTables","GMMresultsCS.pdf"), width = 4, height = 4)
      contour(Psi_grid[,1],Psi_grid[,2],S_CS,xlab="$\\beta_{1}$",ylab="$\\beta_{2}$");
      dev.off()
    }

  }
  return(list("Psihat"= Psihat, "Varhat" = Varhat))
}

HorizontalBars <- function(X,sigma,Psihat,Varhat,cutoffs,symmetric,Studynames) {
  #library(RColorBrewer)
  #install.packages('latex2exp')
  #library(latex2exp)

  color<-brewer.pal(8,"Greys")
  GREY=color[3]
  SLIDE_GREY=color[5]
  color<-brewer.pal(8,"Purples")
  BLUE=color[6]

  #normalize estimates from initial study
  Z1=sort(X/sigma)
  Z1=as.matrix(Z1)

  alpha=0.05
  bonf_alpha=0.04
  bonf_beta=0.01

  #Pad Psihat with additional zeros in symmetric specifications to make dimensions conform

  Psihat_use=Psihat

  #Calculate corrected estimates and confidence bounds

  Z1_U=matrix(0,length(Z1),1)
  Z1_L=matrix(0,length(Z1),1)
  Z1_M=matrix(0,length(Z1),1)

  stepsize=10^-3
  for (n in c(1:length(Z1))) {
    g_U <- function(lambda) {
      (alpha/2-Step_function_normal_cdf(Z1[n],lambda,1,cbind(Psihat_use,1),cutoffs,symmetric))^2
    }
    g_L <- function(lambda) {
      (1-alpha/2-Step_function_normal_cdf(Z1[n],lambda,1,cbind(Psihat_use,1),cutoffs,symmetric))^2
    }
    g_M <- function(lambda) {
      (1/2-Step_function_normal_cdf(Z1[n],lambda,1,cbind(Psihat_use,1),cutoffs,symmetric))^2
    }

    min_U<-optim(par=Z1[n],fn=g_U)
    min_U<-optim(par=min_U$par,fn=g_U,method="BFGS")
    min_L<-optim(par=Z1[n],fn=g_L)
    min_L<-optim(par=min_L$par,fn=g_L,method="BFGS")
    min_M<-optim(par=Z1[n],fn=g_M)
    min_M<-optim(par=min_M$par,fn=g_M,method="BFGS")

    Z1_U[n,1]=min_U$par
    Z1_L[n,1]=min_L$par
    Z1_M[n,1]=min_M$par
  }

  #Calculate bonferroni corrected estimates and confidence bounds
  Z1_UB=matrix(0,length(Z1),1)
  Z1_LB=matrix(0,length(Z1),1)
  sigma_U=matrix(0,length(Z1),1)
  sigma_L=matrix(0,length(Z1),1)

  for (n in c(1:length(Z1))) {
    g_UB <- function(lambda) {
      (bonf_alpha/2-Step_function_normal_cdf(Z1[n],lambda,1,cbind(Psihat_use,1),cutoffs,symmetric))^2
    }
    g_LB <- function(lambda) {
      (1-bonf_alpha/2-Step_function_normal_cdf(Z1[n],lambda,1,cbind(Psihat_use,1),cutoffs,symmetric))^2
    }

    min_UB<-optim(par=Z1_U[n,1],fn=g_UB)
    min_UB<-optim(par=min_UB$par,fn=g_UB,method="BFGS")
    min_LB<-optim(par=Z1_L[n,1],fn=g_LB)
    min_LB<-optim(par=min_LB$par,fn=g_LB,method="BFGS")

    Z1_UB[n,1]=min_UB$par
    Z1_LB[n,1]=min_LB$par

    #Calculate derivatives of corrections with respect to parameters via implicit function theorem
    thetaU_plus=Z1_UB[n,1]+stepsize
    thetaU_minus=Z1_UB[n,1]-stepsize
    F_Uplus=Step_function_normal_cdf(Z1[n],thetaU_plus,1,cbind(Psihat_use,1),cutoffs,symmetric)
    F_Uminus=Step_function_normal_cdf(Z1[n],thetaU_minus,1,cbind(Psihat_use,1),cutoffs,symmetric)
    dFUdtheta=(F_Uplus-F_Uminus)/(2*stepsize)

    thetaL_plus=Z1_LB[n,1]+stepsize
    thetaL_minus=Z1_LB[n,1]-stepsize
    F_Lplus=Step_function_normal_cdf(Z1[n],thetaL_plus,1,cbind(Psihat_use,1),cutoffs,symmetric)
    F_Lminus=Step_function_normal_cdf(Z1[n],thetaL_minus,1,cbind(Psihat_use,1),cutoffs,symmetric)
    dFLdtheta=(F_Lplus-F_Lminus)/(2*stepsize)

    dFUdbeta=matrix(0,length(Psihat_use),1)
    dFLdbeta=matrix(0,length(Psihat_use),1)
    for (n1 in c(1:length(Psihat_use))) {
      Psi_plus=Psihat_use
      Psi_plus[n1]=Psi_plus[n1]+stepsize
      Psi_minus=Psihat_use
      Psi_minus[n1]=Psi_minus[n1]-stepsize

      F_Uplus=Step_function_normal_cdf(Z1[n],Z1_UB[n,1],1,cbind(Psi_plus,1),cutoffs,symmetric)
      F_Uminus=Step_function_normal_cdf(Z1[n],Z1_UB[n,1],1,cbind(Psi_minus,1),cutoffs,symmetric)
      dFUdbeta[n1,1]=(F_Uplus-F_Uminus)/(2*stepsize)

      F_Lplus=Step_function_normal_cdf(Z1[n],Z1_LB[n,1],1,cbind(Psi_plus,1),cutoffs,symmetric)
      F_Lminus=Step_function_normal_cdf(Z1[n],Z1_LB[n,1],1,cbind(Psi_minus,1),cutoffs,symmetric)
      dFLdbeta[n1,1]=(F_Lplus-F_Lminus)/(2*stepsize)
    }

    dFUdtheta=as.numeric(dFUdtheta)
    dFLdtheta=as.numeric(dFLdtheta)

    dthetaUdbeta=-dFUdbeta/dFUdtheta
    dthetaLdbeta=-dFLdbeta/dFLdtheta

    if (symmetric==1) {
      sigma_thetaU=(t(dthetaUdbeta)%*%Varhat%*%dthetaUdbeta)^0.5
      sigma_thetaL=(t(dthetaLdbeta)%*%Varhat%*%dthetaLdbeta)^0.5
    } else {
      sigma_thetaU=(t(dthetaUdbeta)%*%Varhat%*%dthetaUdbeta)^0.5
      sigma_thetaL=(t(dthetaLdbeta)%*%Varhat%*%dthetaLdbeta)^0.5
    }
  }

  sigma_U[n,1]=sigma_thetaU
  sigma_L[n,1]=sigma_thetaL

  Z1_UB_store=Z1_UB
  Z1_LB_store=Z1_LB
  Z1_UB=Z1_UB+qnorm(1-bonf_beta)*sigma_U
  Z1_LB=Z1_LB-qnorm(1-bonf_beta)*sigma_L

  count=1
  store=matrix(0,length(seq(-10,10,0.1)),1)
  for (lambda in seq(-10,10,0.1)) {
    store[count]=g_M(lambda)
    count=count+1
  }


  # Plot original and adjusted confidence sets ------------------------------

  Rl=min(min(Z1_L),min(Z1)-2)-0.5
  Ru=max(max(Z1_U),max(Z1)+2)+0.5
  R=Ru-Rl

  W=6
  H=3.5
  ytick<-Studynames
  xtick<-seq(ceiling(Rl/2)*2,floor(Ru/2)*2+10^-8,2)

  pdf(here("FiguresandTables","GMMresultsOriginalAndAdjusted.pdf"), width = 6, height = 3.5)
  par(mgp = c(3,0.5,0), mar = c(1, 13, 2, 1))
  plot(x = Z1/R*W, y = c(1:n)/n*H,
       pch = 4, frame = FALSE, xlab="", ylab="",yaxt='n', xaxt = "n",
       xlim = c(Rl/R*W,Ru/R*W), ylim = c(0,H), col = BLUE,cex = 0.3); #original
  par(new=TRUE)
  lines(c(0,0),c(0,n+1),col=GREY, xlab="", ylab="");
  par(new=TRUE)
  lines(c(1.96/R*W,1.96/R*W),c(0,n+1), xlab="", ylab="",col=GREY);
  par(new=TRUE)
  lines(c(-1.96/R*W,-1.96/R*W),c(0,n+1), xlab="", ylab="",col=GREY);
  par(new=TRUE)
  plot(x = Z1/R*W, y = c(1:n)/n*H,
       pch = 4, frame = FALSE, xlab="", ylab="",yaxt='n', xaxt = "n",
       xlim = c(Rl/R*W,Ru/R*W), ylim = c(0,H), col = BLUE,cex = 0.3,add=TRUE); #original
  par(new=TRUE)
  plot(x = Z1_M/R*W, y = (c(1:n)-0.3)/n*H,yaxt='n',
       pch = 19, frame = FALSE, xlab="", ylab="",xaxt = "n",
       xlim = c(Rl/R*W,Ru/R*W), ylim = c(0,H),cex = 0.3); #adjusted
  par(new=TRUE)
  for (i in c(1:n)) {
    lines(c(Z1[i]-1.96,Z1[i]+1.96)/R*W,c(i,i)/n*H, xlab="", ylab="",col=BLUE,lwd=1.5);
    lines(c(Z1_L[i],Z1_U[i])/R*W,(c(i,i)-0.3)/n*H, xlab="", ylab="",lwd=1.5);
  }

  axis(2, at = seq(1/n,H,H/n), labels = ytick, las = 1, cex.axis=0.75)
  axis(3, at=seq(ceiling(Rl/2),floor(Ru/2),1)/R*W*2, labels = xtick, cex.axis=0.6)

  legend(2,1,legend=c("Original Estimates", "Adjusted Estimates"),col=c(BLUE,"black"), pch=c(4,19), cex=0.5)

  dev.off()


  # Plot Original and adjusted confidence sets, including bonferroni  -------

  Rl=min(min(Z1_LB),min(Z1)-2)-0.5
  Ru=max(max(Z1_UB),max(Z1)+2)+0.5
  R=Ru-Rl

  W=6
  H=3.5
  ytick<-Studynames
  xtick<-seq(ceiling(Rl/2)*2,floor(Ru/2)*2+10^-8,2)

  pdf(here("FiguresandTables","GMMresultsOriginalAndAdjustedBonferroni.pdf"), width = 6, height = 3.5)
  par(mgp = c(3,0.5,0), mar = c(1, 13, 2, 1))
  par(new=TRUE)
  plot(x = Z1/R*W, y = c(1:n)/n*H,
       pch = 4, frame = FALSE, xlab="", ylab="",yaxt='n', xaxt = "n",
       xlim = c(Rl/R*W,Ru/R*W), ylim = c(1,n)*H/n, col = BLUE,cex = 0.3); #original
  lines(c(0,0),c(0,n+1),col=GREY, xlab="", ylab="");
  par(new=TRUE)
  lines(c(1.96/R*W,1.96/R*W),c(0,n+1), xlab="", ylab="",col=GREY);
  par(new=TRUE)
  lines(c(-1.96/R*W,-1.96/R*W),c(0,n+1), xlab="", ylab="",col=GREY);
  par(new=TRUE)
  plot(x = Z1/R*W, y = c(1:n)/n*H,
       pch = 4, frame = FALSE, xlab="", ylab="",yaxt='n', xaxt = "n",
       xlim = c(Rl/R*W,Ru/R*W), ylim = c(1,n)*H/n, col = BLUE,cex = 0.3,add=TRUE); #original
  par(new=TRUE)
  plot(x = Z1_M/R*W, y = (c(1:n)-0.35)/n*H,yaxt='n',
       pch = 19, frame = FALSE, xlab="", ylab="",xaxt = "n",
       xlim = c(Rl/R*W,Ru/R*W), ylim = c(1,n)*H/n ,cex = 0.3); #adjusted
  par(new=TRUE)
  for (i in c(1:n)) {
    lines(c(Z1[i]-1.96,Z1[i]+1.96)/R*W,c(i,i)/n*H, xlab="", ylab="",col=BLUE,lwd=1.5);
    lines(c(Z1_L[i],Z1_U[i])/R*W,(c(i,i)-0.3)/n*H, xlab="", ylab="",lwd=1.5);

    lines(c(1,1)*Z1_LB[i]/R*W,c(i-0.4,i-0.2)/n*H, xlab="", ylab="");
    lines(c(1,1)*Z1_UB[i]/R*W,(c(i-0.4,i-0.2))/n*H, xlab="", ylab="");
  }

  axis(2, at = seq(1/n,H,H/n), labels = ytick, las = 1,cex.axis=0.75)
  axis(3, at=seq(ceiling(Rl/2),floor(Ru/2),1)/R*W*2, labels = xtick,cex.axis=0.6)

  legend(2,1,legend=c("Original Estimates", "Adjusted Estimates"),col=c(BLUE,"black"), pch=c(4,19), cex=0.5)

  dev.off()


  # Plot corrected estimators and confidence sets ---------------------------

  if (symmetric==1) {
    xgrid=seq(0,5,0.01)
    alpha=0.05

    Theta_U_store=matrix(0,length(xgrid),1)
    Theta_L_store=matrix(0,length(xgrid),1)
    Theta_M_store=matrix(0,length(xgrid),1)

    theta_UB=matrix(0,length(xgrid),1)
    theta_LB=matrix(0,length(xgrid),1)

    sigma_thetaU=matrix(0,length(xgrid),1)
    sigma_thetaL=matrix(0,length(xgrid),1)

    for (n in c(1:length(xgrid))) {
      X=xgrid[n]

      g_U <- function(lambda) {
        (alpha/2-Step_function_normal_cdf(X,lambda,1,cbind(Psihat_use,1),cutoffs,symmetric))^2
      }
      g_L <- function(lambda) {
        (1-alpha/2-Step_function_normal_cdf(X,lambda,1,cbind(Psihat_use,1),cutoffs,symmetric))^2
      }
      g_M <- function(lambda) {
        (1/2-Step_function_normal_cdf(X,lambda,1,cbind(Psihat_use,1),cutoffs,symmetric))^2
      }

      min_U<-optim(par=X,fn=g_U)
      min_U<-optim(par=min_U$par,fn=g_U,method="BFGS")
      min_L<-optim(par=X,fn=g_L)
      min_L<-optim(par=min_L$par,fn=g_L,method="BFGS")
      min_M<-optim(par=X,fn=g_M)
      min_M<-optim(par=min_M$par,fn=g_M,method="BFGS")

      theta_U=min_U$par
      theta_L=min_L$par
      theta_M=min_M$par

      Theta_U_store[n,1]=theta_U
      Theta_L_store[n,1]=theta_L
      Theta_M_store[n,1]=theta_M

      #Calculate Bonferroni corrected estimates and confidence bounds
      g_UB <- function(lambda) {
        (bonf_alpha/2-Step_function_normal_cdf(X,lambda,1,cbind(Psihat_use,1),cutoffs,symmetric))^2
      }
      g_LB <- function(lambda) {
        (1-bonf_alpha/2-Step_function_normal_cdf(X,lambda,1,cbind(Psihat_use,1),cutoffs,symmetric))^2
      }

      min_UB<-optim(par=Theta_U_store[n,1],fn=g_UB)
      min_UB<-optim(par=min_UB$par,fn=g_UB,method="BFGS")
      min_LB<-optim(par=Theta_L_store[n,1],fn=g_LB)
      min_LB<-optim(par=min_LB$par,fn=g_LB,method="BFGS")

      theta_UB[n,1]=min_UB$par
      theta_LB[n,1]=min_LB$par

      #Calculate derivatives of corrections with respect to parameters via implicit function theorem
      thetaU_plus=theta_UB[n,1]+stepsize
      thetaU_minus=theta_UB[n,1]-stepsize
      F_Uplus=Step_function_normal_cdf(X,thetaU_plus,1,cbind(Psihat_use,1),cutoffs,symmetric)
      F_Uminus=Step_function_normal_cdf(X,thetaU_minus,1,cbind(Psihat_use,1),cutoffs,symmetric)
      dFUdtheta=(F_Uplus-F_Uminus)/(2*stepsize)

      thetaL_plus=theta_LB[n,1]+stepsize
      thetaL_minus=theta_LB[n,1]-stepsize
      F_Lplus=Step_function_normal_cdf(X,thetaL_plus,1,cbind(Psihat_use,1),cutoffs,symmetric)
      F_Lminus=Step_function_normal_cdf(X,thetaL_minus,1,cbind(Psihat_use,1),cutoffs,symmetric)
      dFLdtheta=(F_Lplus-F_Lminus)/(2*stepsize)

      dFUdbeta=matrix(0,length(Psihat_use),1)
      dFLdbeta=matrix(0,length(Psihat_use),1)

      for (n1 in c(1:length(Psihat_use))) {
        Psi_plus=Psihat_use
        Psi_plus[n1]=Psi_plus[n1]+stepsize
        Psi_minus=Psihat_use
        Psi_minus[n1]=Psi_minus[n1]-stepsize

        F_Uplus=Step_function_normal_cdf(X,theta_UB[n,1],1,cbind(Psi_plus,1),cutoffs,symmetric)
        F_Uminus=Step_function_normal_cdf(X,theta_UB[n,1],1,cbind(Psi_minus,1),cutoffs,symmetric)
        dFUdbeta[n1,1]=(F_Uplus-F_Uminus)/(2*stepsize)

        F_Lplus=Step_function_normal_cdf(X,theta_LB[n,1],1,cbind(Psi_plus,1),cutoffs,symmetric)
        F_Lminus=Step_function_normal_cdf(X,theta_LB[n,1],1,cbind(Psi_minus,1),cutoffs,symmetric)
        dFLdbeta[n1,1]=(F_Lplus-F_Lminus)/(2*stepsize)
      }

      dFUdtheta=as.numeric(dFUdtheta)
      dFLdtheta=as.numeric(dFLdtheta)

      dthetaUdbeta=-dFUdbeta/dFUdtheta
      dthetaLdbeta=-dFLdbeta/dFLdtheta

      sigma_thetaU[n,1]=(t(dthetaUdbeta)%*%Varhat%*%dthetaUdbeta)^0.5
      sigma_thetaL[n,1]=(t(dthetaLdbeta)%*%Varhat%*%dthetaLdbeta)^0.5
    }

    Theta_UB_store=theta_UB+qnorm(1-bonf_beta)*sigma_thetaU
    Theta_LB_store=theta_LB-qnorm(1-bonf_beta)*sigma_thetaL
  } else {
    xgrid=seq(-5,5,0.1)
    alpha=0.05

    Theta_U_store=matrix(0,length(xgrid),1)
    Theta_L_store=matrix(0,length(xgrid),1)
    Theta_M_store=matrix(0,length(xgrid),1)

    theta_UB=matrix(0,length(xgrid),1)
    theta_LB=matrix(0,length(xgrid),1)

    sigma_thetaU=matrix(0,length(xgrid),1)
    sigma_thetaL=matrix(0,length(xgrid),1)

    for (n in c(1:length(xgrid))) {
      X=xgrid[n]

      g_U <- function(lambda) {
        (alpha/2-Step_function_normal_cdf(X,lambda,1,cbind(Psihat_use,1),cutoffs,symmetric))^2
      }
      g_L <- function(lambda) {
        (1-alpha/2-Step_function_normal_cdf(X,lambda,1,cbind(Psihat_use,1),cutoffs,symmetric))^2
      }
      g_M <- function(lambda) {
        (1/2-Step_function_normal_cdf(X,lambda,1,cbind(Psihat_use,1),cutoffs,symmetric))^2
      }

      min_U<-optim(par=X,fn=g_U)
      min_U<-optim(par=min_U$par,fn=g_U,method="BFGS")
      min_L<-optim(par=X,fn=g_L)
      min_L<-optim(par=min_L$par,fn=g_L,method="BFGS")
      min_M<-optim(par=X,fn=g_M)
      min_M<-optim(par=min_M$par,fn=g_M,method="BFGS")

      theta_U=min_U$par
      theta_L=min_L$par
      theta_M=min_M$par

      Theta_U_store[n,1]=theta_U
      Theta_L_store[n,1]=theta_L
      Theta_M_store[n,1]=theta_M

      #Calculate Bonferroni corrected estimates and confidence bounds
      g_UB <- function(lambda) {
        (bonf_alpha/2-Step_function_normal_cdf(X,lambda,1,cbind(Psihat_use,1),cutoffs,symmetric))^2
      }
      g_LB <- function(lambda) {
        (1-bonf_alpha/2-Step_function_normal_cdf(X,lambda,1,cbind(Psihat_use,1),cutoffs,symmetric))^2
      }

      min_UB<-optim(par=Theta_U_store[n,1],fn=g_UB)
      min_UB<-optim(par=min_UB$par,fn=g_UB,method="BFGS")
      min_LB<-optim(par=Theta_L_store[n,1],fn=g_LB)
      min_LB<-optim(par=min_LB$par,fn=g_LB,method="BFGS")

      theta_UB[n,1]=min_UB$par
      theta_LB[n,1]=min_LB$par

      #Calculate derivatives of corrections with respect to parameters via implicit function theorem
      thetaU_plus=theta_UB[n,1]+stepsize
      thetaU_minus=theta_UB[n,1]-stepsize
      F_Uplus=Step_function_normal_cdf(X,thetaU_plus,1,cbind(Psihat_use,1),cutoffs,symmetric)
      F_Uminus=Step_function_normal_cdf(X,thetaU_minus,1,cbind(Psihat_use,1),cutoffs,symmetric)
      dFUdtheta=(F_Uplus-F_Uminus)/(2*stepsize)

      thetaL_plus=theta_LB[n,1]+stepsize
      thetaL_minus=theta_LB[n,1]-stepsize
      F_Lplus=Step_function_normal_cdf(X,thetaL_plus,1,cbind(Psihat_use,1),cutoffs,symmetric)
      F_Lminus=Step_function_normal_cdf(X,thetaL_minus,1,cbind(Psihat_use,1),cutoffs,symmetric)
      dFLdtheta=(F_Lplus-F_Lminus)/(2*stepsize)

      dFUdbeta=matrix(0,length(Psihat_use),1)
      dFLdbeta=matrix(0,length(Psihat_use),1)

      for (n1 in c(1:length(Psihat_use))) {
        Psi_plus=Psihat_use
        Psi_plus[n1]=Psi_plus[n1]+stepsize
        Psi_minus=Psihat_use
        Psi_minus[n1]=Psi_minus[n1]-stepsize

        F_Uplus=Step_function_normal_cdf(X,theta_UB[n,1],1,cbind(Psi_plus,1),cutoffs,symmetric)
        F_Uminus=Step_function_normal_cdf(X,theta_UB[n,1],1,cbind(Psi_minus,1),cutoffs,symmetric)
        dFUdbeta[n1,1]=(F_Uplus-F_Uminus)/(2*stepsize)

        F_Lplus=Step_function_normal_cdf(X,theta_LB[n,1],1,cbind(Psi_plus,1),cutoffs,symmetric)
        F_Lminus=Step_function_normal_cdf(X,theta_LB[n,1],1,cbind(Psi_minus,1),cutoffs,symmetric)
        dFLdbeta[n1,1]=(F_Lplus-F_Lminus)/(2*stepsize)
      }

      dFUdtheta=as.numeric(dFUdtheta)
      dFLdtheta=as.numeric(dFLdtheta)

      dthetaUdbeta=-dFUdbeta/dFUdtheta
      dthetaLdbeta=-dFLdbeta/dFLdtheta

      sigma_thetaU[n,1]=(t(dthetaUdbeta)%*%Varhat%*%dthetaUdbeta)^0.5
      sigma_thetaL[n,1]=(t(dthetaLdbeta)%*%Varhat%*%dthetaLdbeta)^0.5
    }

    Theta_UB_store=theta_UB+qnorm(1-bonf_beta)*sigma_thetaU
    Theta_LB_store=theta_LB-qnorm(1-bonf_beta)*sigma_thetaL

  }


  # The plot ----------------------------------------------------------------

  xmin=min(xgrid)
  xmax=max(xgrid)

  ymin=min(min(Theta_LB_store,xmin-1.96))
  ymax=max(max(Theta_UB_store,xmax+1.96))

  pdf(here("FiguresandTables","GMMresultsCorrection_plot.pdf"), width = 4, height = 4)
  par(mar = c(3, 3, 1, 1))
  plot(x = xgrid, y = Theta_L_store,
       lty = 1, xlab='', ylab='',
       xlim=c(xmin,xmax), ylim=c(ymin,ymax),
       col = BLUE,cex=0.1)
  title(ylab=TeX('Estimation $\\theta$'), xlab=TeX('$X$'), line=1.75, cex.lab=0.8)
  par(new=TRUE)
  for (n in seq(ceiling(xmin)+1,floor(xmax)-1,1)) {
    lines(c(n,n),c(ymin,ymax),lty=3,col=GREY,lwd=1.5)
  }
  par(new=TRUE)
  for (n in seq(ceiling(ymin),floor(ymax),1)) {
    lines(c(xmin,xmax),c(n,n),lty=3,col=GREY,lwd=1.5)
  }
  par(new=TRUE)
  lines(xgrid,xgrid-1.96,col=GREY, xlab="", ylab="",lwd=1.5);
  par(new=TRUE)
  lines(xgrid,xgrid, xlab="", ylab="",col=GREY,lwd=1.5);
  par(new=TRUE)
  lines(xgrid,xgrid+1.96, xlab="", ylab="",col=GREY,lwd=1.5);
  par(new=TRUE)
  lines(x = xgrid, y = Theta_L_store,
        lty = 1, frame = FALSE, xlab=TeX('$X$'), ylab=TeX('Estimation $\\theta$'),
        xlim=c(xmin,xmax), ylim=c(ymin,ymax),
        col = BLUE,lwd=1.5);
  par(new=TRUE)
  lines(x = xgrid, y = Theta_M_store,
        lty = 1, frame = FALSE, xlab="", ylab="",
        col = BLUE,lwd=3);
  par(new=TRUE)
  lines(x = xgrid, y = Theta_LB_store,
        lty = 1, frame = FALSE, xlab="", ylab="",
        col = "black");
  par(new=TRUE)
  lines(x = xgrid, y = Theta_UB_store,
        lty = 1, frame = FALSE, xlab="", ylab="",
        col = "black");
  lines(x = xgrid, y = Theta_U_store,
        lty = 1, frame = FALSE, xlab="", ylab="",
        col = BLUE,lwd=1.5);

  legend(0,7,legend=c("95% confidence Bounds", "Bonferroni Corrected 95% Bounds", "Median Unbiased Estimator"),col=c(BLUE,"black",BLUE), lty=1, lwd=c(1.5,1,3), cex=0.6)

  dev.off()
}

PublicationbiasGMM <- function(X,sigma,cluster_ID,symmetric,cutoffs,Studynames) {
  #PublicationbiasGMM
  #Arguments the user has to define:
  #X: studies estimates
  #sigma: studies standard deviations
  #symmetric: dummy, whether p is symmetric around 0
  #cluster_ID: stuides for clustered standard errors
  #cutoffs: a COLUMN vector of cutoffs to use in step function,
  #should be given in an increasing order
  #Studynames: the names of the studies the estimates and standard errors are drawn from
  #Note that it is possible to just plug in numbers. However, the values must be
  #of mode character. For defult you may use
  #Studynames <- vector(mode="character", length=length(X)). This will create an empty vector
  #of size n with blank characters

  #Starting values for optimization (can be modified by the researcher)
  Psihat0 <<- matrix(1,1,length(cutoffs)); #betap
  Psihat0_theta <<- matrix(0:1,1,2); #[theta, sd(theta)]

  #Primitives
  n=length(X);
  C=matrix(1,n,1);
  name='GMMResults';

  includeinfigure <<- array(matrix(1,n,1));
  includeinestimation <<- array(matrix(1,n,1));

  #If you want to have the figues that specify where the estimates fall with
  #respect to significance in a 5% level
  dofigures <<- 1;

  #If you want to get the step function p values at the points of discontinuity
  doestimates <<- 1;

  #If you want to get robust confidence sets
  dorobust <<- 1;

  #If you want the corrections and figures that accompany the corrections
  docorrections=1;

  # Producing figures -------------------------------------------------------
  if (dofigures==1) {
    DescriptiveStats(X,sigma)
  }


  # Estimating the model ----------------------------------------------------
  if (doestimates==1) {
    Estimates=EstimatingSelection(X,sigma,symmetric,cluster_ID,cutoffs,Studynames)
    Psihat=Estimates$Psihat
    Varhat=Estimates$Varhat
  }


  # Producing bias-corrected estimates and confidence sets ------------------

  if (docorrections==1) {
    HorizontalBars(X,sigma,Psihat,Varhat,cutoffs,symmetric,Studynames)
  }

}
