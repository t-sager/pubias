HorizontalBars <- function(X,sigma,Psihat,Varhat,cutoffs,symmetric,symmetric_p,Studynames,identificationapproach) {

  color<-brewer.pal(8,"Greys")
  GREY=color[3]
  SLIDE_GREY=color[5]
  color<-brewer.pal(8,"Purples")
  BLUE=color[6]

  #normalize estimates from initial study
  if (identificationapproach==1){
    Z1 <- as.matrix(Z[,1])
  } else if (identificationapproach==2){
    Z1 <- as.matrix(sort(X/sigma))
  }

  alpha <- 0.05
  bonf_alpha <- 0.04
  bonf_beta <- 0.01

  #Pad Psihat_use with additional zeros in symmetric specifications to make dimensions conform
  if (symmetric==1){
    Psihat_use <- Psihat
  } else {
    if (symmetric_p==1){
      Psihat_use <- cbind(Psihat[1], Psihat[2], 1, Psihat[3:length(Psihat)], fliplr(Psihat[3:length(Psihat)]))
    } else {
      if (asymmetric_likelihood_spec==2){
        Psihat_use <- Psihat[2:length(Psihat)]
      } else {
        Psihat_use <- Psihat
      }
    }
  }


  #Calculate corrected estimates and confidence bounds

  Z1_U <- zeros(length(Z1),1)
  Z1_L <- zeros(length(Z1),1)
  Z1_M <- zeros(length(Z1),1)

  stepsize <- 10^-3
  for (n in (1:length(Z1))) {
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
  Z1_UB <- zeros(length(Z1),1)
  Z1_LB <- zeros(length(Z1),1)
  sigma_U <- zeros(length(Z1),1)
  sigma_L <- zeros(length(Z1),1)

  for (n in (1:length(Z1))) {
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

    dFUdbeta <- zeros(length(Psihat_use),1)
    dFLdbeta <- zeros(length(Psihat_use),1)
    for (n1 in (1:length(Psihat_use))) {
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
       pch = 4,   xlab="", ylab="",yaxt='n', xaxt = "n",
       xlim = c(Rl/R*W,Ru/R*W), ylim = c(0,H), col = BLUE,cex = 0.3); #original
  par(new=TRUE)
  lines(c(0,0),c(0,n+1),col=GREY, xlab="", ylab="");
  par(new=TRUE)
  lines(c(1.96/R*W,1.96/R*W),c(0,n+1), xlab="", ylab="",col=GREY);
  par(new=TRUE)
  lines(c(-1.96/R*W,-1.96/R*W),c(0,n+1), xlab="", ylab="",col=GREY);
  par(new=TRUE)
  plot(x = Z1/R*W, y = c(1:n)/n*H,
       pch = 4,   xlab="", ylab="",yaxt='n', xaxt = "n",
       xlim = c(Rl/R*W,Ru/R*W), ylim = c(0,H), col = BLUE,cex = 0.3,add=TRUE); #original
  par(new=TRUE)
  plot(x = Z1_M/R*W, y = (c(1:n)-0.3)/n*H,yaxt='n',
       pch = 19,   xlab="", ylab="",xaxt = "n",
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
       pch = 4,   xlab="", ylab="",yaxt='n', xaxt = "n",
       xlim = c(Rl/R*W,Ru/R*W), ylim = c(1,n)*H/n, col = BLUE,cex = 0.3); #original
  lines(c(0,0),c(0,n+1),col=GREY, xlab="", ylab="");
  par(new=TRUE)
  lines(c(1.96/R*W,1.96/R*W),c(0,n+1), xlab="", ylab="",col=GREY);
  par(new=TRUE)
  lines(c(-1.96/R*W,-1.96/R*W),c(0,n+1), xlab="", ylab="",col=GREY);
  par(new=TRUE)
  plot(x = Z1/R*W, y = c(1:n)/n*H,
       pch = 4,   xlab="", ylab="",yaxt='n', xaxt = "n",
       xlim = c(Rl/R*W,Ru/R*W), ylim = c(1,n)*H/n, col = BLUE,cex = 0.3,add=TRUE); #original
  par(new=TRUE)
  plot(x = Z1_M/R*W, y = (c(1:n)-0.35)/n*H,yaxt='n',
       pch = 19,   xlab="", ylab="",xaxt = "n",
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
        lty = 1,   xlab=TeX('$X$'), ylab=TeX('Estimation $\\theta$'),
        xlim=c(xmin,xmax), ylim=c(ymin,ymax),
        col = BLUE,lwd=1.5);
  par(new=TRUE)
  lines(x = xgrid, y = Theta_M_store,
        lty = 1,   xlab="", ylab="",
        col = BLUE,lwd=3);
  par(new=TRUE)
  lines(x = xgrid, y = Theta_LB_store,
        lty = 1,   xlab="", ylab="",
        col = "black");
  par(new=TRUE)
  lines(x = xgrid, y = Theta_UB_store,
        lty = 1,   xlab="", ylab="",
        col = "black");
  lines(x = xgrid, y = Theta_U_store,
        lty = 1,   xlab="", ylab="",
        col = BLUE,lwd=1.5);

  legend(0,7,legend=c("95% confidence Bounds", "Bonferroni Corrected 95% Bounds", "Median Unbiased Estimator"),col=c(BLUE,"black",BLUE), lty=1, lwd=c(1.5,1,3), cex=0.6)

  dev.off()
}
