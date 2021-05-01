#' Title
#'
#' @param X A `n x 1` matrix containing the estimates, where `n` is the number of estimates.
#' @param Z A `n x 2` matrix where the first (second) column contains the standardized original estimates (replication estimates), where `n` is the number of estimates.
#' @param sigma A `n x 1` matrix containing the standard errors of the estimates, where `n` is the number of estimates.
#' @param result List obJect which contains the results (`Psihat`, `Varhat`, `se_robust`) of the publication probability estimation.
#' @param cutoffs A matrix containing the thresholds for the steps of the publication probability. Should be strictly increasing column
#' vector of size `k x 1` where `k` is the number of cutoffs.
#' @param symmetric If set to `1`, the publication probability is assumed to be symmetric around zero. If set to `0`, asymmetry is allowed.
#' @param symmetric_p If set to `1`, the publication probability is assumed to be symmetric around zero and all cutoffs should be positive.
#' If set to `0`, asymmetry is allowed and cutoffs should be specified in increasing order.
#' @param identificationapproach Indication if we are dealing with replication studies (== 1) or a met analysis (== 2).
#' @param GMM If set to TRUE, the publication probability will be estimated via GMM. Setting it to FALSE uses the MLE
#' method for estimation.
#'
#' @return Returns a list containing the original estimates with their 95% confidence bonds (`Z1`, `Z1_L`, `Z1_U`)
#' as well as the corrected estimates (`Z1_M`) and in addition the Bonferroni corrected 95% confidence bonds (`Z1_LB`, `Z1_UB`).
#' There are additional elements which are mainly used for plotting the results.
#' @export
#'
bias_correction <- function(X,Z,sigma,result,cutoffs,symmetric,symmetric_p,identificationapproach, GMM) {

  Psihat <- result$Psihat
  Varhat <- result$Varhat

  #normalize estimates from initial study
  if (identificationapproach==1){
    Z1 <- as.matrix(Z[,1])
  } else if (identificationapproach==2){
    Z1 <- as.matrix(sort(X/sigma))
  }

  n <- length(Z1)
  alpha <- 0.05
  bonf_alpha <- 0.04
  bonf_beta <- 0.01

  # Pad Psihat_use with additional zeros in symmetric specifications to make dimensions conform
  if (symmetric==1){
    Psihat_use <- Psihat
  } else {
    if (symmetric_p==1){
      Psihat_use <- cbind(Psihat[1], Psihat[2], 1, Psihat[3:length(Psihat)], fliplr(Psihat[3:length(Psihat)]))
    } else {
        Psihat_use <- Psihat
      }
    }


  # Diff between MLE and GMM
  if (GMM == TRUE) {
    Psihat_use <- Psihat_use
    Varhat <- Varhat
  } else {
    Psihat_use <- Psihat_use[-c(1,2)]
    Varhat <- Varhat[-c(1,2), -c(1,2)]
  }

  # Calculate corrected estimates and confidence bounds
  Z1_U <- zeros(length(Z1),1)
  Z1_L <- zeros(length(Z1),1)
  Z1_M <- zeros(length(Z1),1)

  stepsize <- 10^-3
  for (n in (1:length(Z1))) {
    g_U <- function(lambda) {
      (alpha/2-step_function_normal_cdf(Z1[n],lambda,1,c(Psihat_use,1),cutoffs,symmetric))^2
    }
    g_L <- function(lambda) {
      (1-alpha/2-step_function_normal_cdf(Z1[n],lambda,1,c(Psihat_use,1),cutoffs,symmetric))^2
    }
    g_M <- function(lambda) {
      (1/2-step_function_normal_cdf(Z1[n],lambda,1,c(Psihat_use,1),cutoffs,symmetric))^2
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

  # Calculate bonferroni corrected estimates and confidence bounds
  Z1_UB <- zeros(length(Z1),1)
  Z1_LB <- zeros(length(Z1),1)
  sigma_U <- zeros(length(Z1),1)
  sigma_L <- zeros(length(Z1),1)

  for (n in (1:length(Z1))) {
    g_UB <- function(lambda) {
      (bonf_alpha/2-step_function_normal_cdf(Z1[n],lambda,1,c(Psihat_use,1),cutoffs,symmetric))^2
    }
    g_LB <- function(lambda) {
      (1-bonf_alpha/2-step_function_normal_cdf(Z1[n],lambda,1,c(Psihat_use,1),cutoffs,symmetric))^2
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
    F_Uplus=step_function_normal_cdf(Z1[n],thetaU_plus,1,c(Psihat_use,1),cutoffs,symmetric)
    F_Uminus=step_function_normal_cdf(Z1[n],thetaU_minus,1,c(Psihat_use,1),cutoffs,symmetric)
    dFUdtheta=(F_Uplus-F_Uminus)/(2*stepsize)

    thetaL_plus=Z1_LB[n,1]+stepsize
    thetaL_minus=Z1_LB[n,1]-stepsize
    F_Lplus=step_function_normal_cdf(Z1[n],thetaL_plus,1,c(Psihat_use,1),cutoffs,symmetric)
    F_Lminus=step_function_normal_cdf(Z1[n],thetaL_minus,1,c(Psihat_use,1),cutoffs,symmetric)
    dFLdtheta=(F_Lplus-F_Lminus)/(2*stepsize)

    dFUdbeta <- zeros(length(Psihat_use),1)
    dFLdbeta <- zeros(length(Psihat_use),1)
    for (n1 in (1:length(Psihat_use))) {
      Psi_plus=Psihat_use
      Psi_plus[n1]=Psi_plus[n1]+stepsize
      Psi_minus=Psihat_use
      Psi_minus[n1]=Psi_minus[n1]-stepsize

      F_Uplus=step_function_normal_cdf(Z1[n],Z1_UB[n,1],1,c(Psi_plus,1),cutoffs,symmetric)
      F_Uminus=step_function_normal_cdf(Z1[n],Z1_UB[n,1],1,c(Psi_minus,1),cutoffs,symmetric)
      dFUdbeta[n1,1]=(F_Uplus-F_Uminus)/(2*stepsize)

      F_Lplus=step_function_normal_cdf(Z1[n],Z1_LB[n,1],1,c(Psi_plus,1),cutoffs,symmetric)
      F_Lminus=step_function_normal_cdf(Z1[n],Z1_LB[n,1],1,c(Psi_minus,1),cutoffs,symmetric)
      dFLdbeta[n1,1]=(F_Lplus-F_Lminus)/(2*stepsize)
    }

    dFUdtheta <- as.numeric(dFUdtheta)
    dFLdtheta <- as.numeric(dFLdtheta)

    dthetaUdbeta <- -dFUdbeta/dFUdtheta
    dthetaLdbeta <- -dFLdbeta/dFLdtheta

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

  Z1_UB_store <- Z1_UB
  Z1_LB_store <- Z1_LB
  Z1_UB <- Z1_UB+qnorm(1-bonf_beta)*sigma_U
  Z1_LB <- Z1_LB-qnorm(1-bonf_beta)*sigma_L

  count=1
  store=matrix(0,length(seq(-10,10,0.1)),1)
  for (lambda in seq(-10,10,0.1)) {
    store[count]=g_M(lambda)
    count=count+1
  }


  #######################


if (symmetric == 1 | symmetric_p == 1) {

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
        (alpha/2-step_function_normal_cdf(X,lambda,1,cbind(Psihat_use,1),cutoffs,symmetric))^2
      }
      g_L <- function(lambda) {
        (1-alpha/2-step_function_normal_cdf(X,lambda,1,cbind(Psihat_use,1),cutoffs,symmetric))^2
      }
      g_M <- function(lambda) {
        (1/2-step_function_normal_cdf(X,lambda,1,cbind(Psihat_use,1),cutoffs,symmetric))^2
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
        (bonf_alpha/2-step_function_normal_cdf(X,lambda,1,cbind(Psihat_use,1),cutoffs,symmetric))^2
      }
      g_LB <- function(lambda) {
        (1-bonf_alpha/2-step_function_normal_cdf(X,lambda,1,cbind(Psihat_use,1),cutoffs,symmetric))^2
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
      F_Uplus=step_function_normal_cdf(X,thetaU_plus,1,cbind(Psihat_use,1),cutoffs,symmetric)
      F_Uminus=step_function_normal_cdf(X,thetaU_minus,1,cbind(Psihat_use,1),cutoffs,symmetric)
      dFUdtheta=(F_Uplus-F_Uminus)/(2*stepsize)

      thetaL_plus=theta_LB[n,1]+stepsize
      thetaL_minus=theta_LB[n,1]-stepsize
      F_Lplus=step_function_normal_cdf(X,thetaL_plus,1,cbind(Psihat_use,1),cutoffs,symmetric)
      F_Lminus=step_function_normal_cdf(X,thetaL_minus,1,cbind(Psihat_use,1),cutoffs,symmetric)
      dFLdtheta=(F_Lplus-F_Lminus)/(2*stepsize)

      dFUdbeta=matrix(0,length(Psihat_use),1)
      dFLdbeta=matrix(0,length(Psihat_use),1)

      for (n1 in c(1:length(Psihat_use))) {
        Psi_plus=Psihat_use
        Psi_plus[n1]=Psi_plus[n1]+stepsize
        Psi_minus=Psihat_use
        Psi_minus[n1]=Psi_minus[n1]-stepsize

        F_Uplus=step_function_normal_cdf(X,theta_UB[n,1],1,cbind(Psi_plus,1),cutoffs,symmetric)
        F_Uminus=step_function_normal_cdf(X,theta_UB[n,1],1,cbind(Psi_minus,1),cutoffs,symmetric)
        dFUdbeta[n1,1]=(F_Uplus-F_Uminus)/(2*stepsize)

        F_Lplus=step_function_normal_cdf(X,theta_LB[n,1],1,cbind(Psi_plus,1),cutoffs,symmetric)
        F_Lminus=step_function_normal_cdf(X,theta_LB[n,1],1,cbind(Psi_minus,1),cutoffs,symmetric)
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
        (alpha/2-step_function_normal_cdf(X,lambda,1,cbind(Psihat_use,1),cutoffs,symmetric))^2
      }
      g_L <- function(lambda) {
        (1-alpha/2-step_function_normal_cdf(X,lambda,1,cbind(Psihat_use,1),cutoffs,symmetric))^2
      }
      g_M <- function(lambda) {
        (1/2-step_function_normal_cdf(X,lambda,1,cbind(Psihat_use,1),cutoffs,symmetric))^2
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
        (bonf_alpha/2-step_function_normal_cdf(X,lambda,1,cbind(Psihat_use,1),cutoffs,symmetric))^2
      }
      g_LB <- function(lambda) {
        (1-bonf_alpha/2-step_function_normal_cdf(X,lambda,1,cbind(Psihat_use,1),cutoffs,symmetric))^2
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
      F_Uplus=step_function_normal_cdf(X,thetaU_plus,1,cbind(Psihat_use,1),cutoffs,symmetric)
      F_Uminus=step_function_normal_cdf(X,thetaU_minus,1,cbind(Psihat_use,1),cutoffs,symmetric)
      dFUdtheta=(F_Uplus-F_Uminus)/(2*stepsize)

      thetaL_plus=theta_LB[n,1]+stepsize
      thetaL_minus=theta_LB[n,1]-stepsize
      F_Lplus=step_function_normal_cdf(X,thetaL_plus,1,cbind(Psihat_use,1),cutoffs,symmetric)
      F_Lminus=step_function_normal_cdf(X,thetaL_minus,1,cbind(Psihat_use,1),cutoffs,symmetric)
      dFLdtheta=(F_Lplus-F_Lminus)/(2*stepsize)

      dFUdbeta=matrix(0,length(Psihat_use),1)
      dFLdbeta=matrix(0,length(Psihat_use),1)

      for (n1 in c(1:length(Psihat_use))) {
        Psi_plus=Psihat_use
        Psi_plus[n1]=Psi_plus[n1]+stepsize
        Psi_minus=Psihat_use
        Psi_minus[n1]=Psi_minus[n1]-stepsize

        F_Uplus=step_function_normal_cdf(X,theta_UB[n,1],1,cbind(Psi_plus,1),cutoffs,symmetric)
        F_Uminus=step_function_normal_cdf(X,theta_UB[n,1],1,cbind(Psi_minus,1),cutoffs,symmetric)
        dFUdbeta[n1,1]=(F_Uplus-F_Uminus)/(2*stepsize)

        F_Lplus=step_function_normal_cdf(X,theta_LB[n,1],1,cbind(Psi_plus,1),cutoffs,symmetric)
        F_Lminus=step_function_normal_cdf(X,theta_LB[n,1],1,cbind(Psi_minus,1),cutoffs,symmetric)
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

  return(list("Z1" = Z1,
              "Z1_U" = Z1_U,
              "Z1_L" = Z1_L,
              "Z1_M "= Z1_M,
              "Z1_UB" = Z1_UB,
              "Z1_LB" = Z1_LB,
              "Theta_U_store" = Theta_U_store,
              "Theta_L_store" = Theta_L_store,
              "Theta_M_store" = Theta_M_store,
              "Theta_UB_store" = Theta_UB_store,
              "Theta_LB_store" = Theta_LB_store,
              "xgrid" = xgrid))
}
