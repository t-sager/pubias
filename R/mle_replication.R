mle_replication <- function(X, sigma, symmetric, cluster_ID, cutoffs, Studynames) {

    # Stepsize
    stepsize <- 10^(-6)

    # Starting Values
    Psihat0<-c(1,1, rep(1,length(cutoffs)))

    LLH <- function(Psi) {
        f <- replication_analytic_llh(
          Psi[1],
          Psi[2],
          c(reshape(t(Psi[-c(1,2)]), c(length(Psi[-c(1,2)]) / length(cutoffs), length(cutoffs))), 1),
          cutoffs,
          symmetric,
          Z,
          sigmaZ2,
          C
        )
        return(f)
      }


    nn = n

    LLH_only<-function (Psi){
      f<-LLH(Psi);
      return(f$LLH)
    }

    #########################
    # find maximum likelihood estimator using just LLH
    mini <- optim(par=Psihat0,fn=LLH_only,method="BFGS",control = list(abstol=10^-8,maxit=10^5))

    Psihat <- mini$par
    Objval <- mini$value

    Varhat <- robust_variance(stepsize, nn, Psihat, LLH,cluster_ID)
    se_robust <- sqrt(diag(Varhat))

    SelectionTable(Psihat,se_robust,name,(identificationapproach==1),symmetric)

    return(list("Psihat"= Psihat, "Varhat" = Varhat))
  }
