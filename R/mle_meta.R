mle_meta <- function(X, sigma, symmetric, cluster_ID, cutoffs, Studynames) {

  # Stepsize
  stepsize<-10^(-6)

  if (symmetric_p == 1) { # symmetric
        Psihat0<-c(0,1,rep(1,length(cutoffs)))
        LLH <-
          function (Psi) {
            f <- variation_variance_llh(
              Psi[1],
              Psi[2],
              c(1, Psi[-c(1,2)], rev(c(Psi[-c(1,2)])), 1),
              cutoffs,
              symmetric,
              X[includeinestimation],
              sigma[includeinestimation],
              C[includeinestimation, ],
              numerical_integration
            )
            return(f)
          }
        nn = sum(includeinestimation)
  } else { # not symmetric
        Psihat0<-c(0,1,rep(1,length(cutoffs)))
        #Run model with normal distribution for latent effects
        LLH <-
          function (Psi) {
            f <- variation_variance_llh(
              Psi[1],
              Psi[2],
              c(reshape(t(Psi[-c(1,2)]), c(length(Psi[-c(1,2)]) / length(cutoffs), length(cutoffs))), 1),
              cutoffs,
              symmetric,
              X[includeinestimation],
              sigma[includeinestimation],
              C[includeinestimation, ],
              numerical_integration
            )
            return(f)
          }
        nn = sum(includeinestimation)
      }

    LLH_only<-function (Psi){
      f<-LLH(Psi);
      return(f$LLH)
  }

    #########################
    # find maximum likelihood estimator using just LLH
    mini<-optim(par=Psihat0,fn=LLH_only,method="BFGS",control = list(abstol=10^-8,maxit=10^5))

    Psihat <- mini$par
    Objval <- mini$value

    Varhat <- robust_variance(stepsize, nn, Psihat, LLH,cluster_ID)
    se_robust <- sqrt(diag(Varhat))

    SelectionTable(Psihat,se_robust,name,(identificationapproach==1),symmetric)

    return(list("Psihat"= Psihat, "Varhat" = Varhat))
}



