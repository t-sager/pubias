#' Title
#'
#' @param X
#' @param sigma
#' @param symmetric
#' @param cluster_ID
#' @param cutoffs
#' @param Studynames
#'
#' @return
#' @export
#'
#' @examples
mle_meta <- function(X, sigma, symmetric, symmetric_p,cluster_ID, cutoffs, studynames, C) {

  # Stepsize
  stepsize <- 10^(-6)
  n <- length(X)

  if (symmetric_p == 1) { # symmetric

        LLH <-
          function (Psi) {
            f <- variation_variance_llh(
              Psi[1],
              Psi[2],
              c(1, Psi[-c(1,2)], rev(c(Psi[-c(1,2)])), 1),
              cutoffs,
              symmetric,
              X,
              sigma,
              C
            )
            return(f)
          }
  } else { # not symmetric

        #Run model with normal distribution for latent effects
        LLH <-
          function (Psi) {
            f <- variation_variance_llh(
              Psi[1],
              Psi[2],
              c(reshape(t(Psi[-c(1,2)]), c(length(Psi[-c(1,2)]) / length(cutoffs), length(cutoffs))), 1),
              cutoffs,
              symmetric,
              X,
              sigma,
              C
            )
            return(f)
          }
  }

  nn <- n

    LLH_only<-function (Psi){
      f<-LLH(Psi);
      return(f$LLH)
  }

    #########################
    # find maximum likelihood estimator using just LLH
    Psihat0<-c(1,1,rep(1,length(cutoffs)))
    mini <- optim(par=Psihat0,fn=LLH_only,method="BFGS",control = list(abstol=10^-8,maxit=10^5))

    Psihat <<- mini$par
    Objval <- mini$value

    Varhat <- robust_variance(stepsize, nn, Psihat, LLH,cluster_ID)
    se_robust <- sqrt(diag(Varhat))

    return(list("Psihat"= Psihat[-c(1,2)], "Varhat" = Varhat[-c(1,2), -c(1,2)], "se_robust" = se_robust[-c(1,2)]))
}



