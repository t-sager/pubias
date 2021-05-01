#' Computing the publication probability in meta analyses
#'
#'`mle_meta()`  calculates the publication probability, its variance and robust standard errors
#' of replication studies by a Maximum Likelihood approach.
#'
#' @param X A `n x 1` matrix containing the estimates, where `n` is the number of estimates.
#' @param sigma A `n x 1` matrix containing the standard errors of the estimates, where `n` is the number of estimates.
#' @param symmetric If set to `1`, the publication probability is assumed to be symmetric around zero. If set to `0`, asymmetry is allowed.
#' @param symmetric_p If set to `1`, the publication probability is assumed to be symmetric around zero and all cutoffs should be positive.
#' If set to `0`, asymmetry is allowed and cutoffs should be specified in increasing order.
#' @param cluster_ID A `n x 1` matrix containing IDs going from 1 to `n`, where `n` is the number of estimates.
#' @param cutoffs A matrix containing the thresholds for the steps of the publication probability. Should be strictly increasing column
#' vector of size `k x 1` where `k` is the number of cutoffs.
#' @param studynames A vector of type `character` containing all the Studynames of size `n` in the same order as the argument `data`.
#' @param C A `n x 1` matrix with all values being 1. Controls.
#'
#' @return Returns a list object with the publication probability (`Psihat`), its variance (`Varhat`) and robust standard errors (`se_robust`).
#' @export
#'
mle_meta <- function(X, sigma, symmetric, symmetric_p, cluster_ID, cutoffs, studynames, C) {

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

    return(list("Psihat"= Psihat, "Varhat" = Varhat, "se_robust" = se_robust))
}



