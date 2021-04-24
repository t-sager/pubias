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
mle_replication <- function(Z, sigmaZ2, symmetric, cluster_ID, cutoffs, studynames) {

    # Stepsize
    stepsize <- 10^(-6)

    # Starting Values
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
    Psihat0<-c(1,1, rep(1,length(cutoffs)))
    mini <- optim(par=Psihat0,fn=LLH_only,method="BFGS",control = list(abstol=10^-8,maxit=10^5))

    Psihat <- mini$par
    Objval <- mini$value

    Varhat <- robust_variance(stepsize, nn, Psihat, LLH,cluster_ID)
    se_robust <- sqrt(diag(Varhat))

    return(list("Psihat"= Psihat, "Varhat" = Varhat, "se_robust" = se_robust))
  }
