MaxEval<-10^5
MaxIter<-10^5
Tol<-10^(-8)
stepsize<-10^(-6)

mle_meta <-
  function(X,
           sigma,
           symmetric,
           cluster_ID,
           cutoffs,
           Studynames) {
    if (asymmetric_likelihood_spec == 2) {
      if (symmetric_p == 1) {
        Psihat0<-cbind(0,1,10, rep(1,length(cutoffs)))
        LLH <-
          function (Psi) {
            VariationVarianceLogLikelihoodtdist(
              Psi[1],
              Psi[2],
              Psi[3],
              c(1, Psi[-c(1,2,3)], rev(c(Psi[-c(1,2,3)])), 1),
              cutoffs,
              symmetric,
              X[includeinestimation],
              sigma[includeinestimation],
              C[includeinestimation, ],
              numerical_integration
            )
          }

        nn = sum(includeinestimation)

      } else {
        # Run model with t distribution for latent effects
        if (time_trend == 1) {
          Psihat0<-cbind(0,1,10, rep(1,length(cutoffs)))
          LLH <-
            function (Psi) {
              VariationVarianceLogLikelihoodtdist_logit(
                Psi[1],
                Psi[2],
                Psi[3],
                c(reshape(t(Psi[-c(1,2,3)]), c(length(Psi[-c(1,2,3)]) / length(cutoffs), length(cutoffs))), 1),
                cutoffs,
                symmetric,
                X[includeinestimation],
                sigma[includeinestimation],
                C[includeinestimation,],
                numerical_integration
              )
            }
          nn = n

        } else {
          Psihat0<-cbind(0,1,10, rep(1,length(cutoffs)))
          LLH <-
            function (Psi) {
              VariationVarianceLogLikelihoodtdist(
                Psi[1],
                Psi[2],
                Psi[3],
                c(reshape(t(Psi[-c(1,2,3)]), c(length(Psi[-c(1,2,3)]) / length(cutoffs), length(cutoffs))), 1),
                cutoffs,
                symmetric,
                X[includeinestimation],
                sigma[includeinestimation],
                C[includeinestimation,],
                numerical_integration
              )
            }
          nn = n

        }
      }
    } else {
      if (symmetric_p == 1) {
        Psihat0<-cbind(0,1, rep(1,length(cutoffs)))
        LLH <-
          function (Psi) {
            VariationVarianceLogLikelihoodControls(
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
          }
        nn = sum(includeinestimation)
      } else {
        Psihat0<-cbind(0,1, rep(1,length(cutoffs)))
        #Run model with normal distribution for latent effects
        LLH <-
          function (Psi) {
            VariationVarianceLogLikelihoodControls(
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
          }
        nn = sum(includeinestimation)
      }
    }

    #########################
    # find maximum likelihood estimator using just LLH
    if (spec_test == 1) {

      lower.b = c(-Inf,rep(0,length(Psihat0)-1))
      upper.b=rep(Inf,length(Psihat0))

      findmin<-nlminb(objective=LLH, start=Psihat0,lower=lower.b,upper=upper.b,control = list(eval.max = MaxEval, iter.max = MaxIter, abs.tol = Tol));
      Psihat<-findmin$par
      LLHmax<-findmin$objective

      Psihat_test <- c(Psihat, rep(0, length(cutoffs)))
      score_mat <- matrix(0, number , length(Psihat_test))

      for (n1 in 1:length(Psihat_test)) {
        theta_plus <- Psihat_test
        theta_plus[n1] <- theta_plus[n1] + stepsize
        plus <- LLH_test(theta_plus)
        LLH_plus <- plus$LLH
        logL_plus <- plus$logL
        theta_minus <- Psihat_test
        theta_minus[n1] <- theta_minus[n1] - stepsize
        minus <- LLH_test(theta_minus)
        LLH_minus <- minus$LLH
        logL_minus <- minus$logL
        score_mat[, n1] <- (logL_plus - logL_minus) / (2 * stepsize)
      }

      mean_score_mat <- as.matrix(apply(score_mat, 2, mean))
      score_stat <- number * t(mean_score_mat) %*% solve(cov(score_mat)) %*% mean_score_mat
      score_pval <- 1 - pchisq(score_stat, length(cutoffs))


    } else {
      lower.b = c(-Inf,rep(0,length(Psihat0)-1))
      upper.b=rep(Inf,length(Psihat0))

      findmin<-nlminb(objective=LLH, start=Psihat0,lower=lower.b,upper=upper.b,control = list(eval.max = MaxEval, iter.max = MaxIter, abs.tol = Tol));
      Psihat<-findmin$par
      LLHmax<-findmin$objective

    }

    return(list("Psihat"= Psihat, "LLHmax" = LLHmax))

}



