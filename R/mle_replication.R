mle_replication <-
  function(X,
           sigma,
           symmetric,
           cluster_ID,
           cutoffs,
           Studynames) {
    if (spec_test == 1) {
      Psihat0<-cbind(0,1, rep(1,length(cutoffs)))
      LLH <- function(Psi) {
        ReplicationpzthetaLogLikelihood(
          Psi[1],
          Psi[2],
          c(Psi[-c(1, 2)],  1),
          c(zeros(1, length(Psi[-c(1, 2)])), 0),
          cutoffs,
          symmetric,
          Z,
          sigmaZ2,
          vcutoff
        )
      }

    } else {
      Psihat0<-cbind(0,1, rep(1,length(cutoffs)))
      LLH <- function(Psi) {
        ReplicationAnalyticLogLikelihoodControls(
          Psi[1],
          Psi[2],
          c(reshape(t(Psi[-c(1,2)]), c(length(Psi[-c(1,2)]) / length(cutoffs), length(cutoffs))), 1),
          cutoffs,
          symmetric,
          Z,
          sigmaZ2,
          C,
          numerical_integration
        )
      }
    }

    nn = n

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
