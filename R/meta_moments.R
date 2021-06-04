meta_moments <- function(betap, cutoffs, symmetric, X, sigma) {

    # Preliminaries
    n <- length(X)
    T <- X / sigma

    # Create step function for publication probability
    Tpowers <- zeros(n, length(cutoffs) + 1)

    if (symmetric == TRUE) {
        Tpowers[, 1] <- ifelse(abs(T) < cutoffs[1], 1, 0)
        if (length(cutoffs) > 1) {
            for (m in (2:length(cutoffs))) {
                Tpowers[, m] <- (ifelse(abs(T) < cutoffs[m], 1, 0)) * (ifelse(abs(T) >= cutoffs[m-1], 1, 0))

            }
        }
        Tpowers[, length(cutoffs) + 1] <- ifelse(abs(T) >= cutoffs[length(cutoffs)], 1, 0)

    } else {
        Tpowers[, 1] <- ifelse(T < cutoffs[1], 1, 0)

        if (length(cutoffs) > 1) {
            for (m in c(2:length(cutoffs))) {
                Tpowers[, m] <- (ifelse(T < cutoffs[m], 1, 0)) * (ifelse(T >= cutoffs[m-1], 1, 0))

            }
        }
        Tpowers[, length(cutoffs) + 1] <- ifelse(T >= cutoffs[length(cutoffs)], 1, 0)
    }

    phat <- Tpowers %*% t(t(betap))

########################
# Create moments
  Xmat1 <- matrix(X, length(X), length(X))
  Xmat2 <- t(Xmat1)
  sigmamat1 <- matrix(sigma, length(X), length(X))
  sigmamat2 <- t(sigmamat1)

  indicator <- ifelse(sigmamat1 >= sigmamat2, 1, 0)

  sigmadiff <- sqrt(indicator * (sigmamat1 ^ 2 - sigmamat2 ^ 2)) +
      sqrt((1 - indicator) * (sigmamat2 ^ 2 - sigmamat1 ^ 2))
  sigmadiff <- Re(sigmadiff)
  pmat1 <- matrix(phat, length(X), length(X))
  pmat2 <- t(pmat1)

  moment_mean <- zeros(1, length(cutoffs))
  rhat <- zeros(n, length(cutoffs))

  if (symmetric == TRUE) {
      for (k in (1:length(cutoffs))) {
          c <- cutoffs[k]
          base_moments <- (pmat2 ^ -1) * (pmat1 ^ -1) * indicator *
              ((pnorm((c * sigmamat1 - Xmat2) / sigmadiff
              ) - pnorm((-c * sigmamat1 - Xmat2) / sigmadiff
              )) - (((
                  ifelse(Xmat1 <= c * sigmamat1, 1, 0)
              ) - (
                  ifelse(Xmat1 <= -c * sigmamat1, 1, 0)
              )))) + (pmat2 ^ -1) * (pmat1 ^ -1) * (1 - indicator) * ((pnorm((c * sigmamat2 -
                                                                                  Xmat1) / sigmadiff
              ) - pnorm((-c * sigmamat2 - Xmat1) / sigmadiff
              )) - (((
                  ifelse(Xmat2 <= c * sigmamat2, 1, 0)
              ) - (
                  ifelse(Xmat2 <= -c * sigmamat2, 1, 0)
              ))))

          base_moments <- base_moments - diag(diag(base_moments))
          base_moments[is.nan(base_moments)] <- 0
          moment_mean[1, k] <- pracma::nchoosek(n,2)^(-1)/2*sum(sum(base_moments))

          rhat[, k] <- 2 * ((n - 1) ^ -1) * as.matrix(rowSums(base_moments, na.rm = FALSE))
      }
  } else {
      for (k in c(1:length(cutoffs))) {
          c <- cutoffs[k]
          base_moments <- (pmat2 ^ -1) * (pmat1 ^ -1) * indicator *
              (pnorm((c * sigmamat1 - Xmat2) / sigmadiff) - (ifelse(Xmat1 <= c * sigmamat1, 1, 0))) +
              (pmat2 ^ -1) * (pmat1 ^ -1) * (1 - indicator) * (pnorm((c * sigmamat2 -
                                                                          Xmat1) / sigmadiff) - (ifelse(Xmat2 <= c * sigmamat2, 1, 0)))
          base_moments <- base_moments - diag(diag(base_moments))
          base_moments[is.nan(base_moments)] <- 0
          moment_mean[1, k] <- pracma::nchoosek(n,2)^(-1)/2*sum(sum(base_moments))
          rhat[, k] <- 2 * ((n - 1) ^ -1) * as.matrix(rowSums(base_moments, na.rm <- FALSE))
      }
  }

  moment_var <- cov(rhat)
  raw_moments <- rhat

  # Return moments results
  return(
      list(
          "moment_mean" = moment_mean,
          "moment_var" = moment_var,
          "raw_moments" = raw_moments
      )
  )
}
