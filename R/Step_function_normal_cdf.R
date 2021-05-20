step_function_normal_cdf <- function(X,theta,sigma,betap,cutoffs,symmetric) {
    # Throw error if betap "too short"
    if (length(betap) != (length(cutoffs) + 1)) {
        stop("length of betap must be one greater than length of cutoffs ")
    }

    # For symmetric case, create symmetric version of cutoffs and coefficients
    if (symmetric == TRUE) {
        cutoffs_u <- zeros(length(cutoffs), 1)
        betap_u <- zeros(1, length(cutoffs))
        for (n in (1:length(cutoffs))) {
            cutoffs_u[n] <- -cutoffs[length(cutoffs) + 1 - n]
            betap_u[n] <- betap[length(betap) + 1 - n]
        }
        for (n in (1:length(cutoffs))) {
            cutoffs_u[length(cutoffs) + n] <- cutoffs[n]
            betap_u[length(cutoffs) + n] <- betap[n]
        }
        betap_u[length(betap_u) + 1] <- 1
        betap_u <- as.matrix(betap_u)
    } else {
        cutoffs_u <- t(cutoffs)
        betap_u <- t(t(betap))
    }

    # Calculate denominator in CDF
    prob_vec <- zeros(length(cutoffs_u) + 1, 1)

    for (m in (1:length(cutoffs_u))) {
        prob_vec[m + 1] = pnorm((cutoffs_u[m] - theta) / sigma)
    }

    prob_vec <- rbind(prob_vec, 1)
    mean_Z1 <- prob_vec[2:length(prob_vec), 1] - prob_vec[1:(length(prob_vec) - 1), 1]

    denominator <- t(mean_Z1) %*% betap_u

    # Calculate numerator in CDF
    cutoffs_u[length(cutoffs_u) + 1] <- Inf

    if (X <= cutoffs_u[1]) {
        numerator <- pnorm((X - theta) / sigma) * betap_u[1]
    } else {
        numerator <- pnorm((cutoffs_u[1] - theta) / sigma) * betap_u[1]
        m <- 1
        while (X > cutoffs_u[m]) {
            Xcap <- min(X, cutoffs_u[m + 1])
            numerator <-
                numerator + (pnorm((Xcap - theta) / sigma) - pnorm((cutoffs_u[m] - theta) /
                                                                       sigma)) * betap_u[m + 1]
            m <- m + 1
        }
    }

    # CDF
    cdf = numerator / denominator

    # Return CDF
    return(cdf)
}
