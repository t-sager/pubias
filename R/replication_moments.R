replication_moments <- function(betap, cutoffs, symmetric, Z, sigmaZ2){

    # Preliminaries
    n <- nrow(Z)
    k <- length(betap)

    # Create step function for publication probability
    Z1_dummies <- matrix(0, n, length(cutoffs) + 1)

    # Throw error if cutoffs unsorted (not need in current package version --> only one cutoff allowed!)
    if (all(sort(cutoffs) != cutoffs)) {
        stop ("Unsorted cutoffs!")
    }

    if (symmetric == TRUE) {
        # Throw error if cutoffs negative
        if (!all(cutoffs >= 0)) {
            stop("Needs positive cutoff!")
        }
        Z1_dummies[, 1] <- abs(Z[, 1]) < cutoffs[1]

        if (length(cutoffs) > 1) {
            for (m in (2:length(cutoffs))) {
                Z1_dummies[, m] <- (abs(Z[, 1]) < cutoffs[m]) * (abs(Z[, 1]) >= cutoffs[m -
                                                                                            1])

            }
        }
        Z1_dummies[, length(cutoffs) + 1] <-
            abs(Z[, 1]) >= cutoffs[length(cutoffs)]

    } else {
        Z1_dummies[, 1] <- Z[, 1] < cutoffs[1]

        if (length(cutoffs) > 1) {
            for (m in 2:length(cutoffs)) {
                Z1_dummies[, m] <- (Z[, 1] < cutoffs[m]) * (Z[, 1] >= cutoffs[m - 1])

            }
        }
        Z1_dummies[, length(cutoffs) + 1] <-
            Z[, 1] >= cutoffs[length(cutoffs)]


    }

    phat <- Z1_dummies %*% t(t(betap))

########################
# Create moments

sigma_max <- max(max(sigmaZ2),1)
sigma_adj1 <- sqrt((sigma_max^2-1))
sigma_adj2 <- sqrt((sigma_max^2-sigmaZ2^2))

# Currently only implemented for symmetric specification!
if (symmetric == TRUE){
    if (length(cutoffs)==1){
        c <- cutoffs
        lmoment1 <- (pnorm((c-Z[,1])/sigma_adj1)-pnorm((-c-Z[,1])/sigma_adj1))
        lmoment2 <- (pnorm((c-Z[,2])/sigma_adj2)-pnorm((-c-Z[,2])/sigma_adj2))
        lmoments <- lmoment1*(1-lmoment2)-lmoment2*(1-lmoment1)
    } else if (length(cutoffs)==2){
        c1 <- cutoffs(1)
        lmoment1 <- (pnorm((c1-Z[,1])/sigma_adj1)-pnorm((-c1-Z[,1])/sigma_adj1))
        lmoment2 <- (pnorm((c1-Z[,2])/sigma_adj2)-pnorm((-c1-Z[,2])/sigma_adj2))

        c2 <- cutoffs(2)
        lmoment3 <- (pnorm((c2-Z[,1])/sigma_adj1)-pnorm((-c2-Z[,1])/sigma_adj1))
        lmoment4 <- (pnorm((c2-Z[,2])/sigma_adj2)-pnorm((-c2-Z[,2])/sigma_adj2))

        lmoments[,1] <- lmoment1*(1-lmoment2)-lmoment2*(1-lmoment1)
        lmoments[,2] <- lmoment3*(1-lmoment4)-lmoment4*(1-lmoment3)
        lmoments[,3] <- lmoment1*(1-lmoment4)-lmoment2*(1-lmoment3)
        lmoments[,4] <- lmoment3*(1-lmoment2)-lmoment4*(1-lmoment1)
    }
}

# Moments
moments <- lmoments/repmat(phat,1,ncol(as.matrix(lmoments)))

return(moments)

}

