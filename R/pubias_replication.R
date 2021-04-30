#' Title
#'
#' @param data
#' @param studynames
#' @param symmertic
#' @param cutoffs
#' @param GMM
#'
#' @return
#' @export
#'
#' @examples
pubias_replication <- function(data, studynames, symmetric = 1, cutoffs = 1.96, GMM = FALSE, print_plots = FALSE){

  # Preliminaries
  Z <<- cbind(data[,1]/data[,2], data[,3]/data[,2])
  sigmaZ2 <- as.matrix(data[,4]/data[,2])
  X <<- as.matrix(data[,1])
  sigma <- as.matrix(data[,2])
  cluster_ID <<- as.matrix(data[,3])
  n <- nrow(X)
  C <- matrix(1,n,1)
  identificationapproach <- 1

  if (GMM == TRUE) {
    name <- 'GMM_Replication'

    if (print_plots == FALSE) {
    result <- gmm_replication(Z, sigmaZ2, symmetric, cluster_ID, cutoffs, studynames)
    corrected_estimates <- bias_correction(X,Z,sigma,result,cutoffs,symmetric,symmetric_p=0,identificationapproach,GMM)
    pubias_result <<- list("GMM Replication Results" = result, "Corrected Estimates" = corrected_estimates)
    } else {
    result <- gmm_replication(Z, sigmaZ2, symmetric, cluster_ID, cutoffs, studynames)
    corrected_estimates <- bias_correction(X,Z,sigma,result,cutoffs,symmetric,symmetric_p=0,identificationapproach,GMM)
    descriptives <- descriptive_stats(Z, sigmaZ2, identificationapproach, name, symmetric, cluster_ID)
    plots <- plot_correction(Z,sigmaZ2,Psihat,Varhat,cutoffs,symmetric,symmetric_p,studynames,identificationapproach, corrected_estimates)
    rmarkdown::render("R/dashboard.Rmd", params = list(plots=plots, descriptives = descriptives, pub_prob = result$Psihat),output_file = paste0(rprojroot::find_rstudio_root_file(), "/dashboard.html"))
    pubias_result <<- list("GMM Replication Results" = result, "Corrected Estimates" = corrected_estimates, "Descriptive Plots" = descriptives, "Correction Plots" = plots)
    }

    } else {
    name <- 'MLE_Replication'

    if (print_plots == FALSE) {
    result <- mle_replication(Z, sigmaZ2, symmetric, cluster_ID, cutoffs, studynames, C)
    corrected_estimates <- bias_correction(X,Z,sigma,result,cutoffs,symmetric,symmetric_p=0,identificationapproach,GMM)
    pubias_result <<- list("MLE Replication Results" = result, "Corrected Estimates" = corrected_estimates)
    } else {
    result <- mle_replication(Z, sigmaZ2, symmetric, cluster_ID, cutoffs, studynames, C)
    corrected_estimates <- bias_correction(X,Z,sigma,result,cutoffs,symmetric,symmetric_p=0,identificationapproach,GMM)
    descriptives <- descriptive_stats(Z, sigmaZ2, identificationapproach, name, symmetric, cluster_ID)
    plots <- plot_correction(Z,sigmaZ2,Psihat,Varhat,cutoffs,symmetric,symmetric_p,studynames,identificationapproach, corrected_estimates)
    rmarkdown::render("R/dashboard.Rmd", params = list(plots=plots, descriptives = descriptives, pub_prob = result$Psihat[-c(1,2)]),output_file = paste0(rprojroot::find_rstudio_root_file(), "/dashboard.html"))
    pubias_result <<- list("MLE Replication Results" = result, "Corrected Estimates" = corrected_estimates, "Descriptive Plots" = descriptives, "Correction Plots" = plots)
    }
    }
  rm(cluster_ID, X, Z, envir = globalenv())
}
