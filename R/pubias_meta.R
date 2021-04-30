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
pubias_meta <- function(data, studynames, symmetric = 1, symmetric_p = 1, cutoffs = 1.96, GMM = FALSE, print_plots = FALSE){

  # Preliminaries
  X <<- as.matrix(data[,1])
  sigma <- as.matrix(data[,2])
  cluster_ID <<- as.matrix(data[,3])
  n <- nrow(X)
  C <- matrix(1,n,1)
  identificationapproach <- 2

  if (GMM == TRUE) {
    name <- 'GMM_Meta'

    if (print_plots == FALSE) {
    result <- gmm_meta(X, sigma, symmetric, cluster_ID, cutoffs, studynames)
    corrected_estimates <- bias_correction(X,Z,sigma,result,cutoffs,symmetric,symmetric_p,identificationapproach,GMM)
    pubias_result <<- list("GMM Meta Results" = result, "Corrected Estimates" = corrected_estimates)
    rm(corrected_estimates, result)
    } else {
    result <- gmm_meta(X, sigma, symmetric, cluster_ID, cutoffs, studynames)
    corrected_estimates <- bias_correction(X,Z,sigma,result,cutoffs,symmetric,symmetric_p,identificationapproach,GMM)
    descriptives <- descriptive_stats(X, sigma, identificationapproach, name, symmetric, cluster_ID)
    plots <- plot_correction(X,sigma,Psihat,Varhat,cutoffs,symmetric,symmetric_p,studynames,identificationapproach, corrected_estimates)
    rmarkdown::render("R/dashboard.Rmd", params = list(plots=plots, descriptives = descriptives, pub_prob = result$Psihat),output_file = paste0(rprojroot::find_rstudio_root_file(), "/dashboard.html"))
    pubias_result <<- list("GMM Meta Results" = result, "Corrected Estimates" = corrected_estimates, "Descriptive Plots" = descriptives, "Correction Plots" = plots)
    }

    } else {
    name <- 'MLE_Meta'

    if (print_plots == FALSE) {
    result <- mle_meta(X, sigma, symmetric, symmetric_p, cluster_ID, cutoffs, studynames, C)
    corrected_estimates <- bias_correction(X,Z,sigma,result,cutoffs,symmetric,symmetric_p,identificationapproach,GMM)
    pubias_result <<- list("MLE Meta Results" = result, "Corrected Estimates" = corrected_estimates)
    rm(corrected_estimates, result)
    } else {
    result <- mle_meta(X, sigma, symmetric, symmetric_p, cluster_ID, cutoffs, studynames, C)
    corrected_estimates <- bias_correction(X,Z,sigma,result,cutoffs,symmetric,symmetric_p,identificationapproach,GMM)
    descriptives <- descriptive_stats(X, sigma, identificationapproach, name, symmetric, cluster_ID)
    plots <- plot_correction(X,sigma,Psihat,Varhat,cutoffs,symmetric,symmetric_p,studynames,identificationapproach, corrected_estimates)
    rmarkdown::render("R/dashboard.Rmd", params = list(plots=plots, descriptives = descriptives, pub_prob = result$Psihat[-c(1,2)]),output_file = paste0(rprojroot::find_rstudio_root_file(), "/dashboard.html"))
    pubias_result <<- list("MLE Meta Results" = result, "Corrected Estimates" = corrected_estimates, "Descriptive Plots" = descriptives, "Correction Plots" = plots)
    rm(corrected_estimates, result, plots, descriptives)
    }
    }

  rm(cluster_ID, X, envir = globalenv())

}
