#' Identification of and correction for publication bias in meta studies
#'
#'`pubias_meta()` calculates the publication probability either by a Maximum Likelihood or a GMM approach and
#'then reports the corrected estimates for the estimates of a meta-analysis.  Check package vignette for further information.
#'
#' @param data A `n x 3` matrix where the first column contains the treatment effects, the second column the associated
#' standard errors of the estimates and in the third column an ID going from 1 to `k`, where `k` is the number of unique studies.
#' @param studynames Optional. A vector of type `character` containing all the Studynames of size `n` in the same order as the argument `data`.
#' @param sign_lvl A value indicating the significance level at which the analysis should be done. Ultimately leads to the threshold (z-score) for the steps of the publication probability.
#' By default, the significance level is set at 5%, hence `0.05`.
#' @param GMM If set to `TRUE`, the publication probability will be estimated via GMM. By default, it is set to `FALSE` which uses the MLE
#' method for estimation.
#' @param symmetric If set to `TRUE`, the publication probability is assumed to be symmetric around zero. If set to `FALSE`, asymmetry is allowed.
#' @param print_plots If set to `TRUE`, descriptive plots as well as correction plots are printed into the working directory in .pdf format.
#' @param print_dashboard If set to `TRUE`, additionally to the .pdf plots, a dashboard with the same charts in dynamic format is produced.
#' The dashboard is saved in the working directory. Only possible if `print_plots` is set to `TRUE`.
#'
#' @return Returns a list object called `pubias_result` with the elements `Results` and `Corrected Estimates`.
#' If specified, the list also includes the descriptive, as well as the correction plots as `ggplot` objects.
#'
#' `Results` contains the publication probability (`Psihat`), its variance (`Varhat`) and robust standard error (`se_robust`).
#'
#' `Corrected Estimates` contains the original estimates and the median unbiased estimates (`original`, `adj_estimates`) as well as the corrected 95% confidence bounds (`adj_L`, `adj_U`). In addition, the Bonferroni corrected 95% confidence bounds are reported (`adj_LB`, `adj_UB`).
#'
#' @export
#'
pubias_meta <-
  function(data,
           studynames = NULL,
           symmetric = TRUE,
           sign_lvl = 0.05,
           GMM = FALSE,
           print_plots = FALSE,
           print_dashboard = FALSE) {

  #- Preliminaries

    ## No Studynames available
    if (is.null(studynames)) {
      studynames <- as.character(data[, 3])
    }

    ## Setting up the data as needed
    cutoffs <<- qnorm(sign_lvl/2, lower.tail = FALSE)
    X <<- as.matrix(data[, 1])
    sigma <- as.matrix(data[, 2])
    cluster_ID <<- as.matrix(data[, 3])
    n <- nrow(X)
    C <- matrix(1, n, 1)
    identificationapproach <- 2
    wd <- getwd()

  #- Running Identification and Correction by calling other functions
    if (GMM == TRUE) {
      name <- 'GMM_Meta'

      if (print_plots == FALSE) {
        result <-
          gmm_meta(X, sigma, symmetric, cluster_ID, cutoffs)
        corrected_estimates <-
          bias_correction(
            X,
            sigma,
            result,
            cutoffs,
            symmetric,
            identificationapproach,
            GMM
          )
        pubias_result <<-
          list("GMM Meta Results" = result,
               "Corrected Estimates" = corrected_estimates[c("original","adj_estimates","adj_L", "adj_U", "adj_LB", "adj_UB")])
        rm(corrected_estimates, result)
      } else {
        result <-
          gmm_meta(X, sigma, symmetric, cluster_ID, cutoffs)
        corrected_estimates <-
          bias_correction(
            X,
            sigma,
            result,
            cutoffs,
            symmetric,
            identificationapproach,
            GMM
          )
        descriptives <-
          descriptive_stats(X,
                            sigma,
                            identificationapproach,
                            name,
                            symmetric,
                            cluster_ID,
                            cutoffs)

        plots <-
          plot_correction(X, corrected_estimates,cutoffs,symmetric,studynames,identificationapproach)
        if (print_dashboard == TRUE) {
          rmarkdown::render(
            system.file("dashboard.Rmd", package = "pubias"),
            params = list(
              plots = plots,
              descriptives = descriptives,
              pub_prob = result$Psihat
            ),
            output_file = paste0(wd, "/", name, "_Dashboard.html")
          )
        }
        pubias_result <<-
          list(
            "GMM Meta Results" = result,
            "Corrected Estimates" = corrected_estimates[c("original", "adj_estimates","adj_L", "adj_U", "adj_LB", "adj_UB")],
            "Descriptive Plots" = descriptives,
            "Correction Plots" = plots
          )
      }

    } else {
      name <- 'MLE_Meta'

      if (print_plots == FALSE) {
        result <-
          mle_meta(X,
                   sigma,
                   symmetric,
                   cluster_ID,
                   cutoffs,
                   C)
        corrected_estimates <-
          bias_correction(
            X,
            sigma,
            result,
            cutoffs,
            symmetric,
            identificationapproach,
            GMM
          )
        result <- list("Psihat" = result$Psihat[-c(1,2)], "Varhat" = result$Varhat[-c(1,2), -c(1,2)], "se_robust" = result$se_robust[-c(1,2)])
        pubias_result <<-
          list("MLE Meta Results" = result,
               "Corrected Estimates" = corrected_estimates[c("original", "adj_estimates","adj_L", "adj_U", "adj_LB", "adj_UB")])
        rm(corrected_estimates, result)
      } else {
        result <-
          mle_meta(X,
                   sigma,
                   symmetric,
                   cluster_ID,
                   cutoffs,
                   C)
        corrected_estimates <-
          bias_correction(
            X,
            sigma,
            result,
            cutoffs,
            symmetric,
            identificationapproach,
            GMM
          )
        descriptives <-
          descriptive_stats(X,
                            sigma,
                            identificationapproach,
                            name,
                            symmetric,
                            cluster_ID,
                            cutoffs)
        plots <-
          plot_correction(X, corrected_estimates,cutoffs,symmetric,studynames,identificationapproach)
        if (print_dashboard == TRUE) {
          rmarkdown::render(
            system.file("dashboard.Rmd", package = "pubias"),
            params = list(
              plots = plots,
              descriptives = descriptives,
              pub_prob = result$Psihat[-c(1,2)]
            ),
            output_file = paste0(wd, "/", name, "_Dashboard.html")
          )
        }
        result <- list("Psihat" = result$Psihat[-c(1,2)], "Varhat" = result$Varhat[-c(1,2), -c(1,2)], "se_robust" = result$se_robust[-c(1,2)])
        pubias_result <<-
          list(
            "MLE Meta Results" = result,
            "Corrected Estimates" = corrected_estimates[c("original", "adj_estimates","adj_L", "adj_U", "adj_LB", "adj_UB")],
            "Descriptive Plots" = descriptives,
            "Correction Plots" = plots
          )
        rm(corrected_estimates, result, plots, descriptives)
      }
    }

    # Throw error if dashboard wants to be printed but plots is FALSE
    if (print_plots == FALSE && print_dashboard == TRUE) {
      stop("print_plots has to be set to TRUE if the dashboard wants to be printed!")
    }
    rm(cluster_ID, X, cutoffs, envir = globalenv())

  }
