#' Identification of and correction for publication bias in meta studies
#'
#'`pubias_meta()` calculates the publication probability either by a Maximum Likelihood or a GMM approach and
#'then reports the corrected estimates for the estimates of a meta-analysis.
#'
#' @param data A `n x 3` matrix where the first column contains the treatment effects, the second column the associated
#' standard errors of the estimates and in the third column an ID going from 1 to `n`, where `n` is the number of estimates.
#' @param studynames A vector of type `character` containing all the Studynames of size `n` in the same order as the argument `data`.
#' @param cutoffs A matrix containing the thresholds for the steps of the publication probability. Should be strictly increasing column
#' vector of size `k x 1` where `k` is the number of cutoffs. By default, the cutoff is set at 1.96.
#' @param GMM If set to TRUE, the publication probability will be estimated via GMM. By default, it is set to FALSE which uses the MLE
#' method for estimation.
#' @param symmetric If set to `1`, the publication probability is assumed to be symmetric around zero. If set to `0`, asymmetry is allowed.
#' @param symmetric_p If set to `1`, the publication probability is assumed to be symmetric around zero and all cutoffs should be positive.
#' If set to `0`, asymmetry is allowed and cutoffs should be specified in increasing order.
#' @param print_plots If set to TRUE, descriptive plots as well as correction plots are printed into the working directory in .pdf format.
#' @param print_dashboard If set to TRUE, additionally to the .pdf plots, a dashboard with the same charts in dynamic format is produced.
#' The dashboard is saved in the working directory. Only possbile if `print_plots` is set to TRUE.
#'
#' @return Returns a list object called `pubias_result` with the elements `Results` and `Corrected Estimates`.
#' If specified, the list also includes the descriptive, as well as the correction plots as `ggplot` objects.
#'
#' - `Results` contains the publication probability `Psihat`, the variance as well as the robust standard error.
#'
#' - `Corrected Estimates` contains the original estimates with their 95% confidence bonds (`Z1`, `Z1_L`, `Z1_U`)
#' as well as the corrected estimates (`Z1_M`) and in addition the Bonferroni corrected 95% confidence bonds (`Z1_LB`, `Z1_UB`).
#' There are additional elements which are mainly used for plotting the results.
#'
#' @export
#'
#' @examples
pubias_meta <-
  function(data,
           studynames,
           symmetric = 1,
           symmetric_p = 1,
           cutoffs = 1.96,
           GMM = FALSE,
           print_plots = FALSE,
           print_dashboard = FALSE) {
    # Preliminaries
    X <<- as.matrix(data[, 1])
    sigma <- as.matrix(data[, 2])
    cluster_ID <<- as.matrix(data[, 3])
    n <- nrow(X)
    C <- matrix(1, n, 1)
    identificationapproach <- 2
    wd <- getwd()

    if (GMM == TRUE) {
      name <- 'GMM_Meta'

      if (print_plots == FALSE) {
        result <-
          gmm_meta(X, sigma, symmetric, cluster_ID, cutoffs, studynames)
        corrected_estimates <-
          bias_correction(
            X,
            Z,
            sigma,
            result,
            cutoffs,
            symmetric,
            symmetric_p,
            identificationapproach,
            GMM
          )
        pubias_result <<-
          list("GMM Meta Results" = result,
               "Corrected Estimates" = corrected_estimates)
        rm(corrected_estimates, result)
      } else {
        result <-
          gmm_meta(X, sigma, symmetric, cluster_ID, cutoffs, studynames)
        corrected_estimates <-
          bias_correction(
            X,
            Z,
            sigma,
            result,
            cutoffs,
            symmetric,
            symmetric_p,
            identificationapproach,
            GMM
          )
        descriptives <-
          descriptive_stats(X,
                            sigma,
                            identificationapproach,
                            name,
                            symmetric,
                            cluster_ID)
        plots <-
          plot_correction(
            X,
            sigma,
            Psihat,
            Varhat,
            cutoffs,
            symmetric,
            symmetric_p,
            studynames,
            identificationapproach,
            corrected_estimates
          )
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
            "Corrected Estimates" = corrected_estimates,
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
                   symmetric_p,
                   cluster_ID,
                   cutoffs,
                   studynames,
                   C)
        corrected_estimates <-
          bias_correction(
            X,
            Z,
            sigma,
            result,
            cutoffs,
            symmetric,
            symmetric_p,
            identificationapproach,
            GMM
          )
        pubias_result <<-
          list("MLE Meta Results" = result,
               "Corrected Estimates" = corrected_estimates)
        rm(corrected_estimates, result)
      } else {
        result <-
          mle_meta(X,
                   sigma,
                   symmetric,
                   symmetric_p,
                   cluster_ID,
                   cutoffs,
                   studynames,
                   C)
        corrected_estimates <-
          bias_correction(
            X,
            Z,
            sigma,
            result,
            cutoffs,
            symmetric,
            symmetric_p,
            identificationapproach,
            GMM
          )
        descriptives <-
          descriptive_stats(X,
                            sigma,
                            identificationapproach,
                            name,
                            symmetric,
                            cluster_ID)
        plots <-
          plot_correction(
            X,
            sigma,
            Psihat,
            Varhat,
            cutoffs,
            symmetric,
            symmetric_p,
            studynames,
            identificationapproach,
            corrected_estimates
          )
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
            "MLE Meta Results" = result,
            "Corrected Estimates" = corrected_estimates,
            "Descriptive Plots" = descriptives,
            "Correction Plots" = plots
          )
        rm(corrected_estimates, result, plots, descriptives)
      }
    }

    if (print_plots == FALSE && print_dashboard == TRUE) {
      stop("print_plots has to be set to TRUE if the dashboard wants to be printed!")
    }
    rm(cluster_ID, X, envir = globalenv())

  }
