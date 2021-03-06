% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mle_meta.R
\name{mle_meta}
\alias{mle_meta}
\title{Computing the publication probability in meta analyses}
\usage{
mle_meta(X, sigma, symmetric, cluster_ID, cutoffs, C)
}
\arguments{
\item{X}{A \verb{n x 1} matrix containing the estimates, where \code{n} is the number of estimates.}

\item{sigma}{A \verb{n x 1} matrix containing the standard errors of the estimates, where \code{n} is the number of estimates.}

\item{symmetric}{If set to TRUE, the publication probability is assumed to be symmetric around zero. If set to FALSE, asymmetry is allowed.}

\item{cluster_ID}{A \verb{n x 1} matrix containing IDs going from 1 to \code{n}, where \code{n} is the number of estimates.}

\item{cutoffs}{A matrix containing the thresholds for the steps of the publication probability. Should be strictly increasing column
vector of size \verb{k x 1} where \code{k} is the number of cutoffs.}

\item{C}{A \verb{n x 1} matrix with all values being 1. Controls.}
}
\value{
Returns a list object with the publication probability (\code{Psihat}), its variance (\code{Varhat}) and robust standard error (\code{se_robust}).
}
\description{
\code{mle_meta()}  calculates the publication probability, its variance and robust standard errors
of replication studies by a Maximum Likelihood approach.  Check package vignette for further information.
}
