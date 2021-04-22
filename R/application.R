# Funktionen laden
source("R/clustered_covariance_estimate.R")
source("R/gmm_meta.R")
source("R/gmm_replication.R")
source("R/descriptive_stats.R")
source("R/bias_correction.R")
source("R/plot_correction.R")
source("R/MetastudyGMMObjective.R")
source("R/MetastudyMoments.R")
source("R/mle_meta.R")
source("R/mle_replication.R")
source("R/replication_analytic_llh.R")
source("R/ReplicationGMMObjective.R")
source("R/ReplicationMoments.R")
source("R/robust_variance.R")
source("R/step_function_normal_cdf.R")
source("R/variation_variance_llh.R")

# Packages laden
pacman::p_load(tidyverse,
               RColorBrewer,
               latex2exp,
               xtable,
               here,
               invgamma,
               matlab,
               PEIP,
               mvtnorm,
               ggplot2,
               gridExtra,
               cowplot)




#Import data
#Deworming --> Meta example
# data <- read.csv(here("data","deworming","cleaned_deworming_data.csv"),header = FALSE)
# data <- as.matrix(data)
#
# Studynames <- read.csv(here("data","deworming","sorted_names.csv"))
# Studynames <- as.character(Studynames[,1])

#EconExperiments --> Replication example
data <- read.csv(here("data","EconExperiments","cleaned_econ_data.csv"),header = FALSE)
data <- as.matrix(data)

Studynames <- read.csv(here("data","EconExperiments","sorted_names.csv"))
Studynames <- as.character(Studynames[,1])

Z <- cbind(data[,1]/data[,2], data[,3]/data[,2])
sigmaZ2 <- data[,4]/data[,2]
X <- as.matrix(data[,1])
sigma <- as.matrix(data[,2])
cluster_ID <- as.matrix(data[,3])
n <- length(X)
C <- matrix(1,n,1)
name <- 'MLEResults'

symmetric <- 1
symmetric_p <- 0 # nur für mle_meta relevant!
cutoffs <- 1.96
# cutoffs <- c(-1.96, 0, 1.96)
vcutoff <- 1.96

asymmetric_likelihood_spec <- 1

includeinfigure <<- as.logical(rep(1,n))
includeinestimation <- as.logical(rep(1,n))

##############################
# Get estimates
# identificationapproach <- 2
# Estimates <- gmm_meta(X,sigma,symmetric,cluster_ID,cutoffs,Studynames)

identificationapproach <- 1
Estimates <- gmm_replication(Z,sigmaZ2,symmetric,cluster_ID,cutoffs,Studynames)

# identificationapproach <- 2
# Estimates <- mle_meta(X,sigma,symmetric,cluster_ID,cutoffs,Studynames)

# identificationapproach <- 1
# Estimates <- mle_replication(Z,sigmaZ2,symmetric,cluster_ID,cutoffs,Studynames)

# Result
Psihat <- Estimates$Psihat
Varhat <- Estimates$Varhat

# Producing bias-corrected estimates and confidence sets
HorizontalBars(X,sigma,Psihat,Varhat,cutoffs,symmetric,symmetric_p,Studynames,identificationapproach)

#Producing figures
if (identificationapproach==1){
DescriptiveStats(Z,sigmaZ2, 1,name,symmetric,cluster_ID)
DescriptiveStats(data[,1],data[,2], 2,paste0(name, "SanityCheck"),symmetric,cluster_ID)
DescriptiveStatsCombined(Z,data[,1],data[,2],name,symmetric)
}else{
  DescriptiveStats(X, sigma, identificationapproach,name,symmetric,cluster_ID);
}



# GMM Meta --> works
# GMM Replication --> works
# MLE Meta --> works
# MLE Replication --> Still Problem with numerical integration! --> take it out! only differ between symmetric==1/0!
# Grafischer Output (Horizontal Bars) überarbeiten!
####################
