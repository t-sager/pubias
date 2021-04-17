# Funktionen laden
source("R/Clustered_covariance_estimate.R")
source("R/gmm_meta.R")
source("R/gmm_replication.R")
source("R/DescriptiveStats.R")
source("R/DescriptiveStatsCombined.R")
source("R/HorizontalBars.R")
source("R/MetaRegressionTable.R")
source("R/MetastudyGMMObjective.R")
source("R/MetastudyMoments.R")
source("R/mle_meta.R")
source("R/mle_replication.R")
source("R/ReplicationAnalyticLogLikelihoodControls.R")
source("R/ReplicationGMMObjective.R")
source("R/ReplicationMoments.R")
source("R/ReplicationpzthetaLogLikelihood.R")
source("R/RobustVariance.R")
source("R/SelectionTable.R")
source("R/SelectionTableGMM.R")
source("R/Step_function_normal_cdf.R")
source("R/VariationVarianceLogLikelihoodControls.R")
source("R/VariationVarianceLogLikelihoodtdist.R")
source("R/VariationVarianceLogLikelihoodtdist_logit.R")

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
#Deworming --> Meta
data <- read.csv(here("data","deworming","cleaned_deworming_data.csv"),header = FALSE)
data <- as.matrix(data)

Studynames <- read.csv(here("data","deworming","sorted_names.csv"))
Studynames <- as.character(Studynames[,1])

#EconExperiments --> Replication
# data <- read.csv(here("data","EconExperiments","cleaned_econ_data.csv"),header = FALSE)
# data <- as.matrix(data)
#
# Studynames <- read.csv(here("data","EconExperiments","sorted_names.csv"))
# Studynames <- as.character(Studynames[,1])

Z <- cbind(data[,1]/data[,2], data[,3]/data[,2])
sigmaZ2 <- data[,4]/data[,2]
X <- as.matrix(data[,1])
sigma <- as.matrix(data[,2])
cluster_ID <- as.matrix(data[,3])

symmetric <- 0
symmetric_p <- 1
cutoffs <- 1.96
#cutoffs <- c(-1.96, 0, 1.96)
vcutoff <- 1.96
time_trend <- 0
spec_test <- 0 # run baseline model not spec test
numerical_integration <- 0
asymmetric_likelihood_spec <- 1
momentchoice <- 1
n <- length(X)
C <- matrix(1,n,1)
name <- 'MLEResults'

# Psihat0 <- c(0,1,length(cutoffs)); #betap
#Psihat0<-c(0,1, rep(1,length(cutoffs))) # Länge ändert sich je nach Appliaktion!
#Psihat0 <- c(0,1,1,1) #betap

includeinfigure <<- as.logical(rep(1,n))
includeinestimation <- as.logical(rep(1,n))



#If you want to get robust confidence sets
# dorobust <- 1




##############################
# Get estimates
# identificationapproach <- 2
# Estimates <- gmm_meta(X,sigma,symmetric,cluster_ID,cutoffs,Studynames)

# identificationapproach <- 1
# Estimates <- gmm_replication(Z,sigmaZ2,symmetric,cluster_ID,cutoffs,Studynames)

identificationapproach <- 2
Estimates <- mle_meta(X,sigma,symmetric,cluster_ID,cutoffs,Studynames)

# identificationapproach <- 1
# Estimates <- mle_replication(Z,sigmaZ2,symmetric,cluster_ID,cutoffs,Studynames)

# Result
Psihat <- Estimates$Psihat
Varhat <- Estimates$Varhat

# Producing bias-corrected estimates and confidence sets
HorizontalBars(X,sigma,Psihat,Varhat,cutoffs,symmetric,symmetric_p,Studynames,identificationapproach)

#Producing figures
if (dofigures==1){
if (identificationapproach==1){
DescriptiveStats(Z,sigmaZ2, 1,name,symmetric,cluster_ID)
DescriptiveStats(data[,1],data[,2], 2,paste0(name, "SanityCheck"),symmetric,cluster_ID)
DescriptiveStatsCombined(Z,data[,1],data[,2],name,symmetric)
}else{
  DescriptiveStats(X, sigma, identificationapproach,name,symmetric,cluster_ID);
}
}


# GMM Meta --> works
# GMM Replication --> works
# MLE Meta
# MLE Replication -
# Grafischer Output (Horizontal Bars) überarbeiten!
####################
