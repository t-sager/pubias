
pacman::p_load(devtools,
tidyverse,
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

library(pubias)

#Import data
#Deworming --> Meta example
data <- read.csv(here("data","deworming","cleaned_deworming_data.csv"),header = FALSE)
data <- as.matrix(data)

studynames <- read.csv(here("data","deworming","sorted_names.csv"))
studynames <- as.character(studynames[,1])

pubias_meta(data, studynames, symmetric = 1, symmetric_p = 0, GMM = FALSE)


#EconExperiments --> Replication example
data <- read.csv(here("data","EconExperiments","cleaned_econ_data.csv"),header = FALSE)
data <- as.matrix(data)

studynames <- read.csv(here("data","EconExperiments","sorted_names.csv"))
studynames <- as.character(studynames[,1])

pubias_replication(data, studynames, symmetric = 1, cutoffs = 1.96, GMM = FALSE)

##############################
# Get estimates
# Estimates <- gmm_meta(X,sigma,symmetric,cluster_ID,cutoffs,Studynames)

# Estimates <- gmm_replication(Z,sigmaZ2,symmetric,cluster_ID,cutoffs,Studynames)

# Estimates <- mle_meta(X,sigma,symmetric,cluster_ID,cutoffs,Studynames)
#
# Estimates <- mle_replication(Z,sigmaZ2,symmetric,cluster_ID,cutoffs,Studynames)

# Result
Psihat <- Estimates$Psihat
Varhat <- Estimates$Varhat

# correcting estimates
bias_correction(X,sigma,Psihat,Varhat,cutoffs,symmetric,symmetric_p,identificationapproach)

# Producing bias-corrected estimates and confidence sets
plot_correction(X,sigma,Psihat,Varhat,cutoffs,symmetric,symmetric_p,Studynames,identificationapproach)

#Producing figures
if (identificationapproach==1){
descriptive_stats(Z,sigmaZ2, 1,name,symmetric,cluster_ID)
}else{
  descriptive_stats(X, sigma, identificationapproach,name,symmetric,cluster_ID);
}



# GMM Meta --> works
# GMM Replication --> works
# MLE Meta --> works
# MLE Replication --> Still Problem with numerical integration! --> take it out! only differ between symmetric==1/0!
# Grafischer Output (Horizontal Bars) überarbeiten!
####################
