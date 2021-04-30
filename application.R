##############################
devtools::load_all()
library(pubias)

##############################
# Example: Metastudies - Deworming
#Import data
data <- read.csv(here("data","deworming","cleaned_deworming_data.csv"),header = FALSE)
data <- as.matrix(data)

studynames <- read.csv(here("data","deworming","sorted_names.csv"))
studynames <- as.character(studynames[,1])

pubias_meta(data, studynames, symmetric = 1, symmetric_p = 1, cutoffs = 1.96, GMM = FALSE, print_plots = TRUE)

##############################
# Example: Replication Studies - EconExperiments
data <- read.csv(here("data","EconExperiments","cleaned_econ_data.csv"),header = FALSE)
data <- as.matrix(data)

studynames <- read.csv(here("data","EconExperiments","sorted_names.csv"))
studynames <- as.character(studynames[,1])

pubias_replication(data, studynames, symmetric = 1, cutoffs = 1.96, GMM = FALSE, print_plots = TRUE)

