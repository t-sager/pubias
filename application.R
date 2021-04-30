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

setwd("C:/Users/tills/Desktop/output")

pubias_meta(data, studynames, GMM = TRUE, print_plots = FALSE, print_dashboard = TRUE)

setwd("~/UniBe/02 Master/03 FS21/Workshop Econometrics/R/pubias")
##############################
# Example: Replication Studies - EconExperiments
data <- read.csv(here("data","EconExperiments","cleaned_econ_data.csv"),header = FALSE)
data <- as.matrix(data)

studynames <- read.csv(here("data","EconExperiments","sorted_names.csv"))
studynames <- as.character(studynames[,1])

setwd("C:/Users/tills/Desktop/output")

pubias_replication(data, studynames, GMM = TRUE, print_plots = TRUE)
setwd("~/UniBe/02 Master/03 FS21/Workshop Econometrics/R/pubias")
