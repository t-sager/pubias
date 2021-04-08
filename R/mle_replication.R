# Only old Code from 2017 approach ---> Rewrite to 2019 approach! (See Matlab Code)


# not tested
source("ReplicationAnalyticLogLikelihood.R")
source("ReplicationpzthetaLogLikelihood.R")
source("VariationVarianceLogLikelihood.R")
source("SelectionTable.R")
source("RobustVariance.R")
#show(number)
MaxEval<-10^5
MaxIter<-10^5;
Tol<-10^(-8);
stepsize<-10^(-6);

if (identificationapproach==1) {
  if (symmetric == 0) {
    LLH <-function(Psi) {
      A<-ReplicationAnalyticLogLikelihood(Psi[1], Psi[2], c(Psi[-c(1,2)], 1),cutoffs,symmetric, Z, sigmaZ2);
      #A<-A$LLH
      return(A)
    }
  } else {
    if (spec_test == 1) {
      LLH <-function(Psi) {

        A<-ReplicationpzthetaLogLikelihood(0,Psi[1], c(Psi[2:(2+length(cutoffs)-1)], 1),
                                           c(rep (0, length(Psi[2:2+length(cutoffs)-1])),0),
                                           cutoffs,symmetric, Z, sigmaZ2,vcutoff, Zdraws, Thetadraws,zdrawsd,thetadrawsd);
        #A<-A$LLH
        return(A)
      }

      LLH_test <-function(Psi) {
        A<-ReplicationpzthetaLogLikelihood(0,Psi[1],  c(Psi[2:(2+length(cutoffs)-1)], 1),
                                           c(Psi[(2+length(cutoffs)):length(Psi)], 0),
                                           cutoffs,symmetric, Z, sigmaZ2,
                                           vcutoff, Zdraws, Thetadraws,zdrawsd,thetadrawsd);
        #A<-A$LLH
        return(A)
      }
    } else {

      LLH <-function (Psi) {
        A<-ReplicationAnalyticLogLikelihood(0,Psi[1], c(Psi[-1], 1),cutoffs,symmetric, Z, sigmaZ2);
        #A<-A$LLH
        return(A)
      }
    }
  }
  nn<<- number;



if (identificationapproach == 1) {

  #date is formatted as [original study X1, sigma1, replication study X2, sigma2]
  filepath<-"/Applications/EconExperiments/cleaned_econ_data.csv";
  my.data<-read.csv(paste0(pathname,filepath),header=FALSE)

  # Make sure Studynames is an nx1 cell array where n is the number of studies
  filepath<-"/Applications/EconExperiments/sorted_names.csv";
  Studynames<-read.csv(paste0(pathname,filepath),header=FALSE)
  Studynames<<-Studynames$V1

  #Name for outputs
  name="ReplicationEcon";

  #Assign these variables based on your dataset and parametric specification but do not change the variable names
  X=my.data[,1];
  sigma<-my.data[,2];
  Z<-cbind(my.data[,1]/my.data[,2],my.data[,3]/my.data[,2])
  sigmaZ2<-my.data[,4]/my.data[,2];

  # Use a step function symmetric around zero
  symmetric=1;
  # Set cutoffs to use in step function: should be given in increasing order;
  cutoffs=1.96;
  # starting values for optimization: thetabar is fixed at zero due to sign normalization
  # [tau, betap(1)]
  Psihat0=c(2,1);

  number=dim(Z)[1]; #number of studies
  cluster_ID=1:number; #for cluster-robust variance estimation
  #for replication studies, there is no clustering i.e. each study is its own cluster

}
