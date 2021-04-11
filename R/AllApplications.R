#this is the masterfile to do all our empirical applications
clear
rng(1)
#create a folder "FiguresandTables" under this path to save figures and tables

#picking which application to run
# 1: econ experiments
# 2: econ experiments, p(z,omega) specification test
# 3: econ experiments, controlling for journal
# 4: econ experiments, GMM
# 5: econ experiments, GMM with alternate moments
# 6: psychology experiments
# 7: psychology experiments, approved replications
# 8: psychology experiments, denominator dof>=30
# 9: psychology experiments, p(z,theta) specification test
# 10: psychology experiments, controlling for journal
# 11: psychology experiments, GMM
# 12: psychology experiments, GMM with alternate moments
# 13: minimum wage effect
# 14: minimum wage effect, published
# 15: minimum wage, no selection
# 16: minimum wage effect, controlling for year
# 17: minimum wage effect, GMM
# 18: deworming
# 19: deworming, asymmetric
# 20: deworming, asymmetric, single discontinuity
# 21: deworming, no selection
# 22: deworming, GMM
# 23: nature/science experiments
# 24: psychology experiments, single cutoff


application <- 1

# whether we want to produce figures, and/or estimates, and/or bias corrections
dofigures <- 1
doestimates <- 1
docorrections <- 1

#simulations draws for p(z,theta) specifications
vcutoff <- 1.96

fixcoeff <- 0#By default, do not hold any coefficients fixed.

##
#loading the data
if (application == 1){

    Studynames <- textscan(fopen('../Applications/EconExperiments/sorted_names.csv','r'), '%s%[^\n\r]', 'Delimiter', ';')
    Studynames <- Studynames{1}(2:19,1)

    data <- csvread('../Applications/EconExperiments/cleaned_econ_data.csv')
    Z <- [data(:,1)./data(:,2) data(:,3)./data(:,2)]
    sigmaZ2 <- data(:,4)./data(:,2)

    n <- size(Z,1)
    cluster_ID <- (1:n)'

    X <- data(:,1)
    sigma <- data(:,2)
    C <- ones(length(X),1)

    identificationapproach <- 1){
    GMMapproach <- 0
    numerical_integration <- 1
    name <- 'ReplicationEcon'

    #Estimate baseline model, rather than running spec test
    spec_test <- 0
    #Set cutoffs to use in step function: should be given in increasing order;
    cutoffs <- 1.96
    #Use a step function symmetric around zero
    symmetric <- 1
    #starting values for optimization
    Psihat0 <- [1,1,1]#[kappa, lambda, betap(1)]

} else if (application == 2){
    Studynames <- textscan(fopen('../Applications/EconExperiments/sorted_names.csv','r'), '%s%[^\n\r]', 'Delimiter', ';')
    Studynames <- Studynames{1}(2:19,1)

    data <- csvread('../Applications/EconExperiments/cleaned_econ_data.csv')
    Z <- [data(:,1)./data(:,2) data(:,3)./data(:,2)]
    sigmaZ2 <- data(:,4)./data(:,2)

    n <- size(Z,1)
    cluster_ID <- (1:n)'
    X <- data(:,1)
    sigma <- data(:,2)
    #normalize mean standard deviation to one
    X <- X/mean(sigma)
    sigma <- sigma/mean(sigma)


    identificationapproach <- 1){
    GMMapproach <- 0
    name <- 'ReplicationEconSpecTest'

    #Run spec test of baseline model which assumes selection p(z)
    spec_test <- 1
    #Set cutoffs to use in step function: should be given in increasing order;
    cutoffs <- 1.96
    #Use a step function symmetric around zero
    symmetric <- 1
    #starting values for optimization
    Psihat0 <- [0.6,1.7,0.04]

} else if (application == 3){
    Studynames <- textscan(fopen('../Applications/EconExperiments/sorted_names.csv','r'), '%s%[^\n\r]', 'Delimiter', ';')
    Studynames <- Studynames{1}(2:19,1)

    data <- csvread('../Applications/EconExperiments/cleaned_econ_data.csv')
    Z <- [data(:,1)./data(:,2) data(:,3)./data(:,2)]
    sigmaZ2 <- data(:,4)./data(:,2)

    n <- size(Z,1)
    cluster_ID <- (1:n)'

    X <- data(:,1)
    sigma <- data(:,2)

    QJE_dummy <- data(:,5)
    C <- [ones(size(Z,1),1) QJE_dummy]

    identificationapproach <- 1){
    GMMapproach <- 0
    numerical_integration <- 1
    name <- 'ReplicationEconControls'

    #Estimate baseline model, rather than running spec test
    spec_test <- 0
    #Set cutoffs to use in step function: should be given in increasing order;
    cutoffs <- 1.96
    #Use a step function symmetric around zero
    symmetric <- 1
    #starting values for optimization
    Psihat0 <- [1,1,1,0]#[kappa, lambda, betap(1), betap(2)]

} else if (application == 4){

    Studynames <- textscan(fopen('../Applications/EconExperiments/sorted_names.csv','r'), '%s%[^\n\r]', 'Delimiter', ';')
    Studynames <- Studynames{1}(2:19,1)

    data <- csvread('../Applications/EconExperiments/cleaned_econ_data.csv')
    Z <- [data(:,1)./data(:,2) data(:,3)./data(:,2)]
    sigmaZ2 <- data(:,4)./data(:,2)

    n <- size(Z,1)
    cluster_ID <- (1:n)'

    X <- data(:,1)
    sigma <- data(:,2)
    C <- ones(length(X),1)

    identificationapproach <- 1){
    GMMapproach <- 1
    momentchoice <- 1
    numerical_integration <- 1
    name <- 'ReplicationEconGMM'

    #Set cutoffs to use in step function: should be given in increasing order;
    cutoffs <- 1.96
    #Use a step function symmetric around zero
    symmetric <- 1
    #starting values for optimization
    Psihat0 <- 0.1#betap(1)

} else if (application == 5){

    Studynames <- textscan(fopen('../Applications/EconExperiments/sorted_names.csv','r'), '%s%[^\n\r]', 'Delimiter', ';')
    Studynames <- Studynames{1}(2:19,1)

    data <- csvread('../Applications/EconExperiments/cleaned_econ_data.csv')
    Z <- [data(:,1)./data(:,2) data(:,3)./data(:,2)]
    sigmaZ2 <- data(:,4)./data(:,2)

    n <- size(Z,1)
    cluster_ID <- (1:n)'

    X <- data(:,1)
    sigma <- data(:,2)
    C <- ones(length(X),1)

    identificationapproach <- 1){
    GMMapproach <- 1
    momentchoice <- 2
    numerical_integration <- 1
    name <- 'ReplicationEconGMMAlt'

    #Set cutoffs to use in step function: should be given in increasing order;
    cutoffs <- 1.96
    #Use a step function symmetric around zero
    symmetric <- 1
    #starting values for optimization
    Psihat0 <- 1#betap(1)

} else if (application == 6){
    #Studynames = textscan(fopen('../Applications/PsychExperiments/sorted_1stAuthor.csv','r'), '#s#[^\n\r]', 'Delimiter', ';');
    #Studynames=Studynames{1}

    data <- csvread('../Applications/PsychExperiments/cleaned_psych_data.csv')
    dofis1=data(:,7)==1
    use <- dofis1
    data <- data(use,:)

    Z <- [data(:,1)./data(:,2) data(:,3)./data(:,2)]
    sigmaZ2 <- data(:,4)./data(:,2)
    n <- size(Z,1)
    cluster_ID <- (1:n)'
    Studynames <- char.empty(0,n)

    X <- data(:,1)
    sigma <- data(:,2)
    C <- ones(length(X),1)

    identificationapproach <- 1){
    GMMapproach <- 0
    numerical_integration <- 1
    name <- 'ReplicationPsych'

    #Estimate baseline model, rather than running spec test
    spec_test <- 0
    #Set cutoffs to use in step function: should be given in increasing order;
    cutoffs <- [1.64 1.96]'
    #Use a step function symmetric around zero
    symmetric <- 1
    #starting values for optimization
    Psihat0 <- [1,1,1,1]#[kappa, lambda, betap(1), betap(2)]

} else if (application == 7){
    data <- csvread('../Applications/PsychExperiments/cleaned_psych_data.csv')
    approved <- logical(data(:,5))
    dofis1=data(:,7)==1
    use <- dofis1&approved
    data <- data(use,:)

    Z <- [data(:,1)./data(:,2) data(:,3)./data(:,2)]
    sigmaZ2 <- data(:,4)./data(:,2)

    n <- size(Z,1)
    Studynames <- char.empty(0,n)
    cluster_ID <- (1:n)'

    X <- data(:,1)
    sigma <- data(:,2)
    C <- ones(length(X),1)

    identificationapproach <- 1){
    GMMapproach <- 0
    numerical_integration <- 1
    name <- 'ReplicationPsychApproved'

    #Estimate baseline model, rather than running spec test
    spec_test <- 0
    #Set cutoffs to use in step function: should be given in increasing order;
    cutoffs <- [1.64 1.96]'
    #Use a step function symmetric around zero
    symmetric <- 1
    #starting values for optimization
    Psihat0 <- [1,1,1,1]#[kappa, lambda, betap(1), betap(2)]

} else if (application == 8){
    data <- csvread('../Applications/PsychExperiments/cleaned_psych_data.csv')
    dofis1=data(:,7)==1
    denominatordof30=data(:,8)==0|data(:,8)>=30
    use <- dofis1&denominatordof30
    data <- data(use,:)

    Z <- [data(:,1)./data(:,2) data(:,3)./data(:,2)]
    sigmaZ2 <- data(:,4)./data(:,2)
    n <- size(Z,1)
    cluster_ID <- (1:n)'
    Studynames <- char.empty(0,n)

    X <- data(:,1)
    sigma <- data(:,2)
    C <- ones(length(X),1)

    identificationapproach <- 1){
    GMMapproach <- 0
    numerical_integration <- 0
    name <- 'ReplicationPsychLargedof'

    #Estimate baseline model, rather than running spec test
    spec_test <- 0
    #Set cutoffs to use in step function: should be given in increasing order;
    cutoffs <- [1.64 1.96]'
    #Use a step function symmetric around zero
    symmetric <- 1
    #starting values for optimization
    Psihat0 <- [.5,1,.01,.2]#[kappa, lambda, betap(1), betap(2)]

} else if (application == 9){
    Studynames <- textscan(fopen('../Applications/PsychExperiments/sorted_1stAuthor.csv','r'), '%s%[^\n\r]', 'Delimiter', ';')
    Studynames <- Studynames{1}
    data <- csvread('../Applications/PsychExperiments/cleaned_psych_data.csv')
    Z <- [data(:,1)./data(:,2) data(:,3)./data(:,2)]
    sigmaZ2 <- data(:,4)./data(:,2)

    n <- size(Z,1)
    cluster_ID <- (1:n)'
    X <- data(:,1)
    sigma <- data(:,2)
    #normalize mean standard deviation to one
    X <- X/mean(sigma)
    sigma <- sigma/mean(sigma)

    identificationapproach <- 1){
    GMMapproach <- 0
    name <- 'ReplicationPsychSpecTest'

    #Run spec test of baseline model which assumes selection p(z)
    spec_test <- 1
    #Set cutoffs to use in step function: should be given in increasing order;
    cutoffs <- [1.64 1.96]'
    #Use a step function symmetric around zero
    symmetric <- 1
    #starting values for optimization
    Psihat0 <- [1,1,1,1]#[kappa, lambda, betap(1), betap(2)]

} else if (application == 10){
    data <- csvread('../Applications/PsychExperiments/cleaned_psych_data.csv')
    dofis1=data(:,7)==1
    use <- dofis1
    data <- data(use,:)

    Z <- [data(:,1)./data(:,2) data(:,3)./data(:,2)]
    sigmaZ2 <- data(:,4)./data(:,2)
    n <- size(Z,1)
    cluster_ID <- (1:n)'
    Studynames <- char.empty(0,n)

    X <- data(:,1)
    sigma <- data(:,2)

    PS_dummy <- data(:,10)
    JPSP_dummy <- data(:,11)
    C <- [ones(size(Z,1),1) PS_dummy JPSP_dummy]

    identificationapproach <- 1){
    GMMapproach <- 0
    numerical_integration <- 1
    docorrections <- 0
    name <- 'ReplicationPsychControls'

    #Estimate baseline model, rather than running spec test
    spec_test <- 0
    #Set cutoffs to use in step function: should be given in increasing order;
    cutoffs <- [1.64 1.96]'
    #Use a step function symmetric around zero
    symmetric <- 1
    #starting values for optimization
    Psihat0 <- [1,1,1,0,0,1,0,0]#[kappa, lambda, betap(1), betap(2)]

} else if (application == 11){

    data <- csvread('../Applications/PsychExperiments/cleaned_psych_data.csv')
    dofis1=data(:,7)==1
    use <- dofis1
    data <- data(use,:)

    Z <- [data(:,1)./data(:,2) data(:,3)./data(:,2)]
    sigmaZ2 <- data(:,4)./data(:,2)
    n <- size(Z,1)
    cluster_ID <- (1:n)'
    Studynames <- char.empty(0,n)

    X <- data(:,1)
    sigma <- data(:,2)
    C <- ones(length(X),1)

    identificationapproach <- 1){
    GMMapproach <- 1
    momentchoice <- 1
    name <- 'ReplicationPsychGMM'

    #Estimate baseline model, rather than running spec test
    spec_test <- 0
    #Set cutoffs to use in step function: should be given in increasing order;
    cutoffs <- [1.64 1.96]'
    #Use a step function symmetric around zero
    symmetric <- 1
    #starting values for optimization
    Psihat0 <- [1,1]#[betap(1), betap(2)]

} else if (application == 12){

    data <- csvread('../Applications/PsychExperiments/cleaned_psych_data.csv')
    dofis1=data(:,7)==1
    use <- dofis1
    data <- data(use,:)

    Z <- [data(:,1)./data(:,2) data(:,3)./data(:,2)]
    sigmaZ2 <- data(:,4)./data(:,2)
    n <- size(Z,1)
    cluster_ID <- (1:n)'
    Studynames <- char.empty(0,n)

    X <- data(:,1)
    sigma <- data(:,2)
    C <- ones(length(X),1)

    identificationapproach <- 1){
    GMMapproach <- 1
    momentchoice <- 2
    name <- 'ReplicationPsychGMMAlt'

    #Estimate baseline model, rather than running spec test
    spec_test <- 0
    #Set cutoffs to use in step function: should be given in increasing order;
    cutoffs <- 1.96
    #Use a step function symmetric around zero
    symmetric <- 1
    #starting values for optimization
    Psihat0 <- 1

} else if (application == 13){
    data <- csvread('../Applications/MinimumWagev2/cleaned_minwage_data.csv')
    X <- data(:,1)
    X <- -X#flipping sign so that negative employment effects mean X>0
    sigma <- data(:,2)
    n <- length(X)
    Studynames <- char.empty(0,n)
    cluster_ID <- data(:,3)
    C <- ones(length(X),1)

    includeinfigure <- logical(ones(n,1))
    includeinestimation <- logical(ones(n,1))

    identificationapproach <- 2){
    GMMapproach <- 0
    controls <- 0
    numerical_integration <- 0
    name <- 'MinimumWage'
    spec_test <- 0
    time_tr} <- 0

    #Set cutoffs to use in step function: should be given in increasing order;
    cutoffs <- [-1.96 0 1.96]'
    #Use a step function which may by asymmetric around zero
    symmetric <- 0
    symmetric_p <- 0
    #starting values for optimization
    asymmetric_likelihood_spec <- 2#Use a t model for latent distribution of true effects
    Psihat0 <- [0,1,5,1,1,1]#[thetabar, tau, nu, betap(1), betap(2), betap(3)]

} else if (application == 14){
    data <- csvread('../Applications/MinimumWagev2/cleaned_minwage_data.csv')
    X <- data(:,1)
    X <- -X#flipping sign so that negative employment effects mean X>0
    sigma <- data(:,2)
    cluster_ID <- data(:,3)

    published <- data(:,4)
    use=published==1
    X <- X(use,:)
    sigma <- sigma(use,:)
    cluster_ID <- cluster_ID(use,:)
    C <- ones(length(X),1)

    n <- length(X)
    Studynames <- char.empty(0,n)

    includeinfigure <- logical(ones(n,1))
    includeinestimation <- logical(ones(n,1))

    identificationapproach <- 2){
    GMMapproach <- 0
    numerical_integration <- 1
    name <- 'MinimumWagePublished'
    spec_test <- 0
    time_tr} <- 0

    #Set cutoffs to use in step function: should be given in increasing order;
    cutoffs <- [-1.96 0 1.96]'
    #Use a step function which may be asymmetric around zero
    symmetric <- 0
    symmetric_p <- 0
    #starting values for optimization
    asymmetric_likelihood_spec <- 2#Use a t model for latent distribution of true effects
    Psihat0 <- [0,1,5,1,1,1]

} else if (application == 15){
    data <- csvread('../Applications/MinimumWagev2/cleaned_minwage_data.csv')
    X <- data(:,1)
    X <- -X#flipping sign so that negative employment effects mean X>0
    sigma <- data(:,2)
    cluster_ID <- data(:,3)

    #     published=data(:,4);
    #     use=published==1;
    #     X=X(use,:);
    #     sigma=sigma(use,:);
    #     cluster_ID=cluster_ID(use,:);
    C <- ones(length(X),1)

    n <- length(X)
    Studynames <- char.empty(0,n)

    includeinfigure <- logical(ones(n,1))
    includeinestimation <- logical(ones(n,1))

    identificationapproach <- 2){
    GMMapproach <- 0
    numerical_integration <- 1
    name <- 'MinimumWageNoSelection'
    spec_test <- 0
    time_tr} <- 0

    #Set cutoffs to use in step function: should be given in increasing order;
    cutoffs <- 0
    #Use a step function which may be asymmetric around zero
    symmetric <- 0
    symmetric_p <- 1
    #starting values for optimization
    asymmetric_likelihood_spec <- 2#Use a t model for latent distribution of true effects
    Psihat0 <- [0,1,5]#[thetabar, tau, nu, betap(1), betap(2), betap(3)]

} else if (application == 16){
    data <- csvread('../Applications/MinimumWagev2/cleaned_minwage_data.csv')
    X <- data(:,1)
    X <- -X#flipping sign so that negative employment effects mean X>0
    sigma <- data(:,2)
    n <- length(X)
    Studynames <- char.empty(0,n)
    cluster_ID <- data(:,3)
    year <- data(:,5)-2013#Measure year relative to 2013, which is median
    C <- [ones(length(X),1) year]

    includeinfigure <- logical(ones(n,1))
    includeinestimation <- logical(ones(n,1))

    identificationapproach <- 2){
    GMMapproach <- 0
    numerical_integration <- 1
    controls <- 1
    name <- 'MinimumWageControls'
    spec_test <- 0
    time_tr} <- 1

    #Set cutoffs to use in step function: should be given in increasing order;
    cutoffs <- [-1.96 0 1.96]'
    #Use a step function which may by asymmetric around zero
    symmetric <- 0
    symmetric_p <- 0
    #starting values for optimization
    asymmetric_likelihood_spec <- 2
    Psihat0 <- [.02,.02,1.3,1,0,1,0,1,0]#[thetabar, tau, betap(1), betap(2), betap(3), betap(4), betap(5), betap(6)]
    #Psihat0=[.02,.02,1.3,.7,0,.27,0,.323,0];    #[thetabar, tau, betap(1), betap(2), betap(3), betap(4), betap(5), betap(6)]

} else if (application ==17){

    data <- csvread('../Applications/MinimumWagev2/cleaned_minwage_data.csv')
    X <- data(:,1)
    X <- -X#flipping sign so that negative employment effects mean X>0
    sigma <- data(:,2)
    n <- length(X)
    Studynames <- char.empty(0,n)
    cluster_ID <- data(:,3)
    C <- ones(length(X),1)

    includeinfigure <- logical(ones(n,1))
    includeinestimation <- logical(ones(n,1))

    identificationapproach <- 2){
    GMMapproach <- 1
    controls <- 0
    numerical_integration <- 0
    name <- 'MinimumWageGMM'
    spec_test <- 0

    #Set cutoffs to use in step function: should be given in increasing order;
    cutoffs <- [-1.96 0 1.96]'
    #Use a step function which may by asymmetric around zero
    symmetric <- 0
    symmetric_p <- 0
    #starting values for optimization
    Psihat0 <- [1,1,1]#[betap(1), betap(2), betap(3)]


} else if (application == 18){
    data <- csvread('../Applications/deworming/cleaned_deworming_data.csv')
    Studynames <- table2array(readtable('../Applications/deworming/DewormingLabels.csv'))
    X <- data(:,1)
    sigma <- data(:,2)
    n <- length(X)
    cluster_ID <- data(:,3)
    C <- ones(length(X),1)

    includeinfigure <- logical(ones(n,1))
    includeinestimation <- logical(ones(n,1))

    identificationapproach <- 2){
    GMMapproach <- 0
    asymmetric_likelihood_spec <- 1#Use a normal model for latent distribution of true effects
    controls <- 0
    numerical_integration <- 0
    name <- 'Deworming'
    spec_test <- 0

    #Set cutoffs to use in step function: should be given in increasing order;
    cutoffs <- [-1.96 0 1.96]
    #Use a step function symmetric around zero for p(z), but don't enforce
    #symmetric distribution of true effects
    symmetric <- 0
    symmetric_p <- 1
    #starting values for optimization
    Psihat0 <- [0,1,1]#[thetabar, tau, betap(1)]

} else if (application == 19){
    data <- csvread('../Applications/deworming/cleaned_deworming_data.csv')
    Studynames <- table2array(readtable('../Applications/deworming/DewormingLabels.csv'))
    X <- data(:,1)
    sigma <- data(:,2)
    n <- length(X)
    cluster_ID <- data(:,3)
    C <- ones(length(X),1)

    includeinfigure <- logical(ones(n,1))
    includeinestimation <- logical(ones(n,1))

    identificationapproach <- 2){
    GMMapproach <- 0
    asymmetric_likelihood_spec <- 1#Use a normal model for latent distribution of true effects
    controls <- 0
    numerical_integration <- 0
    name <- 'DewormingAsymmetric'
    spec_test <- 0

    #Set cutoffs to use in step function: should be given in increasing order;
    cutoffs <- [-1.96 0 1.96]
    #Allow an asymmetric specification of p(z)
    symmetric <- 0
    symmetric_p <- 0
    #starting values for optimization
    Psihat0 <- [0,1,1,1,1]#[thetabar, tau, betap(1), betap(2), betap(3)]

} else if (application == 20){
    data <- csvread('../Applications/deworming/cleaned_deworming_data.csv')
    Studynames <- table2array(readtable('../Applications/deworming/DewormingLabels.csv'))
    X <- data(:,1)
    sigma <- data(:,2)
    n <- length(X)
    cluster_ID <- data(:,3)
    C <- ones(length(X),1)

    includeinfigure <- logical(ones(n,1))
    includeinestimation <- logical(ones(n,1))

    identificationapproach <- 2){
    asymmetric_likelihood_spec <- 1#Use a normal model for latent distribution of true effects
    controls <- 0
    numerical_integration <- 0
    GMMapproach <- 0
    name <- 'DewormingRestrictedAsymmetric'
    spec_test <- 0

    #Set cutoffs to use in step function: should be given in increasing order;
    cutoffs <- 0
    #Allow an asymmetric specification of p(z)
    symmetric <- 0
    symmetric_p <- 0
    #starting values for optimization
    Psihat0 <- [0,1,1]#[thetabar, tau, betap(1)]

} else if (application == 21){
    data <- csvread('../Applications/deworming/cleaned_deworming_data.csv')
    Studynames <- table2array(readtable('../Applications/deworming/DewormingLabels.csv'))
    X <- data(:,1)
    sigma <- data(:,2)
    n <- length(X)
    cluster_ID <- data(:,3)
    C <- ones(length(X),1)

    includeinfigure <- logical(ones(n,1))
    includeinestimation <- logical(ones(n,1))

    identificationapproach <- 2){
    GMMapproach <- 0
    asymmetric_likelihood_spec <- 1#Use a normal model for latent distribution of true effects
    controls <- 0
    numerical_integration <- 0
    name <- 'DewormingNoSelection'
    spec_test <- 0

    #Set cutoffs to use in step function: should be given in increasing order;
    cutoffs <- 0
    #Use a step function symmetric around zero for p(z), but don't enforce
    #symmetric distribution of true effects
    symmetric <- 0
    symmetric_p <- 1
    #starting values for optimization
    Psihat0 <- [0,1]#[thetabar, tau]

} else if (application == 22){

    data <- csvread('../Applications/deworming/cleaned_deworming_data.csv')
    Studynames <- table2array(readtable('../Applications/deworming/DewormingLabels.csv'))
    X <- data(:,1)
    sigma <- data(:,2)
    n <- length(X)
    cluster_ID <- data(:,3)
    C <- ones(length(X),1)

    includeinfigure <- logical(ones(n,1))
    includeinestimation <- logical(ones(n,1))

    identificationapproach <- 2){
    GMMapproach <- 1
    asymmetric_likelihood_spec <- 1#Use a normal model for latent distribution of true effects
    controls <- 0
    numerical_integration <- 0
    name <- 'DewormingGMM'
    spec_test <- 0

    #Set cutoffs to use in step function: should be given in increasing order;
    cutoffs <- 1.96
    #Use a step function symmetric around zero for p(z)
    symmetric <- 1
    #starting values for optimization
    Psihat0 <- 1#[thetabar, tau, betap(1)]

   } else if (application == 23){

    Studynames <- textscan(fopen('../Applications/NatureScienceExperiments/sorted_names.csv','r'), '%s%[^\n\r]', 'Delimiter', ';')
    Studynames <- Studynames{1}(2:19,1)

    data <- csvread('../Applications/NatureScienceExperiments/cleaned_naturescience_data.csv')
    Z <- [data(:,1)./data(:,2) data(:,3)./data(:,2)]
    sigmaZ2 <- data(:,4)./data(:,2)

    n <- size(Z,1)
    cluster_ID <- (1:n)'

    X <- data(:,1)
    sigma <- data(:,2)
    C <- ones(length(X),1)

    identificationapproach <- 1){
    GMMapproach <- 0
    numerical_integration <- 0
    name <- 'ReplicationNatureScience'

    #Estimate baseline model, rather than running spec test
    spec_test <- 0
    #Set cutoffs to use in step function: should be given in increasing order;
    cutoffs <- [1.96 2.58]
    #Use a step function symmetric around zero
    symmetric <- 1
    #starting values for optimization
    Psihat0 <- [1,1,0,1]#[kappa, lambda, betap(1)]

} else if (application==24){
     #Studynames = textscan(fopen('../Applications/PsychExperiments/sorted_1stAuthor.csv','r'), '#s#[^\n\r]', 'Delimiter', ';');
    #Studynames=Studynames{1}

    data <- csvread('../Applications/PsychExperiments/cleaned_psych_data.csv')
    dofis1=data(:,7)==1
    use <- dofis1
    data <- data(use,:)

    Z <- [data(:,1)./data(:,2) data(:,3)./data(:,2)]
    sigmaZ2 <- data(:,4)./data(:,2)
    n <- size(Z,1)
    cluster_ID <- (1:n)'
    Studynames <- char.empty(0,n)

    X <- data(:,1)
    sigma <- data(:,2)
    C <- ones(length(X),1)

    identificationapproach <- 1){
    GMMapproach <- 0
    numerical_integration <- 0
    name <- 'ReplicationPsychSmall'

    #Estimate baseline model, rather than running spec test
    spec_test <- 0
    #Set cutoffs to use in step function: should be given in increasing order;
    cutoffs <- 1.96
    #Use a step function symmetric around zero
    symmetric <- 1
    #starting values for optimization
    Psihat0 <- [1,1,1]#[kappa, lambda, betap]

}


##
#producing figures
if (dofigures==1){
    if (identificationapproach==1){
        DescriptiveStats(Z,sigmaZ2, 1,name,symmetric,cluster_ID)
        DescriptiveStats(data(:,1),data(:,2), 2,[name 'Sanitycheck'],symmetric,cluster_ID)
        DescriptiveStatsCombined(Z,data(:,1),data(:,2),name,symmetric)
    } else {
        DescriptiveStats(X,sigma, 2,name,symmetric,cluster_ID)
    }
}

##
#estimating the model

if (doestimates==1){
    EstimatingSelection
}


##
#producing bias-corrected estimates and confidence sets

if (docorrections ==1){
    HorizontalBars
}
