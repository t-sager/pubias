Step_function_normal_cdf <- function(X,theta,sigma,betap,cutoffs,symmetric) {
    #Arguments:
    #X: point at which to evaluate cdf
    #theta: parameter value under which to evaluate cdf
    #sigma:standarad deviation of (untruncated) normal variable
    #sigma: stdev of distribution pi of mu
    #cutoffs: vector of thresholds for step function p, in increasing order
    #cutoffs are given in terms of X, not z statistics
    #Symmetric: dummy indicating whether publication probability is symmetric around zero.
    #In symmetric case, cutoffs should include only positive values
    #Note: publication probability for largest category (i.e. for those points beyond largest cutoff) normalized to one.

    if (length(betap)!=(length(cutoffs)+1)) {
        stop("length of betap must be one greater than length of cutoffs ")
    }

    #For symmetric case, create symmetrized version of cutoffs and coefficients
    if (symmetric==1) {
        cutoffs_u=matrix(0,2*length(cutoffs),1)
        betap_u=matrix(0,1,2*length(cutoffs))
        for (n in c(1:length(cutoffs))) {
            cutoffs_u[n]=-cutoffs[length(cutoffs)+1-n]
            betap_u[n]=betap[length(betap)+1-n]
        }
        for (n in c(1:length(cutoffs))) {
            cutoffs_u[length(cutoffs)+n]=cutoffs[n]
            betap_u[length(cutoffs)+n]=betap[n]
        }
        betap_u[length(betap_u)+1]=1
        betap_u=as.matrix(betap_u)
    } else {
        cutoffs_u=cutoffs
        betap_u=betap
    }
    betap_u=as.matrix(betap_u)

    #Calculate denominatore in cdf
    prob_vec=matrix(0,length(cutoffs_u)+2,1)
    for (m in c(1:length(cutoffs_u))) {
        prob_vec[m+1]=pnorm((cutoffs_u[m]-theta)/sigma)
    }
    prob_vec[length(prob_vec)]=1
    mean_Z1=as.matrix(prob_vec[2:length(prob_vec),1])-as.matrix(prob_vec[1:(length(prob_vec)-1),1])
    denominator=t(mean_Z1)%*%betap_u

    #Calculate numerator in cdf
    cutoffs_u[length(cutoffs_u)+1]=Inf
    if (X<=cutoffs_u[1]) {
        numerator=pnorm((X-theta)/sigma)*betap_u[1]
    } else {
        numerator=pnorm((cutoffs_u[1]-theta)/sigma)*betap_u[1]
        m=1
        while (X>cutoffs_u[m]) {
            Xcap=min(X,cutoffs_u[m+1])
            numerator=numerator+(pnorm((Xcap-theta)/sigma)-pnorm((cutoffs_u[m]-theta)/sigma))*betap_u[m+1]
            m=m+1
        }
    }
    #evaluate cdf
    cdf=numerator/denominator

    return(cdf)
}
