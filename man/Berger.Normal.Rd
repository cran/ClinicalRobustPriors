\name{Berger.Normal}
\alias{Berger.Normal}
\title{Normal Log-Odds with Normal and Berger's Priors}
\description{
Compute the distributions (prior, likelihood and posterior) and moments
for the Normal/Normal conjugate  model and Berger/Normal robust model. The plots are processed
in Log-Odds and OR Scale.
}
\usage{
Berger.Normal(mu,ym,tau,sigma,n0,n,min.value=NULL,max.value=NULL,OR.xlim=NULL)
}
\arguments{
\item{mu}{the location parameter in Log-Odds Scale for Normal and Berger's Priors.}
\item{ym}{the location parameter in Log-Odds Scale for Normal Likelihood.}
\item{tau}{the scale parameter in Log-Odds for Normal and Berger's Priors.}
\item{sigma}{the scale parameter in Log-Odds Normal Likelihood.}
\item{n0}{number of prior observations.}
\item{n}{sample size.}
\item{min.value}{minimum value in Log-Odds scale for the plots. The default min.value is 5.}
\item{max.value}{maximum value in Log-Odds scale for the plots. The default max.value is -5.}
\item{OR.xlim}{maximum value for the OR scale in the plot of the Berger/Normal model. The default OR.xlim 
is c(0,2).}
}
\references{
Berger, J. O. (1985), Statistical Decision Theory and Bayesian Analysis, second
edn, Springer-Verlag.

Fuquene, J. A., Cook, J. D. & Pericchi, L. R. (2008), A Case for Robust Bayesian priors with
Applications to Binary Clinical Trials. UT MD Anderson Cancer Center Department of Biostatistics 
Working Paper Series. Working Paper 44. 2008. 
http://www.bepress.com/mdandersonbiostat/paper44.

Spiegelhalter, D. J., Abrams, K. R. & Myles, J. P. (2004), Bayesian Approaches
to Clinical Trials and Health-Care Evaluation, Wiley, London.
}
\author{Jairo A. Fuquene P. <jairo.a.fuquene@uprrp.edu>}
\examples{
#######################
# Example 1: 
#######################
Berger.Normal(-1.97,-0.73,0.31,0.15,406,170)
#######################
# Example 2
#######################
Berger.Normal(-1.97,-0.73,0.31,0.15,406,170,min.value=-3.5,max.value=0,OR.xlim=c(0,1))
#######################
# Example 3
#######################
Berger.Normal(-1.97,-0.73,0.15,0.31,406,170)
}


