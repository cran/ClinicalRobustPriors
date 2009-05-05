\name{Cauchy.Binomial}
\alias{Cauchy.Binomial}
\title{The Binomial Likelihood with Beta and Cauchy Priors}
\description{
Compute the distributions (prior, likelihood, posterior predictive and posterior) and moments
for the Beta/Binomial conjugate  model and Cauchy/Binomial robust model. The plots are processed
in Log-Odds and Theta Scale.}
\usage{
Cauchy.Binomial(n,x,a,b,m,min.value=NULL,max.value=NULL,iter=NULL)
}
\arguments{
\item{n}{sample size or number of observed patients.}
\item{x}{number of positive responses in n trials.}
\item{a}{the usual parameter of Beta prior and the number of positive responses in the prior information.}
\item{b}{the usual parameter of Beta prior and the number of negative responses in the prior information.}
\item{m}{number of additional patients for predictions}
\item{min.value}{minimum value in Log-Odds scale for the plots. The default min.value is 5.}
\item{max.value}{maximum value in Log-Odds scale for the plots. The default max.value is -5.}
\item{iter}{number of iterations in rejection sampling for the moments for the Cauchy/Binomial model. The default iter
is 10000.}
}
\references{
Fuquene, J. A., Cook, J. D. & Pericchi, L. R. (2008), A Case for Robust Bayesian priors with
Applications to Binary Clinical Trials. UT MD Anderson Cancer Center Department of Biostatistics 
Working Paper Series. Working Paper 44. 2008. 
http://www.bepress.com/mdandersonbiostat/paper44.

Spiegelhalter, D. J., Abrams, K. R. & Myles, J. P. (2004), Bayesian Approaches to
Clinical Trials and Health-Care Evaluation, Wiley, London.
}
\author{Jairo A. Fuquene P. <jairo.a.fuquene@uprrp.edu>}
\examples{
##############################################
# Example 1: sample and prior are in conflict
##############################################
Cauchy.Binomial(20,16,3,12,40)
##############################################
# Example 2: sample and prior are consistent
##############################################
Cauchy.Binomial(20,16,12,3,50,min.value=-5,max.value=5,iter=5000)
}