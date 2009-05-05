Cauchy.Binomial<-function(n,x,a,b,m,min.value=NULL,max.value=NULL,iter=NULL){

iter <- if (is.null(iter)) iter <- 5000 else iter <- iter
min.value <- if (is.null(min.value)) min.value <- -5 else min.value <- min.value
max.value <- if (is.null(max.value)) max.value <- 5 else max.value <- max.value

lambda<-seq(min.value, max.value, len=10000)
theta<-exp(lambda)/(1+exp(lambda))

#----------------------------------------------------
# Moments of the Beta Prior
#----------------------------------------------------

location.Beta.Cauchy<-digamma(a) - digamma(b) 
scale.Beta.Cauchy_2<-trigamma(a) + trigamma(b)
scale.Beta.Cauchy<-sqrt(scale.Beta.Cauchy_2)

#-----------------------------------------------------
# Beta Prior
#-----------------------------------------------------

prior.Beta.Binomial<-exp(lgamma(a+b)-lgamma(a)-lgamma(b))*(exp(lambda)/(1+exp(lambda)))^a*(1/(1+exp(lambda)))^b
integrand <- function(lambda) {(exp(lambda)/(1+exp(lambda)))^a*(1/(1+exp(lambda)))^b}
integra<-as.vector(unlist(integrate(integrand, lower = min.value, upper = max.value)))
prior.Beta.Binomial<-prior.Beta.Binomial/(exp(lgamma(a+b)-lgamma(a)-lgamma(b))*as.numeric(integra[1]))

#-----------------------------------------------------
# Cauchy Prior
#-----------------------------------------------------

Cauchy.prior<-dcauchy(lambda,location.Beta.Cauchy,scale.Beta.Cauchy)
integrand <- function(lambda) {dcauchy(lambda,location.Beta.Cauchy,scale.Beta.Cauchy)}
integra<-as.vector(unlist(integrate(integrand, lower = min.value, upper = max.value)))
Cauchy.prior<-Cauchy.prior/(as.numeric(integra[1]))

#------------------------------------------------------
# Normalized Likelihood
#------------------------------------------------------

lik.normalized<-(exp(lgamma(n+1)-lgamma(n-x+1)-lgamma(x+1))*x*(n-x)/n)*exp(x*lambda-n*log(1+exp(lambda)))
integrand <- function(lambda) {(x*(n-x)/n)*exp(x*lambda-n*log(1+exp(lambda)))}
integra<-as.vector(unlist(integrate(integrand, lower = min.value, upper = max.value)))
lik.normalized<-lik.normalized/(exp(lgamma(n+1)-lgamma(n-x+1)-lgamma(x+1))*as.numeric(integra[1]))

#------------------------------------------------------
# Beta/Binomial Posterior Model
#------------------------------------------------------

posterior.Beta.Binomial<-exp(lgamma(a+b+n)-lgamma(a+x)-lgamma(n-x+b))*exp(x*lambda + a*log(((exp(lambda))/(1+exp(lambda))))-(n+b)*log(1+exp(lambda)))
integrand <- function(lambda) {exp(x*lambda + a*log(((exp(lambda))/(1+exp(lambda))))-(n+b)*log(1+exp(lambda)))}
integra<-as.vector(unlist(integrate(integrand, lower = min.value, upper = max.value)))
posterior.Beta.Binomial<-posterior.Beta.Binomial/(exp(lgamma(a+b+n)-lgamma(a+x)-lgamma(n-x+b))*as.numeric(integra[1]))

#-----------------------------------------------------------
# Posterior Moments of the Cauchy/Binomial Posterior Model
#-----------------------------------------------------------

E_Beta.Binomial<-digamma(x+a)-digamma(n-x+b) 
v_Beta.Binomial<-trigamma(x+a) +trigamma(n-x+b)
sd_Beta.Binomial<-sqrt(v_Beta.Binomial)

#------------------------------------------------------
# Predictive Cauchy/Binomial Posterior Model
#------------------------------------------------------

integrand <- function(lambda) {exp(x*lambda-n*log(1+exp(lambda)))*dcauchy(lambda,location.Beta.Cauchy,scale.Beta.Cauchy)}
integra<-as.vector(unlist(integrate(integrand, lower = min.value, upper = max.value)))
predictive.Cauchy.Binomial<-as.numeric(integra[1])

#------------------------------------------------------
# Cauchy/Binomial Posterior Model
#------------------------------------------------------

aproximation.Cauchy.Binomial<-(exp(x*lambda-n*log(1+exp(lambda)))*dcauchy(lambda,location.Beta.Cauchy,scale.Beta.Cauchy))/(predictive.Cauchy.Binomial)
integrand <- function(lambda) {dcauchy(lambda,location.Beta.Cauchy,scale.Beta.Cauchy)}
integra<-as.vector(unlist(integrate(integrand, lower = min.value, upper = max.value)))
aproximation.Cauchy.Binomial<-aproximation.Cauchy.Binomial/(as.numeric(integra[1]))

#----------------------------------------------------------
# Posterior Moments of the Cauchy/Binomial Posterior Model
#----------------------------------------------------------


lambda_alea_calcula<-rep(0,iter)
lambda_alea<-rep(0,iter)
samples<-rep(0,iter)
aleatorio<-rep(0,iter)
post_aprox<-rep(0,iter)
theta_s<-rep(0,iter)
predictive_Bin<-rep(0,iter)
x_b<-x/n
max<-log((x_b)/(1-x_b))
M<-exp(x*max-n*log(1+exp(max)))
MLE<-log(x_b/(1-x_b))



for(i in 1:iter){

while(samples[i]==0) {

lambda_alea[i]<-rcauchy(1, location = location.Beta.Cauchy, scale = scale.Beta.Cauchy)
lambda_alea_calcula[i]<-dcauchy(lambda_alea[i], location = location.Beta.Cauchy, scale = scale.Beta.Cauchy) 
aleatorio[i]<-runif(1, min=0, max=1) 
post_aprox[i]<- exp(x*lambda_alea[i]-n*log(1+exp(lambda_alea[i])))*dcauchy(lambda_alea[i],location.Beta.Cauchy,scale.Beta.Cauchy)

if  (M*aleatorio[i]*lambda_alea_calcula[i]<=post_aprox[i]) samples[i] <- lambda_alea[i]

}
}

E_rejection.Cauchy.Binomial<-mean(samples)
var_rejection.Cauchy.Binomial<-var(samples)
sd_rejection.Cauchy.Binomial<-sqrt(var_rejection.Cauchy.Binomial)

#----------------------------------------------------------
# Predictive
#----------------------------------------------------------

theta_s<-exp(samples)/(1+exp(samples))
predictive_Bin<-rbinom(iter,m,theta_s)

tabla_freq_c<-table(predictive_Bin)/iter
pred<-as.vector(tabla_freq_c)
Xm<-sort(unique(predictive_Bin))

predictive.mean.Cauchy.Binomial.Rejection<-round(mean(predictive_Bin))

predictive.var.Cauchy.Binomial.Rejection<-var(predictive_Bin)
predictive.sd.Cauchy.Binomial.Rejection<-sqrt(predictive.var.Cauchy.Binomial.Rejection)

predictive.mean.Beta.Binomial<-round(m*((a+x)/(a+b+n)))
Predictive.var.Beta.Binomial<-((m*(a+x)*(b+n-x))/(a+b+n)^2)*((a+b+2*n)/(a+b+n+1))
Predictive.sd.Beta.Binomial<-sqrt(Predictive.var.Beta.Binomial)

#---------------------------------------------------
#Credible interval for the posterior moments 95%
#---------------------------------------------------

credible.int95.Beta.Binomial.Inf<-E_Beta.Binomial-qnorm(0.975,0,1)*sd_Beta.Binomial
credible.int95.Beta.Binomial.Sup<-E_Beta.Binomial+qnorm(0.975,0,1)*sd_Beta.Binomial
credible.int95.Beta.Binomial.Theta.Inf<-exp(credible.int95.Beta.Binomial.Inf)/(1+exp(credible.int95.Beta.Binomial.Inf))
credible.int95.Beta.Binomial.Theta.Sup<-exp(credible.int95.Beta.Binomial.Sup)/(1+exp(credible.int95.Beta.Binomial.Sup))

credible.int95.Cauchy.Binomial.Inf<-E_rejection.Cauchy.Binomial-qnorm(0.975,0,1)*sd_rejection.Cauchy.Binomial 
credible.int95.Cauchy.Binomial.Sup<-E_rejection.Cauchy.Binomial+qnorm(0.975,0,1)*sd_rejection.Cauchy.Binomial
credible.int95.Cauchy.Binomial.Theta.Inf<-exp(credible.int95.Cauchy.Binomial.Inf)/(1+exp(credible.int95.Cauchy.Binomial.Inf))
credible.int95.Cauchy.Binomial.Theta.Sup<-exp(credible.int95.Cauchy.Binomial.Sup)/(1+exp(credible.int95.Cauchy.Binomial.Sup))

#----------------------------------------------------
# Figures
#----------------------------------------------------

par(mfrow=c(3,2))

distributions1<-cbind( 
          prior.Beta.Binomial,   
          lik.normalized
)

distributions2<-cbind( 
          lik.normalized,
          posterior.Beta.Binomial
)



distributions3<-cbind( 
          Cauchy.prior,   
          lik.normalized
)

distributions4<-cbind( 
          lik.normalized,
          aproximation.Cauchy.Binomial
)

matplot(lambda, distributions1, type="l", lty=1, xlab="Log-Odds Scale", ylab="Distribution",col=c("red", "darkgreen"),lwd=1.5)
mtext("Beta Prior",3,col = "red",line=1.0,cex=0.7) 
mtext("Binomial Likelihood",3,line=0.2,col = "darkgreen",cex=0.7) 

matplot(lambda, distributions2, type="l", lty=1, xlab="Log-Odds Scale",ylab="Distribution", col=c("darkgreen", "blue"))
mtext("Beta/Binomial Model",3,col = "blue",line=1.0,cex=0.7) 
mtext("Binomial Likelihood",3,line=0.2,col = "darkgreen",cex=0.7) 

matplot(lambda, distributions3, type="l", lty=1, xlab="Log-Odds Scale", ylab="Distribution",col=c("red", "darkgreen"))
mtext("Cauchy Prior",3,col = "red",line=1.0,cex=0.7) 
mtext("Binomial Likelihood",3,line=0.2,col = "darkgreen",cex=0.7) 

matplot(lambda, distributions4, type="l", lty=1, xlab="Log-Odds Scale",ylab="Distribution", col=c("darkgreen", "blue"))
mtext("Cauchy/Binomial Model",3,col = "blue",line=1.0,cex=0.7) 
mtext("Binomial Likelihood",3,line=0.2,col = "darkgreen",cex=0.7) 

matplot(theta, aproximation.Cauchy.Binomial, type="l", lty=1, xlab="Theta Scale",ylab="Distribution", main="Cauchy/Binomial Model",col.main="blue",col="blue")

plot(Xm,pred,type="h",ylab = "",main="Posterior Predictive Distribution",lwd=2)

#----------------------------------------------------------
# Results
#----------------------------------------------------------


prior.data <- round(rbind(MLE,x_b,location.Beta.Cauchy,scale.Beta.Cauchy),digits=2)
colnames(prior.data) <- c("Value")
rownames(prior.data) <- c("MLE of the Log Odds","MLE of Theta", "Location of the Beta and Cauchy Prior", "Scale of the Beta and Cauchy Prior")



Location.Log.Odds.scale <- round(c(E_Beta.Binomial,E_rejection.Cauchy.Binomial),digits=2)
Scale.Log.Odds.scale <- round(c(sd_Beta.Binomial,sd_rejection.Cauchy.Binomial),digits=2)
Location.Theta.scale <- round(exp(Location.Log.Odds.scale)/(1+exp(Location.Log.Odds.scale)),digits=2)
Predictive.Mean<- c(predictive.mean.Beta.Binomial,predictive.mean.Cauchy.Binomial.Rejection)
Predictive.sd<- round(c(Predictive.sd.Beta.Binomial,predictive.sd.Cauchy.Binomial.Rejection),digits=2)  
Cred.I.Inf.Theta<-round(c(credible.int95.Beta.Binomial.Theta.Inf,credible.int95.Cauchy.Binomial.Theta.Inf),digits=2)
Cred.I.Sup.Theta<-round(c(credible.int95.Beta.Binomial.Theta.Sup,credible.int95.Cauchy.Binomial.Theta.Sup),digits=2)

Results <- rbind(Location.Log.Odds.scale,Scale.Log.Odds.scale,Location.Theta.scale,Cred.I.Inf.Theta, Cred.I.Sup.Theta,round(Predictive.Mean),Predictive.sd)
colnames(Results) <- c("Beta/Binomial","Cauchy/Binomial")
rownames(Results) <- c("Posterior Expectation in Log Odds Scale","Posterior Standard Deviation in Log Odds ", "Posterior Expectation in Theta Scale","95% Credible Interval Theta Inf.","95% Credible Interval Theta Sup.", paste("Predictive Mean for m =",m), paste("Predictive sd for m =",m))



print("Results in the Log Odds and Theta Scale")

return(list(prior.data,Results))

print("Results of the Cauchy/Binomial and Beta/Binomial
models in the Log Odds and Theta scale")

return(paste(Results,prior.data))


}


