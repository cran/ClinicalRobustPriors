Cauchy.Normal<-function(mu,ym,tau,sigma,n0,n,min.value=NULL,max.value=NULL,OR.xlim=NULL){

min.value <- if (is.null(min.value)) min.value <- -5 else min.value <- min.value
max.value <- if (is.null(max.value)) max.value <- 5 else max.value <- max.value
OR.xlim <- if (is.null(OR.xlim)) OR.xlim <- c(0,2) else OR.xlim <- OR.xlim


lambda<-seq(min.value, max.value, len=10000)
theta<-exp(lambda)/(1+exp(lambda))

#-----------------------------------------------------
# Cauchy prior 
#-----------------------------------------------------

Cauchy.prior<-dcauchy(lambda,mu,tau)
integrand <- function(lambda) {dcauchy(lambda,mu,tau)}
integra<-as.vector(unlist(integrate(integrand, lower = min.value, upper = max.value)))
Cauchy.prior<-Cauchy.prior/(as.numeric(integra[1]))

#------------------------------------------------------
# Normal Prior
#------------------------------------------------------

Normal.prior<-dnorm(lambda, mean=mu, sd=tau) 
integrand <- function(lambda) {dnorm(lambda, mean=mu, sd=tau)}
integra<-as.vector(unlist(integrate(integrand, lower = min.value, upper = max.value)))
Normal.prior<-Normal.prior/(as.numeric(integra[1]))

#-------------------------------------------------------
# Normalized Likelihood
#-------------------------------------------------------

lik.normalized<-dnorm(lambda, mean=ym, sd=sigma) 
integrand <- function(lambda) {dnorm(lambda, mean=ym, sd=sigma)}
integra<-as.vector(unlist(integrate(integrand, lower = min.value, upper = max.value)))
lik.normalized<-lik.normalized/(as.numeric(integra[1]))

#-------------------------------------------------------
# Normal/Normal Posterior Moments and Posterior Model
#-------------------------------------------------------

E_Normal.Normal <- (n0*mu+n*ym)/(n0+n)
V_Normal.Normal <- 4/(n0+n)
sd_Normal.Normal <- sqrt(V_Normal.Normal)

posterior.Normal.Normal<-dnorm(lambda, mean=E_Normal.Normal, sd=sd_Normal.Normal) 
integrand <- function(lambda) {dnorm(lambda, mean=E_Normal.Normal, sd=sd_Normal.Normal)}
integra<-as.vector(unlist(integrate(integrand, lower = min.value, upper = max.value)))
posterior.Normal.Normal<-posterior.Normal.Normal/(as.numeric(integra[1]))

#-----------------------------------------------------
# Cauchy/Normal Posterior Moments and Posterior Model
#----------------------------------------------------

aproximation.Cauchy.Normal<-dnorm(ym, mean=lambda, sd=sigma)*dcauchy(lambda, location = mu, scale = tau)
integrand <- function(lambda) {dnorm(ym, mean=lambda, sd=sigma)*dcauchy(lambda, location = mu, scale = tau)}
integra<-as.vector(unlist(integrate(integrand, lower = min.value, upper = max.value)))
aproximation.Cauchy.Normal<-aproximation.Cauchy.Normal/(as.numeric(integra[1]))

t_0<-gamma(1) / (tau*sqrt(pi)* gamma(1/2))* (1 + (ym-mu)^2/(tau^2))^-(1)

tau1<-tau/sqrt(3)
t_1<-gamma(2) / (tau1*sqrt(3*pi)* gamma(3/2))* (1 + (ym-mu)^2/(3*tau1^2))^-(2)

tau2<-tau/sqrt(5)
t_2<-gamma(3) / (tau2*sqrt(5*pi)* gamma(5/2))* (1 + (ym-mu)^2/(5*tau2^2))^-(3)

E_Cauchy.Normal<-ym - (2*sigma^2*(ym-mu))/(tau^2+(ym-mu)^2)
var_Cauchy.Normal<-sigma^2-sigma^4*(((ym-mu)^2*t_1^2/tau^4)/t_0^2 - (((ym-mu)^2*3*t_2/tau^4)-(t_1/tau^2))/t_0) 
sd_Cauchy.Normal<-sqrt(var_Cauchy.Normal)

#-------------------
# Credible interval 
#-------------------

credible.int95.Normal.Normal.Inf<-E_Normal.Normal-qnorm(0.975,0,1)*sd_Normal.Normal
credible.int95.Normal.Normal.Sup<-E_Normal.Normal+qnorm(0.975,0,1)*sd_Normal.Normal

credible.int95.Cauchy.Normal.Inf<-E_Cauchy.Normal-qnorm(0.975,0,1)*sd_Cauchy.Normal
credible.int95.Cauchy.Normal.Sup<-E_Cauchy.Normal+qnorm(0.975,0,1)*sd_Cauchy.Normal

credible.int95.likelihood.Inf<-ym-qnorm(0.975,0,1)*sigma
credible.int95.likelihood.Sup<-ym+qnorm(0.975,0,1)*sigma

#----------------------------------------------------
# Figures
#----------------------------------------------------

par(mfrow=c(3,2))

distributions1<-cbind(Normal.prior,
lik.normalized)

distributions2<-cbind(lik.normalized,
posterior.Normal.Normal)

distributions3<-cbind(Cauchy.prior,
lik.normalized)

distributions4<-cbind(lik.normalized,
aproximation.Cauchy.Normal)


matplot(lambda, distributions1, type="l", lty=1, xlab="Log-Odds Scale", ylab="Distribution",col=c("red", "darkgreen"),lwd=1.5)
mtext("Normal Prior",3,col = "red",line=1.0,cex=0.7) 
mtext("Normal Likelihood",3,line=0.2,col = "darkgreen",cex=0.7) 

matplot(lambda, distributions2, type="l", lty=1, xlab="Log-Odds Scale",ylab="Distribution", col=c("darkgreen", "blue"))
mtext("Normal/Normal Model",3,col = "blue",line=1.0,cex=0.7) 
mtext("Normal Likelihood",3,line=0.2,col = "darkgreen",cex=0.7) 

matplot(lambda, distributions3, type="l", lty=1, xlab="Log-Odds Scale", ylab="Distribution",col=c("red", "darkgreen"))
mtext("Cauchy Prior",3,col = "red",line=1.0,cex=0.7) 
mtext("Normal Likelihood",3,line=0.2,col = "darkgreen",cex=0.7) 

matplot(lambda, distributions4, type="l", lty=1, xlab="Log-Odds Scale",ylab="Distribution", col=c("darkgreen", "blue"))
mtext("Cauchy/Normal Model",3,col = "blue",line=1.0,cex=0.7) 
mtext("Normal Likelihood",3,line=0.2,col = "darkgreen",cex=0.7) 


matplot(exp(lambda), distributions3, type="l", lty=1, xlab="OR Scale", ylab="Distribution",col=c("red", "darkgreen"),xlim=OR.xlim)
mtext("Cauchy Prior",3,col = "red",line=1.0,cex=0.7) 
mtext("Normal Likelihood",3,line=0.2,col = "darkgreen",cex=0.7) 

matplot(exp(lambda), aproximation.Cauchy.Normal, type="l", lty=1, xlab="OR Scale",ylab="Distribution", main="Cauchy/Normal Model",col.main="blue",col="blue"
,xlim=OR.xlim)

###############################################
# Results
###############################################

Location.Log.Odds.scale <- round(c(mu,ym,E_Normal.Normal,E_Cauchy.Normal),digits=2)
Scale.Log.Odds.scale <- round(c(tau,sigma,sd_Normal.Normal,sd_Cauchy.Normal),digits=2)
Location.OR.scale <- round(exp(Location.Log.Odds.scale),digits=2)
Cred.I.Inf.OR<-round(exp(c(NA,credible.int95.likelihood.Inf,credible.int95.Normal.Normal.Inf,credible.int95.Cauchy.Normal.Inf)),digits=2)
Cred.I.Sup.OR<-round(exp(c(NA,credible.int95.likelihood.Sup,credible.int95.Normal.Normal.Sup,credible.int95.Cauchy.Normal.Sup)),digits=2)

Results <- rbind(Location.Log.Odds.scale,Scale.Log.Odds.scale,Location.OR.scale,Cred.I.Inf.OR, Cred.I.Sup.OR)
colnames(Results) <- c("Prior","Normalized Lik.","Normal/Normal","Cauchy/Normal")
rownames(Results) <- c("Location Log Odds Scale","Scale", "Location OR Scale","95% Credible Interval OR Inf.","95% Credible Interval OR Sup.")

print("Results in the log odds and OR scale")
return(Results)

}
