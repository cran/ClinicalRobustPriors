
Berger.Normal<-function(mu,ym,tau,sigma,n0,n,min.value=NULL,max.value=NULL,OR.xlim=NULL){

b<-tau^2
d<-sigma^2
if (tau==0) return("The scale of the prior should be bigger to 0")
if (sigma==0) return("The scale of the likelihood should be bigger to 0")
if (b<d) return("If b<d we should use the Cauchy.Normal function!")
min.value <- if (is.null(min.value)) min.value <- -5 else min.value <- min.value
max.value <- if (is.null(max.value)) max.value <- 5 else max.value <- max.value
OR.xlim <- if (is.null(OR.xlim)) OR.xlim <- c(0,2) else OR.xlim <- OR.xlim
mu <- if (mu==0) mu <- 0.001 else mu <- mu


sdn<-tau
sd_ym<-sigma


lambda<-seq(min.value, max.value, len=10000)
theta<-exp(lambda)/(1+exp(lambda))

#-----------------------------------------------------
# Berger's Prior
#-----------------------------------------------------

Berger.prior<-rep(0,10000)

for(i in 1:10000){
integrand <- function(v) {(1/(2*sqrt(v)))*dnorm(lambda[i], mean=mu, sd=sqrt((d+b)/(2*v)- d), log = FALSE)} 
integra<-as.vector(unlist(integrate(integrand, lower = 0, upper = 1)))
Berger.prior[i]<-as.numeric(integra[1])
}

#------------------------------------------------------
# Normal Prior
#------------------------------------------------------

Normal.prior<-dnorm(lambda, mean=mu, sd=sdn) 
integrand <- function(lambda) {dnorm(lambda, mean=mu, sd=sdn)}
integra<-as.vector(unlist(integrate(integrand, lower = min.value, upper = max.value)))
Normal.prior<-Normal.prior/(as.numeric(integra[1]))

#-------------------------------------------------------
# Normalized Likelihood
#-------------------------------------------------------

lik.normalized<-dnorm(lambda, mean=ym, sd=sd_ym) 
integrand <- function(lambda) {dnorm(lambda, mean=ym, sd=sd_ym)}
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

#---------------------------------------------------------
# Berger/Normal Posterior Moments and Posterior Model
#---------------------------------------------------------

factor<- (sqrt(d)*sqrt(d+b))/(sqrt(2)*(ym-mu)^2)
posterior.Robust.Normal<-(Berger.prior*exp(-(ym-lambda)^2/(2*d)))/(factor*(1-exp(-(ym-mu)^2/(d+b))))

F_x.Robust<-exp((ym-mu)^2/(d+b))
E_Robust.Normal<- ym + (2*d*(ym-mu))/((d+b)*(F_x.Robust-1)) - 2*d/(ym-mu)
v_Robust.Normal<- d - (d^2)*((4*(ym-mu)^2*F_x.Robust-2*(d+b)*(F_x.Robust-1))/((d+b)^2*(F_x.Robust-1)^2)-2/(ym-mu)^2)
sd_Robust.Normal<-sqrt(v_Robust.Normal)

#-------------------
# Credible interval 
#-------------------

credible.int95.Normal.Normal.Inf<-E_Normal.Normal-qnorm(0.975,0,1)*sd_Normal.Normal
credible.int95.Normal.Normal.Sup<-E_Normal.Normal+qnorm(0.975,0,1)*sd_Normal.Normal
credible.int95.Berger.Normal.Inf<-E_Robust.Normal-qnorm(0.975,0,1)*sd_Robust.Normal
credible.int95.Berger.Normal.Sup<-E_Robust.Normal+qnorm(0.975,0,1)*sd_Robust.Normal
credible.int95.likelihood.Inf<-ym-qnorm(0.975,0,1)*sd_ym
credible.int95.likelihood.Sup<-ym+qnorm(0.975,0,1)*sd_ym

#----------------------------------------------------
# Figures
#----------------------------------------------------

par(mfrow=c(3,2))

distributions1<-cbind(Normal.prior,
lik.normalized)

distributions2<-cbind(lik.normalized,
posterior.Normal.Normal)

distributions3<-cbind(Berger.prior,
lik.normalized)

distributions4<-cbind(lik.normalized,
posterior.Robust.Normal)

matplot(lambda, distributions1, type="l", lty=1, xlab="Log-Odds Scale", ylab="Distribution",col=c("red", "darkgreen"),lwd=1.5)
mtext("Normal Prior",3,col = "red",line=1.0,cex=0.7) 
mtext("Normal Likelihood",3,line=0.2,col = "darkgreen",cex=0.7) 

matplot(lambda, distributions2, type="l", lty=1, xlab="Log-Odds Scale",ylab="Distribution", col=c("darkgreen", "blue"))
mtext("Normal/Normal Model",3,col = "blue",line=1.0,cex=0.7) 
mtext("Normal Likelihood",3,line=0.2,col = "darkgreen",cex=0.7) 

matplot(lambda, distributions3, type="l", lty=1, xlab="Log-Odds Scale", ylab="Distribution",col=c("red", "darkgreen"))
mtext("Berger's Prior",3,col = "red",line=1.0,cex=0.7) 
mtext("Normal Likelihood",3,line=0.2,col = "darkgreen",cex=0.7) 

matplot(lambda, distributions4, type="l", lty=1, xlab="Log-Odds Scale",ylab="Distribution", col=c("darkgreen", "blue"))
mtext("Berger/Normal Model",3,col = "blue",line=1.0,cex=0.7) 
mtext("Normal Likelihood",3,line=0.2,col = "darkgreen",cex=0.7) 


matplot(exp(lambda), distributions3, type="l", lty=1, xlab="OR Scale", ylab="Distribution",col=c("red", "darkgreen"),xlim=OR.xlim)
mtext("Berger's Prior",3,col = "red",line=1.0,cex=0.7) 
mtext("Normal Likelihood",3,line=0.2,col = "darkgreen",cex=0.7) 

matplot(exp(lambda), posterior.Robust.Normal, type="l", lty=1, xlab="OR Scale",ylab="Distribution", main="Berger/Normal Model",col.main="blue",col="blue"
,xlim=OR.xlim)

###############################################
# Results
###############################################

Location.Log.Odds.scale <- round(c(mu,ym,E_Normal.Normal,E_Robust.Normal),digits=2)
Scale.Log.Odds.scale <- round(c(sqrt(b),sqrt(d),sd_Normal.Normal,sd_Robust.Normal),digits=2)
Location.OR.scale <- round(exp(Location.Log.Odds.scale),digits=2)
Cred.I.Inf.OR<-round(exp(c(NA,credible.int95.likelihood.Inf,credible.int95.Normal.Normal.Inf,credible.int95.Berger.Normal.Inf)),digits=2)
Cred.I.Sup.OR<-round(exp(c(NA,credible.int95.likelihood.Sup,credible.int95.Normal.Normal.Sup,credible.int95.Berger.Normal.Sup)),digits=2)

Results <- rbind(Location.Log.Odds.scale,Scale.Log.Odds.scale,Location.OR.scale,Cred.I.Inf.OR, Cred.I.Sup.OR)
colnames(Results) <- c("Prior","Normalized Lik.","Normal/Normal","Berger/Normal")
rownames(Results) <- c("Location Log Odds Scale","Scale", "Location OR Scale","95% Credible Interval OR Inf.","95% Credible Interval OR Sup.")

print("Results in the Log odds and OR scale")
return(Results)

}












