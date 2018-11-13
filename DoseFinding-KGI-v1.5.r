library(mvtnorm)
y=c() #response vector
dose=c() #dose vector
J=11 #number of doses	
patient=61#number of patients
#true_theta=c(10, 15, 20, 25, 30) #true mean responses
true_theta=c(10, 11, 13, 18, 22, 23, 22, 18, 13, 11, 10)
true_sigma=1 #true deviation of responses
z_j=0 #optimal dose
mu_0=c(10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10) #initial hyperparameter
sigma_0=diag(100, nrow=J, ncol=J) #initial hyperparameter
mu_n=lapply(1:patient, function(x) c())
sigma_n=lapply(1:patient, function(x) matrix(0, nrow=J, ncol=J))
theta_estimate=matrix(NA, nrow=patient, ncol=J)
var_target=c()
M=100
T=1000
start_time=0
end_time=0
####
evolution_eq <-function(mu_n, sigma_n, new_y, z_j){
	e_z = rep(0, J)
	e_z[z_j] <- 1
	sigma_tilde <- (sigma_n %*% e_z)*(1/sqrt(true_sigma+sigma_n[z_j, z_j]))
	new_sigma <- sigma_n - sigma_tilde %*% t(sigma_tilde)
	Var_yF <- true_sigma + sigma_n[z_j, z_j]
	new_X <- (new_y - mu_n[z_j])/sqrt(Var_yF)
	new_mu= mu_n + sigma_tilde * new_X
	return (list("mu"=new_mu, "sigma"=new_sigma))
}
dose_allocation <-function(dose){
	var_j=c()
	if(length(dose)==0){
		K=1
	} else{
		K=length(dose)
	}
	#create a sample of M simulated thetas
	theta_sample<-rmvnorm(M, mu_n[[K]], sigma_n[[K]])
	#print(theta_estimate)
	theta_est<-apply(theta_sample, 2, mean)
	for (j in 1:J){
		theta_j_sample<-theta_sample[,j]
		var_jm=c()
		for (m in 1:M){
			y_jm=rnorm(1, theta_j_sample[m], true_sigma)
			temp_res=evolution_eq(mu_n[[K]], sigma_n[[K]], y_jm, j)
			temp_mu=temp_res$mu
			temp_sigma=temp_res$sigma
			temp_theta_sample<-rmvnorm(T, temp_mu, temp_sigma)
			ED95=c()
			for (t in 1:T){
				for (w in 1:J){
					if (temp_theta_sample[t, w]>=0.95*max(temp_theta_sample[t,])){
						ED95=c(ED95, w)
						break
					}
				}
			}
			var_jm=c(var_jm, var(ED95))
		}
		var_j=c(var_j, mean(var_jm))
	}
	return(list("variance"=var_j, "theta"=theta_est))
}
start_time=Sys.time()
for (K in 1:patient){
	if(K==1){
		sigma_n[[K]]=sigma_0
		mu_n[[K]]=mu_0
	}
	#after 30 patients
	else{
		res<-dose_allocation(dose)
		temp_var<-res$variance #call dose allocation to variance vector for every dose
		#theta_estimate<-res$theta
		var_estimate[K-22,]<-temp_var
		theta_estimate[K-22,]<-res$theta
		z_j=max(which(temp_var==min(temp_var)))
		y=c(y, rnorm(1, true_theta[z_j], true_sigma)) #observing and add the true response of the optimal dose
		dose[K-1]<-z_j #add optimal dose to dose vector
		#call a function of update equations for calculating posterior moments, i.e., mu_n, sigma_n
		res=evolution_eq(mu_n[[K-1]], sigma_n[[K-1]], y[length(y)], dose[length(dose)])
		mu_n[[K]]<-res$mu
		sigma_n[[K]]<-res$sigma
		var_target=c(var_target, sigma_n[[K]][6,6]) #for curve1: 6,6 for curve2: 11,11, curve3: 9,9
	}
}
end_time=Sys.time()
start_time-end_time
write(var_target, file="C:/Users/snasrol/Google Drive/Research-Dose finding/Codes/Results-10dose/08.21.2017/KGI-curve1target_var.txt", append=FALSE, sep="\n")
write(dose, file="C:/Users/snasrol/Google Drive/Research-Dose finding/Codes/Results-10dose/08.21.2017/KGI-curve1doses.txt", append=FALSE, sep="\n")
write.table(var_estimate, file="C:/Users/snasrol/Google Drive/Research-Dose finding/Codes/Results-10dose/08.21.2017/KGI-curve1var.txt", sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)
write.table(theta_estimate, file="C:/Users/snasrol/Google Drive/Research-Dose finding/Codes/Results-10dose/08.21.2017/KGI-curve1thetas.txt", sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)
plot(dose[1:patient])
length(which(dose[1:patient]==1))
length(which(dose[1:patient]==2))
length(which(dose[1:patient]==3))
length(which(dose[1:patient]==4))
length(which(dose[1:patient]==5))
length(which(dose[1:patient]==6))
length(which(dose[1:patient]==7))
length(which(dose[1:patient]==8))
length(which(dose[1:patient]==9))
length(which(dose[1:patient]==10))
length(which(dose[1:patient]==11))
