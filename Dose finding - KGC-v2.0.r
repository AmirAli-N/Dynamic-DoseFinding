library(mvtnorm)
##########################################################################
y=c() #response vector
dose=c() #dose vector
J=11 #number of doses	
patient=1000#number of patients
true_sigma=sqrt(10) #true deviation of responses
mu_0=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) #initial hyperparameter
sigma_0=diag(100, nrow=J, ncol=J) #initial hyperparameter
mu_n=lapply(1:patient, function(x) c())
sigma_n=lapply(1:patient, function(x) matrix(0, nrow=J, ncol=J))
theta_estimate=matrix(NA, nrow=patient, ncol=J)
var_estimate=matrix(NA, nrow=patient, ncol=J)
var_target=c()
M=1000
T=1000
n_simulation=30
start_time=0
end_time=0
alpha=1
############################################################################
true_theta=c(0.0, 0.07, 0.18, 0.47, 1.19, 2.69, 5, 7.31, 8.81, 9.53, 9.82)
curve_st="sigmoid-significant"
target_dose=10
############################################################################
evolution_eq <-function(mu, sigma, new_y, z_j){
	e_z = rep(0, J)
	e_z[z_j] <- 1
	sigma_tilde <- (sigma %*% e_z)*(1/sqrt(true_sigma^2+sigma[z_j, z_j]))
	new_sigma <- sigma - sigma_tilde %*% t(sigma_tilde)
	Var_yF <- true_sigma^2 + sigma[z_j, z_j]
	new_X <- (new_y - mu[z_j])/sqrt(Var_yF)
	new_mu= mu + sigma_tilde * new_X
	return (list("mu"=new_mu, "sigma"=new_sigma))
}

dose_allocation <-function(mu, sigma){
	var_j=c()
	#create a sample of M simulated thetas
	theta_sample<-rmvnorm(M, mu, sigma)
	theta_est<-apply(theta_sample, 2, mean)
	var_j=sapply(1:ncol(theta_sample), function(i) {
		var_jm=unlist(lapply(theta_sample[, i], function(y) {	
			y_jm=rnorm(1, y, true_sigma)
			temp_res=evolution_eq(mu, sigma, y_jm, i)
			temp_mu=temp_res$mu
			temp_sigma=temp_res$sigma
			temp_theta_sample<-rmvnorm(T, temp_mu, temp_sigma)
			ED95=c()
			ED95=apply(temp_theta_sample, 1,function(z) {
				if (all(z<=0)){
					return(NA)
				}else{
					return (min(which(z>=0.95*max(z))))
				}
			})
			return(var(ED95, na.rm=TRUE))
		}))
		return(mean(var_jm))
	})
	return(list("variance"=var_j, "theta"=theta_est))
}

start_time=Sys.time()
for (reps in 1:n_simulation){
	set.seed(reps)
	print(paste("reps=", reps))
	for (K in 1:patient){
		if(K==1){
			for (i in 1:J){
				for(j in 1:J){
					sigma_n[[K]][i,j]<-100*exp(-alpha*(i-j)^2)
				}
			}
			mu_n[[K]]=mu_0
		}
	#after 30 patients
		else{
			res<-dose_allocation(mu_n[[K-1]], sigma_n[[K-1]])
			temp_var<-res$variance #call dose allocation to variance vector for every dose
			#theta_estimate<-res$theta
			var_estimate[K-1,]<-temp_var
			theta_estimate[K-1,]<-res$theta
			z_j=min(which(temp_var==min(temp_var)))
			y=c(y, rnorm(1, true_theta[z_j], true_sigma)) #observing and add the true response of the optimal dose
			dose[K-1]<-z_j #add optimal dose to dose vector
			#call a function of update equations for calculating posterior moments, i.e., mu_n, sigma_n
			res=evolution_eq(mu_n[[K-1]], sigma_n[[K-1]], y[length(y)], dose[length(dose)])
			mu_n[[K]]<-res$mu
			sigma_n[[K]]<-res$sigma
			var_target=c(var_target, sigma_n[[K]][target_dose,target_dose])
		}
	}
	write(var_target, file=paste("C:/Results/",curve_st,"/KGC-target_var-",toString(reps),".txt", sep=""), append=FALSE, sep="\n")
	write(dose, file=paste("C:/Results/",curve_st,"/KGC-doses-",toString(reps),".txt", sep=""), append=FALSE, sep="\n")
	write.table(var_estimate, file=paste("C:/Results/",curve_st,"/KGC-var-",toString(reps),".txt", sep=""), sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)
	write.table(theta_estimate, file=paste("C:/Results/",curve_st,"/KGC-thetas-",toString(reps),".txt", sep=""), sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)
	y=c() #response vector
	dose=c() #dose vector
	mu_n=lapply(1:patient, function(x) c())
	sigma_n=lapply(1:patient, function(x) matrix(0, nrow=J, ncol=J))
	theta_estimate=matrix(NA, nrow=patient, ncol=J)
	var_estimate=matrix(NA, nrow=patient, ncol=J)
	var_target=c()
}
end_time=Sys.time()
start_time-end_time
