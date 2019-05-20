library(mvtnorm)
library(mvnmle)
library(MASS)
library(MHadaptive)
library(coda)

y=c() #response vector
dose=c() #dose vector
J=11 #number of doses	
first_stage_patients=J*2
patient=1000#number of patients
true_sigma=10 #true deviation of responses

updated_thetas=matrix(NA, nrow=patient, ncol=2)
updated_sigma=c()

n_simulation=1
start_time=0
end_time=0

par_estimate=matrix(NA, nrow=patient, ncol=3)
var_estimate=matrix(NA, nrow=patient, ncol=J)
############################################################################
true_theta=c(0.0, 0.07, 0.18, 0.47, 1.19, 2.69, 5, 7.31, 8.81, 9.53, 9.82)
curve_st="sigmoid-significant"
target_dose=10
#true_theta=c(0, 0.01, 0.02, 0.05, 0.12, 0.27, 0.5, 0.73, 0.88, 0.95, 0.98)
#curve_st="sigmoid-not-significant"
#target_dose=10
#true_theta=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#curve_st="flat"
#target_dose=5
############################################################################

theta_max=10
log_density<-function(pars, data){
	theta_0=pars[1]
	theta_1=pars[2]
	sigma=pars[3]
	hill_fun<-theta_max/(1+(theta_0/data[,1])^theta_1)
	log_hill_den=sum(dnorm(data[,2], hill_fun, sigma, log=TRUE))
	prior_theta_0=dunif(theta_0, 0, 10, log=TRUE)
	prior_theta_1=dunif(theta_1, -5, 20, log=TRUE)
	prior_sigma=dunif(sigma, 0, 20, log=TRUE)
	prior_den=prior_theta_0+prior_theta_1+prior_sigma
	return(log_hill_den+prior_den)
}

dose_allocation<-function(updated_thetas, updated_sigma, dose, y){
	var_j=c()
	mcmc_r<-Metro_Hastings(li_func=log_density, pars=c(updated_thetas, updated_sigma), par_names=c('theta_0', 'theta_1', 'sigma'), iterations=1000, quiet=TRUE, data=cbind(dose,y))
	mcmc_r<-mcmc_thin(mcmc_r, thin=100)
	theta_est=mcmc_r$trace[,1:2]
	sigma_est=mcmc_r$trace[,3]
	var_j=sapply(1:J, function(i)	{
										var_jm=unlist(sapply(1:length(sigma_est), function(j)	{
																									y_jm=rnorm(1, theta_max/(1+(theta_est[j,1]/i)^theta_est[j,2]), sigma_est[j])
																									y_temp=c(y, y_jm)
																									dose_temp=c(dose, i)
																									mcmc_temp<-Metro_Hastings(li_func=log_density, pars=c(theta_est[j,], sigma_est[j]), par_names=c('theta_0', 'theta_1', 'sigma'), iterations=1000, quiet=TRUE, data=cbind(dose_temp,y_temp))
																									mcmc_temp<-mcmc_thin(mcmc_temp, thin=100)
																									theta_temp=mcmc_temp$trace[,1:2]
																									sigma_temp=mcmc_temp$trace[,3]
																									ED95=c()
																									ED95=sapply(1:length(sigma_temp), function(z)	{
																																						f_z_theta=theta_max/(1+(theta_temp[z,1]/seq.int(1,J,1))^theta_temp[z,2])
																																						if(all(f_z_theta<=0)){
																																							return(NA)
																																						} else{
																																							return(min(which(f_z_theta>=0.95*max(f_z_theta))))
																																						}
																																					})
																									return(var(ED95, na.rm=TRUE))
																								}))
										return(mean(var_jm))
									})
	return(list("variance"=var_j, "par_est"=colMeans(mcmc_r$trace)))							
}

start_time=Sys.time()
for (reps in 1:n_simulation){
	set.seed(reps)
	print(paste("reps=", reps))
	for(i in 1:J){
		y=c(y, rmvnorm(1, true_theta, diag(true_sigma, nrow=J, ncol=J)))
	}
	dose=rep(seq.int(1,J,1), J)
	for (k in 1:patient){
		if (k==1){
			updated_thetas[1,]=c(0,0)
			updated_sigma[1]=10
			res<-dose_allocation(updated_thetas[k,], updated_sigma[k], dose, y)
		} else{
			res<-dose_allocation(updated_thetas[k-1,], updated_sigma[k-1], dose, y)
			if (k%%10==0){
			  print(k)
			}
		}
		temp_var=res$variance
		par_estimate[k,]=res$par_est
		var_estimate[k,]=temp_var
		z_j=min(which(temp_var==min(temp_var)))
		y=c(y, rnorm(1, true_theta[z_j], true_sigma))
		dose=c(dose, z_j)
		updated_thetas[k,]=par_estimate[k, 1:2]
		updated_sigma[k]=par_estimate[k,3]
	}
	write.table(var_estimate, file=paste("C:/Users/snasrol/Google Drive/Research-Dynamic programming to dose-finding clinical trials/Codes/Results-10dose/12.11.2018/",curve_st,"/1000patients-MCMC-est_var-",toString(reps),".txt", sep=""), sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)
	write.table(par_estimate, file=paste("C:/Users/snasrol/Google Drive/Research-Dynamic programming to dose-finding clinical trials/Codes/Results-10dose/12.11.2018/",curve_st,"/1000patients-MCMC-params-",toString(reps),".txt", sep=""), sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)
	y=c()
	dose=c()
	par_estimate=matrix(NA, nrow=patient, ncol=J)
	var_estimate=matrix(NA, nrow=patient, ncol=J)
}
end_time=Sys.time()
start_time-end_time

