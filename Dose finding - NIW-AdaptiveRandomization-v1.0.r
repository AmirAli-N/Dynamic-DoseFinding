library(mvtnorm)
library(mvnmle)
library(MASS)
##########################################################################
y=c() #response vector
dose=c() #dose vector
J=11 #number of doses	
first_stage_patients=J*2
patient=1000#number of patients
true_sigma=c(10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10) #true deviation of responses

q_0=0
theta_0=c() #initial hyperparameter
b_0=0
B_0=matrix(NA, nrow=J, ncol=J)

q_n=c()
theta_n=matrix(NA, nrow=patient, ncol=J)
b_n=c()
B_n=lapply(1:patient, function(x) matrix(NA, nrow=J, ncol=J))

theta_estimate=matrix(NA, nrow=patient, ncol=J)
prob_estimate=matrix(NA, nrow=patient, ncol=J)
var_obs=matrix(NA, nrow=patient, ncol=J)
M=1000
T=100
n_simulation=1
start_time=0
end_time=0
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
evolution_eq <-function(q_n, theta_n, b_n, B_n, new_y, z_j){
	e_z = rep(0, J)
	e_z[z_j] <- 1
	new_q=q_n+(1/J)
	new_b=b_n+(1/J)
	new_theta=theta_n+B_n%*%e_z*(new_y-theta_n[z_j])/(B_n[z_j, z_j]*(q_n*new_b)/(new_b-J+1)+B_n[z_j, z_j])
	new_B=(new_b/b_n)*B_n+(new_b/(b_n+1))*(q_n*(new_y-theta_n[z_j])^2/(q_n*new_b/(new_b-J+1)+1)-B_n[z_j, z_j]/b_n)*(B_n%*%e_z%*%t(e_z)%*%B_n/B_n[z_j, z_j]^2)
	new_B=nearPD(new_B, keepDiag=TRUE, doSym=TRUE, ensureSymmetry=TRUE, conv.tol=0.000001, conv.norm.type="F")
	#new_B=(new_B+t(new_B))/2
	return(list("que"=q_n, "theta"=new_theta, "b"=new_b, "B"=new_B))
}
dose_allocation <-function(q_n, theta_n, b_n, B_n){
	var_j=c()
	#create a sample of M simulated thetas
	R=rWishart(M, b_n, B_n)
	mu=matrix(unlist(lapply(seq(dim(R)[3]), function(x)	{
															sigma=ginv(R[,,x])
															sigma=(sigma+t(sigma))/2
															rmvnorm(1, theta_n, (1/q_n)*sigma)
														})), nrow=dim(R)[3], ncol=J, byrow=TRUE)
	theta_est<-apply(mu, 2, mean)
	#var_j=sapply(1:ncol(mu), function(i){
	p_j=sapply(1:ncol(mu), function(i){
											#var_jm=unlist(lapply(mu[ ,i], function(mu_ni)	{
											p_jm=unlist(lapply(mu[ ,i], function(mu_ni)	{
																							y_jm=rt(1, b_n-J+1)*(q_n*(b_n-J+1))/((q_n+1)*B_n[i, i])+mu_ni
																							temp_res<-evolution_eq(q_n, theta_n, b_n, B_n, y_jm, i)
																							temp_q=temp_res$que
																							temp_theta=temp_res$theta
																							temp_b=temp_res$b
																							temp_B=temp_res$B$mat #if comming from the nearPD library
																							#temp_R=rWishart(T, temp_b, round(temp_B, digits=6))
																							temp_R=rWishart(T, temp_b, as.matrix(temp_B))
																							temp_mu=matrix(unlist(lapply(seq(dim(temp_R)[3]), function(j)	{	
																																								sigma=ginv(temp_R[,,j])
																																								sigma=(sigma+t(sigma))/2
																																								rmvnorm(1, temp_theta, (1/q_n)*sigma)
																																							})), nrow=dim(temp_R)[3], ncol=J, byrow=TRUE)
																							ED95=c()
																							ED95=apply(temp_mu, 1, function(z)	{
																																	if (all(z<=0)){
																																		return(NA)
																																	} else{
																																		return(min(which(z>=0.95*max(z))))
																																	}
																																})
																							P=length(which(ED95==i))/length(ED95[!is.na(ED95)])
																							return(P)
																							#return(var(ED95, na.rm=TRUE))
																						}))
											#return(mean(var_jm))
											return(mean(p_jm))
										})
	#return(list("variance"=var_j, "theta"=theta_est))
	return(list("probability"=p_j, "theta"=theta_est))
}
start_time=Sys.time()
for (reps in 1:n_simulation){
	set.seed(reps)
	print(paste("reps=", reps))
	for(i in 1:J){
		y=c(y, rmvnorm(1, true_theta, diag(true_sigma, nrow=J, ncol=J)))
	}
	mle.res<-mlest(matrix(y, nrow=J, ncol=J, byrow=TRUE))
	q_0=J
	theta_0=mle.res$muhat
	b_0=J
	B_0=mle.res$sigmahat
	y=c()
	for (k in 1:patient){
		res=list()
		if (k==1){
			res<-dose_allocation(q_0, theta_0, b_0, B_0)
		} else{
			res<-dose_allocation(q_n[k-1], theta_n[k-1,], b_n[k-1], B_n[[k-1]])
		}
		#temp_var<-res$variance #call dose allocation to variance vector for every dose
		temp_prob<-res$probability
		theta_estimate[k,]<-res$theta
		prob_estimate[k,]<-temp_prob
		z_j=min(which(temp_prob==min(temp_prob)))
		y=c(y, rnorm(1, true_theta[z_j], true_sigma[z_j])) #observing and add the true response of the optimal dose
		dose[k]<-z_j
		if (k==1){
			res<-evolution_eq(q_0, theta_0, b_0, B_0, y[k], dose[k])
		} else{
			res<-evolution_eq(q_n[k-1], theta_n[k-1,], b_n[k-1], B_n[[k-1]], y[k], dose[k])
		}
		q_n[k]=res$que
		theta_n[k,]=res$theta
		b_n[k]=res$b
		B_n[[k]]=as.matrix(res$B$mat)
		
		R=rWishart(M, b_n[k], B_n[[k]])
		var_obs[k,]=diag(rowMeans(array(unlist(lapply(seq(dim(R)[3]), function(x)	{
																						sigma=ginv(R[,,x])
																						simga=(sigma+t(sigma))/2
																						(1/q_n[k])*sigma
																					})), c(J,J,M)), dims=2))
	}

	write.table(var_obs, file=paste("C:/Users/snasrol/Google Drive/Research-Dynamic programming to dose-finding clinical trials/Codes/Results-10dose/12.11.2018/",curve_st,"/1000patients-NIW-obs_var-",toString(reps),".txt", sep=""), sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)
	write(dose, file=paste("C:/Users/snasrol/Google Drive/Research-Dynamic programming to dose-finding clinical trials/Codes/Results-10dose/12.11.2018/",curve_st,"/1000patients-NIW-doses-",toString(reps),".txt", sep=""), append=FALSE, sep="\n")
	write.table(var_estimate, file=paste("C:/Users/snasrol/Google Drive/Research-Dynamic programming to dose-finding clinical trials/Codes/Results-10dose/12.11.2018/",curve_st,"/1000patients-NIW-est_var-",toString(reps),".txt", sep=""), sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)
	write.table(theta_estimate, file=paste("C:/Users/snasrol/Google Drive/Research-Dynamic programming to dose-finding clinical trials/Codes/Results-10dose/12.11.2018/",curve_st,"/1000patients-NIW-thetas-",toString(reps),".txt", sep=""), sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)
	y=c() #response vector
	dose=c() #dose vector
	theta_estimate=matrix(NA, nrow=patient, ncol=J)
	var_estimate=matrix(NA, nrow=patient, ncol=J)
	var_obs=matrix(NA, nrow=patient, ncol=J)
}
end_time=Sys.time()
start_time-end_time