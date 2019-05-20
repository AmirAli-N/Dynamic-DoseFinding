library(MASS)
library(mvtnorm)
library(coda)
library(dlm)
library(foreach)
library(doSNOW)
#Model specification#######
#Y(t)=F'(t)*mu(t)+e #######
#mu(t)=G(t)*mu(t-1)+ w(t)##
#mu(t)=(theta_t, delta_t)##
###########################
#Variable initialization###
J=11 #number of doses
patient=61 #number of patients
true_theta=c(10, 11, 13, 18, 22, 23, 22, 18, 13, 11, 10)
true_sigma=2
#Parameter initialization##
y=lapply(1:J, function(x) c())
dose=c()
theta_estimate=matrix(NA, nrow=patient, ncol=J)
var_estimate=matrix(NA, nrow=patient, ncol=J)
start_time=0
end_time=0
#MCMC Gibbs sampler
FFBS_sampler <- function(y, iter, thin){
	true_sigma=2
	lambda=0.5
	params=lapply(1:J, function(x) matrix(NA, nrow=1, ncol=2))
	Params=lapply(1:J, function(x) matrix(NA, nrow=iter/thin, ncol=2))
	theta_0=0
	delta_0=0
	m_0=c(theta_0, delta_0) #constructing a vector of parameters
	C_0=diag(100, nrow=2, ncol=2) #constructing a matrix for variances
	G_j=matrix(c(1, 1, 1, 0), nrow=2, ncol=2, byrow=TRUE)
	W=lapply(1:J, function(x) matrix(NA, nrow=2, ncol=2))
	m=lapply(1:J, function(x) matrix(NA, nrow=2, ncol=1))
	C=lapply(1:J, function(x) matrix(NA, nrow=2, ncol=2))
	m[[1]]=m_0
	C[[1]]=C_0
	d=lapply(1:J, function(x) matrix(NA, nrow=2, ncol=1))
	R=lapply(1:J, function(x) matrix(NA, nrow=2, ncol=2))
	B=lapply(1:J, function(x) matrix(NA, nrow=2, ncol=2))
	H=lapply(1:J, function(x) matrix(NA, nrow=2, ncol=2))
	F_j=list()
	D=list()
	Q=list()
	f=list()
	E=list()
	for (j in 1:J){
			F_j[[j]] = matrix(c(1, 0), nrow=2, ncol=length(y[[j]]), byrow=FALSE)
			E_j = diag(1*as.numeric(true_sigma), length(y[[j]]))
			#E[[j]]=matrix(NA, nrow=length(y[[j]]), ncol=length(y[[j]]))
			#W[[j]]=diag(0.5*as.numeric(true_sigma), 2)
			D[[j]]=matrix(NA, nrow=2, ncol=length(y[[j]]))
			Q[[j]]=matrix(NA, nrow=length(y[[j]]), ncol=length(y[[j]]))
			f[[j]]=matrix(NA, nrow=length(y[[j]]), ncol=1)
			#################################################################
			if (j==1){
				d[[j]]=G_j%*%m_0
				W[[j]]=C_0*((1-lambda)/lambda)
				#E[[j]]=W[[j]]
				R[[j]]=G_j%*%C_0%*%t(G_j)+W[[j]]
			} else{
				d[[j]]=G_j%*%m[[j-1]]
				W[[j]]=C[[j-1]]*((1-lambda)/lambda)
				#E[[j]]=W[[j]]
				R[[j]]=G_j%*%C[[j-1]]%*%t(G_j)+W[[j]]
			}
			Q[[j]]=t(F_j[[j]])%*%R[[j]]%*%F_j[[j]]+E_j
			D[[j]]=R[[j]]%*%F_j[[j]]%*%solve(Q[[j]])
			C[[j]]=R[[j]]-D[[j]]%*%Q[[j]]%*%t(D[[j]])
			f[[j]]=t(F_j[[j]])%*%d[[j]]
			m[[j]] = d[[j]]+D[[j]]%*%(y[[j]]-f[[j]])
	}
	for(j in (J-1):1){
		B[[j]]=C[[j]]%*%t(G_j)%*%solve(R[[j+1]])
		H[[j]]=C[[j]]-B[[j]]%*%R[[j+1]]%*%t(B[[j]])
	}
	for (i in 1:iter){
		params[[J]]=rbind(rmvnorm(1, m[[J]], C[[J]], method="svd"))
		if(i%%thin==0){
				Params[[J]][i/thin,]<-params[[J]]
		}
		for (j in (J-1):1){
			#B[[j]]=C[[j]]%*%t(G_j)%*%solve(R[[j+1]])
			h=m[[j]]+B[[j]]%*%(t(params[[j+1]])-d[[j+1]])
			#H=C[[j]]-B[[j]]%*%R[[j+1]]%*%t(B[[j]])
			params[[j]]=rbind(rmvnorm(1, h, H[[j]], method="svd"))
			if(i%%thin==0){
				Params[[j]][i/thin,]<-params[[j]]
			}
		}
	}
	return (list("params"=Params))
}
#dose allocation simulation
dose_allocation <- function(y){
	var_j=list()
	res=FFBS_sampler(y, 100, 1)#100
	res_theta=lapply(res$params, '[',,1)
	K=length(unlist(y))
	theta_estimate[K+1-11,]<-unlist(lapply(res_theta, mean))
	myCluster<-makeCluster(7)#assign six core to the cluster
	registerDoSNOW(myCluster)
	var_j<- foreach(j=1:J, .packages="mvtnorm", .export=c("J", "FFBS_sampler", "true_sigma"), .verbose=FALSE)%:% foreach(m=1:length(res$params[[j]][,1]), .packages="mvtnorm", .combine="c", .export=c("J", "FFBS_sampler", "true_sigma"), .verbose=FALSE)%dopar%{
		y_col=y[[j]]
		y[[j]]=c(y[[j]],rnorm(1, res$params[[j]][,1][m], true_sigma))
		temp_res=FFBS_sampler(y, 1000, 1)#1000
		theta_updated=lapply(temp_res$params,'[',,1)
		sample_size=length(theta_updated[[j]])
		ED95=c()
		theta_matrix=matrix(unlist(theta_updated), nrow=sample_size, ncol=J, byrow=FALSE)
		theta_max_resp<-apply(theta_matrix, 1, max)
		for(t in 1:sample_size){
			if(length(which((theta_matrix[t,]>=0.95*theta_max_resp[t])==TRUE))>=1){
				 min_val=min(theta_matrix[t,][(theta_matrix[t,]>=0.95*theta_max_resp[t])])
				 ED95=c(ED95, which(theta_matrix[t,]==min_val))
			}
		}
		y[[j]]=y_col
		var(ED95)
	}
	stopCluster(myCluster)
	return (list("variance"=as.vector(unlist(lapply(var_j,mean))), "theta"=theta_estimate))
}
start_time=Sys.time()
for (k in 1:patient){
	if (k==1){
		y[[1]]=c(y[[1]], 10)
		y[[2]]=c(y[[2]], 10)
		y[[3]]=c(y[[3]], 10)
		y[[4]]=c(y[[4]], 10)
		y[[5]]=c(y[[5]], 10)
		y[[6]]=c(y[[6]], 10)
		y[[7]]=c(y[[7]], 10)
		y[[8]]=c(y[[8]], 10)
		y[[9]]=c(y[[9]], 10)
		y[[10]]=c(y[[10]], 10)
		y[[11]]=c(y[[11]], 10)
	}
	#other patients
	else{
		####select the best does with minimum variance
		##saving current data before dose allocation
		##dose allocation
		res<-dose_allocation(y)
		temp_var<-res$variance #call dose allocation to variance vector for every dose
		theta_estimate<-res$theta
		var_estimate[k-1,]<-temp_var
		#z_j=which.min(temp_var) #selecting the oprimal dose with minimum variance
		z_j=max(which(temp_var==min(temp_var)))
		##restor data after dose allocation
		##observe the true response of the optimal dose and add it to original data
		y[[z_j]]=c(y[[z_j]], rnorm(1, true_theta[z_j], true_sigma)) #observing and add the true response of the optimal dose
		dose[k-1]<-z_j #add optimal dose to dose vector
	}
}
end_time=Sys.time()
start_time-end_time
write(dose, file="C:/Users/snasrol/Google Drive/Research-Dose finding/Codes/Results-10dose/09.12.2017/NDLM-curve1-2doses.txt", append=FALSE, sep="\n")
write.table(var_estimate, file="C:/Users/snasrol/Google Drive/Research-Dose finding/Codes/Results-10dose/09.12.2017/NDLM-curve1-2var.txt", sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)
write.table(theta_estimate, file="C:/Users/snasrol/Google Drive/Research-Dose finding/Codes/Results-10dose/09.12.2017/NDLM-curve1-2thetas.txt", sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)
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