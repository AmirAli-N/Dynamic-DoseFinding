library(mvtnorm)
library(mvnmle)
library(MASS)
library(doParallel)
library(Matrix)
library(nicheROVER)
library(matrixcalc)
##########################################################################
J=11 #number of doses	
first_stage_patients=J^2
patient=400#number of patients
true_sigma=c(10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10) #true deviation of responses

M=100
T=100
n_simulation=30
start_time=0
end_time=0
############################################################################
true_theta=c(0.0, 0.07, 0.18, 0.47, 1.19, 2.69, 5, 7.31, 8.81, 9.53, 9.82)
target_dose=10
curve_st="sigmoid"
#true_theta=c(0, 0.01, 0.02, 0.05, 0.12, 0.27, 0.5, 0.73, 0.88, 0.95, 0.98)
#curve_st="sigmoid-not-significant"
#target_dose=10
#true_theta=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#curve_st="flat"
#target_dose=5
############################################################################
evolution_eq <-function(q_n, mu_n, b_n, Beta_n, new_y, z_j){
	e_z=rep(0, J)
	e_z[z_j]=1
	new_q=q_n+(1/J)
	new_b=b_n+(1/J)
	new_mu=mu_n+Beta_n%*%e_z*(new_y-mu_n[z_j])/(Beta_n[z_j, z_j]*(q_n*new_b/(new_b-J+1))+Beta_n[z_j, z_j])
	new_Beta=(new_b/b_n)*Beta_n+(new_b/(b_n+1))*(q_n*(new_y-mu_n[z_j])^2/(q_n*new_b/(new_b-J+1)+1)-Beta_n[z_j, z_j]/b_n)*(Beta_n%*%e_z%*%t(e_z)%*%Beta_n/Beta_n[z_j, z_j]^2)
	return(list("mu"=new_mu, "Beta"=new_Beta))
}
find_ED95 <-function(theta){
	if(all(theta<=0)){
		return(NA)
	} else{
		return(min(which(theta>=0.95*max(theta))))
	}
}
dose_allocation <-function(q_n, mu_n, b_n, Beta_n){
	var_j=c()
	#create a sample of M simulated thetas
	sigma=rwish(M, Beta_n, b_n, inv=TRUE)
	sigma_j=matrix(unlist(lapply(seq(dim(sigma)[3]), function(x) diag(sigma[,,x]))), nrow=M, ncol=J)
	theta=matrix(unlist(lapply(seq(dim(sigma)[3]), function(x) rmvnorm(1, mu_n, (1/q_n)*sigma[,,x], method="chol"))), nrow=J, ncol=J)
	y_jm=matrix(rnorm(M*J, theta, sigma_j), nrow=M, ncol=J, byrow=TRUE)
	temp_res=apply(y_jm, c(1,2), function(x) evolution_eq(q_n, mu_n, b_n, Beta_n, x, which(y_jm==x, arr.ind=TRUE)[2]))
	temp_q=q_n+(1/J)
	temp_mu=matrix(unlist(lapply(temp_res, '[[', 'mu')), nrow=M*J, ncol=J, byrow=TRUE)
	temp_b=b_n+(1/J)
	temp_Beta=lapply(temp_res, '[[', 'Beta')
	temp_sigma=lapply(temp_Beta, function(x) return(tryCatch(rwish(T, x, temp_b, inv=TRUE), error=function(e) NULL)))
	#temp_sigma=compact(temp_sigma)
	temp_theta=lapply(seq(dim(temp_mu)[1]), function(x) lapply(seq(dim(temp_sigma[[x]])[3]), function(y) return(tryCatch(rmvnorm(1, temp_mu[x,], (1/temp_q)*temp_sigma[[x]][,,y]), error=function(e) NULL))))
	ED95=lapply(seq(1:length(temp_theta)), function(x) unlist(lapply(seq(1:length(temp_theta[[x]])), function(y) return(tryCatch(find_ED95(unlist(temp_theta[[x]][[y]])), error=function(e) NULL)))))
	var_ED95=lapply(ED95, function(x) return(tryCatch(var(x, na.rm=TRUE), error=function(x) NULL)))
	var_ED95=as.numeric(as.character(var_ED95))
	var_ED95=matrix(var_ED95, nrow=M, ncol=J, byrow=TRUE)
	var_j=colMeans(var_ED95, na.rm=TRUE)
	return(var_j)
}
############################################################################
############################################################################
rwish<-nicheROVER::rwish
rmvnorm<-mvtnorm::rmvnorm
#compact<-plyr::compact
start_time=Sys.time()
############################################################################
my_cl=makeCluster(detectCores())
registerDoParallel(my_cl)
result=foreach(reps=1:n_simulation, .packages=c("mvtnorm", "mvnmle", "MASS", "Matrix", "nicheROVER", "matrixcalc"), .export=c("evolution_eq", "find_ED95", "dose_allocation", "first_stage_patients", "J", "true_theta", "true_sigma", "patient", "M", "T"))%dopar%{
	set.seed(reps)
	y=c() #response vector
	dose=c() #dose vector
	theta_estimate=matrix(NA, nrow=patient, ncol=J)
	var_estimate=matrix(NA, nrow=patient, ncol=J)
	sigma_estimate=matrix(NA, nrow=patient, ncol=J)
	##initializing with J^2 patients: It seems it is the minimum number to produce positive definite initial matrice B_0
	for(i in 1:(first_stage_patients/J)){
		y=c(y, rmvnorm(1, true_theta, diag(true_sigma, nrow=J, ncol=J), method="chol"))
	}
	M_y=matrix(y, nrow=first_stage_patients/J, ncol=J, byrow=TRUE)
	##maximum likelihood estimate of initial data
	q_0=first_stage_patients/J ##analogous to initial sample size
	mu_0=colMeans(M_y)
	b_0=first_stage_patients/J ##analogous to initial sample size
	Beta_0=matrix(0, nrow=J, ncol=J)
	for(i in 1:(first_stage_patients/J)){
		Beta_0=Beta_0+(M_y[i,]-mu_0)%*%t(M_y[i,]-mu_0)
	}
	if(is.positive.definite(Beta_0)==FALSE){
		Beta_0=nearPD(Beta_0, keepDiag=TRUE, doSym=TRUE, ensureSymmetry=TRUE, conv.tol=0.000001, conv.norm.type="F")
		Beta_0=as.matrix(Beta_0$mat)
	}
	q_n=c()
	mu_n=matrix(NA, nrow=patient, ncol=J)
	b_n=c()
	Beta_n=lapply(1:patient, function(x) matrix(NA, nrow=J, ncol=J))
	#Beta_0=diag(sapply(seq(dim(M_y)[2]), function(x) sum((M_y[,x]-mean(M_y[,x]))^2)))
	for (k in 1:patient){
		ED95.var=0
		if (k==1){
			ED95.var=dose_allocation(q_0, mu_0, b_0, Beta_0)
			q_n[k]=q_0
			b_n[k]=b_0
		} else{
			ED95.var=dose_allocation(q_n[k-1], mu_n[k-1,], b_n[k-1], Beta_n[[k-1]])
		}
		z_j=min(which(ED95.var==min(ED95.var)))
		y=c(y, rnorm(1, true_theta[z_j], true_sigma[z_j])) #observing and add the true response of the optimal dose
		res=list()
		if (k==1){
			res=evolution_eq(q_0, mu_0, b_0, Beta_0, y[k], z_j)
		} else{
			res=evolution_eq(q_n[k-1], mu_n[k-1,], b_n[k-1], Beta_n[[k-1]], y[k], z_j)
		}
		if(k==1){
			q_n[k]=q_0+1/J
			b_n[k]=b_0+1/J
		}else{
			q_n[k]=q_n[k-1]+1/J
			b_n[k]=b_n[k-1]+1/J
		}
		mu_n[k,]=res$mu
		Beta_n[[k]]=as.matrix(res$Beta)
		
		dose[k]=z_j
		var_estimate[k,]=ED95.var
		theta_estimate[k,]<-mu_n[k,]
		sigma_estimate[k,]=diag(Beta_n[[k]])/q_n[k]
	}
	res=list("dose"=dose, "exp_var"=var_estimate, "theta_est"=theta_estimate, "sigma_est"=sigma_estimate)
}
stopCluster(my_cl)
for(i in 1:length(result)){
	write.table(result[[i]]$sigma_est, file=paste("/home/snasrol/R_codes/",curve_st,"-",toString(patient),"-NIW-exp_var-",toString(i),".txt", sep=""), sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)
	write(result[[i]]$dose, file=paste("/home/snasrol/R_codes/",curve_st,"-",toString(patient),"-NIW-doses-",toString(i),".txt", sep=""), append=FALSE, sep="\n")
	write.table(result[[i]]$exp_var, file=paste("/home/snasrol/R_codes/",curve_st,"-",toString(patient),"-NIW-sigma_est-",toString(i),".txt", sep=""), sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)
	write.table(result[[i]]$theta_est, file=paste("/home/snasrol/R_codes/",curve_st,"-",toString(patient),"-NIW-theta_est-",toString(i),".txt", sep=""), sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)
}
############################################################################
end_time=Sys.time()
start_time-end_time