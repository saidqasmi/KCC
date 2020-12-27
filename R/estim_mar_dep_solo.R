
estim_mar_dep_link = function(y1,y2) {
	n1 = length(y1)
	n2 = length(y2)
	
	if (n1 != n2) {
		print("Error in estim_mar_dep_link(): the two time series must have the same length")
		return()
	} else {
		n = n1
	}
	
	theta0 = linkf_mar_dep(c(var1_ar1=var(y1)/2, var2_ar1=var(y2)/2, alpha1_ar1=.4, alpha2_ar1=.8 ))
							    
	log_mar_dep = function(theta,y1,y2) {
		Sigma = Sigma_mar_dep(ilinkf_mar_dep(theta),length(y1))
		return(gauss_log_like_dep(y1,mean(y1),Sigma))
	}
	
	theta = optim(theta0,log_mar_dep,method="BFGS",y1=y1,y2=y2,control=list(trace=4,REPORT=1))
	return(ilinkf_mar_dep(theta$par))
}


Sigma_mar_dep = function(theta,n) {

  Sigma1_ar1 = array(NA,dim=c(n,n))
  for (i in 1:n) {
    for (j in 1:n) {
    	if (i >= j) {
    		Sigma1_ar1[i,j] = 0
    	} else {
	      Sigma1_ar1[i,j] = theta["alpha1_ar1"]^(abs(i-j))
	     }
    }
  }
  
  Sigma2_ar1 = array(NA,dim=c(n,n))
  for (i in 1:n) {
    for (j in 1:n) {
    	if (i <= j) {
    		Sigma2_ar1[i,j] = 0
    	} else {
				Sigma2_ar1[i,j] = theta["alpha2_ar1"]^(abs(i-j))
	    }
    }
  }
  
  Sigma = (theta["var1_ar1"]*theta["var2_ar1"])*(Sigma1_ar1 + Sigma2_ar1 + diag(n))
  return(Sigma)
}

gauss_log_like_dep = function(y,mu,Sigma) {
	n = length(y)
	
	if (!prod(dim(Sigma) == c(n,n))) {
		print("Error in gauss_log_like_dep(): wrong dimension of Sigma")
		return()
	}
	ll = n*log(Sigma[1,1]) + determinant(Sigma/Sigma[1,1],logarithm=T)$modulus + (y-mu) %*% solve(Sigma) %*% (y-mu)
	return(ll)
}

linkf_mar_dep = function(theta) {
	return(c(log(theta["var1_ar1"]),#
				log(theta["var2_ar1"]),#
				tan(theta["alpha1_ar1"]*pi/2),#
				tan(theta["alpha2_ar1"]*pi/2)))
}
ilinkf_mar_dep = function(theta_mod) {
	return(c(exp(theta_mod["var1_ar1"]),#
				exp(theta_mod["var2_ar1"]),#
				atan(theta_mod["alpha1_ar1"])/pi*2,#
				atan(theta_mod["alpha2_ar1"])/pi*2))
}

