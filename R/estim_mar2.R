#' Detrend time series
#'
#' More detailed description
#'
#' @param y A piControl time series
#'
#' @return A detrended time series
#'
#' @examples
#'
#' @export
estim_mar2_link = function(y) {
	n = length(y)
	theta0 = linkf_mar2(c(var1_ar1=var(y)/2,#
								 var2_ar1=var(y)/2,#
								 alpha1_ar1=.4,#
							    alpha2_ar1=.8 ))
	log_mar = function(theta,y) {
		Sigma = Sigma_mar2(ilinkf_mar2(theta),length(y))
		return(gauss_log_like(y,mean(y),Sigma))
	}
	theta = optim(theta0,log_mar,method="BFGS",y=y,control=list(trace=4,REPORT=1))
	return(ilinkf_mar2(theta$par))
}

#' Detrend time series
#'
#' More detailed description
#'
#' @param y A piControl time series
#'
#' @return A detrended time series
#'
#' @examples
#'
#' @export
Sigma_mar2 = function(theta,n) {
  Sigma1_ar1 = array(NA,dim=c(n,n))
  for (i in 1:n) {
    for (j in 1:n) {
      Sigma1_ar1[i,j] = theta["alpha1_ar1"]^(abs(i-j))
    }
  }
  Sigma2_ar1 = array(NA,dim=c(n,n))
  for (i in 1:n) {
    for (j in 1:n) {
      Sigma2_ar1[i,j] = theta["alpha2_ar1"]^(abs(i-j))
    }
  }
  Sigma = theta["var1_ar1"] * Sigma1_ar1 + theta["var2_ar1"] * Sigma2_ar1
  return(Sigma)
}

gauss_log_like = function(y,mu,Sigma) {
	n = length(y)
	if (!prod(dim(Sigma) == c(n,n))) {print("Error in gauss_log_like: dim(Sigma)"); return()}
	ll = n*log(Sigma[1,1]) + determinant(Sigma/Sigma[1,1],logarithm=T)$modulus + (y-mu) %*% solve(Sigma) %*% (y-mu)
	return(ll)
}

linkf_mar2 = function(theta) {
	return(c(log(theta["var1_ar1"]),#
				log(theta["var2_ar1"]),#
				tan(theta["alpha1_ar1"]*pi/2),#
				tan(theta["alpha2_ar1"]*pi/2)))
}
ilinkf_mar2 = function(theta_mod) {
	return(c(exp(theta_mod["var1_ar1"]),#
				exp(theta_mod["var2_ar1"]),#
				atan(theta_mod["alpha1_ar1"])/pi*2,#
				atan(theta_mod["alpha2_ar1"])/pi*2))
}

