#' Fit of a mixture of two autoregressive models of order 1 to a time series
#'
#' \code{estim_mar2_link} computes the parameters relative to the mixture of the
#' two autoregressive processes of order 1 (MAR, see equation 9 in Qasmi and
#' Ribes, 2021) to capture fast and slow components within a given time series.
#' Parameters are estimated by maximum likelihood via the \code{optim} function.
#'
#' @param y a vector time series corresponding to an estimate of observed
#'     internal variability.
#'
#' @return a vector containing the MAR1 parameters fitted to \code{y}.
#'     \code{var1_ar1} (\code{var2_ar1}) is the variance of the white noise
#'     associated with the first (second) AR1. \code{alpha1_ar1}
#'     (\code{alpha2_ar1}) is the coefficient at lag 1 associated with the first
#'     (second) AR1.
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

  skip_to_algo <- FALSE
  theta = tryCatch(optim(theta0,log_mar,method="BFGS",y=y,control=list(trace=4,REPORT=1)), error = function(e) { skip_to_algo <<- TRUE})
  if(skip_to_algo) { 	theta = optim(theta0,log_mar,method="L-BFGS-B",y=y,control=list(trace=4,REPORT=1)) }
	#theta = optim(theta0,log_mar,method="BFGS",y=y,control=list(trace=4,REPORT=1))
	return(ilinkf_mar2(theta$par))
}

#' Compute the covariance matrix corresponding to a mixture of two
#' autoregressive (MAR) processes of order 1
#'
#' \code{Sigma_mar2} takes the parameters returned by \code{estim_mar2_link}
#' (coefficient at lag 1 and variance) to calculate the matrix coefficients
#' whose formulation is given in equations 2 and 3 in Supplementary Material of
#' Qasmi and Ribes (2021).
#'
#' @param theta a vector containing the MAR parameters returned by
#'     \code{estim_mar2_link}
#' @param n the whished number of lines/columns in the covariance matrix,
#'     usually corresponds to the number of years of observations.
#'
#' @return a symmetric matrix corresponding to the MAR model.
#'
#' @examples
#'
#' @export
Sigma_mar2 = function(theta,n) {
  Sigma1_ar1 = array(NA,dim=c(n,n))
  Sigma2_ar1 = array(NA,dim=c(n,n))
  for (i in 1:n) {
    for (j in 1:n) {
      Sigma1_ar1[i,j] = theta["alpha1_ar1"]^(abs(i-j))
      Sigma2_ar1[i,j] = theta["alpha2_ar1"]^(abs(i-j))
    }
  }
  Sigma = theta["var1_ar1"] * Sigma1_ar1 + theta["var2_ar1"] * Sigma2_ar1
  return(Sigma)
}

Sigma_ar = function(alpha,n) {
  Sigma = array(NA,dim=c(n,n))
  for (i in 1:n) {
    for (j in 1:n) {
      Sigma[i,j] = alpha^(abs(i-j))
    }
  }
  return(Sigma)
}

gauss_log_like = function(y,mu,Sigma) {
	n = length(y)

	if (!prod(dim(Sigma) == c(n,n))) {
		print("Error in gauss_log_like: dim(Sigma)")
		return()
	}

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

linkf_dep = function(theta) {
	return(c(log(theta["v1_slow"]),#
				log(theta["v2_slow"]),#
				tan(theta["a1_slow"]*pi/2),#
				tan(theta["a2_slow"]*pi/2),#
				log(theta["v1_fast"]),#
				log(theta["v2_fast"]),#
				tan(theta["a1_fast"]*pi/2),#
				tan(theta["a2_fast"]*pi/2),#
				tan(theta["lambda"]*pi/2))
			)
}

ilinkf_dep = function(theta_mod) {
	return(c(exp(theta_mod["v1_slow"]),#
				exp(theta_mod["v2_slow"]),#
				atan(theta_mod["a1_slow"])/pi*2,#
				atan(theta_mod["a2_slow"])/pi*2,#
				exp(theta_mod["v1_fast"]),#
				exp(theta_mod["v2_fast"]),#
				atan(theta_mod["a1_fast"])/pi*2,#
				atan(theta_mod["a2_fast"])/pi*2,#
				atan(theta_mod["lambda"])/pi*2)
			)
}

#' Compute the parameter lambda between two MAR by the method of moments
#'
#' Following equation 10 in Supplementary Material of Qasmi and Ribes (2021),
#' \code{estim_mar_dep} computes the parameter lambda by the method of moments
#' given two different MAR parameters.
#'
#' @param theta1 a vector containing the MAR parameters returned by
#'     \code{estim_mar2_link} associated with a first time series (eg residuals
#'     related to global internal variability)
#' @param theta2 a vector containing the MAR parameters returned by
#'     \code{estim_mar2_link} associated with a second time series (eg residuals
#'     related to local internal variability)
#' @param y a vector of the concatenated time series associated with
#'     \code{theta1} and \code{theta2}
#' @param nvar the number of residuals (currently set to 2)
#'
#' @return a vector of parameters containing \code{theta1} and \code{theta2},
#'     and the parameter lambda.
#'
#' @examples
#'
#' @export
estim_mar_dep = function(theta1, theta2, y, nvar=2) {
	# nvar: number of blocks (ie spatial dimension)
	ny = length(y)
	if (nvar != 2) stop("Error in estim_mar_dep: nb is not 2")
	n = ny/nvar
	if (n != round(n)) stop("Error in estim_mar_dep: wrong size of y")
	y1 = y[1:n]
	y2 = y[1:n+n]

	theta = rep(NA,9)
	names(theta) = c("v1_slow","v2_slow","a1_slow","a2_slow","v1_fast","v2_fast","a1_fast","a2_fast","lambda")

	# 1st MAR
	if (theta1["alpha1_ar1"]<theta1["alpha2_ar1"]) {
		# This is the expected case
		theta[c("v1_fast","v1_slow","a1_fast","a1_slow")] = theta1[c("var1_ar1","var2_ar1","alpha1_ar1","alpha2_ar1")]
	} else {
		theta[c("v1_slow","v1_fast","a1_slow","a1_fast")] = theta1[c("var1_ar1","var2_ar1","alpha1_ar1","alpha2_ar1")]
	}

	# 2nd MAR
	if (theta2["alpha1_ar1"]<theta2["alpha2_ar1"]) {
		# This is the expected case
		theta[c("v2_fast","v2_slow","a2_fast","a2_slow")] = theta2[c("var1_ar1","var2_ar1","alpha1_ar1","alpha2_ar1")]
	} else {
		theta[c("v2_slow","v2_fast","a2_slow","a2_fast")] = theta2[c("var1_ar1","var2_ar1","alpha1_ar1","alpha2_ar1")]
	}

	# Fit dependence
	v1_tot = theta["v1_fast"] + theta["v1_slow"]
	v2_tot = theta["v2_fast"] + theta["v2_slow"]
	cov_max = sqrt(theta["v1_fast"]) * sqrt(theta["v2_fast"]) * sqrt(1-theta["a1_fast"]^2) * sqrt(1-theta["a2_fast"]^2) / (1-theta["a1_fast"]*theta["a2_fast"]) + #
				 sqrt(theta["v1_slow"]) * sqrt(theta["v2_slow"]) * sqrt(1-theta["a1_slow"]^2) * sqrt(1-theta["a2_slow"]^2) / (1-theta["a1_slow"]*theta["a2_slow"])
	corr_max = cov_max / sqrt(v1_tot * v2_tot)
	cc = cor(y1,y2)
	theta["lambda"] = max(min(cc/corr_max,1),-1)

	return(theta)

}

#' Compute the parameter lambda between two MAR by maximum likelihood
#'
#' \code{estim_mar_dep_full} computes all the parameters, including lambda,
#' relative to the mixture of the two autoregressive processes of order 1 (MAR,
#' see equation 10 in Supplementary Material of Qasmi and Ribes, 2021) to
#' capture fast and slow components within two given time series. Parameters are
#'  estimated by maximum likelihood via the \code{optim} function.
#'
#' @param y a vector of two concatenated time series, typically corresponding to
#'     global and local internal variability
#'
#' @return a vector containing the 9 MAR parameters fitted to \code{y}.
#'     \code{v1_slow} (\code{v1_fast}) is the variance of the white noise
#'     associated with the slow (fast) component of the first time series.
#'     \code{a1_ar1} (\code{a2_ar1}) is the coefficient at lag 1 associated with
#'     the slow (fast) component of the first time series. Similarly
#'     \code{v2_slow}, \code{v2_fast}, \code{a1_ar1}, \code{a2_ar1}) are the
#'     parameters associated with the second time series. \code{lambda} is the
#'     parameter linking the two time series.
#'
#' @examples
#'
#' @export
estim_mar_dep_full = function(y) {
	# nb: number of blocks (ie spatial dimension)
	ny = length(y)
	nb = 2
	n = ny/nb
	if (n != round(n)) stop("Error in estim_mar_dep: wrong size of y")
	y1 = y[1:n]
	y2 = y[1:n+n]

	theta0 = linkf_dep(c(v1_slow=var(y1)/2,
								v2_slow=var(y2)/2,
								a1_slow=.8,
								a2_slow=.8,
								v1_fast=var(y1)/2,
								v2_fast=var(y2)/2,
								a1_fast=.4,
								a2_fast=.4,
								lambda=.5)
							)

	log_mar_dep = function(theta,y) {
		Sigma = Sigma_mar_dep(ilinkf_dep(theta),n)
		return(gauss_log_like(y,0,Sigma))
	}

	#S = Sigma_mar_dep(ilinkf_dep(theta0),n)
	#save(y,S,theta0, theta00,file="tmp.Rdata")
	#print(gauss_log_like(y,0,S))
	#print(log_mar_dep(theta0,y))
	theta = optim(theta0,log_mar_dep,method="L-BFGS-B",y=y,control=list(trace=4,REPORT=1))
	return(ilinkf_dep(theta$par))
}

#' Compute the covariance matrix corresponding to a mixture of four
#' autoregressive (MAR) processes of order 1 and their link.
#'
#' \code{Sigma_mar_dep} takes the parameters returned by \code{estim_mar_dep} or
#' \code{estim_mar_dep_full} to calculate the matrix coefficients whose
#' formulation is given in equation 10 in Qasmi and Ribes (2021).
#'
#' @param theta a vector containing the MAR parameters returned by
#'     \code{estim_mar_dep} or \code{estim_mar_dep_full}
#'
#' @return a symmetric matrix corresponding to the combination of the two MAR
#'     models. The number of lines/columns in the covariance matrix, usually
#'     corresponds to the length of \code{y}.
#'
#' @examples
#'
#' @export
Sigma_mar_dep = function(theta,n) {
	Sigma = array(NA,dim=c(2*n,2*n))
	# (1,1)-block
	theta_11 = theta[c("v1_slow","a1_slow","v1_fast","a1_fast")]
	names(theta_11) = c("var1_ar1","alpha1_ar1","var2_ar1","alpha2_ar1")
	Sigma[1:n,1:n] = Sigma_mar2(theta_11,n)
	# (2,2)-block
	theta_22 = theta[c("v2_slow","a2_slow","v2_fast","a2_fast")]
	names(theta_22) = c("var1_ar1","alpha1_ar1","var2_ar1","alpha2_ar1")
	Sigma[1:n+n,1:n+n] = Sigma_mar2(theta_22,n)
	# (1,2) and (2,1) blocks
	Sigma [ 1:n, 1:n+n ] = theta["lambda"] * #
		( cov_ar_full_dep(theta["v1_slow"],theta["a1_slow"],theta["v2_slow"],theta["a2_slow"],n) + #
		  cov_ar_full_dep(theta["v1_fast"],theta["a1_fast"],theta["v2_fast"],theta["a2_fast"],n) )
	Sigma [ 1:n+n, 1:n ] = t( Sigma [1:n, 1:n+n] )
	return(Sigma)
}

cov_ar_full_dep = function(v1,a1,v2,a2,n) {
	Cov = array(NA,dim=c(n,n))
	for (i in 1:n) {
		for (j in 1:n) {
			if (i>=j)	Cov[i,j] = a1^(i-j)
			else			Cov[i,j] = a2^(j-i)
		}
	}
	Covfd = sqrt(v1*v2)*sqrt(1-a1^2)*sqrt(1-a2^2)/(1-a1*a2) * Cov
	return(Covfd)
}


