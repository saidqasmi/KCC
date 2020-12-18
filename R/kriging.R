#' Calculate the prior distribution from X_fit, then applies the constrain, and returns a 4D array, including the "cons" and "uncons" vesrions of X.
#'
#' More detailed description
#'
#' @param x A piControl time series
#'
#' @return A detrended time series
#'
#' @importFrom abind abind
#'
#' @examples
#'
#' @export
prior2posterior = function(X_fit,Xo,Sigma_obs,Nres=NULL,centering_CX=T,ref_CX="year_obs",S_mean=NULL,Sigma_mod=NULL,X_add=NULL) {

	year = dimnames(X_fit)$year

	if (is.null(S_mean) | is.null(Sigma_mod)) {
		# Calculates prior
		S = abind(X_fit[,"be","all",],#
				 X_fit[,"be","nat",],#
				 along=1,use.dnns=T,
				 new.names=list(year=c(paste0(as.character(year),"_all"),#
											  paste0(as.character(year),"_nat")),#
									 NULL))

		if (!is.null(X_add)) {
			names(dimnames(X_add))[1] = "year"
			S = abind(S,X_add,along=1,use.dnns=T)
		}

		S_mean = apply(S,1,mean)
		Sigma_mod = var(t(S))

	}

	# constraining
	X_cons = constrain(S_mean,Sigma_mod,Xo,Sigma_obs,Nres,centering_CX=centering_CX,ref_CX=ref_CX)

	return(X_cons)
}

constrain = function(S_mean,Sigma_mod,Xo,Sigma_obs,Nres,centering_CX=T,ref_CX=NULL) {

	# Time axis
	year_obs = dimnames(Xo)$year
	nyox  = length(year_obs)

	# Matrix H: selecting "all", centering, extracting relevant years
	S_mean_array = as.array(S_mean)
	dimnames(S_mean_array) = list(year=names(S_mean))
	H_extract = H_extract(S_mean_array,Xo)

	if (centering_CX) {

		if (is.null(ref_CX) | ref_CX == "year_obs")	{

			ref_CX = year_obs

		} else if (!prod(ref_CX %in% year_obs)) {

			stop("Error in constrain(): Reference period is inconsistency with obs")
		}

		Center	= diag(rep(1,nyox)) - rep(1,nyox)%o%(year_obs %in% ref_CX)/length(ref_CX)

	} else {

		Center	= diag(rep(1,nyox))
	}

	H = Center %*% H_extract			# Centering * extracting

	# Other inputs : x, SX, y, SY
	x	= S_mean
	SX = Sigma_mod
	y  = Center %*% apply(Xo,1,median)#Xo[,"median"]
	SY = Center %*% Sigma_obs %*% t(Center)

	# List of prior params
	prior_dist = list(mean=S_mean, var=Sigma_mod)

	# Apply constraint
	post_dist = kriging(x,SX,y,SY,H)

	# Put prior_dist and post_dist together
	dist = list(uncons = prior_dist, cons = post_dist)

	return(dist)

}


H_extract = function(S,Xo) {

	# Time axis
	year_S_str = dimnames(S)$year
	year_obs_str = dimnames(Xo)$year
	# Dimensions
	ny = length(year_S_str)
	ny_obs  = length(year_obs_str)

	# Obs from S
	year_obs_all = paste0(year_obs_str,"_all")
	is_obs_in_year = year_S_str %in% year_obs_all

	if (sum(is_obs_in_year)!=ny_obs) {
		message("Error in H_extract.R: observed years not available in models");
		return()
	}

	H_extract = array(0,dim=c(ny_obs,ny))
	H_extract[,is_obs_in_year] = diag(rep(1,ny_obs))
	return(H_extract)
}

#' Computes a set of realisations from a multivariate Gaussian distribution, and organizes the output as a X-array with "usual" dimensions, ie year, sample, forcing
#'
#' More detailed description
#'
#' @param x	: a priori mean				background
#' @param SX	: cov matrix of x				background error
#' @param y	: observation					obs
#' @param SY	: cov matrix of y				observational error
#' @param H	: observation operator		observation operator
#'
#' @importFrom MASS ginv
#'
#' @return X a (ny,Nres+1,nf) array.
#'
#' @examples
#'
#' @export
kriging = function(x, SX, y, SY, H) {

	Sinv		= ginv( H%*%SX%*%t(H) + SY )
	K			= SX %*% t(H) %*% Sinv
	x_post	= x + K %*% (y-H%*%x)
	SX_post	= SX - SX %*% t(H) %*% Sinv %*% H %*% SX

	x_post_out = as.numeric(x_post)
	names(x_post_out) = names(x)
	out = list(mean=x_post_out, var=SX_post)
	return(out)
}


