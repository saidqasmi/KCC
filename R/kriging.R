#' Derive the posterior distribution of time series from a prior given
#' observations
#'
#' \code{prior2posterior} applies the gaussian conditioning theorem to an array
#' of time series corresponding to the simulated responses to several external
#' forcings (natural, anthropogenic or both) given observations (see equations
#' 12 and 13 in Supplementary Material of Qasmi and Ribes, 2021).
#'
#' @param X_fit a 4-D array of dimension
#'     \code{[length(year), Nres, length(forcing), length(model)]} returned by
#'     \code{x_fit}.
#' @param Xo a vector or a matrix. If a vector, \code{Xo} is a time series of
#'     observations over a given period, and must have names corresponding to
#'     the years of observations (eg \code{1850:2020}). If a matrix, lines
#'     correspond to the years of observations and columns to different type of
#'     observations, which sample measurement uncertainty.
#' @param Sigma_obs a sum of a matrix sampling observed internal variability
#'     (tipically returned by \code{Sigma_mar2}) + a matrix sampling error
#'     measurements in observations (see Methods in Qasmi and Ribes, 2021).
#' @param Nres the whished number of realisations in the posterior gaussian
#'     sample
#' @param centering_CX a logical value indicating whether the constrained time
#'     series must be in anomalies relative to a given period
#'     (see \code{ref_CX})
#' @param ref_CX a vector containing the years corresponding to the reference
#'     period if \code{centering_CX = TRUE}.
#' @param weights a vector of weights of the same length as the number of
#'    models to account for dependencies between models if needed.
#'
#' @return a list of two lists containing the parameters of the unconstrained
#'     (prior) and constrained (posterior) gaussian distributions for the
#'     responses to the forcings in \code{X_fit}. The first (second) list named
#'     \code{uncons} (\code{cons}) contains two other lists, namely \code{mean}
#'     and \code{var}. \code{mean} is a concatenation of time series
#'     corresponding to the mean of the distribution for the different forcings.
#'     The name of each element follows the pattern: \code{year_forcing}, eg
#'     \code{1850_nat} for the mean natural response in 1850. \code{var} is the
#'     covariance matrix associated with \code{mu}, sampling the model
#'     uncertainty.
#'
#' @importFrom abind abind
#'
#' @examples
#'
#' @export
prior2posterior = function(X_fit,Xo,Sigma_obs,Nres=NULL,centering_CX=T,ref_CX="year_obs",S_mean=NULL,Sigma_mod=NULL,weights=NULL) {

	year = dimnames(X_fit)$year

	if (is.null(S_mean) | is.null(Sigma_mod)) {
		# Calculates prior
		S = abind(X_fit[,"be","all",],#
				 X_fit[,"be","nat",],#
				 along=1,use.dnns=T,
				 new.names=list(year=c(paste0(as.character(year),"_all"),#
											  paste0(as.character(year),"_nat")),#
									 NULL))

		if (!is.null(weights)) {
		  S_mean = apply(S, 1, function(x) weighted.mean(x,as.vector(weights)))
		  Sigma_mod = weighted.cov(t(S),as.vector(weights))
		} else {
		  S_mean = apply(S,1,mean)
		  Sigma_mod = var(t(S))
		}

	}

	# constraining
	X_cons = constrain(S_mean,Sigma_mod,Xo,Sigma_obs,Nres,centering_CX=centering_CX,ref_CX=ref_CX)

	return(X_cons)
}

#' Derive the posterior distribution of time series from a prior given
#' observations
#'
#' \code{constrain} applies the gaussian conditioning theorem to a prior
#' distribution for several external forcings (natural, anthropogenic or both)
#' given observations (see equations 12 and 13 in Supplementary Material of
#' Qasmi and Ribes, 2021).
#'
#' @param S_mean a vector accounting for the multi-model mean.
#' @param S_mod a covariance matrix accounting for model uncertainty.
#' @param Xo a vector or a matrix. If a vector, \code{Xo} is a time series of
#'     observations over a given period, and must have names corresponding to
#'     the years of observations (eg \code{1850:2020}). If a matrix, lines
#'     correspond to the years of observations and columns to different type of
#'     observations, which sample measurement uncertainty.
#' @param Sigma_obs a sum of a matrix sampling observed internal variability
#'     (tipically returned by \code{Sigma_mar2}) + a matrix sampling error
#'     measurements in observations (see Methods in Qasmi and Ribes, 2021).
#' @param Nres the whished number of realisations in the posterior gaussian
#'     sample
#' @param centering_CX a logical value indicating whether the constrained time
#'     series must be in anomalies relative to a given period
#'     (see \code{ref_CX})
#' @param ref_CX a vector containing the years corresponding to the reference
#'     period if \code{centering_CX = TRUE}.
#'    models to account for dependencies between models if needed.
#'
#' @return a list of two lists containing the parameters of the unconstrained
#'     (prior) and constrained (posterior) gaussian distributions for the
#'     responses to the forcings in \code{X_fit}. The first (second) list named
#'     \code{uncons} (\code{cons}) contains two other lists, namely \code{mean}
#'     and \code{var}. \code{mean} is a concatenation of time series
#'     corresponding to the mean of the distribution for the different forcings.
#'     The name of each element follows the pattern: \code{year_forcing}, eg
#'     \code{1850_nat} for the mean natural response in 1850. \code{var} is the
#'     covariance matrix associated with \code{mu}, sampling the model
#'     uncertainty.
#'
#' @examples
#'
#' @export
constrain = function(S_mean,Sigma_mod,Xo,Sigma_obs,Nres,centering_CX=T,ref_CX=NULL) {

	# Time axis
	year_obs = dimnames(Xo)$year
	nyox  = length(year_obs)
	locs_in_obs = sub("[0-9][0-9][0-9][0-9]_?","",year_obs)
	locs = unique(locs_in_obs)
	nl = length(locs)

	# Matrix H: selecting "all", centering, extracting relevant years
	S_mean_array = as.array(S_mean)
	dimnames(S_mean_array) = list(year=names(S_mean))
	H_extract = H_extract(S_mean_array,Xo)

	if (centering_CX) {
	  # Initialize ref_CX (if needed)
	  if (is.null(ref_CX) | identical(ref_CX,"year_obs"))	{
	    ref_CX = NULL
	    for (iloc in locs) {
	      is_iloc = ( locs_in_obs == iloc )
	      ref_CX = c(ref_CX,year_obs[is_iloc])
	    }
	  }
	  locs_in_ref_CX = sub("[0-9][0-9][0-9][0-9]_?","",ref_CX)
	  if (!prod(ref_CX %in% year_obs)) {
	    stop("Error in constrain(): Reference period is inconsistent with obs")
	  }
	  if ( !identical(sort(locs),sort(unique(locs_in_ref_CX))) ) {
	    stop("Error in constrain(): Missing ref_CX for some location")
	  }
	  Center = diag(rep(1,nyox))	# a (nyox,nyox) identity matrix
	  for (iloc in locs) {
	    is_iloc = ( locs_in_obs == iloc )
	    nyox_iloc = sum(is_iloc)
	    year_obs_iloc = year_obs[is_iloc]
	    ref_CX_iloc = ref_CX[ locs_in_ref_CX == iloc]
	    Center[is_iloc,is_iloc] = diag(rep(1,nyox_iloc)) - rep(1,nyox_iloc)%o%(year_obs_iloc %in% ref_CX_iloc)/length(ref_CX_iloc)
	  }

	} else {

	  Center	= diag(rep(1,nyox))
	}

	H = Center %*% H_extract			# Centering * extracting

	# Other inputs : x, SX, y, SY
	x	= S_mean
	SX = Sigma_mod
	y  = Center %*% apply(Xo,1,median)
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
	if (sum(grep("_all",year_S_str))) {
	  year_obs_str = paste0(year_obs_str,"_all")
	}
	is_obs_in_year = year_S_str %in% year_obs_str

	if (sum(is_obs_in_year)!=ny_obs) {
		message("Error in H_extract.R: observed years not available in models");
		return()
	}

	H_extract = array(0,dim=c(ny_obs,ny))
	H_extract[,is_obs_in_year] = diag(rep(1,ny_obs))
	return(H_extract)
}

#' Apply the gaussian conditioning theorem to derive a posterior distribution
#' from a prior and observations
#'
#' @param x a vector of time series corresponding to the mean of the prior
#' @param SX the covariance matrix corresponding to \code{x}. The number of
#'     lines and columns must be equal to the length of \code{x}.
#' @param y a vector of time series corresponding to the observations
#' @param SY the covariance matrix corresponding to \code{y}. The number of
#'     lines and columns must be equal to the length of \code{y}.
#' @param H a matrix corresponding to an observation operator. The number of
#'     lines (columns) must be equal to the length of \code{y}(\code{x}).
#'
#' @importFrom MASS ginv
#'
#' @return a list of two lists. The first list \code{mean} contains the mean of
#'     the posterior. The second list \code{var} contains the associated
#'     covariance matrix.
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

weighted.cov = function(x, wt, na.rm = FALSE)
{
  x <- as.matrix(x)
  if (na.rm) {
    x <- na.omit(x)
    wt <- wt[-attr(x, "na.action")]
  }
  wt <- wt/sum(wt)
  mean.x <- colSums(wt * x)
  x <- sqrt(wt) * sweep(x, 2, mean.x, FUN = "-", check.margin = FALSE)
  res <- crossprod(x)/sum(wt)
  return(res)
}


