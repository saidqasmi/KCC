#' Load a single object from a Rdata file into a specified variable
#'
#' @param fileName a character string giving the name of the file to load
#'
#' @examples
#'
#' @export
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#' Computes the square root of a matrix
#'
#' \code{matrix_sqrt} calculates the square root of a matrix by diagonalisation.
#'
#' @param M a square matrix
#'
#' @return the square root of \code{M}
#'
#' @export
matrix_sqrt = function(M) {
	Md = eigen(M)
	return(Re(Md$vectors) %*% diag(sqrt(pmax(Re(Md$values),0))) %*% t(Re(Md$vectors)))
}


nnames = function(X) {
	return(names(dimnames(X)))
}

#' Computes an ensemble of realisations from a multivariate gaussian
#' distribution parameters
#'
#' For multiple types of external forcings (natural, all in the current
#' version), \code{mvgauss_to_Xarray} generates a sample of time series
#' corresponding to the response to each forcing from the gaussian parameters.
#'
#' @param mu a vector of time series corresponding to the mean of the
#'     multivariate gaussian distribution. The names of each element must
#'     correspond to the folowing pattern : \code{year_forcing}, eg
#'     \code{1850_nat} for the mean natural response in 1850. The current
#'     allowed forcings are \code{nat} and \code{all}. \code{mu} is tipically a
#'     vector of concatenated time series relative to the response to the
#'     natural and all external forcings for a given period.
#' @param Sigma a matrix of size \code{[length(mu),length(mu)]}
#'     corresponding to the variance-covariance matrix of the multivariate
#'     gaussian distribution
#' @param Nres the size of the wished sample
#'
#' @return returns a 3-D array of dimension
#'     \code{[Nres,length(year),length(forcing)]} where \code{year} is a time
#'     series of years and \code{forcing} is the type of forcing
#'
#' @examples
#'
#' @export
mvgauss_to_Xarray = function(mu,Sigma,Nres) {

	# Dimension and names
	ns = length(mu)
	yf_array = t(simplify2array(strsplit(names(mu),"_")))
	#print(yf_array)
	year_str = unique(yf_array[,1])
	year = year_str
	ny = length(year)
	forcings = unique(yf_array[,2])
	nf = length(forcings)
	sample_str = c("be",paste0("nres",1:Nres))

	Sigma_sqrt = matrix_sqrt(Sigma)

	# Building X
	X = array(NA,dim=c(ny,Nres+1,nf),#
						dimnames=list(year=year,#
										  sample=sample_str,#
										  forc=forcings))
	epsil_cons = array(0,dim = c(ns,Nres+1))
	epsil_cons[,-1] = rnorm(ns*Nres)
	X_tmp = mu%o%rep(1,Nres+1) + Sigma_sqrt %*% epsil_cons

	for (iforc in forcings) {
		indices_iforc = (yf_array[,2]==iforc)
		year_iforc_str = yf_array[indices_iforc,1]
		X[year_iforc_str,,iforc] = X_tmp[indices_iforc,]
	}

	return(X)

}

#' Extract model name from a pseudo observation chain of character
#'
#' For a given character folowing the pattern \code{model_member}, eg
#' \code{CanESM5_r1i1p1f1}, \code{full_names_to_std} returns the character
#' \code{model}, name, i.e. CanESM5.
#'
#' @param Models_fullnames a character, or a vector of character following the
#'     pattern \code{model_member}
#'
#' @return a character, or a vector of character containing the \code{model}
#'     character.
#'
#' @examples
#'
#' @export
full_names_to_std = function(Models_fullnames,sep="_") {
	Nmod = length(Models_fullnames)
	std_names = rep("",Nmod)
	for (i in 1:Nmod) {
		std_names[i] = strsplit(Models_fullnames[i],sep)[[1]][1]
	}
	return(std_names)
}


