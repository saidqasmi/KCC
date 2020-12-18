#' Loading Rdata with a specific name
#'
#' More detailed description
#'
#' @param fileName a
#'
#' @return As
#'
#' @examples
#'
#' @export
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#' Loading Rdata with a specific name
#'
#' creates an array (/vector/matrix) with only 1 [and appropriate names / dimnames]
#'
#' @param x a
#'
#' @return As
#'
#' @examples
#'
#' @export
ones = function(x) {
  z=x-x+1
  z[is.na(z)]=1
  return(z)
}


matrix_sqrt = function(M) {
	Md = eigen(M)
	return(Re(Md$vectors) %*% diag(sqrt(pmax(Re(Md$values),0))) %*% t(Re(Md$vectors)))
}


nnames = function(X) {
	return(names(dimnames(X)))
}

#' Computes a set of realisations from a multivariate Gaussian distribution, and organizes the output as a X-array with "usual" dimensions, ie year, sample, forcing
#'
#' More detailed description
#'
#' @param mu the mean of MV-Gauss distribution, a vector of size ns, mu needs appropriate names
#' @param Sigma the var-covariance matrix of MV-Gauss distribution, a (ns,ns) matrix
#' @param Nres the size of wished sample
#'
#' @return X a (ny,Nres+1,nf) array.
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


