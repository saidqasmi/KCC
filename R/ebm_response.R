#' Calculate nat response for each set of params
#'
#' ANT NAT DECOMPOSITION PART
#' Calculate EBM responses --> Enat
#'
#' @param ebm_params A piControl time series
#' @param FF A piControl time series
#' @param year A piControl time series
#' @param Nres A piControl time series
#'
#' @examples
#'
#' @export
ebm_response = function(FF,ebm_params,year,Nres) {

	n_params_ebm = length(ebm_params[[1]])
	ny = length(year)
	# Calculate nat response for each set of params
	Enat_multi = matrix(NA,ny,n_params_ebm)
	for (i in 1:n_params_ebm){
		Enat_multi_1750 = hmodel(FF$FF$nat,ebm_params$c[i],ebm_params$c0[i],ebm_params$lamb[i],ebm_params$gamm[i])
		Enat_multi[,i] = Enat_multi_1750[FF$FF$year %in% year]
	}
	# Enat -- Best-estimate and random resampling of responses
	Enat = array(NA,dim=c(ny,Nres+1),#
					 dimnames=list(year=year,sample=c("be",paste0("nres",1:Nres))))
	iparam = sample(1:n_params_ebm,Nres,replace=T)
	Enat[,1]  = apply(Enat_multi,1,"mean")
	Enat[,-1] = Enat_multi[,iparam]

	return(Enat)
}

#' Calculate EBM model function
#'
#' More detailed description
#'
#' @param FF A piControl time series
#' @param c A piControl time series
#' @param c0 A piControl time series
#' @param lamb A piControl time series
#' @param gamm A piControl time series
#'
#' @return A detrended time series
#'
#' @examples
#'
#' @export
hmodel <- function(FF, c, c0, lamb, gamm){
  N <- length(FF)
  dt <- 1; #-- timestep (year)
  #-1- Numerical solutions (explicit scheme)
  T <- numeric(N+1);
  To <- numeric(N+1);
  # H <- numeric(N+1);
  T[1] <- 0;
  To[1] <- 0;
  for(n in 2:(N+1)){
    T[n]  <-  (T[n-1] + dt/c*(FF[n-1] - lamb*T[n-1] - gamm*(T[n-1]-To[n-1])));
    To[n] <-  (To[n-1] + dt/c0*(gamm*(T[n-1]-To[n-1])));
    # H[n]  <-  gamm*(T[n] - To[n]);
  }
  # ans <- cbind(T, To, H)
  # colnames(ans) <- c("T", "To", "H")
  # ans
  T[-1]
}
